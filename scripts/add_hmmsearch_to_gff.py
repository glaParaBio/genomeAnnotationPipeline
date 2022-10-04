#!/usr/bin/env python3

import sys
import re
import argparse
import pandas

def read_tbldom(domfile):
    """Read domain file from hmmer (typically from `--domtblout`) and return a
    DataFrame of selected columns
    """
    if domfile == '-':
        fin = sys.stdin
    else:
        fin = open(domfile)
    tbldom = []
    USECOLS = pandas.DataFrame({'colidx': [0, 3, 4, 12, 17, 18], 
        'colname': ['protein_id', 'pfam_name', 'pfam_id', 'dom_evalue', 'aln_start', 'aln_end'], 
        'coltype': ['category', 'category', 'category', 'float64', 'int64', 'int64']})
    for line in fin:
        if line.startswith('#'):
            continue
        line = line.strip()
        line = re.sub(' +', '\t', line).split('\t')
        if len(line) != 23:
            raise Exception('Error parsing line: "%s"' % ' '.join(line))
        line = [line[i] for i in USECOLS.colidx]
        tbldom.append(line)
    if fin != '-':
        fin.close()
    tbldom = pandas.DataFrame(tbldom, columns=USECOLS.colname)
    tbldom = tbldom.astype(dict(zip(USECOLS.colname, USECOLS.coltype)))
    tbldom.pfam_id = [x.split('.')[0] for x in tbldom.pfam_id]
    if len(tbldom) > 0:
        assert len(tbldom[tbldom.dom_evalue < 0]) == 0
    return tbldom

def get_desc_from_hmm(hmmfile):
    """Read profile file (typically Pfam-A.hmm) to extract PFAM accession and
    description. Return a DataFrame
    """
    accdf = pandas.DataFrame({'pfam_id': [], 'pfam_desc': []})
    accdf = accdf.astype({'pfam_id': 'category', 'pfam_desc': 'str'})
    if hmmfile is None:
        return accdf
    elif hmmfile == '-':
        fin = sys.stdin
    else:
        fin = open(hmmfile)
    acc2desc = {}
    for line in fin:
        if line.startswith('ACC '):
            acc = line.strip().split()
            assert len(acc) == 2
            acc = acc[1]
            acc = acc.split('.')[0]
            assert acc not in acc2desc
            acc2desc[acc] = None
        if line.startswith('DESC '):
            desc = re.sub('^DESC +', '', line.strip())
            assert acc2desc[acc] is None
            acc2desc[acc] = desc
    if hmmfile != '-':
        fin.close()
    accdf.pfam_id = list(acc2desc.keys())
    accdf.pfam_desc = list(acc2desc.values())
    return accdf

def gff_encode(x):
    """Replace reserved characters in `x` with URL encoding to make `x`
    gff-compatible
    """
    x = re.sub(';', '%3B', x)
    x = re.sub('=', '%3D', x)
    x = re.sub(',', '%2C', x)
    return x

def read_pfam2go(pfam2go):
    """Read pfam2go file (from sequenceontolgy.org) and return a DataFrame
    """
    godf= pandas.DataFrame({'pfam_id': [], 'go_name': [], 'go_term': []}, dtype='category')
    if  pfam2go is None:
        return godf
    elif pfam2go == '-':
        fin = sys.stdin
    else:
        fin = open(pfam2go)

    pfam_ids = []
    go_names = []
    go_terms = []
    for line in fin:
        if line.startswith('!'):
            continue
        line = line.strip()
        pfam_id = line.split()[0]
        assert pfam_id.startswith('Pfam:')
        pfam_id = pfam_id.split(':')[1]
        pfam_ids.append(pfam_id)

        go_name = line.split(' > ')[1].split(' ; ')[0].strip()
        assert go_name.startswith('GO:')
        go_name = re.sub('^GO:', '', go_name)
        go_name = gff_encode(go_name)
        go_names.append(go_name)

        go_term = line.split(' ; ')[1]
        assert go_term.startswith('GO:')
        go_terms.append(go_term)
    if pfam2go != '-':
        fin.close()
    data = pandas.DataFrame({'pfam_id': pfam_ids, 'go_name': go_names, 'go_term': go_terms})

    godf = pandas.concat([godf, data])
    godf = godf.groupby('pfam_id').agg({'go_term': lambda x: ','.join(x), 
                                        'go_name': lambda x: ','.join(x)}).reset_index()
    return godf

def get_gff_key(gff_line, key):
    """Extract attribute `key` from gff line
    """
    attr_str = gff_line.split('\t')[8]
    attr = attr_str.rstrip(';').split(';')
    for x in attr:
        if x.startswith(key + '='):
            x = x.replace(key + '=', '')
            return x

def dom2gff(tbldom, chrom, parent_start, parent_end, parent_strand, source='Pfam', feature='protein_match'):
    """Convert the rows in tbldom to GFF. tbldom is a DataFrame of hmmer hits
    for a single protein. Return these hits as a list of strings (lines) in GFF format. 
    
    chrom, parent_start, parent_end, and parent_strand: Genomic coordinates of
    the transcript producing the protein in question.
    
    source and feature: Strings for columns 2 and 3 of output GFF. For
    "feature" use a valid sequence ontology term
    """
    tbldom = tbldom.reset_index(drop=True)
    assert len(set(tbldom.protein_id)) in [0, 1]
    gff = []
    for idx, row in tbldom.iterrows():
        if parent_strand == '+':
            aln_start = parent_start + (row.aln_start - 1) * 3
            aln_end = (parent_start + (row.aln_end - 1) * 3) - 1
        elif parent_strand == '-':
            aln_start = parent_end - (row.aln_end - 1) * 3 + 1
            aln_end = parent_end - (row.aln_start - 1) * 3
        assert aln_start < aln_end
        assert aln_start >= parent_start
        assert aln_end <= parent_end
        pfam_desc = row.pfam_desc
        if pandas.isna(pfam_desc):
            pfam_desc = 'NA'
        pfam_desc = gff_encode(pfam_desc)
        pfam_name = gff_encode(row.pfam_name)
        if pandas.isna(row.go_term):
            go_ann = ''
        else:
            go_ann = f';ontology_name={row.go_name};Ontology_term={row.go_term}'
        attr = f'ID={row.protein_id}.d{idx+1};Parent={row.protein_id};Name={pfam_name};pfam_id={row.pfam_id};dom_evalue={row.dom_evalue};signature_desc={pfam_desc}{go_ann}'
        gff_row = [chrom, source, feature, aln_start, aln_end, row.dom_evalue, parent_strand, '.', attr]
        gff_row = '\t'.join([str(x) for x in gff_row])
        gff.append(gff_row)
    return gff


parser = argparse.ArgumentParser(description='Add the output of hmmsearch to gff. Typical usage: Get proteins from GFF, find protein domains (e.g. from Pfam) using hmmsearch, use this script to annotate the GFF with these domains')
parser.add_argument('--gff', '-gff', help='GFF to annotate', required=True)
parser.add_argument('--domtbl', '-d', help='Domain from `hmmsearch --domtblout <out>` [%(default)s]', required=True, default='-')
parser.add_argument('--hmm', '-H', help='Optional file of hmm profiles to extract domain descriptions, typically "Pfam-A.hmm". Use - to read from stdin [%(default)s]', default=None)
parser.add_argument('--pfam2go', '-g', help='Optional file mapping Pfam domains to GO terms, typically "pfam2go" from http://www.sequenceontology.org/. Use - to read from stdin [%(default)s]', default=None)
parser.add_argument('--evalue', '-e', help='Include domains with e-value below this cutoff [%(default)s]', default=0.01, type=float)
parser.add_argument('--version', '-v', action='version', version='%(prog)s v0.2.0')

if __name__ == '__main__':
    args = parser.parse_args()
    
    pfam2go = read_pfam2go(args.pfam2go)
    acc2desc = get_desc_from_hmm(args.hmm)
    tbldom = read_tbldom(args.domtblout)

    tbldom = tbldom[tbldom.dom_evalue < args.evalue]
    n = len(tbldom)
    tbldom = tbldom.merge(acc2desc, how='left', on='pfam_id')
    tbldom = tbldom.merge(pfam2go, how='left', on='pfam_id')
    assert n == len(tbldom)
    
    if args.gff == '-':
        fin = sys.stdin
    else:
        fin = open(args.gff)

    protein_ids = set(tbldom.protein_id)
    first = True
    for line in fin:
        if first and not line.startswith('##gff-version'):
            print('##gff-version 3')
            first = False
        if line.startswith('#'):
            print(line)
            continue
        print(line.strip())
        gff_id = get_gff_key(line, 'ID')
        if gff_id in protein_ids:
            hmm_hits = tbldom[tbldom.protein_id.isin([gff_id])]
            gff_line = line.split('\t')
            pfam_out = dom2gff(hmm_hits, chrom=gff_line[0], parent_start=int(gff_line[3]), parent_end=int(gff_line[4]), parent_strand=gff_line[6])
            for out in pfam_out:
                print(out)
