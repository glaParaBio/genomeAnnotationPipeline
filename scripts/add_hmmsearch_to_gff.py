#!/usr/bin/env python3

import sys
import re
import argparse
import pandas
import gffutils
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE,SIG_DFL)

def transcriptToGenome(db, txid, tx_start, tx_end, featuretype, parent_feature_type, child_feature_type):
    features = list(db.children(txid, featuretype=featuretype, order_by='start'))
    if len(features) == 0:
        raise Exception(f'No feature of type "{featuretype}" found for ID "{txid}"')

    tx = db[txid]
    txToGenomeMap = {}
    out = []
    tx_pos = 1
    for feature in features:
        genome_pos = range(feature.start, feature.end + 1)
        for i in genome_pos:
            # We map tx coordinates to genome coordinates in a dictionary like:
            # {1:11, 2:12, 3:100, 4:101, ...}
            # That's horrible but easier later to query
            txToGenomeMap[tx_pos] = i
            tx_pos += 1
        
        if tx_start < tx_pos:
            if txToGenomeMap[tx_start] < feature.start:
                gstart = feature.start
            else:
                gstart = txToGenomeMap[tx_start] 
            
            if tx_end >= tx_pos:
                gend = feature.end
            else:
                gend = txToGenomeMap[tx_end]
            out.append(gffutils.Feature(tx.chrom, '.', child_feature_type, gstart, gend))
            if tx_end < tx_pos:
                break
    parent = gffutils.Feature(tx.chrom, '.', parent_feature_type, out[0].start, out[-1].end)
    out.insert(0, parent)
    return(out)

def read_tbldom(domfile):
    """Read domain file from hmmer (typically from `--domtblout`) and return a
    DataFrame of selected columns
    """
    if domfile == '-':
        fin = sys.stdin
    else:
        fin = open(domfile)
    dom = fin.readlines()
    if fin != '-':
        fin.close()
    
    mode = [x for x in dom if x.startswith('# Pipeline mode:')][0]
    if 'SEARCH' in mode:
        mode = 'SEARCH'
        colidx = [0, 3, 4]
    elif 'SCAN' in mode:
        mode = 'SCAN'
        colidx = [3, 0, 1]
    else:
        raise Exception('Cannot determine type of hmmer output')
    colidx.extend([12, 17, 18])

    tbldom = []
    USECOLS = pandas.DataFrame({'colidx': colidx, 
        'colname': ['protein_id', 'pfam_name', 'pfam_id', 'dom_evalue', 'aln_start', 'aln_end'], 
        'coltype': ['category', 'category', 'category', 'float64', 'int64', 'int64']})
    for line in dom:
        if line.startswith('#'):
            continue
        line = line.strip()
        line = re.sub(' +', '\t', line).split('\t')
        if len(line) < 23:
            raise Exception('Error parsing line: "%s"' % ' '.join(line))
        line = [line[i] for i in USECOLS.colidx]
        tbldom.append(line)
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

def read_pfam2go(pfam2go, sep='|'):
    """Read pfam2go file (from sequenceontolgy.org) and return a DataFrame
    sep: When a PFAM has multiple GOP terms, join terms using this separator
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
        go_names.append(go_name)

        go_term = line.split(' ; ')[1]
        assert go_term.startswith('GO:')
        go_terms.append(go_term)
    if pfam2go != '-':
        fin.close()
    data = pandas.DataFrame({'pfam_id': pfam_ids, 'go_name': go_names, 'go_term': go_terms})

    godf = pandas.concat([godf, data])
    godf = godf.groupby('pfam_id').agg({'go_term': lambda x: sep.join(x), 
                                        'go_name': lambda x: sep.join(x)}).reset_index()
    return godf

def dom2gff(tbldom, gff_db, source, featuretype, parent_feature_type, child_feature_type):
    """Convert the rows in tbldom to GFF. tbldom is a DataFrame of hmmer hits
    for a single protein. Return these hits as a list of strings (lines) in GFF format. 
    
    gff_db: GFF database from gffutils.create_db()

    featuretype: Splice domains spanning multiple features across this feature
    type. Most likely you want this to be CDS

    source and feature: Strings for columns 2 and 3 of output GFF. For
    "feature" use a valid sequence ontology term
    """
    tbldom = tbldom.reset_index(drop=True)
    assert len(set(tbldom.protein_id)) in [0, 1]
    gff = []
    for idx, row in tbldom.iterrows():
        pfam_desc = row.pfam_desc
        if pandas.isna(pfam_desc):
            pfam_desc = 'NA'
        pfam_desc = pfam_desc
        pfam_name = row.pfam_name
        if pandas.isna(row.go_term):
            ontology_name = 'NA'
            Ontology_term = 'NA'
        else:
            ontology_name = row.go_name
            Ontology_term = row.go_term
        attr = [('ID', row.protein_id + '.d' + str(idx+1)),
                ('Parent', row.protein_id),
                ('Name', pfam_name),
                ('pfam_id', row.pfam_id),
                ('dom_evalue', row.dom_evalue),
                ('signature_desc', pfam_desc),
                ('ontology_name', ontology_name),
                ('Ontology_term', Ontology_term)]
        
        hmm = transcriptToGenome(gff_db, row.protein_id, row.tx_start, row.tx_end, featuretype=featuretype, parent_feature_type=parent_feature_type, child_feature_type=child_feature_type)

        i = 1
        pid = gff_db[row.protein_id]
        for line in hmm:
            line.source = source
            line.strand = pid.strand
            line.score = str(row.dom_evalue)
            for k, v in attr:
                line.attributes[k] = str(v)
            if line.featuretype == parent_feature_type:
                parent_id = line['ID']
            if line.featuretype == child_feature_type:
                line['ID'] = line['ID'][0] + '.' + str(i)
                line['Parent'] = parent_id
                i += 1
            gff.append(line)
    return gff


parser = argparse.ArgumentParser(description='Add the output of hmmsearch or hmmscan to gff. Typical usage: Get proteins from GFF, find protein domains (e.g. from Pfam) using hmmsearch, use this script to annotate the GFF with these domains')
parser.add_argument('--gff', '-gff', help='GFF to annotate', required=True)
parser.add_argument('--domtblout', '-d', help='Domain file from `hmmsearch or hmmscan --domtblout <out>` [%(default)s]', required=True, default='-')
parser.add_argument('--hmm', '-H', help='Optional file of hmm profiles to extract domain descriptions, typically "Pfam-A.hmm". Use - to read from stdin [%(default)s]', default=None)
parser.add_argument('--pfam2go', '-g', help='Optional file mapping Pfam domains to GO terms, typically "pfam2go" from http://www.sequenceontology.org/. Use - to read from stdin [%(default)s]', default=None)
parser.add_argument('--evalue', '-e', help='Include domains with e-value below this cutoff [%(default)s]', default=0.01, type=float)
parser.add_argument('--featuretype', '-f', help='Assign and splice domains to this feature type (column 3 in gff) [%(default)s]', default='CDS')
parser.add_argument('--parent-feature-type', '-p', help='Set this feature type (column 3 in gff) for grouping spliced matches from the same hit [%(default)s]', default='protein_match')
parser.add_argument('--child-feature-type', '-c', help='Set this feature type (column 3 in gff) for spliced matches [%(default)s]', default='protein_hmm_match')
parser.add_argument('--version', '-v', action='version', version='%(prog)s v0.4.0')

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

    # Convert protein coordinates to transcript coordinates 
    tbldom['tx_start'] = tbldom['aln_start'] * 3
    tbldom['tx_end'] = tbldom['aln_end'] * 3

    if args.gff == '-':
        fin = sys.stdin
    else:
        fin = open(args.gff)
    gff_db = gffutils.create_db(''.join(fin.readlines()), ":memory:", from_string=True)

    protein_ids = set(tbldom.protein_id)
    
    print('##gff-version 3')
    for feature in gff_db.all_features():
        print(feature)
        if feature.id in protein_ids:
            hmm_hits = tbldom[tbldom.protein_id.isin([feature.id])]
            pfam_out = dom2gff(hmm_hits, gff_db, featuretype=args.featuretype, source='Pfam', parent_feature_type=args.parent_feature_type, child_feature_type=args.child_feature_type)
            for out in pfam_out:
                print(out)
