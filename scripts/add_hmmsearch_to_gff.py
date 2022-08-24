#!/usr/bin/env python3

import sys
import re
import argparse

class Hmmer:
    def __init__(self, line):
        self.raw = line.strip()
        xline = line.split()
        self.protein_id = xline[0]
        self.domain_name = xline[2]
        self.domain_accession = xline[3]
        self.best_domain_evalue = float(xline[7])
        self.domain_description = '-'

    def format_hit(self, sep):
        dom_name = self.domain_name.replace('=', '%3D').replace(';', '%3B')
        dom_desc = self.domain_description.replace('=', '%3D').replace(';', '%3B')

        if sep in dom_name or sep in self.domain_accession:
            raise Exception('Value separator "%s" found in output. Choose a different one' % sep)
        out = sep.join([self.domain_accession, dom_name, 'best_dom_evalue:' + str(self.best_domain_evalue), dom_desc])
        return out

def get_desc_from_hmm(hmmfile):
    if hmmfile == '-':
        fin = sys.stdin
    else:
        fin = open(hmmfile)
    acc2desc = {}
    for line in fin:
        if line.startswith('ACC '):
            acc = line.strip().split()
            assert len(acc) == 2
            acc = acc[1]
            assert acc not in acc2desc
            acc2desc[acc] = None
        if line.startswith('DESC '):
            desc = re.sub('^DESC +', '', line.strip())
            assert acc2desc[acc] is None
            acc2desc[acc] = desc
    if hmmfile != '-':
        fin.close()
    return acc2desc

def parse_tblout(tblout, acc2desc, evalue):
    if tblout == '-':
        fin = sys.stdin
    else:
        fin = open(tblout)
    hmmhits = {}
    for line in fin:
        if line.startswith('#'):
            continue
        hit = Hmmer(line)
        if hit.best_domain_evalue < evalue:
            if acc2desc is not None:
                hit.domain_description = acc2desc[hit.domain_accession]
            if hit.protein_id not in hmmhits:
                hmmhits[hit.protein_id] = []
            hmmhits[hit.protein_id].append(hit)
    if tblout != '-':
        fin.close()
    return hmmhits

def get_gff_key(attr_str, key):
    attr = attr_str.rstrip(';').split(';')
    for x in attr:
        if x.startswith(key + '='):
            x = x.replace(key + '=', '')
            return x

def format_hmm_hits(hmmer_list, sep='//'):
    out = []
    for hmmer in hmmer_list:
        fmt = hmmer.format_hit(sep='|')
        if sep in fmt:
            raise Exception('Record %s: Hit separator "%s" found in %s. Choose a different separator' % (hmmer.protein_id, sep, fmt))
        out.append(fmt)
    return sep.join(out)

def add_hmm_to_gff(gff, tblout, newkey):
    if gff == '-':
        fin = sys.stdin
    else:
        fin = open(gff)
    for line in fin:
        line = line.strip()
        if line.startswith('#'):
            print(line)
        else:
            line = line.split('\t')
            attr_str = line[8]
            rec_id = get_gff_key(attr_str,  'ID')
            if rec_id in tblout:
                if get_gff_key(attr_str, newkey) is not None:
                    raise Exception('Key %s already present in gff record' % newkey)
                hmm_hits_fmt = format_hmm_hits(tblout[rec_id])
                attr_str = attr_str.rstrip(';')
                attr_str = attr_str + ';' + newkey + '=' + hmm_hits_fmt + ';'
                line[8] = attr_str
            print('\t'.join(line))

parser = argparse.ArgumentParser(description='Add the output of hmmsearch to gff')
parser.add_argument('--gff', '-gff', help='GFF to annotate', required=True)
parser.add_argument('--key', '-k', help='Attribute key to add [%(default)s]', default='pfam')
parser.add_argument('--tblout', '-t', help='Output of `hmmsearch --tblout <out>` with the annotation', required=True)
parser.add_argument('--hmm', '-H', help='Hmm profile file to extract the domain description. Use - to read from stdin')
parser.add_argument('--evalue', '-e', help='Include domains with e-value below this cutoff [%(default)s]', default=0.01, type=float)
parser.add_argument('--version', '-v', action='version', version='%(prog)s v0.1.0')

if __name__ == '__main__':
    args = parser.parse_args()
    
    if args.hmm is None:
        acc2desc = None
    else:
        acc2desc = get_desc_from_hmm(args.hmm)

    tblout = parse_tblout(args.tblout, acc2desc, args.evalue)
    add_hmm_to_gff(args.gff, tblout, args.key) 
