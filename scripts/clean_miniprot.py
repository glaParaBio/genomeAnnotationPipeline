#!/usr/bin/env python3

import sys
import argparse
import gffutils
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE,SIG_DFL)

def make_id(feature, id_registry):
    p = feature['Parent']
    assert len(p) == 1
    p = p[0]
    if p not in id_registry:
        id_registry[p] = 0
    id_registry[p] += 1
    return p + '.' + str(id_registry[p])

def reorder_keys(feature):
    attr = feature.attributes
    kv = dict(zip(attr.keys(), attr.values()))
    attr.clear()
    attr['ID'] = kv.pop('ID')
    if 'Parent' in kv:
        attr['Parent'] = kv.pop('Parent')
    for k in kv:
        attr[k] = kv[k]

parser = argparse.ArgumentParser(description='Clean gff from miniprot')
parser.add_argument('--gff', '-g', help='Input gff from miniprot [%(default)s]', default='-')
parser.add_argument('--version', '-v', action='version', version='%(prog)s v0.1.0')


MINIPROT_ATTR = """\
# Attributes:
# Rank	int	Rank among all hits of the query
# Identity	real	Fraction of exact amino acid matches
# Positive	real	Fraction of positive amino acid matches
# Donor	str	2bp at the donor site if not GT
# Acceptor	str	2bp at the acceptor site if not AG
# Frameshift	int	Number of frameshift events in alignment
# StopCodon	int	Number of in-frame stop codons
# Target str	Protein coordinate in alignment"""

if __name__ == '__main__':
    args = parser.parse_args()

    if args.gff == '-':
        fin = sys.stdin
    else:
        fin = open(args.gff)
    gff_db = gffutils.create_db(''.join(fin.readlines()), ":memory:", from_string=True)

    for d in gff_db.directives:
        print('##%s' % d)
    print(MINIPROT_ATTR)

    id_registry = {}
    ids = set()

    for feature in gff_db.all_features():
        if 'ID' in feature.attributes:
            assert len(feature['ID']) == 1
            ids.add(feature['ID'][0])
        else: 
            xid = make_id(feature, id_registry)
            if xid in ids:
                sys.stderr.write('ID %s already found in gff. Skipping adding ID\n' % xid)
            else:
                feature['ID'] = xid
                ids.add(xid)
        if 'locus' in feature.attributes.keys():
            parent = feature.attributes.pop('locus')
            feature['Parent'] = parent
        if feature.featuretype == 'locus':
            feature.featuretype = 'gene'

        #if 'Target' in feature.attributes:
        #     target = feature.attributes.pop('Target')
        #     assert len(target) == 1
        #     target = target[0].split()
        #     feature['reference_protein'] = target[0]
        #     feature['reference_aln'] = ' '.join(target[1:])
        reorder_keys(feature)
        print(feature)
