#!/usr/bin/env python3

import sys
import argparse
import gffutils
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description='Clean gff from spaln')
parser.add_argument('--gff', '-g', help='Input gff from spaln [%(default)s]', default='-')
parser.add_argument('--map', '-m', help='Rename "Target" names using this old-to-new file map')
parser.add_argument('--set-source', '-s', help='Rename source column to [%(default)s]', default='SPALN')
parser.add_argument('--version', '-v', action='version', version='%(prog)s v0.1.0')

def reorder_keys(feature):
    attr = feature.attributes
    kv = dict(zip(attr.keys(), attr.values()))
    attr.clear()
    attr['ID'] = kv.pop('ID')
    if 'Parent' in kv:
        attr['Parent'] = kv.pop('Parent')
    for k in kv:
        attr[k] = kv[k]

if __name__ == '__main__':
    args = parser.parse_args()

    if args.gff == '-':
        fin = sys.stdin
    else:
        fin = open(args.gff)
    gff_db = gffutils.create_db(''.join(fin.readlines()), ":memory:", from_string=True)

    if args.map is not None:
        with open(args.map) as xin:
            xmap = {}
            for line in xin:
                line = line.strip().split()
                assert len(line) == 2
                if line[1] in xmap:
                    sys.stderr.write('Duplicate key in old-to-new file\n')
                    sys.exit(1)
                xmap[line[1]] = line[0]

    print('##gff-version 3')
    for feature in gff_db.all_features():
        if feature.start > feature.end:
            sys.stderr.write('Start greater than end for feature:\n%s\nResetting to start = end\n' % feature)
            feature.start = feature.end
        if feature.featuretype == 'cds':
            feature.featuretype = 'CDS'
        if args.set_source is not None:
            feature.source = args.set_source
        if args.map is not None and 'Target' in feature.attributes:
            target = feature.attributes.pop('Target')
            assert len(target) == 1
            target = target[0].split()
            if target[0] not in xmap:
                raise Exception('%s not found in old-to-new map file' % target[0])
            feature['reference_protein'] = xmap[target[0]]
            feature['reference_aln'] = ' '.join(target[1:])
        # Clean after passing through gffread
        for x in ['geneIDs', 'transcripts', 'geneID', 'Name']:
            if x in feature.attributes.keys():
                feature.attributes.pop(x)
        if 'locus' in feature.attributes.keys():
            parent = feature.attributes.pop('locus')
            feature['Parent'] = parent
        if feature.featuretype == 'locus':
            feature.featuretype = 'gene'
            feature.source = 'gffcl'
        reorder_keys(feature)
        print(feature)
