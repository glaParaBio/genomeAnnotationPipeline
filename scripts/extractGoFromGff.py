#!/usr/bin/env python3

import sys
import argparse
import gffutils
import pandas

from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

EVIDENCE_CODE = 'ILS' # See http://geneontology.org/docs/guide-go-evidence-codes/

parser = argparse.ArgumentParser(description='Extract GO nnotation from GFF. Output is suitable for Bioconductor::AnnotationDbi')

parser.add_argument('gff', type= str, help= 'Input GFF file [%(default)s]', default= '-', nargs= '?')

parser.add_argument('--go-key', '-g', type= str, help= 'Attribute key holding the GO terms [%(default)s]', default='Ontology_term')
parser.add_argument('--gene-id', '-id', type= str, help= 'Attribute key of the gene identifier [%(default)s]', default= 'gene_id')
parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

args = parser.parse_args()

if args.gff == '-':
    gff = sys.stdin.readlines()
else:
    gff = open(args.gff) # .readlines()

db = gffutils.create_db(''.join(gff), ':memory:', from_string= True, merge_strategy= 'create_unique')

header = ['#gene_id', 'go_term', 'evidence_code']
data = []
for feature in db.all_features():
    if args.go_key in feature.attributes.keys():
        if args.gene_id not in feature.attributes.keys():
            sys.stderr.write('Gene ID key "%s" not found in record:\n%s\n' % (args.gene_id, feature))
            sys.exit(1)
        gene_id = feature.attributes[args.gene_id]
        assert len(gene_id) == 1
        gene_id = gene_id[0]
        go = feature.attributes[args.go_key]
        for g in go:
            data.append([gene_id, g, EVIDENCE_CODE])

data = pandas.DataFrame(data, columns=header).drop_duplicates()
data.to_csv(sys.stdout, sep='\t', index=False)
