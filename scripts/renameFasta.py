#!/usr/bin/env python3

import sys
import argparse
from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE,SIG_DFL)

parser = argparse.ArgumentParser(description='Rename sequences in fasta file and create a map of old-to-new names')
parser.add_argument('--input', '-i', help='Input fasta [%(default)s]', default='-')
parser.add_argument('--map', '-m', help='Output file for old-to-new names [%(default)s]', required=True)
parser.add_argument('--version', '-v', action='version', version='%(prog)s v0.1.0')

if __name__ == '__main__':
    args = parser.parse_args()

    if args.input == '-':
        fin = sys.stdin
    else:
        fin = open(args.input)

    fmap = open(args.map, 'w')
    i = 1
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            old = line.lstrip('>').split()[0]
            new = '%06d' % i
            print('>%s' % new)
            fmap.write('%s\t%s\n' % (old, new))
            i += 1
        else:
            print(line)

    fmap.close()
