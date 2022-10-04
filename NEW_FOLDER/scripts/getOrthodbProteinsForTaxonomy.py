#!/usr/bin/env python3

import argparse
from itertools import groupby
import sys
import pandas
import numpy

def fasta_iter(fasta_name):
    """
    From https://www.biostars.org/p/710/
    Given a fasta file. yield tuples of header, sequence
    """
    if fasta_name == '-':
        fh = sys.stdin
    else:
        fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (headerStr, seq)

def get_orthodb(level2species_file, levels_file, species):
    levels = pandas.read_csv(levels_file, sep= '\t', names= ['taxid', 'tax_name'], usecols= [0, 1], dtype= {'taxid': int, 'tax_name': str})
    l2s = pandas.read_csv(level2species_file, sep= '\t', names= ['odbid', 'taxid_list'], usecols= [1, 3])
    odbid = []
    taxid = []
    for idx,row in l2s.iterrows():
        taxids = row.taxid_list.strip()
        assert taxids.startswith('{') and taxids.endswith('}')
        taxids = taxids.lstrip('{').rstrip('}').split(',')
        taxids = [int(x) for x in taxids]
        assert "_" in row.odbid
        for x in taxids:
            taxid.append(x)
            odbid.append(row.odbid)
    l2s = pandas.DataFrame({'odbid': odbid, 'taxid': taxid})
    N = len(l2s)
    l2s = pandas.merge(l2s, levels, how= 'left', on= 'taxid')
    assert len(l2s) == N
    l2s = pandas.merge(l2s, species[['taxid', 'species_name']].drop_duplicates(), how= 'left', on= 'taxid')
    assert len(l2s) == N
    l2s['tax_name'] = numpy.where(l2s.tax_name.isna(), l2s['species_name'], l2s['tax_name'])
    l2s.drop(columns= 'species_name', inplace= True)

    tax_name = [str(x).lower() for x in list(l2s.tax_name)]
    l2s['tax_name'] = tax_name
    return l2s

parser = argparse.ArgumentParser(description= r"""Query OrthoDB files to extract proteins in taxonomic groups.
For input files see https://www.orthodb.org/?page=filelist""", formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('--odb-fasta', '-f', required= True, default= '-', 
    help= "Fasta file of proteins, typically odb10v1_all_fasta.tab [%(default)s]")

parser.add_argument('--odb-level2species', '-l2s', required= True, 
    help= "File mapping Ortho DB organism ids based on NCBI taxonomy ids, tipically odb10v1_level2species.tab")

parser.add_argument('--odb-levels', '-l', required= True, 
    help= "File of NCBI taxonomy nodes (levels), tipically odb10v1_levels.tab")

parser.add_argument('--odb-species', '-s', required= True, 
    help= "File of orthodb species and names, tipically odb10v1_species.tab")

parser.add_argument('--include', '-i', required= True, nargs= '+', 
    help= "One or more taxonomic names or IDs to get proteins for. Names are case insensitive must match completely E.g. 'plasmdium berghei anka'")

parser.add_argument('--exclude', '-x', nargs= '*', default= [],
    help= "Optional list of taxonomic names or IDs to drop from the 'include' set")

parser.add_argument('--version', '-v', action= 'version', version= '%(prog)s v0.1.0')

# ------------ #

args = parser.parse_args()

species = pandas.read_csv(args.odb_species, sep= '\t', names= ['taxid', 'odbid', 'species_name'], usecols= [0, 1, 2]).drop_duplicates()
orthodb = get_orthodb(args.odb_level2species, args.odb_levels, species)

keep_odbid = set()
for x in set(args.include):
    if x.isdigit():
        x = int(x)
        keep = set(orthodb[orthodb.taxid == x].odbid)
    else:
        x = x.lower()
        keep = set(orthodb[orthodb.tax_name == x].odbid)
    if len(keep) == 0:
        sys.stderr.write('Taxonomy ID "%s" not found\n' % x)
        sys.exit(1)
    keep_odbid.update(keep)

drop_odbid = set()
for x in set(args.exclude):
    if x.isdigit():
        x = int(x)
        drop = set(orthodb[orthodb.taxid == x].odbid)
    else:
        x = x.lower()
        drop = set(orthodb[orthodb.tax_name == x].odbid)
    if len(drop) == 0:
        sys.stderr.write('Taxonomy ID "%s" not found\n' % x)
        sys.exit(1)
    drop_odbid.update(drop)

keep_odbid = keep_odbid - drop_odbid
assert all([x in set(species.odbid) for x in keep_odbid])

keep = species[species.odbid.isin(keep_odbid)].sort_values('species_name')
sys.stderr.write(keep.to_markdown(index= False) + '\n')
sys.stderr.write('%s species\n' % len(keep))

fiter = fasta_iter(args.odb_fasta)

found = set()
for ff in fiter:
    header, seq = ff
    odbid = header.split(':')
    assert len(odbid) == 2
    odbid = odbid[0]
    if odbid in keep_odbid:
        found.add(odbid)
        print('>%s' % header)
        print(seq)

missing = keep_odbid - found
if len(missing) > 0:
    sys.stderr.write('Error: Could not get sequence for %s OrthoDB:\n%s\n' % (len(missing), ', '.join(missing)))
