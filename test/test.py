#!/usr/bin/env python3

import unittest
import sys
import shutil
import os
import subprocess as sp
import pandas
import gffutils

class Test(unittest.TestCase):

    def setUp(self):
        sys.stderr.write('\n' + self.id().split('.')[-1] + ' ') # Print test name
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')
        os.mkdir('test_out')

    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def testShowHelp(self):
        p = sp.Popen("../scripts/add_hmmsearch_to_gff.py --help", shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)

    def testCanSetFeatureType(self):
        p = sp.Popen("zcat data/augustus.hints.gff3.gz | ../scripts/add_hmmsearch_to_gff.py -gff - -d data/augustus.hints.dom -p PROTEIN_MATCH -c PROTEIN_HMM_MATCH", shell=True, stdout=sp.PIPE, stderr=sp.PIPE, executable='/bin/bash')
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        stdout= stdout.decode()
        self.assertTrue('PROTEIN_MATCH' in stdout)
        self.assertTrue('PROTEIN_HMM_MATCH' in stdout)

    def testCanSetParentChildInSplicedMatch(self):
        p = sp.Popen("zcat data/augustus.hints.gff3.gz | ../scripts/add_hmmsearch_to_gff.py -gff - -d data/augustus.hints.dom", shell=True, stdout=sp.PIPE, stderr=sp.PIPE, executable='/bin/bash')
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        stdout= stdout.decode()
        self.assertTrue('ID=g5.t1.d1;Parent=g5.t1;' in stdout)
        self.assertTrue('ID=g5.t1.d1.1;Parent=g5.t1.d1;' in stdout)
        self.assertTrue('ID=g5.t1.d1.2;Parent=g5.t1.d1;' in stdout)

    def testValidGff(self):
        p = sp.Popen("zcat data/augustus.hints.gff3.gz | ../scripts/add_hmmsearch_to_gff.py -gff - -d data/augustus.hints.dom -g <(zcat data/pfam2go.gz)", shell=True, stdout=sp.PIPE, stderr=sp.PIPE, executable='/bin/bash')
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        stdout= stdout.decode()
        gff = gffutils.create_db(stdout, ':memory:', from_string=True, merge_strategy='error', force_gff=True) 

    def testCanAddOntology(self):
        p = sp.Popen("zcat data/augustus.hints.gff3.gz | ../scripts/add_hmmsearch_to_gff.py -gff - -d data/augustus.hints.dom -g <(zcat data/pfam2go.gz)", shell=True, stdout=sp.PIPE, stderr=sp.PIPE, executable='/bin/bash')
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        stdout= stdout.decode()
        self.assertTrue('Name=TBCA;pfam_id=PF02970;dom_evalue=2.4e-18;signature_desc=NA;ontology_name=beta-tubulin binding|tubulin complex assembly|post-chaperonin tubulin folding pathway;Ontology_term=GO:0048487|GO:0007021|GO:0007023' in stdout)

if __name__ == '__main__':
    unittest.main()
