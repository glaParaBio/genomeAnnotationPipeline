import pandas
import gzip
import os
import shutil
from urllib.parse import urlparse

INSTALL_PATH = os.environ['CONDA_PREFIX']
GENEMARK_PATH = os.path.abspath('genemark')

def is_gzip(filename):
    with gzip.open(filename, 'r') as fh:
        try:
            fh.read(1)
            return True
        except OSError:
            return False

def uri_validator(x):
    # https://stackoverflow.com/questions/7160737/how-to-validate-a-url-in-python-malformed-or-not
    try:
        result = urlparse(x)
        return all([result.scheme, result.netloc])
    except:
        return False

def get_annotation_output(genome_id, source, ext):
    if source == 'augustus.hints':
        return f'{genome_id}/braker/augustus.hints.{ext}'
    elif source == 'spaln':
        return f'{genome_id}/spaln/spaln.{ext}'
    elif source == 'miniprot':
        return f'{genome_id}/miniprot/miniprot.{ext}'
    elif source == 'galba':
        return f'{genome_id}/galba/augustus.hints.{ext}'
    elif source == 'merge':
        return f'{genome_id}/merge/merge.{ext}'
    else:
        raise Exception('Invalid source: %s' % source)

ss = pandas.read_csv(config['sample_sheet'], sep= '\t', comment= '#').dropna(how='all').drop_duplicates()

if len(ss.genome_id) != len(set(ss.genome_id)):
    raise Exception('Duplicate genome IDs found')

if 'generic' in list(ss.genome_id):
    raise Exception('"generic" is a reserved keyword, choose a different genome_id')

wildcard_constraints:
    genome_id = '|'.join([re.escape(x) for x in ss.genome_id]),
    source = '|'.join([re.escape(x) for x in ['augustus.hints', 'spaln', 'miniprot', 'galba', 'merge']]),

rule all:
    input:
        expand('{genome_id}/hmmer/merge.gff3', genome_id= ss.genome_id),

include: 'workflows/installation.smk'
include: 'workflows/galba.smk'

rule clean_sequence_names:
    input:
        genome= lambda wc: ss[ss.genome_id == wc.genome_id].genome_fasta.iloc[0],
    output:
        genome= temp('{genome_id}/repeatmasker/{genome_id}.tmp.fa'),
    run:
        if is_gzip(input.genome):
            xcat = 'gzip -cd'
        else:
            xcat = 'cat'
        cmd = r"""
        {xcat} {input.genome} \
        | awk '{{if($0 ~ "^>") {{
            sub(" .*", "", $0)
        }}
        print $0}}' > {output.genome}
        """
        shell(cmd)


rule build_repeatmasker_db:
    # We build the species databases for RepeatMasker serially by running
    # RepeatMasker on a dummy fasta file. This is to prevent race conditions in
    # case you run the same species on multiple genomes in parallel.
    output:
        rm_database= touch('repeatmasker.db.done'),
    params:
        rm_species= set(ss[pandas.isna(ss.repeatmasker_species) == False].repeatmasker_species)
    run:
        dummy_fa='repeatmasker.dummy.fasta'

        with open(dummy_fa, 'w') as fa:
            fa.write('>dummy\nACTG\n')

        for sp in params.rm_species:
            os.makedirs('tmp_build_repeatmasker_db', exist_ok=True)
            shell(f'RepeatMasker -dir tmp_build_repeatmasker_db -species {sp} {dummy_fa}')
            shutil.rmtree('tmp_build_repeatmasker_db')
            
        os.remove(dummy_fa)


rule mask_genome:
    input:
        genome= '{genome_id}/repeatmasker/{genome_id}.tmp.fa',
        rm_database= 'repeatmasker.db.done',
    output:
        genome= '{genome_id}/repeatmasker/{genome_id}.masked.fasta',
    params:
        rm_species= lambda wc: ss[ss.genome_id == wc.genome_id].repeatmasker_species.iloc[0],
    run:
        outdir = os.path.dirname(output.genome)
        if pandas.isna(params.rm_species):
            shell("cp {input.genome} {output.genome}")
        else:
            shell(r"""
                rm -f {outdir}/{wildcards.genome_id}.tmp.fa.masked
                RepeatMasker -pa 20 -dir {outdir} -species {params.rm_species} -xsmall {input.genome}
                mv {outdir}/{wildcards.genome_id}.tmp.fa.masked {output.genome}
                """)


rule prepare_augustus_config:
    input:
        cfg= os.path.join(INSTALL_PATH, 'augustus/config'),
    output:
        cfg= directory(os.path.abspath('{genome_id}/augustus_config')),
    shell:
        # We copy from the conda env the entire config dir but excluding all
        # species except the "generic". Take care of using correct forward
        # slashes!
        r"""
        rsync -ar --exclude='/species/' {input.cfg}/ {output.cfg}
        rsync -ar {input.cfg}/species/generic {output.cfg}/species/
        """


rule dummy_orthodb_download:
    # If you use a fasta file for training augustus, do not download orthodb
    # and just use this dummy file to signal that downloading is done or not
    # necessary
    output:
        touch('ref/dummy.orthodb'),


rule download_orthodb:
    output:
        db='ref/{db}',
    shell:
        r"""
        curl -L -s -o {output.db} https://data.orthodb.org/download/`basename {output.db}`
        """


def protein_database_input(ss, genome_id):
    # Query the sample sheet to establish whether we need to trigger the
    # download of orthodb. 
    pdb = ss[ss.genome_id == genome_id].protein_database.iloc[0]
    if os.path.isfile(pdb) or uri_validator(pdb):
        return 'ref/dummy.orthodb'
    else:
        orthodb= ['ref/odb11v0_all_fasta.tab.gz', \
                  'ref/odb11v0_level2species.tab.gz', \
                  'ref/odb11v0_levels.tab.gz', \
                  'ref/odb11v0_species.tab.gz']
        return orthodb


rule prepare_protein_database:
    input:
        orthodb= lambda wc: protein_database_input(ss, wc.genome_id),
    output:
        pdb= '{genome_id}/protein_database.fasta',
    params:
        taxonomy= lambda wc: ss[ss.genome_id == wc.genome_id].protein_database.iloc[0],
    run:
        if os.path.isfile(params.taxonomy):
            shell('cp {params.taxonomy} {output.pdb}')
        elif uri_validator(params.taxonomy):
            shell("curl -s -L {params.taxonomy} -o {output.pdb}")
        else:
            shell(r"""
            pigz -cd {input.orthodb[0]} \
            | {workflow.basedir}/scripts/getOrthodbProteinsForTaxonomy.py -f - \
                -l2s {input.orthodb[1]} \
                -l {input.orthodb[2]} \
                -s {input.orthodb[3]} \
                -i {params.taxonomy} > {output.pdb}
                    """)


def get_braker_args(ss, genome_id):
    if 'braker_args' not in ss.columns:
        return ''
    args = ss[ss.genome_id == genome_id].braker_args.iloc[0]
    if pandas.isna(args) or args == '':
        return ''
    else:
        return args

rule braker:
    input:
        augustus= os.path.join(INSTALL_PATH, 'augustus/bin/augustus'),
        braker= os.path.join(INSTALL_PATH, 'bin/braker.pl'),
        genome= '{genome_id}/repeatmasker/{genome_id}.masked.fasta',
        prot_seq= '{genome_id}/protein_database.fasta',
        gm= 'genemark_installation.done',
        cfg= os.path.abspath('{genome_id}/augustus_config'),
    params:
        gmpath= GENEMARK_PATH,
        extra_args= lambda wc: get_braker_args(ss, wc.genome_id),
    output:
        gff= '{genome_id}/braker/augustus.hints.gff3',
        aa= '{genome_id}/braker/augustus.hints.aa',
    shell:
        r"""
        outdir=`dirname {output.gff}`
        rm -rf $outdir
        rm -rf {input.cfg}/species/{wildcards.genome_id}

        braker.pl {params.extra_args} \
            --genome={input.genome} --prot_seq={input.prot_seq} --softmasking --gff3 --cores 16 --species={wildcards.genome_id} \
            --AUGUSTUS_CONFIG_PATH={input.cfg} \
            --AUGUSTUS_BIN_PATH=`dirname {input.augustus}` \
            --AUGUSTUS_SCRIPTS_PATH=`dirname {input.augustus}` \
            --workingdir=$outdir \
            --PROTHINT_PATH={params.gmpath}/ProtHint/bin \
            --GENEMARK_PATH={params.gmpath}
        """


rule merge_annotation:
    input:
        braker='{genome_id}/braker/augustus.hints.gff3',
        miniprot=lambda wc: '{genome_id}/galba/augustus.hints.gff3' if not ss[ss.genome_id == wc.genome_id].reference_proteome.isna().values.any() else [],
    output:
        gff=temp('{genome_id}/merge/merge.gff3'),
    run:
        if input.miniprot == []:
            shell("cp {input.braker} {output.gff}")
        else:
            shell(r"""
            agat_sp_merge_annotations.pl --gff {input.braker} --gff {input.miniprot} --out {output.gff}
            """)


rule extract_proteins:
    input:
        genome=lambda wc: ss[ss.genome_id == wc.genome_id].genome_fasta,
        gff='{genome_id}/merge/merge.gff3',
    output:
        aa='{genome_id}/merge/merge.aa',
    shell:
        r"""
        gffread -y {output.aa} -g {input.genome} {input.gff}
        """


rule download_pfam:
    output:
        pfam= 'ref/Pfam-A.hmm',
        idx= multiext('ref/Pfam-A.hmm', '.h3f', '.h3i', '.h3m', '.h3p'),
    shell:
        r"""
        curl -s -L https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.hmm.gz | gunzip > {output.pfam}
        hmmpress -f {output.pfam}
        """


rule hmmsearch:
    # We want to annotate each protein with Pfam domains. The tools for this is
    # hmmscan. However we use hmmsearch instead since it is much faster and
    # gives the same hits as hmmscan. However, with default settings the
    # evalues of the two programs are different since the databases are
    # different (for hmmscan the database is Pfam and the query is our
    # proteins, vice-versa for hmmsearch). We are interested in the evalue
    # considering the size of Pfam so we set -Z to the size of Pfam.
    input:
        aa= lambda wc: get_annotation_output(wc.genome_id, wc.source, 'aa'),
        pfam= 'ref/Pfam-A.hmm',
        idx= multiext('ref/Pfam-A.hmm', '.h3f', '.h3i', '.h3m', '.h3p'),
    output:
        tab= '{genome_id}/hmmer/{source}.tab',
        dom= '{genome_id}/hmmer/{source}.dom',
    shell:
        r"""
        db_size=`grep '^ACC ' {input.pfam} | wc -l`
        hmmsearch -Z $db_size --cpu 4 --tblout {output.tab} --domtblout {output.dom} --cut_ga {input.pfam} {input.aa} > /dev/null
        """


rule download_pfam2go:
    output:
        pfam2go= 'ref/pfam2go',
    shell:
        r"""
        curl -s -L http://release.geneontology.org/2022-07-01/ontology/external2go/pfam2go > {output.pfam2go}
        """


rule annotate_pfam:
    input:
        gff= lambda wc: get_annotation_output(wc.genome_id, wc.source, 'gff3'),
        dom= '{genome_id}/hmmer/{source}.dom',
        pfam= 'ref/Pfam-A.hmm',
        pfam2go= 'ref/pfam2go'
    output:
        gff= '{genome_id}/hmmer/{source}.gff3',
    shell:
        r"""
        {workflow.basedir}/scripts/add_hmmsearch_to_gff.py --gff {input.gff} \
            --domtblout {input.dom} --hmm {input.pfam} --pfam2go {input.pfam2go} \
        | {workflow.basedir}/scripts/addGeneIdToGff.py -tss '' > {output.gff}
        """
