import pandas
import gzip
import os

INSTALL_PATH = os.environ['CONDA_PREFIX']

def is_gzip(filename):
    with gzip.open(filename, 'r') as fh:
        try:
            fh.read(1)
            return True
        except OSError:
            return False

GENEMARK_PATH = os.path.abspath('genemark')

ss = pandas.read_csv(config['sample_sheet'], sep= '\t', comment= '#').dropna(how='all').drop_duplicates()

if len(ss.genome_id) != len(set(ss.genome_id)):
    raise Exception('Duplicate genome IDs found')

if 'generic' in list(ss.genome_id):
    raise Exception('"generic" is a reserved keyword, choose a different genome_id')

wildcard_constraints:
    genome_id = '|'.join([re.escape(x) for x in ss.genome_id])


rule all:
    input:
        expand('{genome_id}/hmmer/augustus.hints.gff3', genome_id= ss.genome_id),

include: 'workflows/installation.smk'

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


rule mask_genome:
    input:
        genome= '{genome_id}/repeatmasker/{genome_id}.tmp.fa',
    output:
        genome= '{genome_id}/repeatmasker/{genome_id}.masked.fasta',
    params:
        species= lambda wc: ss[ss.genome_id == wc.genome_id].repeatmasker_species.iloc[0],
    run:
        outdir = os.path.dirname(output.genome)
        if pandas.isna(params.species):
            shell("cp {input.genome} {output.genome}")
        else:
            shell(r"""
                rm -f {outdir}/{wildcards.genome_id}.tmp.fa.masked
                RepeatMasker -pa 20 -dir {outdir} -species {params.species} -xsmall {input.genome}
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


rule download_orthodb:
    output:
        fasta= 'ref/odb10v1_all_fasta.tab.gz',
        l2s= 'ref/odb10v1_level2species.tab.gz',
        level= 'ref/odb10v1_levels.tab.gz',
        species= 'ref/odb10v1_species.tab.gz',
    shell:
        r"""
        cd `dirname {output.fasta}`
        curl -L -s -O https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz
        curl -L -s -O https://v101.orthodb.org/download/odb10v1_levels.tab.gz
        curl -L -s -O https://v101.orthodb.org/download/odb10v1_level2species.tab.gz
        curl -L -s -O https://v101.orthodb.org/download/odb10v1_species.tab.gz
        """


rule dummy_orthodb_download:
    output:
        touch('ref/dummy.orthodb'),


def protein_database_input(ss, genome_id):
    pdb = ss[ss.genome_id == genome_id].protein_database.iloc[0]
    if os.path.isfile(pdb):
        return 'ref/dummy.orthodb'
    else:
        orthodb= ['ref/odb10v1_all_fasta.tab.gz', \
                  'ref/odb10v1_level2species.tab.gz', \
                  'ref/odb10v1_levels.tab.gz', \
                  'ref/odb10v1_species.tab.gz']
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
            shell('cp {params.pdb} {output.pdb}')
        else:
            shell(r"""
            pigz -cd {input.orthodb[0]} \
            | {workflow.basedir}/scripts/getOrthodbProteinsForTaxonomy.py -f - \
                -l2s {input.orthodb[1]} \
                -l {input.orthodb[2]} \
                -s {input.orthodb[3]} \
                -i {params.taxonomy} > {output.pdb}
                    """)


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
    output:
        gff= '{genome_id}/braker/augustus.hints.gff3',
        aa= '{genome_id}/braker/augustus.hints.aa',
    shell:
        r"""
        outdir=`dirname {output.gff}`
        rm -rf $outdir
        rm -rf {input.cfg}/species/{wildcards.genome_id}

        braker.pl --genome={input.genome} --prot_seq={input.prot_seq} --softmasking --gff3 --cores 16 --species={wildcards.genome_id} \
            --AUGUSTUS_CONFIG_PATH={input.cfg} \
            --AUGUSTUS_BIN_PATH=`dirname {input.augustus}` \
            --AUGUSTUS_SCRIPTS_PATH=`dirname {input.augustus}` \
            --workingdir=$outdir \
            --PROTHINT_PATH={params.gmpath}/ProtHint/bin \
            --GENEMARK_PATH={params.gmpath}
        """


rule download_pfam:
    output:
        pfam= temp('ref/Pfam-A.hmm'),
        idx= temp(multiext('ref/Pfam-A.hmm', '.h3f', '.h3i', '.h3m', '.h3p')),
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
        aa= '{genome_id}/braker/augustus.hints.aa',
        pfam= 'ref/Pfam-A.hmm',
        idx= multiext('ref/Pfam-A.hmm', '.h3f', '.h3i', '.h3m', '.h3p'),
    output:
        tab= '{genome_id}/hmmer/augustus.hints.tab',
        dom= '{genome_id}/hmmer/augustus.hints.dom',
    shell:
        r"""
        db_size=`grep '^ACC ' {input.pfam} | wc -l`
        hmmsearch -Z $db_size --cpu 4 --tblout {output.tab} --domtblout {output.dom} --cut_ga {input.pfam} {input.aa} > /dev/null
        """


rule annotate_pfam:
    input:
        tab= '{genome_id}/hmmer/augustus.hints.tab',
        gff= '{genome_id}/braker/augustus.hints.gff3',
        pfam= 'ref/Pfam-A.hmm',
    output:
        gff= '{genome_id}/hmmer/augustus.hints.gff3',
    shell:
        r"""
        {workflow.basedir}/scripts/add_hmmsearch_to_gff.py --gff {input.gff} --tblout {input.tab} --hmm {input.pfam} > {output.gff}
        """
