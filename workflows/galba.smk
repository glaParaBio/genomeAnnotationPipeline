GALBA_VERSION = '1.0.0'

rule galba_install:
    output:
        f'GALBA-{GALBA_VERSION}/scripts/galba.pl',
    params:
        v=GALBA_VERSION,
    shell:
        r"""
        mkdir -p GALBA-{params.v}
        cd GALBA-{params.v}
        rm -rf ./*
        curl -s -L -O https://github.com/Gaius-Augustus/GALBA/archive/refs/tags/v{params.v}.tar.gz
        tar zxf v{params.v}.tar.gz
        rm v{params.v}.tar.gz
        mv GALBA-{params.v}/* ./
        rm -r GALBA-{params.v}
        chmod a+x scripts/*
        """

rule prepare_galba_augustus_config:
    input:
        cfg= os.path.join(INSTALL_PATH, 'augustus/config'),
    output:
        cfg= directory(os.path.abspath('{genome_id}/galba_config')),
    shell:
        # We copy from the conda env the entire config dir but excluding all
        # species except the "generic". Take care of using correct forward
        # slashes!
        r"""
        rsync -ar --exclude='/species/' {input.cfg}/ {output.cfg}
        rsync -ar {input.cfg}/species/generic {output.cfg}/species/
        """


rule galba_run:
    input:
        galba=f'GALBA-{GALBA_VERSION}/scripts/galba.pl',
        augustus= os.path.join(INSTALL_PATH, 'augustus/bin/augustus'),
        genome= '{genome_id}/repeatmasker/{genome_id}.masked.fasta',
        proteins=lambda wc: os.path.abspath(ss[ss.genome_id == wc.genome_id].reference_proteome.iloc[0]),
        cfg= os.path.abspath('{genome_id}/galba_config'),
    output:
        gff= '{genome_id}/galba/augustus.hints.gff3',
        aa= '{genome_id}/galba/augustus.hints.aa',
    shell:
        r"""
        galba_path=`dirname {input.galba}`
        galba_path=`realpath $galba_path`
        PATH=$galba_path:$PATH

        augustus_path=`dirname {input.augustus}`
        PATH=$augustus_path:$PATH
        export PATH

        outdir=`dirname {output.gff}`
        rm -rf $outdir

        rm -rf {input.cfg}/species/{wildcards.genome_id}
        
        galba.pl --species={wildcards.genome_id}_galba --genome={input.genome} \
            --prot_seq={input.proteins} --softmasking --gff3 --threads 8 \
            --workingdir=$outdir \
            --AUGUSTUS_CONFIG_PATH={input.cfg} \
            --AUGUSTUS_SCRIPTS_PATH=`dirname {input.augustus}` \
            --AUGUSTUS_BIN_PATH=`dirname {input.augustus}`
        """
