rule spaln_index:
    input:
        fasta=lambda wc: ss[ss.genome_id == wc.genome_id].genome_fasta,
    output:
        fasta='{genome_id}/spaln/genome.fasta',
        idx=touch('{genome_id}/spaln/genome.idx'),
    shell:
        r"""
        ln -s -f `realpath {input.fasta}` {output.fasta}
        spaln -W -KP -o `dirname {output.fasta}` {output.fasta} 
        """


rule spaln_prepare_query:
    input:
        fasta=lambda wc: os.path.abspath(ss[ss.genome_id == wc.genome_id].reference_proteome.iloc[0]),
    output:
        fasta=temp('{genome_id}/spaln/query.faa'),
        idmap=temp('{genome_id}/spaln/query.map'),
    shell:
        r"""
        {workflow.basedir}/scripts/renameFasta.py -i {input.fasta} -m {output.idmap} > {output.fasta}
        """


rule spaln_aln:
    input:
        fasta='{genome_id}/spaln/query.faa',
        idx='{genome_id}/spaln/genome.idx',
    output:
        gff=temp('{genome_id}/spaln/spaln.tmp'),
    shell:
        r"""
        gff=`realpath {output.gff}`
        qry=`realpath {input.fasta}`
        cd `dirname {input.idx}`
        spaln -t 8 -O:0 -Q7 -d `basename {input.idx}` $qry 2> err.log > $gff
        """


rule spaln_clean_gff:
    input:
        gff='{genome_id}/spaln/spaln.tmp',
        idmap='{genome_id}/spaln/query.map',
    output:
        gff='{genome_id}/spaln/spaln.gff3',
    shell:
        r"""
        awk -v FS='\t' -v OFS='\t' '{{if($3 == "cds") {{sub("cds", "CDS", $3); sub("^ID=", "id=", $9)}} print}}' {input.gff} \
        | gffread -F --cluster-only \
        | sed 's/;id=/;ID=/' \
        | {workflow.basedir}/scripts/clean_spaln.py -m {input.idmap} > {output.gff}
        """


rule spaln_proteins:
    input:
        genome=lambda wc: ss[ss.genome_id == wc.genome_id].genome_fasta,
        gff='{genome_id}/spaln/spaln.gff3',
    output:
        aa=temp('{genome_id}/spaln/spaln.aa'),
    shell:
        r"""
        gffread -y {output.aa} -g {input.genome} {input.gff}
        """

