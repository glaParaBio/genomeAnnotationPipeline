rule miniprot_aln:
    input:
        ref=lambda wc: ss[ss.genome_id == wc.genome_id].genome_fasta,
        query=lambda wc: os.path.abspath(ss[ss.genome_id == wc.genome_id].reference_proteome.iloc[0]),
    output:
        gff=temp('{genome_id}/miniprot/miniprot.gff3'),
        paf='{genome_id}/miniprot/miniprot.paf',
    shell:
        r"""
        echo '#prot_name prot_len prot_start prot_end prot_strand ref_name ref_len ref_start ref_end n_match_nt n_match_nt_excl_introns mapq aln_score aln_score_excl_introns n_aa_pos_score dist_start_codon dist_stop_codon cigar diff_string' \
        | tr ' ' '\t' > {output.paf}

        miniprot -u -t 8 --gff -N 1 {input.ref} {input.query} \
        | awk -v FS='\t' -v OFS='\t' '{{
            if($0 ~ "^##PAF\t") {{
                sub("^##PAF\t", "", $0)
                print >> "{output.paf}"
            }} else {{
                print $0
            }}
        }}' \
        | gffread -F --cluster-only \
        | {workflow.basedir}/scripts/clean_miniprot.py > {output.gff}
        """


rule miniprot_proteins:
    input:
        genome=lambda wc: ss[ss.genome_id == wc.genome_id].genome_fasta,
        gff='{genome_id}/miniprot/miniprot.gff3',
    output:
        aa=temp('{genome_id}/miniprot/miniprot.aa'),
    shell:
        r"""
        gffread -y {output.aa} -g {input.genome} {input.gff}
        """

