rule collapse_reads:
    group: lambda wildcards: wildcards.sample
    input:
        expand(f"{cwd}/reads/{{sample}}_trimmed_final.fq.gz", sample=sample_names.keys()),
        trim_read=lambda wildcards: f"{cwd}/reads/{wildcards.sample}_trimmed_final.fq.gz"
    output:
        collapsed=f"{cwd}/reads/{{sample}}_trimmed_collapsed.fa"
    params:
        umi_p=umi_param
    log:
        err_log = f"{cwd}/logs/collapsing_reads/{{sample}}.err"
    shell:
        """
        seqtk seq -a {input.trim_read} | fastx_collapser -i - {params.umi_p} | sed 's/>/>Read/g' > {output.collapsed} 2>> {log.err_log}
        status=$?
        if [ $status -ne 0 ]; then
            echo "Collapsing command failed with status $status" >> {log.err_log}
            exit $status
        fi
        """
