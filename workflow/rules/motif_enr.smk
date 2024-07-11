rule motif_enrichment:
    group: lambda wildcards: wildcards.sample
    input:
        expand(f"{cwd}/genome_viewer_files/{{sample}}_report.txt", sample=sample_names.keys()),
        filt_bed=lambda wildcards: f"{cwd}/genome_viewer_files/{wildcards.sample}_clusters_filtered.bed"
    output:
        motif_enr_report=f"{cwd}/motif_enrichment/{{sample}}_report.txt"
    params:
        genome=get_genome_path(),
        method=params["motif_enrichment_method"],
        params=params["motif_enrichment_params"]
    log:
        out_log = f"{cwd}/logs/motif_enrichment/{{sample}}.log",
        err_log = f"{cwd}/logs/motif_enrichment/{{sample}}.err"
    shell:
        """
        echo "Running motif_enrichment"

        line_count_filt=$(wc -l < {input.filt_bed} | tr -d ' ')

        if [ $line_count_filt -gt 0 ]; then
            method="{params.method}"
            params="{params.params}"
            bedtools getfasta -fi {params.genome} -bed {input.filt_bed} -s -fo motif_enrichment/{wildcards.sample}_clusters_filtered.fa >> {log.out_log} 2>> {log.err_log}
            $method motif_enrichment/{wildcards.sample}_clusters_filtered.fa -oc motif_enrichment/{wildcards.sample}_motif_enr $params >> {log.out_log} 2>> {log.err_log}
            cat {log.out_log} >> {output.motif_enr_report}
            cat {log.err_log} >> {output.motif_enr_report}
        else
            echo "{input.filt_bed} is empty" >> {log.out_log}
            cat {log.out_log} >> {output.motif_enr_report}
        fi
        """
