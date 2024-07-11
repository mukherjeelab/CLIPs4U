rule create_bigwigs:
    group: lambda wildcards: wildcards.sample
    input:
        expand(f"{cwd}/genome_viewer_files/{{sample}}_clusters_unfiltered.bed", sample=sample_names.keys()),
        bam=lambda wildcards: f"{cwd}/PARalyzer/{wildcards.sample}_PARalyzer_Utilized.bam",
        filt_bed=lambda wildcards: f"{cwd}/genome_viewer_files/{wildcards.sample}_clusters_filtered.bed",
        unfilt_bed=lambda wildcards: f"{cwd}/genome_viewer_files/{wildcards.sample}_clusters_unfiltered.bed"
    output:
        bw_report=f"{cwd}/genome_viewer_files/{{sample}}_report.txt"
    params:
        t=num_cores,
        eff_size=eff_size
    log:
        out_log = f"{cwd}/logs/bed_to_bigwig/{{sample}}.log",
        err_log = f"{cwd}/logs/bed_to_bigwig/{{sample}}.err"
    shell:
        """
        echo "Running create_bigwigs"

        line_count_filt=$(wc -l < {input.filt_bed} | tr -d ' ')
        line_count_unfilt=$(wc -l < {input.unfilt_bed} | tr -d ' ')

        if [ $line_count_filt -gt 0 ]; then
            bedtools intersect -abam {input.bam} -b {input.filt_bed} -s > genome_viewer_files/{wildcards.sample}_clusters_filtered.bam 2>> {log.err_log}
            samtools index genome_viewer_files/{wildcards.sample}_clusters_filtered.bam
            bamCoverage -b genome_viewer_files/{wildcards.sample}_clusters_filtered.bam -o genome_viewer_files/{wildcards.sample}_clusters_filtered_fwd.bw --binSize 1 --effectiveGenomeSize {params.eff_size} -p {params.t} --normalizeUsing CPM --filterRNAstrand forward >> {log.out_log} 2>> {log.err_log}
            bamCoverage -b genome_viewer_files/{wildcards.sample}_clusters_filtered.bam -o genome_viewer_files/{wildcards.sample}_clusters_filtered_rev.bw --binSize 1 --effectiveGenomeSize {params.eff_size} -p {params.t} --normalizeUsing CPM --filterRNAstrand reverse >> {log.out_log} 2>> {log.err_log}
        else
            echo "{input.filt_bed} is empty" >> {log.out_log}
        fi

        if [ $line_count_unfilt -gt 0 ]; then
            bedtools intersect -abam {input.bam} -b {input.unfilt_bed} -s > genome_viewer_files/{wildcards.sample}_clusters_unfiltered.bam 2>> {log.err_log}
            samtools index genome_viewer_files/{wildcards.sample}_clusters_unfiltered.bam
            bamCoverage -b genome_viewer_files/{wildcards.sample}_clusters_unfiltered.bam -o genome_viewer_files/{wildcards.sample}_clusters_unfiltered_fwd.bw --binSize 1 --effectiveGenomeSize {params.eff_size} -p {params.t} --normalizeUsing CPM --filterRNAstrand forward >> {log.out_log} 2>> {log.err_log}
            bamCoverage -b genome_viewer_files/{wildcards.sample}_clusters_unfiltered.bam -o genome_viewer_files/{wildcards.sample}_clusters_unfiltered_rev.bw --binSize 1 --effectiveGenomeSize {params.eff_size} -p {params.t} --normalizeUsing CPM --filterRNAstrand reverse >> {log.out_log} 2>> {log.err_log}
            cat {log.out_log} >> {output.bw_report}
            cat {log.err_log} >> {output.bw_report}
        else
            echo "{input.unfilt_bed} is empty" >> {log.out_log}
            cat {log.out_log} >> {output.bw_report}
        fi
        """
