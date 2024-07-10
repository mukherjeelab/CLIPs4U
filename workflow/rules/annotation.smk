rule annotate_groups_and_clusters:
    input:
        groups=lambda wildcards: f"{cwd}/stats/{wildcards.sample}_groups_conv_stats.tsv",
        clusters=lambda wildcards: f"{cwd}/stats/{wildcards.sample}_clusters_conv_stats.tsv"
    output:
        stats_groups_annot=f"{cwd}/stats/{{sample}}_groups_conv_stats_annotated.tsv",
        groups_annot=f"{cwd}/annot/{{sample}}_groups_annotated.tsv",
        stats_clust_annot=f"{cwd}/stats/{{sample}}_clusters_conv_stats_annotated.tsv",
        clust_annot=f"{cwd}/annot/{{sample}}_clusters_annotated.tsv",
        clust_filtered_bed=f"{cwd}/genome_viewer_files/{{sample}}_clusters_filtered.bed",
        clust_unfiltered_bed=f"{cwd}/genome_viewer_files/{{sample}}_clusters_unfiltered.bed"
    params:
        pp2_dir=parpipe2_dir,
        script=os.path.join(parpipe2_dir, "workflow", "scripts", "annot_groups_and_clusters.R"),
        annot_rank=params["annot_rank"],
        main_annot=main_annot_path(),
        trna_annot=os.path.join(parpipe2_dir, "annotation", params["organism"], "rds_files", "trna.rds"),
        repeats_annot=os.path.join(parpipe2_dir, "annotation", params["organism"], "rds_files", "repeats.rds"),
        organism=params["organism"]
    log:
        out_log = f"{cwd}/logs/annot/{{sample}}.log",
        err_log = f"{cwd}/logs/annot/{{sample}}.err"
    shell:
        """
        Rscript {params.script} $(pwd) {params.annot_rank} {params.main_annot} {params.trna_annot} {params.repeats_annot} {input.groups} {input.clusters} {params.organism} 2>> {log.err_log}
        """
