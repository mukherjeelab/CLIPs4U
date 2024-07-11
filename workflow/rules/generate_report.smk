rule generate_report:
    input:
        expand(f"{cwd}/motif_enrichment/{{sample}}_report.txt", sample=sample_names.keys())
    output:
        "final_report.html"
    params:
        script=os.path.join(parpipe2_dir, "workflow", "scripts", "generate_report.R"),
        main_annot=main_annot_path(),
        conv=params["CONVERSION"],
        method=params["motif_enrichment_method"],
        rmd_file=os.path.join(parpipe2_dir, "workflow", "scripts", "report.Rmd")
    shell:
        """
        Rscript {params.script} $(pwd) {params.main_annot} "{params.conv}" {params.method} {params.rmd_file}
        rm gene_types.html
        rm datatable_summary.html
        rm summary_stats.html
        """

