
rule get_genomic_files_and_index:
    input:
        script=os.path.join(parpipe2_dir, "workflow", "scripts", "genomic_files_and_indexes.py")
    output:
        "metadata.json"
    shell:
        prepare_genomes_and_indexes(params, parpipe2_dir, workflow)
