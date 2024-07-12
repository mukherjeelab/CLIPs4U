rule fastqc:
    input:
        input_filenames
    output:
        html=f"{cwd}/qc/fastqc/{{sample}}_fastqc.html",
        zip=f"{cwd}/qc/fastqc/{{sample}}_fastqc.zip"
    log:
        f"{cwd}/logs/fastqc/{{sample}}.log"
    params:
        t=num_cores
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p qc/fastqc
        echo "Running FastQC on {input}" > {log}
        fastqc --quiet -t {params.t} {input} --outdir qc/fastqc >> {log} 2>&1
        """
