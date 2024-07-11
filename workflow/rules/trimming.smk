rule trim_adapters_r1:
    group: lambda wildcards: wildcards.sample
    input:
        expand(f"{cwd}/qc/fastqc/{{sample}}_fastqc.html", sample=sample_names.keys()),
        expand(f"{cwd}/qc/fastqc/{{sample}}_fastqc.zip", sample=sample_names.keys()),
        file=lambda wildcards: sample_names[wildcards.sample]
    output:
        trim_reads1 = f"{cwd}/reads/{{sample}}_trimmed1.fq.gz"
    params:
        three_prime_adapter=params['three_prime_adapter'],
        five_prime_adapter=params['five_prime_adapter'],
        cutadapt_params=params['cutadapt_params1'],
        nextseq=params['nextseq']
    log:
        out_log = f"{cwd}/logs/cutadapt_r1/{{sample}}.log",
        err_log = f"{cwd}/logs/cutadapt_r1/{{sample}}.err"
    shell:
        """
        mkdir -p reads
        
        nextseq_trim=""
        if [ "{params.nextseq}" = "True" ]; then
            nextseq_trim="--nextseq-trim=15"
        fi

        adapters=""
        if [ "{params.three_prime_adapter}" != "" ] && [ "{params.five_prime_adapter}" != "" ]; then
            adapters="-a {params.three_prime_adapter} -g {params.five_prime_adapter}"
        elif [ "{params.three_prime_adapter}" != "" ]; then
            adapters="-a {params.three_prime_adapter}"
        elif [ "{params.five_prime_adapter}" != "" ]; then
            adapters="-g {params.five_prime_adapter}"
        fi

        cutadapt_params="{params.cutadapt_params}"

        echo "Command: cutadapt $adapters $nextseq_trim $cutadapt_params -o {output.trim_reads1} {input.file}" >> {log.out_log}
        cutadapt $adapters $nextseq_trim $cutadapt_params -o {output.trim_reads1} {input.file} >> {log.out_log} 2>> {log.err_log}
        status=$?
        if [ $status -ne 0 ]; then
            echo "cutadapt failed with status $status" >> {log.out_log}
            cat {log.err_log} >> {log.out_log}
            exit $status
        fi
        """
        
if config["cutadapt_params2"] in ["", None]:
    rule move_trimmed_files:
        input:
            expand(f"{cwd}/reads/{{sample}}_trimmed1.fq.gz", sample=sample_names.keys()),
            trim_read=lambda wildcards: f"{cwd}/reads/{wildcards.sample}_trimmed1.fq.gz"
        output:
            trim_final=f"{cwd}/reads/{{sample}}_trimmed_final.fq.gz"
        shell:
            """
            echo "Skipping cutadapt round 2"
            mv {input.trim_read} {output.trim_final}
            """
else:
    rule trim_adapters_r2:
        input:
            expand(f"{cwd}/reads/{{sample}}_trimmed1.fq.gz", sample=sample_names.keys()),
            trim_read=lambda wildcards: f"reads/{wildcards.sample}_trimmed1.fq.gz"
        output:
            trim_final="{cwd}/reads/{{sample}}_trimmed_final.fq.gz"
        params:
            cutadapt_params=params['cutadapt_params2']
        log:
           out_log = f"{cwd}/logs/cutadapt_r2/{{sample}}.log",
            err_log = f"{cwd}/logs/cutadapt_r2/{{sample}}.err"
        shell:
            """
            cutadapt_params="{params.cutadapt_params}"
            echo "Command: cutadapt $cutadapt_params -o {output.trim_final} {input.trim_read}" >> {log.out_log}
            cutadapt $cutadapt_params -o {output.trim_final} {input.trim_read} >> {log.out_log} 2>> {log.err_log}
            status=$?
            if [ $status -ne 0 ]; then
                echo "cutadapt failed with status $status" >> {log.out_log}
                cat {log.err_log} >> {log.out_log}
                exit $status
            fi
            """
