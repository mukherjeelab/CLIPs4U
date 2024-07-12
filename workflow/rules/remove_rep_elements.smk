if config["rem_rep"] in ["True", True, 1, "1"]:
    group: lambda wildcards: wildcards.sample
    rule remove_rep:
        input:
            expand(f"{cwd}/reads/{{sample}}_trimmed_collapsed.fa", sample=sample_names.keys()),
            collapsed_read=lambda wildcards: f"{cwd}/reads/{wildcards.sample}_trimmed_collapsed.fa"
        output:
            rm_rep_read=f"{cwd}/reads/{{sample}}_rm_rep.fa",
            unaligned_bam=f"{cwd}/ali/{{sample}}_reps.bam"
        params:
            star_prefix=f"{cwd}/ali/{{sample}}_rm_rep",
            index_path=params["rep_idx"],
            params=params["star_rem_reps_params"],
            pp2_dir=parpipe2_dir,
            t=num_cores
        log:
            out_log = f"{cwd}/logs/rm_reps/{{sample}}.log",
            err_log = f"{cwd}/logs/rm_reps/{{sample}}.err"
        shell:
            """
            threads="{params.t}"
            prefix="{params.star_prefix}"
            bam="{params.star_prefix}Aligned.out.bam"
            fa="{params.star_prefix}Unmapped.out.mate1"
            rmrep_params="{params.params}"
            idx_path="{params.index_path}"

            STAR $rmrep_params --runThreadN $threads --genomeDir $idx_path --outFileNamePrefix $prefix --readFilesIn {input.collapsed_read} >> {log.out_log} 2>> {log.err_log}

            mv $fa {output.rm_rep_read}
            mv $bam {output.unaligned_bam}

            status=$?
            if [ $status -ne 0 ]; then
                echo "Removing repetitive elaments failed with status $status" >> {log.out_log}
                cat {log.err_log} >> {log.out_log}
                exit $status
            fi
            """
else:
    rule copy_collapsed_files:
        group: lambda wildcards: wildcards.sample
        input:
            expand(f"{cwd}/reads/{{sample}}_trimmed_collapsed.fa", sample=sample_names.keys()),
            collapsed_read=lambda wildcards: f"{cwd}/reads/{wildcards.sample}_trimmed_collapsed.fa"
        output:
            rm_rep_read=f"{cwd}/reads/{{sample}}_rm_rep.fa"
        shell:
            """
            cp {input.collapsed_read} {output.rm_rep_read}
            """
