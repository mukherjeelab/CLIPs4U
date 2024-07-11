if params["aligner"] in ['S', 's', 'star', 'Star', 'STAR']:
     rule genome_alignment_star:
         group: lambda wildcards: wildcards.sample
         input:
             expand(f"{cwd}/reads/{{sample}}_rm_rep.fa", sample=sample_names.keys()),
             rm_rep_read=lambda wildcards: f"{cwd}/reads/{wildcards.sample}_rm_rep.fa"
         output:
             aligned_sam= f"{cwd}/ali/{{sample}}.aligned.mapped.sorted.sam"
             #add all other star output files here (not to rule all)
         params:
             star_prefix=f"{cwd}/ali/{{sample}}",
             index_path=get_index_path('star'),
             params=params["star_map_params"],
             t=num_cores
         log:
             out_log = f"{cwd}logs/star_ali/{{sample}}.log",
             err_log = f"{cwd}logs/star_ali/{{sample}}.err"
         shell:
             """
             threads="{params.t}"
             prefix="{params.star_prefix}"
             bam="{params.star_prefix}Aligned.out.bam"
             sam="{params.star_prefix}.aligned.sam"
             star_params="{params.params}"

             STAR $star_params --runThreadN $threads --genomeDir {params.index_path} --outFileNamePrefix $prefix --readFilesIn {input.rm_rep_read} >> {log.out_log} 2>> {log.err_log}

             samtools view -F 4 -h -o $sam $bam >> {log.out_log} 2>> {log.err_log}
             samtools sort -o {output.aligned_sam} $sam >> {log.out_log} 2>> {log.err_log}

             status=$?
             if [ $status -ne 0 ]; then
                 echo "Mapping reads with STAR failed with status $status" >> {log.out_log}
                 cat {log.err_log} >> {log.out_log}
                 exit $status
             fi
             """
else:
     rule genome_alignment_bowtie:
         group: lambda wildcards: wildcards.sample
         input:
             expand(f"{cwd}/reads/{{sample}}_rm_rep.fa", sample=sample_names.keys()),
             rm_rep_read=lambda wildcards: f"{cwd}/reads/{wildcards.sample}_rm_rep.fa"
         output:
             aligned_sam= f"{cwd}/ali/{{sample}}.aligned.mapped.sorted.sam"
         params:
             index_path=get_index_path('bowtie'),
             params=params["bowtie_map_params"],
             t=num_cores
         log:
             out_log = f"{cwd}/logs/bowtie_ali/{{sample}}.log",
             err_log = f"{cwd}/logs/bowtie_ali/{{sample}}.err"
         shell:
             """
             threads="{params.t}"
             bowtie_params="{params.params}"

             bowtie {params.index_path}/bwt $bowtie_params -p $threads -f {input.rm_rep_read} 2>> {log.err_log} | samtools view -hS -F 4 - | samtools sort -O SAM - -o {output.aligned_sam} >> {log.out_log} 2>> {log.err_log}

             status=$?
             if [ $status -ne 0 ]; then
                 echo "Mapping reads with bowtie failed with status $status" >> {log.out_log}
                 cat {log.err_log} >> {log.out_log}
                 exit $status
             fi
             """
