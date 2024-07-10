rule prepare_paralyzer_ini:
    input:
        sam=lambda wildcards: f"{cwd}/ali/{wildcards.sample}.aligned.mapped.sorted.sam"
    output:
        paralyzer_ini=f"{cwd}/PARalyzer/{{sample}}.ini"
    params:
        genome_2bit=get_2bit_path(),
        bandwidth=params["BANDWIDTH"],
        conv=params["CONVERSION"],
        min_read_cnt_group=params["MINIMUM_READ_COUNT_PER_GROUP"],
        min_read_cnt_clust=params["MINIMUM_READ_COUNT_PER_CLUSTER"],
        min_read_cnt_kde=params["MINIMUM_READ_COUNT_FOR_KDE"],
        min_clust_size=params["MINIMUM_CLUSTER_SIZE"],
        min_conv_loc_clust=params["MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER"],
        min_conv_cnt_clust=params["MINIMUM_CONVERSION_COUNT_FOR_CLUSTER"],
        min_read_cnt_clust_incl=params["MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION"],
        min_read_l=params["MINIMUM_READ_LENGTH"],
        max_non_conv_mism=params["MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES"],
        signal_ext=params["signal_extension"],
        add_nt_signal=params["ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL"],
        pcr_dups=params["use_pcr_duplicates"]
    shell:
        """
        in_sam="{input.sam}"
        if [ "{params.pcr_dups}" = "True" ]; then
            in_sam="{input.sam}=COLLAPSED"
        fi
        bandwidth="{params.bandwidth}"
        conv="{params.conv}"
        min_read_cnt_group="{params.min_read_cnt_group}"
        min_read_cnt_clust="{params.min_read_cnt_clust}"
        min_read_cnt_kde="{params.min_read_cnt_kde}"
        min_clust_size="{params.min_clust_size}"
        min_conv_loc_clust="{params.min_conv_loc_clust}"
        min_conv_cnt_clust="{params.min_conv_cnt_clust}"
        min_read_cnt_clust_incl="{params.min_read_cnt_clust_incl}"
        min_read_l="{params.min_read_l}"
        max_non_conv_mism="{params.max_non_conv_mism}"
        signal_ext="{params.signal_ext}"
        add_nt_signal=""
        if [ "{params.signal_ext}" = "ADDITIONAL_NUCLEOTIDES_BEYOND_SIGNAL" ]; then
            signal_ext="={params.add_nt_signal}"
            if [ ! -z "{params.add_nt_signal}" ] && [ "{params.add_nt_signal}" != "['']" ]; then
                signal_ext="={params.add_nt_signal}"
            fi
        fi
        cat <<EOL > {output.paralyzer_ini}

BANDWIDTH=$bandwidth
CONVERSION=$conv
MINIMUM_READ_COUNT_PER_GROUP=$min_read_cnt_group
MINIMUM_READ_COUNT_PER_CLUSTER=$min_read_cnt_clust
MINIMUM_READ_COUNT_FOR_KDE=$min_read_cnt_kde
MINIMUM_CLUSTER_SIZE=$min_clust_size
MINIMUM_CONVERSION_LOCATIONS_FOR_CLUSTER=$min_conv_loc_clust
MINIMUM_CONVERSION_COUNT_FOR_CLUSTER=$min_conv_cnt_clust
MINIMUM_READ_COUNT_FOR_CLUSTER_INCLUSION=$min_read_cnt_clust_incl
MINIMUM_READ_LENGTH=$min_read_l
MAXIMUM_NUMBER_OF_NON_CONVERSION_MISMATCHES=$max_non_conv_mism

$signal_ext$add_nt_signal

GENOME_2BIT_FILE={params.genome_2bit}
SAM_FILE=$in_sam

OUTPUT_DISTRIBUTIONS_FILE=$(pwd)/PARalyzer/{wildcards.sample}.distribution
OUTPUT_GROUPS_FILE=$(pwd)/PARalyzer/{wildcards.sample}.groups
OUTPUT_CLUSTERS_FILE=$(pwd)/PARalyzer/{wildcards.sample}.clusters
OUTPUT_READS_FILE=$(pwd)/PARalyzer/{wildcards.sample}_PARalyzer_Utilized.sam

EOL
        """
        
rule run_paralyzer:
    input:
        ini_file=lambda wildcards: f"{cwd}/PARalyzer/{wildcards.sample}.ini"
    output:
        bam=f"{cwd}/PARalyzer/{{sample}}_PARalyzer_Utilized.bam",
        bam_idx=f"{cwd}/PARalyzer/{{sample}}_PARalyzer_Utilized.bam.bai",
        distr=f"{cwd}/PARalyzer/{{sample}}.distribution",
        groups=f"{cwd}/PARalyzer/{{sample}}.groups",
        clusters=f"{cwd}/PARalyzer/{{sample}}.clusters"
    params:
        mem=params["paralyzer_memory"]
    log:
        out_log = f"{cwd}/logs/PARalyzer/{{sample}}.log",
        err_log = f"{cwd}/logs/PARalyzer/{{sample}}.err"
    shell:
        """
        mem="{params.mem}"

        PARalyzer $mem {input.ini_file} >> {log.out_log} 2>> {log.err_log}
        samtools sort -o {output.bam} PARalyzer/{wildcards.sample}_PARalyzer_Utilized.sam  >> {log.out_log} 2>> {log.err_log}
        samtools index {output.bam} >> {log.out_log} 2>> {log.err_log}

        status=$?
            if [ $status -ne 0 ]; then
                echo "Running PARalyzer failed with status $status" >> {log.out_log}
                cat {log.err_log} >> {log.out_log}
                exit $status
            fi
        """
        
rule calc_conv_stats:
    input:
        bam=lambda wildcards: f"{cwd}/PARalyzer/{wildcards.sample}_PARalyzer_Utilized.bam",
        groups=lambda wildcards: f"{cwd}/PARalyzer/{wildcards.sample}.groups",
        clusters=lambda wildcards: f"{cwd}/PARalyzer/{wildcards.sample}.clusters",
    output:
        groups_conv=f"{cwd}/stats/{{sample}}_groups_conv_stats.tsv",
        clusters_conv=f"{cwd}/stats/{{sample}}_clusters_conv_stats.tsv"
    params:
        script=os.path.join(parpipe2_dir, "workflow", "scripts", "add_conv_stats.py"),
        genome=get_genome_path(),
        t=num_cores,
        conv=params["CONVERSION"]
    log:
        out_log = f"{cwd}/logs/conv_stats/{{sample}}.log",
        err_log = f"{cwd}/logs/conv_stats/{{sample}}.err"
    shell:
        """
        python {params.script} --bamfile {input.bam} --groups {input.groups} --clusters {input.clusters} --genome {params.genome} --threads {params.t} --conv "{params.conv}" 2>> {log.err_log}
        """
