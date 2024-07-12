import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
import time
import wget
import yaml


def unpack_gzip_file(filename: str, output_filename: str) -> None:
    """Unpacks a gzipped file and save it to specified output"""
    with gzip.open(filename, 'rb') as f_in:
        with open(output_filename, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


def check_checksum(filename: str, checksum_filename: str, checksum_algorithm:str='md5') -> bool:
    """
    Check checksum against the specified file
    Available checksum algorithms: sha1, sha224, sha256, sha384, sha512,
    sha3_224, sha3_256, sha3_384, sha3_512, shake_128, shake_256, blake2b, md5
    """
    if not hasattr(hashlib, checksum_algorithm):
        raise ValueError('No such algorithm')
    if not os.path.isfile(filename):
        raise FileNotFoundError(filename)
    if not os.path.isfile(checksum_filename):
        raise FileNotFoundError(checksum_filename)
    base_filename = os.path.basename(filename)
    with open(checksum_filename, 'r') as f:
        for line in f:
            elements = line.split()
            if elements[1] == base_filename:
                checksum = elements[0]
    algorithm = getattr(hashlib, checksum_algorithm)
    return algorithm(open(filename,'rb').read()).hexdigest() == checksum


def download_and_validate_checksum(url: str, outdir: str, checksum_url: str, max_repeats:int=3) -> bool:
    """
    Try do download and validate file with specified maximum number of repeats
    Args:
        url - url to download
        outdir - where to save files
        checksum_url - url with checksums
        max_repeats - maximum number of repets for downloading file
    Returns:
        True if download was successful, else False
    """
    max_tries = max_repeats
    while max_tries:
        wget.download(url, out=outdir)
        filename = url.split("/")[-1]
        wget.download(checksum_url, out=outdir)
        checksum_filename = checksum_url.split("/")[-1]
        filename_fullpath = os.path.join(outdir, filename)
        checksum_fullpath = os.path.join(outdir, checksum_filename)
        is_ok = check_checksum(filename_fullpath, checksum_fullpath, checksum_algorithm='md5')
        if os.path.exists(checksum_fullpath): # need to delete checksum to prevent problems with other downloads
            os.remove(checksum_fullpath)
        if is_ok:
            print('Download successful, file validated')
            return True
        # cleanup in case of error
        if os.path.exists(filename_fullpath):
            os.remove(filename_fullpath)
    return False


def update_metadata_file(filename="metadata.json", **kwargs):
    with open(filename, "r") as metafile:
        metadata = json.load(metafile)
    metadata.update(kwargs)
    with open(filename, "w") as metafile:
        json.dump(metadata, metafile, indent=4)


def prepare_genome_and_annotation(configuration, parpipe2_dir, cwd, workflow_dir):
    outdir = os.path.join(parpipe2_dir, 'genome')
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    logdir = os.path.join(cwd, 'logs', 'genome_files')
    if not os.path.exists(logdir):
        os.makedirs(logdir, exist_ok=True)
  
    if configuration["genome_fasta"] and configuration["genome_fasta"] != "":
        genome = configuration["genome_fasta"]
        annotation = configuration["gtf"]
        print("Genome and annotation files are provided by user")
    else:
        genome_file = "GRCh38.primary_assembly.genome.fa"
        genome_ver = configuration.get('genome_version', '45')
        if not genome_ver:
            genome_ver = '45'
        if int(genome_ver) > 41:
            gtf_file = f"gencode.v{genome_ver}.primary_assembly.basic.annotation.gtf"
        else:
            gtf_file = f"gencode.v{genome_ver}.primary_assembly.annotation.gtf"
        
        print(f"Selected genome version is {genome_ver}")
        genome = os.path.join(outdir, genome_file)
        annotation = os.path.join(outdir, gtf_file)

        url_genome = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{genome_ver}/{genome_file}.gz"
        url_gtf = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{genome_ver}/{gtf_file}.gz"
        url_checksum = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{genome_ver}/MD5SUMS"

        # checking and downloading separately for genome and annotation
        for full_filepath, filename, url, genome_or_annotation in [
                (genome, genome_file, url_genome, "genome"), (annotation, gtf_file, url_gtf, "annotation")]:
            if not os.path.isfile(full_filepath):
                file_gz = os.path.join(outdir, f"{filename}.gz")
                print(f"Downloading GRCh38 v.{genome_ver} {genome_or_annotation} from Gencode...\n")
                if not download_and_validate_checksum(url, outdir, url_checksum):
                    error = f"Cannot download and validate file: {url}"
                    print(error)
                    raise Exception(error)
                print(f"Unpacking {genome_or_annotation}...")
                try:
                    unpack_gzip_file(file_gz, file_gz.replace('.gz', ''))
                except Exception as e:
                    print(f"An ERROR occurred during {genome_or_annotation} unpacking: {e}")
                    raise
                else:
                    print(f"{genome_or_annotation.capitalize()} unpacked successfully")
            else:
                print(f'{genome_or_annotation.capitalize()} files already exist')
    main_annot_rds = os.path.join(outdir, 'main_annotation.rds')
    
    if os.path.isfile(main_annot_rds) and os.path.getmtime(main_annot_rds) > os.path.getmtime(annotation):
        print('main_annotation.rds file from current gtf already exists')
    else:
        print("Preparing main_annotation.rds file")
        logfile = os.path.join(logdir, 'prep_main_rds_from_gtf.log')
        errfile = os.path.join(logdir, 'prep_main_rds_from_gtf.err')
        ann_script = os.path.join(workflow_dir, "scripts", "prepareAnnotFilesFromGtf.R")
        ann_cmd = f"Rscript {ann_script} {outdir} {annotation} main > {logfile} 2> {errfile}"
        outcode = os.WEXITSTATUS(os.system(ann_cmd))
    
        if outcode == 0:     
            print(f"main_annotation.rds file created successfully\ncheck {logfile} for details")
        else:
            print(f"An ERROR occurred during creating main_annotation.rds file :( check {errfile} for details")
    
    genome_2bit = os.path.normpath(str(genome).replace('fasta', '2bit').replace('fa', '2bit'))
    
    if os.path.isfile(genome_2bit) and os.path.isfile(f'{genome}.fai'):
        print('Index and 2bit files exist. Moving to the next step')
    else:
        logfile = os.path.join(logdir, f'genome_indexing_and_2bit_preparation.log')
        errfile = os.path.join(logdir, f'genome_indexing_and_2bit_preparation.err')
        index_genome_cmd = f"samtools faidx {genome} > {logfile} 2> {errfile}"
        command_2bit = f"faToTwoBit {genome} {genome_2bit} >> {logfile} 2>> {errfile}"
        outcode = os.WEXITSTATUS(os.system(index_genome_cmd))
        if outcode == 0:     
            print(f"Genome file indexed successfully\ncheck {logfile} for details")
        else:
            print(f"An ERROR occurred during genome indexing :( check {errfile} for details")
        print("Preparing genome 2bit file...")
        outcode = os.WEXITSTATUS(os.system(command_2bit))
        if outcode == 0:     
            print(f"Genome 2bit file prepared successfully\ncheck {logfile} for details")
        else:
            print(f"An ERROR occurred during genome 2bit preparation :( check {errfile} for details")

    # Write the paths to a metadata file
    metadata = {
        "genome": genome,
        "annotation": annotation,
        "genome_2bit": genome_2bit,
        "main_annot_rds": main_annot_rds
    }
    with open(os.path.join(cwd, "metadata.json"), "w") as metafile:
        json.dump(metadata, metafile, indent=4)
    return genome, annotation, genome_2bit, main_annot_rds


def prepare_star_index(configuration, genome, annotation, threads, parpipe2_dir, cwd):
    logdir = os.path.join(cwd, 'logs', 'star_idx')
    if not os.path.exists(logdir):
        os.makedirs(logdir, exist_ok=True)

    if not configuration.get('star_index_dir'):
        genome_version = configuration.get('genome_version', '45')
        if not genome_version:
            genome_version = '45'
        star_idx_outdir = os.path.join(parpipe2_dir, 'star_index', os.path.basename(genome).replace('.fa', ''), genome_version)
        if not os.path.exists(star_idx_outdir):
            os.makedirs(star_idx_outdir, exist_ok=True)
        if os.path.exists(os.path.join(star_idx_outdir, 'chrName.txt')):
            print('STAR index found, moving to the next step')
        else:
            logfile = os.path.join(logdir, 'star_index.log')
            errfile = os.path.join(logdir, 'star_index.err')

            default_params = {
                "--runMode": "genomeGenerate",
                "--runThreadN": threads,
                "--genomeDir": star_idx_outdir,
                "--genomeFastaFiles": genome,
                "--sjdbGTFfile": annotation,
                "--sjdbOverhang": "100"
            }
            if configuration.get('star_index_params'):
                params = configuration['star_index_params'].replace(",", " ").split(" ")
                params_dict = {params[2*i]: params[2*i+1] for i in range(len(params)//2)}
                default_params.update(params_dict)
            params_str = " ".join([f"{key} {value}" for key, value in default_params.items()])
            command = f"STAR {params_str} > {logfile} 2> {errfile}"
            print("Generating STAR index...")
            outcode = os.WEXITSTATUS(os.system(command))
            if outcode == 0:
                print(f"STAR index generated successfully\ncheck {logfile} for details")
            else:
                print(f"An ERROR occurred during STAR index generation :(\ncheck {errfile} for details")
    else:
        print("STAR index provided by user")
        star_idx_outdir = configuration['star_index_dir']

    update_metadata_file(star_index=star_idx_outdir)
    return star_idx_outdir

def prepare_bowtie_index(configuration, genome, threads, parpipe2_dir, cwd):
    logdir = os.path.join(cwd, 'logs', 'bowtie_idx')
    if not os.path.exists(logdir):
        os.makedirs(logdir, exist_ok=True)

    if not configuration.get('bowtie_index_dir'):
        genome_version = configuration.get('genome_version', '45')
        if not genome_version:
            genome_version = '45'
        bowtie_idx_outdir = os.path.join(parpipe2_dir, 'bowtie_index', os.path.basename(genome).replace('.fa', ''), genome_version)
        if not os.path.exists(bowtie_idx_outdir):
            print(bowtie_idx_outdir)
            os.makedirs(bowtie_idx_outdir, exist_ok=True)

        if os.path.exists(os.path.join(bowtie_idx_outdir, 'bwt.1.ebwt')):
            print('Bowtie index found, moving to the next step')
        else:
            logfile = os.path.join(logdir, 'bowtie_index.log')
            errfile = os.path.join(logdir, 'bowtie_index.err')
            if not configuration.get('bowtie_index_params'):
                command = f"bowtie-build --threads {threads} -f {genome} {bowtie_idx_outdir}/bwt > {logfile} 2> {errfile}"
            else:
                params = " ".join(configuration['bowtie_index_params'].split(","))
                command = f"bowtie-build {params} -f {genome} {bowtie_idx_outdir}/bwt > {logfile} 2> {errfile}"

            print("Generating bowtie index...")
            outcode = os.WEXITSTATUS(os.system(command))
            if outcode == 0:
                print(f"Bowtie index generated successfully\ncheck {logfile} for details")
            else:
                print(f"An ERROR occurred during bowtie index generation :(\ncheck {errfile} for details")
    else:
        print("Bowtie index provided by user")
        bowtie_idx_outdir = configuration['bowtie_index_dir']

    update_metadata_file(bowtie_index=bowtie_idx_outdir)
    return bowtie_idx_outdir


def prepare_genomes_and_indexes(params: dict, parpipe2_dir: str, workflow) -> None:
    """
    Wrapper function that prepares the genomes and indexes for analysis
    Args:
        params: configuration
        parpipe2_dir: directory containging the parpipe2
        workflow: Workflow object, used to extract current directory and workflow_dir
    Returns: None
    """
    cwd = workflow.workdir_init
    workflow_dir = workflow.basedir
    threads = params['threads']
    aligner = params['aligner']
    genome, annotation, genome_2bit, main_annot_rds  = prepare_genome_and_annotation(params, parpipe2_dir, cwd, workflow_dir)
    if aligner.lower() in ['s', 'star']:
        prepare_star_index(params, genome, annotation, threads, parpipe2_dir, cwd)
    else :
        prepare_bowtie_index(params, genome, threads, parpipe2_dir, cwd)
