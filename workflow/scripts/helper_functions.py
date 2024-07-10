import json


def get_index_path(aligner):
    with open("metadata.json") as f:
        data = json.load(f)
    return data[f"{aligner}_index"]


def get_2bit_path():
    with open("metadata.json") as f:
        data = json.load(f)
    return data["genome_2bit"]


def get_genome_path():
    with open("metadata.json") as f:
        data = json.load(f)
    return data["genome"]


def main_annot_path():
    with open("metadata.json") as f:
        data = json.load(f)
    return data["main_annot_rds"]