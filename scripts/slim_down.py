#! /usr/bin/env python3


import os
import pandas as pd
from Bio import SeqIO
import subprocess


METADATA_PATH = '/mnt/scratch/brownlab/rej110/projects/phylofisher-current-protocols/small_database/metadata.tsv'
ORTHOLOG_DIRECTORY = '/mnt/scratch/brownlab/rej110/projects/phylofisher-current-protocols/small_database/orthologs/'
PARALOG_DIRECTORY = '/mnt/scratch/brownlab/rej110/projects/phylofisher-current-protocols/small_database/paralogs/'
SMALL_PROTEOMES_DIRECTORY = '/mnt/scratch/brownlab/rej110/projects/phylofisher-current-protocols/small_database/proteomes/'
BIG_PROTEOMES_DIRECTORY = '/mnt/scratch/brownlab/rej110/projects/phylofisher-current-protocols/PhyloFisherDatabase_v1.0/database/proteomes'
TAXA_TO_KEEP = [
    'Gregniph', 'Protadhe', 'Stensten', 'Idiovort', 'Vermverm', 'Ancysigm', 'Prascaps', 'Physpate', 'Diphrota',
    'Gemicryo', 'Trypbruc', 'Pharkirb', 'Cyanpara', 'ChoanoFB', 'Prymparv', 'Hemikukw', 'Gefiokel', 'Carpmemb',
    'Thectrah', 'Allomacr', 'Homosapi', 'PB58411a', 'Bigenata', 'Leptvora', 'Rhodlimn', 'Galdphle', 'Rhodmari',
    'Deveeleg', 'PoterBG1', 'Phytpara', 'TeloneP2'
    ]


def parse_metadata(metadata_path):
    meta_dict = {}
    with open(metadata_path, 'r') as metadata_file:
        metadata_file.readline()
        for line in metadata_file:
            metadata = line.strip().split('\t')
            meta_dict[metadata[0]] = {
                'name': metadata[1],
                'higher': metadata[2],
                'lower': metadata[3],
                'type': metadata[4],
                'notes': metadata[5],
                'orthologs': 0,
                'paralogs': 0,
            }

    return meta_dict

def get_homologs(meta_dict, directory, hom_type):

    for filename in os.listdir(directory):
        if filename.endswith('.fas'):
            with open(directory + filename, 'r') as file:
                for line in file:
                    for taxon in meta_dict:
                        if line.startswith('>') and taxon in line:
                            meta_dict[taxon][hom_type] += 1

    return meta_dict


def slim_down():
    for directory in [ORTHOLOG_DIRECTORY, PARALOG_DIRECTORY]:

        for filename in os.listdir(directory):
            if filename.endswith('.fas'):
                records = []
                with open(directory + filename, 'r') as file:
                    for record in SeqIO.parse(file, 'fasta'):
                        for taxon in TAXA_TO_KEEP:
                            if taxon in record.description:
                                records.append(record)
                
                SeqIO.write(records, directory + filename, 'fasta')

def cp_proteomes():
    for filename in os.listdir(BIG_PROTEOMES_DIRECTORY):
        if filename.split('.')[0] in TAXA_TO_KEEP:
            cmd = f'cp {BIG_PROTEOMES_DIRECTORY}/{filename} {SMALL_PROTEOMES_DIRECTORY}/{filename}'
            subprocess.run(cmd, shell=True)

if __name__ == '__main__':
    # Get ortholog an paralog counts
    # metadata_dict = parse_metadata(METADATA_PATH)
    # metadata_dict = get_homologs(metadata_dict, ORTHOLOG_DIRECTORY, 'orthologs')
    # metadata_dict = get_homologs(metadata_dict, PARALOG_DIRECTORY, 'paralogs')

    # Remove unwanted taxa from orthologs and paralogs
    # slim_down()

    # Copy over proteomes
    cp_proteomes()
    