"""
Computer Science 225 Final Project - Target Finder 
- written by Harry Yang
"""
import sys
import os 
import numpy as np 
import scipy 
import pandas as pd 
from collections import defaultdict



####################################
#                                  #
#             Functions            #
#                                  #
####################################


#################
# I/O Functions #
#################

def read_roadmap_bed(input_file):
    # TSV
	# chromosome, start (0-based), stop (1-based), state_label_mnemonic for that region
    index_dict = defaultdict(list)
    with open(input_file, 'r') as f:
        for line in f:
            # access by chromosome, search by start 
            # if ordered, the position we want to find must be within the flanking start positions
            chromosome = line.split()[0]
            start      = line.split()[1]
            end        = line.split()[2]
            state      = line.split()[3]
            index_dict[chromosome].append((start, state))
    return index_dict

def read_whole_sets(input_file_list):
    sample_dict = {}
    with open(input_file_list, 'r') as f:
        for line in f:
            # line = E023_15_coreMarks_mnemonics.bed
            sample_id = line.split('/')[-1].split('_')[0]
            state_dict = read_roadmap_bed(line)
            sample_dict[sample_id] = state_dict

    return sample_dict

def read_metadata(input_metadata):
    metadata_dict = {}
    with open(input_metadata, 'r' ) as f:
        for line in f:
            # sample_id (e***) sample ID
            sample_id = line.split()[0]
            sample_description = "_".join(line.split()[1:])
            metadata_dict[sample_id] = sample_description
    return metadata_dict

def read_gwas_SNP_hits_T1D(input_t1d_snp_hits):
    snp_hits = {}
    with open(input_tid1_snp_hits, 'r') as f:
        for line in f:
            SNP_name = line.split()[0]
            chromosome = line.split()[1]
            position = line.split()[2]
            snp_hits[SNP_name] = (chromosome, position)
    return snp_hits



def test_pipeline_T1D():
    # read whole set 
    # TODO - make file 
    # input files 
    encode_bed_list = ''
    T1D_SNP_list = ''
    metadata_file = ''

    # Read tissue files

    epigenetic_marks = read_whole_sets(input_file)
    metadata_list = read_metadata(metadata_file)

    snp_list = read_gwas_SNP_hits_T1D(T1D_SNP_list)


    # process
    #1. find block