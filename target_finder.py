"""
Computer Science 225 Final Project - Target Finder 
- written by Harry Yang
"""
import sys
import os 
import numpy as np 
import scipy 
import pandas as pd 
from collections import defaultdict, Counter 
from operator import itemgetter

####################################
#                                  #
#            Information           #
#                                  #
####################################


"""

STATE NO.   MNEMONIC    DESCRIPTION COLOR NAME  COLOR CODE
1   TssA    Active TSS  Red 255,0,0
2   TssAFlnk    Flanking Active TSS Orange Red  255,69,0
3   TxFlnk  Transcr. at gene 5' and 3'  LimeGreen   50,205,50
4   Tx  Strong transcription    Green   0,128,0
5   TxWk    Weak transcription  DarkGreen   0,100,0
6   EnhG    Genic enhancers GreenYellow 194,225,5
7   Enh Enhancers   Yellow  255,255,0
8   ZNF/Rpts    ZNF genes & repeats Medium Aquamarine   102,205,170
9   Het Heterochromatin PaleTurquoise   138,145,208
10  TssBiv  Bivalent/Poised TSS IndianRed   205,92,92
11  BivFlnk Flanking Bivalent TSS/Enh   DarkSalmon  233,150,122
12  EnhBiv  Bivalent Enhancer   DarkKhaki   189,183,107
13  ReprPC  Repressed PolyComb  Silver  128,128,128
14  ReprPCWk    Weak Repressed PolyComb Gainsboro   192,192,192
15  Quies   Quiescent/Low   White   255,255,255

"""

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
    count = 0
    with open(input_file_list, 'r') as f:
        for line in f:
            # line = E023_15_coreMarks_mnemonics.bed
            count += 1
            sample_id = line.split('/')[-1].split('_')[0]
            state_dict = read_roadmap_bed(line.rstrip())
            sample_dict[sample_id] = state_dict
            if count % 10 == 0:
                print "Processed: %i tissues" % (count) 
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
    snp_hits = []
    with open(input_t1d_snp_hits, 'r') as f:
        for line in f:
            SNP_name = line.split()[0]
            chromosome = "chr" + line.split()[1]
            position = line.split()[2]
            snp_hits.append((chromosome, position))
    return snp_hits

######################
# Analysis Functions #
######################

def t1d_match(input_roadmap_location_dict, input_snp_hits_list):
    # INPUT - Roadmap_dict - key = tissue ID, value = dict with chromosome -> list of positions
    # INPUT - SNP Hits list - (chromosome, position) tuple

    # Make list SNP first sorted by chromosome
    snp_list = sorted(input_snp_hits_list, key = itemgetter(0))
    # print snp_list

    # look for state in the samples
    SNP_search_result_dict = {}
    for snp in snp_list:
        tissue_result_list = []
        chromosome = snp[0] 
        position = snp[1]
        state = ""
        for tissue, tissue_specific_marker_dict in input_roadmap_location_dict.iteritems():
            previous_position = -1
            previous_chromatin_state = ""
            for gene_block_tuple in tissue_specific_marker_dict[chromosome]:
                # gene_block_tuple = (start_pos, state)
                gene_position = gene_block_tuple[0]
                gene_state = gene_block_tuple[1]
                if gene_position >= position:
                    state = previous_chromatin_state
                    break
                else:
                    previous_position = gene_position
                    previous_chromatin_state = gene_state
            tissue_result_list.append((tissue, state))
        SNP_search_result_dict[snp] = tissue_result_list
    return SNP_search_result_dict


def find_good_match(metadata, snp_tissue_dict):
    # INPUT - metadata - tissue ID -> description 
    # INPUT - SNP tissue dict SNP(chr, pos) -> [(tissue, state)]
    active_state_list = [1,2,4,5,6,7]
    total_count = Counter()
    score_dict = defaultdict(float)
    for snp, tissue_result_list in snp_tissue_dict.iteritems():
        count = Counter()
        active_tissue_list = []
        for tissue_state_tuple in tissue_result_list:
            tissue = tissue_state_tuple[0]
            state = tissue_state_tuple[1]
            state_number = int(state.split('_')[0])
            if state_number in active_state_list:
                tissue_info = metadata[tissue]
                # print snp,tissue_info
                # total_count[tissue] += 1
                active_tissue_list.append(tissue)
            count[state] += 1

        for tissue in active_tissue_list:
            score_dict[tissue] += 1.0/len(active_tissue_list)

        # print count 
    # print total_count


    # threshold > 
    print score_dict
    sorted_score_list = sorted(score_dict.iteritems(), key = lambda (k,v): v, reverse = True)
    for tissue_score_tuple in sorted_score_list:
        tissue = tissue_score_tuple[0]
        score = tissue_score_tuple[1]

        if score >= 1.0:
            print tissue, metadata[tissue], score 


def test_pipeline_T1D():
    # read whole set 
    # TODO - make file 
    # input files 
    encode_bed_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_15core_marks_list.txt'
    T1D_SNP_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/T1D_snp_list.txt'
    metadata_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_metadata.txt'

    # Read tissue files

    epigenetic_marks = read_whole_sets(encode_bed_list)
    metadata_list = read_metadata(metadata_file)

    snp_list = read_gwas_SNP_hits_T1D(T1D_SNP_list)

    # print glance(epigenetic_marks)
    # print snp_list
    print metadata_list

    # process
    #1. find block

    snp_result = t1d_match(epigenetic_marks, snp_list)
    find_good_match(metadata_list, snp_result)

if __name__ == "__main__":
    test_pipeline_T1D()