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
    state_dict = defaultdict(int)
    with open(input_file, 'r') as f:
        for line in f:
            # access by chromosome, search by start 
            # if ordered, the position we want to find must be within the flanking start positions
            chromosome = line.split()[0]
            start      = line.split()[1]
            end        = line.split()[2]
            state      = line.split()[3]
            index_dict[chromosome].append((start, state))
            state_dict[state] += (int(start) - int(end))
    return index_dict, state_dict

def read_whole_sets(input_file_list):
    sample_dict = {}
    sample_dict_state = {}
    count = 0
    with open(input_file_list, 'r') as f:
        for line in f:
            # line = E023_15_coreMarks_mnemonics.bed
            count += 1
            sample_id = line.split('/')[-1].split('_')[0]
            index_dict, state_dict = read_roadmap_bed(line.rstrip())
            sample_dict[sample_id] = index_dict
            sample_dict_state[sample_id] = state_dict
            if count % 10 == 0:
                print "Processed: %i tissues" % (count) 
    return sample_dict, sample_dict_state

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

def import_emission(input_emission_file):
    # key = emission state_num 
    # value = dict with histone mark - key = histone mark, value = emission probability

    emission_dict = {}
    with open(input_emission_file, 'r') as f:
        # state   H3K4me3 H3K4me1 H3K36me3    H3K9me3 H3K27me3
        
        for line in f:
            state_histone_emission_dict = {}
            state = int(line.split()[0].split('_')[0])
            state_histone_emission_dict["H3K4me3"] = float(line.split()[1])
            state_histone_emission_dict["H3K4me1"] = float(line.split()[2])
            state_histone_emission_dict["H3K36me6"] = float(line.split()[3])
            state_histone_emission_dict["H3K9me3"] = float(line.split()[4])
            state_histone_emission_dict["H3K27me3"] = float(line.split()[5].rstrip())
            emission_dict[state] = state_histone_emission_dict
    return emission_dict

######################
# Analysis Functions #
######################

def calculate_whole_sample_state_proportion(input_multiple_sample_state_dict):
    # INPUT ==> Key = sample, value = state_dict for state proportion calculatinon
    new_proportion_dict = {}
    for key, value in input_multiple_sample_state_dict.iteritems():
        sample = key 
        input_dict = value 
        proportion_state_dict = calculate_state_proportions(input_dict)
        new_proportion_dict[sample] = proportion_state_dict
    return new_proportion_dict



def calculate_state_proportions(input_state_dict):
    # INPUT = state_dict - key = state, value = num of bases in that state 
    total_number = 0
    for key, value in input_state_dict.iteritems():
        total_number += value 
    state_proportion_dict = {} 
    # print total_number
    for key, value in input_state_dict.iteritems():
        state_proportion_dict[key] = 1.0 * value / total_number
    return state_proportion_dict



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
    # active_state_list = [1,2,3,4,5,6,7,8] # all the active sites 
    active_state_list = [2,3,6,7,11,12] # H3K4me1 associated states 
    # this is noted by the consortium to be most tissue specific
    total_count = Counter()
    score_dict = defaultdict(float)
    master_list = []
    for snp, tissue_result_list in snp_tissue_dict.iteritems():
        count = Counter()
        active_tissue_list = []
        score_matrix = defaultdict(int)
        tissue_matrix_list = []
        for tissue_state_tuple in tissue_result_list:
            tissue = tissue_state_tuple[0]
            state = tissue_state_tuple[1]
            state_number = int(state.split('_')[0])
            if state_number in active_state_list:
                tissue_info = metadata[tissue]
                # print snp,tissue_info
                # total_count[tissue] += 1
                active_tissue_list.append(tissue)
                score_matrix[state_number] += 1 
            count[state] += 1

        for tissue in active_tissue_list:
            score_dict[tissue] += 1.0/len(active_tissue_list)

        tissue_matrix_list.append(snp)
        for tissue in active_tissue_list:
            tissue_matrix_list.append(score_matrix[tissue])
        master_list.append(tissue_matrix_list) 

        # print count 
    # print total_count


    # threshold > 
    print score_dict
    sorted_score_list = sorted(score_dict.iteritems(), key = lambda (k,v): v, reverse = True)
    for tissue_score_tuple in sorted_score_list:
        tissue = tissue_score_tuple[0]
        score = tissue_score_tuple[1]

        if score >= 0.9:
            print tissue, metadata[tissue], score 
    print master_list
    return sorted_score_list

def emission_matrix_SNP_by_tissue(metadata, snp_tissue_dict, emission_dict):
    # INPUT - metadata - tissue ID -> description 
    # INPUT - SNP tissue dict SNP(chr, pos) -> [(tissue, state)]
    # INPUT - emission dict - key = int(state), value = dict(key = histone modification, value = emission probability)

    histon_marker_list = ["H3K4me3", "H3K4me1", "H3K36me3","H3K9me3", "H3K27me3"]

    # get state_matrix 
    SNP_list = []


    state_matrix = pd.DataFrame()
    
    for SNP, tissue_result_list in snp_tissue_dict.iteritems():
        state_vector = pd.DataFrame(columns = [SNP])
        for tissue_state_tuple in tissue_result_list:
            tissue = tissue_state_tuple[0]
            state = tissue_state_tuple[1]
            temp_state = pd.DataFrame([state], index = [tissue], columns = [SNP])
            print temp_state
            state_vector = state_vector.append(temp_state)
        
        print state_vector
        state_matrix = pd.concat([state_matrix, state_vector], axis = 1) 

    print state_matrix
    sys.exit(33)






    total_count = Counter()
    score_dict = defaultdict(float)
    master_list = []
    histone_mark = "H3K4me1"
    for snp, tissue_result_list in snp_tissue_dict.iteritems():
        count = defaultdict(float)
        active_tissue_list = []
        score_matrix = defaultdict(float)
        tissue_matrix_list = []
        for tissue_state_tuple in tissue_result_list:
            tissue = tissue_state_tuple[0]
            state = tissue_state_tuple[1]
            state_number = int(state.split('_')[0])
            if tissue == "E044":
                print tissue, state, emission_dict[state_number][histone_mark] 
            tissue_info = metadata[tissue]
            score_dict[tissue] += emission_dict[state_number][histone_mark]
            # if emission_dict[state_number][histone_mark] > 0.5:
            #     print "+++++++", tissue, emission_dict[state_number][histone_mark]

        tissue_matrix_list.append(snp)
        for tissue in active_tissue_list:
            tissue_matrix_list.append(score_matrix[tissue])
        master_list.append(tissue_matrix_list) 

        # print count 
    # print total_count


    # threshold > 
    print score_dict
    sorted_score_list = sorted(score_dict.iteritems(), key = lambda (k,v): v, reverse = True)
    for tissue_score_tuple in sorted_score_list:
        tissue = tissue_score_tuple[0]
        score = tissue_score_tuple[1]

        if score >= 0.9:
            print tissue, metadata[tissue], score 
    # print master_list
    return sorted_score_list



def test_pipeline_T1D_simple():
    # read whole set 
    # TODO - make file 
    # input files 
    encode_bed_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_15core_marks_list.txt'
    T1D_SNP_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/T1D_snp_list_new.txt'
    metadata_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_metadata.txt'
    emission_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_emission.txt'
    # Read tissue files

    epigenetic_marks, state_dict = read_whole_sets(encode_bed_list)
    metadata_list = read_metadata(metadata_file)
    snp_list = read_gwas_SNP_hits_T1D(T1D_SNP_list)
    emission_dict = import_emission(emission_file)
    # print glance(epigenetic_marks)
    # print snp_list
    print metadata_list
    print emission_dict
    # sys.exit(2)
    # process
    #1. find block

    snp_result = t1d_match(epigenetic_marks, snp_list)
    # score_sorted_list = find_good_match(metadata_list, snp_result)
    proportion_state_dict = calculate_whole_sample_state_proportion(state_dict)
    
    score_sorted_list_with_emission = emission_matrix_SNP_by_tissue(metadata_list, snp_result, emission_dict)

    simple_output_file_score = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/simple_score_distribution_T1D.txt'
    with open(simple_output_file_score, 'w') as simple_outfile:
        for tissue_score_tuple in score_sorted_list_with_emission:
            # out_str = ""
            tissue = tissue_score_tuple[0]
            score = tissue_score_tuple[1]
            out_str = tissue + '\t' + str(score) + '\n'
            simple_outfile.write(out_str)
            print out_str
            out_str = ""

    print "simple enrichment finished"




def test_pipeline_T1D_novel():
    # read whole set 
    # TODO - make file 
    # input files 
    encode_bed_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_15core_marks_list.txt'
    T1D_SNP_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/T1D_snp_list.txt'
    metadata_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_metadata.txt'

    # Read tissue files

    epigenetic_marks, state_dict = read_whole_sets(encode_bed_list)
    metadata_list = read_metadata(metadata_file)
    snp_list = read_gwas_SNP_hits_T1D(T1D_SNP_list)

    # print glance(epigenetic_marks)
    # print snp_list
    print metadata_list

    # process
    #1. find proportion of the block 

    snp_result = t1d_match(epigenetic_marks, snp_list)
    #find_good_match(metadata_list, snp_result)
    proportion_state_dict = calculate_whole_sample_state_proportion(state_dict)
    print proportion_state_dict
    # make distribution for plots
    tissue_state_proportion_tsv = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_15core_marks_proportion_per_tissue.txt'
    tissue_count = 0
    with open(tissue_state_proportion_tsv, 'w') as outfile:
        for key, value in proportion_state_dict.iteritems():
            tissue_count += 1
            tissue_proportion = str(key) + "\t" + str(metadata_list[key])
            for i in range(1,16):
                tissue_proportion += '\t' + str(value[str(i)])
            if tissue_count % 10 == 0:
                print "exported proportion for %i tissues! " % tissue_count
            outfile.write(tissue_proportion)
    print "DONE! "


def test_pipeline_AD():
    # read whole set 
    # TODO - make file 
    # input files 
    encode_bed_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_15core_marks_list.txt'
    AD_SNP_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/AD_snp_list.txt'
    metadata_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_metadata.txt'

    # Read tissue files

    epigenetic_marks = read_whole_sets(encode_bed_list)
    metadata_list = read_metadata(metadata_file)

    snp_list = read_gwas_SNP_hits_T1D(AD_SNP_list)

    # print glance(epigenetic_marks)
    # print snp_list
    print metadata_list

    # process
    #1. find block

    snp_result = t1d_match(epigenetic_marks, snp_list)
    find_good_match(metadata_list, snp_result)

if __name__ == "__main__":
    test_pipeline_T1D_simple()
    # test_pipeline_T1D_novel()
    # test_pipeline_AD()