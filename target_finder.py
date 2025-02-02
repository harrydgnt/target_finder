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


def read_LD(input_LD_file):
    # INPUT LINE : SNP  Proxy   RSquared    Chromosome  Coordinate_HG18
    SNP_LD_Dict = {}
    SNP_hits = []
    previous_name = ""
    with open(input_LD_file, 'r') as f:
        SNP_hits = []
        for line in f:
            SNP_name = line.split()[0]
            if previous_name != SNP_name:
                print "SNP ChANGED! ", SNP_name 
                SNP_LD_Dict[previous_name] = SNP_hits
                SNP_hits = []
                previous_name = SNP_name
            proxy_name = line.split()[1]
            Rsquared = line.split()[2]
            chromosome = line.split()[3]
            position = line.split()[4].rstrip()
            SNP_hits.append((chromosome, position, Rsquared))
    return SNP_LD_Dict


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

def t1d_LD_match(input_roadmap_location_dict, input_LD_dict):
    # INPUT - Roadmap_dict - key = tissue ID, value = dict with chromosome -> list of positions
    # INPUT - SNP Hits list - (chromosome, position) tuple

    # Make list SNP first sorted by chromosome
    LD_matched_dict = {}
    num_SNP_counter = 0
    for index_SNP, LD_snp_list in input_LD_dict.iteritems():
        snp_list = sorted(LD_snp_list, key = itemgetter(0)) # sort by chromosome 
        # print snp_list
        # chr pos r2 
        # look for state in the samples
        num_SNP_counter += 1
        num_proxy_counter = 0
        SNP_search_result_dict = {}
        for snp in snp_list:
            num_proxy_counter += 1
            if num_proxy_counter % 10 == 0:
                print "processing index snp %i proxy snp %i" %(num_SNP_counter, num_proxy_counter)
            tissue_result_list = []
            chromosome = snp[0] 
            position = snp[1]
            rsquared = snp[2]
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
        LD_matched_dict[index_SNP] = SNP_search_result_dict
    return LD_matched_dict


def find_good_match(metadata, snp_tissue_dict, active):
    # INPUT - metadata - tissue ID -> description 
    # INPUT - SNP tissue dict SNP(chr, pos) -> [(tissue, state)]
    if active:
        active_state_list = [1,2,3,4,5,6,7,8] # all the active sites 
    else:
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
    num_state = 15
    # get state_matrix 
    SNP_list = []


    state_matrix = pd.DataFrame()
    
    for SNP, tissue_result_list in snp_tissue_dict.iteritems():
        state_vector = pd.DataFrame(columns = [SNP])
        for tissue_state_tuple in tissue_result_list:
            tissue = tissue_state_tuple[0]
            state = int(tissue_state_tuple[1].split('_')[0])
            temp_state = pd.DataFrame([state], index = [tissue], columns = [SNP])
            # print temp_state
            state_vector = state_vector.append(temp_state)
        
        # print state_vector
        state_matrix = pd.concat([state_matrix, state_vector], axis = 1) 

    print state_matrix
    # sys.exit(33)

    # make summary_matrix for each tissue
    summary_matrix = pd.DataFrame(columns = range(1,num_state + 1))
    for row_index, row in state_matrix.iterrows():
        in_list = [0] * num_state
        row_vector = pd.DataFrame(0, index = [row_index], columns = range(1,num_state + 1))
        for pos_state in row:
            row_vector[pos_state] += 1 
        print row_vector
        summary_matrix = summary_matrix.append(row_vector)
    print summary_matrix

    return state_matrix, summary_matrix




def emission_matrix_SNP_LD (metadata, emission_dict, LD_SNP_dict):
    # INPUT - metadata - tissue ID -> description 
    # INPUT - SNP tissue dict SNP(chr, pos) -> [(tissue, state)]
    # INPUT - emission dict - key = int(state), value = dict(key = histone modification, value = emission probability)

    histon_marker_list = ["H3K4me3", "H3K4me1", "H3K36me3","H3K9me3", "H3K27me3"]
    num_state = 15
    # get state_matrix 
    SNP_list = []

    # TODO - rewrite 

    state_matrix = pd.DataFrame()
    sum_matrix = pd.DataFrame()
    num_SNP = 0
    snp_dict = {}
    for index_SNP, snp_tissue_dict in LD_SNP_dict.iteritems():
        num_SNP += 1
        state_info_matrix = pd.DataFrame()
        for snp_chr_pos_rsqured_tuple, tissue_state_tuple_list in snp_tissue_dict.iteritems():
            print snp_chr_pos_rsqured_tuple
            snp = snp_chr_pos_rsqured_tuple
            snp_tissue_matrix = pd.DataFrame(columns = [snp])
            chromosome = snp_chr_pos_rsqured_tuple[0]
            position = snp_chr_pos_rsqured_tuple[1]
            rsquared_value = snp_chr_pos_rsqured_tuple[2]
            # new_snp_tuple = (index_SNP, chromosome, position, rsquared_value)
            count = 0
            for tissue_state_tuple in tissue_state_tuple_list:
                count += 1
                tissue = tissue_state_tuple[0]
                state = int(tissue_state_tuple[1].split('_')[0]) 
                new_vector = pd.DataFrame([state],index = [tissue], columns = [snp] )
                # if count % 100 == 0:
                #     print tissue, state, new_vector
                snp_tissue_matrix = snp_tissue_matrix.append(new_vector)
                # print "NEWVECTOR",new_vector
            # print snp_tissue_matrix
            state_info_matrix = pd.concat([state_info_matrix, snp_tissue_matrix], axis = 1)
        print "STATE_MATRIX", state_info_matrix


        print "processing LD calculation for snp %s" % num_SNP
        summary_matrix = pd.DataFrame(columns = range(1,num_state + 1))
        snp_specific_dict = {}
        for row_index, row in state_info_matrix.iterrows():
            in_list = [0] * num_state
            row_vector = pd.DataFrame(0, index = [row_index], columns = range(1,num_state + 1))
            for current_snp in list(state_info_matrix):
                pos_state = row[current_snp]
                rsq = current_snp[2]
                row_vector[pos_state] += 1.0 * float(rsq)
            print row_vector
            summary_matrix = summary_matrix.append(row_vector) 
        snp_dict[index_SNP] = summary_matrix
        sum_matrix = pd.concat([sum_matrix, summary_matrix], axis = 1)



        state_matrix = pd.concat([state_matrix, state_info_matrix], axis = 1)
    print sum_matrix, snp_dict
    return state_matrix, snp_dict



    # for index_SNP, snp_tissue_dict in LD_SNP_dict:
    #     for SNP, tissue_result_list in snp_tissue_dict.iteritems():
    #         state_vector = pd.DataFrame(columns = [SNP])
    #         for tissue_state_tuple in tissue_result_list:
    #             tissue = tissue_state_tuple[0]
    #             state = int(tissue_state_tuple[1].split('_')[0])
    #             temp_state = pd.DataFrame([state], index = [tissue], columns = [SNP])
    #             # print temp_state
    #             state_vector = state_vector.append(temp_state)
            
    #         # print state_vector
    #         state_matrix = pd.concat([state_matrix, state_vector], axis = 1) 

    #     print state_matrix
    #     # sys.exit(33)

    #     # make summary_matrix for each tissue
    #     summary_matrix = pd.DataFrame(columns = range(1,num_state + 1))
    #     for row_index, row in state_matrix.iterrows():
    #         in_list = [0] * num_state
    #         row_vector = pd.DataFrame(0, index = [row_index], columns = range(1,num_state + 1))
    #         for pos_state in row:
    #             row_vector[pos_state] += 1 
    #         print row_vector
    #         summary_matrix = summary_matrix.append(row_vector)
    #     print summary_matrix

    # return state_matrix, summary_matrix


def process_matrix(input_tissue_pos_state_matrix, simple):
    if simple:
        active_state_list = [1,2,3,4,5,6,7,8] # all the active sites 
    else:
        active_state_list = [3,6,7,12] # H3K4me1 associated states 
    new_matrix = pd.DataFrame()
    for tissue, pos_state_vector in input_tissue_pos_state_matrix.iterrows():
        print pos_state_vector.transpose()
        binarized_tissue = [ 1 if x in active_state_list else 0 for x in pos_state_vector]
        print list(input_tissue_pos_state_matrix)
        print binarized_tissue
        print tissue 
        new_vector = pd.DataFrame([binarized_tissue], index = [tissue], columns = list(input_tissue_pos_state_matrix))
        new_matrix = pd.concat([new_matrix, new_vector], axis = 0)
    return new_matrix 

def normalize_matrix(input_binary_matrix):
    new_matrix = input_binary_matrix.copy().astype(float)
    for i in range(len(list(input_binary_matrix))):
        column = input_binary_matrix[list(input_binary_matrix)[i]]
        total_score = 0
        for entry in column:
            total_score += entry
        for j in range(len(column)):
            if total_score == 0.0:
                total_score = 1.0

            print i,j, new_matrix.index[j], list(new_matrix)[i], new_matrix[list(new_matrix)[i]][new_matrix.index[j]]*1.0/total_score
            new_matrix.set_value(new_matrix.index[j], list(new_matrix)[i], new_matrix[list(new_matrix)[i]][new_matrix.index[j]]*1.0/total_score)
    return new_matrix

def sum_rows_in_matrix(input_normalized_matrix):
    sum_dict = {}
    for tissue, score_vector in input_normalized_matrix.iterrows():
        tissue_score = 0.0
        for score in score_vector:
            tissue_score += score
        sum_dict[tissue] = tissue_score
    return sum_dict


def test_pipeline_T1D_simple():
    # read whole set 
    # TODO - make file 
    # input files 
    encode_bed_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_15core_marks_list.txt'
    T1D_SNP_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/T1D_snp_list_new.txt'
    metadata_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_metadata.txt'
    emission_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_emission.txt'
    # Read tissue files
    wrdir = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/'



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
    score_sorted_list = find_good_match(metadata_list, snp_result, True)
    # proportion_state_dict = calculate_whole_sample_state_proportion(state_dict)
    
    tissue_pos_state_matrix, summary_matrix = emission_matrix_SNP_by_tissue(metadata_list, snp_result, emission_dict)
    simple_processed_matrix = process_matrix(tissue_pos_state_matrix, True)
    simple_processed_matrix.to_csv(wrdir+"T1D_binarized_tissue_pos_state_matrix.txt", sep = '\t')

    marker_processed_matrix = process_matrix(tissue_pos_state_matrix, False)
    marker_processed_matrix.to_csv(wrdir + "T1D_binarized_H3K4Me1_pos_state_matrix.txt", sep = '\t')

    overlap_processed_matricis = simple_processed_matrix.add(marker_processed_matrix, fill_value = 0)
    overlap_processed_matricis.to_csv(wrdir + "T1D_binary_active_H3K4Me1_matrix.txt", sep = '\t')

    simple_normalized_matrix = normalize_matrix(simple_processed_matrix)
    print simple_normalized_matrix
    simple_sum_dict = sum_rows_in_matrix(simple_normalized_matrix)

    marker_normalized_matrix = normalize_matrix(marker_processed_matrix)
    marker_sum_dict = sum_rows_in_matrix(marker_normalized_matrix)
    
    simple_sum_dict_sorted = sorted(simple_sum_dict.iteritems(), key = lambda (k,v): v, reverse = True)
    marker_sum_dict_sorted = sorted(marker_sum_dict.iteritems(), key = lambda (k,v): v, reverse = True)

    for item in simple_sum_dict_sorted[0:10]:
        tissue = item[0]
        score = item[1]
        tissue_info = metadata_list[tissue]
        print tissue, score, tissue_info
    for item in marker_sum_dict_sorted[0:10]:
        tissue = item[0]
        score = item[1]
        tissue_info = metadata_list[tissue]
        print tissue, score, tissue_info
    # print simple_sum_dict_sorted
    # print marker_sum_dict_sorted

    # print simple_processed_matrix
    sys.exit(24)
    simple_output_file_score = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/simple_score_distribution_T1D.txt'
    with open(simple_output_file_score, 'w') as simple_outfile:
        for tissue_score_tuple in score_sorted_list:
            # out_str = ""
            tissue = tissue_score_tuple[0]
            score = tissue_score_tuple[1]
            out_str = tissue + '\t' + str(score) + '\n'
            simple_outfile.write(out_str)
            print out_str
            out_str = ""


    tissue_pos_state_matrix.to_csv(wrdir + "T1D_tissue_pos_state_matrix.txt", sep = '\t')
    summary_matrix.to_csv(wrdir + "T1D_summary_matrix.txt", sep = '\t')


    print "simple enrichment finished"


def find_max_state(input_idxsnp_matrix_dict):
    new_snp_tissue_state_dict = {}
    for idx_snp, summary_matrix in input_idxsnp_matrix_dict.iteritems():
        tissue_state_dict = {}
        for tissue, row in summary_matrix.iterrows():
            tissue_state_dict[tissue] = row.idxmax(axis = 1)
        new_snp_tissue_state_dict[idx_snp] = tissue_state_dict
    return new_snp_tissue_state_dict

def process_new_score_dict(input_snp_state_tissue_score_dict):

    heatmap_matrix = pd.DataFrame()
    active_state_list = [1,2,3,4,5,6,7,8]
    for snp, state_dict in input_snp_state_tissue_score_dict.iteritems():
        new_vector = pd.DataFrame(columns = [snp])
        for tissue, state in state_dict.iteritems():
            # binarize
            dummy = -1
            if int(state) in active_state_list:
                dummy = 1
            else:
                dummy = 0
            new_entry = pd.DataFrame([dummy], columns = [snp], index = [tissue])
            new_vector = new_vector.append(new_entry)
        heatmap_matrix = pd.concat([heatmap_matrix,new_vector], axis = 1)




    return heatmap_matrix

def test_pipeline_T1D_novel():
    # read whole set 
    # TODO - make file 
   # input files 
    wrdir = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/'

    encode_bed_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_15core_marks_list.txt'
    T1D_SNP_list = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/T1D_snp_list_new.txt'
    metadata_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_metadata.txt'
    emission_file = '/Users/harryyang/Documents/Research/Class/Com Sci 225/target_finder/roadmap_emission.txt'
    LD_file = wrdir + 'T1D_LD_snps.txt'

    # Read tissue files

    epigenetic_marks, state_dict = read_whole_sets(encode_bed_list)
    metadata_list = read_metadata(metadata_file)
    # snp_list = read_gwas_SNP_hits_T1D(T1D_SNP_list)
    LD_SNP_dict = read_LD(LD_file)
    # for key, value in LD_SNP_dict.iteritems():
    #     print key, value 
    # sys.exit(22)
    emission_dict = import_emission(emission_file)

    # print glance(epigenetic_marks)
    # print snp_list
    # print metadata_list


    # tissue_pos_state_matrix, summary_matrix = emission_matrix_SNP_by_tissue(metadata_list, snp_result, emission_dict)

    # process
    #1. match block <-> positions
    # 
    snp_result = t1d_LD_match(epigenetic_marks, LD_SNP_dict)
    print snp_result

    # 2. Process them for LD aware states
    # 2.1. calculate score_matrix
    # def emission_matrix_SNP_LD (metadata, emission_dict, LD_SNP_dict):

    total_matrix, result_snp_dict = emission_matrix_SNP_LD(metadata_list, emission_dict, snp_result)
    new_score_dict = find_max_state(result_snp_dict)
    print new_score_dict
    new_score_dict = dict((k,v) for k, v in new_score_dict.iteritems() if v)
    heatmap_matrix= process_new_score_dict(new_score_dict)
    print heatmap_matrix
    heatmap_matrix.to_csv(wrdir+'/T1D_LD_aware_heatmap_matrix.txt', sep = 't')


    normalized_ld_matrix = normalize_matrix(heatmap_matrix)
    score = sum_rows_in_matrix(normalized_ld_matrix)

    sorted_list = sorted(score.iteritems(), key = lambda (k,v): v, reverse = True)
    print sorted_list
    for item in sorted_list[0:10]:
        tissue = item[0]
        score = item[1]
        tissue_info = metadata_list[tissue]
        print tissue, score, tissue_info
    sys.exit(23)
    #find_good_match(metadata_list, snp_result)
    # proportion_state_dict = calculate_whole_sample_state_proportion(state_dict)
    # print proportion_state_dict
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
    # test_pipeline_T1D_simple()
    test_pipeline_T1D_novel()
    # test_pipeline_AD()