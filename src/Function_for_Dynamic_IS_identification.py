###Header################################################
header = """

+------------------------------------------------------+
 Module: Function_for_Dynamic_IS_identification
 Author: Stefano Brasca
 Date:  September 23th, 2014
 Contact: brasca.stefano@hsr.it
 Version: inDevelopment
+------------------------------------------------------+

 Description:
  - This module contains functions used in Dynamic IS
    identification framework
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import######
import sys
import os
import random
from subprocess import call
from numpy.random import multinomial
import math
import multiprocessing
from operator import itemgetter
from operator import attrgetter
import scipy
####################################

###Import Module(s)###########################
import Function_for_Gaussian_IS_identification
import DB_connection
import Classes_for_Integration_Analysis
##############################################




##########################################################################################################
### RANKING HISTOGRAMS ###################################################################################
##########################################################################################################

def get_alpha (N_sigma, interaction_limit):
    '''
    given a window of 2*interaction_limit+1, this function yields the proper
    alpha parameter in order to create a gaussian histogram with N_sigma falling
    inside the window (es N_sigma=2 -> 95% of distr inside the window)
    [for Function_for_Gaussian_IS_identification.gaussian_histogram_generator callings]
    '''
    
    span = 2*interaction_limit+1   # def    
    alpha = N_sigma / float(span)   # resulted from exact calculations
    
    return alpha
    

def get_sigma_path (MAX_interaction_limit, MIN_interaction_limit, N_sigma):
    '''
    given the interaction limit range to explore, this function yields the list
    of (interaction_limit, alpha) points that respect N_sigma constraint
    [for Function_for_Gaussian_IS_identification.gaussian_histogram_generator callings]
    
    NOTE: MAX_interaction_limit must be >= MIN_interaction_limit; if =, only one 
    point will be returned
    '''
    
    sigma_path_point_list = []
    interaction_limit_range = range(MIN_interaction_limit, MAX_interaction_limit+1)
    
    for interaction_limit in interaction_limit_range:
        alpha = get_alpha(N_sigma, interaction_limit)
        sigma_path_point = (interaction_limit, alpha)
        sigma_path_point_list.append(sigma_path_point)
        
    return sigma_path_point_list   # list of (interaction_limit, alpha) tuples
    
    
def get_ranking_histograms (MAX_interaction_limit, sigma_paths, MIN_interaction_limit, adaptive_SW):
    '''
    INPUT: MAX_interaction_limit, MIN_interaction_limit - int
           sigma_paths - LIST of kind [2,3,...] (sigma path to explore)
           adaptive_SW - boolean
    OUTPUT: ranking_histogram_dict_list - list of dictionary entries, descripted below
    NOTE: - MAX_interaction_limit must be >= MIN_interaction_limit;
          - sigma_paths must be always a list, even if only one sigma value has to be exploited
    
    DESCRIPTION:
    this function implement 'get_sigma_path' above and gaussian_histogram_generator and 
    normalize_histogram_to_the_peak from Function_for_Gaussian_IS_identification module.
    it yields a list of dictionary of kind:
        
        dict_entry = {'type': "gaussian selection",
                      'N_sigma': N_sigma,
                      'interaction_limit': interaction_limit,
                      'alpha': alpha,
                      'hist_gauss': bin_areas,
                      'bin_boundaries': bin_boundaries,
                      'excluded_perc': diagnostic*100,
                      'hist_normalized_to_peak': hist_gauss_normalized_to_peak}
                      
    if adaptive_SW is True (this name was chosen thinkig at the expected behaviuor of the algorithm),
    further dictionaries will be added to the returned list, derived from limit cases (alpha ->0 or 
    alpha >> 1)
    
        dict_entry = {'type': "adaptive SW",
                              'interaction_limit': interaction_limit,
                              'bin_boundaries': bin_boundaries,
                              'hist_normalized_to_peak': bin_areas}
                              
    NOTE: limit cases of alpha -> 0 are one for each different interaction_limit while for alpha >> 1
    the case is only one (and can be assimilated to the (alpha -> 0, interaction_limit=3) case)
    '''

    ###############################################################################################
    def limit_histogram_generator (interaction_limit):
        # Bin Boundaries Half
        first_bin_half = (0.0, 0.5)
        bin_boundaries_half = [first_bin_half]        
        for i in range(0, interaction_limit):
            bin_boundaries_half.append((bin_boundaries_half[-1][1], bin_boundaries_half[-1][1]+1))        
        # Filling Bins
        bin_areas_half = []
        for x,y in bin_boundaries_half:
            bin_areas_half.append(1.0)
        # Symmetrize Bins
        n_bins = int((2*interaction_limit)+1)
        bin_boundaries = [None]*n_bins
        bin_areas = [None]*n_bins
        bin_boundaries[interaction_limit] = (-0.5, 0.5)
        bin_areas[interaction_limit] = bin_areas_half[0]
        for i in range(1, interaction_limit +1):
            negative_left = -1.0*bin_boundaries_half[interaction_limit+1-i][1]
            negative_right = -1.0*bin_boundaries_half[interaction_limit+1-i][0]
            bin_boundaries[i-1] = (negative_left, negative_right)
            bin_boundaries[i+interaction_limit] = bin_boundaries_half[i]
            bin_areas[i-1] = bin_areas_half[interaction_limit+1-i]
            bin_areas[i+interaction_limit] = bin_areas_half[i]            
        # Remember that bin_boundaries[interaction_limit] is the peak        
        return bin_boundaries, bin_areas
    ###############################################################################################
    
    
    if ((sigma_paths is None) and (adaptive_SW is False)):
        print "\n\n\t get_ranking_histograms function can't run, due to sigma_paths and adaptive_SW variables both set to 'None' and 'False'."
        sys.exit("\n\n\t[ERROR]\tQuit.\n\n")
        
    else:
        
        ranking_histogram_dict_list = []
        
        if (adaptive_SW is True):
            for interaction_limit in range(MIN_interaction_limit, MAX_interaction_limit+1):
                bin_boundaries, bin_areas = limit_histogram_generator (interaction_limit)
                dict_entry = {'type': "adaptive SW",
                              'interaction_limit': interaction_limit,
                              'bin_boundaries': bin_boundaries,
                              'hist_normalized_to_peak': bin_areas}
                ranking_histogram_dict_list.append(dict_entry)
        
        if (sigma_paths is not None):            
            for N_sigma in sigma_paths:
                sigma_path_point_list = get_sigma_path(MAX_interaction_limit, MIN_interaction_limit, N_sigma)
                for interaction_limit, alpha in sigma_path_point_list:
                    bin_boundaries, bin_areas, diagnostic = Function_for_Gaussian_IS_identification.gaussian_histogram_generator(interaction_limit, alpha)
                    hist_gauss_normalized_to_peak = Function_for_Gaussian_IS_identification.normalize_histogram_to_the_peak(bin_areas, interaction_limit)
                    dict_entry = {'type': "gaussian selection",
                                  'N_sigma': N_sigma,
                                  'interaction_limit': interaction_limit,
                                  'alpha': alpha,
                                  'hist_gauss': bin_areas,
                                  'bin_boundaries': bin_boundaries,
                                  'excluded_perc': diagnostic*100,
                                  'hist_normalized_to_peak': hist_gauss_normalized_to_peak}
                    ranking_histogram_dict_list.append(dict_entry)
                    
        return ranking_histogram_dict_list
                                  
                


##########################################################################################################
### EXPLORE AND CLUSTER SOLUTIONS ########################################################################
##########################################################################################################

def check_consistency (ISs_and_configDict_couple_list):
    
    # Check: if all different results are consistent... skip!
    # DEF: consistency <-> same N_IS, same CB 	allocation
    consistency = False
    # Print for Devel
    #print "MULTIPLE SOLUTIONS -> ",
    # same_N_IS
    temp_set = set()
    for IS_list, configDict in ISs_and_configDict_couple_list:
        temp_set.add(len(IS_list))
    if (len(temp_set) == 1):
        consistency = True
        # Print for Devel
        #print "consistency test 1 passed -> ",
    # same CB allocation
    temp_superlist = list()
    if (consistency is True):
        temp_superlist = [IS_list for IS_list, configDict in ISs_and_configDict_couple_list]
        n_superlist = len(temp_superlist)
        for i in range(0,n_superlist-1):
            if (consistency is False):
                break
            else:
                for j in range(i+1, n_superlist):
                    if (consistency is False):
                        break
                    else:
                        IS_list_i = temp_superlist[i]
                        IS_list_j = temp_superlist[j]
                        for IS_i, IS_j in zip(IS_list_i, IS_list_j):  # list of same len due to former control
                            cb_list_i = IS_i.Covered_bases_list
                            cb_list_j = IS_j.Covered_bases_list
                            if (len(cb_list_i) != len(cb_list_j)):
                                consistency = False
                                break
                            else:
                                for cb_i, cb_j in zip(cb_list_i, cb_list_j):
                                    if (cb_i is not cb_j):
                                        consistency = False
                                        break
    return consistency  # Boolean


def solution_clustering (ISs_and_configDict_couple_list):
    
    #Result collector
    putative_unique_solution_list = []
    
    #Loop control ('row' suppression)
    i_indexes = range(0, len(ISs_and_configDict_couple_list))
    #N of different types
    n_of_distinct_putative_unique_solution = 0
#    # Dev check var
#    check_var = 0
    
    #Clustering loop
    for i in range(0, len(ISs_and_configDict_couple_list)):
        if i in i_indexes:
            putative_unique_solution_list.append(Classes_for_Integration_Analysis.Putative_unique_solution(ISs_and_configDict_couple_list[i]))
            n_of_distinct_putative_unique_solution += 1
            for j in range(i+1, len(ISs_and_configDict_couple_list)):
                pair = [ISs_and_configDict_couple_list[i], ISs_and_configDict_couple_list[j]]
                consistency = check_consistency (pair)
                if consistency is True:
                    putative_unique_solution_list[n_of_distinct_putative_unique_solution-1].join(ISs_and_configDict_couple_list[j])
                    i_indexes.remove(j)
#                    check_var += 1
                    
#    # Dev Checks
#    if n_of_distinct_putative_unique_solution > len(ISs_and_configDict_couple_list):
#        sys.exit("\n\n\t[ERROR A]\tQuit.\n\n")
#    if (len(ISs_and_configDict_couple_list) != (n_of_distinct_putative_unique_solution + check_var)):
#        print "len(ISs_and_configDict_couple_list) = ", len(ISs_and_configDict_couple_list)
#        print "n_of_distinct_putative_unique_solution = ", n_of_distinct_putative_unique_solution
#        print "check_var = ", check_var
#        sys.exit("\n\n\t[ERROR B]\tQuit.\n\n")
#    cardinality_sum = 0
#    for item in putative_unique_solution_list:
#        cardinality_sum += item.cardinality
#    if cardinality_sum != len(ISs_and_configDict_couple_list):
#        sys.exit("\n\n\t[ERROR C]\tQuit.\n\n")

    # Return results
    return putative_unique_solution_list
            




##########################################################################################################
### SIMULATIONS - Preliminary Tasks#######################################################################
##########################################################################################################

def get_ID (integration_analysis_object):
    
    object_id = str(integration_analysis_object)
    
    object_id = object_id.split('.')[1]
    object_id = object_id[:-1]
    object_id = object_id.split(' instance at ')
    
    my_id = '_'.join(object_id)
    
    return my_id
    
def reverse_complement (sequence):
    
    rev_seq = sequence[::-1]
    rev_comp_seq = ""
    for char in rev_seq:
        if char == 'A' or char == 'a':
            if char == 'A':
                rev_comp_seq += 'T'
            else:
                rev_comp_seq += 't'
            continue
        if char == 'T' or char == 't':
            if char == 'T':
                rev_comp_seq += 'A'
            else:
                rev_comp_seq += 'a'
            continue
        if char == 'C' or char == 'c':
            if char == 'C':
                rev_comp_seq += 'G'
            else:
                rev_comp_seq += 'g'
            continue
        if char == 'G' or char == 'g':
            if char == 'G':
                rev_comp_seq += 'C'
            else:
                rev_comp_seq += 'c'
            continue
        else:
            rev_comp_seq += char
            
    return rev_comp_seq
            

def analyze_sequences (header_list, conn, seqTracker_conn_dict):
    
    ### DB DATA RETRIEVAL ###
    
    # Retrieve sequences and metadata from DB - iss_table
    table_kind = 'IS'
    final_read_and_metadata_dictionary = DB_connection.retrieve_sequences_and_metadata_from_DB (conn, seqTracker_conn_dict['iss_table'], table_kind, header_list)
    # final_read_and_metadata_dictionary is a dictionary with entries like:
        # {read_header: sub_dictionary}
        # each sub_dictionary is a dictionary with entries like:
            # {'isread_cigar': CIGAR string, 'isread_MD': MD string, 'isread_nasequence': IS (or 'final') sequence, 'seq_len': lenght of the sequence, 'isread_strand': strand}
        
    # Retrieve raw sequences from DB - raw_table
    # NOTE: first 20 nucleotides are removed inside retrieve_sequences_and_metadata_from_DB func.
    table_kind = 'RAW'
    raw_read_dictionary = DB_connection.retrieve_sequences_and_metadata_from_DB (conn, seqTracker_conn_dict['raw_table'], table_kind, header_list)
    # raw_read_dictionary has entries like:
        # {read_header: raw sequence[20:]}
        
    # SUMMARY OF DATA ARRANGEMENT
    # Now header_list is the list of key for all the dictionaries
    # raw_read_dictionary[key] -> raw sequence[20:]
    # final_read_and_metadata_dictionary[key] -> sub_dict with keys: isread_cigar, isread_MD, isread_nasequence, seq_len, isread_strand
    
    
    ### LTR / untrimmedLC COMPUTATION  #  LTR_LC_dictionary ###  --> Returned
    
    # Debug variables
    bad_seq_headers = []
    rev_comp_headers = []
    
    # LTR_LC_dictionary
    LTR_LC_dictionary = {}
    for header in header_list:
        LTR_sequence = None
        untrimmed_LC_sequence = None
        # LTR_sequence, untrimmed_LC_sequence = raw_read_dictionary[header].split(final_read_and_metadata_dictionary[header]['isread_nasequence'])
        sequence_tupla = raw_read_dictionary[header].split(final_read_and_metadata_dictionary[header]['isread_nasequence'])
        if len(sequence_tupla) > 1:
            LTR_sequence = sequence_tupla[0]
            untrimmed_LC_sequence = sequence_tupla[1]
        else:  # for data from old platform
            RC_raw = reverse_complement (raw_read_dictionary[header])
            sequence_tupla = RC_raw.split(final_read_and_metadata_dictionary[header]['isread_nasequence'])
            if len(sequence_tupla) > 1:
                LTR_sequence = sequence_tupla[0]
                untrimmed_LC_sequence = sequence_tupla[1]
                rev_comp_headers.append(header)
            else:
                LTR_sequence = 'ACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA'
                untrimmed_LC_sequence = ''
                bad_seq_headers.append(header)
            
        LTR_len = len(LTR_sequence)
        untrimmed_LC_len = len(untrimmed_LC_sequence)
        if untrimmed_LC_sequence == '':
            untrimmed_LC_sequence = None
        sub_dictionary = {'LTR_sequence': LTR_sequence, 'LTR_len': LTR_len, 'untrimmed_LC_sequence': untrimmed_LC_sequence, 'untrimmed_LC_len': untrimmed_LC_len, 'strand': final_read_and_metadata_dictionary[header]['isread_strand']}
        LTR_LC_dictionary[header] = sub_dictionary
    # LTR_LC_dictionary has entries like:
        # {read_header: sub_dictionary}
        # each sub_dictionary is a dictionary with entries like:
            # {'LTR_sequence': LTR sequence, 'LTR_len': lenght of LTR sequence, 'untrimmed_LC_sequence': sequence of untrimmed LC (may be None), 'untrimmed_LC_len': lenght of LC sequence (always a number, 0 if None), 'strand': strand}
            
    ### CREATE DICTIONARY_FOR_SEQUENCE_SIMULATIONS ###  --> Returned
    dictionary_for_sequence_simulations = {}
    for header in header_list:
        final_read_and_metadata_subdictionary = final_read_and_metadata_dictionary[header]
        # Parse CIGAR for number of Ins and Del
        I = 0  #insertions
        D = 0  #deletions
        for char in final_read_and_metadata_subdictionary['isread_cigar']:
            if ((char == 'I') or (char == 'D')):
                if char == 'I':
                    I += 1
                else:
                    D += 1
        # Parse MD tag for number of Mut
        M = 0  #mutations
        prev_char = None
        for char in final_read_and_metadata_subdictionary['isread_MD']:
            if prev_char != '^':
                prev_char = char
                continue
            else:
                prev_char = char
                if (char in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']):
                    M += 1
        # Correct seq_len for Ins and Del
        len_var = D - I  #NB it's correct: i need a I-shorter sequence; the opposite for Del.
        seq_len = final_read_and_metadata_subdictionary['seq_len'] + len_var
        
        # Prepare sub_dict
        strand = final_read_and_metadata_subdictionary['isread_strand']
        sub_dict = {'strand': strand, 'seq_len': seq_len, 'Mut': M, 'Ins': I, 'Del': D}
        # Fill dictionary_for_sequence_simulations
        dictionary_for_sequence_simulations[header] = sub_dict
        
    # dictionary_for_sequence_simulations is a dictionary with entries like:
        # {read_header: sub_dictionary}
        # each sub_dictionary is a dictionary with entries like:
            # {'seq_len': lenght of the sequence TO RETRIEVE FROM REFERENCE GENOME, 'strand': strand, 'Mut': #Mutations, 'Ins': #Insertions, 'Del': #Deletions}
    
    # Debug print
    if len(rev_comp_headers) != 0 or len(bad_seq_headers) != 0:
        print "  \t+++++++++++++++++++++++++++++++++++++++++ "
        print "  \tBAD SEQUENCES FOUND: {}".format(str(len(rev_comp_headers)+len(bad_seq_headers)))
        print "  \tSolved through RC: {}".format(str(len(rev_comp_headers)))
        if len(bad_seq_headers) != 0:
            print "  \tUNSOLVED: {}".format(str(len(bad_seq_headers)))
            print "  \tUNSOLVED HEADER LIST: {}".format(str(bad_seq_headers))
    
    # Return dictionary_for_sequence_simulations, LTR_LC_dictionary
    return dictionary_for_sequence_simulations, LTR_LC_dictionary
    

def get_assembly_path (reference_genome):
    '''
    input: assembly label like 'hg19'
    output: assembly path if available or sys.exit
    '''
    assembly_path = None
    if reference_genome == 'hg19':
        assembly_path = '/opt/genome/human/hg19/index/bwa_7/hg19.fa'
    elif reference_genome == 'mm9':
        assembly_path = '/opt/genome/mouse/mm9/index/bwa_7/mm9.fa'
    else:
        return False
    return assembly_path
    

    
def get_seq_from_ref (Putative_unique_solution_object, dictionary_for_sequence_simulations, reference_genome, perfect_sequence_folder_path):
    '''
    Putative_unique_solution_object.perfect_sequence_dict = perfect_sequence_dict
    where:
        perfect_sequence_dict = dict with entries {'header': seq}
        perfect_sequence_strandness_dict = dict with entries {'header': strand}
    '''
    
    # Get IS_list
    IS_list = Putative_unique_solution_object.IS_list
    
    perfect_sequence_strandness_dict = {}
    # Prepare BED file for retrieve sequences from reference_genome_assembly
    bed_file_lines = []
    for IS in IS_list:
        chromosome = "chr{}".format(str(IS.chromosome))
        integration_locus = IS.integration_locus
        for header in IS.reads_key_list:
            sub_dict = dictionary_for_sequence_simulations[header]
            line = [chromosome]
            strand = sub_dict['strand']
            perfect_sequence_strandness_dict[header] = strand
            seq_len = sub_dict['seq_len']
            if ((strand == '+') or (strand == '1')):
                strand = '+'
                line.append(str(integration_locus))
                line.append(str(integration_locus+seq_len))
            elif ((strand == '-') or (strand == '2')):
                strand = '-'
                line.append(str(integration_locus-seq_len))
                line.append(str(integration_locus))
            line.append(header)  # Line ID --> fasta file line header!!
            line.append("0") #fake score
            line.append(strand) #strand
            tsv_line = '\t'.join(line)
            bed_file_lines.append(tsv_line)
    
    # Write BED file
    bed_file_name = "ISs_bedfile.bed"
    bedfile_complete_path = os.path.normpath(os.path.join(perfect_sequence_folder_path, bed_file_name))
    with open(bedfile_complete_path, 'w') as bed_file_pointer:
        bed_file_pointer.write('\n'.join(bed_file_lines)+'\n')
        
    # Call bedtools getfasta
    assembly_path = get_assembly_path (reference_genome)
    if assembly_path == False:
        sys.exit("\n\n\t[ERROR] Can't find assembly anymore for genome '{}'\tQuit.\n\n".format(reference_genome))
    fasta_file_name = "SeqFromRef.fa"
    fastafile_complete_path = os.path.normpath(os.path.join(perfect_sequence_folder_path, fasta_file_name))
    command = ['fastaFromBed', '-fi', assembly_path, '-bed', bedfile_complete_path, '-s', '-fo', fastafile_complete_path, '-name']
    exit_status = call(command)
    
    # Read FASTA file --> perfect_sequence_list
    if exit_status != 0:
        sys.exit("\n\n\t[ERROR] Some troubles occured in 'fastaFromBed' call.\tQuit.\n\n")
    
    perfect_sequence_dict = {}
    with open(fastafile_complete_path, 'r') as source_file:
        
        # Load file as a list --> source_file_lines
        source_file_lines = source_file.readlines()
        # Get file features / quick check
        n_of_seq = len(source_file_lines)
        if n_of_seq == 0:
            sys.exit("\n\n\t[ERROR] FASTA file turned out to be empty!\tQuit.\n\n")
        elif n_of_seq%2 != 0:
            sys.exit("\n\n\t[ERROR] FASTA file has an odd number of lines!\tQuit.\n\n")
        else:
            n_of_seq = int(n_of_seq/2.0)
        # Loop over source file lines
        header = None
        for source_file_line in source_file_lines:
            if '>' in source_file_line:
                header = source_file_line.rstrip()[1:]
            else:
                sequence = source_file_line.rstrip()
                perfect_sequence_dict[header] = sequence
                
    # Store perfect_sequence_list in Putative_unique_solution_object
    Putative_unique_solution_object.perfect_sequence_dict = perfect_sequence_dict
    Putative_unique_solution_object.perfect_sequence_strandness_dict = perfect_sequence_strandness_dict
    
    
def get_seq_MID_dict_list (Putative_unique_solution_object, dictionary_for_sequence_simulations):
    """
    Putative_unique_solution_object.seq_MID_dict_list = seq_MID_dict_list
    where:
        seq_MID_dict_list = [{'M':numM, 'I':numI, 'D':numD}, {...}, ... ]
        paired with Putative_unique_solution_object.IS_list
    """
    
    # Get IS_list
    IS_list = Putative_unique_solution_object.IS_list
    # Prepare seq_MID_dict_list
    seq_MID_dict_list = []
    
    for IS in IS_list:
        MID_dict = {'M': 0, 'I': 0, 'D': 0}
        for header in IS.reads_key_list:
            sub_dict = dictionary_for_sequence_simulations[header]
            if sub_dict['Mut'] != 0:
                MID_dict['M'] += sub_dict['Mut']
            if sub_dict['Ins'] != 0:
                MID_dict['I'] += sub_dict['Ins']
            if sub_dict['Del'] != 0:
                MID_dict['D'] += sub_dict['Del']
        seq_MID_dict_list.append(MID_dict)
        
    # Store seq_MID_dict_list in Putative_unique_solution_object
    Putative_unique_solution_object.seq_MID_dict_list = seq_MID_dict_list
    

##########################################################################################################
### SIMULATIONS ##########################################################################################
##########################################################################################################

def RandFloats(Size):
    Scalar = 1.0
    VectorSize = Size
    RandomVector = [random.random() for i in range(VectorSize)]
    RandomVectorSum = sum(RandomVector)
    RandomVector = [Scalar*i/RandomVectorSum for i in RandomVector]
    return RandomVector
    
def RandIntVec(ListSize, ListSumValue, Distribution='normal'):
    """
    Inputs:
    ListSize = the size of the list to return
    ListSumValue = The sum of list values
    Distribution = can be 'uniform' for uniform distribution, 'normal' for a normal distribution ~ N(0,1) with +/- 5 sigma  (default), or a list of size 'ListSize' or 'ListSize - 1' for an empirical (arbitrary) distribution. Probabilities of each of the p different outcomes. These should sum to 1 (however, the last element is always assumed to account for the remaining probability, as long as sum(pvals[:-1]) <= 1).  
    Output:
    A list of random integers of length 'ListSize' whose sum is 'ListSumValue'.
    
    USAGE: result = RandIntVec(ListSize, ListSumValue, Distribution=RandFloats(ListSize))
    """
    if type(Distribution) == list:
        DistributionSize = len(Distribution)
        if ListSize == DistributionSize or (ListSize-1) == DistributionSize:
            Values = multinomial(ListSumValue,Distribution,size=1)
            OutputValue = Values[0]
    elif Distribution.lower() == 'uniform': #I do not recommend this!!!! I see that it is not as random (at least on my computer) as I had hoped
        UniformDistro = [1/ListSize for i in range(ListSize)]
        Values = multinomial(ListSumValue,UniformDistro,size=1)
        OutputValue = Values[0]
    elif Distribution.lower() == 'normal':
        
        """
        Normal Distribution Construction....It's very flexible and hideous
        Assume a +-3 sigma range.  Warning, this may or may not be a suitable range for your implementation!
        If one wishes to explore a different range, then changes the LowSigma and HighSigma values
        """
        LowSigma = -3  #-3 sigma
        HighSigma = 3  #+3 sigma
        StepSize = 1/(float(ListSize) - 1)
        ZValues = [(LowSigma * (1-i*StepSize) +(i*StepSize)*HighSigma) for i in range(int(ListSize))]
        #Construction parameters for N(Mean,Variance) - Default is N(0,1)
        Mean = 0
        Var = 1
        #NormalDistro= [self.NormalDistributionFunction(Mean, Var, x) for x in ZValues]
        NormalDistro= list()
        for i in range(len(ZValues)):
            if i==0:
                ERFCVAL = 0.5 * math.erfc(-ZValues[i]/math.sqrt(2))
                NormalDistro.append(ERFCVAL)
            elif i ==  len(ZValues) - 1:
                ERFCVAL = NormalDistro[0]
                NormalDistro.append(ERFCVAL)
            else:
                ERFCVAL1 = 0.5 * math.erfc(-ZValues[i]/math.sqrt(2))
                ERFCVAL2 = 0.5 * math.erfc(-ZValues[i-1]/math.sqrt(2))
                ERFCVAL = ERFCVAL1 - ERFCVAL2
                NormalDistro.append(ERFCVAL)  
        #print "Normal Distribution sum = %f"%sum(NormalDistro)
        Values = multinomial(ListSumValue,NormalDistro,size=1)
        OutputValue = Values[0]
    else:
        raise ValueError ('Cannot create desired vector')
    return OutputValue
    


def doDeletion(input_string, del_type, del_parameters):
    """
    deletions: <del_type> | <del_parameters>:
        from_last | len >0
        from_first | len >0
        Nbp | [N size {1..seqlen}, starting position 0-based >0 && <len(input_string)]
    """
    instring_list = list(input_string)
    output_string = ""
    if del_type == "from_last":
        output_string = ''.join( instring_list[:len(instring_list)-del_parameters] )
    elif del_type == "from_first":
        output_string = ''.join( instring_list[del_parameters:] )
    elif del_type == "Nbp":
        output_string = ''.join( instring_list[:del_parameters[1]] + instring_list[del_parameters[1]+del_parameters[0]:] )
    else:
        print "[AP]\tError, deletion type not defined."
        sys.exit()
    if ''.join(instring_list) == output_string:
        print "it's a bad story, man!"
    return output_string

def doInsertion(input_string, ins_type, ins_parameters, random_seq = True):
    """
    if random(ACGTN) do:
        insertions: <ins_type> | <ins_parameters>:
            from_last | len >0 
            from_first | len >0
            Nbp | [N size {1..seqlen}, starting position 0-based >0 && <len(input_string)]
    if defined dstring do:
        insertions: <ins_type> | <ins_parameters>:
            from_Nbp (not used because so far it is the only one) | [user defined string, starting position 0-based >0 && <len(input_string)]
    """
    nucleotides = list("ACGTN")
    instring_list = list(input_string)
    output_string = ""
    if random_seq:
        if ins_type == "from_last":
            output_string = ''.join(instring_list) + ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters))
        elif ins_type == "from_first":
            output_string = ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters)) + ''.join(instring_list)
        elif ins_type == "Nbp":
            output_string = ''.join( instring_list[:ins_parameters[1]] ) + ''.join(str(random.choice(nucleotides)) for x in range(0,ins_parameters[0])) + ''.join(instring_list[ins_parameters[1]:] )
        else:
            print "[AP]\tError, insertion type not defined."
            sys.exit()
    else: # if not random
        output_string = ''.join( instring_list[:ins_parameters[1]] ) + ins_parameters[0] + ''.join(instring_list[ins_parameters[1]:])
    if ''.join(instring_list) == output_string:
        print "it's a bad story, man!"
    return output_string

def doMutation(input_string, mut_type, mut_parameters):
    """
    (random) mutation excluding same string to replace: <mut_type> | <mut_parameters>:
        Nbp | [start bp, end bp 0-based] where this is a closed interval thus you are mutating from the staring the the ending bp included.

    NB: each base in the interval MUST be different, else you could obtain an unexpected mutation in 2 positions in the same span: from AAAA, span 2 from start, you may obtain CCAA but also ACAA -> the first one is ok but the second one corresponds to the case of SINGLE mutation at base 2.
    """
    nucleotides = list("ACGTN")
    instring_list = list(input_string)
    output_string = ''.join( instring_list[:mut_parameters[0]] )
    if mut_type == "Nbp":
        if mut_parameters[1]-mut_parameters[0]<=0:
            print "[AP]\tWarning: your mutated string has an input interval NON positive! mut_parameters[1] - mut_parameters[0] <= 0::", mut_parameters[1], mut_parameters[0]
        else:
            for bp_index in range(mut_parameters[0],mut_parameters[1]):
                char_tomutate = str(instring_list[bp_index]).capitalize()
                other_nucleotides = [x for x in nucleotides if x != char_tomutate]
                random_string = str(random.choice(other_nucleotides))
                output_string += random_string
            output_string += ''.join(instring_list[mut_parameters[1]:])
    else:
        print "[AP]\tError, mutation type not defined."
        sys.exit()
    if ''.join(instring_list) == output_string:
        print "it's a bad story, man!"
    return output_string


def simulate_seq (Putative_unique_solution_object, LTR_LC_dictionary_plus, LTR_LC_dictionary_minus, out_q):
    
    # This should fix the UNIX issue about random numbers in parallel processes:
    # "what happens is that on Unix every worker process inherits the same state of the random number generator from the parent process. This is why they generate identical pseudo-random sequences."
    scipy.random.seed()
    
    # Take perfect_sequence_dict, perfect_sequence_strandness_dict
    perfect_sequence_dict = Putative_unique_solution_object.perfect_sequence_dict
    perfect_sequence_strandness_dict = Putative_unique_solution_object.perfect_sequence_strandness_dict
    # Prepare random strand-wise headers for LTR attachment
    shuffled_plus_headers = LTR_LC_dictionary_plus.keys()
    random.shuffle(shuffled_plus_headers)
    shuffled_minus_headers = LTR_LC_dictionary_minus.keys()
    random.shuffle(shuffled_minus_headers)
    
    # Prepare result collector - simulated_sequence_dict - {'header': random-picked LTR (strand-wise) + perfect sequence with MID events random distributed}
    simulated_sequence_dict = {}
    
    # Counter for paired lists
    i = 0
    for IS in Putative_unique_solution_object.IS_list:
        headers = IS.reads_key_list
        # Prepare seq_MID_list for the current IS, the exploded version of seq_MID_dict_list for current IS (es. [M,M,D,M,I,I,D,M,...])
        seq_MID_dict = Putative_unique_solution_object.seq_MID_dict_list[i]
        seq_MID_list = []
        for variant, number in seq_MID_dict.items():
            seq_MID_list += [variant]*number
        random.shuffle(seq_MID_list)
        # Prepare list_of_MID_distr
        # len(list_of_MID_distr) = n of sequence (or headers!) for the current IS
        # list_of_MID_distr[j] = n of MID events to pick from seq_MID_list for the j-th sequence
        n_of_seq = len(headers)
        n_of_MID_events = len(seq_MID_list)
        list_of_MID_distr = RandIntVec(n_of_seq, n_of_MID_events, Distribution=RandFloats(n_of_seq))
        j = 0
        
        for header in headers:
            # Simulate sequence
            simulated_sequence = perfect_sequence_dict[header]
            n_MID_to_do = list_of_MID_distr[j]
            MID_to_do = []  # list of 'M', 'D' or 'I'
            for t in range(0, n_MID_to_do):
                MID_to_do.append(seq_MID_list.pop())
            MID_to_do.sort(key=lambda x: ['D', 'M', 'I'].index(x))
            mut_index_control = []
            for MID in MID_to_do:
                seq_len = len(simulated_sequence)
                if MID == 'D':
                    if seq_len > 2:  # prevent sequence from vanishing
                        index = random.randint(1, seq_len-2) # No Del at first or last bp (equivalent to un-read nucleotide)
                        simulated_sequence = doDeletion(simulated_sequence, 'Nbp', (1, index))
                elif MID == 'M':
                    index = random.randint(0, seq_len-1)
                    if len(mut_index_control) < seq_len:  # Control: can't mutate again nucleotides already mutated (unless you can't help it)
                        while index in mut_index_control:
                            index = random.randint(0, seq_len-1)
                    interval = [index, index+1]
                    simulated_sequence = doMutation(simulated_sequence, 'Nbp', interval)
                    mut_index_control.append(index)
                elif MID == 'I':
                    index = random.randint(1, seq_len-2)
                    interval = [index, index+1] # No Ins at first or last bp (equivalent to Mut!!)
                    simulated_sequence = doInsertion(simulated_sequence, 'Nbp', interval, random_seq = True)
            # Attach LTR
            LTR_random_sequence = None
            if ((perfect_sequence_strandness_dict[header] == '+') or (perfect_sequence_strandness_dict[header] == '1')):
                random_header_plus = shuffled_plus_headers.pop()
                sub_dict = LTR_LC_dictionary_plus[random_header_plus]
                LTR_random_sequence = sub_dict['LTR_sequence']
            elif ((perfect_sequence_strandness_dict[header] == '-') or (perfect_sequence_strandness_dict[header] == '2')):
                random_header_minus = shuffled_minus_headers.pop()
                sub_dict = LTR_LC_dictionary_minus[random_header_minus]
                LTR_random_sequence = sub_dict['LTR_sequence']
            simulated_sequence = LTR_random_sequence + simulated_sequence            
            # Collect final simulated sequence
            simulated_sequence_dict[header] = simulated_sequence
            # Next sequence (header)
            j += 1
            
        # Next IS in IS_list
        i += 1
        
    # Put simulated_sequence_dict (simulation outcomes) in out_q
    out_q.put(simulated_sequence_dict)


def parallelized_simulations (Putative_unique_solution_object, LTR_LC_dictionary_plus, LTR_LC_dictionary_minus, nprocs):
    
    # Queue: here simulate_seq will put simulation results (simulated_sequence_dict)
    out_q = multiprocessing.Queue()
    # Process collector
    procs = []
    
    # Start all the nprocs simulate_seq processes: each process puts result (simulated_sequence_dict) in out_q
    for i in range(nprocs):
        p = multiprocessing.Process(
                target=simulate_seq,
                args=(Putative_unique_solution_object,
                      LTR_LC_dictionary_plus,
                      LTR_LC_dictionary_minus,
                      out_q))
        procs.append(p)
        p.start()
    
    # Join all results (simulated_sequence_dict(s)) into a single list (simulated_sequence_dict_list), taking them from out_q
    simulated_sequence_dict_list = []
    for i in range(nprocs):
        simulated_sequence_dict_list.append(out_q.get())
    
    # Wait for all simulate_seq processes to finish    
    for p in procs:
        p.join()
        
    # Put simulated_sequence_dict_list as Putative_unique_solution_object.simulated_sequence_dict_list attribute
    Putative_unique_solution_object.simulated_sequence_dict_list += simulated_sequence_dict_list


##########################################################################################################
### SIMULATED DATA RETRIEVAL #############################################################################
##########################################################################################################

def simulated_data_retrieval (sim_conn_dict, bp_rule, strand_specific_choice):
    
    ### NB: it returns 'None' instead of list_of_Covered_bases_ensambles, if table is EMPTY. ###
    
    #Initialize output data dictionary
    lam_data_dictionay = None
    reads_data_dictionary = None
    
    # Open DB connection
    connection = DB_connection.dbOpenConnection (sim_conn_dict['host'], sim_conn_dict['user'], sim_conn_dict['passwd'], sim_conn_dict['port'], sim_conn_dict['db']) # init connection to DB for importing data
    # Check table row count
    n_row = DB_connection.getTableRowCount(connection, sim_conn_dict['db_table'])
    if n_row < 1:
        return None
    #reads_data_dictionary
    reads_data_dictionary = DB_connection.import_reads_data_from_DB(connection, sim_conn_dict['db_table'], sim_conn_dict['query_step'], sim_conn_dict['reference_genome'])
    #lam_data_dictionay
    lam_data_dictionay  = DB_connection.import_lam_data_from_DB(connection, sim_conn_dict['db_table'], sim_conn_dict['query_step'], sim_conn_dict['reference_genome'])
    # close connection to DB
    DB_connection.dbCloseConnection(connection)
    
    #ordering keys
    reads_data_dictionary_list = reads_data_dictionary.items() #From reads_data_dictionary to a list of kind [(key1,(value1, value2,...)), (key2,(value1, value2,...)), ...]
    reads_data_dictionary_tuple_list=[]
    reads_data_dictionary_tuple_list[:] = [(reads_data_dictionary_list[i][0],) + reads_data_dictionary_list[i][1] for i in range(len(reads_data_dictionary_list))] #From reads_data_dictionary_list to a list of kind [(key1, value1, value2,...), (key2, value1, value2,...), ...]    
    del reads_data_dictionary_list #now useless, substituted by reads_data_dictionary_tuple_list
    reads_data_dictionary_tuple_list_ordered = sorted(reads_data_dictionary_tuple_list, key=itemgetter(2,4,3)) #reads_data_dictionary_tuple_list_ordered is a list of tuple like reads_data_dictionary_tuple_list but MULTIPLE-SORTED by chromosome first (second element of tuple), integration_locus (fourth element of tuple) second and then STRAND.      
    ordered_keys_for_reads_data_dictionary=[]
    ordered_keys_for_reads_data_dictionary[:] = [reads_data_dictionary_tuple_list_ordered[i][0] for i in range(len(reads_data_dictionary_tuple_list_ordered))] #ordered_keys_for_reads_data_dictionary is an ORDERED-LIST-OF-KEY (by chromosome first, then integration_locus) for reads_data_dictionary. "ORDERED" means "STRING ORDERING" (1 is followed by 11, then 2)
    del reads_data_dictionary_tuple_list_ordered
    
    #CB ---> list_of_Covered_Bases = []
    parameters_list = sim_conn_dict['parameters_list']
    seqTracker = False
    raw_read_dictionary = None
    final_read_dictionary = None
    list_of_Covered_Bases = []
    #First read (retrieved by means of 'ordered_keys_for_reads_data_dictionary[0]') is used to create first Covered_base object, then appended into list_of_Covered_Bases
    list_of_Covered_Bases.append(Classes_for_Integration_Analysis.Covered_base(ordered_keys_for_reads_data_dictionary[0], reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific=strand_specific_choice))
    i=0
    for key in ordered_keys_for_reads_data_dictionary[1:]:
        condition = list_of_Covered_Bases[i].add(key, reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific=strand_specific_choice)
        if (condition == -1):
            Classes_for_Integration_Analysis.Covered_base.collapse(list_of_Covered_Bases[i], seqTracker, raw_read_dictionary, final_read_dictionary) #there, list_of_Covered_Bases[i] is completed, so it has to be 'collapsed' to update and freeze its attributes
            list_of_Covered_Bases.append(Classes_for_Integration_Analysis.Covered_base(key, reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific=strand_specific_choice))
            i+=1
    if (type(list_of_Covered_Bases[-1].selective_reads_count) is not dict):
        Classes_for_Integration_Analysis.Covered_base.collapse(list_of_Covered_Bases[-1], seqTracker, raw_read_dictionary, final_read_dictionary)
        
    # CBE
    list_of_Covered_bases_ensambles = []
    #Get strand types and put in strand_list (strands should be indicated in different ways: +/-, 0/1, 1/2... so this is the only way)
    strand_list = []
    strand_list.append(list_of_Covered_Bases[0].strand)
    for covered_base in list_of_Covered_Bases:
        if (covered_base.strand != strand_list[0]):
            strand_list.append(covered_base.strand)
            break
    # Splitting list_of_Covered_Bases in two, strand-wise
    dic_of_list_of_Covered_Bases = {} #{strand_kind:[ordered list of covered base where strand = strand_kind]}    
    for strand_kind in strand_list:
        dic_of_list_of_Covered_Bases.update({strand_kind:[]})
        for covered_base in list_of_Covered_Bases:
            if (covered_base.strand == strand_kind):
                dic_of_list_of_Covered_Bases[strand_kind].append(covered_base)
    #Results temporally appended here, then ordered and put in list_of_Covered_bases_ensambles
    list_of_Covered_bases_ensambles_temp = []
    #for each (both) strands
    for current_strand in strand_list:
        #List of strand-specific results
        list_of_Covered_bases_ensambles_current_strand = []
        #Creating first covered_bases_ensemble with first_covered_base     
        current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(dic_of_list_of_Covered_Bases[current_strand][0])
        for covered_base in dic_of_list_of_Covered_Bases[current_strand][1:]:
            dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
            if ((dist == "undef") or (dist > bp_rule)):
                list_of_Covered_bases_ensambles_current_strand.append(current_covered_bases_ensemble)
                current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base)
            else:
                current_covered_bases_ensemble.push_in(covered_base)
        #APPEND LAST ENSEMBLE            
        list_of_Covered_bases_ensambles_current_strand.append(current_covered_bases_ensemble)
        #APPEND RESULTS FOR THIS STRAND
        list_of_Covered_bases_ensambles_temp = list_of_Covered_bases_ensambles_temp + list_of_Covered_bases_ensambles_current_strand           
        del list_of_Covered_bases_ensambles_current_strand
    # FOR LOOP OVER STRANDS IS OVER, NOW COVERED BASES ENSEMBLES ARE IN AN UN-ORDERED LIST: list_of_Covered_bases_ensambles_temp
    # Ordering list_of_Covered_bases_ensambles_temp by chr then locus then strand and put results in list_of_Covered_bases_ensambles
    list_of_Covered_bases_ensambles = sorted(list_of_Covered_bases_ensambles_temp, key=attrgetter('chromosome', 'starting_base_locus', 'strand'))
    
    return list_of_Covered_bases_ensambles ### THIS IS THE LIST OF CBE RETRIEVED FROM ONE SIMULATION (THE ONE RELATED TO sim_conn_dict)
    
            
















        