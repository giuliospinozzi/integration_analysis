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

###Requested Package(s) Import#
import sys
###############################

###Import Module(s)###########################
import Function_for_Gaussian_IS_identification
import DB_connection
import Classes_for_Integration_Analysis
from subprocess import call
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
    
    # LTR_LC_dictionary
    LTR_LC_dictionary = {}
    for header in header_list:
        LTR_sequence, untrimmed_LC_sequence = raw_read_dictionary[header].split(final_read_and_metadata_dictionary[header]['isread_nasequence'])
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
            if prev_char == '^':
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
    
    
    # Return dictionary_for_sequence_simulations, LTR_LC_dictionary
    return dictionary_for_sequence_simulations, LTR_LC_dictionary
    

def get_assembly_path (reference_genome):
    '''
    input: assembly label like 'hg19'
    output: assembly path if available or sys.exit
    '''
    assembly_path = None
    if reference_genome == 'hg19':
        assembly_path = '/opt/genome/human/hg19/index/hg19.fa'
    elif reference_genome == 'mm9':
        assembly_path = '/opt/genome/mouse/mm9/index/mm9.fa'
    else:
        sys.exit("\n\n\t[ERROR] Genome Assembly not available for {}\tQuit.\n\n".format(reference_genome))
    return assembly_path
    

    
def get_seq_from_ref (Putative_unique_solution_object, dictionary_for_sequence_simulations, reference_genome):
    '''
    Putative_unique_solution_object.perfect_sequence_dict = perfect_sequence_dict
    where:
        perfect_sequence_dict = dict with entries {'header': seq}
    '''
    
    # Get IS_list
    IS_list = Putative_unique_solution_object.IS_list
    
    # Prepare BED file for retrieve sequences from reference_genome_assembly
    bed_file_lines = []
    for IS in IS_list:
        chromosome = "chr{}".format(str(IS.chromosome))
        integration_locus = IS.integration_locus
        for header in IS.reads_key_list:
            sub_dict = dictionary_for_sequence_simulations[header]
            line = [chromosome]
            strand = sub_dict['strand']
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
    bed_file_name = "temp.bed"
    with open(bed_file_name, 'w') as bed_file_pointer:
        bed_file_pointer.write('\n'.join(bed_file_lines))
        
    # Call bedtools getfasta
    assembly_path = get_assembly_path (reference_genome)
    fasta_file_name = "temp.fa"
    command = ['fastaFromBed', '-fi', assembly_path, '-bed', bed_file_name, '-s', '-fo', fasta_file_name, '-name']
    exit_status = call(command)
    
    # Read FASTA file --> perfect_sequence_list
    if exit_status != 0:
        sys.exit("\n\n\t[ERROR] Some troubles occured in 'fastaFromBed' call.\tQuit.\n\n")
    
    perfect_sequence_dict = {}
    with open(fasta_file_name, 'r') as source_file:
        
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
    
    
def get_seq_MID_dict_list(Putative_unique_solution_object, dictionary_for_sequence_simulations):
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
        MID_dict = {'Mut': 0, 'Ins': 0, 'Del': 0}
        for header in IS.reads_key_list:
            sub_dict = dictionary_for_sequence_simulations[header]
            if sub_dict['Mut'] != 0:
                MID_dict['Mut'] += sub_dict['Mut']
            if sub_dict['Ins'] != 0:
                MID_dict['Ins'] += sub_dict['Ins']
            if sub_dict['Del'] != 0:
                MID_dict['Del'] += sub_dict['Del']
        seq_MID_dict_list.append(MID_dict)
        
    # Store seq_MID_dict_list in Putative_unique_solution_object
    Putative_unique_solution_object.seq_MID_dict_list = seq_MID_dict_list
    

    
    
            
            
            
















        