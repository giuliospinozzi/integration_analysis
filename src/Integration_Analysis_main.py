#!/usr/bin/python

###Header###############################################
header = """

+------------------------------------------------------+
             ***INTEGRATION ANALYSIS***
             
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 11th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - [...]
  
 Note:
  - [...]

 Steps
  - [...]
  
-------------------------------------------------------- 
""" 
########################################################



########################################################
###BEGIN################################################
########################################################

###Print Header###
print header
##################

###Requested Packages##########
from operator import itemgetter
###############################

###Import Module(s)###
import DB_connection_light as DB_connection
import Classes_for_Integration_Analysis
#import Common_Functions #already called by Classes_for_Integration_Analysis
######################

###Input Parameters for DB_connection#################
host = "127.0.0.1"
user = "root"
passwd = ''
db = "sequence_mld01"
db_table = "`redundant_mld01_freeze_18m_separatedcfc`"
reference_genome = "hg19"
######################################################

###Retrieving Reads Data from DB#################################################################################################
#reads_data_dictionary: ["Read Header" => ("reference_genome", "chr", "strand", integration_locusL, read_endL, spanL, "lam_id")]
#lam_data_dictionay: ["lam_id" => ("n_LAM", "tag", "pool", "tissue", "sample", "treatment", "group_name", "enzyme")]
reads_data_dictionary, lam_data_dictionay  = DB_connection.import_data_from_DB(host, user, passwd, db, db_table, reference_genome)
##################################################################################################################################

###Creating ordered_keys_for_reads_data_dictionary################################################################################
###a list of key for reads_data_dictionary, useful for retrieving reads ordered by chromosome and then integration_locus##########
###'ORDER' IS THE STRING'S ONE, so alphabetical and not 'along genome'############################################################
reads_data_dictionary_list = reads_data_dictionary.items() #From reads_data_dictionary to a list of kind [(key1,(value1, value2,...)), (key2,(value1, value2,...)), ...]
reads_data_dictionary_tuple_list=[]
reads_data_dictionary_tuple_list[:] = [(reads_data_dictionary_list[i][0],) + reads_data_dictionary_list[i][1] for i in range(len(reads_data_dictionary_list))] #From reads_data_dictionary_list to a list of kind [(key1, value1, value2,...), (key2, value1, value2,...), ...]    
del reads_data_dictionary_list #now useless, substituted by reads_data_dictionary_tuple_list              
reads_data_dictionary_tuple_list_ordered = sorted(reads_data_dictionary_tuple_list, key=itemgetter(2,4)) #reads_data_dictionary_tuple_list_ordered is a list of tuple like reads_data_dictionary_tuple_list but MULTIPLE-SORTED by chromosome first (second element of tuple) and then integration_locus (fourth element of tuple) 
ordered_keys_for_reads_data_dictionary=[]
ordered_keys_for_reads_data_dictionary[:] = [reads_data_dictionary_tuple_list_ordered[i][0] for i in range(len(reads_data_dictionary_tuple_list_ordered))] #ordered_keys_for_reads_data_dictionary is an ORDERED-LIST-OF-KEY (by chromosome first, then integration_locus) for reads_data_dictionary. "ORDERED" means "STRING ORDERING" (1 is followed by 11, then 2)
del reads_data_dictionary_tuple_list_ordered #now useless, substituted by ordered_keys_for_reads_data_dictionary
###################################################################################################################################

###Creating list of 'Covered Bases' objects ###########################################

list_of_Covered_Bases = []
list_of_Covered_Bases.append(Classes_for_Integration_Analysis.Covered_base(ordered_keys_for_reads_data_dictionary[0], reads_data_dictionary, lam_data_dictionay))

i=0
for key in ordered_keys_for_reads_data_dictionary[1:]:
    condition = list_of_Covered_Bases[i].add(key, reads_data_dictionary, lam_data_dictionay)
    if (condition == -1):
        Classes_for_Integration_Analysis.Covered_base.collapse(list_of_Covered_Bases[i])
        list_of_Covered_Bases.append(Classes_for_Integration_Analysis.Covered_base(key, reads_data_dictionary, lam_data_dictionay))
        i+=1
        print list_of_Covered_Bases[i-1].chromosome, " ", list_of_Covered_Bases[i-1].strand, " ", list_of_Covered_Bases[i-1].locus, list_of_Covered_Bases[i-1].reads_count #########
        print list_of_Covered_Bases[i-1].selective_reads_count ################
        print "\n" ###########
if (type(list_of_Covered_Bases[-1].selective_reads_count) is not dict):
    Classes_for_Integration_Analysis.Covered_base.collapse(list_of_Covered_Bases[-1])
    print list_of_Covered_Bases[-1].chromosome, " ", list_of_Covered_Bases[-1].strand, " ", list_of_Covered_Bases[-1].locus, list_of_Covered_Bases[-1].reads_count #########
    print list_of_Covered_Bases[-1].selective_reads_count #################
    
#######################################################################################



#===============================================================================
# ###Test Creating list of 'Covered Bases' objects ######################################
# i=0
# print "***********************************"
# print "Print for development:\n"
# print "\nList of Covered Bases"
# for base in list_of_Covered_Bases:
#     print "BASE ", i, ": object ", base
#     print "-> Chromosome: ", base.chromosome
#     print "-> Strand: ", base.strand
#     print "-> Locus: ", base.locus
#     print "-> List of keys: ", base.list_of_reads
#     print "-> Reads Count: ", base.reads_count
#     print "-> Selective Reads Count: ", base.selective_reads_count
#     print "---------------\n"
#     i+=1
#     if (i>10):
#         break
# print "***********************************"
# #######################################################################################
#===============================================================================


#===============================================================================
# ###Test get_lam method############
# i=1
# for key in ordered_keys_for_reads_data_dictionary:
#     print i, ") read header: ", key, "; read attribute: ", reads_data_dictionary[key], "."
#     print "...related lam data: ", Common_Functions.get_lam(key, reads_data_dictionary, lam_data_dictionay), " !"
#     i+=1
#     if (i>100):
#         break
# ##################################    
#===============================================================================













