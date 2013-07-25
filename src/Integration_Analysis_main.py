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



###Import Module(s)#########################################################
import DB_connection
import Classes_for_Integration_Analysis
import Matrix_creation
#import Common_Functions #already called by Classes_for_Integration_Analysis
############################################################################



###Input Parameters for DB_connection#################
#Requested by DB_connection.import_data_from_DB###
host = "127.0.0.1"
user = "root"
passwd = ''
db = "sequence_mld01"
db_table = "`redundant_mld01_freeze_18m_separatedcfc`"
query_for_columns="`sample`, `tissue`, `treatment`, `n_LAM`"
reference_genome = "hg19"
######################################################

#Retrieving parameters list from query_for_columns string########################
parameters_list = query_for_columns.split(", ")
parameters_list[:] = [parameter.replace("`","") for parameter in parameters_list]
#################################################################################

###Retrieving Reads Data from DB#################################################################################################
#reads_data_dictionary: ["Read Header" => ("reference_genome", "chr", "strand", integration_locusL, read_endL, spanL, "lam_id")]
#lam_data_dictionay: ["lam_id" => ("n_LAM", "tag", "pool", "tissue", "sample", "treatment", "group_name", "enzyme")]
reads_data_dictionary, lam_data_dictionay  = DB_connection.import_data_from_DB(host, user, passwd, db, db_table, reference_genome)
##################################################################################################################################



###Creating ordered_keys_for_reads_data_dictionary####################################################################################################################
###ordered_keys_for_reads_data_dictionary is a list of keys for reads_data_dictionary, useful for retrieving reads ordered by chromosome and then integration_locus###
###'ORDER' IS THE STRING'S ONE, so alphabetical and not 'along genome'###
reads_data_dictionary_list = reads_data_dictionary.items() #From reads_data_dictionary to a list of kind [(key1,(value1, value2,...)), (key2,(value1, value2,...)), ...]
reads_data_dictionary_tuple_list=[]
reads_data_dictionary_tuple_list[:] = [(reads_data_dictionary_list[i][0],) + reads_data_dictionary_list[i][1] for i in range(len(reads_data_dictionary_list))] #From reads_data_dictionary_list to a list of kind [(key1, value1, value2,...), (key2, value1, value2,...), ...]    
del reads_data_dictionary_list #now useless, substituted by reads_data_dictionary_tuple_list              
reads_data_dictionary_tuple_list_ordered = sorted(reads_data_dictionary_tuple_list, key=itemgetter(2,4)) #reads_data_dictionary_tuple_list_ordered is a list of tuple like reads_data_dictionary_tuple_list but MULTIPLE-SORTED by chromosome first (second element of tuple) and then integration_locus (fourth element of tuple) 
ordered_keys_for_reads_data_dictionary=[]
ordered_keys_for_reads_data_dictionary[:] = [reads_data_dictionary_tuple_list_ordered[i][0] for i in range(len(reads_data_dictionary_tuple_list_ordered))] #ordered_keys_for_reads_data_dictionary is an ORDERED-LIST-OF-KEY (by chromosome first, then integration_locus) for reads_data_dictionary. "ORDERED" means "STRING ORDERING" (1 is followed by 11, then 2)
del reads_data_dictionary_tuple_list_ordered #now useless, substituted by ordered_keys_for_reads_data_dictionary
########################################################################################################################################################################



###Creating list of 'Covered Bases' objects ############################################################################################################################

#Initialize list_of_Covered_Bases
list_of_Covered_Bases = []

#First read (retrieved by means of 'ordered_keys_for_reads_data_dictionary[0]') is used to create first Covered_base object, then appended into list_of_Covered_Bases
list_of_Covered_Bases.append(Classes_for_Integration_Analysis.Covered_base(ordered_keys_for_reads_data_dictionary[0], reads_data_dictionary, lam_data_dictionay, parameters_list))

#This cycle creates the whole list_of_Covered_Bases
#It iterates over ORDERED_keys_for_reads_data_dictionary (mandatory)
#For each read it tries to 'add' it (see Classes_for_Integration_Analysis) into current element of list_of_Covered_Bases: if allowed by 'add' controls, it's done. Else, nothing happens ('add' returns -1) and following
#if clause creates a new  Covered_Bases object, suddenly appended to list_of_Covered_Bases. 
i=0
for key in ordered_keys_for_reads_data_dictionary[1:]:
    condition = list_of_Covered_Bases[i].add(key, reads_data_dictionary, lam_data_dictionay, parameters_list)
    if (condition == -1):
        Classes_for_Integration_Analysis.Covered_base.collapse(list_of_Covered_Bases[i]) #there, list_of_Covered_Bases[i] is completed, so it has to be 'collapsed' to update and freeze its attributes
        list_of_Covered_Bases.append(Classes_for_Integration_Analysis.Covered_base(key, reads_data_dictionary, lam_data_dictionay, parameters_list))
        i+=1
        #Print for development
        #print list_of_Covered_Bases[i-1].chromosome, " ", list_of_Covered_Bases[i-1].strand, " ", list_of_Covered_Bases[i-1].locus, list_of_Covered_Bases[i-1].reads_count
        #print list_of_Covered_Bases[i-1].selective_reads_count
        #print "\n"

#In some cases, necessary to 'collapse' the last Covered_Base object in list_of_Covered_Bases
if (type(list_of_Covered_Bases[-1].selective_reads_count) is not dict):
    Classes_for_Integration_Analysis.Covered_base.collapse(list_of_Covered_Bases[-1])
    #Print for development
    #print list_of_Covered_Bases[-1].chromosome, " ", list_of_Covered_Bases[-1].strand, " ", list_of_Covered_Bases[-1].locus, list_of_Covered_Bases[-1].reads_count
    #print list_of_Covered_Bases[-1].selective_reads_count
    
########################################################################################################################################################################



#Matrix Creation ################################################################################################################

#Retrieving labels for matrix columns
column_labels, merged_column_labels = DB_connection.get_extra_columns_from_DB(host, user, passwd, db, db_table, parameters_list, query_for_columns, reference_genome)

#Create and print final matrix in an output file
Matrix_creation.matrix_output(list_of_Covered_Bases, column_labels, merged_column_labels, db_table)

##################################################################################################################################


