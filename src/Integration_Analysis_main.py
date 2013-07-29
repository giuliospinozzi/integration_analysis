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

description = "This application will create detailed matrix of integration sites per LAM or tissue/sample/timepoint"

usage_example = """
Examples of usage:
    APP --dbschema sequence_mld01 --dbtable redundant_mld01_freeze_18m_separatedcfc --columnsToGroup 'sample,tissue,treatment' -o matrix_redundant_mld01_freeze_18m_separatedcfc.tsv
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
import argparse
###############################

###Import Module(s)#########################################################
import DB_connection
import Classes_for_Integration_Analysis
import Matrix_creation
import Common_Functions #already called by Classes_for_Integration_Analysis
############################################################################

###Parsing Arguments############################################################################################################################################################
parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--dbschema', dest="dbschema", help="The input databse schema", action="store", required=True)
parser.add_argument('--dbtable', dest="dbtable", help="The table of the db schema. No default option.", action="store", required=True)
parser.add_argument('--columnsToGroup', dest="columnsToGroup", help="The columns to group in the final matrix in output. No default option.", action="store", required=True)
parser.add_argument('-o', '--outfilename', dest="outfilename", help="summary columns start and end with _ chars (tsv). No default option.", action="store", required=True)

args = parser.parse_args()
#################################################################################################################################################################################

##########MAIN##################################################################################################################################################################################################################################################################################################################################


def main():
    """
    Main part of the program.
    """
        
    ###Input Parameters for DB_connection#################
    #Requested by DB_connection.import_data_from_DB###
    host = "127.0.0.1"
    user = "root"
    passwd = ''
    db = args.dbschema #db = "sequence_mld01"
    db_table = args.dbtable #db_table = "`redundant_mld01_freeze_18m_separatedcfc`"
    query_for_columns=Common_Functions.prepareSELECT(args.columnsToGroup)   #"`sample`,`tissue`,`treatment`"
    reference_genome = "hg19"
    ######################################################
    
    #Output file name####################
    file_output_name = args.outfilename
    #####################################
        
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
    
    #Retrieving parameters list from query_for_columns string
    parameters_list = query_for_columns.split(",")
    parameters_list[:] = [parameter.replace("`","") for parameter in parameters_list]
    
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
    Matrix_creation.matrix_output(list_of_Covered_Bases, column_labels, merged_column_labels, file_output_name)
    
    ##################################################################################################################################
    
    
    
    #Grouping Covered Bases in ENSEMBLES#############################################################################################
    
    #Defining maximum distance between Correlated Covered Bases
    bushamn_bp_rule = 6
    
    #Dictionary to collect result {'label1': list_of_Covered_bases_ensambles_for_label1, 'label2': list_of_Covered_bases_ensambles_for_label1, ...}
    selective_Covered_bases_ensambles = {}
    
    #Loop over column_labels
    for label in column_labels:
        
        # Refresh over each loop
        current_list_of_Covered_bases_ensambles = []
        selective_Covered_bases_ensambles.update({label:current_list_of_Covered_bases_ensambles})
        first_covered_base = None
        dist = None
        i = 0
        
        #Catch the first covered_base with non-zero read count for "label" of this loop
        for covered_base in list_of_Covered_Bases:
            i+=1
            if (label in covered_base.selective_reads_count.keys()):
                first_covered_base = covered_base
                break #first covered_base catched or list_of_Covered_Bases finished, so stop and keep first_covered_base and i
            
        if (first_covered_base == None):
            continue #if first_covered_base is still "None", there are no covered_base with non-zero read count for "label" of this loop... skip the loop over this label and start with another one
        
        #Creating first covered_bases_ensemble with first_covered_base     
        current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(first_covered_base, label_selection = label)
        
        #Loop over left covered bases, starting from the i-th
        for covered_base in list_of_Covered_Bases[i:]:
            dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base, label_selection = label) #retrieving distance between last element in current_covered_bases_ensemble.Covered_bases_list and current covered_base (looping) 
            
            if (dist == "undef"): #Two possible reasons
                if (label in covered_base.selective_reads_count.keys()): #current covered_base has non-zero read count for "label" but is too much far (e.g. on a different chromosome): append current_covered_bases_ensemble to current_list_of_Covered_bases_ensambles and create another one with current covered_base
                    current_list_of_Covered_bases_ensambles.append(current_covered_bases_ensemble)
                    current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base, label_selection = label)
                else: # current covered_base has zero read count for "label": just skip to the next covered_base over this loop
                    continue
            
            else: #...at this point covered_base has non-zero read count for "label" and dist has to be a positive number.
                if (dist <= bushamn_bp_rule): #dist is small: push covered_base in current_covered_bases_ensemble
                    current_covered_bases_ensemble.push_in(covered_base, label_selection = label)
                else: #dist is large: append current_covered_bases_ensemble to current_list_of_Covered_bases_ensambles and create another one with current covered_base
                    current_list_of_Covered_bases_ensambles.append(current_covered_bases_ensemble)
                    current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base, label_selection = label)
        
        #Loop over left covered bases has finished, store results updating selective_Covered_bases_ensambles dictionary           
        selective_Covered_bases_ensambles.update({label:current_list_of_Covered_bases_ensambles})
    
    #Print for development
    print "\n*** About Covered_bases_ensambles ***"
    print "List of Labels: ", column_labels
    print "\n\nDictionary retrieved: "
    for key, item in selective_Covered_bases_ensambles:
        print key, item
        print "Details: "
        for element in selective_Covered_bases_ensambles[key]:
            print "label: ", element.label, "; chr: ", element.chromosome, "; srd: ", element.strand, "; start: ", element.starting_base_locus, "; end: ", element.ending_base_locus, "; span: ", element.spanned_bases, "; CB: ", element.n_covered_bases, "; total_r: ", element.n_total_reads
        print "\n"
    
    ##################################################################################################################################
    
    
    #Final print#################################
    print "\n[AP]\tTask Finished, closing.\n"
    #############################################
    
################################################################################################################################################################################################################################################################################################################################################


# sentinel
if __name__ == "__main__":
    main()