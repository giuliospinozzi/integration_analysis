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
    APP --dbschema sequence_mld01 --dbtable redundant_mld01_freeze_18m_separatedcfc (--query_steps 1000000) --columnsToGroup 'sample,tissue,treatment' (--IS_method classic) (--bushman_bp_rule 6) (--strand_specific) -o matrix_redundant_mld01_freeze_18m_separatedcfc.tsv
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
from operator import attrgetter
import argparse
###############################

###Import Module(s)#########################################################
import DB_connection
import Classes_for_Integration_Analysis
import Matrix_creation
import Common_Functions #already called by Classes_for_Integration_Analysis
import Integration_Sites_retrieving_methods
############################################################################

###Parsing Arguments############################################################################################################################################################
parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--dbschema', dest="dbschema", help="The input databse schema", action="store", required=True)
parser.add_argument('--dbtable', dest="dbtable", help="The table of the db schema. No default option.", action="store", required=True)
parser.add_argument('--query_steps', dest="query_steps", help="Number of row simultaneously retrieved by a single query. Keep this number low in case of memory leak. Default option is one million row a time", action="store", default = 1000000, required=False)
parser.add_argument('--columnsToGroup', dest="columnsToGroup", help="The columns to group in the final matrix in output. No default option.", action="store", required=True)
parser.add_argument('--IS_method', dest="IS_method", help="Specify which method run to retrieve Integration Sistes. Default option is 'classic'.", action="store", default='classic', required=False)
parser.add_argument('--bushman_bp_rule', dest="bushman_bp_rule", help="If you chose 'classic' method to retrieve IS, here you can set bp number which separate two independent reads cluster. Default option is '3'", action="store", default=3, required=False)
parser.add_argument('--strand_specific', dest="strand_specific", help="If enabled, strands will be treated separately instead of merged together", action="store_true", default=False, required=False)
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
    db = args.dbschema #such as "sequence_mld01"
    db_table = args.dbtable #such as "`redundant_mld01_freeze_18m_separatedcfc`"
    query_for_columns=Common_Functions.prepareSELECT(args.columnsToGroup)   #such as "`sample`,`tissue`,`treatment`"
    reference_genome = "hg19"
    query_step = long(args.query_steps)
    ######################################################
    
    #Output file name####################
    file_output_name = args.outfilename
    #####################################
        
    #Retrieving IS method choice####################################
    strand_specific_choice = args.strand_specific
    IS_method = args.IS_method
    bushamn_bp_rule = 6
    ################################################################
    
    #Defining maximum distance between Correlated Covered Bases#####
    if (IS_method == "classic"):
        bushamn_bp_rule = int(args.bushman_bp_rule)
    ################################################################
        
    ###Retrieving Reads Data from DB#################################################################################################
    #reads_data_dictionary: ["Read Header" => ("reference_genome", "chr", "strand", integration_locusL, read_endL, spanL, "lam_id")]
    #lam_data_dictionay: ["lam_id" => ("n_LAM", "tag", "pool", "tissue", "sample", "treatment", "group_name", "enzyme")]
    reads_data_dictionary, lam_data_dictionay  = DB_connection.import_data_from_DB(host, user, passwd, db, db_table, query_step, reference_genome)
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
    list_of_Covered_Bases.append(Classes_for_Integration_Analysis.Covered_base(ordered_keys_for_reads_data_dictionary[0], reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific=strand_specific_choice))
    
    #This cycle creates the whole list_of_Covered_Bases
    #It iterates over ORDERED_keys_for_reads_data_dictionary (mandatory)
    #For each read it tries to 'add' it (see Classes_for_Integration_Analysis) into current element of list_of_Covered_Bases: if allowed by 'add' controls, it's done. Else, nothing happens ('add' returns -1) and following
    #if clause creates a new  Covered_Bases object, suddenly appended to list_of_Covered_Bases. 
    i=0
    for key in ordered_keys_for_reads_data_dictionary[1:]:
        condition = list_of_Covered_Bases[i].add(key, reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific=strand_specific_choice)
        if (condition == -1):
            Classes_for_Integration_Analysis.Covered_base.collapse(list_of_Covered_Bases[i]) #there, list_of_Covered_Bases[i] is completed, so it has to be 'collapsed' to update and freeze its attributes
            list_of_Covered_Bases.append(Classes_for_Integration_Analysis.Covered_base(key, reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific=strand_specific_choice))
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
    #column_labels and merged_column_labels are used also below
    column_labels, merged_column_labels = DB_connection.get_extra_columns_from_DB(host, user, passwd, db, db_table, parameters_list, query_for_columns, reference_genome)
    
    #Create and print final matrix in an output file
    Matrix_creation.matrix_output(list_of_Covered_Bases, column_labels, merged_column_labels, file_output_name, strand_specific = strand_specific_choice)
    
    ##################################################################################################################################
    
    
    
    
    #Grouping Covered Bases in ENSEMBLES, ALL-label############################################################################################################################
    
    #List of results: list_of_Covered_bases_ensambles
    all_labels_Covered_bases_ensambles = []
    
    
    #Ensemble grouping
    
    #Creating first covered_bases_ensemble with first_covered_base     
    current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(list_of_Covered_Bases[0], strand_specific = strand_specific_choice)
    
    #If strand_specific_choice == False, algorithm goes straight
    if (strand_specific_choice == False):
    
        #different if user chooses "classic" method to retrieving IS (also bushamn_bp_rule value is automatically changed on top)
        if (IS_method == "classic"):
            for covered_base in list_of_Covered_Bases[1:]:
                dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
                if ((dist == "undef") or (dist > bushamn_bp_rule) or (current_covered_bases_ensemble.spanned_bases == (bushamn_bp_rule + 1)) or ((current_covered_bases_ensemble.spanned_bases + dist)>(bushamn_bp_rule + 1))):
                    all_labels_Covered_bases_ensambles.append(current_covered_bases_ensemble)
                    current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base, strand_specific = strand_specific_choice)
                else:
                    current_covered_bases_ensemble.push_in(covered_base)
                    
        #for the others possible user's choices, groping method goes straight
        else:
            for covered_base in list_of_Covered_Bases[1:]:
                dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
                if ((dist == "undef") or (dist > bushamn_bp_rule)):
                    all_labels_Covered_bases_ensambles.append(current_covered_bases_ensemble)
                    current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base, strand_specific = strand_specific_choice)
                else:
                    current_covered_bases_ensemble.push_in(covered_base)
                
        all_labels_Covered_bases_ensambles.append(current_covered_bases_ensemble) #APPEND LAST ENSEMBLE
        
        # NOW COVERED BASES ENSEMBLES ARE IN AN ORDERED LIST: all_labels_Covered_bases_ensambles
    
                    
    #If strand_specific_choice == True, algorithm has to cycle over list_of_Covered_Bases two times, creating two strand-specific list, then merged (preserving ordering)
    #Such a control prevents from "false ensembles splitting"
    elif (strand_specific_choice == True):
        
        #Get strand types and put in strand_list (strands should be indicated in different ways: +/-, 0/1, 1/2... so this is the only way)
        strand_list = []
        strand_list.append(list_of_Covered_Bases[0].strand)
        for covered_base in list_of_Covered_Bases:
            if (covered_base.strand != strand_list[0]):
                strand_list.append(covered_base.strand)
                break
        
        #Results temporally appended here, then ordered and put in all_labels_Covered_bases_ensambles
        all_labels_Covered_bases_ensambles_temp = []
        
        #for both strands
        check = False #necessary to avoid duplicate append
        for current_strand in strand_list:
            
            #List of strand-specific results
            all_labels_Covered_bases_ensambles_current_strand = []
            
            #different if user chooses "classic" method to retrieving IS (also bushamn_bp_rule value is automatically changed on top)
            if (IS_method == "classic"):
                for covered_base in list_of_Covered_Bases[1:]:
                    check = False
                    #If covered_base's strand doesn't match with current_strand, this covered_base is skipped
                    if (covered_base.strand != current_strand):
                        check = True
                        continue
                    dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
                    if ((dist == "undef") or (dist > bushamn_bp_rule) or (current_covered_bases_ensemble.spanned_bases == (bushamn_bp_rule + 1)) or ((current_covered_bases_ensemble.spanned_bases + dist)>(bushamn_bp_rule + 1))):
                        all_labels_Covered_bases_ensambles_current_strand.append(current_covered_bases_ensemble)
                        current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base)
                    else:
                        current_covered_bases_ensemble.push_in(covered_base)
                        
            #for the others possible user's choices, groping method goes straight
            else:
                for covered_base in list_of_Covered_Bases[1:]:
                    check = False
                    #If covered_base's strand doesn't match with current_strand, this covered_base is skipped
                    if (covered_base.strand != current_strand):
                        check = True
                        continue
                    dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
                    if ((dist == "undef") or (dist > bushamn_bp_rule)):
                        all_labels_Covered_bases_ensambles_current_strand.append(current_covered_bases_ensemble)
                        current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base)
                    else:
                        current_covered_bases_ensemble.push_in(covered_base)
            
            if (check == False):        
                all_labels_Covered_bases_ensambles_current_strand.append(current_covered_bases_ensemble) #APPEND LAST ENSEMBLE
            
            all_labels_Covered_bases_ensambles_temp = all_labels_Covered_bases_ensambles_temp + all_labels_Covered_bases_ensambles_current_strand #APPEND RESULTS FOR THIS STRAND
            
            del all_labels_Covered_bases_ensambles_current_strand
            
        # FOR LOOP OVER STRANDS IS OVER, NOW COVERED BASES ENSEMBLES ARE IN AN UN-ORDERED LIST: all_labels_Covered_bases_ensambles_temp
        
        # Ordering all_labels_Covered_bases_ensambles_temp by chr then locus then strand and put results in all_labels_Covered_bases_ensambles
        all_labels_Covered_bases_ensambles = sorted(all_labels_Covered_bases_ensambles_temp, key=attrgetter('chromosome', 'starting_base_locus', 'strand'))
        
        #=======================================================================
        # #Print for development
        # log_file = open('dev_log_file', 'w')
        # for row in all_labels_Covered_bases_ensambles:
        #     log_file.write("{0}\t{1}\t{2}\n".format(str(row.chromosome), str(row.starting_base_locus), str(row.strand)))
        # log_file.close()
        #=======================================================================
    ###########################################################################################################################################################################        
    



    #Integration Sites Retrieving##############################################################################################################################################        
    
    #Initialize list of results:
    IS_list = []
    
    #Classic method
    if (IS_method == "classic"):
        for Covered_bases_ensamble in all_labels_Covered_bases_ensambles:
            IS_list.append(Integration_Sites_retrieving_methods.classic(Covered_bases_ensamble, strand_specific = strand_specific_choice))
    
    #NOW INTEGRATION SITES RETRIEVED THROUGH "CLASSIC" METHOD ARE IN IS_LIST

    
    #Whatever method    
    if (IS_method == "whatever"):
        ###Here the code, when "whatever" new method will be available
        pass
        
    #NOW INTEGRATION SITES RETRIEVED THROUGH "WHATEVER" METHOD ARE IN IS_LIST
    
    ###########################################################################################################################################################################        
    
    
    
   
    #IS matrix creation##############################################################
    Matrix_creation.IS_matrix_output(IS_list, column_labels, merged_column_labels, file_output_name, IS_method, strand_specific=strand_specific_choice)
    #################################################################################
    
    
    
    #Final print#################################
    print "\n[AP]\tTask Finished, closing.\n"
    #############################################
    
################################################################################################################################################################################################################################################################################################################################################


# sentinel
if __name__ == "__main__":
    main()