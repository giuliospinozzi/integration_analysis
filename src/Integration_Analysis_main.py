#!/usr/bin/python

###Header################################################
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
  
 Note for users:
  - If you have spaces in arguments, please use 
    DOUBLE-quoting 

 Note for Devs:
  - Presently, import_data_from_DB function in 
    DB_connection module has been modified to resolve
    an error due to span=NULL in thalassemia dataset
    (line 73, [...] 100 as `span` [...]
  - Some little temporary changes to work on win8
    (search for tmpfile var and #temporary mode to work on win8
    comments)

 Steps
  - [...]

 Requirements:
  - MySQL client installed in the local machine (where this program is called) and globally callable.
  - MySQL port IS NOT specified! -> if it changes from default you must add this variable as option.
  
-------------------------------------------------------- 
""" 

description = "This application will create detailed matrix of integration sites per LAM or tissue/sample/timepoint"

usage_example = """
Examples of usage:
    APP (--host 127.0.0.1) (--user root) (--pw "") --dbDataset "sequence_mld01.redundant_MLD01_FREEZE_18m_separatedCFC,sequence_thalassemia.pool1_tmp" (--reference_genome hg19) (--query_steps 1000000) --columns sample,tissue,treatment (--columnsToGroup sample) (--IS_method classic) (--bushman_bp_rule 3) (--strand_specific) (--collision) (--rowthreshold 100000)
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
from time import gmtime, strftime
import argparse
import os #temporary mode to work on win8
###############################

###Import Module(s)#########################################################
import DB_connection
import Classes_for_Integration_Analysis
import Matrix_creation
import Common_Functions #already called by Classes_for_Integration_Analysis
import Integration_Sites_retrieving_methods
import DB_filedumpparser
import Collision
############################################################################

###Parsing Arguments############################################################################################################################################################
parser = argparse.ArgumentParser(usage = usage_example, epilog = "[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n", description = description)
parser.add_argument('--host', dest="host", help="IP address to establish a connection with the server that hosts DB. Default is '172.25.39.2' - Alien", action="store", default='172.25.39.2', required=False)
parser.add_argument('--user', dest="user", help="Username to log into the server you just chosen through --host argument. Default is a generic read-only user for Alien", action="store", default='readonly', required=False)
parser.add_argument('--pw', dest="pw", help="Password for the user you choose to log through. Default is the password for the generic read-only user for Alien", action="store", default='readonlypswd', required=False)
parser.add_argument('--dbDataset', dest="dbDataset", help='''Here you have to indicate which database(s) you want to query to retrieve dataset(s). The synatx is, e.g. : "dbschema.dbtable" for one only, "dbschema1.dbtable1,dbschema2.dbtable2,dbschema3.dbtable3" for three. Double quote are generally optional, unless you have spaces or key-characters in names. No default option.''', action="store", required=True)
parser.add_argument('--reference_genome', dest="reference_genome", help="Specify reference genome. Default is 'hg19'", action="store", default="hg19", required=False)
parser.add_argument('--query_steps', dest="query_steps", help="Number of row simultaneously retrieved by a single query. Keep this number low in case of memory leak. If you are going to require --collision, choose thinking to the largest DB you are about to call. Default option is one million row a time", action="store", default = 1000000, required=False)
parser.add_argument('--columns', dest="columns", help="The columns in the final matrix in output. No default option. Available fields: {n_LAM, tag, pool, tissue, sample, treatment, group_name, enzyme}. Example: sample,tissue,treatment.", action="store", required=True)
parser.add_argument('--columnsToGroup', dest="columnsToGroup", help="Among categories given as --columns argument, indicate here with the same syntax the ones you want to merge over, if you desire additional merged columns in output.", action="store", default = None, required=False)
parser.add_argument('--IS_method', dest="IS_method", help="Specify which method run to retrieve Integration Sites. Default option is 'classic'.", action="store", default='classic', required=False)
parser.add_argument('--bushman_bp_rule', dest="bushman_bp_rule", help="If you chose 'classic' method to retrieve IS, here you can set bp number which separate two independent reads cluster. Default option is '3'", action="store", default=3, required=False)
parser.add_argument('--strand_specific', dest="strand_specific", help="If enabled, strands will be treated separately instead of merged together", action="store_true", default=False, required=False)
parser.add_argument('--collision', dest="collision", help="For each dataset given in input to --dbDataset, perform collisions with all the others", action="store_true", default=False, required=False)
parser.add_argument('--rowthreshold', dest="rowthreshold", help="Maximum number of rows allowed to use DB connection. Otherwise, the program will use file dump. Default = 10 millions", action="store", default=10000000, type=int)

args = parser.parse_args()
#################################################################################################################################################################################





##########MAIN##################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################

def main():
    """
    Main part of the program.
    """
        
    ###Input Parameters for DB_connection#################
    #Requested by DB_connection.import_data_from_DB###
    host = args.host    #"172.25.39.2" #Alien        #"127.0.0.1" XAMPP for Devolopment
    user = args.user    #"readonly" #Alien, generic user      #"root" XAMPP for Devolopment
    passwd = args.pw    #'readonlypswd' #Alien        #'' XAMPP for Devolopment
    #db = ##given in sentinel loop##  #such as "sequence_mld01"
    #db_table = ##given in sentinel loop##  #such as "`redundant_mld01_freeze_18m_separatedcfc`"
    query_for_columns=Common_Functions.prepareSELECT(args.columns)   #such as "`sample`,`tissue`,`treatment`"
    reference_genome = args.reference_genome
    query_step = long(args.query_steps)
    ######################################################
    
    #Output file name####################
    file_output_name = db + "_" + db_table + ".tsv"
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
        
#     ###Retrieving Reads Data from DB#################################################################################################
#     #reads_data_dictionary: ["Read Header" => ("reference_genome", "chr", "strand", integration_locusL, read_endL, spanL, "lam_id")]
#     #lam_data_dictionay: ["lam_id" => ("n_LAM", "tag", "pool", "tissue", "sample", "treatment", "group_name", "enzyme")]
#     print "\n{0}\tRetrieving data from DB...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
#     reads_data_dictionary, lam_data_dictionay  = DB_connection.import_data_from_DB(host, user, passwd, db, db_table, query_step, reference_genome)
#     print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
#     ##################################################################################################################################
    
    # check table rows. If table rows > threshold, then use file dump and not DB access
    connection = DB_connection.dbOpenConnection (host, user, passwd, db, db_table)
    # init output data dictionary
    lam_data_dictionay = None
    reads_data_dictionary = None
    if DB_connection.getTableRowCount (connection, db_table) < args.rowthreshold:
        print "\n{0}\tRetrieving data from DB...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        reads_data_dictionary, lam_data_dictionay  = DB_connection.import_data_from_DB(host, user, passwd, db, db_table, query_step, reference_genome)
        print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    else:
        print "\n{0}\tRetrieving data from DB, converting into file and parsing data...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        # reads query and dictionary
        query_select_statement_reads = "'hg19' as reference_genome, header, chr, strand, integration_locus, integration_locus + 100 as integration_locus_end, 100 as span, complete_name as lam_id" #%(reference_genome)
        #tmpfile = DB_filedumpparser.dbTableDump(host, user, passwd, db, db_table, "/dev/shm", query_select_statement_reads)
        tmpfile = DB_filedumpparser.dbTableDump(host, user, passwd, db, db_table, os.getcwd(), query_select_statement_reads) #temporary mode to work on win8
        array_field_reads = ['reference_genome', 'chr', 'strand', 'integration_locus', 'integration_locus_end', 'span', 'lam_id']
        reads_data_dictionary = DB_filedumpparser.parseCSVdumpFile (tmpfile, "header", array_field_reads)
        os.remove(tmpfile)
        # lam query and dictionary
        query_select_statement_lam = "DISTINCT complete_name as lam_id, n_LAM, tag, pool, tissue, sample, treatment, group_name, enzyme"
        #tmpfile = DB_filedumpparser.dbTableDump(host, user, passwd, db, db_table, "/dev/shm", query_select_statement_lam)
        tmpfile = DB_filedumpparser.dbTableDump(host, user, passwd, db, db_table, os.getcwd(), query_select_statement_lam) #temporary mode to work on win8
        array_field_lam = ['n_LAM', 'tag', 'pool', 'tissue', 'sample', 'treatment', 'group_name', 'enzyme']
        lam_data_dictionay = DB_filedumpparser.parseCSVdumpFile (tmpfile, "lam_id", array_field_lam)
        os.remove(tmpfile)
        print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    connection.close()
    
    ###Creating ordered_keys_for_reads_data_dictionary####################################################################################################################
    ###ordered_keys_for_reads_data_dictionary is a list of keys for reads_data_dictionary, useful for retrieving reads ordered by chromosome and then integration_locus###
    ###'ORDER' IS THE STRING'S ONE, so alphabetical and not 'along genome'###
    print "\n{0}\tOrdering retrieved data...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    reads_data_dictionary_list = reads_data_dictionary.items() #From reads_data_dictionary to a list of kind [(key1,(value1, value2,...)), (key2,(value1, value2,...)), ...]
    reads_data_dictionary_tuple_list=[]
    reads_data_dictionary_tuple_list[:] = [(reads_data_dictionary_list[i][0],) + reads_data_dictionary_list[i][1] for i in range(len(reads_data_dictionary_list))] #From reads_data_dictionary_list to a list of kind [(key1, value1, value2,...), (key2, value1, value2,...), ...]    
    del reads_data_dictionary_list #now useless, substituted by reads_data_dictionary_tuple_list              
    reads_data_dictionary_tuple_list_ordered = sorted(reads_data_dictionary_tuple_list, key=itemgetter(2,4)) #reads_data_dictionary_tuple_list_ordered is a list of tuple like reads_data_dictionary_tuple_list but MULTIPLE-SORTED by chromosome first (second element of tuple) and then integration_locus (fourth element of tuple) 
    ordered_keys_for_reads_data_dictionary=[]
    ordered_keys_for_reads_data_dictionary[:] = [reads_data_dictionary_tuple_list_ordered[i][0] for i in range(len(reads_data_dictionary_tuple_list_ordered))] #ordered_keys_for_reads_data_dictionary is an ORDERED-LIST-OF-KEY (by chromosome first, then integration_locus) for reads_data_dictionary. "ORDERED" means "STRING ORDERING" (1 is followed by 11, then 2)
    del reads_data_dictionary_tuple_list_ordered #now useless, substituted by ordered_keys_for_reads_data_dictionary
    print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    ########################################################################################################################################################################
    
    
      
    ###Creating list of 'Covered Bases' objects ############################################################################################################################
    
    print "\n{0}\tArranging data structure...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Initialize list_of_Covered_Bases
    list_of_Covered_Bases = []
    
    #Retrieving parameters list from query_for_columns string
    parameters_list = args.columns.split(",")

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
    
    print "{0}\tCovered bases built!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
       
    ########################################################################################################################################################################
    
    
    
    #Preliminary step to elaborate data and generate output according to user's will###################################################################################################################   
    
    print "\n{0}\tCreating data schema according to user request...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Retrieving labels for matrix columns used for computing data, labels as user wishes and a dictionary to relate them (dict details in DB_connection.get_extra_columns_from_DB)
    #column_labels and merged_column_labels are used also below
    column_labels, merged_column_labels, user_label_dictionary = DB_connection.get_extra_columns_from_DB(host, user, passwd, db, db_table, parameters_list, query_for_columns, reference_genome)
    
    #Declare dictionary of user merged labels
    user_merged_labels_dictionary ={}
    list_of_user_merged_labels =[]
    
    if (args.columnsToGroup != None):
        
        ###Instruction to create merged labels
        user_request = parameters_list # for example ["tissue", "sample", "treatment"] if user give 'tissue,sample,treatment' as --columns argument
        position_to_skip_in_merged_labels = []
        category_to_merge = args.columnsToGroup    
        category_to_merge = category_to_merge.split(",")
        #category_to_merge[:] = [category.replace("'","") for category in category_to_merge]
        position = 0
        for category in user_request:
            if (category in category_to_merge):
                position_to_skip_in_merged_labels.append(position)
            position+=1
            
        
        #Create dictionary of user merged labels: key = user merged label; item = list of user related labels
        for key in column_labels:
            user_label_as_tuple = user_label_dictionary[key][1]
            i=-1
            user_merged_label = ""
            for category in user_label_as_tuple:
                i+=1
                if (i in position_to_skip_in_merged_labels):
                    continue
                else:
                    user_merged_label = user_merged_label + "_{0}".format(category)
            if (user_merged_labels_dictionary.has_key(user_merged_label)):
                #user_merged_labels_dictionary[user_merged_label].append(key)
                user_merged_labels_dictionary[user_merged_label].append(user_label_dictionary[key][0])
            else:
                list_of_user_merged_labels.append(user_merged_label)
                user_merged_labels_dictionary[user_merged_label] = [user_label_dictionary[key][0]]
    
    #Now we have: 
    #user_label_dictionary (key in column_labels)
    #user_merged_labels_dictionary (key in list_of_user_merged_labels) #maybe void
    
    print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    ###################################################################################################################################################################################################
    
    
    
     
    #Redundant Reads Matrix Creation ################################################################################################################
    
    print "\n{0}\tProcessing redundant Reads...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Create redundant reads matrix as list and prepare name for output
    redundant_matrix_file_name, redundant_matrix_as_line_list = Matrix_creation.matrix_output(list_of_Covered_Bases, column_labels, merged_column_labels, file_output_name, strand_specific = strand_specific_choice)
    
    #Convert matrix according to user's requests
    redundant_matrix_as_line_list = Common_Functions.convert_matrix(redundant_matrix_as_line_list, user_label_dictionary, user_merged_labels_dictionary)
    
    #Create output
    file_output = open(redundant_matrix_file_name, 'w')
    for line in redundant_matrix_as_line_list:
        file_output.write(line)

    #Close output file    
    file_output.close()
        
    #Tell user this task finished
    print "{0}\t*Redundant Reads Matrix Created --> {1}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())),redundant_matrix_file_name)
    #################################################################################################################################################
    
    
    
    
    #Grouping Covered Bases in ENSEMBLES, ALL-label############################################################################################################################
    
    print "\n{0}\tGrouping Covered Bases in Ensembles...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
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
        
    print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        
    ###########################################################################################################################################################################        
    



    #Integration Sites Retrieving##############################################################################################################################################        
    
    print "\n{0}\tComputing Integration Sites over Covered Bases Ensembles...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
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
    
    print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    ###########################################################################################################################################################################        
    
    
    
   
    #IS matrix creation ###############################################################################################
    
    print "\n{0}\tProcessing Integration Sites...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Create IS matrix as list and prepare output file name
    IS_matrix_file_name, IS_matrix_as_line_list = Matrix_creation.IS_matrix_output(IS_list, column_labels, merged_column_labels, file_output_name, IS_method, strand_specific=strand_specific_choice)
    
    #Convert matrix according to user's requests
    IS_matrix_as_line_list = Common_Functions.convert_matrix(IS_matrix_as_line_list, user_label_dictionary, user_merged_labels_dictionary)
    
    #Create output, only if needed now (no collision. For one only dataset in input, no-collision is mandatory and superimposed)
    if (args.collision == False):
        file_output = open(IS_matrix_file_name, 'w')
        for line in IS_matrix_as_line_list:
            file_output.write(line)
        #Close output file    
        file_output.close()
        #Tell user this task finished
        print "{0}\t*IS Matrix Created --> {1}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())),IS_matrix_file_name)
        
    else:
        print "{0}\t*IS Matrix computed --> Output file will be generated when all datasets will have been processed".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    ####################################################################################################################
    
    #Return IS_matrix_file_name, IS_matrix_as_line_list #####
    return IS_matrix_file_name, IS_matrix_as_line_list
    #########################################################
    
################################################################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################################################################



### SENTINEL, STARTING CONTROLS, LOOP FOR MULTIPLE DATASET ANALYSES AND COLLISION MANAGMENT 
if __name__ == "__main__":
    """
    INPUT
    OUTPUT
    LOGICS
        1. check variables: columns, db args, ...
    """
    
    #Control to verify user's requests make sense
    check = True ## internal variable for checking variables/controls
    reason = " unexpected error. Try to check syntax and DB connection availability."
    
    ### check if the argument "column" is in the list of the argument "columnToGroup" (sub-set)
    #--columnsToGroup argument
    if (args.columnsToGroup != None):
        selected_category = args.columns.split(",")
        merge_over = args.columnsToGroup.split(",")
        for word in merge_over:
            if (word not in selected_category):
                check = False
                reason =  "each category given as --columnsToGroup argument must be given as --columns argument too. Your input: Columns-> {0}; ColumnsToGroup-> {1}".format(selected_category,merge_over)
                break
   
    #--dbDataset and --collision: controls and management
    dbDataset_tuple_list = [] # [('dbtable1', 'dbschema1'), ('dbtable2', 'dbschema2'), ...]
    dbDataset_split = args.dbDataset.split(",")
    if (("" or "'" or '''"''') in dbDataset_split):
        check = False
        reason = "check syntax in --dbDataset argument"
    if ((args.collision == True) and (len(dbDataset_split)<2)):
        check = False
        reason = "can't perform collision with only one input dataset (see --dbDataset argument)"
    if (check == True):
        for db_string in dbDataset_split:
            db_tupla = None
            db_split = db_string.split(".")
            if (len(db_split)!=2):
                check = False
                reason = "check syntax in --dbDataset argument"
                break
            db_tupla = (db_split[0],db_split[1])
            dbDataset_tuple_list.append(db_tupla)
    
    #Here you can put further controls: when control fails just let check=False and reason="explain the reason why"
    
    #CHECK AND MAIN CALLS    
    if (check == True):
        list_of_IS_results_tuple = [] # [(IS_matrix_file_name1, IS_matrix_as_line_list1),(IS_matrix_file_name2, IS_matrix_as_line_list2), ...]
                                    #not empty only if user perform collision. If collision are not requested, output for IS is created in MAIN
                                    # list of tuple from main function
        print "\n{0}\t***Start***".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        loop_to_do = len(dbDataset_tuple_list)
        i=1
        for db_tupla in dbDataset_tuple_list: # for each tuple in the list -> schema.table
            db = db_tupla[0]
            db_table = db_tupla[1]
            print "\n\n\n{0}\t[START]\tTask {1} of {2}: {3} - {4}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())), i, loop_to_do, db, db_table)
            if (args.collision == True): #collect results in list_of_IS_results_tuple to produce output at the end
                IS_matrix_file_name, IS_matrix_as_line_list = main()
                list_of_IS_results_tuple.append((IS_matrix_file_name, IS_matrix_as_line_list))
                del IS_matrix_file_name, IS_matrix_as_line_list
            else: #nothing needed, main does all task, IS output generation too
                main()
            print "\n{0}\t[SUCCESSFUL COMPLETE]\tTask {1} of {2}: {3} - {4}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())), i, loop_to_do, db, db_table)
            i+=1
            
        #Here you have finished, if collision = False. Else, you have IS results for each dataset in list_of_IS_results_tuple
        
        #COLLISION
        if (args.collision == True):
            
            ###Collision DELTA ####
            #delta = 3
            delta = 4
            #######################
            
            # Print for user
            print "\n\n\n{0}\tCOLLISIONS COMPUTING and IS MATRIX OUTPUT GENERATION".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
            
            # Principal loop over each dataset
            #list_of_IS_COLLIDED_results_tuple = []
            for current_dataset_tuple in list_of_IS_results_tuple:
                print "\n{0}\tComputing data for {1} ... ".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())),current_dataset_tuple[0])
                #list_of_IS_COLLIDED_results_tuple.append(Collision.multiple_collision(current_dataset_tuple, list_of_IS_results_tuple, delta))
                current_dataset_IS_matrix_file_name, current_dataset_IS_matrix_as_line_list_collided = Collision.multiple_collision(current_dataset_tuple, list_of_IS_results_tuple, delta)
                print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
                print "\n{0}\tCreating {1} output file... ".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())),current_dataset_IS_matrix_file_name)
                file_output = open(current_dataset_IS_matrix_file_name, 'w')
                for line in current_dataset_IS_matrix_as_line_list_collided:
                    file_output.write(line)
                file_output.close()
                print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
            
            
            
            # Print for user
            print "\n{0}\tCOLLISIONS COMPLETED".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))    
                            
            ################
            

            
            print "\n\n{0}\t***Tasks Finished***\n\n\tQuit.\n".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        else:
            print "\n\n{0}\t***Tasks Finished***\n\n\tQuit.\n".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    else:
        print "\nYour request can't be processed: {0}".format(reason)
        print "Quit.\n"
        