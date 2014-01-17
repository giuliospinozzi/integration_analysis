#!/usr/bin/python

###Header##################################################################
header = """

+--------------------------------------------------------+
             ***INTEGRATION ANALYSIS***
             
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  January 8th, 2014
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 2.0
+--------------------------------------------------------+

 Description:
  - Retrieving data from a network DB, this application
    creates detailed matrixes of Redundant Reads 
    and Integration Sites. User can choose to separate
    (--columns) and partially-aggregate (--columnsToGroup)
    results according to different categories
    (e.g. sample, tissue, treatment...) and to perform
    collisions to compare different input datasets.
  
 Note for users:
  - If you have spaces in arguments, please use 
    DOUBLE-quoting. Likewise, use "" in place of an empty
    argument (e.g. in --pw argument, if chosen account
    doesn't have it)

 Note for Devs:
  - Presently, import_data_from_DB function in 
    DB_connection module has been modified to resolve
    an error due to span=NULL in thalassemia dataset
    (line 73, [...] 100 as `span` [...]
  - Some little temporary changes to work on win8
    (search for tmpfile var and #temporary mode to work 
    on win8 comments, in source code)

 Steps
  - [...]

 Requirements:
  - MySQL client installed in the local machine (where 
    this program is called) and globally callable.
  
+--------------------------------------------------------+ 
"""

print header
###########################################################################


###Requested Packages####################
from operator import itemgetter
from operator import attrgetter
from time import gmtime, strftime
import argparse
import os #temporary mode to work on win8
#########################################


###Import Module(s)#########################################################
import Preliminary_controls
import DB_connection
import Classes_for_Integration_Analysis
import Matrix_creation
import Common_Functions #already called by Classes_for_Integration_Analysis
import Integration_Sites_retrieving_methods
import DB_filedumpparser
import Collision
import Function_for_Gaussian_IS_identification
import output_module
############################################################################


###Parsing Arguments############################################################################################################################################################
usage_example = '''\n\nExamples of usage for 'classic' IS-retrieving-method: python Integration_Analysis.py (--host 172.25.39.57) (--user readonly) (--pw readonlypswd) (--port 3306) --dbDataset "sequence_mld01.fu18m,sequence_mld02.fu18m" (--query_steps 1000000) (--rowthreshold 10000000) (--reference_genome hg19) --columns sample,tissue,treatment (--columnsToGroup sample) --IS_method classic (--bushman_bp_rule 3) (--strand_specific) (--collision (--set_radius 3)) (--tsv)
\nExamples of usage for 'gauss' IS-retrieving-method: python Integration_Analysis.py (--host 172.25.39.57) (--user readonly) (--pw readonlypswd) (--port 3306) --dbDataset "sequence_mld01.fu18m,sequence_mld02.fu18m" (--query_steps 1000000) (--rowthreshold 10000000) (--reference_genome hg19) --columns sample,tissue,treatment (--columnsToGroup sample) --IS_method gauss (--strand_specific) --interaction_limit 3 --alpha 0.6 (--collision (--set_radius 3)) (--tsv)
\nRound brackets highlight settings/arguments that are optional or supplied with defaults\n\n\ndescription:\n'''
description = "This application creates detailed matrixes of Redundant Reads and Integration Sites retrieving data from a network DB. User can choose to separate (--columns) and partially-aggregate (--columnsToGroup) results according to different categories (e.g. sample, tissue, treatment...) and to perform collisions to compare different input datasets"


parser = argparse.ArgumentParser(usage = usage_example, epilog = "\n[ hSR-TIGET - Vector Integration Core - Bioinformatics ] \n\n", description = description)

parser.add_argument('--host', dest="host", help="IP address to establish a connection with the server that hosts DB.\nDefault is 'localhost'.\nTips: '172.25.39.2'-Alien, '172.25.39.57'-Gemini", action="store", default='localhost', required=False)
parser.add_argument('--user', dest="user", help="Username to log into the server you just chosen through --host argument.\nDefault is a generic read-only user for Gemini/Alienware", action="store", default='readonly', required=False)
parser.add_argument('--pw', dest="pw", help="Password for the user you choose to log through.\nDefault is the password for the generic read-only user for Gemini/Alienware", action="store", default='readonlypswd', required=False)
parser.add_argument('--port', dest="dbport", help="Database port.\nDefault is 3306", action="store", default=3306, required=False, type=int)
parser.add_argument('--dbDataset', dest="dbDataset", help='''Here you have to indicate which database(s) you want to query to retrieve dataset(s). The synatx is, e.g. : "dbschema.dbtable" for one only, "dbschema1.dbtable1,dbschema2.dbtable2,dbschema3.dbtable3" for three. Double quote are generally optional, unless you have spaces or key-characters in names.\nNo default option. Required.''', action="store", required=True)
parser.add_argument('--query_steps', dest="query_steps", help="Number of row simultaneously retrieved by a single query. Keep this number low in case of memory leak (choose thinking to the largest DB you are about to call, in case of multiple datasets).\nDefault option is one million row a time", action="store", default = 1000000, required=False, type=int)
parser.add_argument('--rowthreshold', dest="rowthreshold", help="Maximum number of rows allowed to use direct DB connection. Otherwise, the program will use file dump (slower but saves a lot of memory).\nDefault = 10 millions", action="store", default=10000000, type=int)
parser.add_argument('--reference_genome', dest="reference_genome", help="Specify reference genome. Default is 'hg19'", action="store", default="hg19", required=False)
parser.add_argument('--columns', dest="columns", help="Indicate the columns for the final matrix output. Available fields: {n_LAM, tag, pool, tissue, sample, treatment, group_name, enzyme}.\nNo default option. Required.\nExample: sample,tissue,treatment.", action="store", required=True)
parser.add_argument('--columnsToGroup', dest="columnsToGroup", help="Among categories given as --columns argument, indicate here the ones you want to merge over (same syntax) if you desire additional merged columns in final output.\nNo default option. Example: sample", action="store", default = None, required=False)
parser.add_argument('--IS_method', dest="IS_method", help="Specify which method run to retrieve Integration Sites: 'classic' or 'gauss'. You'll be able to tune 'classic' through --bushman_bp_rule, if you don't like defaults, while you'll have to properly set-up 'gauss' through --interaction_limit and --alpha (no defaults provided for it).\nNo default option. Required", action="store", default=None, required=True)
parser.add_argument('--bushman_bp_rule', dest="bushman_bp_rule", help="Minimum number of empty base-pairs between reads belonging to different cluster (also called Covered Bases Ensembles). If you chose 'classic' method to retrieve IS, this number also set the maximum number of Covered Bases allowed for a Covered Bases Ensemble and is totally customizable.\nDefault option is '3', i.e. 'minimum 3 empty-bp between indipendent ensembles, maximum 3 covered-bases per ensemble'. Conversely, if you chose 'gauss' method, it will be automatically set as 2*interaction_limit + 1 (overriding your setting) and no limit of dimension will be set for ensembles construction.", action="store", default=3, required=False, type=int)
parser.add_argument('--strand_specific', dest="strand_specific", help="If called, strands will be treated separately instead of merged together", action="store_true", default=False, required=False)
parser.add_argument('--interaction_limit', dest="interaction_limit", help="Only in case of '--IS_method gauss', here you have to set the 'action radius' of a peak (N of bases flanking the peak: so 'int' and '>=1'); this choice will affect --bushman_bp_rule and --delta arguments, overriding defaults / user's choices with optimal settings.\nNo default option. Tip: if you have no idea, try 3 (stringent) or 4 .", action="store", default=None, required=False, type=int)
parser.add_argument('--alpha', dest="alpha", help="Only in case of '--IS_method gauss', here you have to set 'HOW MANY SIGMAS are equal to HALF-BASEPAIR'. This choice should be made wisely, together with --interaction_limit. Some controls will be performed and you'll be warned in case of bad settings.\nNo default option. Tip: if you have no idea, try 0.6 (stringent) or 0.3.", action="store", default=None, required=False)
parser.add_argument('--collision', dest="collision", help="If called, collisions will be performed for each dataset with all the others.\nCollision radius is set by default equal to bushman_bp_rule+1 if --IS_method classic and equal to 2*interaction_limit + 1 if --IS_method gauss, however you can override it through --set_radius option.", action="store_true", default=False, required=False)
parser.add_argument('--set_radius', dest='collision_radius', help="Along with --collision option, here you can set the maximum number of empty loci separating two covered bases regarded as 'colliding'. If not present, you accept to perform collision with default collision radius. You can change it with an int you like at your own risk.", action="store", default=None, required=False, type=int)
parser.add_argument('--tsv', dest='tsv', help="This option produces *.tsv output files too, as soon as allowed (standard matrixes: Redundant and IS); recommended in development or if an highly compatible output was needed.", action="store_true", default=False, required=False)
parser.add_argument('--no_xlsx', dest='no_xlsx', help="This option prevent from generating *.xlsx output file (Excel Workbook). Sometimes it should be useful, e.g. if you are interested only in *.tsv output (using --tsv option) and you want to save as much time as you can.", action="store_true", default=False, required=False)
args = parser.parse_args()
#################################################################################################################################################################################

###Input Parameters for DB_connection###########################
host = args.host    #'172.25.39.2' Alien; #'172.25.39.57' Gemini
user = args.user    #'readonly' #generic user
passwd = args.pw    #'readonlypswd' #generic user pw
port = args.dbport  # 3306 #standard port
################################################################

###Set Up Variables####################################################################
IS_method = args.IS_method
strand_specific_choice = args.strand_specific
#List of available IS methods    
IS_methods_list = ["classic", "gauss"]  #See check_method in Preliminary_controls
                                        #Choose short name!! (see workbook_output
                                        #in output_module - worksheet name 32char)
#######################################################################################




###MAIN###

def main():
        
    ###Set Up Variables##################################################################
    interaction_limit = args.interaction_limit
    alpha = args.alpha 
    bushman_bp_rule = args.bushman_bp_rule # #see Setting-up parameters section below
    delta = args.collision_radius #see Setting-up parameters below
    #####################################################################################
    
    #Print for user                                                                
    print "\n{0}\t***Start***".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Variables to verify user's requests make sense
    check = True ## internal variable for checking variables/controls
    reason = " unexpected error. Try checking syntax / DB connection availability."
    
    #Print for user                                                                
    print "\n{0}\t[INPUT CHECKING] ... ".format((strftime("%Y-%m-%d %H:%M:%S", gmtime()))),    
    
    #Calling functions from Preliminary_controls module, to verify user's requests make sense
    check, reason = Preliminary_controls.smart_check (args.dbDataset, args.collision, args.collision_radius, host, user, passwd, port, args.columns, args.columnsToGroup, IS_method, bushman_bp_rule, IS_methods_list, interaction_limit, alpha, strand_specific_choice, check, reason)
           
    #CHECK AND Preliminary Operations to PROGRAM CORE CALLS    
    if (check == True):
        
        #Print for user                                                                
        print "OK!"
            
               
        #Setting-up parameters
    
        # Bushman bp Rule
        if (IS_method == "gauss"):
            bushman_bp_rule = int((2*int(interaction_limit)) + 1) # Gauss IS mode: override user's bushman_bp_rule
        else:
            bushman_bp_rule = int(bushman_bp_rule) #Good for classic and for general purpose
        
        # Delta (collision_radius)    
        if (delta == None):
            if (IS_method == "gauss"):
                delta = bushman_bp_rule #Set Defaults
            if (IS_method == "classic"):
                delta = bushman_bp_rule + 1 #Set Defaults, as SciencePaper 
        else:
            delta = int(delta) # User input
        
                                
        #Preparing dbDataset_tuple_list
        dbDataset_tuple_list = [] # [('dbschema1', 'dbtable1'), ('dbschema2', 'dbtable2'), ...]        
        
        dbDataset_split = args.dbDataset.split(",")
        for db_string in dbDataset_split:
            db_tupla = None
            db_split = db_string.split(".")
            db_tupla = (db_split[0],db_split[1])
            dbDataset_tuple_list.append(db_tupla)
        
        #Initialize list_of_IS_results_tuple_for_collision
        list_of_IS_results_tuple_for_collision = []   # [(IS_matrix_file_name1, IS_matrix_as_line_list1),(IS_matrix_file_name2, IS_matrix_as_line_list2), ...], list of tuple for PROGRAM CORE
                                                        # empty only if user asks for no collision and no xlsx output. 
        
        #Initialize list_of_result_dictionaries
        list_of_result_dictionaries = []
        # [{
        #    'dataset_name':dbschema.dbtable,
        #    'redundant_matrix':redundant_matrix_as_line_list,
        #    'IS_matrix':IS_matrix_as_line_list
        #    'IS_matrix_collided':IS_matrix_as_line_list_collided / None
        #    'list_of_Covered_Bases':list_of_Covered_Bases
        #    'list_of_Covered_bases_ensambles':list_of_Covered_bases_ensambles
        #    'IS_list':IS_list
        #    'IS_method': IS_method
        #    'strand_specific_choice':strand_specific_choice
        # }, {...}, ...]
        
        #Loop for PROGRAM CORE over tuple in dbDataset_tuple_list
        loop_to_do = len(dbDataset_tuple_list)
        i=1
        
        for db_tupla in dbDataset_tuple_list: # for each tuple in the list -> schema.table
            
            db = db_tupla[0]
            db_table = db_tupla[1]
            #Print for user
            print "\n\n\n{0}\t[START]\tDataset {1} of {2}: {3} - {4}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())), i, loop_to_do, db, db_table)
            
            #PROGRAM_CORE CALLINGS########################################################################################################################
            
            #Case of xlsx output (default) OR collision request: collect results in list_of_IS_results_tuple_for_collision to produce output at the end
            #(TSV output request is handled by PROGRAM_CORE itself, through args.tsv global argument
            if ((args.collision == True) or (args.no_xlsx == False)): 
                IS_matrix_file_name, IS_matrix_as_line_list, result_dictionary = PROGRAM_CORE(db, db_table, bushman_bp_rule, interaction_limit, alpha)
                list_of_IS_results_tuple_for_collision.append((IS_matrix_file_name, IS_matrix_as_line_list, result_dictionary['dataset_name']))
                list_of_result_dictionaries.append(result_dictionary)
                del IS_matrix_file_name, IS_matrix_as_line_list, result_dictionary
            
            #Case of no-collision AND no-xlsx output: nothing needed, PROGRAM_CORE() can handle all task by its own, also TSV output generation on request
            else:
                PROGRAM_CORE(db, db_table, bushman_bp_rule, interaction_limit, alpha)
            
            ##############################################################################################################################################
            
            #Print for user
            print "\n{0}\t[SUCCESSFUL COMPLETED]\tDataset {1} of {2}: {3} - {4}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())), i, loop_to_do, db, db_table)
            #Cycle counter
            i+=1

        #Here you have finished, if collision = False AND no_xlsx == True. Else, you find IS results for each dataset in list_of_IS_results_tuple_for_collision
                        
        ### COLLISION step ########################################################################
        if (args.collision == True):
            
            #Cast
            delta = int(delta)
                        
            #Print for user
            print "\n\n\n{0}\t[COLLISIONS STEP]".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
            
            #Loop over each dataset tupla in list_of_IS_results_tuple_for_collision: [(IS_matrix_file_name1, IS_matrix_as_line_list1, result_dictionary1['dataset_name']),(IS_matrix_file_name2, IS_matrix_as_line_list2, result_dictionary2['dataset_name']), ...]
            i=0 #counter (used to: 1) choose the right result_dictionary in list_of_result_dictionaries 2) choose a better name for screen printing)
            for current_dataset_tuple in list_of_IS_results_tuple_for_collision:
                
                #prepare name_to_print for screen printing
                name_to_print = dbDataset_tuple_list[i][0] + " - " + dbDataset_tuple_list[i][1]
                
                #Computing collision -> results in current_dataset_IS_matrix_file_name, current_dataset_IS_matrix_as_line_list_collided
                print "\n{0}\tComputing data for {1} ... ".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())), name_to_print)
                current_dataset_IS_matrix_file_name, current_dataset_IS_matrix_as_line_list_collided = Collision.multiple_collision(current_dataset_tuple, list_of_IS_results_tuple_for_collision, delta)
                print "{0}\tDone!".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
                #Collisions completed
                
                #Updating related result_dictionary (key = 'IS_matrix_collided')
                list_of_result_dictionaries[i]['IS_matrix_collided'] = current_dataset_IS_matrix_as_line_list_collided
                                
                #Create *.tsv output file, on request (--tsv option)
                if (args.tsv == True):
                    print "\n{0}\tCreating TSV output file in place ...".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
                    output_module.tsv_output(current_dataset_IS_matrix_file_name, current_dataset_IS_matrix_as_line_list_collided)
                    print "{0}\tDone -> {1}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())), current_dataset_tuple[0])
                    
                #update i counter
                i+=1
            
            #Print for user
            print "\n{0}\t[COLLISIONS COMPLETED]".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
            
        # Here you have finished, if no_xlsx == True. Else output creation starts 
            
        ### EXCEL OUTPUT STEP#######################################
        if (args.no_xlsx == False):
            #Print for user
            print "\n\n\n{0}\t[EXCEL OUTPUT CREATION]".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
            
            #Loop over list_of_result_dictionaries
            for result_dictionary in list_of_result_dictionaries:
                print "\n{0}\tCreating {1} ...".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()), "Integration_Analysis_" +  result_dictionary['dataset_name'].replace(".", "_") + ".xlsx")
                output_module.workbook_output(result_dictionary, host, user, passwd, port)
                print "{0}\tDone!".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
                
            print "\n{0}\t[OUTPUT CREATED]".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))

        #################################################################
        # Here you have finished. Place here the code for further tasks #
        #################################################################
        
        # All tasks finished, goodbye message
        print "\n\n{0}\t***All tasks have finished***\n\n\t[QUIT].\n".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        
    
    else: #Check == False from Preliminary_controls.smart_check - The program doesn't start, de facto.
        print "\n\nYour request can't be processed: {0}".format(reason)
        print "\n\t[QUIT].\n\n"






###PROGRAM CORE###

def PROGRAM_CORE(db, db_table, bushman_bp_rule, interaction_limit, alpha):
    
    #Output file name template
    file_output_name = db + "_" + db_table + ".tsv"
    
    #Preparing queries to DB
    query_for_columns=Common_Functions.prepareSELECT(args.columns)   #such as "`sample`,`tissue`,`treatment`"
    reference_genome = args.reference_genome
    query_step = long(args.query_steps)
    
    
    ###Retrieving data from DB: reads_data_dictionary and lam_data_dictionay ###############################################################################################
    
    # Check n_table rows.
    connection = DB_connection.dbOpenConnection (host, user, passwd, port, db)
    n_table_rows = DB_connection.getTableRowCount (connection, db_table)
    DB_connection.dbCloseConnection(connection)
    
    # Initialize output data dictionary
    lam_data_dictionay = None
    reads_data_dictionary = None
    
    # If n_table_rows > rowthreshold, then use file dump and not DB access
    
    if (n_table_rows < args.rowthreshold): # Retrieving data DIRECTLY from DB
        print "\n{0}\tRetrieving data from DB ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
        
        #reads_data_dictionary 
        connection = DB_connection.dbOpenConnection (host, user, passwd, port, db) # init connection to DB for importing data
        reads_data_dictionary = DB_connection.import_reads_data_from_DB(connection, db_table, query_step, reference_genome)
        DB_connection.dbCloseConnection(connection) # close connection to DB
        
        #lam_data_dictionay
        connection = DB_connection.dbOpenConnection (host, user, passwd, port, db) # init connection to DB for importing data
        lam_data_dictionay  = DB_connection.import_lam_data_from_DB(connection, db_table, query_step, reference_genome)
        DB_connection.dbCloseConnection(connection) # close connection to DB
   
    else: # Retrieving data from DB passing through a tmpfile
        print "\n{0}\tRetrieving data from DB, converting into file and parsing data ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
       
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
    ########################################################################################################################################################################    
    
       
    ###Creating ordered_keys_for_reads_data_dictionary####################################################################################################################
    ###ordered_keys_for_reads_data_dictionary is a list of keys for reads_data_dictionary, useful for retrieving reads ordered by chromosome and then integration_locus###
    ###'ORDER' IS THE STRING'S ONE, so alphabetical and not 'along genome'###
    print "\n{0}\tOrdering retrieved data ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
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
    
    print "\n{0}\tArranging data structure ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
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
    
    print "\n{0}\tCreating data schema according to user request ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Retrieving labels for matrix columns used for computing data, labels as user wishes and a dictionary to relate them
    column_labels, user_label_dictionary = DB_connection.get_column_labels_from_DB(host, user, passwd, port, db, db_table, parameters_list, query_for_columns)
    
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
    
    #Tell user this task has started
    print "\n{0}\tProcessing redundant Reads as Matrix ...".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    
    #Create redundant reads matrix as list and prepare name for output
    redundant_matrix_file_name, redundant_matrix_as_line_list = Matrix_creation.simple_redundant_matrix(list_of_Covered_Bases, column_labels, file_output_name, strand_specific = strand_specific_choice)
        
    #Convert matrix according to user's requests
    redundant_matrix_as_line_list = Common_Functions.convert_matrix(redundant_matrix_as_line_list, user_label_dictionary, user_merged_labels_dictionary)
    
    #Tell user this task has finished
    print "{0}\tDone!".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
    
    #Create *.tsv output file, on request (--tsv option)
    if (args.tsv == True):
        
        #Tell user this task has started
        print "\n{0}\tCreating TSV output file in place ...".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        
        #TSV output file creation
        output_module.tsv_output(redundant_matrix_file_name, redundant_matrix_as_line_list)
               
        #Tell user this task has finished
        print "{0}\t*Redundant Reads Matrix file has been created --> {1}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())),redundant_matrix_file_name)
    
    #################################################################################################################################################
    
     
    #Grouping Covered Bases in ENSEMBLES#######################################################################################################################################
    
    #Tell user this task has started
    print "\n{0}\tGrouping Covered Bases in Ensembles ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #List of results: list_of_Covered_bases_ensambles
    list_of_Covered_bases_ensambles = []
    
    
    #Ensemble grouping
    
    #Creating first covered_bases_ensemble with first_covered_base     
    current_covered_bases_ensemble = None
    
    
    #If strand_specific_choice == False, algorithm goes straight
    if (strand_specific_choice == False):
        
        #Creating first covered_bases_ensemble with first_covered_base     
        current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(list_of_Covered_Bases[0], strand_specific = strand_specific_choice)
    
        #different if user chooses "classic" method to retrieving IS
        if (IS_method == "classic"):
            for covered_base in list_of_Covered_Bases[1:]:
                dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
                if ((dist == "undef") or (dist > bushman_bp_rule) or (current_covered_bases_ensemble.spanned_bases == (bushman_bp_rule + 1)) or ((current_covered_bases_ensemble.spanned_bases + dist)>(bushman_bp_rule + 1))):
                    list_of_Covered_bases_ensambles.append(current_covered_bases_ensemble)
                    current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base, strand_specific = strand_specific_choice)
                else:
                    current_covered_bases_ensemble.push_in(covered_base)
                    
        #for the others possible user's choices, groping method goes straight
        else:
            for covered_base in list_of_Covered_Bases[1:]:
                dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
                if ((dist == "undef") or (dist > bushman_bp_rule)):
                    list_of_Covered_bases_ensambles.append(current_covered_bases_ensemble)
                    current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base, strand_specific = strand_specific_choice)
                else:
                    current_covered_bases_ensemble.push_in(covered_base)
                
        list_of_Covered_bases_ensambles.append(current_covered_bases_ensemble) #APPEND LAST ENSEMBLE
        
        # NOW COVERED BASES ENSEMBLES ARE IN AN ORDERED LIST: list_of_Covered_bases_ensambles
    
                    
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
            
            #different if user chooses "classic" method to retrieving IS
            if (IS_method == "classic"):
                for covered_base in dic_of_list_of_Covered_Bases[current_strand][1:]:
                    dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
                    if ((dist == "undef") or (dist > bushman_bp_rule) or (current_covered_bases_ensemble.spanned_bases == (bushman_bp_rule + 1)) or ((current_covered_bases_ensemble.spanned_bases + dist)>(bushman_bp_rule + 1))):
                        list_of_Covered_bases_ensambles_current_strand.append(current_covered_bases_ensemble)
                        current_covered_bases_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base)
                    else:
                        current_covered_bases_ensemble.push_in(covered_base)
            
            #for the others possible user's choices, groping method goes straight
            else:
                for covered_base in dic_of_list_of_Covered_Bases[current_strand][1:]:
                    dist = current_covered_bases_ensemble.Covered_bases_list[-1].distance(covered_base)
                    if ((dist == "undef") or (dist > bushman_bp_rule)):
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
        
        # NOW COVERED BASES ENSEMBLES ARE IN AN ORDERED LIST: list_of_Covered_bases_ensambles
        
    #Tell user this task has finished    
    print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
        
    ###########################################################################################################################################################################        
    
    
    #Integration Sites Retrieval###############################################################################################################################################        
    
    #Tell user this task has started
    print "\n{0}\tComputing Integration Sites over Covered Bases Ensembles ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Initialize list of results:
    IS_list = []
    
    #Classic method
    if (IS_method == "classic"):
        for Covered_bases_ensamble in list_of_Covered_bases_ensambles:
            IS_list.append(Integration_Sites_retrieving_methods.classic(Covered_bases_ensamble, strand_specific = strand_specific_choice))
            
    
    #NOW INTEGRATION SITES RETRIEVED THROUGH "CLASSIC" METHOD ARE IN IS_LIST

    #Gaussian_IS_identification method    
    if (IS_method == "gauss"):
        bin_boundaries, bin_areas, diagnostic = Function_for_Gaussian_IS_identification.gaussian_histogram_generator(interaction_limit, alpha)
        del bin_boundaries, diagnostic
        hist_gauss_normalized_to_peak = Function_for_Gaussian_IS_identification.normalize_histogram_to_the_peak(bin_areas, interaction_limit)
        for Covered_bases_ensamble in list_of_Covered_bases_ensambles:
            # IS_list = IS_list + Integration_Sites_retrieving_methods.Gaussian_IS_identification(Covered_bases_ensamble, hist_gauss_normalized_to_peak, interaction_limit, strand_specific_choice)
            IS_list = IS_list + Integration_Sites_retrieving_methods.refined_Gaussian_IS_identification(Covered_bases_ensamble, hist_gauss_normalized_to_peak, interaction_limit, strand_specific_choice)
    
    #NOW INTEGRATION SITES RETRIEVED THROUGH "GAUSS" METHOD ARE IN IS_LIST
        
    #Whatever method    
    if (IS_method == "whatever"):
        ###Here the code, when "whatever" new method will be available
        pass
        
    #NOW INTEGRATION SITES RETRIEVED THROUGH "WHATEVER" METHOD ARE IN IS_LIST
    
    print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    ###########################################################################################################################################################################        
    
   
    #IS matrix creation ###############################################################################################
    
    print "\n{0}\tProcessing Integration Sites as Matrix ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Create IS matrix as list and prepare output file name
    IS_matrix_file_name, IS_matrix_as_line_list = Matrix_creation.simple_IS_matrix(IS_list, column_labels, file_output_name, IS_method, strand_specific=strand_specific_choice)
    
    #Convert matrix according to user's requests
    IS_matrix_as_line_list = Common_Functions.convert_matrix(IS_matrix_as_line_list, user_label_dictionary, user_merged_labels_dictionary)
    
    #Tell user this task has finished
    print "{0}\tDone!".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    
    #Create *.tsv output file, on request (--tsv option), if possible (NO --collision option)
    if (args.tsv == True):
        
        if (args.collision == False): #NO collision is necessary to create IS matrix in place
        
            #Tell user this task has started
            print "\n{0}\tCreating TSV output file in place ...".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())))
            
            #TSV output file creation
            output_module.tsv_output(IS_matrix_file_name, IS_matrix_as_line_list)
            
            #Tell user this task has finished
            print "{0}\t*IS Matrix file has been created --> {1}".format((strftime("%Y-%m-%d %H:%M:%S", gmtime())),IS_matrix_file_name)
            
        else:
            #Tell user he has to wait, due to --collision option
            print "\n{0}\t* IS Matrix TSV file won't be created until all datasets have been processed, due to --collision request.".format(strftime("%Y-%m-%d %H:%M:%S", gmtime()))
        
    ####################################################################################################################
    
    
    #Result Dictionary #################################################################################################
    
    # result_dictionary =
    # {
    #    'dataset_name':dbschema.dbtable,
    #    'redundant_matrix':redundant_matrix_as_line_list,
    #    'IS_matrix':IS_matrix_as_line_list
    #    'IS_matrix_collided':IS_matrix_as_line_list_collided / None
    #    'list_of_Covered_Bases':list_of_Covered_Bases
    #    'list_of_Covered_bases_ensambles':list_of_Covered_bases_ensambles
    #    'IS_list':IS_list
    #    'IS_method': IS_method
    #    'strand_specific_choice':strand_specific_choice
    # }
    
    result_dictionary = {'dataset_name':db+"."+db_table, 'redundant_matrix':redundant_matrix_as_line_list, 'IS_matrix':IS_matrix_as_line_list, 'IS_matrix_collided': None, 'list_of_Covered_Bases':list_of_Covered_Bases, 'list_of_Covered_bases_ensambles':list_of_Covered_bases_ensambles, 'IS_list':IS_list, 'IS_method': IS_method, 'strand_specific_choice':strand_specific_choice}
    
    #####################################################################################################################
    
        
    #Return Results #####################################################
    return IS_matrix_file_name, IS_matrix_as_line_list, result_dictionary
    #####################################################################
    




# SENINEL
if __name__ == '__main__':
    main()
    
    