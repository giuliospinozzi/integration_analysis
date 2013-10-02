###Header##################################################
header = """

+--------------------------------------------------------+
 Module: DB_connection
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 11th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+--------------------------------------------------------+

 Description:
  - This module contains functions to manage DB connection
    and provide quick solutions to common tasks involving
    DB data retrieval
  
 Note:
  - None

---------------------------------------------------------- 
""" 
##########################################################


###Requested Package(s) Import###
import MySQLdb
#################################





def dbOpenConnection (host, user, passwd, port, db):
    """
    *** Open DB connection ***
    Input: DB data connection details
    Output: MySQLdb connection object
    """
    conn = MySQLdb.connect( 
     host = host,
     user = user,
     passwd = passwd,
     port = port,
     db = db,
     )
    
    # Return result
    return conn




def dbCloseConnection (conn):
    """
    *** Close the connection to DB ***
    Input: MySQLdb connection object
    Output: None   
    """
    conn.close()




def getTableRowCount (conn, db_table):
    """
    *** Count rows of 'db_table' throgh 'conn' DB connection
    Input: MySQLdb connection object, target table
    Output: integer (target table row count)
    """
    cursor = conn.cursor (MySQLdb.cursors.Cursor)  
    cursor.execute ("SELECT count( * ) FROM %s WHERE 1" %(db_table))
    n_rows = cursor.fetchall()[0][0]
    cursor.close()
    
    # Return result
    return n_rows
    



def import_reads_data_from_DB (conn, db_table, query_step=1000000, reference_genome="unspecified genome"):
    """
    *** Get Reads data in the form of Dictionary, directly from DB ***
    
    INPUT: conn - MySQLdb connection object (you might use 'dbOpenConnection' function)
           db_table - String containing table you want to interrogate (schema was set in 'conn')
    
    OPTIONAL INPUT: query_step - Integer. It fixes the number of rows fetched at a time (useful to reduce memory usage peaks)
                    reference_genome - String. Something like 'hg19'.
        
    OUTPUT: reads_data_dictionary - Dictionary of the form: Key = header; Item = (reference_genome, chr, strand, integration_locus, read end, span, lam_id)
    
    lOGICS: Given a MySQLdb connection object (it already contains DB\DB_schema target), this function perform queries to return reads data in form of a dictionary
            (details explained in 'OUTPUT' section). Queries can be splitted to return only 'query_step'-results-a-time: this can be useful to reduce memory
            usage peaks
            
    !! WARNINGS !! : 1) This function perform a query like the following:
                        "SELECT `header`, `chr`, `strand`, `integration_locus`, 100 as `span`, `complete_name` as lam_id  FROM [...]"
                        PLEASE NOTE *** 100 as `span` ***, necessary to avoid errors when span is NULL: on the other side, this fact could generate
                        WRONG VALUES for READ END element in reads_data_dictionary tuple items
                        
                     2) 'read end' is not retrieved but calculated as 'start + span': conversely 'start' and 'span' are directly retrieved from DB
    """    
    # Initialize reads_data_dictionary to collect results
    reads_data_dictionary={}
    
    # Prepare queries-splitting to reduce memory usage peak
    n_rows = getTableRowCount (conn, db_table)
    splitting = range(0, n_rows, query_step) # by default one million row a time (query_step=1000000)
    
    # Cycle to perform splitted queries
    for n in splitting:
            
        # Create dictionary cursor
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)
        
        # Query for Reads Data and cursor closing
        start = str(n)
        end = str((n+query_step))
        cursor.execute("SELECT `header`, `chr`, `strand`, `integration_locus`, 100 as `span`, `complete_name` as lam_id  FROM {0} WHERE 1 LIMIT {1}, {2}".format(db_table, start, end)) # 100 as `span` to prevent errors due to "NULL" span
        reads_data = cursor.fetchall()
        cursor.close()
        
        #Build reads data dictionary ('reads_data_dictionary' -> return)
        for dat in reads_data:
            reads_data_dictionary[dat['header']]=(reference_genome, dat['chr'], dat['strand'], long(dat['integration_locus']), dat['integration_locus'] + dat['span'], dat['span'], dat['lam_id'])
        del reads_data
    
    # Return result
    return reads_data_dictionary




def import_lam_data_from_DB_lam (conn, db_table, query_step=1000000, reference_genome="unspecified genome"):
    """
    *** Get LAM data in the form of Dictionary, directly from DB ***
    
    INPUT: conn - MySQLdb connection object (you might use 'dbOpenConnection' function)
           db_table - String containing table you want to interrogate (schema was set in 'conn')
    
    OPTIONAL INPUT: query_step - Integer. It fixes the number of rows fetched at a time (useful to reduce memory usage peaks)
                    reference_genome - String. Something like 'hg19'.
        
    OUTPUT: lam_data_dictionary - Dictionary of the form: Key = lam_id; Item = (n_LAM, tag, pool, tissue, sample, treatment, group_name, enzyme)
    
    lOGICS: Given a MySQLdb connection object (it already contains DB\DB_schema target), this function perform queries to return LAM data in form of a dictionary
            (details explained in 'OUTPUT' section). Queries can be splitted to return only 'query_step'-results-a-time: this can be useful to reduce memory
            usage peaks
            
    !! WARNINGS !! :  This function perform a query like the following:
                      "SELECT DISTINCT `complete_name` as lam_id, [...] FROM [...]"
                      PLEASE NOTE *** `complete_name` as lam_id *** because up to now 'lam_id' data are labelled as 'complete_name' in DB
    """
    #Initialize lam_data_dictionay to collect results ('lam_data_dictionay' -> return)
    lam_data_dictionay={}
    
    # Query for Lam Data
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    cursor.execute("SELECT DISTINCT `complete_name` as lam_id,`n_LAM`, `tag`, `pool`, `tissue`, `sample`, `treatment`, `group_name`, `enzyme`  FROM {0} WHERE 1".format(db_table))
    lam_data = cursor.fetchall()
    cursor.close()
    
    #Preparing 'treatment' padding, MANDATORY to preserve alphabetical order
    len_max = 0
    for dat in lam_data:
        if (len(dat['treatment']) > len_max):
            len_max = len(dat['treatment'])
    
    # Build lam data dictionary ('lam_data_dictionary' -> return)
    for dat in lam_data:
        dat['treatment'] = str(dat['treatment']).zfill(len_max) #Padding
        lam_data_dictionay[dat['lam_id']]=(dat['n_LAM'], dat['tag'], dat['pool'], dat['tissue'], dat['sample'], dat['treatment'], dat['group_name'], dat['enzyme'])
    del lam_data
           
    # Return result
    return lam_data_dictionay
    



def get_column_labels_from_DB (host, user, passwd, port, db, db_table, parameters_list, query_for_columns):
    '''
    *** Get all possible column labels from DB, according to categories given as --columns argument ***
    *** Build a dictionary to translate my label (needed for computation) into user labels (change category order according to user wishes, --columns argument order) ***
    
    INPUT: host, user, passwd, port, db, db_table - user input to set up DB connection and queries (args.host, args.user, args.pw, args.dbport;  db, db_table from args.dbDataset.split(","))
           parameters_list - a list of desired category to account for in label construction (from query_for_columns.split(","))
           query_for_columns - String given in input as --columns argument (such as "sample, tissue, treatment", query_for_columns = args.columns)
           
    OUTPUT: column_labels_list - List of all distinct (my) label buildable through categories in parameters_list/query_for_columns
            user_label_dictionary - Dictionary needed to translate my just retrieved 'rigid labels' (here called 'my labels', 'column labels' or simply 'labels') into user labels, according to his wishes (--columns input order)
                                    Key = my label; Item = [user label, user label as tuple, my label as tuple]
    
    LOGIC: 1) Imports and returns column labels for final matrix, that is e.g. all the possible 'sample_tissue_treatment'-like labels buildable with data in db, db_table (column_labels_list)
           2) Creates and returns user_label_dictionary, needed to translate my just retrieved 'rigid labels' (here called 'my labels', 'column labels' or simply 'labels' - strictly needed as they are due to computational reasons)
              into user labels, according to his wishes (refplecting --columns input order). Label as tupla contained in user_label_dictionary items are present only for computational reasons. 
    
    '''    
    # Setting Up Connection to DB and creating cursor 
    conn = dbOpenConnection (host, user, passwd, port, db, db_table, )
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    
    # Query for column labels
    cursor.execute("SELECT DISTINCT {0} FROM {1} WHERE 1".format(query_for_columns, db_table))
    column_labels = cursor.fetchall()
    cursor.close()
    
    # Close DB Connection
    dbCloseConnection (conn)
    
    # User_label_template (list)
    user_label_template = parameters_list
    
    # user_label_dictionary
    user_label_dictionary = {} #of kind: key-> my label; item -> [user label, user label as tuple, my label as tuple]
       
    # Initialize column labels list, ordered at the end ('column_labels_list' -> return)
    column_labels_list=[]
    
    # Preparing padding, MANDATORY to preserve alphabetical order
    len_max = 0
    for dat in column_labels:
        if (len(dat['treatment']) > len_max):
            len_max = len(dat['treatment'])
    
    # Build column labels list
    for dat in column_labels:
        if ('treatment' in parameters_list):
            dat['treatment'] = str(dat['treatment']).zfill(len_max) #Padding
        # Create label: labels building has to keep freezed like this because they need to match with labels builded in "Classes_for_Integration_Analysis" module (class Covered_base)
        label = ""
        label_as_tupla = ()
        if ("group_name" in parameters_list):
            label = label + "_" + dat['group_name']
            label_as_tupla = label_as_tupla + (dat['group_name'],)
        if ("n_LAM" in parameters_list):
            label = label + "_" + dat['n_LAM']
            label_as_tupla = label_as_tupla + (dat['n_LAM'],)
        if ("pool" in parameters_list):
            label = label + "_" + dat['pool']
            label_as_tupla = label_as_tupla + (dat['pool'],)
        if ("tag" in parameters_list):
            label = label + "_" + dat['tag']
            label_as_tupla = label_as_tupla + (dat['tag'],)
        if ("enzyme" in parameters_list):
            label = label + "_" + dat['enzyme']
            label_as_tupla = label_as_tupla + (dat['enzyme'],)
        if ("sample" in parameters_list):
            label = label + "_" + dat['sample']
            label_as_tupla = label_as_tupla + (dat['sample'],)
        if ("tissue" in parameters_list):
            label = label + "_" + dat['tissue']
            label_as_tupla = label_as_tupla + (dat['tissue'],)
        if ("treatment" in parameters_list):
            label = label + "_" + dat['treatment']
            label_as_tupla = label_as_tupla + (dat['treatment'],)
        label = label[1:]
        # Append label to column_labels_list
        column_labels_list.append(label)
        
        # Create user labels (the same as mine but reflecting --colums argument order given in input)
        user_label = ""
        user_label_as_tupla = ()
        for category in user_label_template:
            user_label = user_label + "_" + dat[category]
            user_label_as_tupla = user_label_as_tupla + (dat[category],)
        user_label = user_label[1:]
        
        # Update user_label_dictionary (-> return): this dictionary will be used to "translate" my labels into user labels
        user_label_dictionary.update({label:[user_label, user_label_as_tupla, label_as_tupla]})
            
    # column_labels_list built, column_labels from cursor.fetchall() is became useless
    del column_labels
    
    # Order column_labels_list
    column_labels_list.sort()

    # Return results
    return column_labels_list, user_label_dictionary