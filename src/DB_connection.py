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
import sys
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
    *** Count rows of 'db_table' through 'conn' DB connection ***
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
            # Remove '/1', '/2' from headers
            clean_header = dat['header']
            if "/1" in clean_header:
                clean_header = clean_header.replace("/1", "")
            elif "/2" in clean_header:
                clean_header = clean_header.replace("/2", "")
            reads_data_dictionary[clean_header]=(reference_genome, dat['chr'], dat['strand'], long(dat['integration_locus']), dat['integration_locus'] + dat['span'], dat['span'], dat['lam_id'])
        del reads_data
    
    # Return result
    return reads_data_dictionary




def import_lam_data_from_DB (conn, db_table, query_step=1000000, reference_genome="unspecified genome"):
    """
    *** Get LAM data in the form of Dictionary, directly from DB ***
    
    INPUT: conn - MySQLdb connection object (you might use 'dbOpenConnection' function)
           db_table - String containing table you want to interrogate (schema was set in 'conn')
    
    OPTIONAL INPUT: query_step - Integer. It fixes the number of rows fetched at a time (useful to reduce memory usage peaks)
                    reference_genome - String. Something like 'hg19'.
        
    OUTPUT: lam_data_dictionary - Dictionary of the form: Key = lam_id; Item = (n_LAM, tag, pool, tissue, sample, treatment, group_name, enzyme, vector)
    
    lOGICS: Given a MySQLdb connection object (it already contains DB\DB_schema target), this function perform queries to return LAM data in form of a dictionary
            (details explained in 'OUTPUT' section). Queries can be split to return only 'query_step'-results-a-time: this can be useful to reduce memory
            usage peaks
            
    !! WARNINGS !! :  This function perform a query like the following:
                      "SELECT DISTINCT `complete_name` as lam_id, [...] FROM [...]"
                      PLEASE NOTE *** `complete_name` as lam_id *** because up to now 'lam_id' data are labeled as 'complete_name' in DB
    """
    #Initialize lam_data_dictionay to collect results ('lam_data_dictionay' -> return)
    lam_data_dictionay={}
    
    # Query for Lam Data
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    cursor.execute("SELECT DISTINCT `complete_name` as lam_id,`n_LAM`, `tag`, `pool`, `tissue`, `sample`, `treatment`, `group_name`, `enzyme`, `vector`  FROM {0} WHERE 1".format(db_table))
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
        lam_data_dictionay[dat['lam_id']]=(dat['n_LAM'], dat['tag'], dat['pool'], dat['tissue'], dat['sample'], dat['treatment'], dat['group_name'], dat['enzyme'], dat['vector'])
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
              into user labels, according to his wishes (reflecting --columns input order). Label as tupla contained in user_label_dictionary items are present only for computational reasons. 
    
    '''    
    # Setting Up Connection to DB and creating cursor 
    conn = dbOpenConnection (host, user, passwd, port, db)
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
    if ('treatment' in parameters_list):
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
        if ("vector" in parameters_list):
            label = label + "_" + dat['vector']
            label_as_tupla = label_as_tupla + (dat['vector'],)
            
        label = label[1:]
        # Append label to column_labels_list
        column_labels_list.append(label)
        
        # Create user labels (the same as mine but reflecting --columns argument order given in input)
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




def retrieve_sequences_from_DB (conn, db_table_for_tracking_raw, db_table_for_tracking_final, query_step=1000000):
    """
    *** Get SEQUENCE data in the form of Dictionary, directly from DB ***
                            FOR SEQTRACKER USE
    
    INPUT: conn - MySQLdb connection object (you might use 'dbOpenConnection' function)
           db_table_for_tracking_raw - String containing table you want to interrogate (schema was set in 'conn')
           db_table_for_tracking_final - String containing table you want to interrogate (schema was set in 'conn')
           (note that the schema is supposed to be the same for both tables)
    
    OPTIONAL INPUT: query_step - Integer. It fixes the number of rows fetched at a time (useful to reduce memory usage peaks)
                    SEE NOTES
    
    OUTPUT: raw/final _read_dictionary - Dictionaries of sequence data (raw reads and final reads)
                                         Key = read header; Item = sequence
    
    NOTES: up to now query_step splitting is not yet implemented!!!
    """
    
    # Dictionaries of results
    raw_read_dictionary = {}
    final_read_dictionary = {}
    
    # Query for raw data
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    cursor.execute("SELECT `HEADER`, `SEQUENCE` FROM {0} WHERE 1".format(db_table_for_tracking_raw))
    raw_read_data = cursor.fetchall()
    cursor.close()
    # Fill raw_read_dictionary
    for dat in raw_read_data:
        raw_read_dictionary[dat['HEADER']] = dat['SEQUENCE']
    del raw_read_data
    
    # Query for final data
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    cursor.execute("SELECT `prod_header`, `isread_nasequence` FROM {0} WHERE 1".format(db_table_for_tracking_final))
    final_read_data = cursor.fetchall()
    cursor.close()
    # Fill raw_read_dictionary
    for dat in final_read_data:
        final_read_dictionary[dat['prod_header']] = dat['isread_nasequence']
    del final_read_data
    
    # Return results
    return raw_read_dictionary, final_read_dictionary




def retrieve_sequences_and_metadata_from_DB (conn, table_to_query, table_kind, header_list, query_split=100000):
    """
    To Do
            
    NOTE: query_split=100000 due to the big piece of data retrieved
    """
    
    # Control
    if table_kind not in ['IS', 'RAW']:
        sys.exit("\n\n\t[ERROR] Bad Function call: table_kind arg not recognized in retrieve_sequences_and_metadata_from_DB (conn, table_to_query, table_kind, header_list, query_split=100000).\tQuit.\n\n")
    
    # Dictionary of results
    dictionary_to_return = {}
    
    # Get n_headers
    n_headers = len(header_list)
    
    # Get n_split
    n_split = n_headers / query_split
    reminder = n_headers - (n_split*query_split)
    
    # Prepare list_of_header_lists
    list_of_header_lists = [None]*n_split
    start=0
    end=query_split
    for n in range(0, n_split):
        list_of_header_lists[n] = header_list[start:end]
        start+=query_split
        end+=query_split
    if (reminder != 0):
        list_of_header_lists.append(header_list[start:])
    
    # metadata and sequences - 'IS' case
    if table_kind == 'IS':
        # Loop for query
        for header_list_chunck in list_of_header_lists:
            
            # Query for metadata and sequences - 'IS'
            cursor = conn.cursor (MySQLdb.cursors.DictCursor)
            headers_for_query = "', '".join(header_list_chunck)
            headers_for_query = "'" + headers_for_query + "'"
            cursor.execute("SELECT `prod_header`, `isread_start`, `isread_end`, `isread_strand`, `isread_cigar`, `isread_MD`, `isread_nasequence` FROM {0} WHERE `prod_header` IN ({1})".format(table_to_query, headers_for_query))
            IS_data = cursor.fetchall()
            cursor.close()
    
            # Fill chunck_of_dictionary_to_return
            chunck_of_dictionary_to_return = {}
            for dat in IS_data:
                seq_len = abs(dat['isread_start'] - dat['isread_end'])
                item_dict = {'isread_cigar': dat['isread_cigar'], 'isread_MD': dat['isread_MD'], 'isread_nasequence': dat['isread_nasequence'], 'seq_len': seq_len, 'isread_strand': dat['isread_strand']}
                chunck_of_dictionary_to_return[dat['prod_header']] = item_dict
            
            # Update dictionary_to_return with the current chunck
            dictionary_to_return.update(chunck_of_dictionary_to_return)
            
    if table_kind == 'RAW':
        # Loop for query
        for header_list_chunck in list_of_header_lists:
            
            # Query for metadata and sequences - 'IS'
            cursor = conn.cursor (MySQLdb.cursors.DictCursor)
            headers_for_query = "', '".join(header_list_chunck)
            headers_for_query = "'" + headers_for_query + "'"
            cursor.execute("SELECT `HEADER`, `SEQUENCE` FROM {0} WHERE `HEADER` IN ({1})".format(table_to_query, headers_for_query))
            RAW_data = cursor.fetchall()
            cursor.close()
    
            # Fill chunck_of_dictionary_to_return
            chunck_of_dictionary_to_return = {}
            for dat in RAW_data:
                chunck_of_dictionary_to_return[dat['HEADER']] = dat['SEQUENCE'][20:]  # first 20 nucleotides removed!!!
            
            # Update dictionary_to_return with the current chunck
            dictionary_to_return.update(chunck_of_dictionary_to_return)
        
    
    # Return results
    return dictionary_to_return




def dropTable (host, user, passwd, port, db, db_table_list):
    
    # Open database connection
    conn = dbOpenConnection (host, user, passwd, port, db)
    # prepare a cursor object using cursor() method
    cursor = conn.cursor()
    # Drop table if it already exist using execute() method.
    db_table_list[:] = ["`"+tablename+"`" for tablename in db_table_list]
    db_table_list = ', '.join(db_table_list)
    sql = "DROP TABLE {}".format(db_table_list)
    cursor.execute(sql)
    # Close database connection
    dbCloseConnection (conn)
    
    
    










