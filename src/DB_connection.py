###Header################################################
header = """

+------------------------------------------------------+
 Module: DB_connection
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 11th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - [...]
  
 Note:
  - [...]

-------------------------------------------------------- 
""" 
########################################################


###Requested Package(s) Import###
import MySQLdb
#################################




def dbOpenConnection (host, user, passwd, port, db, db_table, ):
    """
    Input: DB data connection details
    Output: connection object
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
    Close the connection to DB
    """
    conn.close()




def getTableRowCount (conn, db_table):
    """
    Input: connection object, target table
    Output: integer with row count
    """
    # Splitting query to reduce memory usage peak
    cursor = conn.cursor (MySQLdb.cursors.Cursor)  
    cursor.execute ("SELECT count( * ) FROM %s WHERE 1" %(db_table))
    n_rows = cursor.fetchall()[0][0]
    cursor.close()
    
    # Return result
    return n_rows
    



def import_data_from_DB_reads (conn, db_table, query_step=1000000, reference_genome="unspecified genome"):
    """
    [...]
    """    
    # Reads_query dictionary to collect results
    reads_query={}
    
    # Splitting query to reduce memory usage peak
    cursor = conn.cursor (MySQLdb.cursors.Cursor)
    cursor.execute ("SELECT count( * ) FROM {0} WHERE 1".format(db_table))
    n_rows = cursor.fetchall()[0][0]
    cursor.close()
    splitting = range(0, n_rows, query_step) # by default one million row a time
    
    # Cycle to perform queries
    for n in splitting:
            
        # Open DB Connection
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)
        
        # Query for Reads Data
        start = str(n)
        end = str((n+query_step))
        cursor.execute("SELECT `header`, `chr`, `strand`, `integration_locus`, 100 as `span`, `complete_name` as lam_id  FROM {0} WHERE 1 LIMIT {1}, {2}".format(db_table, start, end)) # 100 as `span` to prevent errors due to "NULL" span
        reads_data = cursor.fetchall()
        cursor.close()
        
        #Build reads data dictionary ('reads_query' -> return)
        for dat in reads_data:
            reads_query[dat['header']]=(reference_genome, dat['chr'], dat['strand'], long(dat['integration_locus']), dat['integration_locus'] + dat['span'], dat['span'], dat['lam_id'])
        del reads_data
    
    # Return result
    return reads_query




def import_data_from_DB_lam (conn, db_table, query_step=1000000, reference_genome="unspecified genome"):
    """
    [...]
    """
    # Query for Lam Data
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    cursor.execute("SELECT DISTINCT `complete_name` as lam_id,`n_LAM`, `tag`, `pool`, `tissue`, `sample`, `treatment`, `group_name`, `enzyme`  FROM {0} WHERE 1".format(db_table))
    lam_data = cursor.fetchall()
    cursor.close()
    
    #Build lam data dictionary ('lam_query' -> return)
    lam_query={}
    
    #Preparing padding, MANDATORY to preserve alphabetical order
    len_max = 0
    for dat in lam_data:
        if (len(dat['treatment']) > len_max):
            len_max = len(dat['treatment'])
    
    for dat in lam_data:
        dat['treatment'] = str(dat['treatment']).zfill(len_max) #Padding
        lam_query[dat['lam_id']]=(dat['n_LAM'], dat['tag'], dat['pool'], dat['tissue'], dat['sample'], dat['treatment'], dat['group_name'], dat['enzyme'])
    del lam_data
           
    # Return result
    return lam_query
    



def get_column_labels_from_DB (host, user, passwd, port, db, db_table, parameters_list, query_for_columns, reference_genome):
    '''
    Imports and returns column labels for final matrix e.g. 'sample_tissue_treatment' -> column_labels_list
    Creates and returns user_label_dictionary, needed to translate my 'rigid labels' just retrieved as user wishes (input format)
    
    '''    
    # Setting Up Connection to DB
    conn = dbOpenConnection (host, user, passwd, port, db, db_table, )
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    
    # Query for column labels
    cursor.execute("SELECT DISTINCT {0} FROM {1} WHERE 1".format(query_for_columns, db_table))
    column_labels = cursor.fetchall()
    cursor.close()
    
    # User_label_template
    user_label_template = parameters_list
    
    # user_label_dictionary
    user_label_dictionary = {} #of kind: key-> my label, item -> [user label, user label as tupla, my label as tupla]
       
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
            
    # column_labels_list built, column_labels useless
    del column_labels
    # Order column_labels_list
    column_labels_list.sort()

    # Close DB Connection
    conn.close()

    # Return results
    return column_labels_list, user_label_dictionary