# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 15:13:00 2014

@author: stefano
"""

###Requested Package(s) Import###
import sqlite3
#################################


def dbOpenConnection (db):
    """
    *** Open DB connection ***
    Input: DB complete path
    Output: sqlite3 connection object
    """
    conn = sqlite3.connect(db, timeout=1440)
    return conn
    
def dbCloseConnection (conn):
    """
    *** Close the connection to DB ***
    Input: sqlite3 connection object
    Output: None   
    """
    conn.close()
    
def getTableRowCount (conn, db_table):
    """
    *** Count rows of 'db_table' through 'conn' DB connection ***
    Input: sqlite3 connection object, target table
    Output: integer (target table row count)
    """
    cursor = conn.cursor()  
    cursor.execute ("SELECT count( * ) FROM %s WHERE 1" %(db_table))
    n_rows = cursor.fetchall()[0][0]
    cursor.close()
    return n_rows
    
def import_reads_data_from_DB (conn, db_table, reference_genome="unspecified genome"):
    """
    *** Get Reads data in the form of Dictionary, directly from DB ***
    
    INPUT: conn
           db_table
    
    OPTIONAL INPUT: reference_genome - String. Something like 'hg19'.
        
    OUTPUT: reads_data_dictionary - Dictionary of the form: Key = header; Item = (reference_genome, chr, strand, integration_locus, read end, span, lam_id)
    
    lOGICS: Given a connection object (it already contains DB\DB_schema target), this function perform queries to return reads data in form of a dictionary
            (details explained in 'OUTPUT' section).
            
    !! WARNINGS !! : 1) This function perform a query like the following:
                        "SELECT `header`, `chr`, `strand`, `integration_locus`, 100 as `span`, `complete_name` as lam_id  FROM [...]"
                        PLEASE NOTE *** 100 as `span` ***, necessary to avoid errors when span is NULL: on the other side, this fact could generate
                        WRONG VALUES for READ END element in reads_data_dictionary tuple items
                        
                     2) 'read end' is not retrieved but calculated as 'start + span': conversely 'start' and 'span' are directly retrieved from DB
    """    
    # Initialize reads_data_dictionary to collect results
    reads_data_dictionary={}
    
    # Create dictionary cursor
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    
    # Query for Reads Data and cursor closing
    cursor.execute("SELECT `header`, `chr`, `strand`, `integration_locus`, 100 as `span`, `complete_name` as lam_id  FROM {0} WHERE 1".format(db_table)) # 100 as `span` to prevent errors due to "NULL" span
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
    
    # Return result
    return reads_data_dictionary
    
def import_lam_data_from_DB (conn, db_table, reference_genome="unspecified genome"):
    """
    *** Get LAM data in the form of Dictionary, directly from DB ***
    
    INPUT: conn
           db_table
    
    OPTIONAL INPUT: reference_genome - String. Something like 'hg19'.
        
    OUTPUT: lam_data_dictionary - Dictionary of the form: Key = lam_id; Item = (n_LAM, tag, pool, tissue, sample, treatment, group_name, enzyme, vector)
    
    lOGICS: Given a connection object (it already contains DB\DB_schema target), this function perform queries to return LAM data in form of a dictionary
            (details explained in 'OUTPUT' section).
            
    !! WARNINGS !! :  This function perform a query like the following:
                      "SELECT DISTINCT `complete_name` as lam_id, [...] FROM [...]"
                      PLEASE NOTE *** `complete_name` as lam_id *** because up to now 'lam_id' data are labeled as 'complete_name' in DB
    """
    #Initialize lam_data_dictionay to collect results ('lam_data_dictionay' -> return)
    lam_data_dictionay={}
    
    # Create dictionary cursor
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()
    # Query for Lam Data
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
        lam_data_dictionay[dat['lam_id']]=(dat['n_LAM'], dat['tag'], dat['pool'], dat['tissue'], dat['sample'], str(dat['treatment']).zfill(len_max), dat['group_name'], dat['enzyme'], dat['vector'])
           
    # Return result
    return lam_data_dictionay