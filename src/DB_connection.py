###Header###############################################
header = """

+------------------------------------------------------+
 Module: DB_connection_light
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
    

###Import input data from DB###############################################################################################################
def import_data_from_DB (host, user, passwd, db, db_table, reference_genome="unspecified genome"):
    '''
    [...]
    '''
     
    ###Setting Up Connection to DB ########################
    #(Retrieved from input variables)######################
    userdb = {
     'host': host,
     'user': user,
     'passwd': passwd,
     'db': db
     }
    
    conn = MySQLdb.connect( 
     host = userdb["host"],
     user = userdb["user"],
     passwd = userdb["passwd"],
     db = userdb["db"]
     )
    #######################################################
    
    
    ###QUERIES######################################################################################################
    
    #Open DB Connection###
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    
    #Query for Reads Data###
    cursor.execute("SELECT `header`, `chr`, `strand`, `integration_locus`, `span`, `lam_id`  FROM {0} WHERE 1".format(db_table))
    reads_data = cursor.fetchall()
    cursor.close()
    
    #Build reads data dictionary ('reads_query' -> return)
    reads_query={}
    for dat in reads_data:
        reads_query[dat['header']]=(reference_genome, dat['chr'], dat['strand'], dat['integration_locus'], dat['integration_locus'] + dat['span'], dat['span'], dat['lam_id'])
    del reads_data
    
    #Query for Lam Data###
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    cursor.execute("SELECT DISTINCT `lam_id`,`n_LAM`, `tag`, `pool`, `tissue`, `sample`, `treatment`, `group_name`, `enzyme`  FROM {0} WHERE 1".format(db_table))
    lam_data = cursor.fetchall()
    cursor.close()
    
    #Close DB Connection###
    conn.close()
    
    #Build lam data dictionary ('lam_query' -> return)
    lam_query={}
    for dat in lam_data:
        if (len(dat['treatment'])==1):
            dat.update(treatment="0{0}".format(dat['treatment']))
        lam_query[dat['lam_id']]=(dat['n_LAM'], dat['tag'], dat['pool'], dat['tissue'], dat['sample'], dat['treatment'], dat['group_name'], dat['enzyme'])
    del lam_data
    #################################################################################################################

        
    ###Return results##################################
    return reads_query, lam_query
    ###################################################
    
#################################################################################################################################################


###Import extra-columns for final matrix 'sample_tissue_treatment'###############################################################################

def get_extra_columns_from_DB (host, user, passwd, db, db_table, parameters_list, query_for_columns, reference_genome):
    '''
    [...]
    '''
     
    ###Setting Up Connection to DB ########################
    #(Retrieved from input variables)######################
    userdb = {
     'host': host,
     'user': user,
     'passwd': passwd,
     'db': db
     }
    
    conn = MySQLdb.connect( 
     host = userdb["host"],
     user = userdb["user"],
     passwd = userdb["passwd"],
     db = userdb["db"]
     )
    #######################################################

    ###QUERIES######################################################################################################
    
    #Open DB Connection###
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    
    #Query for column labels###
    cursor.execute("SELECT DISTINCT {0} FROM {1} WHERE 1".format(query_for_columns, db_table))
    column_labels = cursor.fetchall()
    cursor.close()
       
    #Build column labels list, ordered ('column_labels_list' -> return)
    column_labels_list=[]
    
    for dat in column_labels:
        if ('treatment' in parameters_list):
            if (len(dat['treatment'])==1):
                dat.update(treatment="0{0}".format(dat['treatment']))
        #create label
        label = ""
        for parameter in parameters_list:
            label = label + "{0}_".format(dat[parameter])
        label = label[:-1]
        #append label to column_labels_list
        column_labels_list.append(label)
        
    del column_labels
    column_labels_list.sort()
    
    #Query for merged column labels '_tissue_treatment' ### IF POSSIBLE
    merged_column_labels_list=[]
    if ((len(parameters_list) >= 3) and ('tissue' in parameters_list) and ('treatment' in parameters_list)):
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)
        cursor.execute("SELECT DISTINCT `tissue` , `treatment` FROM {0} WHERE 1".format(db_table))
        merged_column_labels = cursor.fetchall()
        cursor.close()
        
        #Build merged column labels list, ordered ('merged_column_labels_list' -> return)
        for dat in merged_column_labels:
            if (len(dat['treatment'])==1):
                dat.update(treatment="0{0}".format(dat['treatment']))
            merged_column_labels_list.append("_{0}_{1}".format(dat['tissue'], dat['treatment']))
        del merged_column_labels
        merged_column_labels_list.sort()
    #################################################################################################################

    #Close DB Connection###
    conn.close()

    ###Return results##################################
    return column_labels_list, merged_column_labels_list
    ###################################################

#################################################################################################################################################




