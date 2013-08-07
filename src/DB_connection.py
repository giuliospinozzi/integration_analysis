###Header###############################################
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
    

###Import input data from DB###############################################################################################################
def import_data_from_DB (host, user, passwd, db, db_table, query_step=1000000, reference_genome="unspecified genome"):
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
    
    #reads_query dictionary to collect results
    reads_query={}
    
    #splitting query to reduce memory usage peak
    cursor = conn.cursor (MySQLdb.cursors.Cursor) #  
    cursor.execute ("SELECT count( * ) FROM {0} WHERE 1".format(db_table))
    n_rows = cursor.fetchall()[0][0]
    cursor.close()
    splitting = range(0, n_rows, query_step) #by default one million row a time
    
    for n in splitting:
            
        #Open DB Connection###
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)
        
        #Query for Reads Data###
        start = str(n)
        end = str((n+query_step))
        cursor.execute("SELECT `header`, `chr`, `strand`, `integration_locus`, `span`, `lam_id`  FROM {0} WHERE 1 LIMIT {1}, {2}".format(db_table, start, end))
        reads_data = cursor.fetchall()
        cursor.close()
        
        #Build reads data dictionary ('reads_query' -> return)
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
    
    #User_label_template
    user_label_template = parameters_list
    
    #user_label_dictionary
    user_label_dictionary = {} #of kind: key-> my label, item -> [user label, user label as tupla, my label as tupla]
       
    #Build column labels list, ordered ('column_labels_list' -> return)
    column_labels_list=[]
    
    for dat in column_labels:
        if ('treatment' in parameters_list):
            if (len(dat['treatment'])==1):
                dat.update(treatment="0{0}".format(dat['treatment']))
        #create label #labels building has to keep freezed like this because they need to match with labels builded in "Classes_for_Integration_Analysis" module (class Covered_base)
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
        #append label to column_labels_list
        column_labels_list.append(label)
        
        #create user labels
        user_label = ""
        user_label_as_tupla = ()
        for category in user_label_template:
            user_label = user_label + "_" + dat[category]
            user_label_as_tupla = user_label_as_tupla + (dat[category],)
        user_label = user_label[1:]
        
        #update user_label_dictionary
        user_label_dictionary.update({label:[user_label, user_label_as_tupla, label_as_tupla]})
            
            
        
    del column_labels
    column_labels_list.sort()
    
    ###POSSIBLE IMPROVEMENTS HERE###
    #Here you can implement code to produce "default column merging", directly in "redundant" ad "IS" output
    #This is a template: if user perform a query (--columns argument) such as 'sample, tissue, treatment', output files show also "merged columns"
    #with labels such as '_tissue_treatment' (sum over / merged over sample)
    #By default this option is commented and merged_column_labels_list is returned void; if uncommented, all outputs will be created correctly automatically
    #TO MAKE POSSIBLE MORE KIND OF "MERGED COLUMNS" --> See Classes_for_Integration_Analysis for merged_column_labels_list.append("_{0}_{1}..... labels building)
    merged_column_labels_list=[]
    #===========================================================================
    # if ((len(parameters_list) >= 3) and ('tissue' in parameters_list) and ('treatment' in parameters_list)):
    #     cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    #     cursor.execute("SELECT DISTINCT `tissue` , `treatment` FROM {0} WHERE 1".format(db_table))
    #     merged_column_labels = cursor.fetchall()
    #     cursor.close()
    #     
    #     #Build merged column labels list, ordered ('merged_column_labels_list' -> return)
    #     for dat in merged_column_labels:
    #         if (len(dat['treatment'])==1):
    #             dat.update(treatment="0{0}".format(dat['treatment']))
    #         merged_column_labels_list.append("_{0}_{1}".format(dat['tissue'], dat['treatment']))
    #     del merged_column_labels
    #     merged_column_labels_list.sort()
    #===========================================================================
    #################################################################################################################

    #Close DB Connection###
    conn.close()

    ###Return results##########################################################
    return column_labels_list, merged_column_labels_list, user_label_dictionary
    ###########################################################################

#################################################################################################################################################




