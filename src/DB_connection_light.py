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


def import_data_from_DB (host, user, passwd, db, db_table, reference_genome="unspecified genome"):
    '''
    [...]
    '''

    ###Requested Package(s) Import###
    import MySQLdb
    #################################
     
    ###Setting Up Connection to DB ########################
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

    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    cursor.execute("SELECT `header`, `chr`, `strand`, `integration_locus`, `span`, `lam_id`  FROM {0} WHERE 1".format(db_table))
    reads_data = cursor.fetchall()
    cursor.close()
    
    reads_query={}
    for dat in reads_data:
        reads_query[dat['header']]=(reference_genome, dat['chr'], dat['strand'], dat['integration_locus'], dat['integration_locus'] + dat['span'], dat['span'], dat['lam_id'])
    del reads_data
    
    cursor = conn.cursor (MySQLdb.cursors.DictCursor)
    cursor.execute("SELECT `lam_id`,`n_LAM`, `tag`, `pool`, `tissue`, `sample`, `treatment`, `group_name`, `enzyme`  FROM {0} WHERE 1".format(db_table))
    lam_data = cursor.fetchall()
    cursor.close()
    
    conn.close()
    
    lam_query={}
    for dat in lam_data:
        lam_query[dat['lam_id']]=(dat['n_LAM'], dat['tag'], dat['pool'], dat['tissue'], dat['sample'], dat['treatment'], dat['group_name'], dat['enzyme'])
    del lam_data
    
    #################################################################################################################

         
    #===========================================================================
    # #DEVELOPMENT
    # #Remove comment-block to print the 10 former element in "query" dictionary, as control.
    # i=0
    # print "**********************"
    # print "Print for development:\n"
    # print "Reads Dictionary"
    # for key, value in reads_query.iteritems():
    #     i+=1
    #     print key, value
    #     if (i>=10):
    #         break
    # i=0
    # print "\nLAM Dictionary"
    # for key, value in lam_query.iteritems():
    #     i+=1
    #     print key, value
    #     if (i>=10):
    #         break
    # print "**********************\n"
    #===========================================================================
        
    ###Return results by means of "query" dictionary###
    return reads_query, lam_query
    ###################################################

    
