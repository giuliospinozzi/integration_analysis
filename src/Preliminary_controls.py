'''
Created on 26/set/2013

@author: Stefano
'''

###Requested Package(s) Import#
import MySQLdb
###############################

###Import Module(s)#
import DB_connection
####################





def check_syntax (args_dbDataset, args_collision, check, reason):
    '''
    *** This function controls for common syntax mistake ***
    
    INPUT: args_dbDataset - user input, a string such as 'dbschema.dbtable' for one only, 'dbschema1.dbtable1,dbschema2.dbtable2,dbschema3.dbtable3' for three (args.dbDataset)
           args_collision - user input Boolean (args.collision)
           check - Boolean
           reason - String
           
    OUTPUT: check - Boolean, set 'False' only if input variables don't pass controls, otherwise left as given in input
            reason - String, modified only if input variables don't pass controls, otherwise left as given in input
            
    LOGIC: if 'check' is given True this function controls "dbDataset" string structure first, then collision feasibility,
           switching 'check' to False and explaining why in 'reason', if necessary.
    '''
    if (check == True):
        dbDataset_split = args_dbDataset.split(",")
        
        # Check for errors in dbDataset choice syntax (args_dbDataset), such as commas and quotes...
        if (("" or "'" or '''"''') in dbDataset_split):
            check = False
            reason = "check syntax in --dbDataset argument, there should be a mistake using commas or quotes"
            return check, reason
        
        # Check for errors in dbDataset names
        for db_string in dbDataset_split:
            db_split = db_string.split(".")
            if (len(db_split)!=2):
                check = False
                reason = "check syntax in --dbDataset argument, there should be a mistake related to dots (can't recognize dbschema.dbtable structure)"
                return check, reason
                
        # Check feasibility of collision request     
        if ((args_collision == True) and (len(dbDataset_split)<2)):
            check = False
            reason = "can't perform collision with only one input dataset (see --dbDataset argument)"
            return check, reason
        
    return check, reason





def check_DB_for_data (host, user, passwd, port, args_dbDataset, check, reason):
    '''
    *** This function controls if desired DB schema(s) and DB table(s) are available at selected Host ***
    
    INPUT: host, user, passwd, port - user input to set up DB connection (args.host, args.user, args.pw, args.dbport)
           args_dbDataset - user input, a string such as 'dbschema.dbtable' for one only, 'dbschema1.dbtable1,dbschema2.dbtable2,dbschema3.dbtable3' for three (args.dbDataset)
           check - Boolean
           reason - String
           
    OUTPUT: check - Boolean, set 'False' only if input variables don't pass controls, otherwise left as given in input
            reason - String, modified only if input variables don't pass controls, otherwise left as given in input
            
    LOGIC: if 'check' is given True this function asks to host for table_schema(s) (dbschema) and table_name(s) (dbtable) availability,
           switching 'check' to False and explaining why in 'reason', if necessary.
           
    
    WARNING: to test (see n = cursor.fetchall()[0][0], dunno if it's ok!!)
    '''
    if (check == True):
        
        # Preparing data for queries
        dbDataset_tuple_list = [] # [('dbschema1', 'dbtable1'), ('dbschema2', 'dbtable2'), ...]
        dbDataset_split = args_dbDataset.split(",")
        for db_string in dbDataset_split:
            db_tupla = None
            db_split = db_string.split(".")
            db_tupla = (db_split[0],db_split[1])
            dbDataset_tuple_list.append(db_tupla)
        
        # Loop for queries
        for db_tupla in dbDataset_tuple_list:
            conn = DB_connection.dbOpenConnection (host, user, passwd, port, db_split[0], db_split[1])
            cursor = conn.cursor (MySQLdb.cursors.Cursor)
            cursor.execute ("SELECT count(*) FROM information_schema.tables WHERE table_schema = '{0}' AND table_name = '{1}'".format(db_split[0], db_split[1]))
            n = cursor.fetchall()[0][0] #...una prova! bisogna lanciare
            if (int(n)==0):
                check = False
                reason = "verify --dbDataset argument: can't find db_schema = '{0}' and db_table = '{1}' for user '{2}' on host '{3}'".format(db_split[0], db_split[1], user, host)
                cursor.close()
                DB_connection.dbCloseConnection(conn)
                return check, reason        
            # Close cursor and connection
            cursor.close()
            DB_connection.dbCloseConnection(conn)
        
    return check, reason
                
                


###





def check_columnsToGroup (args_columnsToGroup, args_columns, check, reason):
    '''
    *** This function controls if "columnsToGroup" choice is compatible with selected "columns" ***
    
    INPUT: args_columnsToGroup - user input, a string such as 'sample,tissue,treatment' (args.columnsToGroup)
           args_columns - user input, a string such as 'tissue,treatment' (args.columns)
           check - Boolean
           reason - String
           
    OUTPUT: check - Boolean, set 'False' only if input variables don't pass controls, otherwise left as given in input
            reason - String, modified only if input variables don't pass controls, otherwise left as given in input
            
    LOGIC: if 'check' is given True and 'args_columnsToGroup' is not None, this function controls if "columnsToGroup" categories (args_columnsToGroup) are a sub-set of "columns" categories (args_columns),
           switching 'check' to False and explaining why in 'reason', if necessary.
        
    '''
    if (check == True):
        
        if (args_columnsToGroup != None):
            selected_category = args_columns.split(",")
            merge_over = args_columnsToGroup.split(",")
            
            for word in merge_over:
                
                if (word not in selected_category):
                    check = False
                    reason =  "each category given as --columnsToGroup argument must be given as --columns argument too. Your input: Columns-> {0}; ColumnsToGroup-> {1}".format(selected_category,merge_over)
                    
    return check, reason




def check_method (IS_method, bushamn_bp_rule, IS_methods_list, check, reason):
    '''
    INPUT:  IS_method - user input, a string such as 'classic', reflecting user choice of IS retrieving method (args.IS_method)
            bushamn_bp_rule - int number (args.bushman_bp_rule, suddenly put in bushamn_bp_rule in order to modify 
            IS_method_list - a list of strings, such as ['classic', 'whatever', ... ], collecting all available IS retrieving methods
            args_columns - user input, a string such as 'sample,tissue,treatment' (args.columns)
    
    OUTPUT:
    
    LOGIC: if 'check' is given True this function perform queries to DB(s) (args_dbDataset) to check availability of data desired by user (args_columns),
           switching 'check' to False and explaining why in 'reason', if necessary.
    '''

    
    
    
    
    
    
    
    
    
    
    
    