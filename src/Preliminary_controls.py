###Header################################################
header = """

+------------------------------------------------------+
 Module: Preliminary_controls
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  September 27th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - This module contains functions to control launch 
    command and input arguments. They are all collected
    in a unique function (smart_check) suitable for 
    typical usage
  
 Note:
  - check_method should be updated (inside its code and
    through IS_methods_list variable, see IS method tuning
    box in Integration_Analysis.py file)

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import MySQLdb
###############################

###Import Module(s)#
import DB_connection
####################





def smart_check (args_dbDataset, args_collision, host, user, passwd, port, args_columns, args_columnsToGroup, IS_method, bushamn_bp_rule, IS_methods_list, check, reason):
    '''
    *** This function controls for user's input ***
    
    LOGIC: This function collects all the following ones in the typical 'usage in chain'.
           See docs of each function called inside for details.
           
    INPUT: check - Boolean, HAS TO BE GIVEN AS 'TRUE' IN INPUT in order to work properly and performing controls
           reason - String.
                    Choose wisely: this string will be the returned as-it-is in case of 
                    any control will fail (check becomes 'False') but an explanation of the reason why
                    is not foreseen ('unexpected error. Try to check syntax and DB connection availability.'
                    fits the case).
                    Likewise, this string will be returned also if 'check' is given 'False' in input and any
                    control will be skipped.
                    
    OUTPUT: check - Boolean, set 'False' only if input variables don't pass controls, otherwise left as given in input (typically 'True')
            reason - String, modified only if input variables don't pass controls and a specific reason is foreseen, otherwise left as given in input   
    '''

    check, reason = check_syntax (args_dbDataset, args_collision, check, reason)
    check, reason = check_DB_for_data (host, user, passwd, port, args_dbDataset, check, reason)
    check, reason = check_DB_for_columns (host, user, passwd, port, args_dbDataset, args_columns, check, reason)
    check, reason = check_columnsToGroup (args_columnsToGroup, args_columns, check, reason)
    check, reason = check_method (IS_method, bushamn_bp_rule, IS_methods_list, check, reason)
    
    return check, reason






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
                        
        #Check for repeated/duplicated dataset (problems during collisions)
        n_dataset = len(dbDataset_split)
        n_unique_dataset = len(set(dbDataset_split))
        if (n_dataset != n_unique_dataset):
            check = False
            reason = "check syntax in --dbDataset argument: you give the same dataset as argument twice."
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
        
        # Open Connection to DB
        conn = MySQLdb.connect(host = host, user = user, passwd = passwd, port = port)
        
        # Loop for queries
        for db_tupla in dbDataset_tuple_list:
            
            # Check for existence
            cursor = conn.cursor (MySQLdb.cursors.Cursor)
            cursor.execute ("SELECT count(*) FROM information_schema.tables WHERE table_schema = '{0}' AND table_name = '{1}'".format(db_tupla[0], db_tupla[1]))
            n = cursor.fetchall()[0][0]

            if (int(n)==0):
                check = False
                reason = "can't find db_schema = '{0}' and db_table = '{1}' for user '{2}' on host '{3}', please verify --dbDataset argument".format(db_tupla[0], db_tupla[1], user, host)
                cursor.close()
                DB_connection.dbCloseConnection(conn)
                return check, reason
            
            cursor.close()
            
            # Check for data inside
            cursor = conn.cursor (MySQLdb.cursors.Cursor)
            cursor.execute ("SELECT count(*) FROM {0}.{1} WHERE 1".format(db_tupla[0], db_tupla[1]))
            n_row = cursor.fetchall()[0][0]
            
            if (int(n_row)==0):
                check = False
                reason = "[db_schema = '{0}', db_table = '{1}'] exists but is EMPTY on host '{2}'".format(db_tupla[0], db_tupla[1], host)
                cursor.close()
                DB_connection.dbCloseConnection(conn)
                return check, reason
            
            cursor.close()
            
        # Close Connection to DB    
        DB_connection.dbCloseConnection(conn)
        
    return check, reason
                
                



def check_DB_for_columns (host, user, passwd, port, args_dbDataset, args_columns, check, reason):
    '''
    *** This function controls if categories (columns) chosen by user are effectively in each dbschema.dbtable couple ***
    
    INPUT: host, user, passwd, port - user input to set up DB connection (args.host, args.user, args.pw, args.dbport)
           args_dbDataset - user input, a string such as 'dbschema.dbtable' for one only, 'dbschema1.dbtable1,dbschema2.dbtable2,dbschema3.dbtable3' for three (args.dbDataset)
           args_columns - user input, a string such as 'sample,tissue,treatment' (args.columns)
           check - Boolean
           reason - String
           
    OUTPUT: check - Boolean, set 'False' only if input variables don't pass controls, otherwise left as given in input
            reason - String, modified only if input variables don't pass controls, otherwise left as given in input
            
    LOGIC: if 'check' is given True this function asks to host if categories desired by user (--columns argument, args_columns here) are available in each given 
           'dbschema.dbtable' couple (--dbDataset argument, args_dbDataset here), switching 'check' to False and explaining why in 'reason', if necessary.    
    '''
    if (check == True):
        
        # Preparing selected_categories list
        selected_categories = args_columns.split(",")
        
        # Preparing problems_list
        problems_list =[]
        
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
            #conn = DB_connection.dbOpenConnection (host, user, passwd, port, db_tupla[0], db_tupla[1])
            conn = DB_connection.dbOpenConnection (host, user, passwd, port, db_tupla[0])
            cursor = conn.cursor (MySQLdb.cursors.Cursor)
            cursor.execute ("show columns from {1} from {0}".format(db_tupla[0], db_tupla[1]))
            data_retrieved = cursor.fetchall() 
            
            #Storing available categories in retrieved_categories
            retrieved_categories =[]
            for item in data_retrieved:
                retrieved_categories.append(item[0])
            
            #Check for matching between retrieved_categories and selected_categories
            for column in selected_categories:
                if (column not in retrieved_categories):
                    check = False
                    problems_list.append("'{0}' in '{1} - {2}'".format(column, db_tupla[0], db_tupla[1]))
        
            # Close cursor and connection
            cursor.close()
            DB_connection.dbCloseConnection(conn)
            
        # Format reason
        if (check == False):
            problems_string = ", ".join(problems_list)
            reason = "can't find " + problems_string + ". Please verify you are querying the correct DB and/or --columns argument"
        
    return check, reason





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
    *** This function controls if "IS_method" user's choice is available ***
    
    INPUT:  IS_method - user input, a string such as 'classic', reflecting user choice of IS retrieving method (args.IS_method)
            bushamn_bp_rule - int number (args.bushman_bp_rule, suddenly put in bushamn_bp_rule in order to modify 
            IS_method_list - a list of strings, such as ['classic', 'whatever', ... ], collecting all available IS retrieving methods
            check - Boolean
            reason - String
            
    OUTPUT: check - Boolean, set 'False' only if input variables don't pass controls, otherwise left as given in input
            reason - String, modified only if input variables don't pass controls, otherwise left as given in input
            
    LOGIC: if 'check' is given True, this function controls if IS retrieving method selected by user (IS_method) exists (IS_method_list),
           switching 'check' to False and explaining why in 'reason', if necessary.
           If user also specified a bushamn_bp_rule by hand, besides a method having its mandatory default, it produce a *warning* informing
           user that his bushamn_bp_rule choice will be ignored.
           
    WARNING: this function has to be updated every time a new IS retrieving method will have been added or bushamn_bp_rule standard (now 6)
             will be changed. Check it!
    '''
    if (check == True):
        
        if (IS_method not in IS_methods_list):
            check = False
            reason = "selected IS retrieving methods doesn't exist: your input -> '{0}'. Please check --IS_method argument and retry with an available one, within {1}  .".format(IS_method, str(IS_methods_list))
            return check, reason
        
        else:
            if ((IS_method != "classic") and (bushamn_bp_rule != 6)):
                check = True
                print "\n\t  *WARNING*\n\t  *{0} method has its default for bushamn_bp_rule, that is '{1}'*\n\t  *Your bushamn_bp_rule setting will be ignored!!!*\n".format(IS_method, str(bushamn_bp_rule))
                
    return check, reason       
    