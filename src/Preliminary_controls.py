'''
Created on 26/set/2013

@author: Stefano
'''

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
           switching 'check' in False and explaining why in 'reason', if necessary.
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
           switching 'check' in False and explaining why in 'reason', if necessary.
        
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