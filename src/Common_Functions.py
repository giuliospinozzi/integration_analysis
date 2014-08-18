###Header################################################
header = """

+------------------------------------------------------+
 Module: Common_functions
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 22th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.2
+------------------------------------------------------+

 Description:
  - [...] Set of functions of common usage
  
 Note:
  - Description still to do!!

-------------------------------------------------------- 
""" 
########################################################

### temp mod ###
import copy
################

###For parsing#####################################
def prepareSELECT(columnsToGroup):
    #select_temp = columnsToGroup[1:-1].split(",")
    select_temp = columnsToGroup.split(",")
    select = "`" + "`,`".join(str(x) for x in select_temp) + "`" 
    #select = ""
    #for s in select_temp:
    #    select = select + "`" + s + "`" + ","
    #select = select[:-1]
    return select
####################################################



###Given a key for read dictionary, it returns related lam data###################
def get_lam (reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay):
    lam_data_dictionay_Key = reads_data_dictionary[reads_data_dictionary_Key][-1]
    lam_data = lam_data_dictionay[lam_data_dictionay_Key]
    return lam_data #tuple
##################################################################################



###Find header of the longest read and related sequences ###############################
def find_longest_read (list_of_reads_key, raw_read_dictionary, final_read_dictionary):
    ### temp mod ###
    old_list_of_reads_key = copy.deepcopy(list_of_reads_key)
    list_of_reads_key = []
    for old_header in old_list_of_reads_key:
        if "/1" in old_header:
            header = old_header.replace("/1", "")
        elif "/2" in old_header:
            header = old_header.replace("/2", "")
        list_of_reads_key.append(header)
    ################
    selected_header=list_of_reads_key[0]
    for header in list_of_reads_key[1:]:
        if (len(raw_read_dictionary[header])>len(raw_read_dictionary[selected_header])):
            selected_header = header
    longest_raw_sequence = raw_read_dictionary[selected_header]
    longest_final_sequence = final_read_dictionary[selected_header]
    return selected_header, longest_raw_sequence, longest_final_sequence

########################################################################################



###From my matrix (the ones returned by functions in Matrix_creation module - somewhat_matrix_as_line_list) to user's matrix (returns matrix_as_line_list) ##############
def convert_matrix (matrix_as_line_list, user_label_dictionary, user_merged_labels_dictionary):
    #no problem if user_merged_labels_dictionary is void!
    
    #switch to user labels changing first matrix line, preparing dictionary_of_columns and its bunch of key (list_of_columns_keys)
    dictionary_of_columns = {} #key: user_label (new column label), item: list of value for that column
    list_of_columns_keys =[] #list of keys for dictionary_of_columns
    first_line_split = matrix_as_line_list[0].split("\t")
    i=0
    for cell in first_line_split:
        if (user_label_dictionary.has_key(cell)):
            first_line_split[i] = user_label_dictionary[cell][0]
            dictionary_of_columns.update({user_label_dictionary[cell][0]:[]})
            list_of_columns_keys.append(user_label_dictionary[cell][0])
        i+=1
    matrix_as_line_list[0] = "\t".join(first_line_split)
    
    #Fill dictionary of columns
    j=0 #used below, to know list length
    for line in matrix_as_line_list[1:]:
        j+=1
        line_split = line.split("\t")
        line_split = line_split[3:-1]
        i=3
        for cell in line_split:
            dictionary_of_columns[first_line_split[i]].append(cell)
            i+=1
    
    #Creating merged columns and updating dictionary of columns
    void_list = [0]*j
    for key in user_merged_labels_dictionary.keys():
        # dictionary_of_columns.update({key:void_list})
        dictionary_of_columns[key] = void_list
        i=0
        for item in user_merged_labels_dictionary[key]:
            i+=1
            tmplist = [int(k) for k in dictionary_of_columns[item]]
            dictionary_of_columns[key] = map(lambda t,z:t+z, dictionary_of_columns[key], tmplist)
                        
    #Updating list_of_columns_keys and sort
        
    #list_of_columns_keys = list_of_columns_keys + user_merged_labels_dictionary.keys()
    #list_of_columns_keys.sort()
    #This (deprecated) way of sorting shows some trouble because capital cases < "_"
    #while lower case > "_". This may cause a wrong columns placement among 'straight'
    #and 'merged' ones, in matrixes.
            
    list_of_merged_column_keys = user_merged_labels_dictionary.keys()
    list_of_columns_keys.sort()
    list_of_merged_column_keys.sort()
    list_of_columns_keys.extend(list_of_merged_column_keys)
    
    #Switch to the new matrix_as_line_list
    first_line_list = first_line_split[:3] + list_of_columns_keys + [first_line_split[-1]]
    matrix_as_line_list[0] = "\t".join(first_line_list)
    for i in range(len(matrix_as_line_list[1:])):
        line_split = matrix_as_line_list[i+1].split("\t")
        current_line = "\t".join(line_split[:3])
        for label in list_of_columns_keys:
            current_line = current_line + "\t" + str(dictionary_of_columns[label][i])
        current_line = current_line + "\t" + line_split[-1]
        matrix_as_line_list[i+1] = current_line
    
        
    #Return results###########    
    return matrix_as_line_list
    ##########################
        
#########################################################################################################################################################################

