###Header################################################
header = """

+------------------------------------------------------+
 Module: output_module
 Author: Stefano Brasca
 Date:  January 8th, 2014
 Contact: brasca.stefano@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - This module contains functions to create output in 
    *.tsv format and in 'Excel Workbook' format (*.xlsx)
    
 Requirements:
  - XlsxWriter, a Python module for creating Excel XLSX
    files (see http://xlsxwriter.readthedocs.org/ )
  
 Note:
  - Developed under XlsxWriter version 0.5.2

-------------------------------------------------------- 
""" 
#########################################################

###Requested Package(s) Import#
import xlsxwriter
###############################





#########################################################
def tsv_output (matrix_file_name, matrix_as_line_list):
    
    #Open a file named as 'redundant_matrix_file_name'
    file_output = open(matrix_file_name, 'w')
    
    #Fill file line by line
    for line in matrix_as_line_list:
        file_output.write(line)
    
    #Close file    
    file_output.close()
#########################################################

    


###########################################################################################################################################################################
def workbook_output (result_dictionary, mode = 'feature_rich'): # or mode = 'basic'
    '''
    *** This function generates an output summary file of kind 'Excel Workbook' ***
    
    INPUT: - result_dictionary: a dictionary like:  {
                                                    'dataset_name':dbschema.dbtable,
                                                    'redundant_matrix':redundant_matrix_as_line_list,
                                                    'IS_matrix':IS_matrix_as_line_list
                                                    'IS_matrix_collided':IS_matrix_as_line_list_collided / None
                                                    'list_of_Covered_Bases':list_of_Covered_Bases
                                                    'list_of_Covered_bases_ensambles':list_of_Covered_bases_ensambles
                                                    'IS_list':IS_list
                                                    'IS_method': IS_method
                                                    'strand_specific_choice':strand_specific_choice
                                                    }
                                                    
    [...]                                                
                                                    
    '''
    
    ### CREATE WORKBOOK ###########################################################################
    
    # Create Workbook name
    file_name_part = result_dictionary['dataset_name'].replace(".", "_")
    workbook_file_name = "Integration_Analysis_" +  file_name_part + ".xlsx"
     
    # Create Workbook instance and set policy
    workbook_output = xlsxwriter.Workbook(workbook_file_name,
        {'in_memory': True, # Memory usage policy: True means 'work fast at the expense of memory usage'
        'strings_to_formulas': False, # String conversion policy: False means 'Don't try to auto-convert strings to formulas'
        'strings_to_urls': False, # String conversion policy: False means 'Don't try to auto-convert strings to URL links'
        'strings_to_numbers': True}) # String conversion policy: False means 'Don't try to auto-convert strings to numbers'
     
    # Set Workbook metadata
    title = 'Integration Analysis'
    dataset = result_dictionary['dataset_name'].replace(".", " - ")
    author = 'Stefano Brasca'
    manager = 'Eugenio Montini'
    company = 'TIGET - Safety of Gene Therapy and Insertional Mutagenesis Research Unit'
    comments = '''Created by Montini's Bioinfo Team: Andrea Calabria - calabria.andrea@hsr.it; Stefano Brasca - brasca.stefano@hsr.it; Giulio Spinozzi - spinozzi.giulio@hsr.it'''
    workbook_output.set_properties({
        'title':    title,
        'subject':  dataset,
        'author':   author,
        'manager':  manager,
        'company':  company,
        'category': '',
        'keywords': '',
        'comments': comments})
    
    # Initialize Variables for mode = 'feature_rich'
    standard_data_rows_indexes = None 
    standard_data_columns_indexes = None 
    merged_data_columns_indexes = None
    collision_data_columns_indexes = None
            
    
    ### REDUNDANT WORKSHEET #######################################################################
    
    # Create Worksheet name    #Note: must be less than 32char
    redundant_worksheet_name = "RedundantReads"
    if (result_dictionary['strand_specific_choice'] == True):
        redundant_worksheet_name = redundant_worksheet_name + "_StrandSpecific"
    
    # Create Worksheet instance
    redundant_worksheet = workbook_output.add_worksheet(redundant_worksheet_name)
    
    # Case of mode = 'feature_rich':
    if (mode == 'feature_rich'):
        standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes = find_indexes(result_dictionary['redundant_matrix'])
    
    # Fill Worksheet with Redundant Reads data
    write_matrix_in_worksheet(workbook_output, redundant_worksheet, result_dictionary['redundant_matrix'], mode, standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes)
    
        
    ### IS WORKSHEET ##############################################################################
    
    # Create Worksheet name    #Note: must be less than 32char
    IS_worksheet_name = "ISmatrix" + "_{0}_".format(result_dictionary['IS_method'])
    if (result_dictionary['IS_matrix_collided'] != None):
        IS_worksheet_name = IS_worksheet_name + "&collisions"
    
    # Create Worksheet instance
    IS_worksheet = workbook_output.add_worksheet(IS_worksheet_name)
    
    # Case of mode = 'feature_rich':
    if (mode == 'feature_rich'):
        standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes = find_indexes(result_dictionary['IS_matrix'], result_dictionary['IS_matrix_collided'])
    
    # Select IS data to use
    selected_matrix_as_line_list = None
    if (result_dictionary['IS_matrix_collided'] != None):
        selected_matrix_as_line_list = result_dictionary['IS_matrix_collided']
    else:
        selected_matrix_as_line_list = result_dictionary['IS_matrix']
            
    # Fill Worksheet with IS data
    write_matrix_in_worksheet(workbook_output, IS_worksheet, selected_matrix_as_line_list, mode, standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes)
    
    
    ### CONCLUSIVE ACTIONS ########################################################################
        
    # Closing
    workbook_output.close()
        
    # Debug find_indexes function below
    
    #===========================================================================
    # #Create log file
    # file_output = open('log.txt', 'a')
    # file_output.write('Dataset: ' + result_dictionary['dataset_name'])
    # 
    # #Redundant
    # standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes = find_indexes(result_dictionary['redundant_matrix'])
    # file_output.write('\n\nRedundant Indexes:')
    # file_output.write('\nstandard_data_rows_indexes = ' + str(standard_data_rows_indexes))
    # file_output.write('\nstandard_data_columns_indexes = ' + str(standard_data_columns_indexes))
    # file_output.write('\nmerged_data_columns_indexes = ' + str(merged_data_columns_indexes))
    # file_output.write('\ncollision_data_columns_indexes = ' + str(collision_data_columns_indexes))
    # 
    # #IS
    # standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes = find_indexes(result_dictionary['IS_matrix'], result_dictionary['IS_matrix_collided'])
    # file_output.write('\n\nIS Indexes:')
    # file_output.write('\nstandard_data_rows_indexes = ' + str(standard_data_rows_indexes))
    # file_output.write('\nstandard_data_columns_indexes = ' + str(standard_data_columns_indexes))
    # file_output.write('\nmerged_data_columns_indexes = ' + str(merged_data_columns_indexes))
    # file_output.write('\ncollision_data_columns_indexes = ' + str(collision_data_columns_indexes))
    # 
    # #Close file
    # file_output.write('\n\n\n**********************\n\n')    
    # file_output.close()
    #===========================================================================
###########################################################################################################################################################################    




###########################################################################
def write_matrix_in_worksheet(workbook_object, worksheet_object, matrix_as_line_list, mode, standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes):
    
    ### Basic Mode ###############################################
    if (mode == 'basic'):                
        row=0
        col=0
        for line in matrix_as_line_list:
            line_as_cells = line.strip().split('\t')
            worksheet_object.write_row(row, col, line_as_cells)
            row+=1
    ###############################################################
    
    ### Feature Rich mode #########################################
    elif (mode == 'feature_rich'):
        
        # Prepare matrix_as_cells: a list of list of kind matrix_as_cells[here_the_row][here_the_column]
        matrix_as_cells = []
        for line in matrix_as_line_list:
            line_as_cells = line.strip().split('\t')
            matrix_as_cells.append(line_as_cells)
            
        # Prepare further lists of indexes
        genome_location_columns_indexes = [0,1,2]
        
        # Define different formats
        
        #Format for 'genome location' data label
        genome_location_label_format = workbook_object.add_format()
        genome_location_label_format.set_font_name('Arial')
        genome_location_label_format.set_font_size(14)
        #genome_location_label_format.set_bg_color('#C0C0C0') #Silver
        
        #Format for 'genome location' data type
        genome_location_format = workbook_object.add_format()
        genome_location_format.set_font_name('Arial')
        genome_location_format.set_font_size(12)
        #genome_location_format.set_bg_color('#C0C0C0') #Silver
        
        #Format for 'standard' data label
        pass        
        #Format for 'standard' data type
        standard_data_format = workbook_object.add_format()
        
        #Format for 'merged' data label
        pass        
        #Format for 'merged' data type
        merged_data_format = workbook_object.add_format()
        
        #Format for 'total' label
        pass        
        #Format for 'total' data type        
        total_format = workbook_object.add_format()
        
        #Format for 'collided' data label
        pass        
        #Format for 'collided' data type
        collision_data_format = workbook_object.add_format()
        
        
        
        # WRITING LOOPS: Label worksheet column
        i = 0 #Set row as first
        
        # genome location labels
        for j in genome_location_columns_indexes:
            worksheet_object.write(i, j, matrix_as_cells[i][j], genome_location_label_format)
        # standard data labels
        for j in standard_data_columns_indexes:
            worksheet_object.write(i, j, matrix_as_cells[i][j])
        # merged data labels
        for j in merged_data_columns_indexes:
            worksheet_object.write(i, j, matrix_as_cells[i][j])
        # total data labels
        j+=1
        worksheet_object.write(i, j, matrix_as_cells[i][j])
        # collision data labels
        for j in collision_data_columns_indexes:
            worksheet_object.write(i, j, matrix_as_cells[i][j])
        
        

        
        # WRITING DOUBLE LOOP: Fill the worksheet formatting cells in place
        for i in standard_data_rows_indexes:
            
            # (i,j) Genome location data cell
            # format = genome_location_format
            for j in genome_location_columns_indexes:
                worksheet_object.write(i, j, matrix_as_cells[i][j], genome_location_format)
            
            # (i,j) Standard data cell
            # format = standard_data_format
            for j in standard_data_columns_indexes:
                if (matrix_as_cells[i][j] == '0'):
                    pass
                else:
                    worksheet_object.write(i, j, matrix_as_cells[i][j])
            
            # (i,j) Merged data cell
            # format = merged_data_format
            for j in merged_data_columns_indexes:
                if (matrix_as_cells[i][j] == '0'):
                    pass
                else:
                    worksheet_object.write(i, j, matrix_as_cells[i][j])
            
            # (i,j) Total cell
            # format = total_format
            j+=1
            worksheet_object.write(i, j, matrix_as_cells[i][j])
            
            # (i,j) Collision cell
            # format = collision_data_format
            for j in collision_data_columns_indexes:
                if (matrix_as_cells[i][j] == '0'):
                    pass
                else:
                    worksheet_object.write(i, j, matrix_as_cells[i][j])


    ###############################################################

###########################################################################




#################################################################################################################################
def find_indexes(matrix_as_line_list, matrix_as_line_list_collided = None):
    
    #########################################
    # For redundant give only first argument
    # For IS give both, in ANY case
    # It returns lists of indexes 
    #########################################

    # Row indexes about standard data
    # By default first row is supposed to be the header
    # Other rows are supposed to be 'standard data' till the end
    standard_data_rows_indexes = range(1, len(matrix_as_line_list))    
    
    # Columns indexes
    # By default first 3 columns are supposed to be 'as usual': (chromosome/locus/strand)
    standard_data_columns_indexes = [] # Column indexes about standard data    
    merged_data_columns_indexes = [] # Column indexes about merged data    
    collision_data_columns_indexes = [] # Column indexes about collision data

    # Prepare first line
    first_row_as_list_of_cells = matrix_as_line_list[0].strip().split('\t')

    # Merged indexes
    i = 0
    for cell in first_row_as_list_of_cells:
        if (cell[0] == '_'):
            merged_data_columns_indexes.append(i)       
        i+=1
        
    # Standard indexes
    if (len(merged_data_columns_indexes) != 0):
        standard_data_columns_indexes = range(3, merged_data_columns_indexes[0])
    else:
        standard_data_columns_indexes = range(3, (len(first_row_as_list_of_cells)-1))
        
    # Collision indexes
    if (matrix_as_line_list_collided != None):
        first_collided_row_as_list_of_cells = matrix_as_line_list_collided[0].strip().split('\t')
        collision_data_columns_indexes = range(len(first_row_as_list_of_cells), len(first_collided_row_as_list_of_cells))
        
    # Return Results
    return standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes

#################################################################################################################################        
