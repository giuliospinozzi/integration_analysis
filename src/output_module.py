###Header################################################
from DB_connection import dbOpenConnection
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

###Requested Package(s) Import############################
import xlsxwriter
from xlsxwriter.utility import xl_range, xl_rowcol_to_cell
##########################################################

###Import Module(s)###
import DB_connection
######################




#########################################################
def tsv_output (matrix_file_name, matrix_as_line_list):
    '''
    A function to produce a tsv file (named by
    matrix_file_name string) from matrix_as_line_list.
    
    TIPICAL USAGES:
    tsv_output (redundant_matrix_file_name, redundant_matrix_as_line_list)
    tsv_output (IS_matrix_file_name, IS_matrix_as_line_list)
    (e.g. the one in PROGRAM_CORE function)
    '''
    
    #Open a file named as 'redundant_matrix_file_name'
    file_output = open(matrix_file_name, 'w')
    
    #Fill file line by line
    for line in matrix_as_line_list:
        file_output.write(line)
    
    #Close file    
    file_output.close()
#########################################################

    


#####################################################################################################################################################################################################################################
def workbook_output (result_dictionary, host, user, passwd, port, mode = 'feature_rich'): # or mode = 'basic'
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
            (e.g. the one returned by PROGRAM_CORE function)
            
    OUTPUT: it produces an output summary file of kind 'Excel Workbook' from result_dictionary.
            The workbook file (*.xlsx) has two sheet, the former for redundant matrix and the
            latter for IS matrix / IS matrix collided
                                                    
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
        'strings_to_numbers': True}) # String conversion policy: True means 'auto-convert strings to numbers, if possible'
     
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
    
    # Add coherence controls
    add_coherence_controls(workbook_output, redundant_worksheet, result_dictionary['redundant_matrix'], mode, standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes, result_dictionary['dataset_name'], host, user, passwd, port)
    
        
    ### IS WORKSHEET ##############################################################################
    
    # Create Worksheet name    #Note: must be less than 32char
    IS_worksheet_name = "ISmatrix" + "_{0}".format(result_dictionary['IS_method'])
    if (result_dictionary['IS_matrix_collided'] != None):
        IS_worksheet_name = IS_worksheet_name + "_&collisions"
    
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
    
    # Add coherence controls
    add_coherence_controls(workbook_output, IS_worksheet, selected_matrix_as_line_list, mode, standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes, result_dictionary['dataset_name'], host, user, passwd, port)
    
    
    ### CONCLUSIVE ACTIONS ########################################################################
        
    # Closing
    workbook_output.close()

#####################################################################################################################################################################################################################################




#####################################################################################################################################################################################################################################
def write_matrix_in_worksheet(workbook_object, worksheet_object, matrix_as_line_list, mode, standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes):
    '''
    *** This function write matrix_as_line_list in the worksheet_object of a workbook_object ***
    
    INPUT NOTE: this function can be called in two mode: 'basic' and 'feature_rich': the latter
                requires further inputs (see TIPICAL CONTEXT in find_indexes function) than the former;
                below the input required by both mode.
                
    INPUT: - workbook_object: The XlsxWriter Workbook Object containing the worksheet_object
           - worksheet_object: the sheet of the workbook_object where you want to 
                               write the matrix(_as_line_list)
           - matrix_as_line_list: ...as name says; E.g. one of those returned by PROGRAM_CORE function
                                (result_dictionary['redundant_matrix'], result_dictionary['IS_matrix']
                                or result_dictionary['IS_matrix_collided'])
           - mode: a string among 'basic' and 'feature_rich'
           [...]
           
    OUTPUT: nothing
    
    GENERAL NOTE: actually this function was written just to lighten workbook_output and make the
                  code more readable
                  
    LOGIC: basic mode is straight ('tsv-like'), just see the code below. 
           feature_rich mode is complex; briefly: 
           - the matrix is splitted in cells (array of array, so [i][j] notation is allowed)
           - many format object are defined: the intention is to create formats to couple with
             the different kind of data
           - further lists of indexes, if needed in writing loop, are given/computed; 
           - taking advantage of lists of indexes and formats the sheet is filled cell by cell:
             single loop for the first row of label, double loop for the whole. The code is 
             quite clear, check it!    
    '''
    
    ### Basic Mode ###############################################
    if (mode == 'basic'):                
        row=0
        col=0
        for line in matrix_as_line_list:
            line_as_cells = line.strip().split('\t')
            worksheet_object.write_row(row, col, line_as_cells)
            row+=1
    ###############################################################
    
    ### Feature Rich mode ###################################################################################
    elif (mode == 'feature_rich'):
        
        # Prepare matrix_as_cells: a list of list of kind matrix_as_cells[here_the_row][here_the_column]
        matrix_as_cells = []
        for line in matrix_as_line_list:
            line_as_cells = line.strip().split('\t')
            matrix_as_cells.append(line_as_cells)
        
        
        ### Define different formats ###
        
        #Format for 'genome location' data label
        genome_location_label_format = workbook_object.add_format()
        genome_location_label_format.set_font_name('Arial')
        genome_location_label_format.set_font_size(14)       
        #Format for 'genome location' data type
        genome_location_format = workbook_object.add_format()
        genome_location_format.set_font_name('Arial')
        genome_location_format.set_font_size(12)
        
        #Format for 'standard' data label
        standard_label_format = workbook_object.add_format()
        standard_label_format.set_font_name('Arial')
        standard_label_format.set_italic()
        standard_label_format.set_font_size(14)
        standard_label_format.set_rotation(45)        
        #Format for 'standard' data type
        standard_data_format = workbook_object.add_format()
        standard_data_format.set_font_name('Arial')
        standard_data_format.set_italic()
        standard_data_format.set_font_size(12)
        
        #Format for 'merged' data label
        merged_label_format = workbook_object.add_format()
        merged_label_format.set_font_name('Arial')
        merged_label_format.set_italic()
        merged_label_format.set_bold()
        merged_label_format.set_font_size(14)
        merged_label_format.set_rotation(45)
        #Format for 'merged' data type
        merged_data_format = workbook_object.add_format()
        merged_data_format.set_font_name('Arial')
        merged_data_format.set_italic()
        merged_data_format.set_bold()
        merged_data_format.set_font_size(12)
        
        #Format for 'total' label
        total_label_format = workbook_object.add_format()
        total_label_format.set_font_name('Arial')
        total_label_format.set_bold()
        total_label_format.set_font_size(14)
        #Format for 'total' data type        
        total_format = workbook_object.add_format()
        total_format.set_font_name('Arial')
        total_format.set_bold()
        total_format.set_font_size(12)
        
        #Format for 'collided' data label
        collision_label_format = workbook_object.add_format()
        collision_label_format.set_font_name('Arial')
        collision_label_format.set_bold()
        collision_label_format.set_font_size(14)
        collision_label_format.set_font_color('#800080')
        #Format for 'collided' data type
        collision_data_format = workbook_object.add_format()
        collision_data_format.set_font_name('Arial')
        collision_data_format.set_bold()
        collision_data_format.set_font_size(12)
        collision_data_format.set_font_color('#800080')
        
        
        # Prepare further lists of indexes, needed in writing loop
        genome_location_columns_indexes = [0,1,2]
                
        
        ### LABELLING LOOPS: Label worksheet columns ###
        ########### (writing first row) ################
        i = 0 #Set row as first
        
        # genome location labels
        for j in genome_location_columns_indexes:
            worksheet_object.write(i, j, matrix_as_cells[i][j], genome_location_label_format)
        # standard data labels
        for j in standard_data_columns_indexes:
            worksheet_object.write(i, j, matrix_as_cells[i][j], standard_label_format)
        # merged data labels
        for j in merged_data_columns_indexes:
            worksheet_object.write(i, j, matrix_as_cells[i][j], merged_label_format)
        # total data labels
        j+=1
        worksheet_object.write(i, j, matrix_as_cells[i][j], total_label_format)
        # collision data labels
        for j in collision_data_columns_indexes:
            worksheet_object.write(i, j, matrix_as_cells[i][j], collision_label_format)
        
                
        ### WRITING DOUBLE LOOP: Fill the worksheet formatting cells in place ###
        #################### (writing all the data row) #########################
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
                    worksheet_object.write(i, j, matrix_as_cells[i][j], standard_data_format)
            
            # (i,j) Merged data cell
            # format = merged_data_format
            for j in merged_data_columns_indexes:
                if (matrix_as_cells[i][j] == '0'):
                    pass
                else:
                    worksheet_object.write(i, j, matrix_as_cells[i][j], merged_data_format)
            
            # (i,j) Total cell
            # format = total_format
            j+=1
            worksheet_object.write(i, j, matrix_as_cells[i][j], total_format)
            
            # (i,j) Collision cell
            # format = collision_data_format
            for j in collision_data_columns_indexes:
                if (matrix_as_cells[i][j] == '0'):
                    pass
                else:
                    worksheet_object.write(i, j, matrix_as_cells[i][j], collision_data_format)


    #########################################################################################################

#####################################################################################################################################################################################################################################




#################################################################################################################################
def add_coherence_controls(workbook_object, worksheet_object, matrix_as_line_list, mode, standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes, db_dataset_name, host, user, passwd, port):
    
    # Prepare matrix_as_cells: a list of list of kind matrix_as_cells[here_the_row][here_the_column]
    matrix_as_cells = []
    for line in matrix_as_line_list:
        line_as_cells = line.strip().split('\t')
        matrix_as_cells.append(line_as_cells)
        
    # Define formats
    ################
    #### ToDo ######
    ################


    
    ####### COLUMN CONTROLS #######
    
    
    # Get suitable i (row index)
    i = standard_data_rows_indexes[-1] + 3 # leave 2 rows void
    
    # Write row title
    worksheet_object.write(i, 0, "Columns sums:")
   
    ### Columns sums
    
    # standard columns
    for j in standard_data_columns_indexes:
        cells_range = xl_range(standard_data_rows_indexes[0], j, standard_data_rows_indexes[-1], j)
        worksheet_object.write_formula(i, j, '=SUM({0})'.format(cells_range))
    # merged columns
    for j in merged_data_columns_indexes:
        cells_range = xl_range(standard_data_rows_indexes[0], j, standard_data_rows_indexes[-1], j)
        worksheet_object.write_formula(i, j, '=SUM({0})'.format(cells_range))
    # total column
    j+=1
    cells_range = xl_range(standard_data_rows_indexes[0], j, standard_data_rows_indexes[-1], j)
    worksheet_object.write_formula(i, j, '=SUM({0})'.format(cells_range))

    
    ### DB comparison
    
    # Get suitable i (row index)
    i+=1
    
    # Write row title
    worksheet_object.write(i, 0, "DB data:")
    
    # Get total from DB    
    db, db_table = db_dataset_name.split('.')
    conn = DB_connection.dbOpenConnection(host, user, passwd, port, db)
    #db_total = DB_connection.getTableRowCount(conn, db_table)
    cursor = conn.cursor()  
    cursor.execute ("SELECT count(DISTINCT header) FROM %s WHERE 1" %(db_table))
    db_total = cursor.fetchall()[0][0]
    cursor.close()
    DB_connection.dbCloseConnection(conn)
        
    # Write total
    worksheet_object.write(i, j, db_total)
    
    # Write row title
    worksheet_object.write(i+1, 0, "Check")
    
    # Write formula
    total_from_sum = xl_rowcol_to_cell(i-1, j)
    total_from_db = xl_rowcol_to_cell(i, j)
    worksheet_object.write_formula(i+1, j, '''IF({0}={1},TRUE,FALSE)'''.format(total_from_sum, total_from_db))
    
    
        
    ####### ROW CONTROLS #######
    
    
    # Get suitable j (column index)
    j = None
    if (len(collision_data_columns_indexes) != 0):
        j = collision_data_columns_indexes[-1] + 2 # leave a column void
    else:
        if (len(merged_data_columns_indexes) != 0):
            j = merged_data_columns_indexes[-1] + 3 # leave a column void
        else:
            j = standard_data_columns_indexes[-1] + 3 # leave a column void
        
    # Write column titles
    temp_j = j
    worksheet_object.write(0, temp_j, "Standard data rows sums:")
    temp_j+=1
    worksheet_object.write(0, temp_j, "Check:")
    if (len(merged_data_columns_indexes) != 0):
        temp_j+=1
        worksheet_object.write(0, temp_j, "Merged data rows sums:")
        temp_j+=1
        worksheet_object.write(0, temp_j, "Check:")
        
    # Rows checks
    temp_standard_data_rows_indexes = standard_data_rows_indexes
    temp_standard_data_rows_indexes.append(standard_data_rows_indexes[-1]+3)
    for i in temp_standard_data_rows_indexes:
        
        #Get total cell
        total_cell = xl_rowcol_to_cell(i, j-2-len(collision_data_columns_indexes))
        
        # Standard Sum
        cells_range = xl_range(i, standard_data_columns_indexes[0], i, standard_data_columns_indexes[-1])
        worksheet_object.write_formula(i, j, '=SUM({0})'.format(cells_range))
        
        # Merged Sum
        if (len(merged_data_columns_indexes) != 0):
            cells_range = xl_range(i, merged_data_columns_indexes[0], i, merged_data_columns_indexes[-1])
            worksheet_object.write_formula(i, j+2, '=SUM({0})'.format(cells_range))
        
        # Standard Check
        standard_sum_cell = xl_rowcol_to_cell(i,j)
        worksheet_object.write_formula(i, j+1, '''IF({0}={1},TRUE,FALSE)'''.format(total_cell, standard_sum_cell))
        
        # Merged Check
        if (len(merged_data_columns_indexes) != 0):
            merged_sum_cell = xl_rowcol_to_cell(i,j+2)
            worksheet_object.write_formula(i, j+3, '''IF({0}={1},TRUE,FALSE)'''.format(total_cell, merged_sum_cell))
        
    
        
    
#################################################################################################################################





#################################################################################################################################
def find_indexes(matrix_as_line_list, matrix_as_line_list_collided = None):
    '''
    *** Given matrix(es) as line list, this function provides arrays of indexes ***
           for each category of data - 'straight', merged, (collided), ...
           
    INPUT: matrix_as_line_list - ...as name says; E.g. one of those returned by PROGRAM_CORE function
                                (result_dictionary['redundant_matrix'] or result_dictionary['IS_matrix'])
           matrix_as_line_list_collided -  'collided version' of the matrix above, if exists (e.g. 
                                           result_dictionary['IS_matrix'] and result_dictionary['IS_matrix_collided'])
                                           
    OUTPUT: standard_data_rows_indexes, standard_data_columns_indexes, merged_data_columns_indexes, collision_data_columns_indexes
            (all list of integers)
            
    LOGIC: 
    - about rows indexes:
      standard_data_rows_indexes = range(1, len(matrix_as_line_list)), so it should be always ok, unless you 
      want to add more row beside the default 'label row'
      
    - about column indexes:
      this function starts looking at the label-row, supposed to be matrix_as_line_list[0], in order to catch
      'merged columns', identified by initial '_' (-> merged_data_columns_indexes). Taking advantage of results
      of this search, the function then find standard_data_columns_indexes (the first of which is supposed to
      be indexed by '3' -> standard_data_columns_indexes = range(3, ...) ). In the end, 
      if (matrix_as_line_list_collided != None), collision_data_columns_indexes are computed through a comparison
      between 'matrix_as_line_list' and 'matrix_as_line_list_collided': this way the column of total, for example,
      is automatically accounted for.
      
      TIPICAL USAGES:
      - find_indexes(result_dictionary['redundant_matrix'])
      - find_indexes(result_dictionary['IS_matrix'], result_dictionary['IS_matrix_collided'])
      
      TIPICAL CONTEXT:
      Written to provide inputs necessary for write_matrix_in_worksheet function, when called in
      'feature_rich' mode  
    '''
    
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
