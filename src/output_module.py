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
def workbook_output (result_dictionary):
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
    
    
    ### REDUNDANT WORKSHEET #######################################################################
    
    # Create Worksheet name    #Note: must be less than 32char
    redundant_worksheet_name = "RedundantReads"
    if (result_dictionary['strand_specific_choice'] == True):
        redundant_worksheet_name = redundant_worksheet_name + "_StrandSpecific"
    
    # Create Worksheet instance
    redundant_worksheet = workbook_output.add_worksheet(redundant_worksheet_name)
    
    # Fill Worksheet with Redundant Reads data
    write_matrix_in_worksheet(redundant_worksheet, result_dictionary['redundant_matrix'], mode = 'basic')
    
        
    ### IS WORKSHEET ##############################################################################
    
    # Create Worksheet name    #Note: must be less than 32char
    IS_worksheet_name = "ISmatrix" + "_{0}_".format(result_dictionary['IS_method'])
    if (result_dictionary['IS_matrix_collided'] != None):
        IS_worksheet_name = IS_worksheet_name + "&collisions"
    
    # Create Worksheet instance
    IS_worksheet = workbook_output.add_worksheet(IS_worksheet_name)
    
    # Select IS data to use
    selected_matrix_as_line_list = None
    if (result_dictionary['IS_matrix_collided'] != None):
        selected_matrix_as_line_list = result_dictionary['IS_matrix_collided']
    else:
        selected_matrix_as_line_list = result_dictionary['IS_matrix']
            
    # Fill Worksheet with IS data
    write_matrix_in_worksheet(IS_worksheet,selected_matrix_as_line_list, mode = 'basic')
    
    
    ### CONCLUSIVE ACTIONS ########################################################################
        
    # Closing
    workbook_output.close()
###########################################################################################################################################################################    




###########################################################################
def write_matrix_in_worksheet(worksheet_object, matrix_as_line_list, mode):
    
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
        
        ###FUTURE RELEASE ###        
        pass#################
        #####################
    ###############################################################

###########################################################################
