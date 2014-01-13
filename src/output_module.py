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




def tsv_output (matrix_file_name, matrix_as_line_list):
    
    #Open a file named as 'redundant_matrix_file_name'
    file_output = open(matrix_file_name, 'w')
    
    #Fill file line by line
    for line in matrix_as_line_list:
        file_output.write(line)
    
    #Close file    
    file_output.close()

    

###### IN DEVELOPMENT ############################################################################################################

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
 
    # Create Workbook object
    file_name_part = result_dictionary['dataset_name'].replace(".", "_")
    workbook_file_name = "Integration_Analysis_" +  file_name_part + ".xlsx"
     
    # Set Workbook policy
    workbook_output = xlsxwriter.Workbook(workbook_file_name,
        {'in_memory': True, # Memory usage policy: True means 'work fast at the expense of memory usage'
        'strings_to_formulas': False, # String conversion policy: False means 'Don't try to auto-convert strings to formulas'
        'strings_to_urls': False, # String conversion policy: False means 'Don't try to auto-convert strings to URL links'
        'strings_to_numbers': True}) # String conversion policy: False means 'Don't try to auto-convert strings to numbers'
     
    # Set Workbook properties
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
    
    # REDUNDANT WORKSHEET #################################
    redundant_worksheet = workbook_output.add_worksheet("Redundant Reads Matrix")
    
    row=0
    col=0
    for line in result_dictionary['redundant_matrix']:
        line_as_cells = line.strip().split('\t')
        redundant_worksheet.write_row(row, col, line_as_cells)
        row+=1
    #######################################################
    
    # Closing
    workbook_output.close()
    


##################################################################################################################################    

