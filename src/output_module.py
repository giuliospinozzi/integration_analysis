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
# import xlsxwriter
###############################




def tsv_output (matrix_file_name, matrix_as_line_list):
    
    #Open a file named as 'redundant_matrix_file_name'
    file_output = open(matrix_file_name, 'w')
    
    #Fill file line by line
    for line in matrix_as_line_list:
        file_output.write(line)
    
    #Close file    
    file_output.close()

    

#===============================================================================
# ###### IN DEVELOPMENT ###### (It'll become a function)
# 
# # Create Workbook object
# workbook_file_name = "Temporary_file_name.xlsx" # Tipo Integration_Analysis_dbschema_dbtable.xlsx
# 
# # Set Workbook policy
# workbook_output = xlsxwriter.Workbook(workbook_file_name,
#     {'in_memory': True}, # Memory usage policy: True means 'work fast at the expense of memory usage'
#     {'strings_to_formulas': False}, # String conversion policy: False means 'Don't try to auto-convert strings to formulas'
#     {'strings_to_urls': False}, # String conversion policy: False means 'Don't try to auto-convert strings to URL links'
#     {'strings_to_numbers': False}) # String conversion policy: False means 'Don't try to auto-convert strings to numbers'
# 
# # Set Workbook properties
# title = 'Integration Analysis'
# dataset = 'Temporary_dataset_name' # Tipo dbschema - dbtable
# author = 'Stefano Brasca'
# manager = 'Eugenio Montini'
# company = 'TIGET - Safety of Gene Therapy and Insertional Mutagenesis Research Unit'
# comments = '''Created by Montini's Bioinfo Team: Andrea Calabria - calabria.andrea@hsr.it; Stefano Brasca - brasca.stefano@hsr.it; Giulio Spinozzi - spinozzi.giulio@hsr.it'''
# workbook_output.set_properties({
#     'title':    title,
#     'subject':  dataset,
#     'author':   author,
#     'manager':  manager,
#     'company':  company,
#     'category': '',
#     'keywords': '',
#     'comments': comments})
#===============================================================================

