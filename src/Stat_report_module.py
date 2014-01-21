###Header################################################
header = """

+------------------------------------------------------+
 Module: Stat_report_module
 Author: Stefano Brasca
 Date:  January 21th, 2014
 Contact: brasca.stefano@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - This module contains functions to create some
    statistical report in *.tsv format and in
    'Excel Workbook' format (*.xlsx)
    
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
#from xlsxwriter.utility import xl_range, xl_rowcol_to_cell
##########################################################

###Import Module(s)###
#import DB_connection
######################


def stat_report (result_dictionary, bushman_bp_rule, interaction_limit, alpha, args_tsv, args_no_xlsx):
    
    ### DEFINE VARIABLES ###########################################################################
    
    # IS column labels
    IS_column_labels = ["Ensemble ID",
                        "IS ID",
                        "Chromosome",
                        "Strand",
                        "Starting base locus",
                        "Ending base locus",
                        "# Spanned bases",
                        "# Covered bases",
                        "Integration locus",
                        "# Reads",
                        "Locus of peak",
                        "Peak height",
                        "%Reads in peak (over IS reads)",
                        "%Reads in IS (over Ensemble reads)",
                        "Landscape",
                        "ReadCount per CB"]
    
    # Ensemble column labels
    Ensemble_column_labels = ["Ensemble ID",
                              "Chromosome",
                              "Strand",
                              "Starting base locus",
                              "Ending base locus",
                              "# Spanned bases",
                              "# Covered bases",
                              "# Reads",
                              "# IS derived",
                              "Landscape",
                              "ReadCount per CB"]
    
    
    ### CREATE WORKBOOK ###########################################################################
    
    # Create Workbook name
    file_name_part = result_dictionary['dataset_name'].replace(".", "_") + "_StatREPORT"
    workbook_file_name = "Integration_Analysis_" + file_name_part + ".xlsx"
     
    # Create Workbook instance and set policy
    workbook_output = xlsxwriter.Workbook(workbook_file_name,
        {'in_memory': True, # Memory usage policy: True means 'work fast at the expense of memory usage'
        'strings_to_formulas': False, # String conversion policy: False means 'Don't try to auto-convert strings to formulas'
        'strings_to_urls': False, # String conversion policy: False means 'Don't try to auto-convert strings to URL links'
        'strings_to_numbers': True}) # String conversion policy: True means 'auto-convert strings to numbers, if possible'
     
    # Set Workbook metadata
    title = 'Integration Analysis [STATISTICAL REPORT]'
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
    
    
    ### IS WORKSHEET ##############################################################################
    
    # Create Worksheet name    #Note: must be less than 32char
    IS_worksheet_name = "ISs" + "_" + result_dictionary['IS_method']
    
    if (result_dictionary['IS_method'] == 'gauss'):
        IS_worksheet_name = IS_worksheet_name + "_intLim" + str(interaction_limit) + "_alpha" + str(alpha)
           
    # Create Worksheet instance
    IS_worksheet = workbook_output.add_worksheet(IS_worksheet_name)
    
    
    ### ENSEMBLES WORKSHEET #######################################################################
    
    # Create Worksheet name    #Note: must be less than 32char
    ensembles_worksheet_name = "Ensembles" 
    
    if (result_dictionary['strand_specific_choice'] == True):
        ensembles_worksheet_name = ensembles_worksheet_name + "_StrandSpecific"
    ensembles_worksheet_name = ensembles_worksheet_name + "_bpRule" + str(bushman_bp_rule)
    
    # Create Worksheet instance
    ensembles_worksheet = workbook_output.add_worksheet(ensembles_worksheet_name)
    
    
    
    
    
    
    ### CONCLUSIVE ACTIONS ########################################################################
        
    # Closing
    workbook_output.close()




###########################################################

def write_IS_labels ():
    
    pass
###########################################################





###########################################################
def write_IS_row ():
    
    pass
###########################################################




###########################################################

def write_ensemble_labels ():
    
    pass
###########################################################





###########################################################
def write_ensemble_row ():
    
    pass
###########################################################










