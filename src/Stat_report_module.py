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
import copy
import xlsxwriter
#from xlsxwriter.utility import xl_range, xl_rowcol_to_cell
##########################################################

###Import Module(s)###
#import DB_connection
######################




####################################################################################################################################################
def stat_report (result_dictionary, bushman_bp_rule, interaction_limit, alpha, args_tsv, args_no_xlsx):
    '''
    *** This function generates a STAT REPORT file of kind 'Excel Workbook' ***
    
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
    
    ### DEFINE VARIABLES ###########################################################################
    
    # CB column label
    CB_column_labels = ["Ensemble ID",
                        "IS ID",
                        "Chromosome",
                        "Strand",
                        "Locus",
                        "# Reads"]
                        
    
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
    
    
    ### DECLARE VARIABLES ########################################################################
    
    # Workbook variable
    workbook_output = None
    
    # Worksheet variables
    CB_worksheet = None
    ensembles_worksheet = None
    IS_worksheet = None
    
    # Line List variables
    CB_stat_as_line_list = []
    Ensembles_stat_as_line_list = []
    IS_stat_as_line_list = []
    
    
    
    ### CREATE WORKBOOK ###########################################################################
    
    if (args_no_xlsx == False):
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

    
        ### COVERED BASES WORKSHEET ###################################################################    
        
        # Create Worksheet name    #Note: must be less than 32char
        CB_worksheet_name = "CoveredBases" 
        
        if (result_dictionary['strand_specific_choice'] == True):
            CB_worksheet_name = CB_worksheet_name + "_StrandSpecific"
        else:
            CB_worksheet_name = CB_worksheet_name + "_AspecificStrand"
        
        # Create Worksheet instance
        CB_worksheet = workbook_output.add_worksheet(CB_worksheet_name)
        
            
        ### ENSEMBLES WORKSHEET #######################################################################
        
        # Create Worksheet name    #Note: must be less than 32char
        ensembles_worksheet_name = "Ensembles" + "_bpRule" + str(bushman_bp_rule)
           
        # Create Worksheet instance
        ensembles_worksheet = workbook_output.add_worksheet(ensembles_worksheet_name)
        
        
        ### IS WORKSHEET ##############################################################################
        
        # Create Worksheet name    #Note: must be less than 32char
        IS_worksheet_name = "ISs" + "_" + result_dictionary['IS_method']
        
        if (result_dictionary['IS_method'] == 'gauss'):
            IS_worksheet_name = IS_worksheet_name + "_intLim" + str(interaction_limit) + "_alpha" + str(alpha)
               
        # Create Worksheet instance
        IS_worksheet = workbook_output.add_worksheet(IS_worksheet_name)
    
    
    
    ### LABELLING #################################################################################
    write_labels (CB_worksheet, CB_stat_as_line_list, CB_column_labels, args_tsv, args_no_xlsx)
    write_labels (ensembles_worksheet, Ensembles_stat_as_line_list, Ensemble_column_labels, args_tsv, args_no_xlsx)
    write_labels (IS_worksheet, IS_stat_as_line_list, IS_column_labels, args_tsv, args_no_xlsx)

        
    
    
    
    
    ###############################################################################################
    ### LOOP OVER ENSEMBLES #######################################################################
    ###############################################################################################
    
    ### TO ###
    
        
    
    ### CONCLUSIVE ACTIONS ########################################################################
    
    if (args_no_xlsx == False):            
        # Closing Workbook
        workbook_output.close()
        
    # Dev control
    print "\n\n\t\t### DEV CONTROL ###"
    print "\t", CB_stat_as_line_list
    print "\t", Ensembles_stat_as_line_list
    print "\t", IS_stat_as_line_list

####################################################################################################################################################





###########################################################

def write_labels (worksheet, line_list, column_labels_list, args_tsv, args_no_xlsx):
        
    if (args_no_xlsx == False):
        worksheet.write_row(0, 0, column_labels_list)
        
    if (args_tsv == True):
        if (column_labels_list[-2] == "Landscape"):
            deepcopy_column_labels_list = copy.deepcopy(column_labels_list)
            deepcopy_column_labels_list.pop(-2)
            line_list.append("\t".join(deepcopy_column_labels_list))
        else:
            line_list.append("\t".join(column_labels_list))
###########################################################





###########################################################
def write_row (worksheet, line_list, args_tsv, args_no_xlsx):
    
    pass
###########################################################











