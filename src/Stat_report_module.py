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
from xlsxwriter.utility import xl_range, xl_rowcol_to_cell
##########################################################

###Import Module(s)###
import output_module
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
    ensembles_stat_as_line_list = []
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
    write_labels (ensembles_worksheet, ensembles_stat_as_line_list, Ensemble_column_labels, args_tsv, args_no_xlsx)
    write_labels (IS_worksheet, IS_stat_as_line_list, IS_column_labels, args_tsv, args_no_xlsx)

    
    ###############################################################################################
    ### TRIPLE LOOP -> CBE -> IS -> CB ############################################################
    ###############################################################################################
    
    #Counters
    cbe_row = 0
    is_row = 0
    cb_row = 0
    
    for CBE in result_dictionary['list_of_Covered_bases_ensambles']:
        
        # ENSEMBLES
        cbe_row += 1 # labels already present
        ensemble_line_as_cells = []
        ensemble_line_as_cells.append(str(CBE)) # Ensemble ID
        ensemble_line_as_cells.append(CBE.chromosome)
        ensemble_line_as_cells.append(CBE.strand) # Good in my opinion, even if is 'None' in case of 'aspecific strand'
        ensemble_line_as_cells.append(CBE.starting_base_locus)
        ensemble_line_as_cells.append(CBE.ending_base_locus)
        ensemble_line_as_cells.append(CBE.spanned_bases)
        ensemble_line_as_cells.append(CBE.n_covered_bases)
        ensemble_line_as_cells.append(CBE.n_total_reads)
        ensemble_line_as_cells.append(len(CBE.IS_derived))
        # Prepare variables for last two lines
        cbe_column = 8 # 9 cells appended above
        CBE_loci_range = range(CBE.starting_base_locus, CBE.ending_base_locus + 1)
        CBE_reads_count_per_CB = []
        
        
        for IS in CBE.IS_derived:
            
            # IS
            is_row += 1 # labels already present
            IS_line_as_cells = []
            IS_line_as_cells.append(str(CBE)) # Ensemble ID
            IS_line_as_cells.append(str(IS)) # IS ID
            IS_line_as_cells.append(IS.chromosome)
            IS_line_as_cells.append(IS.strand) # Good in my opinion, even if is 'None' in case of 'aspecific strand'
            IS_line_as_cells.append(IS.starting_base_locus)
            IS_line_as_cells.append(IS.ending_base_locus)
            IS_line_as_cells.append(IS.spanned_bases)
            IS_line_as_cells.append(IS.n_covered_bases)
            IS_line_as_cells.append(IS.integration_locus)
            IS_line_as_cells.append(IS.reads_count)
            IS_line_as_cells.append(max([covered_base.locus for covered_base in IS.Covered_bases_list])) # Locus of peak. Generally different for integration locus
            IS_line_as_cells.append(IS.peak_height)
            IS_line_as_cells.append((float(IS.peak_height)/float(IS.reads_count))*100.0)
            IS_line_as_cells.append((float(IS.reads_count)/float(CBE.n_total_reads))*100.0)
            # Prepare variables for last two lines
            is_column = 13 # 14 cells appended above
            IS_loci_range = range(IS.starting_base_locus, IS.ending_base_locus + 1)
            IS_reads_count_per_CB = []
            
                        
            for CB in IS.Covered_bases_list:
        
                # COVERED BASES
                cb_row += 1 # labels already present
                CB_line_as_cells = []
                CB_line_as_cells.append(str(CBE)) # Ensemble ID
                CB_line_as_cells.append(str(IS)) # IS ID
                CB_line_as_cells.append(CB.chromosome)
                CB_line_as_cells.append(CB.strand) # Good in my opinion, even if is 'None' in case of 'aspecific strand'
                CB_line_as_cells.append(CB.locus)
                CB_line_as_cells.append(CB.reads_count)
                
                # CBE operations involving CB
                if (CB.locus in CBE_loci_range):
                    CBE_reads_count_per_CB.append(CB.reads_count)
                else:
                    CBE_reads_count_per_CB.append(0)
                    
                # IS operations involving CB
                if (CB.locus in IS_loci_range):
                    IS_reads_count_per_CB.append(CB.reads_count)
                else:
                    IS_reads_count_per_CB.append(0)
                                   
                #Write (on workbook and/or in CB_stat_as_line_list)
                cb_column = None
                write_row (CB_worksheet, cb_row, cb_column, CB_line_as_cells, CB_stat_as_line_list, args_tsv, args_no_xlsx)
            
                
            # IS again: set up last two lines
            IS_sparkline_cell = xl_rowcol_to_cell(is_row, is_column+1)
            IS_sparkline_parameters = {}
            IS_sparkline_range = xl_range(is_row, is_column+2, is_row, is_column+1+len(IS_reads_count_per_CB))
            IS_sparkline_parameters.update({'range': IS_sparkline_range})
            IS_sparkline_parameters.update({'type': 'column'})
            IS_sparkline_parameters.update({'style': 5})
            
            #Write (on IS_worksheet and/or in IS_stat_as_line_list)    
            write_row (IS_worksheet, is_row, is_column, IS_line_as_cells, IS_stat_as_line_list, args_tsv, args_no_xlsx, IS_sparkline_cell, IS_sparkline_parameters, IS_reads_count_per_CB)
        
            
        # CBE again: set up last two lines
        CBE_sparkline_cell = xl_rowcol_to_cell(cbe_row, cbe_column+1)
        CBE_sparkline_parameters = {}
        CBE_sparkline_range = xl_range(cbe_row, cbe_column+2, cbe_row, cbe_column+1+len(CBE_reads_count_per_CB))
        CBE_sparkline_parameters.update({'range': CBE_sparkline_range})
        CBE_sparkline_parameters.update({'type': 'column'})
        CBE_sparkline_parameters.update({'style': 3})
            
        
        ### Here write ensemble_line_as_cells ###
        write_row (ensembles_worksheet, cbe_row, cbe_column, ensemble_line_as_cells, ensembles_stat_as_line_list, args_tsv, args_no_xlsx, CBE_sparkline_cell, CBE_sparkline_parameters, CBE_reads_count_per_CB)    

    
    ### CONCLUSIVE ACTIONS ########################################################################
    
    if (args_no_xlsx == False):            
        ### Closing Workbook
        workbook_output.close()
        
    if (args_tsv == True):
        ### Produce CB *.tsv files
        
        #Name
        CB_tsv_file_name = "Integration_Analysis_" + file_name_part + "_CBs_file"
        if (result_dictionary['strand_specific_choice'] == True):
            CB_tsv_file_name = CB_tsv_file_name + "_StrandSpecific"
        else:
            CB_tsv_file_name = CB_tsv_file_name + "_AspecificStrand"
        CB_tsv_file_name = CB_tsv_file_name + ".tsv" 
        #Writing       
        output_module.tsv_output(CB_tsv_file_name, '\n'.join(CB_stat_as_line_list))
        
        ### Produce ensembles *.tsv files
        
        #Name
        ensembles_tsv_file_name = "Integration_Analysis_" + file_name_part  + "_Ensembles_file_bpRule" + str(bushman_bp_rule) + ".tsv"
        #Writing
        output_module.tsv_output(ensembles_tsv_file_name, '\n'.join(ensembles_stat_as_line_list))
        
        ### Produce IS *.tsv files
        
        #Name
        IS_tsv_file_name = "Integration_Analysis_" + file_name_part + "_ISs_file_" + result_dictionary['IS_method']
        if (result_dictionary['IS_method'] == 'gauss'):
            IS_tsv_file_name = IS_tsv_file_name + "_intLim" + str(interaction_limit) + "_alpha" + str(alpha)
        IS_tsv_file_name = IS_tsv_file_name + ".tsv"
        #Writing
        output_module.tsv_output(IS_tsv_file_name, '\n'.join(IS_stat_as_line_list))

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
def write_row (worksheet, row, column, line_as_cells, line_list, args_tsv, args_no_xlsx, sparkline_cell = None, sparkline_parameters = None, reads_count_per_CB = None):
    
    if (args_no_xlsx == False):
        worksheet.write_row(row, 0, line_as_cells)
        if (reads_count_per_CB != None):
            worksheet.add_sparkline(sparkline_cell, sparkline_parameters) #column+1
            worksheet.write_row(row, column+2, reads_count_per_CB)
    
    if (args_tsv == True):
        if (reads_count_per_CB != None):
            line_list.append('\t'.join(str(cell) for cell in line_as_cells+reads_count_per_CB))
        else:
            line_list.append('\t'.join(str(cell) for cell in line_as_cells))
###########################################################





###########################################################



###########################################################






