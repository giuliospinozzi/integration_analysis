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
def stat_report (result_dictionary, bp_rule, interaction_limit, alpha, scale, shape, args_tsv, args_no_xlsx, args_seqTracker):
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
    
    # Additional labels for seqTracking
    seqTracking_labels = ["Header", "Raw Read (longest)", "Trimmed Read (longest)"]
    
    # CB column label
    CB_column_labels = ["Ensemble ID",
                        "IS ID",
                        "Chromosome",
                        "Strand",
                        "Locus",
                        "# Reads"]
    if (args_seqTracker is True):
        CB_column_labels = CB_column_labels + seqTracking_labels
        
    adapt_labels (CB_column_labels, result_dictionary)
                        
    
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
    IS_max_row_len = None
    if (args_seqTracker is True):
        IS_column_labels, IS_max_row_len = append_seqTracking_labels (IS_column_labels, seqTracking_labels, result_dictionary['IS_list'])
    
    adapt_labels (IS_column_labels, result_dictionary)
    
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
    adapt_labels (Ensemble_column_labels, result_dictionary)
    
    # dataset and file kind
    file_name_part = "Integration_Analysis_" + result_dictionary['dataset_name'].replace(".", "_")[9:] + "_StatREPORT"
    
    
    
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
        
        workbook_file_name =  file_name_part + ".xlsx"
         
        # Create Workbook instance and set policy
        workbook_output = xlsxwriter.Workbook(workbook_file_name,
            {'in_memory': True, # Memory usage policy: True means 'work fast at the expense of memory usage'
            'strings_to_formulas': False, # String conversion policy: False means 'Don't try to auto-convert strings to formulas'
            'strings_to_urls': False, # String conversion policy: False means 'Don't try to auto-convert strings to URL links'
            'strings_to_numbers': True}) # String conversion policy: True means 'auto-convert strings to numbers, if possible'
         
        # Set Workbook metadata
        title = 'Integration Analysis [STATISTICAL REPORT]'
        dataset = result_dictionary['dataset_name'].replace(".", " - ")[9:]
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
        ensembles_worksheet_name = "Ensembles" + "_bpRule" + str(bp_rule)
           
        # Create Worksheet instance
        ensembles_worksheet = workbook_output.add_worksheet(ensembles_worksheet_name)
        
        
        ### IS WORKSHEET ##############################################################################
        
        # Create Worksheet name    #Note: must be less than 32char
        IS_worksheet_name = "ISs" + "_" + result_dictionary['IS_method']
        
        if (result_dictionary['IS_method'] == 'gauss'):
            IS_worksheet_name = IS_worksheet_name + "_intLim" + str(interaction_limit) + "_alpha" + str(alpha)
        if (result_dictionary['IS_method'] == 'skewedG'):
            shape = str(shape)
            shape = shape[1:]
            if (shape[-2:] == ".0"):
                shape = shape[:-2]
            IS_worksheet_name = IS_worksheet_name + "_iLim" + str(interaction_limit) + "_sc" + str(scale) + "_sh" + shape   
                       
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
    
    #Sparkline variables
    IS_sparkline_locations = []
    IS_sparkline_ranges = []
    CBE_sparkline_locations = []
    CBE_sparkline_ranges = []
    IS_max_len=0
    CBE_max_len=0
    
    for CBE in result_dictionary['list_of_Covered_bases_ensambles']:
        
        # ENSEMBLES
        cbe_row += 1 # labels already present
        ensemble_line_as_cells = []
        ensemble_line_as_cells.append(get_ID(CBE)) # Ensemble ID
        ensemble_line_as_cells.append(CBE.chromosome)
        if (result_dictionary['strand_specific_choice'] == True):
            ensemble_line_as_cells.append(str(CBE.strand))
        else:
            ensemble_line_as_cells.append(str(CBE.strand_aspecific))
        ensemble_line_as_cells.append(CBE.starting_base_locus)
        ensemble_line_as_cells.append(CBE.ending_base_locus)
        ensemble_line_as_cells.append(CBE.spanned_bases)
        ensemble_line_as_cells.append(CBE.n_covered_bases)
        ensemble_line_as_cells.append(CBE.n_total_reads)
        ensemble_line_as_cells.append(len(CBE.IS_derived))
        # Prepare variables for last two lines
        cbe_column = 8 # 9 cells appended above
        CBE_loci_range = range(CBE.starting_base_locus, CBE.ending_base_locus + 1)
        CBE_reads_count_per_CB = [0]*len(CBE_loci_range)
        
        
        for IS in CBE.IS_derived:
            
            # IS
            is_row += 1 # labels already present
            IS_line_as_cells = []
            IS_line_as_cells.append(get_ID(CBE)) # Ensemble ID
            IS_line_as_cells.append(get_ID(IS)) # IS ID
            IS_line_as_cells.append(IS.chromosome)
            if (result_dictionary['strand_specific_choice'] == True):
                IS_line_as_cells.append(str(IS.strand))
            else:
                IS_line_as_cells.append(str(IS.strand_aspecific))
            IS_line_as_cells.append(IS.starting_base_locus)
            IS_line_as_cells.append(IS.ending_base_locus)
            IS_line_as_cells.append(IS.spanned_bases)
            IS_line_as_cells.append(IS.n_covered_bases)
            IS_line_as_cells.append(IS.integration_locus)
            IS_line_as_cells.append(IS.reads_count)
            locus_of_peak = None
            for cb in IS.Covered_bases_list:
                if (cb.reads_count == IS.peak_height):
                    locus_of_peak = cb.locus
            IS_line_as_cells.append(locus_of_peak) # Locus of peak. Generally different for integration locus
            IS_line_as_cells.append(IS.peak_height)
            IS_line_as_cells.append((float(IS.peak_height)/float(IS.reads_count))*100.0)
            IS_line_as_cells.append((float(IS.reads_count)/float(CBE.n_total_reads))*100.0)
            # Prepare variables for last two lines
            is_column = 13 # 14 cells appended above
            IS_loci_range = range(IS.starting_base_locus, IS.ending_base_locus + 1)
            IS_reads_count_per_CB = [0]*len(IS_loci_range)
            # Adapt IS_reads_count_per_CB to seqTracking columns
            if (args_seqTracker is True):
                n_of_void_cells = IS_max_row_len - len(IS_reads_count_per_CB) - len(IS_line_as_cells) -1 # -1 accounts for sparkline column
                void_cells = [""]*n_of_void_cells
                IS_reads_count_per_CB = IS_reads_count_per_CB + void_cells
                # Add seqTracking columns in IS_reads_count_per_CB, in order to exploit write_row as it is
                IS_reads_count_per_CB.append(IS.longest_seq_header)
                IS_reads_count_per_CB.append(IS.longest_raw_seq)
                IS_reads_count_per_CB.append(IS.longest_final_seq)
            
                        
            for CB in IS.Covered_bases_list:
        
                # COVERED BASES
                cb_row += 1 # labels already present
                CB_line_as_cells = []
                CB_line_as_cells.append(get_ID(CBE)) # Ensemble ID
                CB_line_as_cells.append(get_ID(IS)) # IS ID
                CB_line_as_cells.append(CB.chromosome)
                if (result_dictionary['strand_specific_choice'] == True):
                    CB_line_as_cells.append(str(CB.strand))
                else:
                    CB_line_as_cells.append(str(CB.strand_aspecific))
                CB_line_as_cells.append(CB.locus)
                CB_line_as_cells.append(CB.reads_count)
                if (args_seqTracker is True):
                    CB_line_as_cells.append(CB.longest_seq_header)
                    CB_line_as_cells.append(CB.longest_raw_seq)
                    CB_line_as_cells.append(CB.longest_final_seq)
                
                # CBE operations involving CB
                CBE_reads_count_per_CB[CBE_loci_range.index(CB.locus)] = CB.reads_count
                    
                # IS operations involving CB
                IS_reads_count_per_CB[IS_loci_range.index(CB.locus)] = CB.reads_count
                                   
                #Write (on workbook and/or in CB_stat_as_line_list)
                cb_column = None
                write_row (CB_worksheet, cb_row, cb_column, CB_line_as_cells, CB_stat_as_line_list, args_tsv, args_no_xlsx)
            
                
            # IS again: set up sparklines variables
            if (IS_max_len < len(IS_reads_count_per_CB)):
                IS_max_len = len(IS_reads_count_per_CB)
            
            #Write (on IS_worksheet and/or in IS_stat_as_line_list)    
            write_row (IS_worksheet, is_row, is_column, IS_line_as_cells, IS_stat_as_line_list, args_tsv, args_no_xlsx, IS_reads_count_per_CB)
        
            
        # CBE again: set up sparklines variables
        if (CBE_max_len < len(CBE_reads_count_per_CB)):
            CBE_max_len = len(CBE_reads_count_per_CB)
        
        #Write (on ensembles_worksheet and/or in ensembles_stat_as_line_list)
        write_row (ensembles_worksheet, cbe_row, cbe_column, ensemble_line_as_cells, ensembles_stat_as_line_list, args_tsv, args_no_xlsx, CBE_reads_count_per_CB)
    
    
    ### Define Sparklines parameters and create them
    
    if (args_no_xlsx == False):
        
        # CBE
        for row in range(1,cbe_row+1):
            CBE_sparkline_locations.append(xl_rowcol_to_cell(row, cbe_column+1))
            CBE_sparkline_ranges.append(xl_range(row, cbe_column+2, row, cbe_column+1+CBE_max_len))            
        CBE_sparkline_parameters = {}
        CBE_sparkline_parameters.update({'location': CBE_sparkline_locations})
        CBE_sparkline_parameters.update({'range': CBE_sparkline_ranges})
        CBE_sparkline_parameters.update({'type': 'column'})
        CBE_sparkline_parameters.update({'style': 3})
        ensembles_worksheet.add_sparkline(CBE_sparkline_locations[0], CBE_sparkline_parameters)
    
        # IS 
        for row in range(1,is_row+1):
            IS_sparkline_locations.append(xl_rowcol_to_cell(row, is_column+1))
            IS_sparkline_ranges.append(xl_range(row, is_column+2, row, is_column+1+IS_max_len))
        IS_sparkline_parameters = {}
        IS_sparkline_parameters.update({'location': IS_sparkline_locations})
        IS_sparkline_parameters.update({'range': IS_sparkline_ranges})
        IS_sparkline_parameters.update({'type': 'column'})
        IS_sparkline_parameters.update({'style': 5})
        IS_worksheet.add_sparkline(IS_sparkline_locations[0], IS_sparkline_parameters)
                    
    
    ### CONCLUSIVE ACTIONS ########################################################################
    
    if (args_no_xlsx == False):            
        ### Closing Workbook
        workbook_output.close()
        
    if (args_tsv == True):
        ### Produce CB *.tsv files
        
        #Name
        CB_tsv_file_name = file_name_part + "_CBs-File"
        if (result_dictionary['strand_specific_choice'] == True):
            CB_tsv_file_name = CB_tsv_file_name + "_StrandSpecific"
        else:
            CB_tsv_file_name = CB_tsv_file_name + "_AspecificStrand"
        CB_tsv_file_name = CB_tsv_file_name + ".tsv" 
        #Writing       
        output_module.tsv_output(CB_tsv_file_name, '\n'.join(CB_stat_as_line_list))
        
        ### Produce ensembles *.tsv files
        
        #Name
        ensembles_tsv_file_name = file_name_part  + "_Ensembles-File_bpRule" + str(bp_rule) + ".tsv"
        #Writing
        output_module.tsv_output(ensembles_tsv_file_name, '\n'.join(ensembles_stat_as_line_list))
        
        ### Produce IS *.tsv files
        
        #Name
        IS_tsv_file_name = file_name_part + "_ISs-File_" + result_dictionary['IS_method']
        if (result_dictionary['IS_method'] == 'gauss'):
            IS_tsv_file_name = IS_tsv_file_name + "_intLim" + str(interaction_limit) + "_alpha" + str(alpha)
        IS_tsv_file_name = IS_tsv_file_name + ".tsv"
        #Writing
        output_module.tsv_output(IS_tsv_file_name, '\n'.join(IS_stat_as_line_list))

####################################################################################################################################################




###########################################################

def append_seqTracking_labels (IS_column_labels, seqTracking_labels, IS_list):
    
    IS_max_len=1
    
    for IS in IS_list:
        if (IS_max_len < IS.spanned_bases):
            IS_max_len = IS.spanned_bases
    if (IS_max_len > 1):
        void_cells = [""]*(IS_max_len-1)
        IS_column_labels = IS_column_labels + void_cells
        
    IS_max_row_len = len(IS_column_labels)
    
    IS_column_labels = IS_column_labels + seqTracking_labels
    
    return IS_column_labels, IS_max_row_len
    
###########################################################





###########################################################

def write_labels (worksheet, line_list, column_labels_list, args_tsv, args_no_xlsx):
        
    if (args_no_xlsx == False):
        worksheet.write_row(0, 0, column_labels_list)
        
    if (args_tsv == True):
        if ("Landscape" in column_labels_list):
            deepcopy_column_labels_list = copy.deepcopy(column_labels_list)
            deepcopy_column_labels_list.remove("Landscape")
            line_list.append("\t".join(deepcopy_column_labels_list))
        else:
            line_list.append("\t".join(column_labels_list))            
###########################################################





###########################################################
def write_row (worksheet, row, column, line_as_cells, line_list, args_tsv, args_no_xlsx, reads_count_per_CB = None):
    
    if (args_no_xlsx == False):
        worksheet.write_row(row, 0, line_as_cells)
        if (reads_count_per_CB != None):
            worksheet.write_row(row, column+2, reads_count_per_CB)
    
    if (args_tsv == True):
        if (reads_count_per_CB != None):
            line_list.append('\t'.join(str(cell) for cell in line_as_cells+reads_count_per_CB))
        else:
            line_list.append('\t'.join(str(cell) for cell in line_as_cells))
###########################################################




###########################################################
def get_ID (integration_analysis_object):
    
    object_id = str(integration_analysis_object)
    
    object_id = object_id.split('.')[1]
    object_id = object_id[:-1]
    object_id = object_id.split(' instance at ')
    
    my_id = '_'.join(object_id)
    
    return my_id
###########################################################




#####################################################################
def adapt_labels (label_list, result_dictionary):
    
    # Strand
    if (result_dictionary['strand_specific_choice'] == False):
        label_list[label_list.index("Strand")] = "Aspecific_Strand"
#####################################################################



    

