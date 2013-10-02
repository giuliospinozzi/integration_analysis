###Header################################################
header = """

+------------------------------------------------------+
 Module: Matrix_creation
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 24th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - [...]
  
 Note:
  - [...]

-------------------------------------------------------- 
""" 
########################################################





def simple_redundant_matrix (list_of_Covered_Bases, column_labels, file_output_name, strand_specific = True):
    '''
    *** This function produces 'Simple' Redundant Reads Matrix ***
        ('simple' stands for: - column labels format is mine
                              - no merged columns
                              - no collision               )
                              
    INPUT: - list_of_Covered_Bases: ORDERED list of Covered Base Objects (ORDER 'along genome' IS MANDATORY, at least 'inside the same chromosome';
                                    then chromosomes could be in any order you like, e.g. alphabetical)
           - column_labels: 
    
    '''    
    # File name creation (file_output_name)
    file_output_name_temp = "redundant_reads_matrix"
    if (strand_specific == True):
        file_output_name_temp = file_output_name_temp + "_strand_specific_"
    file_output_name = file_output_name_temp + file_output_name
    # Uncomment following line and the ones at the end to produce output right there
    #file_output = open(("output_"+file_output_name), 'w')
    
    # Create matrix header (first line of following line_list)
    matrix_header = "chr\tintegration_locus\t"
    if (strand_specific == True):
        matrix_header = matrix_header + "strand\t"
    elif (strand_specific == False):
        matrix_header = matrix_header + "aspecific_strand\t"
    matrix_header = matrix_header + '\t'.join(column_labels)+"\t"
    matrix_header = matrix_header + "total_sequence_count"
    
    # Initialize line_list
    line_list = []
    line_list.append(matrix_header)
    
    
    # Creating Matrix line per line
    if (strand_specific == True):
        for covered_base in list_of_Covered_Bases:
            line = "\n{0}\t{1}\t{2}".format(covered_base.chromosome, covered_base.locus, covered_base.strand)
            tot = 0
            
            for label in column_labels:
                count = "0"
                if (covered_base.selective_reads_count.has_key(label)):
                    count = str(covered_base.selective_reads_count[label])
                line = line + "\t" + count
                
            tot = covered_base.reads_count     
            line = line + "\t" + str(tot)
            line_list.append(line)
            
    elif (strand_specific == False):
        for covered_base in list_of_Covered_Bases:
            line = "\n{0}\t{1}\t{2}".format(covered_base.chromosome, covered_base.locus, covered_base.strand_aspecific)
            tot = 0
            
            for label in column_labels:
                count = "0"
                if (covered_base.selective_reads_count.has_key(label)):
                    count = str(covered_base.selective_reads_count[label])
                line = line + "\t" + count
            
            tot = covered_base.reads_count     
            line = line + "\t" + str(tot)
            line_list.append(line)
    
    # To produce output right there, see above    
    #===========================================================================
    # #Write file
    # for line in line_list:
    #     file_output.write(line)
    # 
    # #Close output file    
    # file_output.close()
    #===========================================================================
    
    # Return file_output_name and line_list
    return file_output_name, line_list






def simple_IS_matrix (IS_list, column_labels, file_output_name, IS_method, strand_specific = True):
    
    # File name creation (file_output_name)
    file_output_name_temp = "IS_matrix_{0}".format(IS_method)
    if (strand_specific == True):
        file_output_name_temp = file_output_name_temp + "_strand_specific"
    file_output_name = file_output_name_temp + "_method_{0}".format(file_output_name)
    # Uncomment following line and the ones at the end to produce output right there
    #file_output = open(("output_"+file_output_name), 'w')
    
    # Create matrix header (first line of following line_list)
    all_labels_list = column_labels + ["all"]
    matrix_header = "chr\tintegration_locus\t"
    if (strand_specific == True):
        matrix_header = matrix_header + "strand\t"
    elif (strand_specific == False):
        matrix_header = matrix_header + "aspecific_strand\t"
    matrix_header = matrix_header + '\t'.join(all_labels_list)
    
    # Initialize line_list
    line_list = []
    line_list.append(matrix_header)
    
    # Creating Matrix line per line
    if (strand_specific == True):
    
        for IS in IS_list:
            current_line = "\n{0}\t{1}\t{2}".format(str(IS.chromosome), str(IS.integration_locus), str(IS.strand))
            for column_label in column_labels:
                if (column_label in IS.selective_reads_count.keys()):
                    current_line = current_line + "\t{0}".format(IS.selective_reads_count[column_label])
                else:
                    current_line = current_line + "\t0"
            current_line = current_line + "\t{0}".format(str(IS.reads_count))
            line_list.append(current_line)
                    
    elif (strand_specific == False):
        
        for IS in IS_list:
            current_line = "\n{0}\t{1}\t{2}".format(str(IS.chromosome), str(IS.integration_locus), str(IS.strand_aspecific))
            for column_label in column_labels:
                if (column_label in IS.selective_reads_count.keys()):
                    current_line = current_line + "\t{0}".format(IS.selective_reads_count[column_label])
                else:
                    current_line = current_line + "\t0"
            current_line = current_line + "\t{0}".format(str(IS.reads_count))
            line_list.append(current_line)
        
        
    # To produce output right there, see above             
    #===========================================================================
    # #Write file
    # for line in line_list:
    #     file_output.write(line)
    #         
    # #Close output file    
    # file_output.close()
    #===========================================================================
    
    # Return file_output_name and line_list
    return file_output_name, line_list
