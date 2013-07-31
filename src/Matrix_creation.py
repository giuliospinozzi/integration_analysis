###Header###############################################
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


### Create and print dataset matrix on output file ##############################################################################################
def matrix_output (list_of_Covered_Bases, column_labels, merged_column_labels, file_output_name):
    
    #Open output file
    file_output = open(file_output_name, 'w')
    
    #Print matrix header, first line
    matrix_header = "chr\tintegration_locus\tstrand\t"+'\t'.join(column_labels)+"\t"+'\t'.join(merged_column_labels)+"\ttotal_sequence_count"
    if (len(merged_column_labels)<1):
        matrix_header = "chr\tintegration_locus\tstrand\t"+'\t'.join(column_labels)+"\ttotal_sequence_count"
    file_output.write(matrix_header)
    
    #Print each line
    for covered_base in list_of_Covered_Bases:
        line = "\n{0}\t{1}\t{2}".format(covered_base.chromosome, covered_base.locus, covered_base.strand)
        tot = 0
        
        for label in column_labels:
            count = "0"
            if (covered_base.selective_reads_count.has_key(label)):
                count = str(covered_base.selective_reads_count[label])
            line = line + "\t" + count
            
        for label in merged_column_labels:
            count = 0
            for key in covered_base.selective_reads_count.keys():
                if (label in key):
                    count = count + long(covered_base.selective_reads_count[key])           
            line = line + "\t" + str(count)
        
        tot = covered_base.reads_count     
        line = line + "\t" + str(tot)
        file_output.write(line)
    
    #Close output file    
    file_output.close()

#################################################################################################################################################



### Create and print IS matrix on output file ###################################################################################################
def IS_matrix_output (list_of_Covered_Bases, IS_Dictionary, Keys_of_Final_Dictionary, file_output_name, IS_method):
    
    #Open output file
    file_output_name = "IS_{0}_method_{1}".format(IS_method, file_output_name)
    file_output = open(file_output_name, 'w')
    
    #Print matrix header, first line
    matrix_header = "chr\tintegration_locus\tstrand\t"+'\t'.join(Keys_of_Final_Dictionary)
    file_output.write(matrix_header)
    
    #Print each line
    for covered_base in list_of_Covered_Bases:
        line = "\n{0}\t{1}\t{2}".format(covered_base.chromosome, covered_base.locus, covered_base.strand)
        
        for label in Keys_of_Final_Dictionary:
            count = "0"
            
            for IS in IS_Dictionary[label]:
                if ((IS.chromosome == covered_base.chromosome) and (IS.integration_locus == covered_base.locus) and (IS.strand == covered_base.strand)):
                    count = str(IS.reads_count)
                    break
            
            line = line + "\t" + count
            
        file_output.write(line)
            
    #Close output file    
    file_output.close()

#################################################################################################################################################


    
#################################################################################################
