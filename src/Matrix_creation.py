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


### Create and print matrix on output file ##############################################################################################
def matrix_output (list_of_Covered_Bases, column_labels, merged_column_labels, file_output_name):
    
    #Open output file
    file_output = open(file_output_name, 'w')
    
    #Print matrix header, first line
    matrix_header = "chr\tintegration_locus\t"+'\t'.join(column_labels)+"\t"+'\t'.join(merged_column_labels)+"\ttotal_sequence_count"
    file_output.write(matrix_header)
    
    #Print each line left
    for covered_base in list_of_Covered_Bases:
        line = "\n{0}\t{1}".format(covered_base.chromosome, covered_base.locus)
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
        
        if (len(merged_column_labels)<1):
            line = line + "\t"
        tot = covered_base.reads_count     
        line = line + "\t" + str(tot)
        file_output.write(line)
    
    #Close output file    
    file_output.close()
    
#################################################################################################
