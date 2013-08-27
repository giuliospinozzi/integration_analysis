###Header################################################
header = """

+------------------------------------------------------+
 Module: Collision
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  August 26th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - [...] Up to now only "simple_collision" method
  
 Note:
  - [...]

-------------------------------------------------------- 
""" 
########################################################



###Simple collision#####################################
def simple_collision (current_dataset_IS_matrix_file_name, current_dataset_IS_matrix_as_line_list, dataset_to_collide_IS_matrix_as_line_list, dataset_to_collide_IS_matrix_file_name, delta):
    
    collision_column = []
    collision_column.append(("\t" + dataset_to_collide_IS_matrix_file_name))# Temporary, needs refinements
    
    for current_dataset_line in current_dataset_IS_matrix_as_line_list[1:]:
        collision_count_for_current_genome_location = 0
        current_dataset_line_split = current_dataset_line.split("\t")
        current_genome_location = (current_dataset_line_split[0], current_dataset_line_split[1]) # ("\nchromosome", "locus") tupla, strand ignored
        
        for dataset_to_collide_line in dataset_to_collide_IS_matrix_as_line_list[1:]:
            dataset_to_collide_line_split = dataset_to_collide_line.split("\t")
            dataset_to_collide_genome_location = (dataset_to_collide_line_split[0], dataset_to_collide_line_split[1]) # ("\nchromosome", "locus") tupla, strand ignored
            dataset_to_collide_sc = int(dataset_to_collide_line_split[-1])
            
            if ((dataset_to_collide_genome_location[0] == current_genome_location[0]) and (abs(int(dataset_to_collide_genome_location[1]) - int(current_genome_location[1])) <= delta)):
                collision_count_for_current_genome_location = collision_count_for_current_genome_location + dataset_to_collide_sc
            
            elif ((dataset_to_collide_genome_location[0] == current_genome_location[0]) and (int(dataset_to_collide_genome_location[1]) > int(current_genome_location[1]) + delta)):
                break
                           
            elif (dataset_to_collide_genome_location[0] > current_genome_location[0]):
                break
            
        collision_column.append(("\t" + str(collision_count_for_current_genome_location)))
    
    #===========================================================================
    # ###DEV LOG###
    # filename = "LOG_"+current_dataset_IS_matrix_file_name+"_COLLISION_WITH_"+dataset_to_collide_IS_matrix_file_name+".txt"
    # file_output = open(filename, 'w')
    # for cell in collision_column:
    #     file_output.write(cell+"\n")
    # file_output.close()
    # #############
    #===========================================================================
    
        
    return collision_column
#########################################################        
        
            
            
            
            
            
            