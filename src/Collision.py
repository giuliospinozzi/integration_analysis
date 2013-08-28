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
        
            
def multiple_collision (current_dataset_tuple, list_of_IS_results_tuple, delta):
    
    # current_ dataset_tuple is like an element of list_of_IS_results_tuple             
    # list_of_IS_results_tuple is a list of tuples, created appending results provided by main
    
    #Take the current dataset
    current_dataset_IS_matrix_file_name = current_dataset_tuple[0]
    current_dataset_IS_matrix_as_line_list = current_dataset_tuple[1]
    
    #Initialize current_dataset_IS_matrix_as_line_list_collided
    current_dataset_IS_matrix_as_line_list_collided = [None]
    
    # Preparing first line of current_dataset_IS_matrix_as_line_list_collided, soon containing columns labels and in the end data to return
    number_of_collision = 0 #take advantage of following loop to enumerate collision to perform
    for dataset in list_of_IS_results_tuple:
        if (dataset[0]!=current_dataset_IS_matrix_file_name):
            if (number_of_collision == 0):  
                current_dataset_IS_matrix_as_line_list_collided[0] = current_dataset_IS_matrix_as_line_list[0] + "\t" + dataset[0] #dataset[0] is the name of a dataset
            else:
                current_dataset_IS_matrix_as_line_list_collided[0] = current_dataset_IS_matrix_as_line_list_collided[0] + "\t" + dataset[0] #dataset[0] is the name of a dataset
            number_of_collision +=1
    
    #Loop over each line of current dataset IS matrix
    for current_dataset_line in current_dataset_IS_matrix_as_line_list[1:]:
        current_dataset_line_split = current_dataset_line.split("\t")
        current_genome_location = (current_dataset_line_split[0], current_dataset_line_split[1]) # ("\nchromosome", "locus") tupla, strand ignored
        
        #Loop over each dataset to collide
        temp_row = ""
        for dataset_to_collide in list_of_IS_results_tuple:
            if (dataset_to_collide[0]!=current_dataset_IS_matrix_file_name):
                collision_count = 0
                for dataset_to_collide_line in dataset_to_collide[1][1:]:
                    dataset_to_collide_line_split = dataset_to_collide_line.split("\t")
                    dataset_to_collide_genome_location = (dataset_to_collide_line_split[0], dataset_to_collide_line_split[1]) # ("\nchromosome", "locus") tupla, strand ignored
                    dataset_to_collide_sc = int(dataset_to_collide_line_split[-1])
    
                    if ((dataset_to_collide_genome_location[0] == current_genome_location[0]) and (abs(int(dataset_to_collide_genome_location[1]) - int(current_genome_location[1])) <= delta)):
                        collision_count = collision_count + dataset_to_collide_sc
            
                    elif ((dataset_to_collide_genome_location[0] == current_genome_location[0]) and (int(dataset_to_collide_genome_location[1]) > int(current_genome_location[1]) + delta)):
                        break
                           
                    elif (dataset_to_collide_genome_location[0] > current_genome_location[0]):
                        break
                temp_row = temp_row+"\t"+str(collision_count)
        current_dataset_IS_matrix_as_line_list_collided.append(current_dataset_line+temp_row)

    
    #Return results
    return current_dataset_IS_matrix_file_name, current_dataset_IS_matrix_as_line_list_collided #returning a tuple like the one given in input
            
            
            
            