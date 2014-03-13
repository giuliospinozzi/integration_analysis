###Header################################################
header = """

+------------------------------------------------------+
 Module: Collision
 Author: Stefano Brasca
 Date:  August 26th, 2013
 Contact: brasca.stefano@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - [...]
  
 Note:
  - [...]

-------------------------------------------------------- 
""" 
########################################################


###Requested Packages###
import copy #
########################
        
            
def multiple_collision (current_dataset_tuple, list_of_IS_results_tuple_for_collision, delta):
    '''
    *** This function produce "collisions" between one dataset (current_dataset_tuple) and a list of some others
        (list_of_IS_results_tuple_for_collision) through delta "collision radius" ***
        
    INPUT: - current_dataset_tuple, i.e. a tuple composed by:
               - IS_matrix_file_name
               - IS_matrix_as_line_list
               - result_dictionary['dataset_name']
               returned by PROGRAM_CORE function; this represent the dataset we want to make bump into others
           
           - list_of_IS_results_tuple_for_collision: a list of tuple of 'current_dataset_tuple' kind; this list collects all
                                       target datasets for current_dataset_tuple
                                         
           - delta: integer, a sort of collision-radius, setting the minimum distance between two different
                    genome locations in order to collide each others
                    NOTE: delta is the threshold difference between two integration locus! Therefore, if you
                          want to consider 2 covered bases separated by 3empty loci as "colliding", delta must
                          be set equal to 4!!!
           
    OUTPUT: - current_dataset_IS_matrix_file_name: like IS_matrix_file_name in current_dataset_tuple given in input 
            - current_dataset_IS_matrix_as_line_list_collided: like IS_matrix_as_line_list in current_dataset_tuple
                                                               given in input, but now containing collision columns
            
    LOGIC: [...]
    
    NOTE: Remember - strand is ignored. Is it a problem for strand_specific colliding matrixes?? ...think at 'duplicate'
                                        genome location... 
    
    '''
    
    # current_dataset_tuple is like an element of list_of_IS_results_tuple_for_collision             
    # list_of_IS_results_tuple_for_collision is a list of tuples (IS_matrix_file_name, IS_matrix_as_line_list), created appending results returned by PROGRAM_CORE function
    
    #Take the current dataset
    current_dataset_IS_matrix_file_name = current_dataset_tuple[0]
    current_dataset_IS_matrix_as_line_list = current_dataset_tuple[1]
    
    #Initialize current_dataset_IS_matrix_as_line_list_collided
    current_dataset_IS_matrix_as_line_list_collided = [None]
    
    # Preparing first line of current_dataset_IS_matrix_as_line_list_collided, soon containing columns labels and in the end data to return
    number_of_collision = 0 #take advantage of following loop to enumerate collision to perform
    for dataset in list_of_IS_results_tuple_for_collision:
        if (dataset[0]!=current_dataset_IS_matrix_file_name):
            colliding_dataset_name = "Collision against " + dataset[2]
            if (number_of_collision == 0):  
                #current_dataset_IS_matrix_as_line_list_collided[0] = current_dataset_IS_matrix_as_line_list[0] + "\t" + dataset[0] #dataset[0] is the name of a dataset
                current_dataset_IS_matrix_as_line_list_collided[0] = current_dataset_IS_matrix_as_line_list[0] + "\t" + colliding_dataset_name #dataset[0] is the name of a dataset
            else:
                #current_dataset_IS_matrix_as_line_list_collided[0] = current_dataset_IS_matrix_as_line_list_collided[0] + "\t" + dataset[0] #dataset[0] is the name of a dataset
                current_dataset_IS_matrix_as_line_list_collided[0] = current_dataset_IS_matrix_as_line_list_collided[0] + "\t" + colliding_dataset_name #dataset[0] is the name of a dataset
            number_of_collision +=1
    
    
    line_start_list = [1]*number_of_collision #position tracking
    
    #Loop over each line of current dataset IS matrix
    for current_dataset_line in current_dataset_IS_matrix_as_line_list[1:]:
        current_dataset_line_split = current_dataset_line.split("\t")
        current_genome_location = (current_dataset_line_split[0], current_dataset_line_split[1]) # ("\nchromosome", "locus") tupla, strand ignored
        
        line_start_list_index = 0
                
        #Loop over each dataset to collide
        temp_row = ""
        for dataset_to_collide in list_of_IS_results_tuple_for_collision:
            if (dataset_to_collide[0]!=current_dataset_IS_matrix_file_name):                   
                collision_count = 0
                line_count = copy.deepcopy(line_start_list[line_start_list_index]) # copy.deepcopy
                for dataset_to_collide_line in dataset_to_collide[1][line_start_list[line_start_list_index]:]:
                    line_count += 1
                    dataset_to_collide_line_split = dataset_to_collide_line.split("\t")
                    dataset_to_collide_genome_location = (dataset_to_collide_line_split[0], dataset_to_collide_line_split[1]) # ("\nchromosome", "locus") tupla, strand ignored
                    dataset_to_collide_sc = int(dataset_to_collide_line_split[-1])
    
                    if ((dataset_to_collide_genome_location[0] == current_genome_location[0]) and (abs(int(dataset_to_collide_genome_location[1]) - int(current_genome_location[1])) <= delta)):
                        collision_count = collision_count + dataset_to_collide_sc
            
                    elif ((dataset_to_collide_genome_location[0] == current_genome_location[0]) and (int(dataset_to_collide_genome_location[1]) > int(current_genome_location[1]) + delta)):
                        line_start_list[line_start_list_index] = line_count - delta - 1 # -1
                        if (line_start_list[line_start_list_index] <= 0):
                            line_start_list[line_start_list_index] = 1  
                        break
                           
                    elif (dataset_to_collide_genome_location[0] > current_genome_location[0]):
                        line_start_list[line_start_list_index] = line_count - delta - 1 # -1
                        if (line_start_list[line_start_list_index] <= 0):
                            line_start_list[line_start_list_index] = 1
                        break
                temp_row = temp_row+"\t"+str(collision_count)
                
                line_start_list_index += 1
                
        current_dataset_IS_matrix_as_line_list_collided.append(current_dataset_line+temp_row)

    
    #Return results
    return current_dataset_IS_matrix_file_name, current_dataset_IS_matrix_as_line_list_collided #returning a tuple like the one given in input
            
                                