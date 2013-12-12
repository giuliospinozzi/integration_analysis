###Header################################################
header = """

+-------------------------------------------------------+
 Module: Integration_Sites_retrieving_methods
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 31th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+-------------------------------------------------------+

 Description:
  - This module contains functions that are able to
    retrieve one (or even more) Integration Site Object(s)
    from an Ensemble of Covered Bases appropriately
    created (e.g. 'classic' method demands ensemble built
    in a specific way - bushamn_bp_rule could be specified
    by user and it sets both maximum distance between each
    covered base in the same ensemble and maximum number
    of bases an ensemble can span) 
  
 Note:
  - For each new method added, you'll surely have to update
    something in Integration_Analysis.py file:
    * IS_methods_list variable (see IS method tuning box)
    * See 'if (IS_method == "whatever"): ' template in #Integration
      Sites Retrieving# part in PROGRAM_CORE function
    
    Maybe you would like to make some minor refinements:
    * See check_method function in Preliminary_controls module
    * Think about delta variable in #COLLISION step# in main
    
    !! Use classic function as input/output template !!

--------------------------------------------------------- 
""" 
########################################################


###Import Module(s)####################
import Classes_for_Integration_Analysis
import Function_for_Gaussian_IS_identification
import copy
from operator import itemgetter
#######################################





def classic (Covered_bases_ensamble_object, strand_specific = True):
    '''
     *** This function creates an IS object from a Covered_bases_ensamble_object ***
     
    INPUT: - Covered_bases_ensamble_object. NOTE: Ensemble of Covered Bases have to be created in a specific way in order to use this function for IS
                                                  retrieving (the 'Science MLD WAS papers' one). When user chooses this retrieval method, it happens
                                                  automatically due to "if (IS_method == "classic"):" controls in PROGRAM_CORE function called in main
           - strand_specific: boolean; it specifies if the matrix computation algorithm had to account for strand: generally, the choice made here should
                              reflect the ones previously made elsewhere. For this purpose, you can find a variable called 'strand_specific_choice' in main,
                              retrieved from user input, so the best usage is strand_specific = strand_specific_choice
                              
    OUTPUT: IS_object.
    
    LOGIC: Classic way to retrieve IS from Covered Bases Ensembles (I.E. the 'Science MLD WAS papers' one). 
           In order to work properly, ensembles need to be built in a specific way (bushamn_bp_rule could be specified by user - DEFAULT IS 3) and it sets both
           maximum distance between each covered base in the same ensemble and maximum number of bases an ensemble can span.
           
           First covered base sets the "integration locus" while the related reads count is fixed as the overall reads count of the whole ensemble.
    '''
    
    #IS object instance
    IS_object = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object, strand_specific = strand_specific)
    
    #IS starting base locus
    IS_object.starting_base_locus = Covered_bases_ensamble_object.starting_base_locus
    
    #IS ending base locus
    IS_object.ending_base_locus = Covered_bases_ensamble_object.ending_base_locus
    
    #Set Integration Locus
    IS_object.integration_locus = Covered_bases_ensamble_object.starting_base_locus
    
    #IS spanned bases
    IS_object.spanned_bases = Covered_bases_ensamble_object.spanned_bases
    
    #IS covered.bases
    IS_object.n_covered_bases = Covered_bases_ensamble_object.n_covered_bases
    
    #Overall Reads Count
    IS_object.reads_count = Covered_bases_ensamble_object.n_total_reads
    
    #Selective Reads Count
    IS_object.selective_reads_count = {}
    for covered_base in Covered_bases_ensamble_object.Covered_bases_list:
        for label in covered_base.selective_reads_count.keys():
            if (label in IS_object.selective_reads_count):
                current_count = IS_object.selective_reads_count[label] + covered_base.selective_reads_count[label]
                IS_object.selective_reads_count.update({label:current_count})
            else:
                IS_object.selective_reads_count.update({label:covered_base.selective_reads_count[label]})
                
    #Peak Height
    IS_object.peak_height = Covered_bases_ensamble_object.covered_base_of_max.reads_count
    
    #Reads Keys List
    IS_object.reads_key_list = []
    for covered_base in Covered_bases_ensamble_object.Covered_bases_list:
        IS_object.reads_key_list = IS_object.reads_key_list + covered_base.list_of_reads_key
    
    return IS_object





def Gaussian_IS_identification (Covered_bases_ensamble_object, hist_gauss_normalized_to_peak, interaction_limit, strand_specific_choice):
    
    #Cast
    interaction_limit = int(interaction_limit)
    
    # Copy of Covered_bases_ensamble_object
    current_ensemble = copy.deepcopy(Covered_bases_ensamble_object)
    
    # N of bases to assign
    bases_to_assign = Covered_bases_ensamble_object.n_covered_bases
    
    # List of IS to return
    IS_list =[]
    
    
    while (bases_to_assign > 0):
        
        # Build histogram for current ensemble
        current_ensemble_bin_areas, current_ensemble_list_of_loci, current_ensemble_max_read_count, current_ensemble_index_of_max = Function_for_Gaussian_IS_identification.CBE__histogram_generator(current_ensemble)
        current_ensemble_bin_areas_normalized = Function_for_Gaussian_IS_identification.normalize_histogram_to_the_peak(current_ensemble_bin_areas, current_ensemble_index_of_max)
        del current_ensemble_max_read_count
        
        # n of allowed step left and right from the current ensemble's peak
        index_last_bin = len(current_ensemble_bin_areas) - 1
        n_step_right = index_last_bin - current_ensemble_index_of_max
        if (n_step_right > interaction_limit):
            n_step_right = interaction_limit
        n_step_left = current_ensemble_index_of_max
        if (n_step_left > interaction_limit):
            n_step_left = interaction_limit
            
        # starting and ending indexes for hist_gauss
        starting_index = interaction_limit - n_step_left
        ending_index = interaction_limit + n_step_right
        
        # list of allowed indexes hist_gauss
        allowed_indexes_gauss = range(starting_index, ending_index+1)
        
        # starting and ending indexes for current_ensemble_bin_areas_normalized
        starting_index = current_ensemble_index_of_max - n_step_left
        ending_index = current_ensemble_index_of_max + n_step_right
        
        # list of allowed indexes for current_ensemble_bin_areas_normalized
        allowed_indexes_CBE = range(starting_index, ending_index+1)
        
        # list of indexes tuples [(allowed_indexes_gauss1, allowed_indexes_CBE1), (allowed_indexes_gauss2, allowed_indexes_CBE2), ... ]
        indexes_tuples = []
        for i in range(0, n_step_left+n_step_right+1):
            indexes_tuples.append((allowed_indexes_gauss[i],allowed_indexes_CBE[i]))
            
        # comparison loop: i from allowed_indexes_gauss and j from allowed_indexes_CBE. IS_indexes is the list of these j that passed if statement
        IS_indexes =[]
        for i,j in indexes_tuples:
            if (hist_gauss_normalized_to_peak[i] >= current_ensemble_bin_areas_normalized[j]):
                IS_indexes.append(j)
        
        # create temp_ensamble: Covered Bases Ensemble made only with gauss-selected covered bases from current_ensemble
        # since temp_ensamble is created adding covered bases by means of push_in method, all its attribute are correct
        # conversely, current_ensemble here is used just like a 'covered bases container, emptied step-by-step
        temp_ensamble = None                   
        for j in IS_indexes:
            for covered_base in current_ensemble.Covered_bases_list:
                if (covered_base.locus == current_ensemble_list_of_loci[j]):
                    if (temp_ensamble == None):
                        temp_ensamble = Classes_for_Integration_Analysis.Covered_bases_ensamble(covered_base, strand_specific = strand_specific_choice)
                        current_ensemble.Covered_bases_list.remove(covered_base)
                        bases_to_assign = bases_to_assign - 1
                        break
                    else:
                        check = temp_ensamble.push_in(covered_base)
                        if (check == 1):
                            current_ensemble.Covered_bases_list.remove(covered_base)
                            bases_to_assign = bases_to_assign - 1
                        break
        
        # retrieved_IS: Covered_bases_ensamble_object for instance, temp_ensamble for attributes
        retrieved_IS = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object, strand_specific = strand_specific_choice)
        retrieved_IS.starting_base_locus = temp_ensamble.starting_base_locus
        retrieved_IS.ending_base_locus = temp_ensamble.ending_base_locus
        retrieved_IS.integration_locus = temp_ensamble.covered_base_of_max.locus
        retrieved_IS.spanned_bases = temp_ensamble.spanned_bases
        retrieved_IS.n_covered_bases = temp_ensamble.n_covered_bases
        retrieved_IS.reads_count = temp_ensamble.n_total_reads
        
        retrieved_IS.selective_reads_count = {}
        for covered_base in temp_ensamble.Covered_bases_list:
            for label in covered_base.selective_reads_count.keys():
                if (label in retrieved_IS.selective_reads_count):
                    current_count = retrieved_IS.selective_reads_count[label] + covered_base.selective_reads_count[label]
                    retrieved_IS.selective_reads_count.update({label:current_count})
                else:
                    retrieved_IS.selective_reads_count.update({label:covered_base.selective_reads_count[label]})

        retrieved_IS.peak_height = temp_ensamble.covered_base_of_max.reads_count
        
        retrieved_IS.reads_key_list = []        
        for covered_base in temp_ensamble.Covered_bases_list:
            retrieved_IS.reads_key_list.append(covered_base.list_of_reads_key)
            
        # append retrieved_IS to IS_list
        IS_list.append(retrieved_IS)
        
    # Return Result
    return IS_list




def refined_Gaussian_IS_identification (Covered_bases_ensamble_object, hist_gauss_normalized_to_peak, interaction_limit, strand_specific_choice):
    
    #Cast
    interaction_limit = int(interaction_limit)
    
    # Partitioning Covered_bases_ensamble_object: CBE_list_of_slices is a list of CBE objects, created from CB object from Covered_bases_ensamble_object, ordered by peak's height in list.
    CBE_list_of_slices = Function_for_Gaussian_IS_identification.explore_and_split_CBE(Covered_bases_ensamble_object, strand_specific_choice)
    
    # Creating global_score_dic, a dictionary of kind : {locus:[(CBE_slice, score), (...), ...]}
    # Each locus (key) of Covered_bases_ensamble_object has as item a list of tuples: each tulpes[0] is a CBE_slice laying claim for above mentioned locus and each tulpes[1] is the 'claim score'
    global_score_dic = Function_for_Gaussian_IS_identification.global_score_dictionary(CBE_list_of_slices, Covered_bases_ensamble_object, hist_gauss_normalized_to_peak, interaction_limit, strand_specific_choice)
    
    # Create list_of_bases_to_assign, a list of kind: [(covered_base to assign, CBE_slice claiming it with highest among positive scores), ...]
    list_of_bases_to_assign = []
    for covered_base in Covered_bases_ensamble_object.Covered_bases_list:
        if (global_score_dic.has_key(covered_base.locus)):
            ordered_score_tuples = sorted(global_score_dic[covered_base.locus], key=itemgetter(1), reverse=True)
            if (ordered_score_tuples[0][1] >= 0):
                list_of_bases_to_assign.append((covered_base, ordered_score_tuples[0][0]))
                
    # Reconstruct CBE slices through list_of_bases_to_assign
    new_CBE_list_of_slices =[] # Since CBE_list_of_slices is ordered by peak's height, also new_CBE_list_of_slices will be
    for CBE_slice in CBE_list_of_slices:
        new_CBE_slice = Function_for_Gaussian_IS_identification.reconstruct_CBE_slice(CBE_slice, global_score_dic, list_of_bases_to_assign)
        new_CBE_list_of_slices.append(new_CBE_slice)
            
    # Now you have to mutually compare new CBE slices in order to manage duplicate
    list_of_new_CBE_slice_to_remove = []
    list_of_bases_to_remove_from_new_CBE_slices =[] # list of tuples, such as: (new_CBE_slices_to_fix, [list of CB to remove])
    list_of_bases_to_reassign = [] # a list of tuples, such as: (new_CBE_slices_to_fix, [list of CB to add]
    temp_list_of_CB = []
    i = 0 
    for subject_new_CBE_slice in new_CBE_list_of_slices:
        i+=1
        for object_new_CBE_slice in new_CBE_list_of_slices[i:]:
            if (subject_new_CBE_slice != object_new_CBE_slice): #they have not to be the same
                if (all(CB in subject_new_CBE_slice.Covered_bases_list for CB in object_new_CBE_slice.Covered_bases_list)): # object_new_CBE_slice inglobato da subject_new_CBE_slice
                    list_of_new_CBE_slice_to_remove.append(object_new_CBE_slice)
                else: #object_new_CBE_slice 'mutilato' ma non inglobato da subject_new_CBE_slice
                    if (object_new_CBE_slice.Covered_bases_list[0] in subject_new_CBE_slice.Covered_bases_list): #il picco di object_new_CBE_slice e' stato preso da subject_new_CBE_slice
                        list_of_new_CBE_slice_to_remove.append(object_new_CBE_slice)
                        list_of_bases_to_reassign.append((subject_new_CBE_slice, object_new_CBE_slice.Covered_bases_list))
                    else: #il picco non e' stato preso
                        temp_list_of_CB = []
                        for CB in object_new_CBE_slice.Covered_bases_list:
                            if (CB in subject_new_CBE_slice.Covered_bases_list):
                                temp_list_of_CB.append(CB)
                        list_of_bases_to_remove_from_new_CBE_slices.append((object_new_CBE_slice, temp_list_of_CB))
    
    
    new_CBE_list_of_slices_from_bottom = new_CBE_list_of_slices[::-1]            
    for new_CBE_slice in new_CBE_list_of_slices_from_bottom:
        
        if (new_CBE_slice in list_of_new_CBE_slice_to_remove):
            reassignment = False
            list_of_bases_to_redirect = []
            new_CBE_list_of_slices.remove(new_CBE_slice)
            for new_CBE_slice_to_fix, list_of_CB_to_add in list_of_bases_to_reassign:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    list_of_bases_to_redirect = list_of_bases_to_redirect + list_of_CB_to_add
                    reassignment = True
            # Find new_CBE_slice to whom redirect
            if (reassignment == True):
                for address_slice in new_CBE_list_of_slices:
                    if (new_CBE_slice.covered_base_of_max in address_slice.Covered_bases_list):
                        list_of_bases_to_reassign.append(address_slice, list_of_bases_to_redirect)
        
        else:
            for new_CBE_slice_to_fix, list_of_CB_to_add in list_of_bases_to_reassign:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    i = new_CBE_list_of_slices.index(new_CBE_slice)
                    for CB in list_of_CB_to_add:
                        new_CBE_list_of_slices[i].push_in(CB)
            for new_CBE_slice_to_fix, list_of_CB_to_remove in list_of_bases_to_remove_from_new_CBE_slices:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    i = new_CBE_list_of_slices.index(new_CBE_slice)
                    for CB in list_of_CB_to_remove:
                        check = new_CBE_list_of_slices[i].push_out(CB)
                        if (check == -1): # Dev
                            print "\n\n\n Vedi questa scritta? Moooolto male!!! Comunicamelo subito \n\n\n" #Dev

# Dev: Schema loop qui sopra                                
        # sei da buttare?
            #SI:
                # ti spettavano basi in piu'?
                    #SI: reindirizzale a chi ti ha inglobato o ti ha sottratto il picco (aggiornare lista list_of_bases_to_reassign) 
                    #NO: bene.
            
            #NO:
                #prenditi le basi in piu' che ti spettano
                #butta le basi che devi buttare
        
            
    # List of IS to return
    IS_list =[]
    
    #Filling IS_list with IS object instantiated with new conveniently-created CBE_slice
    for CBE_slice in new_CBE_list_of_slices:
        # retrieved_IS: Covered_bases_ensamble_object for instance, temp_ensamble for attributes
        retrieved_IS = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object, strand_specific = strand_specific_choice)
        retrieved_IS.starting_base_locus = CBE_slice.starting_base_locus
        retrieved_IS.ending_base_locus = CBE_slice.ending_base_locus
        retrieved_IS.integration_locus = CBE_slice.covered_base_of_max.locus
        retrieved_IS.spanned_bases = CBE_slice.spanned_bases
        retrieved_IS.n_covered_bases = CBE_slice.n_covered_bases
        retrieved_IS.reads_count = CBE_slice.n_total_reads
        
        retrieved_IS.selective_reads_count = {}
        for covered_base in CBE_slice.Covered_bases_list:
            for label in covered_base.selective_reads_count.keys():
                if (label in retrieved_IS.selective_reads_count):
                    current_count = retrieved_IS.selective_reads_count[label] + covered_base.selective_reads_count[label]
                    retrieved_IS.selective_reads_count.update({label:current_count})
                else:
                    retrieved_IS.selective_reads_count.update({label:covered_base.selective_reads_count[label]})

        retrieved_IS.peak_height = CBE_slice.covered_base_of_max.reads_count
        
        retrieved_IS.reads_key_list = []        
        for covered_base in CBE_slice.Covered_bases_list:
            retrieved_IS.reads_key_list.append(covered_base.list_of_reads_key)
            
        # append retrieved_IS to IS_list
        IS_list.append(retrieved_IS)
            
    # Return Result
    return IS_list
    
    
        
        
    
    
                    