###Header################################################
header = """

+-------------------------------------------------------+
 Module: Integration_Sites_retrieving_methods
 Author: Stefano Brasca
 Date:  March 12th, 2014
 Contact: brasca.stefano@hsr.it
 Version: 1.0
+-------------------------------------------------------+

 Description:
  - This module contains functions that are able to
    retrieve one (or even more) Integration Site Object(s)
    from an Ensemble of Covered Bases
    PLEASE NOTE: 'classic' method demands ensemble built
    in a specific way - bushamn_bp_rule could be specified
    by user and it sets both maximum distance between each
    covered base in the same ensemble and maximum number
    of bases an ensemble can span
  
 Note for developers:
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


###Import Module(s)###########################
import Classes_for_Integration_Analysis
import Function_for_Gaussian_IS_identification
import Function_for_SkewedGaussian_IS_identification
##############################################

###Requested Package(s) Import#
import copy
import sys
from operator import itemgetter
from operator import attrgetter
###############################





def classic (Covered_bases_ensamble_object, strand_specific):
    '''
     *** This function creates an IS object from a Covered_bases_ensamble_object ***
     
    INPUT: - Covered_bases_ensamble_object. NOTE: Ensemble of Covered Bases have to be created in a specific way in order to use this function for IS
                                                  retrieving (the 'Science MLD WAS papers' one). When user chooses this retrieval method, it happens
                                                  automatically due to "if (IS_method == "classic"):" controls in PROGRAM_CORE function called in main
           - strand_specific: boolean; it specifies if the matrix computation algorithm had to account for strand: generally, the choice made here should
                              reflect the ones previously made elsewhere. For this purpose, you can find a variable called 'strand_specific_choice' in main,
                              retrieved from user input, so the best usage is strand_specific = strand_specific_choice
                              
    OUTPUT: IS_object.
    
    LOGIC: Classic way to retrieve IS from Covered Bases Ensembles (I.E. the 'Science MLD/WAS papers' one). 
           In order to work properly, ensembles need to be built in a specific way (bushamn_bp_rule could be specified by user - DEFAULT IS 3) and it sets both
           maximum distance between each covered base in the same ensemble and maximum number of bases an ensemble can span.
           
           First covered base sets the "integration locus" while the related reads count is fixed as the overall reads count of the whole ensemble.
    '''
    
    #IS object instance
    IS_object = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object, strand_specific = strand_specific)
    
    #IS Covered_bases_list
    IS_object.Covered_bases_list = Covered_bases_ensamble_object.Covered_bases_list
    
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
        
    #Refreshing 'IS_derived' attribute of Covered_bases_ensamble_object
    Covered_bases_ensamble_object.IS_derived = [IS_object]
    
    return IS_object





def refined_Gaussian_IS_identification (Covered_bases_ensamble_object, hist_gauss_normalized_to_peak, interaction_limit, strand_specific_choice):
    '''
     *** This function creates a list of IS object from a Covered_bases_ensamble_object ***
     
    INPUT: - Covered_bases_ensamble_object: no further specifications needed
           - hist_gauss_normalized_to_peak: list of float representing bin heights( len(hist_gauss_normalized_to_peak) is always odd, peak is in the middle)
                                            coming from main -> PROGRAM_CORE ('if (IS_method == "gauss"):')
                                            See Function_for_Gaussian_IS_identification module for further details (especially gaussian_histogram_generator
                                            and normalize_histogram_to_the_peak functions)
           - interaction_limit: a string from interaction_limit = args.interaction_limit in main, passed to PROGRAM_CORE - It should be an 'int'
           - strand_specific_choice: boolean; it specifies if the matrix computation algorithm had to account for strand: generally, the choice made here should
                                     reflect the ones previously made elsewhere. For this purpose, you can find a variable called 'strand_specific_choice' in main,
                                     retrieved from user input, so the best usage is strand_specific_choice = strand_specific_choice
                              
    OUTPUT: IS_list: a list of IS object.
    
    LOGIC: 
    1) Covered_bases_ensamble_object is incrementally split in 'slices' (objects of Covered bases ensemble class, of kind 'CB of peaks plus its vicinity (side CBs), if present') -> CBE_list_of_slices
    2) Each CBE_slice evaluates all the other bases in Covered_bases_ensamble_object through comparison with hist_gauss_normalized_to_peak -> global_score_dic
    3) global_score_dic is converted to a list_of_bases_to_assign (CBs that are 'vicinity of any peak' are not accounted: they will share the destiny of their 'master' (the peak))
    4) list_of_bases_to_assign is used as a guide to create new and updated new_CBE_slice(s) -> new_CBE_list_of_slices
    5) management of duplicate CBs, of CBE_slice absorbed by others and so on through mutual comparisons
    6) new_CBE_list_of_slices allow creation of IS_list
    A wide use of functions from Function_for_Gaussian_IS_identification module is made
    
    * Coherence control of algorithm behavior are performed: some prints may be produced or even a sys.exit call in worst cases *
    '''
    
    ### Var for Coherence Control
    N_reads_before = Covered_bases_ensamble_object.n_total_reads
    N_cb_before = Covered_bases_ensamble_object.n_covered_bases
    
    ### Cast
    interaction_limit = int(interaction_limit)
    
    ### Partitioning Covered_bases_ensamble_object: CBE_list_of_slices is a list of CBE objects
    # created from CB object from Covered_bases_ensamble_object, these CBE slices will contain original CB, not copy!
    # this list is ordered by peak's height.
    CBE_list_of_slices = Function_for_Gaussian_IS_identification.explore_and_split_CBE(Covered_bases_ensamble_object, strand_specific_choice)
    
    ### Creating global_score_dic, a dictionary of kind : {locus:[(CBE_slice, score), (...), ...]}
    # Each locus (key) of Covered_bases_ensamble_object has a list of tuples as item: each tuple[0] is a CBE_slice laying claim for above mentioned locus and each tuple[1] is the 'claim score'
    global_score_dic = Function_for_Gaussian_IS_identification.global_score_dictionary(CBE_list_of_slices, Covered_bases_ensamble_object, hist_gauss_normalized_to_peak, interaction_limit, strand_specific_choice)
    
    ### Create list_of_bases_to_assign, a list of kind: [(covered_base to assign, CBE_slice claiming it with the highest among non-negative scores), (...), ...]
    list_of_bases_to_assign = []
    for covered_base in Covered_bases_ensamble_object.Covered_bases_list:
        is_adjacent = False    
        if (global_score_dic.has_key(covered_base.locus)):
            ordered_score_tuples = sorted(global_score_dic[covered_base.locus], key=itemgetter(1), reverse=True)
            #If 'the winner of the claiming challange for current covered_base is adjacent to any peak (i.e. is in any CBE.Covered_bases_list[1:])
            #it will be not assigned to anyone (its destiny is tied to its peak's one)
            for CBE in CBE_list_of_slices:
                for CB in CBE.Covered_bases_list[1:]:
                    if (covered_base == CB):
                        is_adjacent = True
                        break
                if (is_adjacent == True):
                    break
            #assignment            
            if ((ordered_score_tuples[0][1] >= 0) and (is_adjacent == False)):
                list_of_bases_to_assign.append((covered_base, ordered_score_tuples[0][0]))
                
    ### Reconstruct CBE slices through just-created list_of_bases_to_assign
    # Each CBE_slice in CBE_list_of_slices gains all the bases it deserved, becoming a 'new_CBE_slice' in new_CBE_list_of_slices:
    # management of duplicate CBs, of CBE_slice absorbed by others and so on will be done at the next stage(s)
    new_CBE_list_of_slices =[] # Since CBE_list_of_slices is ordered by peak's height, also new_CBE_list_of_slices will be
    for CBE_slice in CBE_list_of_slices:
        new_CBE_slice = Function_for_Gaussian_IS_identification.reconstruct_CBE_slice(CBE_slice, list_of_bases_to_assign)
        new_CBE_list_of_slices.append(new_CBE_slice)
        
            
    # Now you have to mutually compare new CBE slices in order to manage duplicate
    # Results (action to perform) will be stored in following list for later usage
    list_of_new_CBE_slice_to_remove = [] # list of new_CBE_slice completely absorbed by others (literally or due to peak's absorption)
    list_of_bases_to_remove_from_new_CBE_slices =[] # list of tuples, such as: (new_CBE_slices_to_fix, [list of CB to remove])
    list_of_bases_to_reassign = [] # a list of tuples, such as: (new_CBE_slices_to_fix, [list of CB to add]
    temp_list_of_CB = [] # support variable in loop
    
    # Loop is performed 'all-against-all' (subject_new_CBE_slice in external loop, object_new_CBE_slice in the nested one)
    # Thanks to order, each subject_new_CBE_slice needs comparison only with following object_new_CBE_slice (that's 'i' and new_CBE_list_of_slices[i:])
    i = 0 
    for subject_new_CBE_slice in new_CBE_list_of_slices:
        i+=1
        for object_new_CBE_slice in new_CBE_list_of_slices[i:]:
            if (subject_new_CBE_slice != object_new_CBE_slice): #they have not to be the same
                
                # object_new_CBE_slice completely absorbed by subject_new_CBE_slice
                # NOTE: due to bug fix, list_of_bases_to_assign does not contain 'adjacent-bases' anymore so this case should be redundant
                if (all(CB in subject_new_CBE_slice.Covered_bases_list for CB in object_new_CBE_slice.Covered_bases_list)):
                    list_of_new_CBE_slice_to_remove.append(object_new_CBE_slice)
                
                #object_new_CBE_slice 'chopped' but not completely absorbed subject_new_CBE_slice    
                else: 
                    
                    # object_new_CBE_slice's peak taken by subject_new_CBE_slice
                    if (object_new_CBE_slice.Covered_bases_list[0] in subject_new_CBE_slice.Covered_bases_list):
                        list_of_new_CBE_slice_to_remove.append(object_new_CBE_slice)
                        list_of_bases_to_reassign.append((subject_new_CBE_slice, object_new_CBE_slice.Covered_bases_list))
                        
                    #object_new_CBE_slice's peak NOT taken by subject_new_CBE_slice
                    else:
                        temp_list_of_CB = []
                        for CB in object_new_CBE_slice.Covered_bases_list:
                            if (CB in subject_new_CBE_slice.Covered_bases_list):
                                temp_list_of_CB.append(CB)
                        list_of_bases_to_remove_from_new_CBE_slices.append((object_new_CBE_slice, temp_list_of_CB))
    
    
    # Revert new_CBE_list_of_slices -> new_CBE_list_of_slices_from_bottom
    # Actions stored in lists above will be performed starting from smaller slices
    new_CBE_list_of_slices_from_bottom = new_CBE_list_of_slices[::-1]            
    for new_CBE_slice in new_CBE_list_of_slices_from_bottom:
        
        # are you (new_CBE_slice) signed as 'to waste'?
        
        if (new_CBE_slice in list_of_new_CBE_slice_to_remove): # YES, new_CBE_slice is to waste
            reassignment = False
            list_of_bases_to_redirect = []
            new_CBE_list_of_slices.remove(new_CBE_slice) #...so remove it!
            #If some CB was signed as 'to be assigned to this new_CBE_slice', it has to be redirected to another new_CBE_slice (The one who make this new_CBE_slice to waste, i.e. the one who absorbed it)
            for new_CBE_slice_to_fix, list_of_CB_to_add in list_of_bases_to_reassign:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    list_of_bases_to_redirect = list_of_bases_to_redirect + list_of_CB_to_add
                    reassignment = True
            # Find new_CBE_slice to whom redirect such bases and update list_of_bases_to_reassign (no matter about removing obsolete entry, adding the new ones is enough)
            if (reassignment == True):
                for address_slice in new_CBE_list_of_slices:
                    if (new_CBE_slice.covered_base_of_max in address_slice.Covered_bases_list):
                        list_of_bases_to_reassign.append((address_slice, list_of_bases_to_redirect))
        
        else: # NO, new_CBE_slice is NOT to waste
            # Let it take the further bases it gained 
            for new_CBE_slice_to_fix, list_of_CB_to_add in list_of_bases_to_reassign:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    i = new_CBE_list_of_slices.index(new_CBE_slice)
                    for CB in list_of_CB_to_add:
                        new_CBE_list_of_slices[i].push_in(CB)
            # Let it back down the bases it has to leave (to another CBE slice)
            for new_CBE_slice_to_fix, list_of_CB_to_remove in list_of_bases_to_remove_from_new_CBE_slices:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    i = new_CBE_list_of_slices.index(new_CBE_slice)
                    for CB in list_of_CB_to_remove:
                        check = new_CBE_list_of_slices[i].push_out(CB)
                        if (check == -1): # Another Coherence Control
                            print "\n\n\t You should have never sought this notice: Some trouble occurs."
                            print "\tPlease let the administrator know this fact!"
                            sys.exit("\n\n\t[ERROR]\tQuit.\n\n")

    
            
    # List of IS to return
    IS_list =[]
    
    #Filling IS_list with IS object instantiated with new conveniently-created CBE_slice
    for CBE_slice in new_CBE_list_of_slices:
        # retrieved_IS: Covered_bases_ensamble_object for instance (allow tracking), CBE_slice for all attributes (CB in Covered_bases_list are the same object in Covered_bases_ensamble_object! Not copies)
        retrieved_IS = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object, strand_specific = strand_specific_choice)
        retrieved_IS.Covered_bases_list = CBE_slice.Covered_bases_list
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
        
    ### Coherence Control
    # N_reads_after and N_cb_after have to be equal to N_reads_before and N_cb_before
    N_reads_after = 0
    N_cb_after = 0
    for CBE_slice in new_CBE_list_of_slices:
        N_reads_after = N_reads_after + CBE_slice.n_total_reads
        N_cb_after = N_cb_after + CBE_slice.n_covered_bases
    if ((N_reads_after != N_reads_before) or (N_cb_after != N_cb_before)):
        print "\n\n\n\tSome troubles found in Refined Gaussian IS retrieving method:"
        print "\tsee CHR {0} from locus {1} to {2}".format(Covered_bases_ensamble_object.chromosome, Covered_bases_ensamble_object.starting_base_locus, Covered_bases_ensamble_object.ending_base_locus)
        sys.exit("\n\n\t[ERROR]\tQuit.\n\n")        
    
    ### Re-order IS_list 'along genome' 
    # Each IS_list returned is joined to others then showed in output as they are
    IS_list = sorted(IS_list, key=attrgetter('chromosome', 'integration_locus', 'strand'))
    
    #Refreshing 'IS_derived' attribute of Covered_bases_ensamble_object
    Covered_bases_ensamble_object.IS_derived = IS_list

            
    ### Return Result
    return IS_list
    
    


def refined_SKEWED_Gaussian_IS_identification (Covered_bases_ensamble_object, two_hist_gauss_normalized_to_peak, strand_specific_choice):
    '''
     *** This function creates a list of IS object from a Covered_bases_ensamble_object ***
     
    INPUT: - Covered_bases_ensamble_object: no further specifications needed
           - two_hist_gauss_normalized_to_peak: dictionary; e.g. {'positive':positive_hist_gauss_normalized_to_peak; 'negative':negative_hist_gauss_normalized_to_peak}
                                                both '_hist_gauss_normalized_to_peak' are list of float representing bin heights,returned by SKEWED_gaussian_histogram_generator
                                                in Function_for_SkewedGaussian_IS_identification module: one as it is, the other is reversed.
                                                In tipical usage such a dictionary comes from main -> PROGRAM_CORE ('if (IS_method == "skewedG"):')
           - strand_specific_choice: boolean; it specifies if the matrix computation algorithm had to account for strand: generally, the choice made here should
                                     reflect the ones previously made elsewhere. For this purpose, you can find a variable called 'strand_specific_choice' in main,
                                     retrieved from user input, so the best usage is strand_specific_choice = strand_specific_choice
                              
    OUTPUT: IS_list: a list of IS object.
    
    LOGIC: 
    1) Covered_bases_ensamble_object is incrementally split in 'slices' (objects of Covered bases ensemble class, of kind 'CB of peaks plus its vicinity (side CBs), if present') -> CBE_list_of_slices
    2) Each CBE_slice evaluates all the other bases in Covered_bases_ensamble_object through comparison with the correct (strand-specifically) _hist_gauss_normalized_to_peak -> global_score_dic
    3) global_score_dic is converted to a list_of_bases_to_assign (CBs that are 'vicinity of any peak' are not accounted: they will share the destiny of their peak)
    4) list_of_bases_to_assign is used as a guide to create new and updated new_CBE_slice(s) -> new_CBE_list_of_slices
    5) management of duplicate CBs, of CBE_slice absorbed by others and so on through mutual comparisons
    6) new_CBE_list_of_slices allow creation of IS_list
    NOTEs:
    In contrast to 'regular' gauss method, here we consider as 'vicinity' of any peak only the CB on the immediate left/right (strand specifically)  
    A wide use of functions from Function_for_SkewedGaussian_IS_identification module is made
    
    * NOTE: strands must be of kind (+,-) or (1,2); else sys.exit is called. 
    * Coherence control of algorithm behavior are performed: some prints may be produced or even a sys.exit call in worst cases *
    
    '''
    
    ### Var for Coherence Control
    N_reads_before = Covered_bases_ensamble_object.n_total_reads
    N_cb_before = Covered_bases_ensamble_object.n_covered_bases
    
       
    ### Choose the correct histogram in two_hist_gauss_normalized_to_peak {dictionary}
    CBE_strand = Covered_bases_ensamble_object.strand
    
    # Strand control
    if ((CBE_strand != '+')and(CBE_strand != '-')and(CBE_strand != '1')and(CBE_strand != '2')):
        print "\n\n\n\tSome troubles found in SKEWED_Gaussian IS retrieving method:"
        print "\tStrand type not recognized. Found: '{0}'.".format(CBE_strand)
        sys.exit("\n\n\t[ERROR]\tQuit.\n\n")
        
    # Find key
    key = None
    if ((CBE_strand == '+')or(CBE_strand == '1')):
        key = 'positive'
    else:
        key = 'negative'
        
    # Fix template hist
    hist_gauss_normalized_to_peak = two_hist_gauss_normalized_to_peak[key]
    
    
    ### Partitioning Covered_bases_ensamble_object: CBE_list_of_slices is a list of CBE objects
    # created from CB object from Covered_bases_ensamble_object, these CBE slices will contain original CB, not copy!
    # this list is ordered by peak's height.
    CBE_list_of_slices = Function_for_SkewedGaussian_IS_identification.explore_and_split_CBE(Covered_bases_ensamble_object, strand_specific_choice)   
    
    ### Creating global_score_dic, a dictionary of kind : {locus:[(CBE_slice, score), (...), ...]}
    # Each locus (key) of Covered_bases_ensamble_object has a list of tuples as item: each tuple[0] is a CBE_slice laying claim for above mentioned locus and each tuple[1] is the 'claim score'
    global_score_dic = Function_for_SkewedGaussian_IS_identification.global_score_dictionary(CBE_list_of_slices, Covered_bases_ensamble_object, hist_gauss_normalized_to_peak, strand_specific_choice)
    
    ### Create list_of_bases_to_assign, a list of kind: [(covered_base to assign, CBE_slice claiming it with the highest among non-negative scores), (...), ...]
    list_of_bases_to_assign = []
    for covered_base in Covered_bases_ensamble_object.Covered_bases_list:
        is_adjacent = False    
        if (global_score_dic.has_key(covered_base.locus)):
            ordered_score_tuples = sorted(global_score_dic[covered_base.locus], key=itemgetter(1), reverse=True)
            #If 'the winner of the claiming challange for current covered_base is adjacent to any peak (i.e. is in any CBE.Covered_bases_list[1:])
            #it will be not assigned to anyone (its destiny is tied to its peak's one)
            for CBE in CBE_list_of_slices:
                for CB in CBE.Covered_bases_list[1:]:
                    if (covered_base == CB):
                        is_adjacent = True
                        break
                if (is_adjacent == True):
                    break
            #assignment            
            if ((ordered_score_tuples[0][1] >= 0) and (is_adjacent == False)):
                list_of_bases_to_assign.append((covered_base, ordered_score_tuples[0][0]))
                
    ### Reconstruct CBE slices through just-created list_of_bases_to_assign
    # Each CBE_slice in CBE_list_of_slices gains all the bases it deserved, becoming a 'new_CBE_slice' in new_CBE_list_of_slices:
    # management of duplicate CBs, of CBE_slice absorbed by others and so on will be done at the next stage(s)
    new_CBE_list_of_slices =[] # Since CBE_list_of_slices is ordered by peak's height, also new_CBE_list_of_slices will be
    for CBE_slice in CBE_list_of_slices:
        new_CBE_slice = Function_for_SkewedGaussian_IS_identification.reconstruct_CBE_slice(CBE_slice, list_of_bases_to_assign)
        new_CBE_list_of_slices.append(new_CBE_slice)
        
            
    # Now you have to mutually compare new CBE slices in order to manage duplicate
    # Results (action to perform) will be stored in following list for later usage
    list_of_new_CBE_slice_to_remove = [] # list of new_CBE_slice completely absorbed by others (literally or due to peak's absorption)
    list_of_bases_to_remove_from_new_CBE_slices =[] # list of tuples, such as: (new_CBE_slices_to_fix, [list of CB to remove])
    list_of_bases_to_reassign = [] # a list of tuples, such as: (new_CBE_slices_to_fix, [list of CB to add]
    temp_list_of_CB = [] # support variable in loop
    
    # Loop is performed 'all-against-all' (subject_new_CBE_slice in external loop, object_new_CBE_slice in the nested one)
    # Thanks to order, each subject_new_CBE_slice needs comparison only with following object_new_CBE_slice (that's 'i' and new_CBE_list_of_slices[i:])
    i = 0 
    for subject_new_CBE_slice in new_CBE_list_of_slices:
        i+=1
        for object_new_CBE_slice in new_CBE_list_of_slices[i:]:
            if (subject_new_CBE_slice != object_new_CBE_slice): #they have not to be the same
                
                # object_new_CBE_slice completely absorbed by subject_new_CBE_slice
                # NOTE: due to bug fix, list_of_bases_to_assign does not contain 'adjacent-bases' anymore so this case should be redundant
                if (all(CB in subject_new_CBE_slice.Covered_bases_list for CB in object_new_CBE_slice.Covered_bases_list)):
                    list_of_new_CBE_slice_to_remove.append(object_new_CBE_slice)
                
                #object_new_CBE_slice 'chopped' but not completely absorbed subject_new_CBE_slice    
                else: 
                    
                    # object_new_CBE_slice's peak taken by subject_new_CBE_slice
                    if (object_new_CBE_slice.Covered_bases_list[0] in subject_new_CBE_slice.Covered_bases_list):
                        list_of_new_CBE_slice_to_remove.append(object_new_CBE_slice)
                        list_of_bases_to_reassign.append((subject_new_CBE_slice, object_new_CBE_slice.Covered_bases_list))
                        
                    #object_new_CBE_slice's peak NOT taken by subject_new_CBE_slice
                    else:
                        temp_list_of_CB = []
                        for CB in object_new_CBE_slice.Covered_bases_list:
                            if (CB in subject_new_CBE_slice.Covered_bases_list):
                                temp_list_of_CB.append(CB)
                        list_of_bases_to_remove_from_new_CBE_slices.append((object_new_CBE_slice, temp_list_of_CB))
    
    
    # Revert new_CBE_list_of_slices -> new_CBE_list_of_slices_from_bottom
    # Actions stored in lists above will be performed starting from smaller slices
    new_CBE_list_of_slices_from_bottom = new_CBE_list_of_slices[::-1]            
    for new_CBE_slice in new_CBE_list_of_slices_from_bottom:
        
        # are you (new_CBE_slice) signed as 'to waste'?
        
        if (new_CBE_slice in list_of_new_CBE_slice_to_remove): # YES, new_CBE_slice is to waste
            reassignment = False
            list_of_bases_to_redirect = []
            new_CBE_list_of_slices.remove(new_CBE_slice) #...so remove it!
            #If some CB was signed as 'to be assigned to this new_CBE_slice', it has to be redirected to another new_CBE_slice (The one who make this new_CBE_slice to waste, i.e. the one who absorbed it)
            for new_CBE_slice_to_fix, list_of_CB_to_add in list_of_bases_to_reassign:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    list_of_bases_to_redirect = list_of_bases_to_redirect + list_of_CB_to_add
                    reassignment = True
            # Find new_CBE_slice to whom redirect such bases and update list_of_bases_to_reassign (no matter about removing obsolete entry, adding the new ones is enough)
            if (reassignment == True):
                for address_slice in new_CBE_list_of_slices:
                    if (new_CBE_slice.covered_base_of_max in address_slice.Covered_bases_list):
                        list_of_bases_to_reassign.append((address_slice, list_of_bases_to_redirect))
        
        else: # NO, new_CBE_slice is NOT to waste
            # Let it take the further bases it gained 
            for new_CBE_slice_to_fix, list_of_CB_to_add in list_of_bases_to_reassign:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    i = new_CBE_list_of_slices.index(new_CBE_slice)
                    for CB in list_of_CB_to_add:
                        new_CBE_list_of_slices[i].push_in(CB)
            # Let it back down the bases it has to leave (to another CBE slice)
            for new_CBE_slice_to_fix, list_of_CB_to_remove in list_of_bases_to_remove_from_new_CBE_slices:
                if (new_CBE_slice_to_fix == new_CBE_slice):
                    i = new_CBE_list_of_slices.index(new_CBE_slice)
                    for CB in list_of_CB_to_remove:
                        check = new_CBE_list_of_slices[i].push_out(CB)
                        if (check == -1): # Another Coherence Control
                            print "\n\n\t You should have never sought this notice: Some trouble occurs."
                            print "\tPlease let the administrator know this fact!"
                            sys.exit("\n\n\t[ERROR]\tQuit.\n\n")

    
            
    # List of IS to return
    IS_list =[]
    
    #Filling IS_list with IS object instantiated with new conveniently-created CBE_slice
    for CBE_slice in new_CBE_list_of_slices:
        # retrieved_IS: Covered_bases_ensamble_object for instance (allow tracking), CBE_slice for all attributes (CB in Covered_bases_list are the same object in Covered_bases_ensamble_object! Not copies)
        retrieved_IS = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object, strand_specific = strand_specific_choice)
        retrieved_IS.Covered_bases_list = CBE_slice.Covered_bases_list
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
        
    ### Coherence Control
    # N_reads_after and N_cb_after have to be equal to N_reads_before and N_cb_before
    N_reads_after = 0
    N_cb_after = 0
    for CBE_slice in new_CBE_list_of_slices:
        N_reads_after = N_reads_after + CBE_slice.n_total_reads
        N_cb_after = N_cb_after + CBE_slice.n_covered_bases
    if ((N_reads_after != N_reads_before) or (N_cb_after != N_cb_before)):
        print "\n\n\n\tSome troubles found in Refined SKEWED Gaussian IS retrieving method:"
        print "\tsee CHR {0} from locus {1} to {2}".format(Covered_bases_ensamble_object.chromosome, Covered_bases_ensamble_object.starting_base_locus, Covered_bases_ensamble_object.ending_base_locus)
        sys.exit("\n\n\t[ERROR]\tQuit.\n\n")        
    
    ### Re-order IS_list 'along genome' 
    # Each IS_list returned is joined to others then showed in output as they are
    IS_list = sorted(IS_list, key=attrgetter('chromosome', 'integration_locus', 'strand'))
    
    #Refreshing 'IS_derived' attribute of Covered_bases_ensamble_object
    Covered_bases_ensamble_object.IS_derived = IS_list

            
    ### Return Result
    return IS_list





        
########################################################################################################################################################################################################################        

# DEPRECATED ############################################
# It works, although poorly tested                      #
# refined_Gaussian_IS_identification has been preferred #
# (this approach's too rigid)                           # 
#########################################################

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
    
                    