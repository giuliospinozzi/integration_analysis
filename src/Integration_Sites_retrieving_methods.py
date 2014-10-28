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
    in a specific way - bp_rule could be specified
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
import Function_for_Dynamic_IS_identification
import DB_connection
##############################################

###Requested Package(s) Import#
import sys
from operator import itemgetter
from operator import attrgetter
###############################





def classic (Covered_bases_ensamble_object, strand_specific, center_on_mode=False):
    '''
     *** This function creates an IS object from a Covered_bases_ensamble_object ***
     
    INPUT: - Covered_bases_ensamble_object. NOTE: Ensemble of Covered Bases have to be created in a specific way in order to use this function for IS
                                                  retrieving (the 'Science MLD WAS papers' one). When user chooses this retrieval method, it happens
                                                  automatically due to "if (IS_method == "classic"):" controls in PROGRAM_CORE function called in main
           - strand_specific: boolean; it specifies if the matrix computation algorithm had to account for strand: generally, the choice made here should
                              reflect the ones previously made elsewhere. For this purpose, you can find a variable called 'strand_specific_choice' in main,
                              retrieved from user input, so the best usage is strand_specific = strand_specific_choice
           - center_on_mode: boolean, FALSE by default. If True, the IS_object.integration_locus will be the locus with the highest SC.
                              
    OUTPUT: IS_object.
    
    LOGIC: Classic way to retrieve IS from Covered Bases Ensembles (I.E. the 'Science MLD/WAS papers' one). 
           In order to work properly, ensembles need to be built in a specific way (bp_rule could be specified by user - DEFAULT IS 3) and it sets both
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
    if center_on_mode is False:
        IS_object.integration_locus = Covered_bases_ensamble_object.starting_base_locus
    elif center_on_mode is True:
        CB_of_mode = Covered_bases_ensamble_object.Covered_bases_list[0]
        for CB in Covered_bases_ensamble_object.Covered_bases_list[1:]:
            if CB.reads_count > CB_of_mode.reads_count:
                CB_of_mode = CB
        IS_object.integration_locus = CB_of_mode.locus
    
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
                            print "\n\n\t You should have never seen this notice: Some trouble occurs."
                            print "\tPlease let the administrator know this fact!"
                            sys.exit("\n\n\t[ERROR]\tQuit.\n\n")

    
            
    # List of IS to return
    IS_list =[]
    
    #Filling IS_list with IS object instantiated with new conveniently-created CBE_slice
    for CBE_slice in new_CBE_list_of_slices:
        # retrieved_IS: Covered_bases_ensamble_object for instance (allow tracking), CBE_slice for all attributes (CB in Covered_bases_list are the same object in Covered_bases_ensamble_object! Not copies)
        retrieved_IS = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object, strand_specific = strand_specific_choice)
        #retrieved_IS.Covered_bases_list = CBE_slice.Covered_bases_list
        retrieved_IS.Covered_bases_list = sorted(CBE_slice.Covered_bases_list, key=attrgetter('chromosome', 'locus', 'strand')) # order CB 'along genome'
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
            retrieved_IS.reads_key_list = retrieved_IS.reads_key_list + covered_base.list_of_reads_key
            
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
        #retrieved_IS.Covered_bases_list = CBE_slice.Covered_bases_list
        retrieved_IS.Covered_bases_list = sorted(CBE_slice.Covered_bases_list, key=attrgetter('chromosome', 'locus', 'strand')) # order CB 'along genome'
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
            retrieved_IS.reads_key_list = retrieved_IS.reads_key_list + covered_base.list_of_reads_key
            
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




def dynamic_IS_identification (list_of_Covered_bases_ensambles, ranking_histogram_dict_list, seqTracker_conn_dict, strand_specific_choice, n_parallel_simulations = 2, N_simulations_per_solution = 10):
    '''
    TO DO
    
    SUPER-ALPHA VERSION
    
    NOTE: the aim is to make raw_read_dictionary, final_read_dictionary useless.
    '''

    ###Result collector
    Global_Final_IS_list = []
    
    ### Open a connection to DB through seqTracker_conn_dict
    conn= None
    try:
        conn = DB_connection.dbOpenConnection (seqTracker_conn_dict['host'], seqTracker_conn_dict['user'], seqTracker_conn_dict['passwd'], seqTracker_conn_dict['port'], seqTracker_conn_dict['db'])
    except:
        sys.exit("\n\n\t[ERROR]\tCan't establish a connection with DB anymore.\n\t[QUIT]\n\n")
        
    ### Prepare progressbar
    import progressbar
    global_counter = 0
    global_maxval = len(list_of_Covered_bases_ensambles)
    bar = progressbar.ProgressBar(maxval=global_maxval, widgets=['                        * Ensambles processed: ', progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    
    
    
    # Loop over ensembles
    for Covered_bases_ensamble_object in list_of_Covered_bases_ensambles:
        
        # Strand control
        CBE_strand = Covered_bases_ensamble_object.strand
        if ((CBE_strand != '+')and(CBE_strand != '-')and(CBE_strand != '1')and(CBE_strand != '2')):
            print "\n\tSome troubles found in Dynamic IS retrieval method:"
            print "\tStrand type not recognized. Found: '{0}'.".format(CBE_strand)
            sys.exit("\n\t[ERROR]\tQuit.\n\n")
        
        # Print for Devel
        #print "\t\t{ens_ID} :  ".format(ens_ID=str(Covered_bases_ensamble_object)),
        
        # Check: skip trivial ensambles (1cb or 2adjacent_cb)
        if (Covered_bases_ensamble_object.spanned_bases < 3):
            IS_found = classic(Covered_bases_ensamble_object, strand_specific_choice, center_on_mode=True)
            ### NOTE - in case, fix/set here IS_found/Covered_bases_ensamble_object attributes, according to the general policy of this method (see bottom) ###
            Covered_bases_ensamble_object.flag = "TRIVIAL"
            Covered_bases_ensamble_object.n_of_putative_unique_solution = ""
            IS_found.flag = "TRIVIAL"
            Global_Final_IS_list.append(IS_found)
            # Print for Devel
            #print "TRIVIAL"
            ### Update progressbar
            global_counter += 1
            bar.update(global_counter)
            continue  # 'classic' method is called with the custom mod 'center_on_mode=True' -> return a single IS -> appended to Global_Final_IS_list -> next loop (next Covered_bases_ensamble_object)
        
        # Loop over gauss method conigurations (configDict(s), also called ranking_histogram_dict(s))
        ISs_and_configDict_couple_list = []
        for configDict in ranking_histogram_dict_list:
            ### NOTE - refined_Gaussian_IS_identification also set 'Covered_bases_ensamble_object.IS_derived = Local_Proposed_IS_list': ATTRIBUTE FIXED AFTER 'CHOICE'
            Local_Proposed_IS_list = refined_Gaussian_IS_identification (Covered_bases_ensamble_object, configDict['hist_normalized_to_peak'], configDict['interaction_limit'], strand_specific_choice)
            ISs_and_configDict_couple = (Local_Proposed_IS_list, configDict)
            ISs_and_configDict_couple_list.append(ISs_and_configDict_couple)
        
        # Check: if only one ranking_histogram was tried, the solution is unique... skip!
        if len(ranking_histogram_dict_list) < 2:
            ISs_found = ISs_and_configDict_couple_list[0][0]
            ### NOTE - in case, fix/set here ISs_found/Covered_bases_ensamble_object attributes, according to the general policy of this method (see bottom) ###
            Covered_bases_ensamble_object.flag = "FORCED UNIQUE SOLUTION"
            Covered_bases_ensamble_object.n_of_putative_unique_solution = ""
            for IS in ISs_found:
                IS.flag = "FORCED UNIQUE SOLUTION"
            Global_Final_IS_list = Global_Final_IS_list + ISs_found
            # Print for Devel
            #print "FORCED UNIQUE SOLUTION"
            ### Update progressbar
            global_counter += 1
            bar.update(global_counter)
            continue  # the only list of ISs_found is concatenated to Global_Final_IS_list -> next loop (next Covered_bases_ensamble_object)
        
        # Check: if all different results are consistent... skip!
        consistency = Function_for_Dynamic_IS_identification.check_consistency(ISs_and_configDict_couple_list)
        if (consistency is True):
            # Print for Devel
            #print "consistency test 2 passed -> CONSISTENT, unique solution, no choice needed"
            ISs_found = ISs_and_configDict_couple_list[0][0]
            ### NOTE - in case, fix/set here ISs_found/Covered_bases_ensamble_object attributes, according to the general policy of this method (see bottom) ###
            Covered_bases_ensamble_object.flag = "CONSISTENT"
            Covered_bases_ensamble_object.n_of_putative_unique_solution = ""
            for IS in ISs_found:
                IS.flag = "CONSISTENT"
            Global_Final_IS_list = Global_Final_IS_list + ISs_found
            ### Update progressbar
            global_counter += 1
            bar.update(global_counter)
            continue  # all the IS_lists found are equivalent -> one is concatenated to Global_Final_IS_list -> next loop (next Covered_bases_ensamble_object)
        #else:
            # Print for Devel
            #print "overall consistency test failed, take a choice! "
        
        # HERE:
        #   1) current Covered_bases_ensamble_object is not trivial
        #   2) more than one configDict (ranking_histogram) was tried
        #   3) retrieved solutions are not equivalent and stored in ISs_and_configDict_couple, a list of tuple of kind (Local_Proposed_IS_list, configDict)
        #   --> LET'S CLUSTERIZE THEM
        
        # Solution clustering -> putative_unique_solution_list, a list of 'Putative_unique_solution' objects
        putative_unique_solution_list = Function_for_Dynamic_IS_identification.solution_clustering (ISs_and_configDict_couple_list)
        # Dev Check
        #if len(putative_unique_solution_list) < 2:
            #sys.exit("\n\n\t[ERROR D]\tQuit.\n\n")
        
        # Simulation workflow
        # 0) Retrieve data from DB
        # For each putative_unique_solution:
        #   1) Characterize LTRs and sequences IS by IS (consensus LTR - overall #ins, del, mut - record lenghts)
        #   2) Simulate *many* alternative read realizations of each IS (sequences from refGenome, mapping point and len already defined, scrabling according to #ins, del, mut; then attach consensus LTR)
        #   3) For each alternative realization:
        #       a) Merge sequences in a FASTA file
        #       b) Launch PIPE
        #       c) Retrive data
        #       d) Arrange data in CB -> make a CBE (something like a 'one among the CBEs that maight be, if the reality was the putative solution')
        #       e) Store CBs, CBE and so on in putative_unique_solution as 'SIMULATIONS'
        
        # Assessment step
        # 1) cross-processing among putative_unique_solutions of SIMULATIONS and CONFIG(S)
        # 2) store success/failure rate
        # 3) compute p-value
        # 4) Choose (p-value) the best putative_unique_solution


        ### Analyze sequences ###
        header_list = Covered_bases_ensamble_object.get_headers()
        dictionary_for_sequence_simulations, LTR_LC_dictionary = Function_for_Dynamic_IS_identification.analyze_sequences (header_list, conn, seqTracker_conn_dict)
        # Split LTR_LC_dictionary strand-wise
        LTR_LC_dictionary_plus = {}
        LTR_LC_dictionary_minus = {}
        for header, sub_dict in LTR_LC_dictionary.items():
            if ((sub_dict['strand'] == '+') or (sub_dict['strand'] == '1')):
                LTR_LC_dictionary_plus[header] = sub_dict
            elif ((sub_dict['strand'] == '-') or (sub_dict['strand'] == '2')):
                LTR_LC_dictionary_minus[header] = sub_dict
        
        ### SIMULATIONS ###
        for putative_unique_solution_object in putative_unique_solution_list:
            
            ### Prepare simulations ### --> Fill putative_unique_solution_object attributes: perfect_sequence_dict, perfect_sequence_strandness_dict, seq_MID_dict_list
            
            ### Add perfect_sequence_dict attribute to putative_unique_solution objects : {'header': sequence}
            ### Add perfect_sequence_strandness_dict attribute to putative_unique_solution objects : {'header': strand}
            Function_for_Dynamic_IS_identification.get_seq_from_ref (putative_unique_solution_object, dictionary_for_sequence_simulations, reference_genome='hg19')
            ### Add seq_MID_dict_list simulated_sequence_dict, a list paired with putative_unique_solution_object.IS_list like [{'M':numM, 'I':numI, 'D':numD}, {...}, ... ]
            Function_for_Dynamic_IS_identification.get_seq_MID_dict_list (putative_unique_solution_object, dictionary_for_sequence_simulations)
            
            ### Simulate! ### --> Fill putative_unique_solution_object attribute: simulated_sequence_dict_list
            
            ### Organize 'parallelized_simulations' calls
            n_loop = N_simulations_per_solution / n_parallel_simulations
            reminder = N_simulations_per_solution - (n_loop*n_parallel_simulations)
            ### Main loop
            for n in range(n_loop):
                #Append simulation(s) to simulated_sequence_dict_list attribute of putative_unique_solution objects : append({'header': sequence})
                Function_for_Dynamic_IS_identification.parallelized_simulations (putative_unique_solution_object, LTR_LC_dictionary_plus, LTR_LC_dictionary_minus, n_parallel_simulations)
            ### Extra call if n_parallel_simulations can't exactly divide N_simulations_per_solution
            if reminder != 0:
                #Append simulation(s) to simulated_sequence_dict_list attribute of putative_unique_solution objects : append({'header': sequence})
                Function_for_Dynamic_IS_identification.parallelized_simulations (putative_unique_solution_object, LTR_LC_dictionary_plus, LTR_LC_dictionary_minus, reminder)
                
            
        ### Fake choice just to conclude - take the putative_unique_solution with the highest cardinality ###
        Local_Selected_IS_list = None
        max_cardinality = 0
        for putative_unique_solution in putative_unique_solution_list:
            if putative_unique_solution.cardinality > max_cardinality:
                Local_Selected_IS_list = putative_unique_solution.IS_list
        ### Join Local_Selected_IS_list with Global_Final_IS_list ###
        for IS in Local_Selected_IS_list:
            IS.flag = "SIMULATIONS"
        Global_Final_IS_list = Global_Final_IS_list + Local_Selected_IS_list
        ### After choice, fix attributes '.IS_derived' for and '.flag' for current Covered_bases_ensamble_object
        Covered_bases_ensamble_object.IS_derived = Local_Selected_IS_list
        Covered_bases_ensamble_object.flag = "SIMULATIONS"
        Covered_bases_ensamble_object.n_of_putative_unique_solution = len(putative_unique_solution_list)
        
        ### Update progressbar
        global_counter += 1
        bar.update(global_counter)
    
    
    
    ### Close DB connection
    DB_connection.dbCloseConnection (conn)
    ### Close progressbar
    bar.finish()
    
    ### Return Result
    return Global_Final_IS_list
