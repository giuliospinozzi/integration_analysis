###Header################################################
header = """

+------------------------------------------------------+
 Module: Function_for_Gaussian_IS_identification
 Author: Stefano Brasca
 Date:  November 27th, 2013
 Contact: brasca.stefano@hsr.it
 Version: 1.0
+------------------------------------------------------+

 Description:
  - This module contains functions used in Gaussian
    IS identification framework
    (refined_Gaussian_IS_identification function and
    Gaussian_IS_identification [deprecated] in
    Integration_Sites_retrieving_methods module)
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import numpy
import copy
###############################

###Import Module(s)#####################
import Classes_for_Integration_Analysis
########################################





def gaussian_histogram_generator (interaction_limit, alpha):
    '''
    *** This function generates a Gaussian-shaped histogram ***
    
    INPUT: interaction_limit - Integer number, forced to be (read description and note CAREFULLY)
           alpha - Float number  (read description and note CAREFULLY)
           
    OUTPUT: bin_boundaries - list of tuples of two elements, setting bins boundaries
            bin_areas - list of float, setting the area of each bin
            diagnostic - float number, contains the fraction of distribution lost (area of cutted tails)
    
            
    DESCRIPTION AND NOTES:
    
    # INTERACTION_LIMIT
    interaction_limit states, de facto, the number of bin of the histogam you get (n_bins = (2*interaction_limit)+1)
    
    examples about interaction_limit:
    interaction_limit = 3 means you'll get an histogram having 7 bins (central one flanked by 3 on each side)
    interaction_limit = 4 means you'll get an histogram having 9 bins (central one flanked by 4 on each side)
    
    # ALPHA
    alpha states HOW MANY SIGMAS are equal to HALF-BASEPAIR
    
    examples about alpha:
    alpha = 1 means that sigma is half-bp long; then 3bp are long 6sigma
    alpha = 0.5 means that sigma is 1-bp long
    
    
    # ! PLEASE NOTE:
    Since we assume the peak has no interaction with bases more than interaction_limit-bp far, you got an histogram having 2*interaction_limit+1 bins
    under the implicit assumption that the excluded area is "negligible". 
    BE AWARE THAT IT'S THE FOLLOWING CHOICE OF ALPHA THAT STASES IF THIS AREA IS REALLY NEGLIGIBLE OR NOT (bigger alphas improve negligibleness).
    The idea is that you should choose a couple that excludes only portions at-least-3-sigma-far from the central bin,
    thus neglecting 0.3% of the area or less. Conversely, keep in mind that if alpha is too big, the histogram you'll get will be very peaky (lower alphas 
    improve gaussian shape).
    You can try to make some guesses but verify them through 'diagnostic' variable that contains the fraction of distribution lost; don't try tricky choices
    and keep in mind that the calculation are approximate, not only because of floating algebra, but also due to Formula (tends to overstimate tails of density)
    [High Accurate Simple Approximation ofNormal Distribution Integral, Hector Vazquez-Leal et al., doi:10.1155/2012/124029]
    
    
    # SOME GUESSES
    
    # super-stringent: low range of interaction, high precision; 7bins and 14sigmas of span (7 each side)
    alpha = 1
    interaction_limit = 3
    
    # an alternative, less stringent (in term of interaction limit and tail loss - nevertheless acceptable) but better shaped;
    # 9 bins and 9 sigmas of span (4.5 each side)
    alpha = 0.5
    interaction_limit = 4
    
    # maybe a compromise (we loose 0.002% of distribution - 2 orders of magnitude smaller than boundaries bins)
    alpha = 0.6 #(good gaussian shape)
    interaction_limit = 3 #(as bushman said)
    
    # validated in silico on simulated data (optimal choice):
    alpha = 0.3
    interaction_limit = 4
    
    '''
    
        # Function for normal  CDF ######################################################################
    def norm_CDF_from_x_to_y (x, y, alpha):
        # Given alpha (the number of sigmas equals to half-bp)
        # and given x and y (coordinates in bp unit), this function
        # provides the normal CDF from x to y, thanks to an approximate Formula
        # Normal distribution is meant to be located in zero with a sigma of 1/(2alpha) bp
    
        # Sigma
        sigma = 1.0 / numpy.sqrt(2.0*numpy.pi)
        alpha = float(alpha)
    
        ### Approximate Formula for CDF of normal distribution (0, (2pi)^-1/2), from zero to x ###
        a = 39.0/2.0
        b = 111.0/2.0
        c = 35.0/111.0
        def Formula (x, a, b, c):
            result = 0.5*numpy.tanh((a*x)-(b*numpy.arctan(c*x)))
            return result
        ##########################################################################################
    
        CDF_from_0_to_x = Formula(x*2*alpha*sigma, a, b, c) # 2 because alpha states the relation
        CDF_from_0_to_y = Formula(y*2*alpha*sigma, a, b, c) # between sigma and HALF-base
    
        result = CDF_from_0_to_y - CDF_from_0_to_x
        return result
        ###############################################################################################
    
    #Cast
    interaction_limit = int(interaction_limit)
    
    # Bin Boundaries Half
    first_bin_half = (0.0, 0.5)
    bin_boundaries_half = [first_bin_half]
    
    for i in range(0, interaction_limit):
        bin_boundaries_half.append((bin_boundaries_half[-1][1], bin_boundaries_half[-1][1]+1))
    
    # Filling Bins
    bin_areas_half = []
    for x,y in bin_boundaries_half:
        bin_areas_half.append(norm_CDF_from_x_to_y (x, y, alpha))
    
    # Symmetrize Bins
    n_bins = int((2*interaction_limit)+1)
    bin_boundaries = [None]*n_bins
    bin_areas = [None]*n_bins
    
    bin_boundaries[interaction_limit] = (-0.5, 0.5)
    bin_areas[interaction_limit] = bin_areas_half[0] * 2.0
    for i in range(1, interaction_limit +1):
        negative_left = -1.0*bin_boundaries_half[interaction_limit+1-i][1]
        negative_right = -1.0*bin_boundaries_half[interaction_limit+1-i][0]
        bin_boundaries[i-1] = (negative_left, negative_right)
        bin_boundaries[i+interaction_limit] = bin_boundaries_half[i]
        bin_areas[i-1] = bin_areas_half[interaction_limit+1-i]
        bin_areas[i+interaction_limit] = bin_areas_half[i]
        
    # Remember that bin_boundaries[interaction_limit] is the peak
    
    # Diagnostic: fraction of CDF lost
    diagnostic = 1.0 - sum(bin_areas)
    
    # Return Results
    return bin_boundaries, bin_areas, diagnostic





def CBE__histogram_generator (Covered_bases_ensamble):
    '''
    *** From a covered bases ensemble objects, it creates an 'histogram' ***
    
    INPUT: Covered_bases_ensamble object. Not modified during the process
    
    OUTPUT: bin_areas - List of Float (ReadCounts), representing the histogram
                        of the Covered_bases_ensamble given in input (bin heights)
            list_of_loci - list of int, reporting the loci for bin_areas
            max_read_count - Float, the highest ReadCount in bin_areas
            index_of_max - Integer, index of max_read_count in bin_areas
            
    NOTE for developers: 
          Code was wrote this way in order to limit usage of Covered_bases_ensamble
          attributes (just 'spanned_bases' and 'Covered_bases_list'): this way you
          can get a usable and reliable histogram even when some covered_bases have
          been removed not properly (acting directly on Covered_bases_list).
          Forthcoming improvements aim at make this approach useless (push_out
          method for Covered_bases_ensamble objects and development of new function
          'refined_Gaussian_IS_identification') but old Gaussian_IS_identification
          method [deprecated and not in use anymore] needs this code
    '''
    n_bins = Covered_bases_ensamble.spanned_bases
    bin_areas = [0.0]*n_bins
    
    list_of_loci =[]
    for Covered_base in Covered_bases_ensamble.Covered_bases_list:
        list_of_loci.append(Covered_base.locus)
    locus_min = min(list_of_loci)
    del list_of_loci
    
    list_of_loci =[]        
    for i in range(0, n_bins):
        current_locus = locus_min + i
        list_of_loci.append(current_locus)
        for Covered_base in Covered_bases_ensamble.Covered_bases_list:
            if (Covered_base.locus == current_locus):
                bin_areas[i] = float(Covered_base.reads_count)
                
    # Now bin_areas is a list of number, reproducing the shape of reads_count histogram for the Covered_bases_ensamble given in input
    
    max_read_count = max(bin_areas)
    index_of_max = bin_areas.index(max_read_count)
    
    # Return Results
    return bin_areas, list_of_loci, max_read_count, index_of_max





def normalize_histogram_to_the_peak (bin_areas, index_of_max):
    '''
    *** 'Normalize' an histogram to make peak's area = 1 ***
    '''
    index_of_max = int(index_of_max)
    max_height = bin_areas[index_of_max]
    bin_areas_normalized = []
    for one_bin in bin_areas:
        bin_areas_normalized.append(one_bin/max_height)
        
    return bin_areas_normalized





def explore_and_split_CBE (Covered_bases_ensamble_object, strand_specific_choice):
    '''
    *** Explore a Covered Base Ensemble in order to split it in 'peaks + vicinity' ***
               
    INPUT: - Covered_bases_ensamble_object: no further specifications needed. Not modified during the process
           - strand_specific_choice: boolean; it specifies if the matrix computation algorithm had to account for strand: generally, the choice made here should
                                     reflect the ones previously made for Covered_bases_ensamble_object. For this purpose, you can find a variable called 
                                     'strand_specific_choice' in main, retrieved from user input, so the best usage is strand_specific_choice = strand_specific_choice
                                     
    OUTPUT: CBE_list_of_slices - list of Covered Base Ensemble objects, ordered (descending) by peak's height; objects in this list constitute a partition of the
                                 Covered_bases_ensamble_object given in input (and covered bases in each 'Covered_bases_list' attribute are also the very ones from
                                 Covered_bases_ensamble_object given in input)
                                 
    LOGIC: Starting from the covered base with highest reads count (highest 'peak') in Covered_bases_ensamble_object given in input, a new Covered Base Ensemble
           is instanced; then, if adjacent loci host covered bases they are 'push(ed)_in'. This product (current_CBE_slice) is appended to CBE_list_of_slices, 
           then removed from Covered_bases_ensamble_object and this logic is ready to be applied again ( while (bases_to_assign > 0) ).
           - NOTE for Developers: Covered_bases_ensamble_object IS NOT REALLY MODIFIED during this process
    '''
    # Covered_bases_ensamble's slices list
    CBE_list_of_slices =[]
    
    #Order Covered_bases_ensamble_object.Covered_bases_list by reads_count (Descending order)
    ordered_covered_bases_list = sorted(Covered_bases_ensamble_object.Covered_bases_list, key=lambda x: x.reads_count, reverse=True)
    
    #Splitting loop
    bases_to_assign = Covered_bases_ensamble_object.n_covered_bases
    while (bases_to_assign > 0):
        current_CBE_slice = Classes_for_Integration_Analysis.Covered_bases_ensamble(ordered_covered_bases_list[0], strand_specific=strand_specific_choice)
        bases_to_assign = bases_to_assign - 1
        ordered_covered_bases_list.remove(ordered_covered_bases_list[0])
        # Check peak's sides
        covered_base_to_remove = []
        for covered_base in ordered_covered_bases_list:
            if ((covered_base.locus == current_CBE_slice.Covered_bases_list[0].locus + 1) or (covered_base.locus == current_CBE_slice.Covered_bases_list[0].locus - 1)):
                current_CBE_slice.push_in(covered_base)
                covered_base_to_remove.append(covered_base)
        for covered_base in covered_base_to_remove:
            bases_to_assign = bases_to_assign - 1
            ordered_covered_bases_list.remove(covered_base)
        # Append current_CBE_slice to CBE_list_of_slices
        CBE_list_of_slices.append(current_CBE_slice)
        
    return CBE_list_of_slices # CBE_list_of_slices is ordered by peak's height; items are CBE object created from CB object coming from Covered_bases_ensamble_object in input




def evaluate_surroundings (CBE_slice, whole_CBE, hist_gauss_normalized_to_peak, interaction_limit):
    '''
    *** Given a CBE_slice, it returns a score dictionary about loci in whole_CBE ***
               [Designed to be used in global_score_dictionary function,
                      input from explore_and_split_CBE function]
                      
    INPUT: - CBE_slice: Covered Base Ensemble object, typically an element of CBE_list_of_slices from explore_and_split_CBE function
           - whole_CBE: Covered Base Ensemble object from which CBE_slice is derived. Not modified during the process
           - hist_gauss_normalized_to_peak: a list of float representing the histogram of the discrete-gaussian used as model (bin heights),
                                            'normalized' to make central bin = 1.
                                            Typically:
                                            bin_boundaries, bin_areas, diagnostic = Function_for_Gaussian_IS_identification.gaussian_histogram_generator(interaction_limit, alpha)
                                            hist_gauss_normalized_to_peak = Function_for_Gaussian_IS_identification.normalize_histogram_to_the_peak(bin_areas, interaction_limit)
           - interaction_limit: int, see gaussian_histogram_generator function for further details
    
    OUTPUT: score_dic - dictionary of kind: {'locus1':(CBE_slice given in input, evaluation score), 'locus2':(CBE_slice given in input, evaluation score), ...}
                        Comment about score: evaluation score is a real number, always comparable with other score of dictionary like this. It allows to states which CBE_slice has more
                                             influence over each locus (once such a dictionary is available for each CBE_slice derived from whole_CBE). A non-negative score asserts that CBE_slice
                                             has influence over related locus, the higher the more.
                        Comment about evaluation: loci belonging to CBE_slice and loci whose bins are 0 (see (*) comment in code) are not evaluated
    
    TYPICAL USAGE IN THIS CONTEXT: CBE_slice has to be the slice of whole_CBE with the highest peak. After first use, if you want to loop, you have to remove CBE_slice content from whole_CBE
                                   then pass to the second CBE_slice in terms of peak height 
    '''

    # Build histogram for whole_CBE
    whole_CBE_bin_areas, whole_CBE_list_of_loci, whole_CBE_max_read_count, whole_CBE_index_of_max = CBE__histogram_generator(whole_CBE)
    whole_CBE_bin_areas_normalized = normalize_histogram_to_the_peak(whole_CBE_bin_areas, whole_CBE_index_of_max)
    del whole_CBE_max_read_count
    
    # n of allowed step left and right from whole_CBE's peak
    index_last_bin = len(whole_CBE_bin_areas) - 1
    n_step_right = index_last_bin - whole_CBE_index_of_max
    if (n_step_right > interaction_limit):
        n_step_right = interaction_limit
    n_step_left = whole_CBE_index_of_max
    if (n_step_left > interaction_limit):
        n_step_left = interaction_limit
        
    # starting and ending indexes for hist_gauss
    starting_index = interaction_limit - n_step_left
    ending_index = interaction_limit + n_step_right
    
    # list of allowed indexes hist_gauss
    allowed_indexes_gauss = range(starting_index, ending_index+1)
    
    # starting and ending indexes for current_ensemble_bin_areas_normalized
    starting_index = whole_CBE_index_of_max - n_step_left
    ending_index = whole_CBE_index_of_max + n_step_right
    
    # list of allowed indexes for current_ensemble_bin_areas_normalized
    allowed_indexes_CBE = range(starting_index, ending_index+1)
    
    # list of indexes tuples [(allowed_indexes_gauss1, allowed_indexes_CBE1), (allowed_indexes_gauss2, allowed_indexes_CBE2), ... ]
    # indexes of peak's locus and of loci beside peak will be excluded
    indexes_tuples = []
    for i in range(0, n_step_left+n_step_right+1):
        if ((allowed_indexes_CBE[i] != whole_CBE_index_of_max) and (allowed_indexes_CBE[i] != whole_CBE_index_of_max + 1) and (allowed_indexes_CBE[i] != whole_CBE_index_of_max - 1)):
            indexes_tuples.append((allowed_indexes_gauss[i],allowed_indexes_CBE[i]))
    
    score_dic = {} #dictionary of kind: {locus:(CBE_slice, score)}         
    for i,j in indexes_tuples:
        if (whole_CBE_bin_areas_normalized[j] != 0): # Added if clause, it seems to fix the whole algorithm!! (*)
            score = hist_gauss_normalized_to_peak[i] - whole_CBE_bin_areas_normalized[j]
            score_dic.update({whole_CBE_list_of_loci[j]:(CBE_slice, score)})
            # (*) indeed, is not useful to evaluate loci whose bins are 0: such loci in fact are really empty or at least
            # previously occupied by an higher peak + its adjacent cb!
        
    return score_dic 
    #this dictionary contains each locus evaluate with respect to CBE_slice as key and as value shows a tuple of kind (which CBE slice gave the mark to that locus, how much the mark is) 
    #about mark: is a real number, positive if it's OK (the higher the better) and negative if it's not (the lower, the worse)

    
        

def global_score_dictionary (CBE_list_of_slices, whole_CBE, hist_gauss_normalized_to_peak, interaction_limit, strand_specific_choice):
    '''
    *** it returns a global score dictionary about loci in whole_CBE and each CBE slice in CBE_list_of_slices ***
     [Designed as a smart 'looping box' for evaluate_surroundings, helping with correct usage and showing merged 
                                       results in a unique dictionary]
                                       
    INPUT: see evaluate_surroundings function above
    
    OUTPUT: global_score_dic - a dictionary of kind: {locus:[(CBE_slice, score), (...), ...]}, obtained 'merging' score_dic dictionary
                               returned by evaluate_surroundings function callings (if the key is the same, items (tuple) are joined in a list 
                                       
    LOGIC: 
    It takes in input CBE_list_of_slices coming from explore_and_split_CBE function and Covered_bases_ensamble_object from which CBE
    slices are derived (here whole_CBE), then it makes a loop usage (over CBE_slice in CBE_list_of_slices) of evaluate_surroundings function:
    single score_dicS returned at each cycle are merged in global_score_dic. A deeper insight is derivable from description of 
    evaluate_surroundings function and a quick look to the following code.
    Note that whole_CBE / Covered_bases_ensamble_object is not modified during the process (the working copy is created inside this function then discarded)
    '''
    
    current_ensemble = copy.deepcopy(whole_CBE)
    
    global_score_dic = {}
    for CBE_slice in CBE_list_of_slices:
        
        current_dic = evaluate_surroundings (CBE_slice, current_ensemble, hist_gauss_normalized_to_peak, interaction_limit)
        
        for key, item in current_dic.iteritems():
            if global_score_dic.has_key(key):
                global_score_dic[key].append(current_dic[key])
            else:
                global_score_dic.update({key:[item]})
                
        list_of_locus_to_remove = []
        for covered_base in CBE_slice.Covered_bases_list:
            list_of_locus_to_remove.append(covered_base.locus)
        
        new_Covered_base_list = []    
        for covered_base in current_ensemble.Covered_bases_list:
            if (covered_base.locus not in list_of_locus_to_remove):
                new_Covered_base_list.append(covered_base)
        
        if (len(new_Covered_base_list) != 0):        
            current_ensemble = Classes_for_Integration_Analysis.Covered_bases_ensamble(new_Covered_base_list[0], strand_specific=strand_specific_choice)
            for covered_base in new_Covered_base_list[1:]:
                current_ensemble.push_in(covered_base)
            
    # global_score_dic is ready
    
    return global_score_dic #dictionary of kind: {locus:[(CBE_slice, score), (...), ...]}




def reconstruct_CBE_slice (CBE_slice, list_of_bases_to_assign):
    '''
    *** Given a CBE_slice, it returns a new one with some bases more, the ones in list_of_bases_to_assign ***
    
    INPUT: CBE_slice - Covered Base Ensemble object, typically an element of CBE_list_of_slices from explore_and_split_CBE function
           list_of_bases_to_assign - a list of kind: [(covered_base (object) to assign, CBE_slice claiming it), (...), ...]
    
    OUTPUT: new_CBE_slice - a new CBE slice containing all the CB in CBE_slice given in input, plus the ones in list_of_bases_to_assign
                            for which there is a CBE_slice matching
    '''
        
    list_of_already_present_bases = []
    for CB in CBE_slice.Covered_bases_list:
        list_of_already_present_bases.append(CB)
        
    list_of_bases_just_assigned = []
    for CB, claiming_CBE_slice in list_of_bases_to_assign:
        if (claiming_CBE_slice.Covered_bases_list == CBE_slice.Covered_bases_list):
            list_of_bases_just_assigned.append(CB)
        
    new_list_of_CB = list_of_already_present_bases + list_of_bases_just_assigned
        
    # Create a new CBE slice
    new_CBE_slice = Classes_for_Integration_Analysis.Covered_bases_ensamble(new_list_of_CB[0])
    for CB in new_list_of_CB[1:]:
        new_CBE_slice.push_in(CB)
       
    return new_CBE_slice    
    
    
    
        