###Header###################################################
header = """

+--------------------------------------------------------+
 Module: Function_for_SkewedGaussian_IS_identification
 Author: Stefano Brasca
 Date:  March 4th, 2013
 Contact: brasca.stefano@hsr.it
 Version: 1.0
+--------------------------------------------------------+

 Description:
  - This module contains functions used in SkewedGaussian
    IS identification framework
    ('refined_SKEWED_Gaussian_IS_identification' in
    Integration_Sites_retrieving_methods module)
  
 Note: [...]

+--------------------------------------------------------+ 
""" 
###########################################################

###Requested Package(s) Import##################################################
from math import sqrt, copysign, pi
from numpy import where, zeros, ones, int64, float64, array, arange, inner, kron
from numpy import exp as np_exp
from numpy import arctan as np_arctan
from scipy.stats import norm
from scipy.special import gamma as sp_gamma
import copy
import sys
################################################################################

###Import Module(s)#####################
import Classes_for_Integration_Analysis
########################################





def SKEWED_gaussian_histogram_generator (interaction_limit, location, scale, shape):
    
    ######################################################################
    ### shape MUST BE ALWAYS NEGATIVE there, for correct hist building ###
    ######################################################################
    
    '''
    *** This function generates a SkewedGaussian-shaped histogram ***
     [info at http://en.wikipedia.org/wiki/Skew_normal_distribution]
     [further stuff at http://azzalini.stat.unipd.it/SN/index.html - code too]
     [R-based alternative at http://cran.r-project.org/web/packages/sn/index.html]
     
    INPUT: interaction_limit - Integer number, better if even
           location - Integer number: MUST BE '0.0' , to let the whole function work properly.
           scale - Float positive number, any.
           shape - Float number: MUST BE GIVEN ALWAYS NEGATIVE, to let the whole function work properly.
           
    OUTPUT: bin_boundaries - list of tuples of two elements, setting bins boundaries
            bin_areas - list of float, setting the area of each bin
            diagnostic - float number, contains the fraction of distribution lost (area of cut tails)
    
            
    DESCRIPTION AND NOTES:
    
    # INTERACTION_LIMIT
    interaction_limit states, de facto, the number of bin of the histogram you get. 1 bin is always occupied by
    the peak, then interaction_limit/2 (or interaction_limit+1/, if odd) bins span the shorter tail of the 
    distribution and three times more span the longest.
    e.g. if interaction_limit is even, you get 2*interaction_limit + 1 bins
    
    # LOCATION, SCALE, SHAPE
    see http://en.wikipedia.org/wiki/Skew_normal_distribution (xi, omega and alpha, respectively)
    In brief, they are analogous of mean, variance and skewness. If you prefer to give the very mean, variance and skewness
    as input, some function are available below, commented. 
    Note that:
    - skewness has an upper limit value (and has to be given negative, like 'shape')
    - location must turn out to be 0.0 in any case!
    
    # ! PLEASE NOTE: 
    It being understood that location = 0.0 is required and shape has to be < 0, it's difficult to find couples of
    scale-shape suitable for all the purposes, such as:
    - minimize diagnostic
    - produce a well-shaped histogram
    - get a final histogram that works good with data
    The idea for the future is to implement a function that can handle this task. By now,
    we empirically found that 
    "SKEWED_gaussian_histogram_generator (interaction_limit, location=0.0, scale=3.0, shape=-4.0)"
    is a quite good solution.     
    '''
    
    ## FUNCTIONS FOR SKEWED GAUSSIAN ####################################################
    
    # Maybe useful in the future
#===============================================================================
#     def skew_max():
#         """
#         Return maximum skewness of a skew normal distribution
#         """
#         beta = 2.0 - pi / 2.0
#         #lim(delta, shape-> inf) = 1.0
#         eps = sqrt(2.0 / pi)
#         return beta * pow(eps, 3.0) / pow(1.0 - eps * eps, 3.0 / 2.0) - 1e-16
# 
#     def skewnormal_parms(mean=0.0, stdev=1.0, skew=0.0):
#         """
#         Return parameters for a skew normal distribution function:
#             location (xi), scale (omega) and shape (alpha)
#         """
#         if abs(skew) > skew_max():
#             print('Skewness must be between %.8f and %.8f' % (-skew_max(), skew_max()))
#             print('None, None, None returned')
#             return None, None, None
#         beta = (2.0 - pi / 2.0)
#         skew_23 = pow(skew * skew, 1.0 / 3.0)
#         beta_23 = pow(beta * beta, 1.0 / 3.0)
#         eps2 = skew_23 / (skew_23 + beta_23)
#         eps = copysign(sqrt(eps2), skew)
#         delta = eps * sqrt(pi / 2.0)
#         alpha = delta / sqrt(1.0 - delta * delta)
#         omega = stdev / sqrt(1.0 - eps * eps)
#         xi = mean - omega * eps
#         return xi, omega, alpha
#===============================================================================
    
    def fui(h, i):
        return (h ** (2 * i)) / ((2 ** i) * sp_gamma(i + 1))
    
    def T_Owen_int(h, a, jmax=50, cut_point=6):
        """
        Return Owens T
        """
        if type(h) in (float, float64):
            h = array([h])
        low = where(h <= cut_point)[0]
        high = where(h > cut_point)[0]
        n_low = low.size
        n_high = high.size
        irange = arange(0, jmax, dtype=int64)
        series = zeros(h.size)
        if n_low > 0:
            h_low = h[low].reshape(n_low, 1)
            b = fui(h_low, irange)
            cumb = b.cumsum(axis=1)
            b1 = np_exp(-0.5 * h_low ** 2) * cumb
            matr = ones((jmax, n_low)) - b1.transpose()
            jk = kron(ones(jmax), [1.0, -1.0])
            jk = jk[0: jmax] / (2 * irange + 1)
            matr = inner((jk.reshape(jmax, 1) * matr).transpose(),
                         a ** (2 * irange + 1))
            series[low] = (np_arctan(a) - matr.flatten(1)) / (2 * pi)
        if n_high > 0:
            h_high = h[high]
            atana = np_arctan(a)
            series[high] = (atana * np_exp(-0.5 * (h_high ** 2) * a / atana) *
                        (1.0 + 0.00868 * (h_high ** 4) * a ** 4) / (2.0 * pi))
        return series
       
    def T_Owen_series(h, a, jmax=50, cut_point=6):
        """
        Return Owens T
        """
        if abs(a) <= 1.0:
            return T_Owen_int(h, a, jmax=jmax, cut_point=cut_point)
        else:
            """D.B. Owen Ann. Math. Stat. Vol 27, #4 (1956), 1075-1090
             eqn 2.3, 2.4 and 2.5
             Available at: http://projecteuclid.org/DPubS/Repository/1.0/
                Disseminate?view=body&id=pdf_1&handle=euclid.aoms/1177728074"""
            signt = copysign(1.0, a)
            a = abs(a)
            h = abs(h)
            ha = a * h
            gh = norm.cdf(h)
            gha = norm.cdf(ha)
            t = 0.5 * gh + 0.5 * gha - gh * gha - T_Owen_int(ha, 1.0 / a, jmax=jmax, cut_point=cut_point)
            return signt * t
       
    def cdf_skewnormal(x, location=0.0, scale=1.0, shape=0.0):
        """
        Return skew normal cdf
        """
        xi = (x - location) / scale
        return norm.cdf(xi) - 2.0 * T_Owen_series(xi, shape)
    

    ## FUNCTIONS FOR HISTOGRAM ###########################################################
    def create_bin_boundaries (interaction_limit):        
        nBin_right = None
        
        if ((interaction_limit/2) == (float(interaction_limit)/2.0)):
            nBin_right = interaction_limit/2
        else:
            nBin_right = interaction_limit/2 + 1
    
        nBin_left = 3*nBin_right
    
        nBin_tot = nBin_left+1+nBin_right
    
        bin_boundaries = [None]*nBin_tot
    
        bin_boundaries[nBin_left] = (-0.5,0.5) #bin of max
    
        for i in range(0,nBin_left):
            bin_boundaries[nBin_left-i-1] = (bin_boundaries[nBin_left][0]-i-1,bin_boundaries[nBin_left][1]-i-1)
        for i in range(0,nBin_right):
            bin_boundaries[nBin_left+i+1] = (bin_boundaries[nBin_left][0]+i+1,bin_boundaries[nBin_left][1]+i+1)
    
        return bin_boundaries
    
    def create_bin_areas (bin_boundaries, location, scale, shape):
        bin_boundaries_extended = []
    
        for i in reversed(range(1,(len(bin_boundaries)-1)/2 +1)):
            tupla = (bin_boundaries[0][0]-i,bin_boundaries[0][1]-i)
            bin_boundaries_extended.append(tupla)
        bin_boundaries_extended = bin_boundaries_extended + bin_boundaries
        for i in range(1,(len(bin_boundaries)-1)/4 +1):
            tupla = (bin_boundaries[-1][0]+i,bin_boundaries[-1][1]+i)
            bin_boundaries_extended.append(tupla)
    
        bin_areas_extendex = []
    
        for i,j in bin_boundaries_extended:
            bin_areas_extendex.append(float(cdf_skewnormal(j, location=location, scale=scale, shape=shape) - cdf_skewnormal(i, location=location, scale=scale, shape=shape)))
    
        index_of_max = bin_areas_extendex.index(max(bin_areas_extendex))
    
        bin_areas = []
    
        for i in range(index_of_max-(len(bin_boundaries)-1)/4*3,index_of_max):
            bin_areas.append(bin_areas_extendex[i])
        bin_areas.append(bin_areas_extendex[index_of_max])
        for i in range(index_of_max+1, index_of_max+(len(bin_boundaries)-1)/4 +1):
            bin_areas.append(bin_areas_extendex[i])
    
        #print bin_boundaries_extended
        #print bin_areas_extendex
    
        return bin_areas
    #####################################################################################
    
    bin_boundaries = create_bin_boundaries (interaction_limit)
    bin_areas = create_bin_areas (bin_boundaries, location, scale, shape)
    diagnostic = 1 - sum(bin_areas)
    
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
          been removed not properly (e.g. acting directly on Covered_bases_list).
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
           is instanced; then, if adjacent locus (left/right according to strand) host a covered base, it's 'push(ed)_in'. This product (current_CBE_slice) is appended to CBE_list_of_slices, 
           then removed from Covered_bases_ensamble_object and this logic is ready to be applied again ( while (bases_to_assign > 0) ).
           - NOTE for Developers: Covered_bases_ensamble_object IS NOT ACTUALLY MODIFIED during this process
    '''
    # Covered_bases_ensamble's slices list
    CBE_list_of_slices =[]
    
    # CBE strand
    CBE_strand = str(Covered_bases_ensamble_object.strand) # Cast should be redundant
       
    #Order Covered_bases_ensamble_object.Covered_bases_list by reads_count (Descending order) and, second, by locus (Ascend/descend according to strand)
    choice = None
    if ((CBE_strand == '+')or(CBE_strand == '1')):
        choice = False
    else:
        choice = True    
    # Sorting made in 2 step (2 then 1) beacause of changes in reverse policy according to strand
    semi_ordered_covered_bases_list = sorted(Covered_bases_ensamble_object.Covered_bases_list, key=lambda x: x.locus, reverse=choice)    
    ordered_covered_bases_list = sorted(semi_ordered_covered_bases_list, key=lambda x: x.reads_count, reverse=True)
    
    
    #Splitting loop
    bases_to_assign = Covered_bases_ensamble_object.n_covered_bases
    while (bases_to_assign > 0):
        current_CBE_slice = Classes_for_Integration_Analysis.Covered_bases_ensamble(ordered_covered_bases_list[0], strand_specific=strand_specific_choice)
        bases_to_assign = bases_to_assign - 1
        ordered_covered_bases_list.remove(ordered_covered_bases_list[0])
        # Check peak's sides
        covered_base_to_remove = []
        for covered_base in ordered_covered_bases_list:
            # Following condition modified to take in account orientation
            if ((((CBE_strand == '+')or(CBE_strand == '1'))and(covered_base.locus == current_CBE_slice.Covered_bases_list[0].locus + 1)) or (((CBE_strand == '-')or(CBE_strand == '2'))and(covered_base.locus == current_CBE_slice.Covered_bases_list[0].locus - 1))):
                current_CBE_slice.push_in(covered_base)
                covered_base_to_remove.append(covered_base)
        for covered_base in covered_base_to_remove:
            bases_to_assign = bases_to_assign - 1
            ordered_covered_bases_list.remove(covered_base)
        # Append current_CBE_slice to CBE_list_of_slices
        CBE_list_of_slices.append(current_CBE_slice)
        
    return CBE_list_of_slices # CBE_list_of_slices is ordered by peak's height; items are CBE object created from CB object coming from Covered_bases_ensamble_object in input




def evaluate_surroundings (CBE_slice, whole_CBE, hist_gauss_normalized_to_peak):
    '''
    *** Given a CBE_slice, it returns a score dictionary about loci in whole_CBE ***
               [Designed to be used in global_score_dictionary function,
                      input from explore_and_split_CBE function]
                      
    INPUT: - CBE_slice: Covered Base Ensemble object, typically an element of CBE_list_of_slices from explore_and_split_CBE function
           - whole_CBE: Covered Base Ensemble object from which CBE_slice is derived. Not modified during the process
           - hist_gauss_normalized_to_peak: a list of float representing the histogram of the skewed-gaussian used as model (list of bin heights), 'normalized' to make peak's bin = 1.
                                            NOTE: IT'S UP TO YOU PASSING THE CORRECT HISTOGRAM, ACCORDING TO CBE's STRAND.
    
    OUTPUT: score_dic - dictionary of kind: {'locus1':(CBE_slice given in input, evaluation score), 'locus2':(CBE_slice given in input, evaluation score), ...}
                        Comment about score: evaluation score is a real number, always comparable with other score of dictionary like this. It allows to states which CBE_slice has more
                                             influence over each locus (once such a dictionary is available for each CBE_slice derived from whole_CBE). A non-negative score asserts that CBE_slice
                                             has influence over related locus, the higher the more.
                        Comment about evaluation: loci belonging to CBE_slice and loci whose bins are 0 (see (*) comment in code) are not evaluated
    
    TYPICAL USAGE IN THIS CONTEXT: CBE_slice has to be the slice of whole_CBE with the highest peak. After first use, if you want to loop, you have to remove CBE_slice content from whole_CBE
                                   then pass to the second CBE_slice in terms of peak height 
    '''
    
    # CBE strand
    CBE_slice_strand = str(CBE_slice.strand) # Cast should be redundant
    
    # Fix boundaries for hist_gauss_normalized_to_peak
    hist_gauss_index_of_max = hist_gauss_normalized_to_peak.index(max(hist_gauss_normalized_to_peak))
    gauss_n_step_right = hist_gauss_normalized_to_peak.index(hist_gauss_normalized_to_peak[-1]) - hist_gauss_index_of_max
    gauss_n_step_left = hist_gauss_index_of_max
    
    # Build histogram for whole_CBE
    whole_CBE_bin_areas, whole_CBE_list_of_loci, whole_CBE_max_read_count, whole_CBE_index_of_max = CBE__histogram_generator(whole_CBE)
    whole_CBE_bin_areas_normalized = normalize_histogram_to_the_peak(whole_CBE_bin_areas, whole_CBE_index_of_max)
    del whole_CBE_max_read_count
    
    # n of allowed step left and right from whole_CBE's peak
    index_last_bin = len(whole_CBE_bin_areas) - 1
    n_step_right = index_last_bin - whole_CBE_index_of_max
    if (n_step_right > gauss_n_step_right):
        n_step_right = gauss_n_step_right
    n_step_left = whole_CBE_index_of_max
    if (n_step_left > gauss_n_step_left):
        n_step_left = gauss_n_step_left
        
    # starting and ending indexes for hist_gauss
    starting_index = hist_gauss_index_of_max - n_step_left
    ending_index = hist_gauss_index_of_max + n_step_right
    
    # list of allowed indexes hist_gauss
    allowed_indexes_gauss = range(starting_index, ending_index+1)
    
    # starting and ending indexes for current_ensemble_bin_areas_normalized
    starting_index = whole_CBE_index_of_max - n_step_left
    ending_index = whole_CBE_index_of_max + n_step_right
    
    # list of allowed indexes for current_ensemble_bin_areas_normalized
    allowed_indexes_CBE = range(starting_index, ending_index+1)
    
    # list of indexes tuples [(allowed_indexes_gauss1, allowed_indexes_CBE1), (allowed_indexes_gauss2, allowed_indexes_CBE2), ... ]
    # indexes of peak's locus and of loci beside peak will be excluded strand-specifically
    indexes_tuples = []
    for i in range(0, n_step_left+n_step_right+1):
        if ((CBE_slice_strand == '+')or(CBE_slice_strand == '1')):
            if ((allowed_indexes_CBE[i] != whole_CBE_index_of_max) and (allowed_indexes_CBE[i] != whole_CBE_index_of_max + 1)):
                indexes_tuples.append((allowed_indexes_gauss[i],allowed_indexes_CBE[i]))
        else:
            if ((allowed_indexes_CBE[i] != whole_CBE_index_of_max) and (allowed_indexes_CBE[i] != whole_CBE_index_of_max - 1)):
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



        
def global_score_dictionary (CBE_list_of_slices, whole_CBE, hist_gauss_normalized_to_peak, strand_specific_choice):
    '''
    *** it returns a global score dictionary about loci in whole_CBE and each CBE slice in CBE_list_of_slices ***
     [Designed as a smart 'looping box' for evaluate_surroundings, helping with correct usage and showing merged 
                                       results in a unique dictionary]
                                       
    INPUT: see evaluate_surroundings function above.
           IMPORTANT NOTE about hist_gauss_normalized_to_peak: IT'S UP TO YOU PASSING THE CORRECT HISTOGRAM, ACCORDING TO CBE's STRAND
    
    OUTPUT: global_score_dic - a dictionary of kind: {locus:[(CBE_slice, score), (...), ...]}, obtained merging 'score_dic' dictionary
                               returned by evaluate_surroundings function calls (if the key is the same, items (tuple) are joined in a list 
                                       
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
        
        current_dic = evaluate_surroundings (CBE_slice, current_ensemble, hist_gauss_normalized_to_peak)
        
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
           list_of_bases_to_assign - a list of kind: [(covered_base object to assign, CBE_slice claiming it), (...), ...]
    
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
    
    
    
        