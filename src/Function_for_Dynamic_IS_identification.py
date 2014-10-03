###Header################################################
header = """

+------------------------------------------------------+
 Module: Function_for_Dynamic_IS_identification
 Author: Stefano Brasca
 Date:  September 23th, 2014
 Contact: brasca.stefano@hsr.it
 Version: inDevelopment
+------------------------------------------------------+

 Description:
  - This module contains functions used in Dynamic IS
    identification framework
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import sys
###############################

###Import Module(s)###########################
import Function_for_Gaussian_IS_identification
##############################################





def get_alpha (N_sigma, interaction_limit):
    '''
    given a window of 2*interaction_limit+1, this function yields the proper
    alpha parameter in order to create a gaussian histogram with N_sigma falling
    inside the window (es N_sigma=2 -> 95% of distr inside the window)
    [for Function_for_Gaussian_IS_identification.gaussian_histogram_generator callings]
    '''
    
    span = 2*interaction_limit+1   # def    
    alpha = N_sigma / float(span)   # resulted from exact calculations
    
    return alpha
    
    



def get_sigma_path (MAX_interaction_limit, MIN_interaction_limit, N_sigma):
    '''
    given the interaction limit range to explore, this function yields the list
    of (interaction_limit, alpha) points that respect N_sigma constraint
    [for Function_for_Gaussian_IS_identification.gaussian_histogram_generator callings]
    
    NOTE: MAX_interaction_limit must be >= MIN_interaction_limit; if =, only one 
    point will be returned
    '''
    
    sigma_path_point_list = []
    interaction_limit_range = range(MIN_interaction_limit, MAX_interaction_limit+1)
    
    for interaction_limit in interaction_limit_range:
        alpha = get_alpha(N_sigma, interaction_limit)
        sigma_path_point = (interaction_limit, alpha)
        sigma_path_point_list.append(sigma_path_point)
        
    return sigma_path_point_list   # list of (interaction_limit, alpha) tuples
    


    
def get_ranking_histograms (MAX_interaction_limit, sigma_paths, MIN_interaction_limit=1, adaptive_SW=True):
    '''
    INPUT: MAX_interaction_limit, MIN_interaction_limit - int
           sigma_paths - LIST of kind [2,3,...] (sigma path to explore)
           adaptive_SW - boolean
    OUTPUT: ranking_histogram_dict_list - list of dictionary entries, descripted below
    NOTE: - MAX_interaction_limit must be >= MIN_interaction_limit;
          - sigma_paths must be always a list, even if only one sigma value has to be exploited
    
    DESCRIPTION:
    this function implement 'get_sigma_path' above and gaussian_histogram_generator and 
    normalize_histogram_to_the_peak from Function_for_Gaussian_IS_identification module.
    it yields a list of dictionary of kind:
        
        dict_entry = {'type': "gaussian selection",
                      'N_sigma': N_sigma,
                      'interaction_limit': interaction_limit,
                      'alpha': alpha,
                      'hist_gauss': bin_areas,
                      'bin_boundaries': bin_boundaries,
                      'excluded_perc': diagnostic*100,
                      'hist_normalized_to_peak': hist_gauss_normalized_to_peak}
                      
    if adaptive_SW is True (this name was chosen thinkig at the expected behaviuor of the algorithm),
    further dictionaries will be added to the returned list, derived from limit cases (alpha ->0 or 
    alpha >> 1)
    
        dict_entry = {'type': "adaptive SW",
                              'interaction_limit': interaction_limit,
                              'bin_boundaries': bin_boundaries,
                              'hist_normalized_to_peak': bin_areas}
                              
    NOTE: limit cases of alpha -> 0 are one for each different interaction_limit while for alpha >> 1
    the case is only one (and can be assimilated to the (alpha -> 0, interaction_limit=3) case)
    '''

    ###############################################################################################
    def limit_histogram_generator (interaction_limit):
        # Bin Boundaries Half
        first_bin_half = (0.0, 0.5)
        bin_boundaries_half = [first_bin_half]        
        for i in range(0, interaction_limit):
            bin_boundaries_half.append((bin_boundaries_half[-1][1], bin_boundaries_half[-1][1]+1))        
        # Filling Bins
        bin_areas_half = []
        for x,y in bin_boundaries_half:
            bin_areas_half.append(1.0)
        # Symmetrize Bins
        n_bins = int((2*interaction_limit)+1)
        bin_boundaries = [None]*n_bins
        bin_areas = [None]*n_bins
        bin_boundaries[interaction_limit] = (-0.5, 0.5)
        bin_areas[interaction_limit] = bin_areas_half[0]
        for i in range(1, interaction_limit +1):
            negative_left = -1.0*bin_boundaries_half[interaction_limit+1-i][1]
            negative_right = -1.0*bin_boundaries_half[interaction_limit+1-i][0]
            bin_boundaries[i-1] = (negative_left, negative_right)
            bin_boundaries[i+interaction_limit] = bin_boundaries_half[i]
            bin_areas[i-1] = bin_areas_half[interaction_limit+1-i]
            bin_areas[i+interaction_limit] = bin_areas_half[i]            
        # Remember that bin_boundaries[interaction_limit] is the peak        
        return bin_boundaries, bin_areas
    ###############################################################################################
    
    
    if ((sigma_paths is None) and (adaptive_SW is False)):
        print "\n\n\t get_ranking_histograms function can't run, due to sigma_paths and adaptive_SW variables both set to 'None' and 'False'."
        sys.exit("\n\n\t[ERROR]\tQuit.\n\n")
        
    else:
        
        ranking_histogram_dict_list = []
        
        if (adaptive_SW is True):
            for interaction_limit in range(MIN_interaction_limit, MAX_interaction_limit+1):
                bin_boundaries, bin_areas = limit_histogram_generator (interaction_limit)
                dict_entry = {'type': "adaptive SW",
                              'interaction_limit': interaction_limit,
                              'bin_boundaries': bin_boundaries,
                              'hist_normalized_to_peak': bin_areas}
                ranking_histogram_dict_list.append(dict_entry)
        
        if (sigma_paths is not None):            
            for N_sigma in sigma_paths:
                sigma_path_point_list = get_sigma_path(MAX_interaction_limit, MIN_interaction_limit, N_sigma)
                for interaction_limit, alpha in sigma_path_point_list:
                    bin_boundaries, bin_areas, diagnostic = Function_for_Gaussian_IS_identification.gaussian_histogram_generator(interaction_limit, alpha)
                    hist_gauss_normalized_to_peak = Function_for_Gaussian_IS_identification.normalize_histogram_to_the_peak(bin_areas, interaction_limit)
                    dict_entry = {'type': "gaussian selection",
                                  'N_sigma': N_sigma,
                                  'interaction_limit': interaction_limit,
                                  'alpha': alpha,
                                  'hist_gauss': bin_areas,
                                  'bin_boundaries': bin_boundaries,
                                  'excluded_perc': diagnostic*100,
                                  'hist_normalized_to_peak': hist_gauss_normalized_to_peak}
                    ranking_histogram_dict_list.append(dict_entry)
                    
        return ranking_histogram_dict_list
                                  
                
                
        
        