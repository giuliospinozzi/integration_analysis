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

###############################

###Import Module(s)###########################
import Function_for_Gaussian_IS_identification
##############################################





def get_alpha (N_sigma, interaction_limit):
    
    span = 2*interaction_limit+1   # def    
    alpha = N_sigma / float(span)   # resulted from exact calculations
    
    return alpha
    
    



def get_sigma_path (MAX_interaction_limit, MIN_interaction_limit, N_sigma):
    
    sigma_path_point_list = []
    interaction_limit_range = range(MIN_interaction_limit, MAX_interaction_limit+1)
    
    for interaction_limit in interaction_limit_range:
        alpha = get_alpha(N_sigma, interaction_limit)
        sigma_path_point = (interaction_limit, alpha)
        sigma_path_point_list.append(sigma_path_point)
        
    return sigma_path_point_list   # list of (interaction_limit, alpha) tuples
    


    
def get_ranking_histograms (MAX_interaction_limit, sigma_paths, MIN_interaction_limit=1, adaptive_SW=True):
    ### sigma_paths is a list of kind [2,3,...] (sigma path to explore)

    ###############################################################################################
    def limit_histogram_generator (interaction_limit):
        ### per alpha che 'tende a zero'
        
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
        bin_areas[interaction_limit] = bin_areas_half[0] * 2.0
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
        pass ### ERRORE!!! Uscire ###
        
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
                                  
                
                
        
        