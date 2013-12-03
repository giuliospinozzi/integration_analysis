###Header################################################
header = """

+------------------------------------------------------+
 Module: Function_for_Gaussian_IS_identification
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  November 27th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - This module contains functions used in Gaussian
    IS identification framework [...]
  
 Note: [...]

-------------------------------------------------------- 
""" 
########################################################

###Requested Package(s) Import#
import numpy
###############################

###Import Module(s)#
import Classes_for_Integration_Analysis
####################





def gaussian_histogram_generator (interaction_limit, alpha):
    '''
    *** This function generates a Gaussian-shaped histogram ***
    
    INPUT: interaction_limit - Integer number (read description and note CAREFULLY)
           alpha - Float number  (read description and note CAREFULLY)
           
    OUTPUT: bin_boundaries - list of tuples of two elements, setting bins boundaries
            bin_areas - list of float, setting the area of each bin
            diagnostic - float number, contains the fraction of distribution lost (area of cutted tails)
    
            
    DESCRIPTION AND NOTES:
    
    # INTERACTION_LIMIT
    interaction_limit states, de facto, the number of bin of the histogam you get (n_bins = (2*interaction_limit)+1)
    ! Note that for theoretical consistency, it should be equal to half of bushamn_bp_rule
    
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
    
    
    # SOME GOOD GUESSES
    
    # stringent: low range of interaction, high precision; 7bins and 14sigmas of span (7 each side)
    alpha = 1
    interaction_limit = 3
    
    # an alternative, less stringent (in term of interaction limit and tail loss - nevertheless acceptable) but better shaped;
    # 9 bins and 9 sigmas of span (4.5 each side)
    alpha = 0.5
    interaction_limit = 4
    
    # maybe a compromise (we loose 0.002% of distribution - 2 orders of magnitude smaller than boundaries bins)
    alpha = 0.6 #(good gaussian shape)
    interaction_limit = 3 #(as bushman said)
    
    '''
    
    ### Function for normal  CDF ######################################################################
    def norm_CDF_from_x_to_y (x, y, alpha):
        # Given alpha (the number of sigmas equals to half-bp)
        # and given x and y (coordinates in bp unit), this funcion
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
    ###################################################################################################
    
    
    # Bin Boundaries Half
    first_bin_half = (0.0, 0.5)
    bin_boundaries_half = [first_bin_half]
    
    for i in range(0, interaction_limit):
        bin_boundaries_half.append((bin_boundaries_half[-1][1], bin_boundaries_half[-1][1]+1))
    
    # Filling Bins
    bin_areas_half = []
    for x,y in bin_boundaries_half:
        bin_areas_half.append(norm_CDF_from_x_to_y (x, y, alpha))
    
    # Simmetrize Bins
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
    *** From a covered bases ensamble objects, it creates an histogram ***
    
    INPUT: Covered_bases_ensamble object
    
    OUTPUT: bin_areas - List of Float, representing the histogram of the
                        Covered_bases_ensamble given in input (ReadCount)
            list_of_loci - list of int, reporting the loci for bin_areas
            max_read_count - Float, the highest ReadCount in bin_areas
            index_of_max - Integer, index of max_read_count in bin_areas 
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
    *** Normalize an histogram to make peak area 1 ***
    [...]
    '''
    max_height = bin_areas[index_of_max]
    bin_areas_normalized = []
    for one_bin in bin_areas:
        bin_areas_normalized.append(one_bin/max_height)
        
    return bin_areas_normalized


