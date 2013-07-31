###Header###############################################
header = """

+------------------------------------------------------+
 Module: Integration_Sites_retrieving_methods
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 31th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - [...]
  
 Note:
  - [...]

-------------------------------------------------------- 
""" 
########################################################


###Import Module(s)####################
import Classes_for_Integration_Analysis
#######################################


###Classic method#############################################################################

def classic (Covered_bases_ensamble_object):
    '''
    Classic way to retrieve IS from Covered Bases Ensembles:
    covered base with the higher number of reads sets the "integration locus"
    while the related reads count is fixed as the overall reads count of the whole
    ensemble.
    '''
    
    #===========================================================================
    # #Template for further methods that will require to distinguish between labels
    # 
    # #For label "all"
    # if (Covered_bases_ensamble_object.label == "all"):
    #     pass
    # 
    # #For merged labels
    # elif (Covered_bases_ensamble_object.label[0] == "_"):
    #     pass
    # 
    # #For standard labels (matrix column labels)
    # else:
    #     pass
    #===========================================================================
    
    #Thanks to class design, classic methods operates without distinguish between labels
    
    #IS object instance
    IS_object = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object)
    
    #Set Locus
    IS_object.integration_locus = Covered_bases_ensamble_object.covered_base_of_max.locus
    
    #Overall Reads Count
    IS_object.reads_count = Covered_bases_ensamble_object.n_total_reads
    
    return IS_object

##############################################################################################


