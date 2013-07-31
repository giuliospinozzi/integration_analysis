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



def classic (Covered_bases_ensamble_object):
    '''
    Classic way to retrieve IS from Covered Bases Ensembles:
    covered base with the higher number of reads sets the "integration locus"
    while the related reads count is fixed as the overall reads count of the whole
    ensemble.
    '''
    
    #IS object instance
    IS_object = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object)
    
    #For label "all"
    if (Covered_bases_ensamble_object.label == "all"):
        pass
    
    #For merged labels
    elif (Covered_bases_ensamble_object.label[0] == "_"):
        pass
    
    #For standard labels (matrix column labels)
    else:
        pass
    
    return IS_object