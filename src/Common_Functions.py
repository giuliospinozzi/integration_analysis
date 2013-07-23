###Header###############################################
header = """

+------------------------------------------------------+
 Module: Common_functions
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 22th, 2013
 Contact: brasca.stefano@hsr.it, spinozzi.giulio@hsr.it
 Version: 0.1
+------------------------------------------------------+

 Description:
  - [...] Set of functions of common usage
  
 Note:
  - [...]

-------------------------------------------------------- 
""" 
########################################################


###Given a key for read dictionary, it returns related lam data###################
def get_lam (reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay):
    lam_data_dictionay_Key = reads_data_dictionary[reads_data_dictionary_Key][-1]
    lam_data = lam_data_dictionay[lam_data_dictionay_Key]
    return lam_data #tupla
##################################################################################
