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
    
    #Thanks to class design, classic methods operates without distinguish between labels
    
    #IS object instance
    IS_object = Classes_for_Integration_Analysis.IS(Covered_bases_ensamble_object)
    
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
    
    return IS_object

##############################################################################################


