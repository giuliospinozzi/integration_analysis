###Header###############################################
header = """

+------------------------------------------------------+
 Module: Classes_for_Integration_Analysis
 Author: Stefano Brasca, Giulio Spinozzi
 Date:  July 12th, 2013
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


###Class of covered bases######################       
class Covered_base:
    '''
    Class of covered bases.
    [...]
    '''
 
 
    def __init__(self, reads_data_dictionary_Key, reads_data_dictionary):
        '''
        [...]
        '''
        self.list_of_reads = [reads_data_dictionary_Key]
        self.chromosome = reads_data_dictionary[reads_data_dictionary_Key][1]
        self.strand = reads_data_dictionary[reads_data_dictionary_Key][2]
        self.locus = reads_data_dictionary[reads_data_dictionary_Key][3]
        self.reads_count = 1
         
    def add (self, reads_data_dictionary_Key, reads_data_dictionary):
        '''
        [...]
        '''
        if ((reads_data_dictionary[reads_data_dictionary_Key][1] == self.chromosome) and (reads_data_dictionary[reads_data_dictionary_Key][2] == self.strand) and (reads_data_dictionary[reads_data_dictionary_Key][3] == self.locus)):
            self.list_of_reads.append(reads_data_dictionary_Key)
            self.reads_count = self.reads_count + 1
            return 1
        else:
            return -1
         
################################################



#===============================================================================
# ###Class of read sequences ensemble######################       
# class Covered_bases_ensamble:
#     '''
#     Class of covered bases ensembles, grouped by mutual distance (Bushman bp rule)
#     [...]
#     '''
# 
# 
#     def __init__(selfparams):
#         '''
#         [...]
#         '''
# ########################################################       
#===============================================================================
        
        
#===============================================================================
# ###Class of Integration Sites###########################        
# class IS:
#     '''
#     Class of Integration Sites
#     [...]
#     '''
# 
# 
#     def __init__(selfparams):
#         '''
#         [...]
#         '''
# ########################################################        
#===============================================================================