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

###Import Module(s)###
import Common_Functions
######################

###Class of covered bases######################       
class Covered_base:
    '''
    Class of covered bases.
    [...]
    '''
 
 
    def __init__(self, reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay):
        '''
        [...]
        '''
        self.list_of_reads = [reads_data_dictionary_Key]
        self.chromosome = reads_data_dictionary[reads_data_dictionary_Key][1]
        self.strand = reads_data_dictionary[reads_data_dictionary_Key][2]
        self.locus = reads_data_dictionary[reads_data_dictionary_Key][3]
        self.reads_count = 1
        self.selective_reads_count = []
        lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay)
        column_label = "{0}_{1}_{2}".format(lam_data[4],lam_data[3],lam_data[5])
        self.selective_reads_count.append(column_label)

         
    def add (self, reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay):
        '''
        [...]
        '''
        if ((reads_data_dictionary[reads_data_dictionary_Key][1] == self.chromosome) and (reads_data_dictionary[reads_data_dictionary_Key][2] == self.strand) and (reads_data_dictionary[reads_data_dictionary_Key][3] == self.locus)):
            self.list_of_reads.append(reads_data_dictionary_Key)
            self.reads_count = self.reads_count + 1
            lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay)
            column_label = "{0}_{1}_{2}".format(lam_data[4],lam_data[3],lam_data[5])
            self.selective_reads_count.append(column_label)
            return 1
        else:
            return -1
        
    def collapse (self):
        i=0
        self.selective_reads_count.sort()
        collapsed_selective_reads_count = {self.selective_reads_count[0]:1}
        for label in self.selective_reads_count[1:]:
            i+=1
            if (label in collapsed_selective_reads_count.keys()):
                collapsed_selective_reads_count[label] = collapsed_selective_reads_count[label] + 1
            else:
                collapsed_selective_reads_count.update({self.selective_reads_count[i]:1})
        self.selective_reads_count = collapsed_selective_reads_count
                
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