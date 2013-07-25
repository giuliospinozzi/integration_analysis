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


###Import Module(s)####
import Common_Functions #requested for collapse method for Covered_base class
#######################


###Class of covered bases#############################################################################################      
class Covered_base:
    '''
    Class of covered bases.
    [...]
    '''
    
    #Constructor######################################################################################################
    def __init__(self, reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay):
        '''
        [...]
        '''
        self.list_of_reads_key = [reads_data_dictionary_Key]
        self.chromosome = reads_data_dictionary[reads_data_dictionary_Key][1]
        self.strand = reads_data_dictionary[reads_data_dictionary_Key][2]
        self.locus = reads_data_dictionary[reads_data_dictionary_Key][3]
        self.reads_count = 1
        self.selective_reads_count = [] # a list of labels, one for each read in this object (the first created below)
                                        # 'collapse' method changes it into a dictionary (see 'collapse')
        #Create first elemnt of selective_reads_count list
        lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay) #retrieve lam data related to the read, by means of "get_lam" function in "Common_Functions" module : stored in lam_data
        column_label = "{0}_{1}_{2}".format(lam_data[4],lam_data[3],lam_data[5]) #from lam_data create a label of kind 'sample'_'tissue'_'treatment'
        self.selective_reads_count.append(column_label) #append column_label in selective_reads_count
    ###################################################################################################################


    #Methods###########################################################################################################

    #Add a read in Covered_base object with controls about chromosome, strand and locus matching###
    #If the read passes controls it will be added, Covered_base attributes will be updated and you get 1 (if clause)###
    #Else nothing happens and you get -1 (else clause)###       
    def add (self, reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay):
        '''
        [...]
        '''
        if ((reads_data_dictionary[reads_data_dictionary_Key][1] == self.chromosome) and (reads_data_dictionary[reads_data_dictionary_Key][2] == self.strand) and (reads_data_dictionary[reads_data_dictionary_Key][3] == self.locus)):
            self.list_of_reads_key.append(reads_data_dictionary_Key)
            self.reads_count = self.reads_count + 1
            lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay)
            column_label = "{0}_{1}_{2}".format(lam_data[4],lam_data[3],lam_data[5])
            self.selective_reads_count.append(column_label)
            return 1
        else:
            return -1
        
    #Collapse method change "selective_reads_count" attribute (originally a list, see above) into a dictionary of kind { 'label1' : #n_of_label1_in_selective_reads_count_list, 'label2' : #n_of_label2_in_selective_reads_count_list, ... } ###
    #You should use this method AFTER having added EVERY READ you need, mainly because of type change for selective_reads_count attribute ###
    #Return nothing###     
    def collapse (self):
        '''
        [...]
        '''
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
        
    #Distance method for Covered_base returns distance from another Covered_base, given in input.
    #If the distance doesn't make sense (e.g. distance between CBs in different chromosome, this method returns 'undef' instead of a number
    def distance (self, another_Covered_base):
        '''
        [...]
        '''
        dist = "undef"
        if ((self.chromosome == another_Covered_base.chromosome) and (self.strand == another_Covered_base.strand) and (self.locus != another_Covered_base.locus) and (self != another_Covered_base)):
            dist = abs(self.locus - another_Covered_base.locus)
        return dist
            
    ###################################################################################################################

#######################################################################################################################




###Class of read sequences ensemble####################################################################################      
class Covered_bases_ensamble:
    '''
    Class of covered bases ensembles, grouped by mutual distance (Bushman bp rule)
    [...]
    '''
 
    #Constructor####################################################################################################### 
    def __init__(self, Covered_base_object):
        '''
        [...]
        '''
        self.Covered_bases_list = [Covered_base_object]
        self.chromosome = Covered_base_object.chromosome
        self.strand = Covered_base_object.strand
        self.starting_base_locus = Covered_base_object.locus
        self.ending_base_locus = Covered_base_object.locus
        self.spanned_bases = 1
        self.n_covered_bases = 1
        self.n_total_reads = Covered_base_object.reads_count
        self.covered_base_of_max = Covered_base_object
    ####################################################################################################################
        
    #Methods############################################################################################################
    def push_in (self, Covered_base_object):
        check = -1
        if ((Covered_base_object not in self.Covered_bases_list) and (self.chromosome == Covered_base_object.chromosome) and (self.strand == Covered_base_object.strand)):
            self.Covered_bases_list.append(Covered_base_object)
            self.starting_base_locus = min(self.starting_base_locus, Covered_base_object.locus)
            self.ending_base_locus = max(self.ending_base_locus, Covered_base_object.locus)
            self.spanned_bases = self.ending_base_locus - self.starting_base_locus + 1
            self.n_covered_bases = self.n_covered_bases + 1
            self.n_total_reads = self.n_total_reads + Covered_base_object.reads_count
            if (Covered_base_object.reads_count > self.covered_base_of_max.reads_count):
                self.covered_base_of_max = Covered_base_object
            check = 1
        return check
    ####################################################################################################################    

    
########################################################################################################################      
        
        
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