###Header################################################
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
    def __init__(self, reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific = True):
        '''
        [...]
        '''
        self.list_of_reads_key = [reads_data_dictionary_Key]
        self.chromosome = reads_data_dictionary[reads_data_dictionary_Key][1]
        self.strand = None
        self.strand_aspecific = None
        if (strand_specific == True):
            self.strand = reads_data_dictionary[reads_data_dictionary_Key][2]
        elif (strand_specific == False):
            self.strand_aspecific = reads_data_dictionary[reads_data_dictionary_Key][2]
        self.locus = reads_data_dictionary[reads_data_dictionary_Key][3]
        self.reads_count = 1
        self.selective_reads_count = [] # a list of labels, one for each read in this object (the first created below)
                                        # 'collapse' method changes it into a dictionary (see 'collapse')
        #Create first element of selective_reads_count list
        lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay) #retrieve lam data related to the read, by means of "get_lam" function in "Common_Functions" module : stored in lam_data
        column_label = ""
        if ("group_name" in parameters_list):
            column_label = column_label + "_" + lam_data[6]
        if ("n_LAM" in parameters_list):
            column_label = column_label + "_" + lam_data[0]
        if ("pool" in parameters_list):
            column_label = column_label + "_" + lam_data[2]
        if ("tag" in parameters_list):
            column_label = column_label + "_" + lam_data[1]
        if ("enzyme" in parameters_list):
            column_label = column_label + "_" + lam_data[7]    
        if ("sample" in parameters_list):
            column_label = column_label + "_" + lam_data[4]    
        if ("tissue" in parameters_list):
            column_label = column_label + "_" + lam_data[3]
        if ("treatment" in parameters_list):
            column_label = column_label + "_" + lam_data[5]
        column_label = column_label[1:]
        self.selective_reads_count.append(column_label) #append column_label in selective_reads_count
    ###################################################################################################################


    #Methods###########################################################################################################

    #Add a read in Covered_base object with controls about chromosome, strand and locus matching###
    #If the read passes controls it will be added, Covered_base attributes will be updated and you get 1 (if clause)###
    #Else nothing happens and you get -1 (else clause)###
    #"strand specific" option allows to choose if consider strand in controls or not ###       
    def add (self, reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific = True):
        '''
        [...]
        '''
        if (strand_specific == True):
            if ((reads_data_dictionary[reads_data_dictionary_Key][1] == self.chromosome) and (reads_data_dictionary[reads_data_dictionary_Key][2] == self.strand) and (reads_data_dictionary[reads_data_dictionary_Key][3] == self.locus)):
                self.list_of_reads_key.append(reads_data_dictionary_Key)
                self.reads_count = self.reads_count + 1
                lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay)
                column_label = ""
                if ("group_name" in parameters_list):
                    column_label = column_label + "_" + lam_data[6]
                if ("n_LAM" in parameters_list):
                    column_label = column_label + "_" + lam_data[0]
                if ("pool" in parameters_list):
                    column_label = column_label + "_" + lam_data[2]
                if ("tag" in parameters_list):
                    column_label = column_label + "_" + lam_data[1]
                if ("enzyme" in parameters_list):
                    column_label = column_label + "_" + lam_data[7]    
                if ("sample" in parameters_list):
                    column_label = column_label + "_" + lam_data[4]    
                if ("tissue" in parameters_list):
                    column_label = column_label + "_" + lam_data[3]
                if ("treatment" in parameters_list):
                    column_label = column_label + "_" + lam_data[5]
                column_label = column_label[1:]
                self.selective_reads_count.append(column_label)
                return 1
            else:
                return -1
        elif (strand_specific == False):
            if ((reads_data_dictionary[reads_data_dictionary_Key][1] == self.chromosome) and (reads_data_dictionary[reads_data_dictionary_Key][3] == self.locus)):
                self.list_of_reads_key.append(reads_data_dictionary_Key)
                self.reads_count = self.reads_count + 1
                lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay)
                column_label = ""
                if ("group_name" in parameters_list):
                    column_label = column_label + "_" + lam_data[6]
                if ("n_LAM" in parameters_list):
                    column_label = column_label + "_" + lam_data[0]
                if ("pool" in parameters_list):
                    column_label = column_label + "_" + lam_data[2]
                if ("tag" in parameters_list):
                    column_label = column_label + "_" + lam_data[1]
                if ("enzyme" in parameters_list):
                    column_label = column_label + "_" + lam_data[7]    
                if ("sample" in parameters_list):
                    column_label = column_label + "_" + lam_data[4]    
                if ("tissue" in parameters_list):
                    column_label = column_label + "_" + lam_data[3]
                if ("treatment" in parameters_list):
                    column_label = column_label + "_" + lam_data[5]
                column_label = column_label[1:]
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
    #If the distance doesn't make sense at all (e.g. distance between CBs in different chromosome) this method returns 'undef' instead of a number
    def distance (self, another_Covered_base, label_selection = "all"):
        '''
        [...]
        '''
        dist = "undef"
        
        if ((label_selection == "all") and (self.chromosome == another_Covered_base.chromosome) and (self.strand == another_Covered_base.strand) and (self.locus != another_Covered_base.locus) and (self != another_Covered_base)):
            dist = abs(self.locus - another_Covered_base.locus)
        #label-selective: dist remains "undef" if, for a given label, both covered bases haven't a non-zero reads count for that label
        elif ((label_selection in another_Covered_base.selective_reads_count.keys()) and (self.chromosome == another_Covered_base.chromosome) and (self.strand == another_Covered_base.strand) and (self.locus != another_Covered_base.locus) and (self != another_Covered_base)):
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
    def __init__(self, Covered_base_object, label_selection = "all", strand_specific = True):
        '''
        [...]
        '''
        self.label = label_selection #you can create Covered_bases_ensamble looking only at one label
        self.Covered_bases_list = [Covered_base_object]
        self.chromosome = Covered_base_object.chromosome
        self.strand = Covered_base_object.strand
        self.starting_base_locus = Covered_base_object.locus
        self.ending_base_locus = Covered_base_object.locus
        self.spanned_bases = 1
        self.n_covered_bases = 1
        self.strand_aspecific = None
        
        if (strand_specific == False):
            self.strand_aspecific = Covered_base_object.strand_aspecific
        
        if (label_selection == "all"):
            self.n_total_reads = Covered_base_object.reads_count
        else: #label-selective
            self.n_total_reads = Covered_base_object.selective_reads_count[label_selection]
                    
        self.covered_base_of_max = Covered_base_object
    ####################################################################################################################
        
    #Methods############################################################################################################
    def push_in (self, Covered_base_object, label_selection = "all"):
        check = -1

        if ((label_selection == "all") and (Covered_base_object not in self.Covered_bases_list) and (self.chromosome == Covered_base_object.chromosome) and (self.strand == Covered_base_object.strand)):
            self.Covered_bases_list.append(Covered_base_object)
            self.starting_base_locus = min(self.starting_base_locus, Covered_base_object.locus)
            self.ending_base_locus = max(self.ending_base_locus, Covered_base_object.locus)
            self.spanned_bases = self.ending_base_locus - self.starting_base_locus + 1
            self.n_covered_bases = self.n_covered_bases + 1
            self.n_total_reads = self.n_total_reads + Covered_base_object.reads_count
            if (Covered_base_object.reads_count > self.covered_base_of_max.reads_count):
                self.covered_base_of_max = Covered_base_object
            check = 1
        #label-selective
        elif ((Covered_base_object not in self.Covered_bases_list) and (self.chromosome == Covered_base_object.chromosome) and (self.strand == Covered_base_object.strand)):
            self.Covered_bases_list.append(Covered_base_object)
            self.starting_base_locus = min(self.starting_base_locus, Covered_base_object.locus)
            self.ending_base_locus = max(self.ending_base_locus, Covered_base_object.locus)
            self.spanned_bases = self.ending_base_locus - self.starting_base_locus + 1
            self.n_covered_bases = self.n_covered_bases + 1
            self.n_total_reads = self.n_total_reads + Covered_base_object.selective_reads_count[label_selection]
            if (Covered_base_object.selective_reads_count[label_selection] > self.covered_base_of_max.selective_reads_count[label_selection]):
                self.covered_base_of_max = Covered_base_object
            check = 1
            
        return check
    ####################################################################################################################    
    
########################################################################################################################      
        


        
###Class of Integration Sites###########################################################################################         
class IS:
    '''
    Class of Integration Sites
    [...]
    '''
 
 
    def __init__(self, Covered_bases_ensamble_object, strand_specific = True):
        '''
        [...]
        '''
        self.related_ensemble = Covered_bases_ensamble_object
        self.chromosome = Covered_bases_ensamble_object.chromosome
        self.strand = Covered_bases_ensamble_object.strand
        self.starting_base_locus = None #to be evaluated by Integration_Sites_retrieving_methods
        self.ending_base_locus = None #to be evaluated by Integration_Sites_retrieving_methods
        self.integration_locus = None #to be evaluated by Integration_Sites_retrieving_methods
        self.spanned_bases = None #to be evaluated by Integration_Sites_retrieving_methods
        self.n_covered_bases = None #to be evaluated by Integration_Sites_retrieving_methods
        self.reads_count = None #to be evaluated by Integration_Sites_retrieving_methods
        self.selective_reads_count = None
        self.peak_height = None  #to be evaluated by Integration_Sites_retrieving_methods
        self.reads_key_list = None #to be evaluated by Integration_Sites_retrieving_methods
        self.strand_aspecific = None 
        if (strand_specific == False):
            self.strand_aspecific = Covered_bases_ensamble_object.strand_aspecific
        
########################################################################################################################         
