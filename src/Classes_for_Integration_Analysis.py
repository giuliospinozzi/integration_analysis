###Header################################################
header = """

+------------------------------------------------------+
 Module: Classes_for_Integration_Analysis
 Author: Stefano Brasca
 Date:  January 20th, 2014
 Contact: brasca.stefano@hsr.it
 Version: 1.0
+------------------------------------------------------+

 Description:
  - This module contains classes and methods used in the
    whole Integration Analysis project
  
 Note:
  - Description of methods are lacking and should be
    improved!

-------------------------------------------------------- 
""" 
########################################################


###Import Module(s)####
import Common_Functions #requested for read tracking
#######################


###Class of covered bases#############################################################################################      
class Covered_base:
    '''
    *** Class of covered bases ***
    
    CONCEPT: a 'covered base' (CB) is a genome location univocally determined by chromosome-locus(-strand) whose reads
             coverage (in terms of starting base position) consists at least of 1 read.
             To get an object of Covered_base class is therefore necessary and sufficient to provide only a single read:
             other potential reads will be added later, add(ing them).
             Such reads may come from different LAM, cell-types, samples, time-points...
             because of this, beside the overall reads-count, it's useful to provide a 'selective reads-count' that splits
             the total according to these categories.
             
    STRUCTURE:
    
        __init__ INPUT: - reads_data_dictionary: see output of import_reads_data_from_DB function in DB_connection module
                                                 (dictionary of reads identified by header (the key) )
                        - reads_data_dictionary_Key: the key of reads_data_dictionary related to the read you want to use
                                                     to create the covered base object
                        - lam_data_dictionay: see output of import_lam_data_from_DB function in DB_connection module
                                              (dictionary of lam data, identified and related to reads through lam_id (the key) )
                        - parameters_list: a list of string retrieved from user input (--columns arg), reflecting categories
                                           of interest for user among the ones in lam_data_dictionay
                        - strand_specific: boolean; if true (default), the covered base is univocally determined also by strand
                                           NOTE: the choice made will condition the following ones. For this purpose, you can find
                                           a variable called 'strand_specific_choice' in main, retrieved from user input, so the
                                           best usage is strand_specific = strand_specific_choice
                        
        ATTRIBUTES: 'list_of_reads_key' - list of headers (reads_data_dictionary keys) of the reads in the covered base
                    'chromosome' - string; the chromosome hosting the covered base
                    'strand' / 'strand_aspecific' - string; the strand hosting the covered base
                                                    NOTE: always both present and set 'None'. Then, if strand_specific = True,
                                                    'strand' is suddenly set as the actual strand; else, if strand_specific = False,
                                                    'aspecific_strand' is set as the strand of the read used to create the covered
                                                    base object.
                    'locus': long int; the locus hosting the covered base
                    'reads_count' - int; the overall number of reads in the covered base
                    'selective_reads_count' - list of labels, one for each read in the covered base, specifying the reads origin 
                                              according to parameters_list / lam_data_dictionay
                                              NOTE: 'collapse' method turn it into a dictionary of kind 
                                                    { 'label1' : #n_of_label1_in_selective_reads_count_list,
                                                    'label2' : #n_of_label2_in_selective_reads_count_list, ... }
                    'longest_seq_header' - for seqTrakcker; filled by collapse method
                    'longest_raw_seq' - for seqTrakcker; filled by collapse method
                    'longest_final_seq' - for seqTrakcker; filled by collapse method
        
        METHODS: add - [...]
                 collapse - [...]
                 distance - [...]
        
        NOTE for developers: - This rigid structure of labels in selective_reads_count is mandatory, due to algorithm conception, and coherent
                               in the whole Integration Analysis python project. DO NOT IMPROVISE CHANGES!
                             - 'Add' method refreshes attributes in real time. 
                             - When you finished adding reads, use collapse!
                                                     
             

    '''
    
    #Constructor######################################################################################################
    def __init__(self, reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific = True):
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
        self.longest_seq_header = None
        self.longest_raw_seq = None
        self.longest_final_seq = None
        #Create first element of selective_reads_count list
        lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay) #retrieve lam data related to the read, by means of "get_lam" function in "Common_Functions" module : stored in lam_data
        column_label = ""
        if ("group_name" in parameters_list):
            column_label = column_label + "_" + str(lam_data[6])
        if ("n_LAM" in parameters_list):
            column_label = column_label + "_" + str(lam_data[0])
        if ("pool" in parameters_list):
            column_label = column_label + "_" + str(lam_data[2])
        if ("tag" in parameters_list):
            column_label = column_label + "_" + str(lam_data[1])
        if ("enzyme" in parameters_list):
            column_label = column_label + "_" + str(lam_data[7])
        if ("sample" in parameters_list):
            column_label = column_label + "_" + str(lam_data[4])
        if ("tissue" in parameters_list):
            column_label = column_label + "_" + str(lam_data[3])
        if ("treatment" in parameters_list):
            column_label = column_label + "_" + str(lam_data[5])
        if ("vector" in parameters_list):
            column_label = column_label + "_" + str(lam_data[8])

        column_label = column_label[1:]
        self.selective_reads_count.append(column_label) #append column_label in selective_reads_count
    ###################################################################################################################


    #Methods###########################################################################################################

    #Add a read in Covered_base object with controls about chromosome, strand and locus matching###
    #If the read passes controls it will be added, Covered_base attributes will be updated and you get 1 (if clause)###
    #Else nothing happens and you get -1 (else clause)###
    #"strand specific" option allows to choose if consider strand in controls or not ###       
    def add (self, reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay, parameters_list, strand_specific = True):
        if (strand_specific == True):
            if ((reads_data_dictionary[reads_data_dictionary_Key][1] == self.chromosome) and (reads_data_dictionary[reads_data_dictionary_Key][2] == self.strand) and (reads_data_dictionary[reads_data_dictionary_Key][3] == self.locus)):
                self.list_of_reads_key.append(reads_data_dictionary_Key)
                self.reads_count = self.reads_count + 1
                lam_data = Common_Functions.get_lam(reads_data_dictionary_Key, reads_data_dictionary, lam_data_dictionay)
                column_label = ""
                if ("group_name" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[6])
                if ("n_LAM" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[0])
                if ("pool" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[2])
                if ("tag" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[1])
                if ("enzyme" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[7]) 
                if ("sample" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[4])
                if ("tissue" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[3])
                if ("treatment" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[5])
                if ("vector" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[8])
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
                    column_label = column_label + "_" + str(lam_data[6])
                if ("n_LAM" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[0])
                if ("pool" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[2])
                if ("tag" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[1])
                if ("enzyme" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[7])
                if ("sample" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[4])
                if ("tissue" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[3])
                if ("treatment" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[5])
                if ("vector" in parameters_list):
                    column_label = column_label + "_" + str(lam_data[8])
                column_label = column_label[1:]
                self.selective_reads_count.append(column_label)
                return 1
            else:
                return -1
            
        
    #Collapse method change "selective_reads_count" attribute (originally a list, see above) into a dictionary of kind { 'label1' : #n_of_label1_in_selective_reads_count_list, 'label2' : #n_of_label2_in_selective_reads_count_list, ... } ###
    #Further, it fills 'longest_seq_header', 'longest_raw_seq' and 'longest_final_seq' attributes if seqTracker is True.
    #You should use this method AFTER having added EVERY READ you need, mainly because of type change for selective_reads_count attribute ###
    #Return nothing###     
    def collapse (self, seqTracker, raw_read_dictionary, final_read_dictionary):
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
        
        # seqTracker
        if (seqTracker == True):
            selected_header, longest_raw_sequence, longest_final_sequence = Common_Functions.find_longest_read (self.list_of_reads_key, raw_read_dictionary, final_read_dictionary)
            self.longest_seq_header=selected_header
            self.longest_raw_seq=longest_raw_sequence
            self.longest_final_seq =longest_final_sequence
        
    #Distance method for Covered_base returns distance from another Covered_base, given in input.
    #If the distance doesn't make sense at all (e.g. distance between CBs in different chromosomes) this method returns 'undef' instead of a number
    def distance (self, another_Covered_base, label_selection = "all"):
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
    *** Class of covered bases ensembles ***
    
    CONCEPT: a 'covered_bases_ensemble' (CBE) is a set of covered_base objects, grouped according to somewhat criteria, in order
             to be easily processed together by Integration Sites Retrieving methods, on a second time; conceptually,
             this is a practical way to partition covered_bases from a dataset into groups: CBs in the same group are supposed
             to be mutually correlated hence they should be processed as a whole during IS retrieval.
             To get an object of Covered_bases_ensamble class is therefore necessary and sufficient to provide just a single CB:
             other potential CBs will be added later, 'push(ing them) in'.
             In order to belong to the same ensemble, covered bases must:
             - not to be the same CB (a kind of errors control, to avoid duplicates)
             - be placed on the same chromosome - mandatory
             - be placed on the same strand, if 'strand_specific option is turned 'True' (this is a real restriction only if CBs
               have been constructed 'strand specifically': otherwise 'strand' attribute is 'none' and strand control doesn't actually act)
             - belong to the same category (eg. sample, tissue, etc.), specified through 'label_selection' option
               (this feature is an old remains from formers versions of this program and should be considered as deprecated; letting  
                label_selection = "all" as defaults, this control doesn't actually act)
             CBs that respect this constraints are allowed to be push(ed)_in the same covered_bases_ensemble even if, usually, further
             condition are superimposed (in PROGRAM_CORE function we set also 'mutual distance between CBs =< Bushman bp rule')
    
    
    STRUCTURE:
    
        __init__ INPUT: - Covered_base_object: object of Covered_base class
                        - label_selection: a string like the ones in 'selective_reads_count' attribute of Covered_base_object, allowing
                                           ensemble creation looking only at one specific label. Default is 'all' (no label selection)
                        - strand_specific: boolean; if true (default), Covered_bases_ensamble construction accounts also for strand.
                                           usually 'strand_specific' should be considered as a kind of analysis so the choice made here
                                           should be inherited from the previous ones (e.g. strand_specific choice in CB construction)
                                           For this purpose, you can find a variable called 'strand_specific_choice' in main, retrieved 
                                           from user input, so the best usage is strand_specific = strand_specific_choice
                        
        ATTRIBUTES: 'label' - string, resulting from label_selection variable given in input
                    'chromosome' - string; the chromosome hosting the Covered_bases_ensamble
                    'strand' / 'strand_aspecific' - string; the strand hosting the Covered_bases_ensamble.
                                                    NOTE: always both present;'strand_aspecific' set 'None' and self.strand = Covered_base_object.strand.
                                                    Then, if strand_specific = False, self.strand_aspecific = Covered_base_object.strand_aspecific
                    'starting_base_locus' - long int; min value among CB's loci 
                    'ending_base_locus' - long int; max value among CB's loci 
                    'spanned_bases' - int; ending_base_locus - starting_base_locus + 1
                    'n_covered_bases' - int; number of CBs in the ensemble (len(Covered_bases_list))
                    'n_total_reads' - int; the sum of reads of all CBs in the ensemble (coherent with label_selection choice)
                    'covered_base_of_max' - Covered_base_object, the one in the ensemble with the highest reads count (coherent with label_selection choice)
                    'Covered_bases_list' - list of Covered_base_object(s) hosted in the ensemble
                    'IS_derived' - initialized as 'None', it will contain derived IS objects as list, after IS retrieval
                    
        METHODS: push_in - [...]
                 push_out - [...]
    
    NOTE for developers: - push_in method refreshes attributes in real time.
                         - push_out method was recently created specifically for gauss IS retrieval purpose: it needs further refinement and improvements
                           to be safely used out of its context
    '''
 
    #Constructor####################################################################################################### 
    def __init__(self, Covered_base_object, label_selection = "all", strand_specific = True):
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
        self.IS_derived = None # Not known 'a priori', it will become a list of IS object
    ####################################################################################################################
    
        
    #Methods############################################################################################################
    
    #Through this method you can push_in a Covered_base_object in an already existing Covered_bases_ensamble object
    #Controls about 'duplicate pushing' (Covered_base_object not in self.Covered_bases_list), about chromosome and about strand
    #are performed. 'label_selection' as in Constructor.
    #If the Covered_base_object is suitable for the ensemble, it will be added, attributes of the ensemble refreshed and 1 returned
    #Else nothing happens and -1 is returned
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
    
    #Recently created specifically for gauss IS retrieval purpose: it needs further refinement and improvements to be safely used out of that context
    #Note: works only with label_selection = "all", doesen't support removal of last CB (simply return -1)
    def push_out (self, Covered_base_object, label_selection = "all"): 
        check = -1
        
        if ((label_selection == "all") and (Covered_base_object in self.Covered_bases_list) and (len(self.Covered_bases_list) - 1 > 0)):
            self.Covered_bases_list.remove(Covered_base_object)
            if ((self.starting_base_locus == Covered_base_object.locus) or (self.ending_base_locus == Covered_base_object.locus)):
                list_of_loci = []
                for covered_base in self.Covered_bases_list:
                    list_of_loci.append(covered_base.locus)
                self.starting_base_locus = min(list_of_loci)
                self.ending_base_locus = max(list_of_loci)
                self.spanned_bases = self.ending_base_locus - self.starting_base_locus + 1
            self.n_covered_bases = self.n_covered_bases - 1
            self.n_total_reads = self.n_total_reads - Covered_base_object.reads_count    
            if (Covered_base_object == self.covered_base_of_max):
                cb_of_max = self.Covered_bases_list[0]
                for covered_base in self.Covered_bases_list:
                    if (covered_base.reads_count > cb_of_max.reads_count):
                        cb_of_max = covered_base
                self.covered_base_of_max = cb_of_max
            check = 1
        
        return check  
            
    ####################################################################################################################    
    
########################################################################################################################      
        


        
###Class of Integration Sites###########################################################################################         
class IS:
    '''
    *** Class of Integration Sites ***

    CONCEPT: avoiding any biological explanation about IS, here should be considered as a 'results collector'.
             ISs are derived from a Covered_bases_ensamble_object (the same CBE may produce one or more IS) hence by defaults
             they are characterized by the 'related_ensemble', therefore 'chromosome' and 'strand'/'strand_aspecific' too.
             All other attributes are fixed to None in order to be properly set during IS retrieval stage, in an easier way
             
    ALL FEATURES IN BRIEF:
    
    Input:
    You need a Covered_bases_ensamble_object to create an IS (the ensemble from which the IS is derived), strand_specific as above
    
    Attributes:
    Attributes names are intuitive hence no further specifications are needed. Characteristic of attributes should be obvious (e.g. type)
    but in any case unpredictable due to their complete dependency from the IS retrieval method used. It's up to you to make a good work
    (give a look to already implemented IS retrieval method before developing a new one!)
    
    NOTE for developers: 'selective_reads_count' should be a dic of kind {label:read_count}, likewise CB, in order to work smooth with functions for matrixes generation!             
    '''
 
    #Constructor####################################################################################################### 
    def __init__(self, Covered_bases_ensamble_object, strand_specific = True):
        '''
        [...]
        '''
        self.related_ensemble = Covered_bases_ensamble_object
        self.Covered_bases_list = None # different from Covered_bases_ensamble_object.Covered_bases_list, e.g. in Gauss IS retrieval method
        self.chromosome = Covered_bases_ensamble_object.chromosome
        self.strand = Covered_bases_ensamble_object.strand
        self.starting_base_locus = None #to be evaluated by Integration_Sites_retrieving_methods
        self.ending_base_locus = None #to be evaluated by Integration_Sites_retrieving_methods
        self.integration_locus = None #to be evaluated by Integration_Sites_retrieving_methods
        self.spanned_bases = None #to be evaluated by Integration_Sites_retrieving_methods
        self.n_covered_bases = None #to be evaluated by Integration_Sites_retrieving_methods
        self.reads_count = None #to be evaluated by Integration_Sites_retrieving_methods
        self.selective_reads_count = None # a dic of kind {label:read_count}, likewise CB, in order to work smooth with functions for matrixes generation!
        self.peak_height = None  #to be evaluated by Integration_Sites_retrieving_methods
        self.reads_key_list = None #to be evaluated by Integration_Sites_retrieving_methods
        self.strand_aspecific = None 
        if (strand_specific == False):
            self.strand_aspecific = Covered_bases_ensamble_object.strand_aspecific
        self.longest_seq_header = None
        self.longest_raw_seq = None
        self.longest_final_seq = None
    ####################################################################################################################
    
    #Methods############################################################################################################
    def track_sequences (self):
        integration_cb = None
        for cb in self.Covered_bases_list:
            if (cb.locus == self.integration_locus):
                integration_cb = cb            
        self.longest_seq_header=cb.longest_seq_header
        self.longest_raw_seq=cb.longest_raw_seq
        self.longest_final_seq =cb.longest_final_seq
    ####################################################################################################################
    
    
        
########################################################################################################################         
