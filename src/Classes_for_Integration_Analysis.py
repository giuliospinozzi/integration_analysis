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


###Class of read sequences###########################################################################
class Read:
    '''
    Class of read sequences
    [...]
    '''

    def __init__(self, read_header, read_attributes):
        '''
        Constructor of Read objects
        [...]
        '''
        #Legenda
        #read_attributes = ()
        #read_attributes[0] = reference_genome
        #read_attributes[1] = chromosome
        #read_attributes[2] = strand
        #read_attributes[3] = integration_locus
        #read_attributes[4] = span
        #read_attributes[5] = lam_id
        
        #Read object attributes
        self.header = read_header
        self.genome = read_attributes[0]
        self.chromosome = read_attributes[1]
        self.strand = read_attributes[2]
        self.start = read_attributes[3]
        self.end = read_attributes[3] + read_attributes[4] #read_end = integration_locus + span
        self.lam_id = read_attributes[5]        
#####################################################################################################

###Class of covered bases######################       
class Covered_base:
    '''
    Class of covered bases.
    [...]
    '''
 
 
    def __init__(self, Read_object):
        '''
        [...]
        '''
        self.list_of_reads = [Read_object]
        self.chromosome = Read_object.chromosome
        self.strand = Read_object.strand
        self.locus = Read_object.start
        self.reads_count = 1
         
    def add (self, Read_object):
        '''
        [...]
        '''
        self.list_of_reads.append(Read_object)
        self.reads_count = self.reads_count + 1
         
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