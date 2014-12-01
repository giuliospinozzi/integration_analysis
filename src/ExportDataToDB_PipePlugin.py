#!/usr/bin/python
# -*- coding: utf-8 -*-
import MySQLdb, sys, os, argparse, csv
from operator import itemgetter
import warnings

parser = argparse.ArgumentParser()

parser.add_argument('-b', '--bedfile', dest="bedfile", help="BED file to process. No default option.", action="store", required=True)
parser.add_argument('-a', '--associationfile', dest="associationfile", help="File of association: for each barcode this file (TSV format) specifies TAG sequence, tissue, sample and time point, lam id, lam name, marker, enzyme, vector. No default option.", action="store", required=True)
#parser.add_argument('--headerformat', dest="headerformat", help="Header format of the filter file. Cases: type 1:: '@M00174:25:000000000-A0T21:1:1:15041:1491 1:N:0:0'; type 2:: '@M00571:5:000000000-A267F:1:1101:14475:1610/1'. Select type number {1, .., N}. No default value assigned because you must explicitely be aware of what you are doing.", action="store", required=True)
parser.add_argument('-r', '--seqfileraw', dest="seqfileraw", help="File of raw sequences, extracted with a custom script of Andrea that formats for each read ID, as reported in the BED file, it stores raw sequences. Default is 'None'.", action="store", default=None)
parser.add_argument('-t', '--seqfiletrimmed', dest="seqfiletrimmed", help="File of trimmed sequences, extracted with a custom script of Andrea that formats for each read ID, as reported in the BED file, it stores trimmed sequences. Default is 'None'.", action="store", default=None)
parser.add_argument('--patient', dest="patient", help="Patient ID. No default option.", action="store", required=True)
parser.add_argument('--pool', dest="pool", help="Pool ID. No default option.", action="store", required=True)
parser.add_argument('--tag', dest="tag", help="TAG ID. No default option.", action="store", required=True)
parser.add_argument('--host', dest="host", action="store", required=True)
parser.add_argument('--user', dest="user", action="store", required=True)
parser.add_argument('--passwd', dest="passwd", action="store", required=True)
parser.add_argument('--port', dest="port", action="store", required=True, type=int)
parser.add_argument('--dbschema', dest="dbschema", help="Target Database schema in which importing data. E.g.: MLD01. Mandatory.", action="store", required=True)
parser.add_argument('--dbtable', dest="dbtable", help="Target Database table in which importing data. E.g.: redundant_MLD01_ALL_NEW_NAMES. Mandatory.", action="store", required=True)

args = parser.parse_args()

########################################################################
####### GLOBAL VARS
########################################################################

# init db values
userdb = {
	'host': args.host,
	'user': args.user,
	'passwd': args.passwd,
	'port': args.port,
	'db': args.dbschema,
	}

conn = MySQLdb.connect( 
	host = userdb["host"],
	user = userdb["user"],
	passwd = userdb["passwd"],
	port = userdb["port"],
	db = userdb["db"]
	)

cursor = conn.cursor (MySQLdb.cursors.DictCursor)


########################################################################
####### FUNCTIONS
########################################################################

def parseAssociationfile(infile):
	"""
	IN: association file
	OUT: dictionary of TAGs (k=tag, v=dict: k={tissue, sample, timepoint}, v=values)
	FILE FORMAT: <barcode_id> <barcode_sequence> <tissue> <sample> <timepoint> <lam_id> <lam_fullname> <marker> <enzyme> <vector>
	"""
	assodict = {}
	with open(infile, 'rb') as inf:
		reader = csv.reader(inf, delimiter = "\t")
		for row in reader:
			if len(row) > 0:
				barcode_id = row[0] 
				barcode_sequence = row[1]
				tissue = row[2]
				sample = row[3]
				timepoint = row[4]
				assodict[barcode_sequence] = {
					"tissue": tissue,
					"sample": sample,
					"timepoint": timepoint,
					"barcode_id": barcode_id,
					"lam_id": row[5],
					"lam_fullname": row[6],
					"marker": row[7],
					"enzyme": row[8],
					"vector": row[9],
					}
			else:
				print "[AP] WARNING: this file contais a blank row!!!! Check it please."
	return assodict

def parseBEDfile(inbedfile):
	"""
	IN: bed file
	OUT: BED tuple array with chr, is, header, strand, score
	LOGICS: 	
		1. recognize start point from BED
		2. create tuple array of BED elements
	"""
	bedarray = [] # return var
	with open(inbedfile, 'rb') as inf:
		reader = csv.reader(inf, delimiter = "\t")
		for row in reader:
			chr = row[0].replace('chr', '')
			if row[5] == '+':
				bedarray.append( (chr, row[1], row[3], row[5], row[4]) )
			else:
				bedarray.append( (chr, row[2], row[3], row[5], row[4]) )
	results = sorted(bedarray, key=itemgetter(0, 1), reverse = False) # sorted by chr, start, end -> now do sort by orientation!!!
	return results
	
def queryIfTableExists(dbschema, dbtable):
	"""
	IN: db options
	OUT: boolean (true/false)
	LOGICS: 	
		1. query info schema to see if table exists
		2. return response
	"""
	exists = False
	query = """
		SELECT *
		FROM information_schema.tables
		WHERE table_schema = '%s'
			AND table_name = '%s';
		""" %(dbschema, dbtable)
	cursor.execute(query)
	results = cursor.fetchall()
	if len(results) == 1:
		exists = True
	return exists

def DropTableIfExists(dbschema, dbtable):
	with warnings.catch_warnings():
		warnings.filterwarnings('ignore', 'unknown table')
		query = "DROP TABLE IF EXISTS `{0}`.`{1}`".format(str(dbschema), str(dbtable))
		cursor.execute(query)
	return True

def createTargetTable(dbschema, dbtable):
	"""
	IN: db options
	OUT: -
	LOGICS: if target table does NOT exist, create table
	"""
	query = """
		CREATE TABLE IF NOT EXISTS `%s`.`%s` (
		  `ID` int(11) NOT NULL AUTO_INCREMENT,
		  `group_name` varchar(200) DEFAULT NULL,
		  `n_LAM` varchar(200) DEFAULT NULL,
		  `pool` varchar(200) DEFAULT NULL,
		  `tag` varchar(20) DEFAULT NULL,
		  `sample` varchar(200) DEFAULT NULL,
		  `vector` varchar(50) DEFAULT NULL,
		  `tissue` varchar(50) DEFAULT NULL,
		  `treatment` varchar(20) DEFAULT NULL,
		  `enzyme` varchar(50) DEFAULT NULL,
		  `complete_name` varchar(200) DEFAULT NULL,
		  `header` varchar(200) DEFAULT NULL,
		  `sequence_raw` varchar(2000) DEFAULT NULL,
		  `sequence_trimmed` varchar(2000) DEFAULT NULL,
		  `chr` varchar(50) DEFAULT NULL,
		  `integration_locus` int(11) DEFAULT NULL,
		  `sequence_count` int(3) NOT NULL DEFAULT '1',
		  `start_seq` int(11) DEFAULT NULL,
		  `end_seq` int(11) DEFAULT NULL,
		  `score` double DEFAULT NULL,
		  `identity` double DEFAULT NULL,
		  `e_value` double DEFAULT NULL,
		  `strand` varchar(2) DEFAULT NULL,
		  `q_size` int(11) DEFAULT NULL,
		  `span` int(11) DEFAULT NULL,
		  `transcript_ucsc` varchar(20) DEFAULT NULL,
		  `orientation_ucsc` varchar(20) DEFAULT NULL,
		  `distance_to_TSS_ucsc` int(11) DEFAULT NULL,
		  `in_gene_ucsc` double DEFAULT NULL,
		  `upstream_end_ucsc` int(11) DEFAULT NULL,
		  `downstream_end_ucsc` int(11) DEFAULT NULL,
		  `gene_lenght_ucsc` int(11) DEFAULT NULL,
		  `ref_gene_ucsc` varchar(200) DEFAULT NULL,
		  `ref_seq_name_ucsc` varchar(200) DEFAULT NULL,
		  `entrez_gene` varchar(200) DEFAULT NULL,
		  `entrez_name` varchar(200) DEFAULT NULL,
		  `ensembl_gene` varchar(200) DEFAULT NULL,
		  `refSeqName` varchar(200) DEFAULT NULL,
		  `refSeqSymbol` varchar(200) DEFAULT NULL,
		  `gene_orientation_refseq` varchar(2) DEFAULT NULL,
		  `ingene_refseq` double DEFAULT NULL,
		  `distance_to_TSS_refseq` int(11) DEFAULT NULL,
		  `upstream_end_refseq` int(11) DEFAULT NULL,
		  `downstream_end_refseq` int(11) DEFAULT NULL,
		  `gene_lenght_refseq` int(11) DEFAULT NULL,
		  `miRNA` varchar(200) DEFAULT NULL,
		  `orientation_mirna` varchar(200) DEFAULT NULL,
		  `inmirna` double DEFAULT NULL,
		  `distance_to_TSS_mirna` int(11) DEFAULT NULL,
		  `upstream_end_mirna` int(11) DEFAULT NULL,
		  `downstream_end_mirna` int(11) DEFAULT NULL,
		  `lenght_mirna` int(11) DEFAULT NULL,
		  PRIMARY KEY (`ID`),
		  KEY `integration_locus` (`integration_locus`),
		  KEY `header` (`header`)
		) ENGINE=MyISAM DEFAULT CHARSET=latin1 AUTO_INCREMENT=1 ;
		""" %(dbschema, dbtable)
	cursor.execute(query)
	return True

def aquireSeqData(infile):
	"""
	IN: seqfile formatted as <id> <sequence>
	OUT: dictionary of data from CSV file
	LOGICS: 
	"""
	if not os.path.isfile(infile):
		print "[ERROR]\tNo valid input files. Check them please then run me again.\n"
		sys.exit()
	else:
		dicseq = {}
		with open(infile, 'rb') as inf:
			reader = csv.reader(inf, delimiter = "\t")
			for row in reader:
				if len(row)>1:
					dicseq[row[0]] = row[1].strip()
		return dicseq

######################################################					


#print "[AP]\tChecking input files"
if not os.path.isfile(args.bedfile) or not os.path.isfile(args.associationfile):
	print "[ERROR]\tNo valid input files. Check them please then run me again.\n"
	sys.exit()

#print "[AP]\tGet Association data"
dict_asso = parseAssociationfile(args.associationfile)

#print "[AP]\tGet BED data"
beddata = parseBEDfile(args.bedfile)

#print "[AP]\tDrop target Table if exists"
DropTableIfExists(args.dbschema, args.dbtable)

#print "[AP]\tCreate target Table if not exists"
createTargetTable(args.dbschema, args.dbtable)

## in case of having or not sequence file, change import method
if args.seqfileraw is not None and args.seqfiletrimmed is not None:
	#print "[AP]\tGet SEQUENCE data (raw and trimmed)"
	dict_seq_raw = aquireSeqData(args.seqfileraw)
	dict_seq_trimmed = aquireSeqData(args.seqfiletrimmed)
	
	#print "[AP]\tImport data into target Table"
	### now format data to import all and insert them all
	query_import = """
		INSERT INTO `%s`.`%s` (`group_name`, `n_LAM`, `pool`, `tag`, `sample`, `tissue`, `treatment`, `enzyme`, `complete_name`, `header`, `chr`, `integration_locus`, `strand`, `score`, `vector`, `sequence_raw`, `sequence_trimmed`) 
		VALUES 
		""" %(args.dbschema, args.dbtable)
	# write array of data in this order:: patient, lamid, pool, tag, sample, tissue, timepoint, enzyme, lam_name, header, chr, is, strand, score
	for (chr, locus, header, strand, score) in beddata:
		query_import += " ('" + "', '".join( str(x) for x in  
			[args.patient, 
			dict_asso[args.tag]["lam_id"], 
			args.pool, 
			args.tag, 
			dict_asso[args.tag]["sample"], 
			dict_asso[args.tag]["tissue"], 
			dict_asso[args.tag]["timepoint"],
			dict_asso[args.tag]["enzyme"],
			dict_asso[args.tag]["lam_fullname"],
			header,
			chr,
			locus,
			strand,
			score,
			dict_asso[args.tag]["vector"],
			dict_seq_raw[header],
			dict_seq_trimmed[header]		 ]
		  ) + "'), "
	#print query_import.rstrip(', ')
	cursor.execute( query_import.rstrip(', ') )

else: # caso of not available sequence file
	#print "[AP]\tImport data into target Table"
	### now format data to import all and insert them all
	step = 5000
	bedstart = 0 # starting point of bed data array for each step
	bedend = 0 # starting point of bed data array for each step
	# check max len of beddata to set bedend, no more than beddata len!
	if len(beddata) > step:
		bedend = step
	else:
		bedend = len(beddata)
	# now process data into steps of a maximum number of blocks (step)
	for elems in range(0, len(beddata), step):
		query_import = """
			INSERT INTO `%s`.`%s` (`group_name`, `n_LAM`, `pool`, `tag`, `sample`, `tissue`, `treatment`, `enzyme`, `complete_name`, `header`, `chr`, `integration_locus`, `strand`, `score`, `vector`) 
			VALUES 
			""" %(args.dbschema, args.dbtable)
		# write array of data in this order:: patient, lamid, pool, tag, sample, tissue, timepoint, enzyme, lam_name, header, chr, is, strand, score
		for (chr, locus, header, strand, score) in beddata[bedstart:bedend]:
			query_import += " ('" + "', '".join( str(x) for x in  
				[args.patient, 
				dict_asso[args.tag]["lam_id"], 
				args.pool, 
				args.tag, 
				dict_asso[args.tag]["sample"], 
				dict_asso[args.tag]["tissue"], 
				dict_asso[args.tag]["timepoint"],
				dict_asso[args.tag]["enzyme"],
				dict_asso[args.tag]["lam_fullname"],
				header,
				chr,
				locus,
				strand,
				score,
				dict_asso[args.tag]["vector"] ]
			  ) + "'), "
		#print query_import.rstrip(', ')
		cursor.execute( query_import.rstrip(', ') )
		
		# increment bed indexes
		bedstart += step
		bedend += step

#print "[AP]\tTask finished\n\tCheck target table\n\n"


######################## finalize objects ##################		
cursor.close()
conn.close()