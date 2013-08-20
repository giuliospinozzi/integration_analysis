'''
Created on Aug 20, 2013

@author: Andrea
@note: this module aims at: (1) dump db data into a file CSV with header using a system call, (2) parse file using CSV module, (3) return same data structure as DB_connection.import_data_from_DB 
'''

import sys, os, csv
from time import gmtime, strftime

def dbTableDump (host, user, passwd, db, db_table, dest_folder, query_select_statement):
    """
    Input: same parameters as DB_connection.py -> db parameters, destination folder
    Output: file name (with absolute path) -> delimiter is by default \t
    Required: mysql client installed and globally named
    """
    tmpdate = strftime("%Y-%m-%d %H:%M:%S", gmtime())
    destfile = os.path.join(dest_folder, "tmp.%s.%s.tsv" %(db_table, tmpdate, ) )
    os.system( "mysql -h %(host)s -u %(user)s --password=%(passwd)s %(db)s -e 'SELECT %(query_select_statement)s FROM %(db_table)s WHERE 1' > %(destfile)s " %{
     'host': host,
     'user': user,
     'passwd': passwd,
     'db': db,
     'db_table': db_table,
     'destfile': destfile,
     'query_select_statement': query_select_statement,
     }
    )
    return destfile


def parseCSVdumpFile (infile, key_field, array_field, delimiter = "\t"):
    """
    Input: CSV/TSV file with delimiter \t
    Output: dictionary where key = key_field, v = tuple of array_field => data structure same as reads_query, lam_query (same as DB_connection.import_data_from_DB)
    Logics: (1) acquire header and create dictionary of positions (this structure will be used to replicate the tuple data structure from the fetchall call to MySQL; (2) acquire data and return data structure of dictionary in which key is key_field and value is a tuple of values ordered by array_field 
    """
    # init vars
    dict_header = {} # dictionary of header fields: k = field, v = number
    filerow_index = 0 # file row index
    dict_filecontent = {} # dictionary of the file content: k = key_field, v = tuple of fields sorted as array_field
    # get header opening the file
    with open(infile, 'rb') as inf: # parse input file line by line
        reader = csv.reader(inf, delimiter = delimiter) # using csv
        try:
            for row in reader:
                if len(row) > 0: # just a control on row len
                    if filerow_index == 0: # you are reading the header line
                        for col_number in range(0, len(row)): # parse all fields of the columns header
                            dict_header[row[col_number]] = col_number
                        for k in array_field: # now check that all array_field are also in the header keys 
                            if k not in dict_header.keys():
                                print "[AP]\tError:\tinput file %s does not contain the field %s" %(infile, k)
                                sys.exit()
                    else: # all other lines, the content of the file
                        if not dict_header.has_key(key_field): # just check key field in dictionary (because later you will use it as key)
                            print "[AP]\tError:\tthe KEY field %s is not in the dictionary of acquired columns from the dump file (all keys: %s)" %(key_field, dict_header.keys())
                            sys.exit()
                        for col_number in range(0, len(row)): # parse all fields of the columns header
                            sorted_row_array = [] # init of array of the row content in the final sort
                            for k in array_field: # parse content in this oder
                                if row[dict_header[k]].isdigit(): # check data type (int, float, string)
                                    sorted_row_array.append(int(row[dict_header[k]]))
                                elif row[dict_header[k]].isfloat():
                                    sorted_row_array.append(float(row[dict_header[k]]))
                                else:
                                    sorted_row_array.append(row[dict_header[k]])
                                ### sorted_row_array = [x for x in row[dict_header[k]]] -> this is good unless you want to change type in the array
                            dict_filecontent[row[dict_header[key_field]]] = tuple(sorted_row_array)
                else:
                    print "[AP]\tWarning:\tinput file %s has got the line number %d blank" %(infile, filerow_index)
                filerow_index += 1
        except csv.Error, e:
            sys.exit("[AP]\tError while reading file %s, line %d: %s" % (infile, reader.line_num, e))
    # return data structure (dictionary)
    return dict_filecontent






