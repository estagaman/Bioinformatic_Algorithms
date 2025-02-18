#given a genus name and two dates in YYYY/MM/DD format 
#find number of nucleotide database entries on that genus published between dates 

import sys
import argparse
#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="prob1.py") #specify the name of the script
    parser.add_argument("-i", "--input", #add input file argument
    help="input file",
    required=True)
    parser.add_argument("-o", "--output", #add output file argument
    help="output file",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments 
arguments = check_arg(sys.argv[1:]) #get the arguments that I inputted from the terminal
infile = arguments.input #assign input to "infile" variable
outfile = arguments.output #assign output location to "outfile" variable for later

from Bio import Entrez #import Entrez so we can search NCBI

file_list = open(infile).readlines() #make a list with each line of the file as one element
genus = file_list[0] #first line is the genus we are searching
date_1 = file_list[1] #second line is the starting date
date_2 = file_list[2] #third line is the ending date

search_term = genus + " [ORGANISM]" #making sure ["ORGANISM"] specification is included when we search

Entrez.email = "elisestagaman@gmail.com" #give NCBI an email to contact

#perform the search: term is the genus, datetype shows we're looking at publication date, and mindate and maxdate are our date range
handle = Entrez.esearch(db="nucleotide", term = search_term, datetype = 'pdat', mindate = date_1, maxdate = date_2) 
record = Entrez.read(handle)

with open(outfile, "w") as output_handle: #write to output file 
    output_handle.write(record["Count"]) #save the number of records returned that match search criteria
