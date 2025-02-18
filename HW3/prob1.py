#argparse 
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

#import packages
from Bio.Seq import Seq
import pandas as pd

with open(infile) as f: #open the file line by line
    lines = [line.rstrip() for line in f]

revlines = [] #start list of reverse complements

for line in lines: #for each line
    line = Seq(line) #convert to sequence object
    revcomp = line.reverse_complement() #find reverse complement
    revlines.append(str(revcomp)) #add rev comp to the sequence list

all_seqs = lines + revlines #add seqs and their reverse complements to make one list

adjac = pd.DataFrame(columns = ["seq1", "seq2"]) #begin a dataframe for adjacent kmers

for seq in all_seqs: #for each sequence in the list
    adjac.loc[len(adjac)] = [seq[0:3], seq[1:4]] #add a row to the dataframe with adjacencies

adjac = adjac.drop_duplicates().sort_values(by = ["seq1", "seq2"]) #remove duplicate rows and sort alphabetically

with open(outfile, "w") as output_handle: #write to output file 
    for i in range(len(adjac)): #for each row in the data frame 
        output_handle.write(adjac["seq1"].iloc[i] + " " + adjac["seq2"].iloc[i] + "\n") #write the adjacent kmers separated by a space
