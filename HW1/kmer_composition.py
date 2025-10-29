#kmer_composition.py: find the kmer composition of a sequence based on a user-specified kmer length
#INPUT: k-mer length as the first line, sequence on second line 
#OUTPUT: txt file with the count of each unique kmer

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

from itertools import product #using this tool to calculate all possible kmer combinations

file_list = open(infile).readlines() #make a list with each line of the file as one element
count = int(file_list[0]) #the kmer length is the first line
seq = file_list[1] #the sequence is the second line

kmer_list = [] #start a list of all possible kmers you could make from that sequence
for i in range(len(seq) - (count - 1)): #iterate through the sequence
    kmer = seq[i:(i + count)] #each kmer starts at i, count from there until kmer length is reached
    kmer_list.append(kmer) #add each kmer to a list

all_combos = product("ACGT", repeat = count) #find all possible combinations of bases for that kmer length
possible_kmers = [] #make it into a list
for combo in all_combos: #for each combination of bases, make it a string (product function returns a list)
    kmer_str = "".join(combo)
    possible_kmers.append(kmer_str)

total_counts = [] #calculating the number of times each kmer appears in the sequence
for i in range(len(possible_kmers) - 1): #iterate through our list of all possibilities
    kmer = possible_kmers[i]
    kmer_count = kmer_list.count(kmer) #count how many times that kmer appears in the sequence
    total_counts.append(kmer_count) #add this to our list of counts

total_counts = [str(x) for x in total_counts] #convert all the counts to strings so we can format them for the outfile

with open(outfile, "w") as output_handle: #write to output file 
    output_handle.write(" ".join(total_counts)) #join with " " so each count separated by a space
