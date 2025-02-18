#problem 2

#input: one DNA-sequencing read per line 
#S is the set of all possible kmers in the reads and their reverse complements (adjacency list but one column)
#for S, there is a value k 
#want De Brujn graph of left and right k-1-mers to generate 2 directed cycles

#write both cyclic superstrings to an output file 

#tips: 
#iterate to find the right value of k 
#when the suffix of one edge matches the prefix of another edge --> add last character of second edge to superstring
#superstrings could rotate example, they are cyclic
#can be present in the cycle, but not the read, make sure to wrap around 
#for each value of k: 
    #start with read length
    #build the k-1 DeBruijn Graph - try to build exactly 2 directed cycles 
    #break when 2 directed cycles are created

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


from Bio.Seq import Seq
import pandas as pd

def make_seq_list(infile):
    with open(infile) as f:
        lines = [line.rstrip() for line in f]

    revlines = []

    for line in lines: 
        line = Seq(line)
        revcomp = line.reverse_complement()
        revlines.append(str(revcomp))

    all_seqs = lines + revlines  
    return all_seqs

def make_adj_list(seq_list, k): #need to adjust this to fit multiple values of k
    adjac = pd.DataFrame(columns = ["seq1", "seq2"])
    for seq in seq_list:
        num_kmers = len(seq) - k
        for j in range(num_kmers):
            adjac.loc[len(adjac)] = [seq[j:j+k], seq[j+1:j+1+k]]

    adjac = adjac.drop_duplicates().sort_values(by = ["seq1", "seq2"])
    return adjac

def test_graph(seq_list):
    k = len(seq_list[0])
    while k != 0:
        adjac = make_adj_list[seq_list, k]
        start_index = 0
        start = adjac["seq1"].iloc[start_index] #starting with the first kmer in top left
        next_try = adjac["seq2"].iloc[start_index]
        list1 = [start, next_try]
        while next_try in adjac["seq1"] and next_try not in list1:
            indices = adjac.index[adjac["seq1"] == next_try].tolist()
            next_try = adjac["seq2"].iloc[indices]
            list1.append(next_try)
        
        row = adjac["seq1"].iloc[0]
        i = 0

        while row in list1 and i != adjac.shape[0]: #iterate through rows in the dataframe to kmer not in first cycle
            row = adjac["seq1"].iloc[i] #look at each row 
            i += 1

        #now we have a row that is not in list 1 with index = i - 1  
        #start our second loop   
        start2 = adjac["seq1"].iloc[i - 1] #starting with the first kmer in top left
        next_try2 = adjac["seq2"].iloc[i - 1]
        list2 = [start2, next_try2]
        while next_try2 in adjac["seq1"] and next_try2 not in list2:
            indices = adjac.index[adjac["seq1"] == next_try2].tolist()
            next_try2 = adjac["seq2"].iloc[indices]
            list2.append(next_try2)

        #created our second loop 

        #now we check if there is anything left in the dataframe not included in the lists
        




def find_superstrings(infile):
k = 3




seqs = make_seq_list(infile)
adj_list = make_adj_list(seqs, k)

with open(outfile, "w") as output_handle: #write to output file 
    for i in range(len(adj_list)): 
        output_handle.write(adj_list["seq1"].iloc[i] + " " + adj_list["seq2"].iloc[i] + "\n")

import graphviz as gv
dbg = gv.Digraph('dBgraph', format='png')
with open("prob2_output_small") as f:
    for line in f:
        e1, e2 = line.strip().split()
        dbg.edge(e1, e2)
dbg.render()
