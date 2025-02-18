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

from Bio.Seq import Seq #import Seq from Biopython
import pandas as pd # import pandas package 

def make_seq_list(infile): #function to make a list of all seqs and reverse complements
    with open(infile) as f:
        lines = [line.rstrip() for line in f]

    revlines = []

    for line in lines:  #for each sequence in the file
        line = Seq(line) #make into a sequence object
        revcomp = line.reverse_complement() #find reverse complement
        revlines.append(str(revcomp)) #add to list of reverse complements

    all_seqs = lines + revlines #adds seqs and complements together into one list
    return all_seqs

def make_adj_list(seq_list, k): #make an adjacency list for specified value of k
    adjac = pd.DataFrame(columns = ["seq1", "seq2"]) #make empty data frame for adjacencies
    for seq in seq_list: #for each read in our list
        num_kmers = len(seq) - k #number of kmers expected from one read
        for j in range(num_kmers): #find each possible kmer starting at index j
            adjac.loc[len(adjac)] = [seq[j:j+k], seq[j+1:j+1+k]] #add that kmer pair to our dataframe as a row

    adjac = adjac.drop_duplicates().sort_values(by = ["seq1", "seq2"]) #drop duplicate pairings and sort alphabetically
    adjac.index = range(0, len(adjac)) #adjust the indices back since we sorted

    return adjac #return the adjacency list for de bruijn graph testing 

def test_graph(seq_list): #tests whether we create two directed cycles for each value of k
    k = len(seq_list[0]) - 1 #highest possible value of k 
    stop = k #stop signal for our while loop
    while stop != 0:
        adjac = make_adj_list(seq_list, k) #make an adjacency list for k
        start_index = 0 #start with first row in our dataframe
        start = adjac["seq1"].iloc[start_index] #starting with the first kmer in top left
        next_try = adjac["seq2"].iloc[start_index] #find the adjacent kmer in top right
        list1 = [start] #add the first kmer to a list
        indices = [start_index] #index we're currently at is 0 or start_index
        while len(indices) > 0 and next_try != start: #as long as we're still finding matches for the adjacent kmer and we haven't completed a cycle:
            list1.append(next_try) #add the adjacent kmer to our list
            indices = adjac.index[adjac["seq1"] == next_try].tolist() #check to see if the column 2 kmer is anywhere in column 1
            if len(indices) > 0: #if we found a match
                next_try = adjac["seq2"].iloc[indices[0]] #check the adjacent kmer for that match

        row = adjac["seq1"].iloc[0] #start over at the beginning of the dataframe
        i = 0
        while row in list1 and i < adjac.shape[0]: #iterate through rows in the dataframe to find kmer that wasn't in first cycle
            row = adjac["seq1"].iloc[i] 
            i += 1

        #now we have a row that is not in list 1 with index = i - 1  
        #start our second loop
        start_index2 = i - 1 
        start2 = adjac["seq1"].iloc[i - 1]
        next_try2 = adjac["seq2"].iloc[i - 1]
        list2 = [start2] #initialize a new list to hold the second loop
        indices2 = [start_index2]
        while len(indices2) > 0 and next_try2 != start2: #functions the same as the first loop
            list2.append(next_try2)
            indices2 = adjac.index[adjac["seq1"] == next_try2].tolist()
            if len(indices2) > 0:
                next_try2 = adjac["seq2"].iloc[indices2[0]]

    #created our second loop
    #now we check if there is anything left in the dataframe not included in the lists
        #if there are pairings left in the dataframe, we don't have two directed cycles, set k = k - 1 and try again
        if len(list1 + list2) < len(adjac): # if not all the kmers have been used
            k = k - 1
            stop = k
        else: 
            stop = 0 #exit the loop if they've all been used

    return(list1, list2, k) 

def find_superstring(my_list, k): #make the superstrings from a kmer list
    ss = my_list[0] #start with first kmer
    stop_point = -k + 1 #correcting for loop: tells when to stop adding letters on
    for word in my_list[1:stop_point]:
        ss = ss + word[-1] #add the last letter of each kmer onto our string
    
    return ss

seqs = make_seq_list(infile) #make the sequence list from our file
list1, list2, k = test_graph(seqs) #make kmer lists

ss1 = find_superstring(list1, k) #make superstrings from those kmer lists
ss2 = find_superstring(list2, k)

with open(outfile, "w") as output_handle: #write to output file 
    output_handle.write(ss1 + "\n" + ss2)
