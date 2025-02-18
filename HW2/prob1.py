#dynamic programming 
#tips: use numPy or make a list of lists 
#calculate the edit distance between two sequences 
#first function to calculate the alignment scores 
#second function starts in the lower corner to tally up the whole score and build the alignment 

#argparse 
import sys
import argparse
#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="prob2.py") #specify the name of the script
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

from Bio import SeqIO #helps take in the fasta file (SeqIO.parse())

#open the fasta file
with open(infile) as handle:
    all_seq = list(SeqIO.parse(handle, "fasta")) #separate out each sequence
    
    seq1 = str(all_seq[0].seq) #mark just sequence 1
    seq2 = str(all_seq[1].seq) #mark just sequence 2

import numpy as np

def diag_left_above(x, y, matrix): #function that finds the values diagonal, left, and above where we currently are
    diag = matrix[x-1, y-1] #what value is diagonal
    above = matrix[x-1, y] #what value is above
    left = matrix[x, y-1] #what value is to the left
    return(diag, left, above)

def create_distance_matrix(seq1, seq2): #create the edit distance matrix

    my_matrix = np.full([len(seq1) + 1, len(seq2) + 1], np.nan) #create an empty numpy matrix full of NAs

    my_matrix[0, :] = np.array(range(len(seq2) + 1)) #make the first column numbered 1,2,3,4 etc. until length of sequence 2 is reached
    my_matrix[:, 0] = np.array(range(len(seq1) + 1)) #make the first row numbered 1,2,3,4,etc. until length of sequence 1 is reached

    for i in range(1,len(seq1) + 1): #for every letter in sequence 1
        for j in range(1,len(seq2) + 1): #compare to every letter in sequence 2
            if seq1[i-1] == seq2[j-1]: #if it's a match
                my_matrix[i, j] = my_matrix[i-1, j-1] #value stays the same as the diagonal point, no adding/subtracting
            else: #if it's not a match
                diag, left, above = diag_left_above(i,j, my_matrix) #find the values diagonal, left, and above

                take = min(diag, above, left) #take the minimum of these 3
                my_matrix[i, j] = take + 1 #add 1 to it as penalty for mismatch or gap, assign to spot in matrix

    return my_matrix.T #returning transposed version because this is more consistent with the paper example

def find_best_alignment(dist_matrix): #find the edit score of a distance matrix 

    total_score = dist_matrix[-1,-1] #bottom left value 

    return total_score #return that score

distance_matrix = create_distance_matrix(seq1, seq2) #create the matrix 

total_score = find_best_alignment(distance_matrix) #find the best score 

with open(outfile, "w") as output_handle: #write to output file 
        output_handle.write(str(int(total_score))) #convert the score to an integer, then a string for saving 
