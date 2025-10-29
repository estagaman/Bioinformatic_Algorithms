#build_alignment.py: build a Needleman-Wunsch global alignment for two sequences 

#INPUT: fasta file with one sequence per line (total sequences)
#OUTPUT: txt file with alignment score on first line, sequence 1 alignment on second line, and sequence 2 alignment on third line

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

def diag_left_above(x, y, matrix): 
    diag = matrix[x-1, y-1] #what value is diagonal
    above = matrix[x-1, y] #what value is above
    left = matrix[x, y-1] #what value is to the left
    return(diag, left, above)

def create_distance_matrix(seq1, seq2): 

    my_matrix = np.full([len(seq1) + 1, len(seq2) + 1], np.nan) #create empty matrix with sequence 1 on one axis and sequence 2 on the other

    my_matrix[0, :] = np.array(range(len(seq2) + 1)) #number the axes 
    my_matrix[:, 0] = np.array(range(len(seq1) + 1))

    for i in range(1,len(seq1) + 1): #for each letter in sequence 1
        for j in range(1,len(seq2) + 1): #match to each letter in sequence 2
            if seq1[i-1] == seq2[j-1]: #if it's a match
                my_matrix[i, j] = my_matrix[i-1, j-1] #value stays the same as the diagonal value
            else: #if it's not a match
                diag, left, above = diag_left_above(i,j, my_matrix) #find the values diagonal, left, and above

                take = min(diag, above, left) #take the minimum of these 3 values
                my_matrix[i, j] = take + 1 #add 1 to it as penalty for mismatch or gap, assign to spot in matrix

    return my_matrix.T #returning transposed version because that follows the paper example 

def find_best_alignment(seq1, seq2, dist_matrix): #function to find the alignment once we have a distance matrix
    seq1 = list(seq1) #make a list of each letter in seq1
    seq2 = list(seq2) #make a list of each letter in seq2
    seq1_align = [] #initialize a list to add the alignment to 
    seq2_align = [] #initialize a list to add the alignment to 

    edit_dist = dist_matrix[-1, -1] #find the edit distance in the bottom left corner 

    nrow = len(dist_matrix) #how many rows are in our matrix? 
    ncol = len(dist_matrix[0]) #how many columns are in our matrix?

    y = ncol - 1 #use these values to start from the bottom left corner 
    x = nrow - 1

    while x > 0 and y > 0: #unti we reach the axes 
        diag, left, above = diag_left_above(x, y, dist_matrix) #find the values diagonal, left, and above our current coordinates 
        min_edit = min(diag, above, left) #when computing distance, we took from the minimum of these 3
        
        if len(seq1) > 0 and len(seq2) > 0: #as long as both sequences still have letter to compare:
            if diag == min_edit and diag == dist_matrix[x, y]: #it's a match
                x = x - 1 #move into the diagonal spot 
                y = y - 1
                seq1_align.insert(0, seq1[-1]) #add that letter to the beginning of each alignment list
                seq2_align.insert(0, seq2[-1])
                seq1.pop() #remove that letter from sequence 1
                seq2.pop() #remove that letter from sequence 2
            elif diag == min_edit:   #it's a mismatch 
                x = x - 1 #move into the diagonal spot 
                y = y - 1
                seq1_align.insert(0, seq1[-1]) #add that letter to the beginning of each alignment list
                seq2_align.insert(0, seq2[-1])
                seq1.pop() #remove that letter from sequence 1
                seq2.pop()  #remove that letter from sequence 2
            elif above == min_edit:  #it's a gap 
                x = x - 1 #move over one cell, using x because I transposed the matrix
                seq2_align.insert(0, seq2[-1]) #insert the letter into sequence 2 alignment 
                seq1_align.insert(0, "-") #insert a gap character into sequence 1 alignment 
                seq2.pop() #remove the letter from sequence 2
            elif left == min_edit:   #it's a gap 
                seq2_align.insert(0, "-") #insert a gap character into sequence 2 alignment
                seq1_align.insert(0, seq1[-1]) #insert the letter into sequence 1 alignment 
                seq1.pop() #remove the letter from sequence 1
                y = y - 1 #move up one cell, using y because I transposed the matrix 
        else:  #once one sequence is down to zero characters
            if len(seq1) > len(seq2): #if there are still some letters left in sequence 1:
                x = 0 #make x = 0 to end the loop
                diff = len(seq1) - len(seq2) #find how many letters are left 
                while diff > 0: #add gap characters to alignment for the letters that are left unaligned 
                    seq2_align.insert(0, "-")
                    seq1_align.insert(0, seq1[-1])
                    diff = diff - 1
            elif len(seq2) > len(seq1): #same process, but if sequence 2 has leftover letters instead 
                x = 0
                diff = len(seq2) - len(seq1)
                while diff > 0:
                    seq1_align.insert(0, "-")
                    seq2_align.insert(0, seq2[-1])
                    diff = diff - 1

    return edit_dist, seq1_align, seq2_align #return the edit distance and the alignments 

distance_matrix = create_distance_matrix(seq1, seq2) #use function to create distance matrix

edit_dist, seq1_align, seq2_align = find_best_alignment(seq1, seq2, distance_matrix) #use function to create alignment 

with open(outfile, "w") as output_handle: #write to output file 
    output_handle.write(str(int(edit_dist)) + '\n') #convert the score to an integer, then a string for saving 
    output_handle.write(str(''.join(seq1_align)) + '\n') #convert the score to an integer, then a string for saving 
    output_handle.write(str(''.join(seq2_align)) + '\n') #convert the score to an integer, then a string for saving 
