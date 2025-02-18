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

from Bio import SeqIO #import to look at phred scores and read in fastq

with open(infile) as f:
    num_min = int(f.readline().strip().split()[0]) #read first line only for minimum mean score allowed per position 
    num_low = 0 #start number of positions that do not meet the minimum at zero 

    all_scores = [] #initialize a list to keep all the scores in 

    for record in SeqIO.parse(f, "fastq"): #look at each record individually
        scores = record.letter_annotations["phred_quality"] #find phred scores 
        all_scores.append(scores) #add the scores from that record to our list of all scores 
    
    for i in range(len(all_scores[0])): #iterate down the length of a score starting with index i 
        scores_at_i = [] #initialize a list to keep all scores at the ith base position 
        for score_list in all_scores: #check each sequence individually 
            scores_at_i.append(score_list[i]) #add the score at the ith base position from each sequence to a list

        mean_at_i = sum(scores_at_i)/len(scores_at_i) #find the mean score at the ith base position 
        if mean_at_i < num_min: #if the mean is below our minimum 
            num_low += 1 #add one to our count

with open(outfile, "w") as output_handle: #write to output file 
    output_handle.write(str(num_low)) #save the number of base positions where mean score does not meet the minimum 
