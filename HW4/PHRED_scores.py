#INPUT: fastq file, where the first line indicates the minimum PHRED score and percentage of a sequence needed to pass that PHRED score to keep it 
#OUTPUT: the number of sequences in the fastq file that meet the quality constraints given by the first line

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

def check_scores(scores, nums): #function to check whether score meets conditions
    mini = int(nums[0]) #minimum score allowed
    perc = int(nums[1]) #percentage in sequence that must be at that score

    #find percentage above that quality score 
    count_above_mini = 0 #start count of good quality scores at 0 
    for score in scores: #check each score separately
        if score >= mini: #if at or above the minimum
            count_above_mini += 1 #add one to our count of good scores 
    
    if count_above_mini/len(scores)*100 >= perc: #if a high enough percent of scores pass the check
        checked = True #return true
    else: 
        checked = False #if the scores don't pass the check, return false

    return checked 


from Bio import SeqIO #import SeqIO for Phred scores 

with open(infile) as f:
    nums = f.readline().split() #read just the first line, split minimum score and percentage into a list
    num_good = 0 #start number of sequences that pass the test at 0 

    for record in SeqIO.parse(f, "fastq"): #check each record separately
        scores = record.letter_annotations["phred_quality"] #find phred scores
        checked = check_scores(scores, nums) #check whether that record meets the conditions for quality 
        if checked == True: #if it meets the conditions
            num_good += 1 #add to our count of good sequences 

with open(outfile, "w") as output_handle: #write to output file 
    output_handle.write(str(num_good)) #save the number of good sequences 
