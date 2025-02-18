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

with open(infile) as f: #open the file
    contigs = f.read().splitlines() #make a list of all the contigs
    
    total_length = 0 
    for contig in contigs: 
        total_length += len(contig) #add all contigs together to get the total length
    
    contigs.sort(key = len, reverse = True) #sort list of contigs, longest to shortest

    length_50 = total_length*0.5 #find 50% of length
    length_75 = total_length*0.75 #find 75% of length 

    acc_length = 0 #start accumulated length at zero 
    i = -1 #start index so we can iterate through the contigs
    while acc_length < length_50: #as long as the contigs we've added together don't reach 50% of the total length 
        i = i + 1 #move to the next contig in our list
        acc_length += len(contigs[i]) #add the length of that contig to our accumulation 
    
    n50 = len(contigs[i]) #the n50 is the shortest contig that had to be added to reach length_50 

    while acc_length < length_75: #same process, but until we reach 75% of the total length 
        i = i + 1 
        acc_length += len(contigs[i])
    
    n75 = len(contigs[i]) #the n75 is the shortest contig that had to be added to reach length_75

with open(outfile, "w") as output_handle: #write to output file 
    output_handle.write(str(n50) + " " + str(n75)) #write n50 and n75 to output
