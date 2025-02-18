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

from Bio import SeqIO #helps take in the fasta file (SeqIO.parse())
from Bio.Seq import Seq #helps us convert the sequences to Seq() object when necessary 

#open the fasta file
with open(infile) as handle:
    all_seq = list(SeqIO.parse(handle, "fasta")) #separate out each sequence
    
    mrna = all_seq[0].seq #mark the first sequence as the mrna, save just the sequence instead of the whole record
    introns = all_seq[1:] #mark the second sequence and all following as the introns

    for intron in introns: #iterate through the introns
        if intron.seq in mrna: #if the intron sequence is contained in the mrna
            mrna = str(mrna) #convert the mrna to a string
            intron = str(intron.seq) #convert the intron to a string
            mrna = mrna.replace(intron, "") #replace the intron sequence with nothing - removes it
            mrna = Seq(mrna) #convert mrna back to a sequence 
    
    translated = mrna.translate(to_stop = True) #translate the mrna to get the protein - have it stop at the stop codon

    with open(outfile, "w") as output_handle: #write to output file 
        output_handle.write(str(translated)) #convert to string because our output is not a fasta 
