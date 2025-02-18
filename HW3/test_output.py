with open("HW3_prob2_output", "r") as f:
    lines = [line.rstrip() for line in f]
    seq1_ex = lines[0]
    seq2_ex = lines[1]

with open("prob2_output_large", "r") as f:
    lines = [line.rstrip() for line in f]
    seq1_me = lines[0]
    seq2_me = lines[1]

print("ex1:", len(seq1_ex), "ex2:", len(seq2_ex))
print("me1:", len(seq1_me), "me2:", len(seq2_me))

from collections import Counter

def nucleotide_percentage(seq):
    seq = seq.upper()  # Ensure the sequence is uppercase
    counts = Counter(seq)  # Count occurrences of each nucleotide
    total = len(seq)  # Total length of the sequence
    
    percentages = {nuc: (counts[nuc] / total) * 100 for nuc in "ACTG"}
    return percentages

perc_ex1 = nucleotide_percentage(seq1_ex)
perc_ex2 = nucleotide_percentage(seq2_ex)
perc_me1 = nucleotide_percentage(seq1_me)
perc_me2 = nucleotide_percentage(seq2_me)

print("ex1:", perc_ex1, "ex2:", perc_ex2)
print("me1:", perc_me1, "me2:", perc_me2)

def is_cyclic_permutation(seq1, seq2):
    if len(seq1) != len(seq2):  # Sequences must be the same length
        return False
    return seq2 in (seq1 + seq1)

# Example usage:

print(is_cyclic_permutation(seq1_ex, seq1_me))  # Output: True
print(is_cyclic_permutation(seq2_ex, seq2_me))  # Output: True
