# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 16:08:48 2024

@author: eblack
"""

# Create a genetic code dictionary

file_path = 'genetic_code_rna.txt'

# Initialize an empty dictionary
genetic_code_rna = {}

# Read the file and populate the dictionary
with open(file_path, 'r') as file:
    for line in file:
        # Split each line into codon and amino acid
        codon, amino_acid = line.strip().split(' : ')
        # Add the entry to the dictionary
        genetic_code_rna[codon] = amino_acid

# A function to translate genetic code into an amino acid sequence

file_path = 'genetic_code_dna.txt'

# Initialize an empty dictionary
genetic_code_dna = {}

# Read the file and populate the dictionary
with open(file_path, 'r') as file:
    for line in file:
        # Split each line into codon and amino acid
        codon, amino_acid = line.strip().split(' : ')
        # Add the entry to the dictionary
        genetic_code_dna[codon] = amino_acid

def ProteinTranslation(seq):
    amino_acids = []
    
    for i in range(0,len(seq),3):
        if len(seq) < i+2:
            print("No stop codon")
            break
        elif genetic_code_rna[seq[i:i+3]] == '*':
            break
        else:
            amino_acids.append(genetic_code_rna[seq[i:i+3]])
    
    amino_acid_string = ''.join(amino_acids)
    
    return amino_acid_string

'''
A given genome fragment has 6 reading frames: 3 different starting positions
the sequence, x2 for being double-stranded
'''

# Functions to produce reverse complement

def Reverse(Pattern):
    rev = ""
    for char in Pattern:
        rev = char + rev
    return rev


def Complement(Pattern):
    comp_pattern = ""
    for char in Pattern:
        if char == "A":
            comp = "T"
        elif char == "T":
            comp = "A"
        elif char == "C":
            comp = "G"
        elif char == "G":
            comp = "C"
        comp_pattern = comp_pattern + comp
    return(comp_pattern)


def ReverseComplement(Pattern):
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return(Pattern)

# A function to find all instances of a peptide being encoded by a string of DNA

def PeptideEncoding(dna, peptide):
    rc_dna = ReverseComplement(dna)
    encoding_patterns = []
    
    def PeptideSearch(frame, peptide, reverse):
        n = len(peptide) * 3
        if reverse:
            frame = ReverseComplement(frame)
        
        for i in range(0, len(frame) - n + 1):
            peptide_check = frame[i:i+n]
            test_peptide = []
            
            for j in range(0, len(peptide_check) - 2, 3):
                test_peptide.append(genetic_code_dna[peptide_check[j:j+3]])
            
            test_peptide = ''.join(test_peptide)
            if test_peptide == peptide:
                if reverse:
                    encoding_patterns.append(ReverseComplement(frame[i:i+n]))
                else:
                    encoding_patterns.append(frame[i:i+n])

    
    for i in range(3):
        frame = dna[i:len(dna)]
        PeptideSearch(frame, peptide, reverse = False)
    
    for i in range(3):
        frame = rc_dna[i:len(rc_dna)]
        PeptideSearch(frame, peptide, reverse = True)
        
    encoding_patterns = list(set(encoding_patterns))

    return encoding_patterns



' '.join(map(str, PeptideEncoding('ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA', 'MA')))