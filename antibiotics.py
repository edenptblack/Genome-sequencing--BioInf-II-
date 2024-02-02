# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 16:08:48 2024

@author: eblack
"""

# Create a genetic code dictionary

file_path = 'genetic_code.txt'

# Initialize an empty dictionary
genetic_code = {}

# Read the file and populate the dictionary
with open(file_path, 'r') as file:
    for line in file:
        # Split each line into codon and amino acid
        codon, amino_acid = line.strip().split(' : ')
        # Add the entry to the dictionary
        genetic_code[codon] = amino_acid

# A function to translate genetic code into an amino acid sequence

def ProteinTranslation(seq):
    amino_acids = []
    
    for i in range(0,len(seq),3):
        if genetic_code[seq[i:i+3]] == '*':
            break
        else:
            amino_acids.append(genetic_code[seq[i:i+3]])
    
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
    
    def PeptideSearch(text, peptide):
        n = len(peptide)
        
        # Need to account for likely error from lack of stop codon!