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
        if len(seq[i:i+3]) < 3:
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

def DnaProteinTranslation(seq):
    amino_acids = []
    
    for i in range(0,len(seq),3):
        if len(seq[i:i+3]) < 3:
            print("No stop codon")
            break
        elif genetic_code_dna[seq[i:i+3]] == '*':
            break
        else:
            amino_acids.append(genetic_code_dna[seq[i:i+3]])
    
    amino_acid_string = ''.join(amino_acids)
    
    return amino_acid_string


def PeptideEncoding(dna, peptide):
    encoding_patterns = []
    n = len(peptide) * 3
    
    for i in range(len(dna)):
        frame = dna[i:i+n]        
        test_peptide = DnaProteinTranslation(frame)
        
        if test_peptide == peptide:
            encoding_patterns.append(frame)
        
        rc_frame = ReverseComplement(dna[i:i+n])
        rc_test_peptide = DnaProteinTranslation(rc_frame)
        
        if rc_test_peptide == peptide:
            encoding_patterns.append(frame)

    return encoding_patterns


'''
Some peptides are cyclic, meaning we need to search the genome for the linear
version of every possible starting point on the circle.

In bacillus, tyrocidines & gramicidins are non-ribosomal peptides - rather
than being DNA encoded they are assembled amino acid by acid by NRP synthetase.
Therefore we can't infer them from the genome.

Mass spectrometry allows us to identify peptides via their molecular weight in
Daltons, where 1Da is approximately the weight of one proton/neutron. Each
amino acid has a known weight, and we therefore know the weight of peptides.

Mass spectrometry works by breaking peptides into fragments, with each copy
breaking it's own way. We need to be able to identify peptides from the
resulting experimental spectrum.

To start, we will assume that a cyclic peptide will be split at every possible
2 bonds, so the spectrum contains all possible subpeptides. The theoretical
spectrum of a peptide is an ordered list of the masses of all possible
subpeptides, starting at 0 and ending at the mass of the whole peptide.
'''

# Dictionary of amino acid masses

file_path = 'integer_mass_table.txt'

mass_table = {}

with open(file_path, 'r') as file:
    for line in file:

        amino_acid, mass = line.strip().split(' ')
        # Add the entry to the dictionary
        mass_table[amino_acid] = mass

print(mass_table['N'])

# Generating the theoretical spectrum of a linear peptide

def LinearSpectrum(peptide):
    PrefixMass = [0]
    
    for i in range(0, len(peptide)):
        PrefixMass.append(PrefixMass[i] + int(mass_table[peptide[i]]))
    
    linear_spectrum = [0]
    
    for i in range(len(peptide) - 1):
        for j in range(i+1:len(peptide)):
            linear_spectrum.append(PrefixMass[j] - PrefixMass[i])
    
    linear_spectrum.sort()      
    
    return linear_spectrum


peptide = 'NQEL'    

LinearSpectrum(peptide)

        