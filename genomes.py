# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 15:38:27 2024

@author: eblack
"""
import sys
sys.path.append("/Users/eblack/OneDrive - Imperial College London/Documents/Python training/Coursera bioinformatics for beginners/Genome sequencing (BioInf II)")

'''
Initially, we will assume a fictional ideal scenario: all reads from the same
strand, with no errors, with perfect coverage (every possible k-mer substring 
generated).
'''

# Find k-mer composition of a read

def StringComposition(Text, k):
    substrings = []
    n = len(Text)
    for i in range(n-k+1):
        substrings.append(Text[i:i+k])
    return substrings

'''
Genome assembly is solved by doing the inverse - constructing a string from its
constituent k-mers. This becomes difficult when k-mers repeat - it is not
possible to immediately determine which other k-mer "fits" in the next position
since more than one will. Need to "look ahead" to deal with this.
'''

#First simply reconstruct a whole string from all of the k-mers in order

def StringSpelledByPath(Path):
    genome = Path[0]
    n = len(Path[0])
    for i in range(1, len(Path)):
        genome = genome + Path[i][n-1]
    return genome

genome_path = open("genome_path.txt", 'r')
Path_raw = genome_path.read()
Path = Path_raw.split(' ')

String = StringSpelledByPath(Path)
genome_path.close()

print(String)