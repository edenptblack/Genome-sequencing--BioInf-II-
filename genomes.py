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

'''
To apply this to an unordered collection of reads, we need to use graph theory.
Graph nodes and directional edges can be represented as a dictionary of lists
with the entry for each node containing a list of all out-connected nodes
'''

# Directional graph connecting overlapping k-1 suffixes to prefixes

def OverlapGraph(Dna):
    graph = {}
    
    for i in Dna:
        if i in graph:
            pass
        else:
            for j in range(len(Dna)):
                if i[1:] == Dna[j][:-1]:
                    if i not in graph:
                        graph[i] = []
                    graph[i].append(Dna[j])
    
    graph = dict(sorted(graph.items()))

    return graph

'''
The k-mers are now edges, with nodes representing the overlaps

To reconstruct a string, we need to find a Hamiltonian path, ie a path which
goes through every node exactly 1 time
'''

# Create a de Bruijn graph by gluing together identical overlap nodes (adjacency list form)

def deBruijnAdjacency(Text, k):
    nodes = StringComposition(Text, k-1)
    graph = {}
    
    for i in range(len(nodes)):
        if i + 1 >= len(nodes):  # Adjust the condition to avoid index out of range
            break
        elif nodes[i] not in graph:
            graph[nodes[i]] = []
            graph[nodes[i]].append(nodes[i + 1])
        elif nodes[i] in graph:
            graph[nodes[i]].append(nodes[i + 1])
    
    graph = dict(sorted(graph.items()))
    
    return graph

'''
Reconstructing the string requires finding a path through the de Bruijn graph
which touches each edge exactly once. Called the Eulerian path.

For this to be a viable solution, we need to be able to construct a de Bruijn
graph without needing the full genome. Can do this from a composition graph.
'''

# Create a de Bruijn graph from a set of k-mers

def deBruijn(Patterns):
    graph = {}
    for i in range(len(Patterns)):
        if Patterns[i][:-1] not in graph:
            graph[Patterns[i][:-1]] = [Patterns[i][1:]]
        else:
            graph[Patterns[i][:-1]].append(Patterns[i][1:])
    
    return graph

'''
The Hamiltonian method and the de Bruijn graph are different approaches to the 
same problem: Every k-mer node once, or every k-mer edge once. With 2 different
approaches, the one which provides the most efficient algorithm should be used

de Bruijn graph approach can be solved with Eulerian algorithms.


For a directed graph to be Eulerian:
    It must be balanced (each node has the same indegree & outdegree)
    It must be strongly connected (all nodes must be reachable from each other)

Euler's Theorem states that any graphs meeting these criteria are Eulerian.

It is impossible to get stuck in a cycle except at the starting node. The only 
way to get stuck is to return to the start before visiting every node.

Therefore: Complete a cycle. If it is not Eulerian, choose a node in cycle1 with
unused edges entering & leaving - complete cycle 1 and travel on unused edges
for cycle 2. If still not Eulerian, iterate until full cycle found.
'''

# Find the Eulerian cycle of a graph from a dictionary of nodes & edges

import random


def EulerianCycleNonLinear(graph):
    cycle = []
    path_travelled = {key:[] for key in graph}
    
    node = random.choice(list(graph.keys()))
    cycle.append(node)
    
    while not all(item in cycle for item in graph.keys()):
        print("New loop")
        
        if len(path_travelled[node]) == 0:
            next_steps = graph[node]
        else:
            next_steps = [item for item in graph[node] if item not in path_travelled[node]]
        
        
        if len(next_steps) > 0:
            next_node = random.choice(next_steps)
            path_travelled[node].append(next_node)
            node = next_node
        
        else:
            new_starts = []
            
            for item in cycle:
                if item != node:
                    if not all(item in graph[node] for item in cycle):
                        new_starts.append(item)
            print("NEW CYCLE")
            cycle = []
            path_travelled = {key:[] for key in graph}
            node = random.choice(new_starts)
         
        cycle.append(node)    
         
    return cycle

def EulerianCycle(graph):
    def dfs(node):
        while graph[node]:
            next_node = graph[node].pop()
            dfs(next_node)
        cycle.append(node)

    cycle = []
    stack = [random.choice(list(graph.keys()))]

    while stack:
        current_node = stack[-1]
        if graph[current_node]:
            next_node = graph[current_node].pop()
            stack.append(next_node)
        else:
            cycle.append(stack.pop())

    # Reverse the cycle to get the correct order
    cycle.reverse()

    return cycle
        
