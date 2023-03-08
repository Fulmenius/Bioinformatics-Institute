# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 18:38:03 2022

@author: Fulmenius User

Three-level Manhattan global sequence alignment
"""

import numpy as np

sequence = input().split() # seq1 seq2 match mismatch gap


# Encode bases as numbers and vice versa
bases_to_numbers = {'_':0, 'A':1, 'T':2, 'G':3, 'C':4, 'U':5}
numbers_to_bases = {0:'_', 1:'A', 2:'T', 3:'G', 4:'C', 5:'U'}

match = float(sequence[2]) # match score
mu = float(sequence[3])    # mismatch penalty
delta = float(sequence[4]) # gap penalty
rho = delta                # gap start penalty
sigma = float(sequence[5]) # gap prolongation penalty


# Transform sequence into numeric form
def btn(s: str) -> np.array:
    return np.array(list(map(lambda x: bases_to_numbers[x], list(s))))

def ntb(s: np.array) -> np.array:
    return np.array(list(map(lambda x: numbers_to_bases[x], list(s))))

weight_matrix_standard = np.array([[-1, delta, delta, delta, delta, delta],   #    gap A T G C U
                                  [delta, match, mu, mu, mu, mu],             # gap
                                  [delta, mu, match, mu, mu, mu],             # A 
                                  [delta, mu, mu, match, mu, mu],             # T
                                  [delta, mu, mu, mu, match, mu],             # G
                                  [delta, mu, mu, mu, mu, match]])            # C
                                                                              # U

# Get optimal alignment from optimal path
def path_to_alignment(seq1: np.array, seq2: np.array, path: tuple) -> (str, str):
    seqA = list(seq1[::-1]) # We move from the beginning of the sequence to the end
    seqB = list(seq2[::-1]) # and .pop() gives us the last element, so we reverse the order
    path = path[::-1]       # path is generated from the end to the beginning, so we reverse it too
    aligned_seq1, aligned_seq2  = [], []
    
    def diff(i): # displacement vector at each step. Tells us what pair of symbols to choose
        return np.array(path[i+1]) - np.array(path[i])
    
    diff_path = list(map(diff, range(len(path)-1))) # rewrite the path as a sequence of displacements
    for i in range(len(diff_path)):                 # if (1, 1, x), append two nucleotides, 
                                                    # if (0, 1, x) or (1, 0, x), append gap+nucl or vice versa
        A = numbers_to_bases[seqA.pop()] if diff_path[i][0] == 1 else '_' 
        B = numbers_to_bases[seqB.pop()] if diff_path[i][1] == 1 else '_'
        aligned_seq1.append(A)
        aligned_seq2.append(B)
        
    return ''.join(aligned_seq1), ''.join(aligned_seq2)

def three_level_manhattan(seq1, seq2, weight_matrix):
    len1, len2 = len(seq1), len(seq2)
    
    paths = {} # For each element of the path scores matrix we will determine and remember 
               # the optimal precedent element
               # at the same time as we determine the value of the element itslef. paths{} is a
               # dictionary of links: path : (i, j) -> optimal_precedent((i, j))
                
    path_scores = np.zeros([len1+1, len2+1, 3]) #Path scores matrix initialization
    for i in range(3):
        path_scores[0, 0, i] = 0
    
    path_scores[1:, 0, 0] = [rho + sigma*i for i in range(1, len1+1)]
    path_scores[0, 1:, 0] = [-100 for j in range(1, len2+1)]
  
    path_scores[1:, 0, 1] = [-100 for i in range(1, len1+1)]
    path_scores[0, 1:, 1] = [-100 for j in range(1, len2+1)]

    path_scores[1:, 0, 2] = [-100 for i in range(1, len1+1)]
    path_scores[0, 1:, 2] = [rho + sigma*j for j in range(1, len2+1)]
    
#-----------------------------------------------------------------------------------------------------------
    
    for i in range(1, len1+1): # Optimal precedent elements on the boundary are determined uniquely
            paths[(i, 0, 0)] = (i-1, 0, 0)
            
    for j in range(1, len2+1):
            paths[(0, j, 2)] = (0, j-1, 2)
    
#-----------------------------------------------------------------------------------------------------------
# T_ij, layer 0
    
    for i in range(1, len1+1):
        for j in range(1, len2+1):
            
            prev_pos = ((i-1, j, 0), (i-1, j, 1), (i-1, j, 2)) # Possible precedent elements
            # Scores of paths that come from the set of optimal precedent elements
            prev_scores = (path_scores[i-1, j, 0] + sigma,\
                           path_scores[i-1, j, 1] + rho + sigma,\
                           path_scores[i-1, j, 2] + rho + sigma)
                                    
            path_scores[i, j, 0] = max(*prev_scores)
            # Add the link to the optimal precedent element into paths{}
            paths[(i, j, 0)] = max(zip(prev_pos, prev_scores), key = lambda x: x[1])[0]
                
#-----------------------------------------------------------------------------------------------------
# S_ij, layer 1

            prev_pos = ((i-1, j-1, 0), (i-1, j-1, 1), (i-1, j-1, 2)) # Possible precedent elements
            # Scores of paths that come from the set of optimal precedent elements
            prev_scores = (path_scores[i-1, j-1, 0] + weight_matrix[seq1[i-1], seq2[j-1]],\
                           path_scores[i-1, j-1, 1] + weight_matrix[seq1[i-1], seq2[j-1]],\
                           path_scores[i-1, j-1, 2] + weight_matrix[seq1[i-1], seq2[j-1]])
                                    
            path_scores[i, j, 1] = max(*prev_scores)
            # Add the link to the optimal precedent element into paths{}
            paths[(i, j, 1)] = max(zip(prev_pos, prev_scores), key = lambda x: x[1])[0]
                
#-----------------------------------------------------------------------------------------------------      
# U_ij, layer 2

            prev_pos = ((i, j-1, 0), (i, j-1, 1), (i, j-1, 2)) # Possible precedent elements
            # Scores of paths that come from the set of optimal precedent elements
            prev_scores = (path_scores[i, j-1, 0] + rho + sigma,\
                           path_scores[i, j-1, 1] + rho + sigma,\
                           path_scores[i, j-1, 2] + sigma)
                                    
            path_scores[i, j, 2] = max(*prev_scores)
            # Add the link to the optimal precedent element into paths{}
            paths[(i, j, 2)] = max(zip(prev_pos, prev_scores), key = lambda x: x[1])[0]
                
    optimal_path = [(len1, len2, max(zip([0, 1, 2], [path_scores[len1, len2, i] for i in range(3)]),\
                                    key = lambda x: x[1])[0])]
    
    # Optimal path must contain the element with the highest score
    optimal_score = max([path_scores[len1, len2, k] for k in range(3)])
    # The highest score in the table is the score of the alignment
    
    while optimal_path[-1] in paths: # Follow the links to the optimal elements backwards 
                                     # until we reach the element (0, 0)
        optimal_path.append(paths[optimal_path[-1]])
    
    return seq1, seq2, path_scores.transpose(), optimal_path, path_to_alignment(seq1, seq2, optimal_path), int(optimal_score)

