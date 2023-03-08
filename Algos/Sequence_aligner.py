# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 23:29:32 2022

@author: Fulmenius User
"""

from abc import ABC, abstractclassmethod
import numpy as np
import traceback
from numba import jit


class Sequence_aligner(ABC):
    
    @property
    @abstractclassmethod
    def algorithm(self):
        pass
    
    @abstractclassmethod
    def align(self, seq1, seq2, weight_matrix):
        pass
    

    