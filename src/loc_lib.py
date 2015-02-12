
import pdb
import random
import doctest
import re
import threading
import time
import itertools


def hamdist(str1, str2):
   """Count the # of differences between equal length strings str1 and str2. If string length unequal returns length of the longer string"""
   diffs = 0
   if len(str1) != len(str2):
      return max(len(str1),len(str2))
   for ch1, ch2 in zip(str1, str2):
      if ch1 != ch2:
	  diffs += 1
   return diffs

def mutations(word, hamming_distance, charset='ATCG'):
    answer=set()
    for indices in itertools.combinations(range(len(word)), hamming_distance):
        for replacements in itertools.product(charset, repeat=hamming_distance):
                mutation = list(word)
                for index, replacement in zip(indices, replacements):
                    mutation[index] = replacement
                #yield "".join(mutation)
                answer.add("".join(mutation))
    return answer
