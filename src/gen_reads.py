'''
Created on Feb 11, 2015

@author: govinda

syntax ./gen_reads id
'id' is an identifier used to distinguish files created due to this run of gen reads. All files created are given a suffix id.
inputs: none
outputs: Writes the following files:
underlying_seq : Contains the underlying sequences that were generated for this experiment:

'''
#!/usr/bin/env python

import sys
import random
import itertools
import profile
import time
import math
from loc_lib import hamdist,mutations
import operator
import numpy as np
from noise_add import add_sub_noise1
import pdb



if __name__ == '__main__':
    assert len(sys.argv) == 2
    ipid=sys.argv[1]
    
    ''' Parameters of reads generated'''
    m=100 #number of sequences
    L=100 #length of each sequence
    N=1000# No of reads
    k=25 # k-mer size
    error_prob=100 # noise rate is 1/error_prob
    ep=float(1)/error_prob #floating point representation of error probability

    
    '''Generate Sequences'''
    sequences=[] #sequences store the underlying sequences in the reads
    for i in xrange(m):
        seq=''
        for j in xrange(L):
            seq=seq+random.choice('AGCT')
        
        sequences.append(seq)

    '''Write sequences to underlying_seq#'''
    g=open('underlying_seq'+ipid,'w')
    for i in xrange(len(sequences)):
        g.write(str(sequences[i]) +'\n' )

    g.close()
    
    '''Generate Abundances'''
    p=np.random.lognormal(0,2,m) # p is the abundance of each sequence
    p=p/sum(p);
    print 'Sum of distribution is '+str(sum(p))+'.'
    g=open('prob_distribution'+ipid,'w')
    for i in xrange(len(p)):
        g.write(str(p[i])+'\n')
        
    g.close()
    
    '''Generate number of reads'''
    no_reads=np.random.multinomial(N,p) #no_reads is the number of reads of each of the sequences
    g=open('num_reads'+ipid,'w')
    for i in xrange(len(no_reads)):
        g.write(str(no_reads[i])+'\n')
        
    g.close()    
    
    '''Generate actual reads'''
    reads=[] #reads stores the actual reads
    for i in xrange(len(no_reads)):
        for j in xrange(no_reads[i]):
            temp=add_sub_noise1(sequences[i], error_prob)
            reads.append(temp)
            
    assert len(reads)==N
    
    g=open('reads_obtained'+ipid,'w')
    for i in xrange(len(reads)):
        g.write(str(reads[i])+'\n')
     
    g.close()   
    
    
    '''Generate true kmers'''
    true_kmers=set() #True kmers stores the actual kmers
    for i in xrange(len(sequences)):
        for j in xrange(L-k+1):
            true_kmers.add(sequences[i][j:j+k])
    
    assert len(true_kmers)==m*(L-k+1)
    
    g=open('true_kmers'+ipid,'w')
    for kmer in true_kmers:
        g.write(kmer+'\n')
        
    g.close()
    
    '''Generate observed kmer counts'''
    kmer_dict={} #kmers generated.
    
    for i in xrange(len(reads)):
        for j in xrange(L-k+1):
            kmer=reads[i][j:j+k]
            if kmer in kmer_dict.keys():
                kmer_dict[kmer]+=1
            else:
                kmer_dict[kmer]=1
                
    
    kmers_sorted=sorted(kmer_dict.items(), key=operator.itemgetter(1),reverse=True)
    g=open('kmer_dict'+ipid,'w')
    for kmer in kmers_sorted:
        g.write(kmer[0]+'\t'+str(kmer[1])+'\n')        
    g.close()
    
    
'''Partition  the  true kmers based on thresholds. Note that this assumes that each kmer only appears in  only 1 read.'''
    
    
