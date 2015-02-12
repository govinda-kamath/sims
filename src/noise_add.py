import pdb
import random


def add_sub_noise(read,noise_level): #level of noise =  one-in-noise_level.
    temp = ''
    k=0
    for i in xrange(len(read)-1):
	  j=random.choice(range(int(noise_level)))
	  if j==1:
	      #print 'added noise'
	      #k=k+1
	      #print k
	      if read[i]=='A':
		  temp=temp+random.choice('GCT')
	      elif read[i]=='G':
		  temp=temp+random.choice('ACT')
	      elif read[i]=='C':
		  temp=temp+random.choice('AGT')
	      elif read[i]=='T':
		  temp=temp+random.choice('ACG')
	  else:
	      temp=temp+read[i]
    temp=temp+'\n'
    return temp
    

def add_sub_noise1(read,noise_level): #level of noise =  one-in-noise_level. No newline at end of read here.
    temp = ''
    k=0
    for i in xrange(len(read)):
	  j=random.choice(range(int(noise_level)))
	  if j==1:
	      #print 'added noise'
	      #k=k+1
	      #print k
	      if read[i]=='A':
		  temp=temp+random.choice('GCT')
	      elif read[i]=='G':
		  temp=temp+random.choice('ACT')
	      elif read[i]=='C':
		  temp=temp+random.choice('AGT')
	      elif read[i]=='T':
		  temp=temp+random.choice('ACG')
	  else:
	      temp=temp+read[i]
    temp=temp
    return temp
