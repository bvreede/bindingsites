#! /usr/bin/python

import dsx_config, re


# Import the genome
genome = open("%s/%s" %(dsx_config.path_box,dsx_config.path_genome))
genomedict = dsx_config.open_genome(genome) #read genome file and save as dictionary

# Generate the sequences that will be searched
erdman = "[AG][ATCG][ATCG]AC[AT]A[AT]GT[ATCG][ATCG][CT]"
yizarkower = "[ACG][ATC][AT]ACAATGT[AT][ATG]C"
murphy =  "[TC][AG][CT][AT]AC[AT][AT][AT]GT[AT][ATG]"
luobaker = "GCAACAATGTTGC"

n_mismatch = 2 #number of elements per sequence that can be a mismatch

dsxbinding = []
for seq in erdman,yizarkower,murphy,luobaker:
	dsxbinding += dsx_config.list_of_re(seq,n_mismatch) #generate all possible sequences with mismatch
	
dsxbinding = list(set(dsxbinding)) #remove duplicates
dsxbinding_compiled = [re.compile(r) for r in dsxbinding] #compile the regular expressions

