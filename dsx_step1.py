#! /usr/bin/python
'''
Description: First step of dsx binding site search: identifying all putative binding site
motifs in a genome.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 23 August 2016
'''


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

dsxbinding,dsxinverse = [],[]
for seq in erdman,yizarkower,murphy,luobaker:
	poss_sequences,poss_inverses = dsx_config.list_of_re(seq,n_mismatch) #generate all possible sequences with mismatch
	dsxbinding += poss_sequences
	dsxinverse += poss_inverses
	
dsxbinding_compiled = [re.compile(r) for r in dsxbinding] #compile the regular expressions
dsxinverse_compiled = [re.compile(r) for r in dsxinverse]

# Search the genome for all possible sequences
# Write the results to an output file
out = open("step1.tsv", "w")
for scaffold in genomedict:
	duplicate_prevention = [] #empty the duplicate prevention for each scaffold
	direction = '+' #first search the plus strand
	for motif in dsxbinding_compiled:
		for i in motif.finditer(genomedict[scaffold]):
			mseq = i.group() # the sequence picked up by the RE
			strt = i.start() # the start site of the identified motif
			if strt not in duplicate_prevention:
				out.write("%s\t%s\t%s\t%s\n" %(scaffold,strt,direction,mseq)) #write results to file
				duplicate_prevention.append(strt)
	duplicate_prevention = [] #empty the duplicate prevention before starting on the inverse!
	direction = '-' #now search the minus strand
	for motif in dsxinverse_compiled: #now search for the inverse motifs! 
		for i in motif.finditer(genomedict[scaffold]):
			mseq = dsx_config.invert_re(i.group()) # the sequence picked up by the RE, inversed!
			strt = i.end() # the END site of the identified motif, which is the START site from the other side!
			if strt not in duplicate_prevention:
				out.write("%s\t%s\t%s\t%s\n" %(scaffold,strt,direction,mseq)) #write results to file
				duplicate_prevention.append(strt)