#! /usr/bin/python

'''
Description: Second step of dsx binding site search: cataloging all identified motifs and
determining the most frequently found putative binding sites.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 24 August 2016
'''

from collections import Counter
import dsx_config

infile = open("step1.tsv")
out1 = open("step2_fasta.fa", "w")
out2 = open("step2_fasta-all.fa", "w")

#collect the information from the tsv file into:
## - a dictionary identifying all possible locations of a 13-mer
dict13mer = {}
## - a list of all 13-mers
list13mer = []
## - a dictionary collecting all hit locations in a given scaffold
dictscaf = {}
for line in infile:
	line = line.split() #line[0] is the scaffold, line[1] is the location, line[2] is the sequence
	list13mer.append(line[2]) #add the identified sequence to the list
	out2.write(">%s_%s\n%s\n\n" %(line[0],line[1],line[2]))
	# sequence data
	if line[2] not in dict13mer: #if the sequence is unknown: make a new entry
		dict13mer[line[2]] = [[line[0],line[1]]]
	else: #the sequence is already in the dictionary, so append the location
		dict13mer[line[2]] += [[line[0],line[1]]]
	# scaffold data
	if line[0] not in dictscaf:
		dictscaf[line[0]] = [int(line[1])]
	else:
		dictscaf[line[0]] += [int(line[1])]

		
#sort the data in the scaffold dictionary so that the locations are in order
for scaffold in dictscaf:
	dictscaf[scaffold] = sorted(dictscaf[scaffold])
	

#count how many times each 13 mer was found
count13mer = Counter(list13mer)
#identify the 15 most common sequences
for mc in count13mer.most_common(15):
	for i in dict13mer[mc[0]]:
		out1.write(">%s_%s\n%s\n\n" %(i[0],i[1],mc[0]))
		
# check for each 13mer whether they appear in clusters
clusterscores = {}
for mer in list13mer:
	hits = dict13mer[mer]
	totalscore = 0
	for h in hits:
		scaffold = h[0]
		loclist = dictscaf[scaffold]
		location = h[1]
		maxdist = 150 #the maximum distance between the hit and other hits
		score = dsx_config.test_scaffold_proximity(loclist,location,maxdist)
		totalscore += score
	clusterscores[mer] = totalscore
	
for mc in count13mer.most_common(25):
	print mc[0], count13mer[mc[0]], clusterscores[mc[0]]
	


# combine the two methods: select the 100 most common 13mers and check for clustering only between them