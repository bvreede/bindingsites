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

infile = open("step1-random.tsv")

#collect the information from the tsv file into:
## - a dictionary identifying all possible locations of a 13-mer
dict13mer = {}
## - a list of all 13-mers, plus lists for the minus and plus strand
list13mer,pluslist, minuslist = [],[],[]
## - a dictionary collecting all hit locations in a given scaffold
dictscaf = {}
for line in infile:
	line = line.split() #line[0] is the scaffold, line[1] is the location, line[2] is the strand, and line[3] is the sequence
	list13mer.append(line[3]) #add the identified sequence to the list
	#add the sequence to the minus or plus strand
	if line[2] == '+':
		pluslist.append(line[3])
	else:
		minuslist.append(line[3])
	# sequence data
	if line[3] not in dict13mer: #if the sequence is unknown: make a new entry
		dict13mer[line[3]] = [[line[0],line[1],line[2]]]
	else: #the sequence is already in the dictionary, so append the location
		dict13mer[line[3]] += [[line[0],line[1],line[2]]]
	# scaffold data
	if line[0] not in dictscaf:
		dictscaf[line[0]] = [int(line[1])]
	else:
		dictscaf[line[0]] += [int(line[1])]

		
#sort the data in the scaffold dictionary so that the locations are in order
for scaffold in dictscaf:
	dictscaf[scaffold] = sorted(dictscaf[scaffold])
	
pluscommon,minuscommon,totalcommon = {},{},{}

#count how many times each 13 mer was found
count13mer = Counter(list13mer)
#identify the 15 most common sequences
print "25 most common sequences"
for k,mc in enumerate(count13mer.most_common(25)):
	print mc
	totalcommon[mc[0]]=[mc[1],k]
	
countplus = Counter(pluslist)
print "25 most common sequences on plus strand"
for k,mc in enumerate(countplus.most_common(25)):
	print mc
	pluscommon[mc[0]]=[mc[1],k]
	
countminus = Counter(minuslist)
print "25 most common sequences on minus strand"
for k,mc in enumerate(countminus.most_common(25)):
	print mc
	minuscommon[mc[0]]=[mc[1],k]
	
	
print "The most common sequences:"
for it in count13mer.most_common(25):
	s = it[0]
	print s, totalcommon[s]
	if s in pluscommon:
		print "on plus:", pluscommon[s]
	else:
		print "not in top-25 on plus strand"
	if s in minuscommon:
		print "on minus:", minuscommon[s]
	else:
		print "not in top-25 on minus strand"
	
"""	
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
	
#for mc in count13mer.most_common(25):
#	print mc[0], count13mer[mc[0]], clusterscores[mc[0]]
	
"""