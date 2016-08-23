'''
Description: Config file for dsx binding site search.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 23 August 2016
'''

#Customize paths here
path_box = "/Users/BarbaraMaria/Box Sync/PROJECT doublesex"
path_genome = "Data/genome/Otaur.scaffolds.fa" #relative path to box folder



def open_genome(fastadb):
	'''
	Function to open large sequence file (in fasta format) and read into
	memory; make dictionary.
	'''
	fastadict = {}
	sequence,header = "",""
	for line in fastadb:
		if line[0] == ">":
			if header != "":
				fastadict[header] = sequence
			header = line[1:].strip()
			sequence = ""
		else:
			sequence += line.strip()
	fastadict[header] = sequence
	return fastadict