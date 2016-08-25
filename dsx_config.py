'''
Description: Config file for dsx binding site search.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 23 August 2016
'''
import itertools

#Customize paths here
path_box = "/Users/BarbaraMaria/Box Sync/PROJECT doublesex"
path_genome = "Data/genome/Otaur.scaffolds.fa" #relative path to box folder


# nucleotide translation: DNA opposite strand
opp = {'A':'T','T':'A','C':'G','G':'C'}

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
	
def re2li(s):
	'''
	Translate a string containing regular expressions to a list
	where each element of the re occupies a single item.
	'''
	reli = []
	flag = 0
	re_ele = ['[',']']
	for i in s:
		if i not in re_ele:
			if flag == 0:
				reli.append(i)
			if flag == 1:
				i_ += i
		else:
			flag = 1
			if i == '[':
				i_ = i
			elif i == ']':
				i_ += i
				reli.append(i_)
				flag = 0
			else:
				i_ += i
	return reli
	
def invert_re(expression):
	'''
	Inverting a regular expression so it can search the opposite strand.
	'''
	#turn expression into list 
	expli = re2li(expression)
	#reverse list
	expli.reverse()
	#list back to string, so expression is restored
	expression = ''.join(expli)
	inverted = ''
	#for each item in the expression, add the opposite nucleotide to the inverted expression
	for i in expression:
		if i in opp:
			inverted += opp[i] #add the nucleotide from the opposite strand to the inverted string
		else:
			inverted += i
	return inverted
				
	
def list_of_re(seq_in,n_mismatch):
	'''
	Turns a regular expression or DNA sequence into a list of regular expressions
	that can be used to search a scaffold. n_mismatch is the number of mismatches
	allowed in the sequence; the larger that number, the more sequences will result
	from this function, obviously.
	'''
	seq_in_li = re2li(seq_in)
	# make a list of the locations that will be replaced by [ATCG]
	replace_locations = [n for n in itertools.combinations(range(len(seq_in_li)),n_mismatch)]
	# for each replacement, make a regular expression and add to list
	relist = []
	for l in replace_locations:
		newre = ""
		for k,i in enumerate(seq_in_li):
			if k in l:
				newre += "[ATCG]"
			else:
				newre += i
		relist.append(newre)
	# turn list of re into set to remove duplicates
	relist = list(set(relist))
	invlist = []
	# obtain the inverse expression for each item in the list
	for r in relist:
		inv = invert_re(r)
		invlist.append(inv)
	return relist,invlist
	
	
def test_scaffold_proximity(loclist,location,maxdist):
	'''
	Checks whether a hit in a given location has hits closeby, from a list of other hits,
	within a given distance.
	'''
	minloc = int(location)-maxdist
	maxloc = int(location)+maxdist
	score = 0
	for hit in loclist:
		hit = int(hit)
		if hit == int(location): #don't count other hits on the exact same sequence
			continue
		elif hit >minloc and hit<maxloc:
			score+= 1
		elif hit > maxloc: #stop iterating once you've passed the max location
			break
	return score

		