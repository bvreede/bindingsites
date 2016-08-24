'''
Description: Config file for dsx binding site search.
Author: Barbara Vreede
Contact: b.vreede@gmail.com
Date: 23 August 2016
'''

#Customize paths here
path_box = "/Users/BarbaraMaria/Box Sync/PROJECT doublesex"
path_genome = "Data/genome/Otaur.scaffolds.fa" #relative path to box folder


#SEQUENCES ASSOCIATE WITH DSX BINDING
#Erdman et al 1996 (G/A)NNAC(A/T)A(T/A)GTNN(C/T)
#Yi and Zarkower 1999 (nT)(nG)(T/A)ACAATGT(A/T)(nC)C
#Murphy et al 2007 (T/C)(G/A)(C/T)(A/T)AC(A/T)(T/A)(T/A)GT(A/T)(nC)
#Luo and Baker 2011 GCAACAATGTTGC


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
	re_ele = ['{','|','}']
	for i in s:
		if i not in re_ele:
			if flag == 0:
				reli.append(i)
			if flag == 1:
				i_ += i
		else:
			flag = 1
			if i == '{':
				i_ = i
			elif i == '}':
				i_ += i
				reli.append(i_)
				flag = 0
			else:
				i_ += i
return reli
	
	
def list_of_re(seq_in,n_mismatch):
	'''
	Turns a regular expression or DNA sequence into a list of regular expressions
	that can be used to search a scaffold. n_mismatch is the number of mismatches
	allowed in the sequence; the larger that number, the more sequences will result
	from this function, obviously.
	'''
	seq_out = ""
	seq_in_li = re2li(seq_in)
	# make a list of the locations that will be replaced by {A|T|C|G}
	replace_locations = [n for n in itertools.combinations(range(len(seq_in_li)),n_mismatch)]


		