#!/user/bin/python/2.7

# /////////////////////////////////////////////////////////////////////////////
# Transpoable Element Annotation 
# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 
# Michelle Hwang
# Baucom Lab
# University of Michigan, Ann Arbor
# 2016
  

from __future__  import print_function
from __future__  import division
from collections import defaultdict as dd
from Bio 		 import SeqIO
import argparse
import operator
import re
import sys


# =============================================================================
# ARGUMENTS
# =============================================================================

parser = argparse.ArgumentParser(description='''
	This script generates an annotated database of transposable elements as a
	FASTA file. Assumes reads are already filtered by chloroplast, mitochondria,
	and host sequences.


	INPUT:
		(1) MCL outfile
		(2) BLAST outfile (outfmt=6 + stitle)
		(3) FASTA database of raw reads.

	OUTPUT:
		(1) FASTA file
		(2) Family info file
		(3) Report to stdout


	HOW TO RUN:
	annotate.py mcl.out blast.out reads.fa


	FOR HELP: 
	python annoate.py -h


	FASTA OUTPUT HEADER:
	>idfromAAARF#class#order#superfam#mcl''')

parser.add_argument('mcl', 
	help='Name of MCL outfile.')

parser.add_argument('blast',
	help='Name of BLAST outfile.')

parser.add_argument('fasta', 
	help='Name of database of raw reads in FASTA format.')

parser.add_argument('-n', '--ncbi',
	help='Name of NCBI BLAST outfile.')

parser.add_argument('-pre', '--prefix', default='result',
	help='''Output prefix. (DEFAULT=result)''')

parser.add_argument('-e', '--evalue', type=str, default=1, 
	help='''E-value threshold. (DEFAULT=1e1) 
	An input of "-2" will be interpreted as "1e-2". 
	All BLAST hits above this threshold will not be considered.''')

parser.add_argument('-l', '--length', type=int, default=0.6, 
	help='''Length threshold (bp). (DEFAULT=None). 
	Any BLAST hits shorter than this threshold 
	will not be considered.''')

parser.add_argument('-b', '--bit', type=float, default=0, 
	help='''Bit score threshold. (DEFAULT=None). 
	Any BLAST hit with a bit score lower than this 
	threshold will not be considered.''')

parser.add_argument('-p', '--pid', type=int, default=0, 
	help='''Percent identity threshold. (DEFAULT=None). 
	Any BLAST hit with a percent identity score 
	lower than this value will not be considered.''')

parser.add_argument('-t', '--thresh', type=float, default=0.2, 
	help='''(DEFAULT=0.2) This is interpreted as a percentage which is the
	percentage of elements in a TE family that must be annotated with the
	BLAST results for the family to be considered an annotated family and
	not an "NA" family.''')

parser.add_argument('-st', '--sthresh', type=float, default=0.5, 
	help='''(DEFAULT=0.5) This threshold is interpreted as a percentage which 
	is the percentage of elements in a *small* TE family that must be annotated 
	with the BLAST results for the family to be considered an annotated family 
	and not an "NA" family. The number of elements to be considered a small 
	family is determined by the "fthresh" argument.''')

parser.add_argument('-ft', '--fthresh', type=int, default=5, 
	help='''(DEFAULT=5) This is the threshold of members at which a 
	family is considered a "small family."''')

args = parser.parse_args()

## Get command line arguments
evalue 	= float(args.evalue)
length 	= args.length
bit 	= args.bit
pid 	= args.pid
thresh 	= args.thresh
sthresh = args.sthresh
fthresh = args.fthresh
PREFIX  = args.prefix


# =============================================================================
# GLOBAL VARIABLES
# =============================================================================

## Open up files
mcl_fh  	= open(args.mcl, 'r' )
mcl 		= mcl_fh.readlines() #Slurp
mcl_fh.close()

blast_fh	= open(args.blast, 'r')
BLAST 		= blast_fh.readlines() #Slurp
blast_fh.close()

NCBI = None
if args.ncbi is not None:
	ncbi_fh = open(args.ncbi, 'r')
	NCBI 	= ncbi_fh.readlines() #Slurp
	ncbi_fh.close()

myfasta		= open(args.fasta, 'r')
fastaiter 	= SeqIO.parse(myfasta, 'fasta')

myout 		= open('ANNOTATE_'+PREFIX+'.fasta', 'w')
fam_file 	= open('FAM_INFO_'+PREFIX+'.txt', 'w')

## Variables to hold TE family info
CLASSES  = ("1", "2", "HOST")
ORDERS1  = ("LTR", "SINE", "LINE", "RETRO")
ORDERS2  = ("TIR", "HELITRON")
SUPFAMS1 = ("GYPSY", "COPIA", "ANGELA", "ATHILA", "RETROLYC1")
SUPFAMS2 = ("ACDS", "ENSPM", "CACTA", "PIF/HARBRINGER", "HAT", "MITE", "HELITRON")


# =============================================================================
# FUNCTIONS
# =============================================================================

def account_for_ties(family, ties, top, i):
	"""If there are ties, then pick the one with the lower evalue.

	ARGS:
		family (list):  all elements within family
		ties (list): 	list of elements in family with top annotation
		top (string): 	top annotation
		i (int): 		class=1, order=2 or superfamily=3
	RETURNS:
		string: top annotation accounting for ties
	"""
	if ties != [] and len(ties) > 1:

		winning = 0
		winner 	= ''

		for tie in ties:
			evals = []
			for member in family:
				if member[i] == tie:
					evals.append(member[4])
			if get_mean(evals) < winning: 
				winning = get_mean(evals)
				winner 	= tie
		top = winner
		if winner is '':
			winner = 'NA'

	return top


def annotate_descr(s):
	"""Annotates TE based on BLAST description.

	ARGS:
		s (string): BLAST description 

	RETURNS:
		list: ID, class, order superfamily
	"""

	# Host:
	if(re.search('\s18s', s, re.I) != None 
		or re.search('\s5s', s, re.I) != None
		or re.search('ribosomal', s, re.I) != None
		or re.search('its1', s, re.I) != None
		or re.search('its2', s, re.I) != None 
		or re.search('internal', s, re.I) != None
		or re.search('\s45s', s, re.I) != None):
		return ['HOST', 'NA', 'NA']

	# Superfamily match:
	if(re.search('gypsy', s, re.I) != None
		or re.search( 'ty3', s, re.I) != None
		or re.search( 'cinful', s, re.I) != None
		or re.search( 'grande', s, re.I) != None
		or re.search( '\ssmilt', s, re.I) != None
		or re.search( '\sshuck', s, re. I) != None):
		return ['1', 'LTR', 'GYPSY']

	elif(re.search('copia', s, re.I) != None 
		or re.search('ty1', s, re.I) != None
		or re.search('\sji', s, re.I) != None
		or re.search('\sopie', s, re.I) != None
		or re.search('hopscotch', s, re.I) != None
		or re.search('\ssire', s, re.I) != None
		or re.search('\ssto', s, re.I) != None):
		return ['1', 'LTR', 'COPIA']

	elif(re.search('angela', s, re.I) != None):
		return ['1', 'LTR', 'ANGELA']

	elif(re.search('athila', s, re.I) != None):
		return ['1', 'LTR', 'ATHILA']

	elif(re.search('retrolyc1', s, re.I ) != None):
		return ['1', 'RETRO', 'RETROLYC1']

	elif( re.search( 'tip100', s, re.I ) != None
		or re.search( 'tip201', s, re.I ) != None
		or re.search( 'dihydroflavonol', s, re.I ) != None
		or re.search( 'activator', s, re.I ) != None
		or re.search( 'Ac-Ds', s, re.I ) != None
		or re.search( 'AcDs', s, re.I ) != None
		or re.search( 'Ac/Ds', s, re.I ) != None ):
		return ['2', 'TIR', 'ACDS']

	elif(re.search('tpn1', s, re.I ) != None
		or re.search('en/spm', s, re.I) != None
		or re.search('en-spm', s, re.I) != None
		or re.search('enspm', s, re.I) != None
		or re.search('\sspm', s, re.I) != None):
		return ['2', 'TIR', 'ENSPM']

	elif(re.search('cacta', s, re.I) != None):
		return ['2', 'TIR', 'CACTA']

	elif(re.search('tourist', s, re.I) != None):
		return ['2', 'TIR', 'PIF/HARBRINGER']
	elif(re.search( '\sshat', s, re.I) != None):
		return ['2', 'TIR', 'HAT']

	elif(re.search( 'mite', s, re.I) != None
		or re.search( 'mite-adh', s, re.I) != None):
		return ['2', 'TIR', 'MITE']

	elif(re.search('helitron', s, re.I) != None):
		return ['2', 'HELITRON', 'HELITRON']

	# Order match:
	if(re.search('endonuclease', s, re.I) != None
		or re.search( '\sline\s', s, re.I) != None ):
		return ['1', 'LINE', 'NA']
	elif(re.search( '\ssine\s', s, re.I) != None ):
		return ['1', 'SINE', 'NA']
	elif(re.search('\sgag\s', s, re.I) != None 
		or re.search('retroposon', s, re.I) != None
		or re.search('retrotransposon', s, re.I) != None
		or re.search('retrotransposon-like', s, re.I) != None
		or re.search('transcriptase-like', s, re.I) != None
		or re.search('polyprotein', s, re.I) != None 
		or re.search('reverse transcriptase', s, re.I) != None
		or re.search('\sbare', s, re.I) != None
		or re.search('sukkula', s, re.I) != None
		or re.search('leviathan', s, re.I) != None):
		return ['1', 'LTR', 'NA']
	elif(re.search('terminal repeat', s, re.I) != None
		or re.search('transposase', s, re.I) != None
		or re.search('tetn02', s, re.I) != None 
		or re.search(' mariner', s, re.I) != None ):
		return ['2', 'TIR', 'NA']

	return ['NA', 'NA', 'NA']

def classify(l):
	"""
	Classifies a transposable element by class, order, 
	and superfamily based on BLAST hits.

	ARGS:
		l (list): ID, fam_num

	RETURNS:
		list: ID, Class, Order, Supfam, Eval
	"""

	seqid 	= l[0]
	fam 	= l[1]

	eval1, length1, bit1, pid1, descr1 = 'NA', 'NA', 'NA', 'NA', 'NA'
	eval2, length2, bit2, pid2, descr2 = 'NA', 'NA', 'NA', 'NA', 'NA'

	for hits in BLAST:
		hit = hits.split('\t')
		if hit[0] == seqid:
			eval2 	= hit[10]	
			length2 = hit[3]
			bit2 	= hit[11]
			pid2 	= hit[2]
			descr2  = hit[12]

			if hit[10] != 'NA': 
				eval2 = float(hit[10])
			if hit[11] != 'NA': 
				bit2 = float(hit[11])						
			break

	if args.ncbi is None or descr2 is 'NA':
		if(eval2 == 'NA'):
			return [seqid, 'NA', 'NA', 'NA', 'NA']
		elif(bit2 < bit 
			or length2 < length 
			or eval2 > evalue 
			or pid2 < pid):
			return [seqid, 'NA', 'NA', 'NA', eval2 ]
		else:
			tedb = annotate_descr(descr2)
			tedb.insert(0, seqid)
			tedb.append(eval2)
			return tedb

	else:
		# 1 = NCBI
		# 2 = TE database
		for hits in NCBI:
			hit = hits.split('\t')
			if hit[0] == seqid:
				eval1	= hit[10]
				bit1 	= hit[11]
				length1 = hit[3]
				pid1 	= hit[2]
				descr1  = hit[12]

				if hit[10] != 'NA': 
					eval1 = float(hit[10])
				if hit[11] != 'NA': 
					bit1 = float(hit[11])

		## If no hits at all:
		if eval1 == 'NA' and eval2 == 'NA':
			return [seqid, 'NA', 'NA', 'NA', 'NA']

		## If BLAST hit only:
		elif eval1 == 'NA':
			if(bit2 < bit 
				or length2 < length 
				or eval2 > evalue 
				or pid2 < pid):
				return [seqid, 'NA', 'NA', 'NA', eval2]
			else:
				tedb = annotate_descr(descr2)
				tedb.insert(0, seqid)
				tedb.append(eval2)
				return tedb

		## If NCBI hit only:
		elif eval2 == 'NA':
			if(bit1 < bit 
				or length1 < length 
				or eval1 > evalue 
				or pid1 < pid):
				return [seqid, 'NA', 'NA', 'NA', eval1 ]
			else:
				ncbi = annotate_descr(descr1);
				ncbi.insert(0, seqid) # Put seqid at front of list
				ncbi.append(eval1) # Put BLAST e-value at end of list
				return ncbi

		## If hits for both, TE DB priority:
		else:
			db1 = annotate_descr(descr1) # NCBI
			db1.insert(0, seqid)
			db1.append(eval1)
			db2 = annotate_descr(descr2) # TE DB
			db2.insert(0, seqid)
			db2.append(eval2)

			## If NCBI scores sucks, use TE DB:
			if(bit1 < bit
				or length1 < length 
				or eval1 > evalue
				or eval1 > eval2
				or pid1 < pid):

				# But if TE DB score also sucks:
				if(bit2 < bit or length2 < length or eval2 > evalue or pid2 < pid):
					return [seqid, 'NA', 'NA', 'NA', 'NA'] 
				else:
					return db2 

			## If NCBI does not have superfam but TE does, use TE DB:
			elif db1[2] == 'NA' and db2[2] != 'NA':
				return db2 

			## If NCBI does not have order but TE does, use TE DB:
			elif db1[1] == 'NA' and db2[1] != 'NA':
				return db2 

			## If NCBI does not have class but TE does, use TE DB:
			elif db1[0] == 'NA' and db2[0] != 'NA':
				return db2 

			else:
				return db1 


def create_dict(l, i):
	"""
	ARGS:
		l (list): List of annotated TE lists within a family
		i (int):  1=class, 2=order, 3=superfamily
	RETURNS:
		dict: Dictionary of TEs and their annotation numbers.
	"""
	d = {}

	try:
		i = int(i)
	except ValueError:
		print('create_dict() function expects an integer')

	for x in l:
		if not d:
			d[x[i]] = 1
			next
		if d.has_key(x[i]) == True:
			d[x[i]] += 1
		else:
			d[x[i]] = 1
	return d


def get_ties(fam_info, i, top):
	"""Gets list of ties.

	ARGS:
		fam_info (list): list of dictionaries of all annotations and their counts
		i (int): 		 0=class, 1=order, 2=superfamily
		top (string): 	 top class/order/superfamily

	RETURNS:
		list: All elements within the list that match the top annotation
	"""
	ties = list()
	for element in fam_info[i]:
		# If both have same count
		if fam_info[i][element] == fam_info[i][top]:
			ties.append(element)
	return ties


def get_top(top, fam_info, i):
	""" Get top annotation. Will exclude superfamily or order values 
	not in determined class.

	ARGS:
		top (string): 		top annotation... class or order
		fam_info (list):	list of dictionaries of annotations and their counts
		i (int): 			order=1, superfamily=2

	RETURNS:
		string: top annotation
	"""

	if top == 'HOST':
		return 'NA'

	else:
		if top == '1' :
			for x in ORDERS2:
				if fam_info[i].has_key(x) is True: 
					del(fam_info[i][x])
			for y in SUPFAMS2:
				if fam_info[i+1].has_key(y) is True: 
					del(fam_info[i+1][y])
		elif top == '2' :
			for x in ORDERS1:
				if fam_info[i].has_key(x) is True: 
					del(fam_info[i][x])
			for y in SUPFAMS1:
				if fam_info[i+1].has_key(y) is True: 
					del(fam_info[i+1][y])

	return max(fam_info[i].iteritems(), key=operator.itemgetter(1))[0]


def get_mean(l):
	"""Returns mean given list of numbers.

	ARGS:
		l (list): list of numbers
	RETURNS:
		float: mean of numbers
	"""
	total = 0
	for n in l:
		if n == 'NA': continue 
		total = total + float(n)
	return(total/len(l))


def get_fasta(fasta_id):
	"""Given ID, get fasta sequence from fasta db.

	ARGS:
		fasta_id (string): ID of TE
	RETURNS:
		string: FASTA Sequence 
	"""
	myfasta.seek(0)
	fastaiter = SeqIO.parse(myfasta, 'fasta')
	return ''.join(str(seq.seq) for seq in fastaiter if seq.id in fasta_id)


def print_annot(l, remove):
	"""
	Prints TE annotations to stdout.

	ARGS:
		l (list): Class, order, superfamily, evalue 
		remove (bool): True or False if element is removed

	RETURNS:
		nothing
	"""
	if remove is False:
		print(l[1], l[2], l[3], l[4], sep='\t')
	else:
		print(l[1], l[2], l[3], l[4], '=REMOVED', sep='\t')


def print_family(l, t, i, host_tie):
	"""This function will print the results of the annotation for an element 
	in the annotation outfile.

	ARGS:
		l (list): current id, class, order, superfamily
		t (list): family top class, order, superfamily
		e (int): evalue
	RETURNS:
		nothing
	"""

	## When element is HOST:
	if l[1] == 'HOST' and host_tie is False:
		if t[0] == 'HOST':
			print_annot(l, False)
			print(">"+l[0]+"#"+l[1]+"#"+l[2]+"#"+l[3]+"#"+str(i)+"#"+PREFIX, 
				get_fasta(l[0]), sep="\n", file=myout)
		else:
			print_annot(l, True) ### WHEN DOES THIS HAPPEN??

	elif ((l[1] == t[0] and l[2] == t[1] and l[3] == t[2])
		 or (l[1] == t[0] and l[2] == t[1] and l[3] is 'NA')
		 or (l[1] == t[0] and l[2] == 'NA' and l[3] is 'NA')
		 or (l[1] == 'NA' and l[2] == 'NA' and l[3] == 'NA')):

		print_annot(l, False)
		print(">"+l[0]+"#"+t[0]+"#"+t[1]+"#"+t[2]+"#"+str(i)+"#"+PREFIX, 
			get_fasta(l[0]), sep="\n", file=myout)

	else:
		print_annot(l, True)


# =============================================================================
# MAIN
# =============================================================================

def main():
	i = 0 #Family index
	family = list()

	for myline in mcl:
		myline_strip = myline.strip() # Remove newline
		line = myline_strip.split('\t') 

		## If first family:
		if(i == 0):
			i = 1
			family.append(classify(line))

		## If not new family:
		elif(i == int(line[1])):
			family.append(classify(line))

		## If new family:
		else:
			n = len(family)

			## Create list of dictionaries of annotations - class, order, superfamily
			fam_info = [create_dict(family, 1), 
						create_dict(family, 2), 
						create_dict(family, 3)] 

			print(">FAM="+str(i), "N="+str(n), sep=", ")

			if fam_info[0].has_key('NA') is False:
				fam_info[0]['NA'] = 0		

			## If at least thresh annotated:
			if ((fam_info[0]['NA']/n <= 1 - thresh and n > fthresh) or 
				(fam_info[0]['NA']/n <= 1 - sthresh and n <= fthresh)):

				print("Passes annotation threshold.")
				fam_info[0]['NA'] = 0
				fam_info[1]['NA'] = 0
				fam_info[2]['NA'] = 0

				## Get top class
				top_class = max(fam_info[0].iteritems(), key=operator.itemgetter(1))[0]
				ties = get_ties(fam_info, 0, top_class) 

				host_tie = False
				if ('HOST' in ties and ('ENSPM' in ties or 'ACDS' in ties) and i==1):
					host_tie = True
					top = '2'
					continue

				top_class  = account_for_ties(family, ties, top_class, 1)

				## Get top order
				top_order  = get_top(top_class, fam_info, 1)
				top_order  = account_for_ties(family, get_ties(fam_info, 1, top_order), top_order, 2)

				## Get top superfamily
				top_supfam = get_top(top_order, fam_info, 2)
				top_supfam = account_for_ties(family, get_ties(fam_info, 2, top_supfam ), top_supfam, 3)

				## Resolves hierarhical issues
				if top_class == 'NA' or top_class == 'HOST':
					top_order, top_supfam = 'NA', 'NA'

				if top_order == 'NA':
					top_supfam = 'NA'

				for m in family:
					print_family(m, [top_class, top_order, top_supfam], i, host_tie)

			## If specified annotation threshold is not passed:
			else:
				print("Does not pass annotation threshold.")
				top_class, top_order, top_supfam = 'NA', 'NA', 'NA'

				for m in family:
					if m[1] == 'NA':
						print_annot(m, False)
						print(">"+m[0]+"#"+'NA'+"#"+'NA'+"#"+'NA'+"#"+str(i)+"#"+PREFIX, 
							get_fasta(m[0]), sep="\n", file=myout)
					else:
						print_annot(m, True)

			print(top_class, top_order, top_supfam, '(Top)', sep="\t")
			print( '------------------------------------------------------')

			## Create fam_file:
			print(i, top_class, top_order, top_supfam, sep='\t', file=fam_file)

			## Classify first member of the NEW family:
			family = []
			family.append(classify(line))

			## Update family index:
			i+=1

main()
myout.close()
fam_file.close()





