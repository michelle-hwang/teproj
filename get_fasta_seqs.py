#!/usr/local/python/2.7.3/bin/python


from Bio import SeqIO
import argparse
from numpy import random as rand
import sys


# -----------------------------------------------------------------------------


parser = argparse.ArgumentParser(description='''
	This script parses a multi-line FASTA file and extracts sequences based
	on instructions.

	AUTHOR: Michelle Hwang''')

parser.add_argument('fasta', help='Name of FASTA db.')

parser.add_argument('-k', '--keep', help='''List of names, 1 per line, 
	of sequences, to keep.''')
parser.add_argument('-i', '--invert', action="store_true", help='''If 
	specified, remove list of sequences instead of grabbing them.''')

parser.add_argument('-d', '--divide', help='''Instead of grabbing/filtering 
	sequences, divide FASTA file into separate files of x sequences each.''')
parser.add_argument('-p', '--prefix', help='''Out prefix for divide flag.''')

parser.add_argument('-r', '--random', type=int,
	help='''Randomly sample this number of sequences from the FASTA file.''')

parser.add_argument('-a', '--append', type=str, 
	help='''Append string to all headers.''')
args = parser.parse_args()


seqiter = SeqIO.parse(open(args.fasta), 'fasta')


# -----------------------------------------------------------------------------

def append(seqiter):
	for seq in seqiter:
		seq.id = seq.id + str(args.append)
		SeqIO.write(seq, out, 'fasta') 		


def divide():
	i 	  = 1
	label 	  = 10
	seqs 	  = list()
	seqiter   = list(seqiter)
	l 	  = len(seqiter)
	d	  = int(args.divide)

	for seq in seqiter:
		seqs.append(seq)

		if ((i == d) or (i % d == 0) or (i == l)):
			out   = open(args.prefix+'-'+str(label)+'.fa', 'w')
			SeqIO.write((s for s in seqs), out, 'fasta') 		
			out.close()
			label = label + 10
			seqs  = list()

		i = i+1 


def random(seqiter):
	n 	   = int(args.random)
	seqs   = list(seqiter)
	SeqIO.write((seq for seq in random.sample(seqs, n)), 
		sys.stdout, 'fasta')	


def get(seqiter):
	wanted 	= [line.strip() for line in open(args.keep)]

	if args.invert:
		SeqIO.write((seq for seq in seqiter if seq.id.split(' ')[0] 
			not in wanted), sys.stdout, 'fasta')
	else:
		#SeqIO.write((seq for seq in seqiter if seq.id.split(' ')[0] in wanted), sys.stdout, 'fasta')
		for seq in seqiter:
			if seq.id.split(' ')[0] in wanted or seq.id.split('#')[0] in wanted:
				SeqIO.write(seq, sys.stdout, 'fasta')

# -----------------------------------------------------------------------------


def main():
	seqiter = SeqIO.parse(open(args.fasta), 'fasta')

	## Divide sequences
	if args.divide is not None:
		divide(seqiter)

	## Randomly sample sequences
	elif args.random is not None:
		random(seqiter)

	elif args.append is not None:
		append(seqiter)

	## Grab or filter sequences
	else:
		get(seqiter)


main()

