#!/usr/local/python2.7

####### UNTESTED

from __future__ import print_function, division
import numpy as np 
import numpy.random as npr
import sys
import argparse

parser = argparse.ArgumentParser(description='''
	Run bootstrap resampling on percent reads masked. Prints output to stdout.

	USE: bootstrap.py <(grep 'COPIA' ANNOTATE.txt | awk '{print $1}') repeatmasker.out

	To get specific element repetitiveness:
		1. With RepeatMasker, first get AAARG ID's for all consensus 
			sequences in families annotated as COPIA.
		2. For each read sequence: bp masked/bp of read sequence = percent masked
		3. Bootstrap 10000 iterations of percent masked
		4. Calculate a 95 CI

	AUTHOR: Michelle Hwang
''')
parser.add_argument('fasta_headers', help='Name of file with FASTA headers.')
parser.add_argument('repeatmasker', help='Name of RepeatMasker outfile.')
#parser.add_argument('bp', type=int, help='Total bp of short reads')

args 		= parser.parse_args()
fasta_headers 	= args.fasta_headers
repeatmasker  	= args.repeatmasker
#total_bp		= args.bp

headers = open(fasta_headers, 'r')
RMfile  = open(repeatmasker, 'r')

#-----------
# Functions |
# ---------------------------------------------------------------------------------

def bootstrap(x, n, s, a):
	"""

	ARGS:
		x (list): percent masked for each read
		n (int): number of iterations
		s (np statistic): statistic to bootstrap
		a (float): alpha significance
	RETURNS:
		list: bootstrapped values
		int: lower CI
		int: upper CI
	"""
	l 		= len(x)
	rand 	= npr.randint(0, l, (n, l))
	samples = x[rand]
	stat 	= np.sort(s(samples, 1))
	return (stat, stat[int((a/2.0)*n)], stat[int((1-a/2.0)*n)])

def parse_headers(headers):
	"""
	ARGS:		
		headers (filehandle)
	RETURNS:
		list: list of ids
	"""
	ids = list()
	for header in headers:
		fields = header[1:].split('#')
		ids.append(fields[0])
	headers.close()
	return ids

def parse_RMfile(RMfile):
	"""
	ARGS:		
		RMfile (filehandle)
	RETURNS:
		list: list of lines in filehandle
	"""
	rms = list()
	for line in RMfile: #skip first 3 lines
		rms.append(line.rstrip())
	return rms

#------
# Main |
# ---------------------------------------------------------------------------------

def main():
	ids = parse_headers(headers) #slurp
	rms = parse_RMfile(RMfile) #slurp

	p = [] #list of % masked for all aaarf ids
	t = 0 #total bp masked

	for rm in rms[3:len(rms)-1]:
		col = rm.split()
		aaarf = col[9]
		left  = col[7]
		if aaarf in ids:
			if left == '(0)':
				p.append((int(col[6])-int(col[5])+1)/(int(col[6])+1))
			else:
				p.append((int(col[6])-int(col[5])+1)/(int(col[6])+1+int(left[1:-1])))
			t+=(int(col[6])-int(col[5])+1)
		else:
			p.append(0)


	# Bootstrap resample vector of % masked values 
	(samples, low, high) = bootstrap(np.array(p), 10000, np.mean, 0.05)

	print('BOOTSTRAP RESULTS')
	print('Statistic: ', np.mean(samples))
	print('95 CI:', low, high, sep=' ')
	print('Total BP Masked: ', t, end='\n\n')


main()












#if left == '(0)':
#p.append((int(col[6])-int(col[5])+1)/(int(col[6])+1))
#else:
#p.append((int(col[6])-int(col[5])+1)/(int(col[6])+1+int(left[1:-1])))