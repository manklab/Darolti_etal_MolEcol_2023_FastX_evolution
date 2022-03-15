#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
''' Identify sites passing minimum coverage threshold
Filters an mpileup file to identify scaffolds with at least one site that passes coverage threshold imposed in 05.filter-snps.py.'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import time
import cPickle as pickle
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
					help="An mpileup file")
parser.add_argument("coverage", type=str,
					help="Minimum coverage (DP >= 20) in order to call a SNP at a given site")
parser.add_argument("coverage_threshold", type=str,
					help="Minimum no. of individuals with DP >= 20 in order to call a SNP at a given site")
parser.add_argument("gene_position", type=str,
					help="A file containing positional information on each gene ie which scaffold it is located on"
					"Format: Geneid,Scaffold,Chromosome,Startwithinscaffold,Startwithingenome")
parser.add_argument("outfile", type=str,
					help="An outfile containing a list of valid scaffolds in the mpileup file and the sites that pass the coverage threshold")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

def scaffold_gene_dictionary(gene_position):
	gene_scaffold_dict = {}
	scaffold_gene_dict = defaultdict(list)
	with open(gene_position, "r") as gene_infile:
		# Get dictionary of gene and which scaffold it is located on
		for line in gene_infile:
			if not line.startswith("Gene"):
				line = line.rstrip()
				gene = line.split(",")[0]
				scaffold = line.split(",")[1]
				gene_scaffold_dict[gene] = scaffold
				scaffold_gene_dict[scaffold].append(gene)
		print "Number of genes with positional information =", len(gene_scaffold_dict)
	return scaffold_gene_dict

def parse_mpileup(source, scaffold_gene_dict, coverage, coverage_threshold):
	#'mpileup format: scaffold site ref SDP1 readbases1 quality1 SDP2 readbases2 quality2...'
	#site refers to position within scaffold ie 1 = 1st bp in scaffold
	count = 0
	count_lines = 0
	scaffold_names = {}
	valid_sites_in_scaffold = defaultdict(list)
	with open(source, "r") as infile:
		for line in infile:
			sys.stdout.write('%d\r' % (count_lines))
			sys.stdout.flush()
			count += 1
			count_lines += 1 
			sample_SDP = 0
			line = line.rstrip()
			line = line.split("\t")
			# Count number of scaffolds
			scaffold = "scaffold"+line[0]
			if scaffold in scaffold_gene_dict:
				scaffold_names[scaffold] = 0
				site = line[1]
				# Extract mapping info for each individual
				sample_info = line[3:]
				sample_number = int(float(len(sample_info))/3)
				# Check last value in sample info is the same as the last value in line
				if line[-1] != sample_info[-1]:
					print "ERROR - info for each sample has been selected incorrectly"
				# range stop: Generate numbers up to, but not including this number.
				# not need to add one to range because want in python counting as these are positions in list
				sample_SDP_position = [3*i for i in range(0,sample_number)]
				# Count no of individuals with SDP >= 20
				for i in sample_SDP_position:
					SDP = float(sample_info[i])
					if SDP >= float(args.coverage):
						sample_SDP += 1
				# Store valid sites
				if sample_SDP >= float(args.coverage_threshold): 
					valid_sites_in_scaffold[scaffold].append(site)

	print "Number of lines in mpileup =", count
	print "Number of scaffolds in mpileup =", len(scaffold_names)
	print "Number of scaffolds with at least one site that passes coverage threshold =", len(valid_sites_in_scaffold)
	return valid_sites_in_scaffold

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	startTime = time.time()

	# Get dictionary of genes with positional information
	scaffold_gene_dict = scaffold_gene_dictionary(args.gene_position)

	# Identify sites in each scaffold passing threshold
	valid_sites_in_scaffold = parse_mpileup(args.infile, scaffold_gene_dict, args.coverage, args.coverage_threshold)

	# save as a pickle file
	print "Starting to pickle..."
	with open(args.outfile, 'w') as outfile:
		pickle.dump(valid_sites_in_scaffold, outfile, -1)

	print ('The script took {0} seconds !'.format(time.time() - startTime))

if __name__ == '__main__':
	main()
