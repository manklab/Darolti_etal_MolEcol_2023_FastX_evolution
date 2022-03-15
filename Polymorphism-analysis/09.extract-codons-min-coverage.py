#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
''' Identify codons where all sites pass minimum coverage threshold imposed in 05.filter-snps.py. 
Reports codons within ORF. Position reported in relation to codons within ORF. i.e. 1 = 1st codon in ORF/protein. Reverse complements are accounted for.'''
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
					help="An mpileup pickle file")
parser.add_argument("gene_position", type=str,
					help="A file containing positional information on each gene ie which scaffold it is located on"
					"Format: Geneid,Scaffold,Chromosome,Startwithinscaffold,Startwithingenome")
parser.add_argument("ORF_sites_in_scaffold", type=str,
					help="A file containing genes and the position of sites within the ORF in the scaffold"
					"Output of 07.orf-position-within-scaffold.py")
parser.add_argument("outfile", type=str,
					help="An outfile file of codons passing coverage threshold")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

def read_pkl(pkl_file):
	"Reads in a pickled file and returns the stored object"
	with open(pkl_file, "rb") as infile:
		return pickle.load(infile)

def get_dictionary_which_gene_scaffold(gene_position):
	gene_scaffold_dict = {}
	with open(gene_position, "r") as gene_infile:
		# Get dictionary of gene and which scaffold it is located on
		next(gene_infile) # skips header
		for line in gene_infile:
			line = line.rstrip()
			gene = line.split(",")[0]
			scaffold = line.split(",")[1]
			gene_scaffold_dict[gene] = scaffold
		print "Number of genes with positional information =", len(gene_scaffold_dict)
	return gene_scaffold_dict

def get_ORF_sites_in_scaffold(ORF_sites_in_scaffold):
	ORF_sites_in_scaffold_dict = {}
	with open(ORF_sites_in_scaffold, "r") as infile:
		for line in infile:
			line = line.rstrip().split("\t")
			gene = line[0].split(".")[0]+"."+line[0].split(".")[1]
			sites = line[1].split("]")[0].split("[")[1].split(", ")
			ORF_sites_in_scaffold_dict[gene] = sites
	print "Number of genes with ORF =", len(ORF_sites_in_scaffold_dict)
	return ORF_sites_in_scaffold_dict

def pass_sites_in_ORF(valid_sites_in_scaffold, gene_scaffold_dict, ORF_sites_in_scaffold):
	#site still refers to position within scaffold ie 1 = 1st bp in scaffold
	valid_sites_in_gene = {}
	for gene in gene_scaffold_dict:
		scaffold = gene_scaffold_dict[gene]
		valid_sites = valid_sites_in_scaffold[scaffold]
		if len(valid_sites) != 0:
			if gene in ORF_sites_in_scaffold:
				ORF_valid_sites = list(set(valid_sites).intersection(ORF_sites_in_scaffold[gene]))
				if len(ORF_valid_sites) > 0:
					valid_sites_in_gene[gene] = ORF_valid_sites
	print "Number of genes with ORF with at least one site within the ORF that passes the coverage threshold =", len(valid_sites_in_gene)
	return valid_sites_in_gene

def pass_codons_in_ORF(valid_sites_in_gene, gene_scaffold_dict, ORF_sites_in_scaffold_dict):
	#codon refers to position within ORF ie 1 = 1st codon in ORF/protein.
	valid_codons_in_ORF = defaultdict(list)
	for gene in gene_scaffold_dict:
		if gene in ORF_sites_in_scaffold_dict:
			if gene in valid_sites_in_gene:
				sites = ORF_sites_in_scaffold_dict[gene]
				# Get dictionary of codon position (within ORF) and nucleotide position (within scaffold)
				position_dictionary = identify_codons(sites)
				# Check all sites within each codon pass minimum coverage
				valid_sites = valid_sites_in_gene[gene]
				for codon in position_dictionary:
					triplet = position_dictionary[codon]
					if len(list(set(triplet).intersection(valid_sites))) == 3:
						valid_codons_in_ORF[gene].append(codon)
	print "Number of ORF with at least one codon where all sites pass coverage threshold =", len(valid_codons_in_ORF)
	return valid_codons_in_ORF

def identify_codons(sites):
	ORF_nucleotides = sites
	# Check length of ORF is a multiple of three
	if len(ORF_nucleotides) % 3:
		print "ERROR - ORF nucleotides not selected correctly"
	# Get list of codon positions - human counting
	ORF_codons = range(1,(len(ORF_nucleotides)/3)+1)
	# Check there are a third the number of codons than nucleotides in ORF
	if len(ORF_nucleotides)/3 != len(ORF_codons):
		print "ERROR - codon position not correct"
	# No need to reverse the nucleotides if reverse complement because already reversed in 07.orf-position-within-scaffold.py
	# Split nucleotides into triplets
	triplets = [ORF_nucleotides[i:i +3] for i in range(0, len(ORF_nucleotides), 3)]
	# Check there are the same number of triplets and codons
	if len(triplets) != len(ORF_codons):
		print "ERROR - codon position not correct"
	# Make dictionary between codon position and triplet position
	position_dictionary = {}
	for codon in ORF_codons:
		triplet = triplets[codon-1] # Codon is in human counting. Need to substract 1 as looking up position in triplets.
		position_dictionary[codon] = triplet # position_dictionary is in human counting
	return position_dictionary

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	startTime = time.time()

	# Read mpileup pickle infile (Extract sites in each scaffold passing threshold)
	valid_sites_in_scaffold = read_pkl(args.infile)
	print "Number of scaffolds with at least one site that passes coverage threshold =", len(valid_sites_in_scaffold)

	# Get dictionary of coordinates of ORF within scaffold
	gene_scaffold_dict = get_dictionary_which_gene_scaffold(args.gene_position)

	# Make dictionary of valid genes and the position in scaffold of all sites within ORF
	ORF_sites_in_scaffold_dict = get_ORF_sites_in_scaffold(args.ORF_sites_in_scaffold)

	# Identify sites in each ORF passing threshold
	valid_sites_in_gene = pass_sites_in_ORF(valid_sites_in_scaffold, gene_scaffold_dict, ORF_sites_in_scaffold_dict)

	# Identify codons in each ORF where all sites pass threshold
	valid_codons_in_ORF = pass_codons_in_ORF(valid_sites_in_gene, gene_scaffold_dict, ORF_sites_in_scaffold_dict)

	with open(args.outfile,"w") as outfile:
		for gene in valid_codons_in_ORF:
			outfile.write(gene)
			outfile.write("\t")
			outfile.write(str(valid_codons_in_ORF[gene]))
			outfile.write("\n")

	print ('The script took {0} seconds !'.format(time.time() - startTime))

if __name__ == '__main__':
	main()