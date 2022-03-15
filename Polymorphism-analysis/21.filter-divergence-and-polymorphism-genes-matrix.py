#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
''' Takes the divergence and polymorphism files and finds the genes that are common. 
	Creates an output file containing only the genes in common together with their divergence and polymorphism data'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
import time
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("divergence_file", type=str,
					help="An infile containing orthogroup name, t, N, S, dN/dS, dN, dS, N*dN, S*dS"
					"The output file of 14.paml-extract-branch-lengths-dSSdSfilter.py")
parser.add_argument("polymorphism_file", type=str,
					help="An infile containing gene name, NpN and SpS"
					"the output file of 16.identify-NpN-SpS.py")
parser.add_argument("coordinates", type=str,
					help="Containing correspondence between target gene name and orthogroup name"
					"The output file of 07.orf-position-within-scaffold.py")
parser.add_argument("outfile", type=str,
					help="An outfile containing target gene name, Gene, Dn, Ds, Pn, Ps")

# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

def get_dictionary_gene_orthogroup(coordinates):
	orthogroup_gene = {}
	with open(coordinates, "r") as coord:
		# Get dictionary of gene name and orthogroup name
		next(coord) # skips header
		for line in coord:
			line = line.rstrip().split("\t")
			target_gene_name = line[0].split("_")[1].split(".")[0] + "." + line[0].split("_")[1].split(".")[1]
			orthogroup_name = line[1].split("_")[1]
			orthogroup_gene[orthogroup_name] = target_gene_name
	return orthogroup_gene

def get_divergence_data(divergence_file, orthogroup_gene):
	gene_divergence = defaultdict(list)
	with open(divergence_file, "r") as divergence:
		next(divergence)
		for line in divergence:
			line = line.rstrip()
			orthogroup = line.split(",")[0].split("/")[-1]
			N_dN = line.split(",")[7]
			S_dS = line.split(",")[8]
			if orthogroup in orthogroup_gene:
				gene = orthogroup_gene[orthogroup]
				gene_divergence[gene].append(N_dN)
				gene_divergence[gene].append(S_dS)
	return gene_divergence

def get_polymorphism_data(polymorphism_file):
	gene_polymorphism = defaultdict(list)
	with open(polymorphism_file, "r") as polymorphism:
		next(polymorphism)
		for line in polymorphism:
			line = line.rstrip()
			gene = line.split("\t")[0]
			NpN = line.split("\t")[1]
			SpS = line.split("\t")[2]
			gene_polymorphism[gene].append(NpN)
			gene_polymorphism[gene].append(SpS)	
	return 	gene_polymorphism

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	startTime = time.time()

	# Make dictionary of gene and orthogroup
	orthogroup_gene = get_dictionary_gene_orthogroup(args.coordinates)
	print "Number of genes with orthogroup =", len(orthogroup_gene)

	# Make dictionary of target gene name and divergence data
	gene_divergence = get_divergence_data(args.divergence_file, orthogroup_gene)
	print "Number of orthogroups with divergence data =", len(gene_divergence)

	# Make dictionary of gene namd and polymorphism information
	gene_polymorphism = get_polymorphism_data(args.polymorphism_file)
	print "No of genes with polymorphsim data =", len(gene_polymorphism)

	with open(args.outfile, "w") as out:
		out.write("Gene,Dn,Ds,Pn,Ps")
		out.write("\n")
		genes = list(set(gene_divergence).intersection(gene_polymorphism))
		print "Number of genes with both divergence and polymorphism data =", len(genes)
		
		for gene in genes:
			out.write(gene)
			out.write(",")
			for info in gene_divergence[gene]:
				out.write(info)
				out.write(",")
			count = 0 
			for info in gene_polymorphism[gene]:
				count += 1
				out.write(info)
				if count != len(gene_polymorphism[gene]):
					out.write(",")
			out.write("\n")

	print "The script too {0} seconds !".format(time.time() - startTime)

if __name__ == '__main__':
    main()