#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
''' Keeps SNPs from a SNP input file that are positioned in codons without gaps, without Ns and that pass the coverage threshold '''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
					help="An snp file")
parser.add_argument("coordinates", type=str,
					help="The output file of 06.orf-position-within-transcript.py")
parser.add_argument("gene_position", type=str,
					help="A file containing positional information on each gene ie which scaffold it is located on"
					"Format: Geneid,Scaffold,Chromosome,Startwithinscaffold,Startwithingenome")
parser.add_argument("codons_without_gaps", type=str,
                    help="A file containing genes with codons where no site is in an alignment gap")
parser.add_argument("codons_without_Ns", type=str,
                    help="A file containing genes with codons where no site is an N")
parser.add_argument("codons_with_coverage", type=str,
                    help="A file containing genes with codons where no site fails the coverage threshold")
parser.add_argument("ORF_sites_in_scaffold", type=str,
					help="A file containing genes and the position of sites within the ORF in the scaffold"
					"Output of 07.orf-position-within-scaffold.py")
parser.add_argument("outfile", type=str,
					help="An outfile file of kept SNPs")

# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

def get_dictionary_which_gene_scaffold(gene_position):
	scaffold_gene_dict = defaultdict(list)
	with open(gene_position, "r") as gene_infile:
		# Get dictionary of gene and which scaffold it is located on
		next(gene_infile) # skips header
		for line in gene_infile:
			line = line.rstrip()
			gene = line.split(",")[0]
			scaffold = line.split(",")[1]
			scaffold_gene_dict[scaffold].append(gene)
	return scaffold_gene_dict

def get_dictionary_gene_orthogroup(coordinates):
	gene_orthogroup = {}
	with open(coordinates, "r") as coord:
		# Get dictionary of gene name and orthogroup name
		next(coord) # skips header
		for line in coord:
			line = line.rstrip().split("\t")
			target_gene_name = line[0].split("_")[1].split(".")[0] + "." + line[0].split("_")[1].split(".")[1]
			orthogroup_name = line[1].split("_")[1]
			gene_orthogroup[target_gene_name] = orthogroup_name
	return gene_orthogroup

def get_ORF_sites_in_scaffold(ORF_sites_in_scaffold):
	ORF_sites_in_scaffold_dict = {}
	with open(ORF_sites_in_scaffold, "r") as infile:
		for line in infile:
			line = line.rstrip().split("\t")
			gene = line[0].split(".")[0]+"."+line[0].split(".")[1]
			sites = line[1].split("]")[0].split("[")[1].split(", ")
			ORF_sites_in_scaffold_dict[gene] = sites
	return ORF_sites_in_scaffold_dict

def get_codons(file):
    # Get dictionary of codons to be kept
    codons_to_keep = defaultdict(list)
    with open(file, "r") as infile:
        for line in infile:
            line = line.rstrip().split("\t")
            name = line[0]
            codons_keep = line[1].split("[")[1].split("]")[0].split(",")
            for codon in codons_keep:
                codon = float(codon)
                codons_to_keep[name].append(codon)
    return codons_to_keep

def get_gene_names(coordinates, codons_with_coverage):
    codons_with_coverage_ortho_gene_name = {}
    with open(coordinates, "r") as coord:
        next(coord) # skips header
        for line in coord:
            line = line.rstrip().split("\t")
            target_gene_name = line[0].split("_")[1].split(".")[0] + "." + line[0].split("_")[1].split(".")[1]
            orthogroup_gene_name = line[1].split("_")[1]
            if target_gene_name in codons_with_coverage:
                codons_with_coverage_ortho_gene_name[orthogroup_gene_name] = codons_with_coverage[target_gene_name]
            else:
                codons_with_coverage_ortho_gene_name[orthogroup_gene_name] = []
    return codons_with_coverage_ortho_gene_name

def get_sites_codons(valid_genes_codons, ORF_sites_in_scaffold_dict):
	gene_sites_in_valid_codon = defaultdict(list)
	gene_sites_in_valid_codon_check = defaultdict(list)
	for gene in valid_genes_codons:
		sites = ORF_sites_in_scaffold_dict[gene]
		position_dictionary = identify_codons(sites)
		for codon in position_dictionary:
			if codon in valid_genes_codons[gene]:
				for site in position_dictionary[codon]:
					gene_sites_in_valid_codon[gene].append(site)
	return gene_sites_in_valid_codon

def identify_codons(sites):
	ORF_nucleotides = sites
	# Check length of ORF is a multiple of three
	if len(ORF_nucleotides) % 3:
		print "ERROR - ORF nucleotides not selected correctly"
	# Get list of codon positions - human counting
	ORF_codons = range(1,(len(ORF_nucleotides)//3)+1)
	# Check there are a third the number of codons than nucleotides in ORF
	if len(ORF_nucleotides)//3 != len(ORF_codons):
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

def remove_invalid_snps(infile, scaffold_gene_dict, gene_sites_in_valid_codon, outfile):
	count_total_SNPs = 0
	count_valid_SNPs = 0
	count_missing_scaffold = 0
	with open(outfile, 'w') as out:
		out.write("CHROM")
		out.write("\t")
		out.write("POS")
		out.write("\t")
		out.write("REF")
		out.write("\t")
		out.write("ALT")
		out.write("\t")
		out.write("GENE_NAME")
		out.write("\n")
		with open(infile, 'r') as infile:
			next(infile)
			for line in infile:
				printed = None
				count_total_SNPs += 1
				line = line.rstrip().split("\t")
				scaffold = "scaffold" + line[0]
				position = line[1]
				ref = line[2]
				alt = line[3]
				if scaffold in scaffold_gene_dict:
					genes = scaffold_gene_dict[scaffold]
					for gene in genes:
						if position in gene_sites_in_valid_codon[gene]:
							if printed is None:
								printed = "yes"
								count_valid_SNPs += 1
								out.write(scaffold)
								out.write("\t")
								out.write(position)
								out.write("\t")
								out.write(ref)
								out.write("\t")
								out.write(alt)
								out.write("\t")
								out.write(gene)
								out.write("\n")	
				else:
					if scaffold.startswith("s"):
						count_missing_scaffold += 1
			print count_missing_scaffold 
	return count_total_SNPs, count_valid_SNPs

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	# Make dictionary of gene and scaffold
	scaffold_gene_dict = get_dictionary_which_gene_scaffold(args.gene_position)
	print "Number of genes with positional information =", len(scaffold_gene_dict) 

	# Make dictionary of gene and orthogroup
	gene_orthogroup = get_dictionary_gene_orthogroup(args.coordinates)
	print "Number of genes with orthogroup =", len(gene_orthogroup) 

	# Make dictionary of valid genes and the position in scaffold of all sites within ORF
	ORF_sites_in_scaffold_dict = get_ORF_sites_in_scaffold(args.ORF_sites_in_scaffold)
	print "Number of genes with ORF =", len(ORF_sites_in_scaffold_dict) 

	# Extract codons without gaps
	# Position refers to ORF of target species
	codons_without_gaps = get_codons(args.codons_without_gaps)
	print "Number of ORF with codons without gaps =", len(codons_without_gaps) 

	# Extract codons without Ns
	# Position refers to ORF of target species
	codons_without_Ns = get_codons(args.codons_without_Ns)  
	print "Number of ORF with codons without Ns =", len(codons_without_Ns)   

	# Extract codons with coverage
	# Position refers to ORF of target species
	# Note: the created dictionary codons_with_coverage contains gene names instead of orthogroup names
	codons_with_coverage = get_codons(args.codons_with_coverage)  
	print "Number of ORF with codons with coverage =", len(codons_with_coverage)  

	# Use dictionary between orthogroup names and target species' gene names
	# and create dictionary between orthogroup name and codons that pass coverage threshold (codons kept in codons_with_coverage)
	# This is needed to be able to make the intersection between the 3 codon dictionaries
	codons_with_coverage_ortho_gene_name = get_gene_names(args.coordinates, codons_with_coverage)
	print "Number of ORF with codons with coverage ortho names =", len(codons_with_coverage_ortho_gene_name) 

	valid_genes_codons = {}
	for gene in ORF_sites_in_scaffold_dict:
		no_gaps = codons_without_gaps[gene_orthogroup[gene]]
		no_Ns = codons_without_Ns[gene_orthogroup[gene]]
		coverage = codons_with_coverage_ortho_gene_name[gene_orthogroup[gene]]
		no_gaps_Ns_and_coverage = list(set(no_gaps).intersection(no_Ns).intersection(coverage))
		if len(no_gaps_Ns_and_coverage) != 0:
			valid_genes_codons[gene] = no_gaps_Ns_and_coverage
	print "Number genes that have valid codons (no gaps, no Ns, pass coverage) =", len(valid_genes_codons)

	# Make dictionary of genes and sites for codons in valid_genes_codons
	gene_sites_in_valid_codon = get_sites_codons(valid_genes_codons, ORF_sites_in_scaffold_dict)

	# Remove invalid SNPs
	count_total_SNPs, count_valid_SNPs = remove_invalid_snps(args.infile, scaffold_gene_dict, gene_sites_in_valid_codon, args.outfile)
	print "Total number of SNPs =", count_total_SNPs 
	print "Number of valid SNPs =", count_valid_SNPs 

if __name__ == '__main__':
    main()