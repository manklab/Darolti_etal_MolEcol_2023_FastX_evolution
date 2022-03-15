#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
''' Takes the filtered SNP file (output of alielw/McDonald-Kreitman-test/09.filter-snps) and identifies the codon corresponding to each SNP.
	It finds the encoding protein of each codon when the reference, and separately when the alternative allele, of each SNP is used
	Identifies proteins that differ when replacing the reference with the alternative allele (nonsynonymous SNPs) and those that stay the same (synonymous SNPs) '''
#==============================================================================
import argparse
import sys
import os
from Bio.Seq import Seq
from collections import defaultdict
import time
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
					help="A filtered snp file"
					"Output of alielw/McDonald-Kreitman-test/09.filter-snps")
parser.add_argument("coordinates", type=str,
					help="The output file of alielw/McDonald-Kreitman-test/01.orf-extract-coordinates")
parser.add_argument("gene_position", type=str,
					help="A file containing positional information on each gene ie which scaffold it is located on"
					"Output of alielw/Quantify-gene-expression-reference-based/12.extract-counts-annotated.py"
					"Format: Geneid,Scaffold,Chromosome,Startwithinscaffold,Startwithingenome")
parser.add_argument("ORF_sites_in_scaffold", type=str,
					help="A file containing genes and the position of sites within the ORF in the scaffold"
					"Output of alielw/McDonald-Kreitman-test/02.orf-position-within-scaffold.py")
parser.add_argument("genome_assembly", type=str,
					help="The 1k filtered genome assembly")
parser.add_argument("outfile", type=str,
					help="An outfile containing gene name, NpN and SpS")

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
	gene_orthogroup = {}
	with open(coordinates, "r") as coord:
		# Get dictionary of gene name and orthogroup name
		next(coord) # skips header
		for line in coord:
			line = line.rstrip().split("\t")
			target_gene_name = line[0].split("_")[1].split(".")[0] + "." + line[0].split("_")[1].split(".")[1]
			orthogroup_name = line[1].split("_")[1]
			# orthogroup_name = line[1].split("_")[1] + "_" + line[1].split("_")[2]
			gene_orthogroup[target_gene_name] = orthogroup_name
	return gene_orthogroup

def get_dictionary_which_gene_scaffold(gene_position, gene_orthogroup):
	gene_scaffold_dict = {}
	scaffolds_to_keep = []
	with open(gene_position, "r") as gene_infile:
		# Get dictionary of gene and which scaffold it is located on
		next(gene_infile) # skips header
		for line in gene_infile:
			line = line.rstrip()
			gene = line.split(",")[0]
			scaffold = line.split(",")[1]
			# Use the gene_orthogroup output to filter out the genes without orthogroup
			if gene in gene_orthogroup:
				gene_scaffold_dict[gene] = scaffold
				# Make list of valid scaffolds to make the get_sequence function faster
				scaffolds_to_keep.append(scaffold)
	return gene_scaffold_dict, scaffolds_to_keep

def get_sequence(genome_assembly, scaffolds_to_keep):
	scaffold_sequence = defaultdict(list)
	with open(genome_assembly, 'r') as assembly:
		for line in assembly:
			line = line.rstrip()
			if line.startswith(">"):
				scaffold = "scaffold" + line.split(">")[1].split(" ")[0]
			else:
				# Only keep scaffolds in valid scaffold list to make the function faster
				if scaffold in scaffolds_to_keep:
					for site in line:
						scaffold_sequence[scaffold].append(site)
	return scaffold_sequence

def get_ORF_sites_in_scaffold(ORF_sites_in_scaffold):
	ORF_sites_in_scaffold_dict = {}
	with open(ORF_sites_in_scaffold, "r") as infile:
		for line in infile:
			line = line.rstrip().split("\t")
			gene = line[0].split(".")[0]+"."+line[0].split(".")[1]
			ORF_sites = line[1].split("]")[0].split("[")[1].split(", ")
			ORF_sites_in_scaffold_dict[gene] = ORF_sites
	return ORF_sites_in_scaffold_dict

def get_ORF_sequence(ORF_sites_in_scaffold_dict, scaffold_sequence, gene_scaffold_dict):
	gene_ORF_nucleotides = defaultdict(list)
	for gene in ORF_sites_in_scaffold_dict:
		scaffold = gene_scaffold_dict[gene]
		for ORF_site in ORF_sites_in_scaffold_dict[gene]:
			# python counting
			nucleotide_pos = int(ORF_site) - 1
			nucleotide = scaffold_sequence[scaffold][nucleotide_pos]
			gene_ORF_nucleotides[gene].append(nucleotide)
	return gene_ORF_nucleotides

def get_dictionary_gene_snps(snps):
	gene_snps = defaultdict(list)
	with open(snps, "r") as snps:
		next(snps)
		for line in snps:
			snp_info = []
			line = line.rstrip().split("\t")
			gene = line[4]
			position = line[1]
			ref = line[2]
			alt = line[3]
			snp_info.append(position)
			snp_info.append(ref)
			snp_info.append(alt)
			gene_snps[gene].append(snp_info)
	return gene_snps

def get_gene_snps_info(ORF_sites_in_scaffold, gene_snps):
	gene_snps_info = defaultdict(list)
	with open(ORF_sites_in_scaffold, "r") as infile:
		for line in infile:
			line = line.rstrip().split("\t")
			gene = line[0].split(".")[0]+"."+line[0].split(".")[1]
			sites = line[1].split("]")[0].split("[")[1].split(", ")
			position_dictionary, site_codon_dict = identify_codons(sites)
			for info in gene_snps[gene]:
				snp_position_in_scaffold = info[0]
				ref = info[1]
				alt = info[2]
				position_in_codon = 0
				# Get the codon in which the snp scaffold position lies
				codon = site_codon_dict[snp_position_in_scaffold]
				snp_info = []
				for site in position_dictionary[codon]:
					if float(site) == float(snp_position_in_scaffold):
						snp_info.append(snp_position_in_scaffold)
						snp_info.append(codon)
						snp_info.append(position_in_codon)
						snp_info.append(ref)
						snp_info.append(alt)
						gene_snps_info[gene].append(snp_info)
					position_in_codon += 1
	return gene_snps_info

def identify_codons(sites):
	ORF_nucleotides = sites
	# Check length of ORF is a multiple of three
	if len(ORF_nucleotides) % 3:
		print("ERROR - ORF nucleotides not selected correctly")
	# Get list of codon positions - human counting
	ORF_codons = range(1,(len(ORF_nucleotides)//3)+1)
	# Check there are a third the number of codons than nucleotides in ORF
	if len(ORF_nucleotides)//3 != len(ORF_codons):
		print("ERROR - codon position not correct")
	# No need to reverse the nucleotides if reverse complement because already reversed in 02.orf-position-within-scaffold.py
	# Split nucleotides into triplets
	triplets = [ORF_nucleotides[i:i +3] for i in range(0, len(ORF_nucleotides), 3)]
	# Check there are the same number of triplets and codons
	if len(triplets) != len(ORF_codons):
		print("ERROR - codon position not correct")
	# Make dictionary between codon position and triplet position
	position_dictionary = {}
	site_codon_dict = {}
	for codon in ORF_codons:
		triplet = triplets[codon-1] # Codon is in human counting. Need to substract 1 as looking up position in triplets.
		position_dictionary[codon] = triplet # position_dictionary is in human counting
		for site in triplet:
			site_codon_dict[site] = codon
	return position_dictionary, site_codon_dict
			

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	startTime = time.time()

	# Make dictionary of gene and orthogroup
	gene_orthogroup = get_dictionary_gene_orthogroup(args.coordinates)
	print("No of genes with orthogroup =", len(gene_orthogroup))

	# Make dictionary of gene and scaffold
	gene_scaffold_dict, scaffolds_to_keep = get_dictionary_which_gene_scaffold(args.gene_position, gene_orthogroup)
	print("No of valid genes with positional information =", len(gene_scaffold_dict))

	# Make dictionary of valid scaffolds and their sequence
	scaffold_sequence = get_sequence(args.genome_assembly, scaffolds_to_keep)
	print("No of valid scaffolds with sequence info =", len(scaffold_sequence))

	# Make dictionary of valid genes and the position in scaffold of all sites within ORF
	ORF_sites_in_scaffold_dict = get_ORF_sites_in_scaffold(args.ORF_sites_in_scaffold)
	print("No of genes with ORF sites before any filtering =", len(ORF_sites_in_scaffold_dict))

	# Make dictionary of gene and nucleotides in ORF
	gene_ORF_nucleotides = get_ORF_sequence(ORF_sites_in_scaffold_dict, scaffold_sequence, gene_scaffold_dict) 
	print("No of genes with ORF nucleotides =", len(gene_ORF_nucleotides))

	# Make dictionary of gene, snp position, reference and alternative alleles
	gene_snps = get_dictionary_gene_snps(args.infile)

	# Make dictionary of gene and snp detailed info. 
	# Format: Gene_name :['SNP_position_in_scaffold', 'Codon_in_ORF', 'Position_in_codon', 'Reference_allele', 'Alternative_allele']
	gene_snps_info = get_gene_snps_info(args.ORF_sites_in_scaffold, gene_snps)
	print("No of genes with SNPs =", len(gene_snps_info))

	with open(args.outfile, "w") as out:
		out.write("Gene")
		out.write("\t")
		out.write("NpN")
		out.write("\t")
		out.write("SpS")
		out.write("\n")
		for gene in gene_snps_info:
			count_synonymous_snps = 0
			count_nonsynonymous_snps = 0
			# Get nucleotides in ORF for each gene with SNPs
			nucleotides = gene_ORF_nucleotides[gene]
			# Split the nucleotides into codons to identify the corresponding nucleotide for each SNP 
			position_dictionary, site_codon_dict = identify_codons(nucleotides)
			for snp in gene_snps_info[gene]:
				snp_position_in_scaffold = snp[0]
				snp_codon = snp[1]
				snp_position_in_codon = snp[2]
				snp_ref = snp[3]
				snp_alt = snp[4]
				# Check that the nucleotide at the position of each SNP corresponds to that SNP's reference allele
				# Use .upper() as some of the nucleotides in the sequence will be lowercase and will cause Error 
				if position_dictionary[snp_codon][snp_position_in_codon].upper() != snp_ref:
					print("Error! - SNP nucleotide does not correspond to SNP reference allele")
					print(gene, "SNP pos in scaffold =", snp_position_in_scaffold, "SNP codon =", snp_codon, "SNP pos in codon =", snp_position_in_codon)
					print(position_dictionary[snp_codon][snp_position_in_codon], snp_ref)
					# sys.exit()

				# Find amino acid sequence containing the SNP reference allele and the protein it encodes
				ref_codon = position_dictionary[snp_codon][0] + position_dictionary[snp_codon][1] + position_dictionary[snp_codon][2]
				ref_protein = Seq(ref_codon).translate()
				
				# Find amino acid sequence containing the SNP alternative allele and the protein it encodes
				counter = 0
				alt_codon = [""]
				for site in position_dictionary[snp_codon]:
					if counter == snp_position_in_codon:
						alt_codon[0] += snp_alt
					else:
						alt_codon[0] += site
					counter += 1
				for codon in alt_codon:
					alt_codon = codon
					alt_protein = Seq(alt_codon).translate()

				# Count synonymous and nonsynonymous snps for each gene
				if ref_protein != alt_protein:
					count_nonsynonymous_snps += 1
				else:
					count_synonymous_snps += 1
			out.write(gene)
			out.write("\t")
			out.write(str(count_nonsynonymous_snps))
			out.write("\t")
			out.write(str(count_synonymous_snps))
			out.write("\n")

	print ('The script too {0} seconds !'.format(time.time() - startTime))

if __name__ == '__main__':
    main()
