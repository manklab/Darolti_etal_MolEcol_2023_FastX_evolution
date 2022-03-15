#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
'''  Get coordinates of open reading frame within scaffold
This script uses the output from 06.orf-position-within-transcript.py and, using the gtf file of
genes' position in scaffolds, calculates the coordinates of the open reading frame within the scaffold.'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
from Bio.Seq import Seq
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("coordinates", type=str,
					help="A file containing position of first and last nucleotide within ORF"
					"Output of 06.orf-position-within-transcript.py")
parser.add_argument("gene_position", type=str,
					help="A file containing positional information on each gene ie which scaffold it is located on"
					"Format: Geneid,Scaffold,Chromosome,Startwithinscaffold,Startwithingenome")
parser.add_argument("gtf", type=str,
					help="A gtf file")
parser.add_argument("outfile", type=str,
					help="Outfile containing coordinates of open reading frame within scaffold")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================

def coordinate_dictionary(coordinates, gene_position):
	coordinate_dict = defaultdict(list)
	gene_scaffold_dict = {}
	with open(coordinates, "r") as coord_infile, open(gene_position, "r") as gene_infile:
		# Get dictionary of gene and which scaffold it is located on
		for line in gene_infile:
			if not line.startswith("Gene"):
				line = line.rstrip()
				gene = line.split(",")[0]
				scaffold = line.split(",")[1]
				gene_scaffold_dict[gene] = scaffold
		print "No of genes with positional information =", len(gene_scaffold_dict)
		# Get dictionary of gene and scaffold, direction, start and stop within transcript
		for line in coord_infile:
			if not line.startswith("Gene"):
				line = line.rstrip()
				transcript = line.split("\t")[0].split("_")[1]
				gene = transcript.split(".")[0]+"."+transcript.split(".")[1]
				direction = line.split("\t")[2]
				start = line.split("\t")[3]
				stop = line.split("\t")[4]
				# Only select genes which are assigned to scaffolds with known chromosomes
				if gene in gene_scaffold_dict:
					scaffold = gene_scaffold_dict[gene]
					coordinate_dict[transcript].append(scaffold)
					coordinate_dict[transcript].append(direction)
					coordinate_dict[transcript].append(start)
					coordinate_dict[transcript].append(stop)
	print "Number of genes with ORF =", len(coordinate_dict)
	return coordinate_dict

def get_exon_start_stop(stringtie):
	exon_start_stop_in_scaffold = defaultdict(list)
	gene_size = {}
	with open(stringtie, 'r') as stringtie:
		for line in stringtie:
			line = line.rstrip()
			if not line.startswith("#"):
				if line.split("\t")[2] == "transcript":
					scaffold = line.split("\t")[0]
					transcript = line.split("\t")[8].split(";")[1].split('"')[1]
					gene_size[transcript] = 0
				else:
					transcript_id = line.split("\t")[8].split(";")[1].split('"')[1]
					if line.split("\t")[2] == "exon" and transcript_id == transcript:
						start_exon = int(line.split("\t")[3])
						stop_exon = int(line.split("\t")[4])
						# Calculate moving sum of the size of the exons
						exon_size = stop_exon - start_exon + 1
						gene_size[transcript] += exon_size
						pos = [start_exon, stop_exon]
						exon_start_stop_in_scaffold[transcript_id].append(pos)
					else:
						print "ERROR!"
						print line
						print "transcript", transcript
						sys.exit()
	return exon_start_stop_in_scaffold, gene_size

def ORF_positions(exon_start_stop_in_scaffold, gene_size, coordinate_dict):
	ORF_pos_in_scaffold = defaultdict(list)
	for transcript in exon_start_stop_in_scaffold:
		if transcript in coordinate_dict:
			sites_in_exons_in_scaffold = []
			# Storing position of every site in all exons into list
			for exon in exon_start_stop_in_scaffold[transcript]:
				transcript = transcript.split("_")[0]
				for i in range(int(exon[0]), int(exon[1]+1)):
					sites_in_exons_in_scaffold.append(i)
			# Checking lengths are correct
			if len(sites_in_exons_in_scaffold) != gene_size[transcript]:
					print "ERROR!"
					sys.exit()
			# Make dictionary between position of site in transcript and sites in exon within scaffold
			dict_sites_ORF_scaffold = {}
			for i, scaffold_site in enumerate(sites_in_exons_in_scaffold):
				i = i+1
				dict_sites_ORF_scaffold[i] = scaffold_site
			# Identify where sites in the ORF are located in the scaffold
			pos_in_ORF = defaultdict(list)
			start = int(coordinate_dict[transcript][2])
			stop = int(coordinate_dict[transcript][3])
			#Reverse order of sites that are printed if reverse complement
			if coordinate_dict[transcript][1] == "reversecomplement":
				while stop >= start:
					site = dict_sites_ORF_scaffold[stop]
					ORF_pos_in_scaffold[transcript].append(site)
					stop -= 1
			else:
				while start <= stop:
					site = dict_sites_ORF_scaffold[start]
					ORF_pos_in_scaffold[transcript].append(site)
					start +=1
	return ORF_pos_in_scaffold

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	# Get dictionary of coordinates of ORF within the gene
	coordinate_dict = coordinate_dictionary(args.coordinates, args.gene_position)

	# Get dictionary of start and stop postions within the scaffold of each exon within a gene
	exon_start_stop_in_scaffold, gene_size = get_exon_start_stop(args.gtf)

	# For every site in the ORF, identify its location within the scaffold
	ORF_pos_in_scaffold = ORF_positions(exon_start_stop_in_scaffold, gene_size, coordinate_dict)

	with open(args.outfile, "w") as outfile:
		for gene in ORF_pos_in_scaffold:
			outfile.write(gene)
			outfile.write("\t")
			outfile.write(str(ORF_pos_in_scaffold[gene]))
			outfile.write("\n")

if __name__ == '__main__':
	main()