#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

''' Extract branch lengths
This script takes an infolder containing folders with the output of branch model test 
and outputs branch lengths for the specified branch. Filters orthologs with branch-specific dS >= 2
and if SdS for that branch is <= 1.'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
					help="Folder containing folders with paml output files")
parser.add_argument("branchofinterest", type=str,
					help="Branch for which branch lengths should be extracted eg 5..6")
parser.add_argument("outfile", type=str,
					help="")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def list_folder(infolder):
	'''Returns a list of all files in a folder with the full path'''
	return [os.path.join(infolder, f) for f in os.listdir(infolder) if not f.endswith(".DS_Store")]

def list_files(current_dir):
	file_list = []
	for path, subdirs, files in os.walk(current_dir): # Walk directory tree
		for name in files:
			if name.endswith(".txt"):
				f = os.path.join(path, name)
				file_list.append(f)
	return file_list

def read_paml_branch(source,branchid):
	# branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS
	''' PAML branch parser
	Reads a PAML output file and extracts the branch of interest.
	'''
	try:
		with open(source, "r") as infile:
			branch_lengths = []
			for line in infile:
				line = line.rstrip()
				if line.startswith(branchid):
					line = line.split(" ")
					for l in line:
						if is_number(l):
							branch_lengths.append(l)
			if len(branch_lengths) != 8:
				print "ERROR - extracted wrong number of lengths"
				print branch_lengths
		return branch_lengths
	except IOError:
		print "File %s does not exit!" (source)
		raise

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
#==============================================================================
#Main==========================================================================
#==============================================================================

def main():
	#get files
	infolders = list_folder(args.infolder)
	print "Number of infolders =", len(infolders)

	#get id
	branchid = "   "+args.branchofinterest
	print "Finding branch ...",branchid

	#extract branch lengths
	branch_length_dict = {}
	for infolder in infolders:
		for infile in list_files(infolder):
			genename = os.path.dirname(infile).split("_")[-1]
			# t       N       S   dN/dS      dN      dS  N*dN  S*dS
			branch_lengths = read_paml_branch(infile,branchid)
			branch_length_dict[genename] = branch_lengths
	print "Number of genes with branch lengths =", len(branch_length_dict)

	#filter dS>2 and S*dS > 1
	branch_length_dict_filtered = {}
	for gene in branch_length_dict:
		branch_lengths = branch_length_dict[gene]
		dS = float(branch_lengths[5])
		SdS = float(branch_lengths[-1])
		if dS < 2:
			if SdS > 1:
				branch_length_dict_filtered[gene] = branch_length_dict[gene]

	print "Number of genes with branch lengths after dS <2 and S*dS >1 filter =", len(branch_length_dict_filtered)

	with open(args.outfile, "w") as outfile:
		count = 0
		outfile.write("Gene,t,N,S,dN/dS,dN,dS,N*dN,S*dS")
		outfile.write("\n")
		for gene in branch_length_dict_filtered:
			count += 1
			sys.stdout.write('%d\r' % (count))
			sys.stdout.flush()
			outfile.write(gene)
			for g in branch_length_dict_filtered[gene]:
				outfile.write(",")
				outfile.write(str(g))
			outfile.write("\n")

	print "Number of files processed =", count

if __name__ == "__main__":
	main()