#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
''' Check that position of codons passing gap and ambiguity filter are correct
Takes folder of folders of sequences aligned by PRANK before gap removal. Identifies codons in ORF which pass both the
gap and ambiguity filter. Filters the PRANK alignment on the basis of these codons and checks each is identical 
to the final masked alignments where both gaps and ambiguity data have been removed. Thereby ensuring codon positions 
have been correctly identified.'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder_prefiltering", type=str,
                    help="An input folder containing folders of fasta sequences of orthogroups aligned by PRANK before gap removal")
parser.add_argument("infolder_postfiltering", type=str,
                    help="An input folder containing folders of phylip sequences of orthogroups aligned after gap removal, masking and N removal"
                    "Same folder than PAML branch site and branch model tests would be run on")
parser.add_argument("target_species", type=str,
                    help="Species for which divergence data will be compared to polymorphism data")
parser.add_argument("codons_without_gaps", type=str,
                    help="A file containing genes with codons where no site is in an alignment gap")
parser.add_argument("codons_without_Ns", type=str,
                    help="A file containing genes with codons where no site is an N")
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

def list_files(current_dir, file_ending):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if not name.endswith(".DS_Store"):
                if name.endswith(file_ending):
                    f = os.path.join(path, name)
                    file_list.append(f)
    return file_list

def extract_alignment_files(infolder, file_ending):
    aln_file_dict = defaultdict(list)
    for orthogroup in list_folder(infolder):
        for aln in list_files(orthogroup,file_ending):
            name = os.path.basename(aln).split(file_ending)[0]
            fasta = read_fasta(aln)
            aln_file_dict[name] = fasta
    return aln_file_dict

def read_fasta(source):
    ''' Read fasta file and outputs a dictionary 
    key [seq name] values [sequence] [header]'''
    SG = {}
    try:
        with open(source, "r") as file:
            for line in file.readlines():
                if line[0] == ">":
                    name = line[1:].rstrip().split()[0]
                    header = line[1:].rstrip()
                    SG[name] = ["", header]
                else:
                    SG[name][0] += line.rstrip()
            return SG
    except IOError:
        print "File does not exit!"

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

def remove_target_species_gaps(fasta):
    for s in fasta:
        if s == args.target_species:
            target_species_gaps = {}
            for i, p in enumerate(fasta[s][0]):
                if p == "-":
                    target_species_gaps[i] = 0
    # Store new alignments with target species gaps removed
    # Alignment now is equvialent to ORF of target species ie starts and ends at same position
    fasta_target_gaps_removed = {}
    for s in fasta:
        fasta_target_gaps_removed[s] = [""]
        for i, fs in enumerate(fasta[s][0]):
            if i not in target_species_gaps:
                fasta_target_gaps_removed[s][0] += fs
    return fasta_target_gaps_removed

#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():

    # Extract alignments without gaps removed
    aln_file_dict_prefiltering = extract_alignment_files(args.infolder_prefiltering,".stripped.prank.best.fas")
    print "Number of alignments before gaps removed =", len(aln_file_dict_prefiltering)

    # Extract alignments after gap removal, masking and ambiguity data removal
    aln_file_dict_postfiltering = extract_alignment_files(args.infolder_postfiltering,"_masked.Nrm.fa")
    print "Number of alignments after alignment filtering =", len(aln_file_dict_postfiltering)

    # Extract codons without gaps
    # Position refers to ORF of target species
    codons_without_gaps = get_codons(args.codons_without_gaps)
    print "Number of ORF with codons without gaps =", len(codons_without_gaps)

    # Extract codons without Ns
    # Position refers to ORF of target species
    codons_without_Ns = get_codons(args.codons_without_Ns)  
    print "Number of ORF with codons without Ns =", len(codons_without_Ns)     

    # Filter alignment to keep only codons passing both filters
    count_same_aln = 0
    for aln in aln_file_dict_prefiltering:
        # Find codons passing both gap and N filter
        no_gaps = codons_without_gaps[aln]
        no_Ns = codons_without_Ns[aln]
        no_gaps_and_Ns = list(set(no_gaps).intersection(no_Ns))
        # First remove gaps within target ORF. These are not listed in no_gaps_and_Ns as they are not in ORF
        fasta_target_gaps_removed = remove_target_species_gaps(aln_file_dict_prefiltering[aln])
        # Keep only codons passing both filters
        for s in fasta_target_gaps_removed:
            no_gaps_and_Ns_filtered_seq = [""]
            # Identify sites in sequence and triplets 
            triplets = [fasta_target_gaps_removed[s][0][i:i +3] for i in range(0, len(fasta_target_gaps_removed[s][0]), 3)]
            for i, fs in enumerate(triplets):
                i = i+1
                if i in no_gaps_and_Ns:
                    no_gaps_and_Ns_filtered_seq[0] += fs
            #Check that filtered alignment is identical to post filtered alignment used in PAML
            if aln in aln_file_dict_postfiltering:
                aln_file_postfiltered = aln_file_dict_postfiltering[aln]
                aln_postfiltered_seq = aln_file_postfiltered[s][0]
                if no_gaps_and_Ns_filtered_seq[0] == aln_postfiltered_seq:
                    count_same_aln += 1
                else:
                    print "ERROR"
                    print aln

if __name__ == '__main__':
    main()
