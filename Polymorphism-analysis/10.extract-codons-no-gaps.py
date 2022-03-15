#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
''' Identifies codons that are not in gaps
Takes folder of folders of sequences aligned by PRANK before gaps have been removed.
Outfile codons_without_gaps reports codons where all sites are not in a gap. i.e. codons to keep
Outfile codons_with_gaps reports codons where any site is in a gap. i.e. codons to remove. 
Gaps before or within target species sequence are not included in either list because they are not in the ORF and do not need to be removed.
Position reported in relation to codons within ORF. i.e. 1 = 1st codon in ORF/protein.'''
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
                    help="An input folder containing folders of sequences of orthogroups aligned by PRANK before gap removal")
parser.add_argument("target_species", type=str,
                    help="Species for which divergence data will be compared to polymorphism data")
parser.add_argument("outfile", type=str,
                    help="An outfile file of codons passing gap filtering")
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
            if not name.endswith(".DS_Store"):
                if name.endswith("prank.best.fas"):
                    f = os.path.join(path, name)
                    file_list.append(f)
    return file_list

def read_fasta(source):
    ''' Read fasta file and outputs a dictionary 
    key [seq name] values [sequence] [header] '''
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

def identify_sites_with_gaps(fasta):
    sites_with_gaps = {}
    for s in fasta:
        for i, p in enumerate(fasta[s][0]):
            if p == "-":
                sites_with_gaps[i] = 0
    return sites_with_gaps

def identify_codons_gaps(name, fasta_target_gaps_removed, sites_with_gaps, codons_without_gaps, codons_with_gaps):
    # Identify position of sites within alignment (same as target species ORF as gaps within target species have been removed)
    ORF_nucleotides_pos = range(0, len(fasta_target_gaps_removed[args.target_species][0]))
    triplets = [ORF_nucleotides_pos[i:i +3] for i in range(0, len(ORF_nucleotides_pos), 3)]
    for codon, sites in enumerate(triplets):
        # Only keep codon if no sites are in gaps
        if len([i for i in sites if i in sites_with_gaps]) == 0:
            # Make human counting of codons
            codon_in_ORF = codon+1
            codons_without_gaps[name].append(codon_in_ORF)
        else:
            # Make human counting of codons
            codon_in_ORF = codon+1
            codons_with_gaps[name].append(codon_in_ORF)
    return codons_without_gaps,codons_with_gaps

#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():
    # Need to account for gaps within target species
    # target species ---TTTAGG---ATT
    # other species  AGT---AGTAGT---
    # Need to ignore gaps at beginning of target species - these are not in ORF
    # Need to ignore gaps within target species - these are not in ORF
    files = 0
    codons_without_gaps = defaultdict(list)
    codons_with_gaps = defaultdict(list)
    
    # Loop through alignments
    for orthogroup in list_folder(args.infolder):
        files +=1
        for aln in list_files(orthogroup):
            name = os.path.basename(aln).split(".stripped.prank.best.fas")[0]
            fasta = read_fasta(aln)
            
            # Remove gaps in target species from the alignment
            # Now can correctly identify position of codons to be removed from target species ORF due to gaps in the rest of the alignment
            # Example filtered alignment:
            # target species TTTAGGATT
            # other species  ---AGT---
            fasta_target_gaps_removed = remove_target_species_gaps(fasta)
            
            # Identify position of gaps in filtered alignment (python counting)
            # Position refers to position of site within filtered alignnment
            sites_with_gaps = identify_sites_with_gaps(fasta_target_gaps_removed)
            
            # Identify position of codons in target species ORF (human counting)
            codons_without_gaps, codons_with_gaps = identify_codons_gaps(name, fasta_target_gaps_removed, sites_with_gaps, codons_without_gaps, codons_with_gaps)

    print "Number of files processed =", files
    print "Number of ORF with at least one codon where all sites are not in gaps =", len(codons_without_gaps)
    print "Number of ORF with at least one codon in a gap =", len(codons_with_gaps)
    
    with open(args.outfile+"codons_without_gaps.txt","w") as outfile_without:
        for gene in codons_without_gaps:
            outfile_without.write(gene)
            outfile_without.write("\t")
            outfile_without.write(str(codons_without_gaps[gene]))
            outfile_without.write("\n")
    
    with open(args.outfile+"codons_with_gaps.txt","w") as outfile_with:
        for gene in codons_with_gaps:
            outfile_with.write(gene)
            outfile_with.write("\t")
            outfile_with.write(str(codons_with_gaps[gene]))
            outfile_with.write("\n")

if __name__ == '__main__':
    main()
