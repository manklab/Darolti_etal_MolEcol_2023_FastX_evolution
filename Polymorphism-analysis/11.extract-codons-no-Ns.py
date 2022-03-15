#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
''' Identifies codons where no sites have ambiguity data
Takes folder of folders of sequences aligned by PRANK after gap removal and masking with SWAMP but before N removal. 
Assumes that gaps have been removed before masking. Accounts for position of gaps when identifying the position of codons to remove due to ambiguity data. 
Reports codons without Ns within ORF. Position reported in relation to codons within ORF. i.e. 1 = 1st codon in ORF/protein.'''
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
                    help="An input folder containing folders of phylip sequences of orthogroups aligned after gap removal and masking but before N removal")
parser.add_argument("target_species", type=str,
                    help="Species for which divergence data will be compared to polymorphism data")
parser.add_argument("species", type=str,
                    help="List of species names e.g., Poeciliareticulata,Poeciliaformosa,Xiphophorusmaculatus,Oryziaslatipes")
parser.add_argument("codons_with_gaps", type=str,
                    help="A file containing genes with codons within alignment gaps")
parser.add_argument("outfile", type=str,
                    help="An outfile file of codons passing ambiguity data filtering")
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
                if name.endswith("_masked.phy"):
                    f = os.path.join(path, name)
                    file_list.append(f)
    return file_list

def read_phylip(source, species):
    ''' Read phylip file and outputs a dictionary 
    key [seq name] values [sequence] [header]'''
    species = species.split(",")
    SG = {}
    try:
        with open(source, "r") as file:
            for line in file.readlines():
                line = line.rstrip()
                if not line.startswith(" "):
                    if line in species:
                        name = line
                        header = line
                        SG[name] = ["", header]
                    else:
                        SG[name][0] += line
            return SG
    except IOError:
        print "File does not exit!"

def get_gaps(file):
    # Get dictionary of codons with gaps in ORF
    codons_with_gaps = defaultdict(list)
    with open(args.codons_with_gaps, "r") as infile:
        for line in infile:
            line = line.rstrip().split("\t")
            name = line[0]
            codons_gaps = line[1].split("[")[1].split("]")[0].split(",")
            for codon in codons_gaps:
                codon = float(codon)
                codons_with_gaps[name].append(codon)
    print "Number of ORF with codons in gaps =", len(codons_with_gaps)
    return codons_with_gaps

def add_gaps(fasta, codons_with_gaps):
    fasta_gaps_added = defaultdict(list)
    length_ORF = len(codons_with_gaps) + float(len(fasta[args.target_species][0]))/3
    for s in fasta:
        fasta_gaps_added[s] = [""]
        ORF_nucleotides = fasta[s][0]
        triplets = [ORF_nucleotides[i:i +3] for i in range(0, len(ORF_nucleotides), 3)]
        counter = 1
        offset = 0
        while counter <= length_ORF:
            if counter in codons_with_gaps:
                fasta_gaps_added[s][0] += "---"
                offset += 1
            else:
                triplet_pos = int(counter-1-offset)
                fasta_gaps_added[s][0] += triplets[triplet_pos]
            counter += 1
    return fasta_gaps_added

def identify_sites_with_Ns(fasta):
    sites_with_Ns = {}
    for s in fasta:
        for i, p in enumerate(fasta[s][0]):
            if p == "N":
                sites_with_Ns[i] = 0
    return sites_with_Ns

def identify_codons_Ns(name, fasta_gaps_added, sites_with_Ns, codons_without_Ns):
    # Identify position of sites within alignment (same as target species ORF as gaps in alignment have been added back in)
    ORF_nucleotides_pos = range(0, len(fasta_gaps_added[args.target_species][0]))
    triplets = [ORF_nucleotides_pos[i:i +3] for i in range(0, len(ORF_nucleotides_pos), 3)]
    for codon, sites in enumerate(triplets):
        # Only keep codon if no sites are Ns
        if len([i for i in sites if i in sites_with_Ns]) == 0:
            # Make human counting of codons
            codon_in_ORF = codon+1
            codons_without_Ns[name].append(codon_in_ORF)
    return codons_without_Ns

#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():

    files = 0
    codons_without_Ns = defaultdict(list)

    # Extract codons in gaps
    codons_with_gaps = get_gaps(args.codons_with_gaps)
    
    # Loop through alignments
    for orthogroup in list_folder(args.infolder):
        files +=1
        for aln in list_files(orthogroup):
            name = os.path.basename(aln).split("_masked.phy")[0]
            fasta = read_phylip(aln, args.species)

            # Add gaps back into masked alignment
            # But not gaps before or within target species ORF - not in codons_with_gaps
            if len(codons_with_gaps[name]) == 0:
                fasta_gaps_added = fasta
            else:
                fasta_gaps_added = add_gaps(fasta, codons_with_gaps[name])

            # Identify position of Ns in alignment (python counting)
            # Position refers to position of site within alignnment where gaps have been added in
            sites_with_Ns = identify_sites_with_Ns(fasta_gaps_added)

            # Identify position of codons in ORF where all sites are not N (human counting)
            codons_without_Ns = identify_codons_Ns(name, fasta_gaps_added, sites_with_Ns, codons_without_Ns)

    print "Number of files processed =", files
    print "Number of ORF with at least one codon where all sites are not Ns =", len(codons_without_Ns)

    with open(args.outfile,"w") as outfile:
        for gene in codons_without_Ns:
            outfile.write(gene)
            outfile.write("\t")
            outfile.write(str(codons_without_Ns[gene]))
            outfile.write("\n")

if __name__ == '__main__':
    main()
