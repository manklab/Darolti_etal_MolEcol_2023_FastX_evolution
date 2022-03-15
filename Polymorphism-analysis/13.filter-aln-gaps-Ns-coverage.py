#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
''' Filters for codons passing gap, ambiguity and coverage filters
Takes a file of genes for which divergence data will be compared to polymorphism data.
Takes folder of folders of sequences aligned by PRANK before gap removal and masking. Identifies codons in ORF which pass both
gap, ambiguity and coverage filters. Filters the PRANK alignment on the basis of these codons. Removes filtered alignments
with a length less than threshold specified. Creates a folder for each orthogroup and writes a new file with the filtered PRANK alignment.'''
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
                    "Same folder that PAML would be run on")
parser.add_argument("target_species", type=str,
                    help="Species for which divergence data will be compared to polymorphism data")
parser.add_argument("codons_without_gaps", type=str,
                    help="A file containing genes with codons where no site is in an alignment gap")
parser.add_argument("codons_without_Ns", type=str,
                    help="A file containing genes with codons where no site is an N")
parser.add_argument("codons_with_coverage", type=str,
                    help="A file containing genes with codons where no site fails the coverage threshold")
parser.add_argument("coordinates", type=str,
                    help="The output file of 06.orf-position-within-transcript.py")
parser.add_argument("paml_branch_lengths", type=str,
                    help="A file of genes for which divergence data will be compared to polymorphism data."
                    "Output of Estimate-divergence-rates/15.paml-extract-branch-lengths-dSSdSfilter.py"
                    "Header=Gene,t,N,S,dN/dS,dN,dS,N*dN,S*dS")
parser.add_argument("outfolder", type=str,
                    help="A folder containing orthogroup infolders")
parser.add_argument("-cutoff", type=str,
                    help="bp")
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
            name = os.path.dirname(aln).split("/")[-1].split("_")[1]
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

def extract_genes_for_paml(source):
    genes_for_paml = []
    with open(source, "r") as infile:
        for line in infile:
            line = line.rstrip()
            if not line.startswith("Gene"):
                ortholog = line.split(",")[0]
                genes_for_paml.append(ortholog)
    return genes_for_paml

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

    # Extract codons with coverage
    # Position refers to ORF of target species
    codons_with_coverage = get_codons(args.codons_with_coverage)  
    print "Number of ORF with codons with coverage =", len(codons_with_coverage)  

    # Make dictionary between orthogroup names and target species' gene names (this is needed as the gene names in the codons_with_coverage file only are taget species' names)
    codons_with_coverage_ortho_gene_name = get_gene_names(args.coordinates, codons_with_coverage)
    print "Number of ORF with codons with coverage ortho names =", len(codons_with_coverage_ortho_gene_name)

    # Extract genes to run PAML on for MK test
    genes_for_paml = extract_genes_for_paml(args.paml_branch_lengths)
    print "Number of genes with PAML data which need alignments filtering =", len(genes_for_paml)   

    # Filter alignment to keep only codons passing both filters
    count_same_aln = 0
    count_genes = 0
    count_passed = 0
    count_failed = 0
    for gene in genes_for_paml:

        # Find codons passing gap, N and coverage filters
        count_genes += 1
        no_gaps = codons_without_gaps[gene]
        no_Ns = codons_without_Ns[gene]
        coverage = codons_with_coverage_ortho_gene_name[gene]
        no_gaps_Ns_and_coverage = list(set(no_gaps).intersection(no_Ns).intersection(coverage))
        
        # First remove gaps within target ORF. These are not listed in no_gaps_and_Ns as they are not in ORF
        fasta_target_gaps_removed = remove_target_species_gaps(aln_file_dict_prefiltering[gene])
        
        # Keep only codons passing all filters
        fasta_no_gaps_Ns_and_coverage = {}
        for s in fasta_target_gaps_removed:
            fasta_no_gaps_Ns_and_coverage[s] = [""]
            # Identify sites in sequence and triplets 
            triplets = [fasta_target_gaps_removed[s][0][i:i +3] for i in range(0, len(fasta_target_gaps_removed[s][0]), 3)]
            for i, fs in enumerate(triplets):
                i = i+1
                if i in no_gaps_Ns_and_coverage:
                    fasta_no_gaps_Ns_and_coverage[s][0] += fs
            
        # Check if length passes cutoff
        if float(len(fasta_no_gaps_Ns_and_coverage[args.target_species][0])) < int(args.cutoff):
            # print gene, float(len(fasta_no_gaps_Ns_and_coverage[args.target_species][0]))
            count_failed += 1
        else:
            count_passed += 1
            # Make a new folder for orthogroup
            folder = args.outfolder+"/"+gene
            outfile_path = folder+"/"+gene+".prank.gaps.Ns.coverage.filtered.fa"
            os.makedirs(folder)
            # Print filtered alignment to new file
            with open(outfile_path,"w") as outfile:
                for s in fasta_no_gaps_Ns_and_coverage:
                    header = ">"+s+"\n"
                    outfile.write(header)
                    outfile.write(fasta_no_gaps_Ns_and_coverage[s][0])
                    outfile.write("\n")
    
    print "Number of genes processed =", count_genes
    print "Number of genes that post filtering pass the length threshold =", count_passed
    print "Number of genes that post filtering failed the length threshold =", count_failed

if __name__ == '__main__':
    main()
