#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

''' Extract Gene Position on Chromosomes
This script uses a GTF file in StringTie format to identify the chromosomal position of genes that remain after filtering (genes in dSSdSfilter.txt file)'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("orthocluster", type=str,
                    help="A txt file containging validated and mapped orthogroups. Output of 03.ortho-cluster.py")
parser.add_argument("outgroupid", type=str,
                    help="Name of species with Ensembl proteome. Must be exactly as in the orthocluster file")
parser.add_argument("speciesname", type=str,
                    help="Name of target species. Must be exactly as in the orthocluster file")
parser.add_argument("gtf", type=str,
                    help="A gtf file for the target species")
parser.add_argument("dSSdS", type=str,
                    help="File with branch lengths information for filtered genes. Output file of 15.paml-extract-branch-lengths-dSSdSfilter.py")
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
#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():

    # The dSSdSfilter file contains only the outgroup gene names
    # We need to correlate the outgroup gene name with the orthologous target species gene name
    dicty_orthologs = {}
    with open(args.orthocluster, "r") as orthocluster:
        for line in orthocluster:
            line = line.split("'")
            sp1 = line[1]
            sp2 = line[3]
            sp3 = line[5]
            sp4 = line[7]

            if sp1.startswith(args.outgroupid):
                outgroupsp = sp1.split("_")[1]
            elif sp2.startswith(args.outgroupid):
                outgroupsp = sp2.split("_")[1]
            elif sp3.startswith(args.outgroupid):
                outgroupsp = sp3.split("_")[1]
            elif sp4.startswith(args.outgroupid):
                outgroupsp = sp4.split("_")[1]

            if sp1.startswith(args.speciesname):
                targetsp = sp1.split("_")[1]
            elif sp2.startswith(args.speciesname):
                targetsp = sp2.split("_")[1]
            elif sp3.startswith(args.speciesname):
                targetsp = sp3.split("_")[1]
            elif sp4.startswith(args.speciesname):
                targetsp = sp4.split("_")[1]

            dicty_orthologs[outgroupsp] = targetsp

    # Extract gene position in chromosomes based on a gtf file in StringTie output format.
    dicty_position = defaultdict(list)
    with open(args.gtf, "r") as gtf:
        for line in gtf:
            if not line.startswith("#"):
                if line.split()[2] == "transcript": # If Ensembl format is used instead, replace "transcript" with "gene"
                    gene_id = line.split('transcript_id "')[1].split('";')[0] # If Ensembl format is used instead, replace "transcript_id" with "gene_id"
                    chromosome = line.split()[0]
                    start_position = line.split()[3]
                    dicty_position[gene_id].append(chromosome)
                    dicty_position[gene_id].append(start_position)

    # Obtain positional information for genes in the dSSdSfilter file
    count_dSSdS = -1
    count_dSSdS_position = 0
    with open(args.outfile, "w") as out:
        with open(args.dSSdS, "r") as dSSdS:
            for line in dSSdS:
                line = line.rstrip()
                if count_dSSdS == -1 :
                    out.write(line)
                    out.write(",TargetGene,LG,StartPos\n")
                else:
                    outgroupgene = line.split(",")[0]
                    targetgene = dicty_orthologs[outgroupgene]
                    if targetgene in dicty_position:
                        count_dSSdS_position += 1
                        out.write(line)
                        out.write(",")
                        out.write(targetgene)
                        out.write(",")
                        out.write(dicty_position[targetgene][0])
                        out.write(",")
                        out.write(dicty_position[targetgene][1])
                        out.write("\n")
                count_dSSdS += 1

    print "Number of genes with orthogroup = ", len(dicty_orthologs)
    print "Number of genes with positional information = ", len(dicty_position)
    print "Number of genes with dNdS values = ", count_dSSdS
    print "Number of genes with dNdS values and positional information = ", count_dSSdS_position
    
if __name__ == '__main__':
    main()
