#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
''' Extract Gene Position on Chromosomes'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("gtf", type=str,
                    help="A gtf file for the target species")
parser.add_argument("infile", type=str,
                    help="File with genes containing both divergence and polymorphismm data. Output file of 17.filter-divergence-and-polymorphism-genes.py")
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

    # Obtain positional information for genes in the infile
    count_dSSdS = -1
    count_dSSdS_position = 0
    with open(args.outfile, "w") as out:
        with open(args.infile, "r") as infile:
            for line in infile:
                line = line.rstrip()
                if count_dSSdS == -1 :
                    out.write(line)
                    out.write(",LG,StartPos\n")
                else:
                    targetgene = line.split(",")[0]
                    if targetgene in dicty_position:
                        count_dSSdS_position += 1
                        out.write(line)
                        out.write(",")
                        out.write(dicty_position[targetgene][0])
                        out.write(",")
                        out.write(dicty_position[targetgene][1])
                        out.write("\n")
                count_dSSdS += 1
                
    print "Number of genes with positional information = ", len(dicty_position)
    print "Number of genes with dNdS values = ", count_dSSdS
    print "Number of genes with dNdS values and positional information = ", count_dSSdS_position
    
if __name__ == '__main__':
    main()
