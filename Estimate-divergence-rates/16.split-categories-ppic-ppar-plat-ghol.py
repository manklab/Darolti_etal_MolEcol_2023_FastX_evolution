#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

''' Split genes into different categories based on their chromosomal position'''
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
                    help="A dSSdSfilter file with positional information of genes into chromosomes")
parser.add_argument("outfile_autosomes", type=str,
                    help="A file containing all autosomal genes")
parser.add_argument("outfile_sexchromosomes_all", type=str,
                    help="A file containing all X-linked genes")
parser.add_argument("outfile_sexchromosomes_noPAR", type=str,
                    help="A file containing X-linked genes excluding those in the pseudoautosomal region")
parser.add_argument("PAR", type=str,
                    help="A file of genes in the pseudoautosomal region")
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

    count_all = -1
    count_autosomes = 0
    count_sexchromosomes_all = 0
    count_sexchromosomes_noPAR = 0
    count_PAR = 0

    with open(args.PAR, "w") as PAR:
        with open(args.outfile_sexchromosomes_noPAR, "w") as sex_noPAR:
            with open(args.outfile_autosomes, "w") as auto:
                with open(args.outfile_sexchromosomes_all, "w") as sex_all:
                    with open(args.infile, "r") as infile:
                        for line in infile:
                            if count_all == -1:
                                auto.write(line)
                                sex_all.write(line)
                                sex_noPAR.write(line)
                                PAR.write(line)
                            else:
                                chromo = line.split(",")[-2]
                                if chromo == "chr8":
                                    sex_all.write(line)
                                    count_sexchromosomes_all += 1
                                    start = float(line.split(",")[-1])
                                    if start >= 20000000.0:
                                        PAR.write(line)
                                        count_PAR += 1
                                    else:
                                        sex_noPAR.write(line)
                                        count_sexchromosomes_noPAR += 1
                                else:
                                    auto.write(line)
                                    count_autosomes += 1
                            count_all += 1

    print "Total number of genes = ", count_all
    print "Number of autosomal genes = ", count_autosomes
    print "Number of sex chromosome genes (all) = ", count_sexchromosomes_all 
    print "Number of sex chromosome genes (no PAR) = ", count_sexchromosomes_noPAR
    print "Number of PAR genes = ", count_PAR
   
if __name__ == '__main__':
    main()
