#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

''' Correct for the X inversion present in the P. reticulata reference geneome '''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
from collections import OrderedDict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="dSSdSfilter file with positional information of genes on chromosomes")
parser.add_argument("outfile", type=str,
                    help="dSSdSfilter file with positional information corrected for the inversion")
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
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    with open(args.outfile, "w") as outfile:
        with open(args.infile, "r") as infile:
            for line in infile:
                line = line.rstrip()
                if line.startswith("Gene"):
                    outfile.write(line)
                    outfile.write(",StartPos_inversion\n")
                else:
                    chromosome = line.split(",")[-2]
                    if chromosome == "chr12":
                        start_pos = float(line.split(",")[-1])
                        if start_pos < 40654:
                            new_start_pos = start_pos
                        if start_pos >= 40654 and start_pos < 9942741:
                            new_start_pos = 9942741 - start_pos + 10805812
                        if start_pos >= 9942741 and start_pos < 20707899:
                            new_start_pos = start_pos - 9942741 + 40654
                        if start_pos >= 20707899:
                            new_start_pos = start_pos
                        outfile.write(line)
                        outfile.write(",")
                        outfile.write(str(new_start_pos))
                        outfile.write("\n")
                    else:
                        start_pos = float(line.split(",")[-1])
                        outfile.write(line)
                        outfile.write(",")
                        outfile.write(str(start_pos))
                        outfile.write("\n")

if __name__ == '__main__':
    main()