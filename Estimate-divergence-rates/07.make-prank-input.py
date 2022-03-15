#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

''' Prepare PRANK
This script prepares sequences for prank alignment of orthologous sequences.
Takes each folder of orthologous sequences and removes all files apart from cdna
fasta file. Makes a new file with species names as the sequence id. Necessary for
consistency with prank tree.'''
#==============================================================================
import argparse
import sys
import os
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="An input folder containing folders of sequences of orthogroups")
parser.add_argument("uniq", type=str,
                    help="Unique identifier which file to align ends with")
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
                f = os.path.join(path, name)
                file_list.append(f)
    return file_list

def read_fasta(source):
    ''' Read Ensembl fasta file and outputs a dictionary 
    key [seq name] values [sequence] [header]
    '''
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
#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():
    count = 0 
    for orthogroup in list_folder(args.infolder):
        for f in list_files(orthogroup):
            #remove everything except alignment file
            if not f.endswith(args.uniq):
                print f
                os.remove(f)
            else:
                #makes new files with species names as seq id
                seq = read_fasta(f)
                with open(f,"wb") as outfile:
                    for s in seq:
                        count += 1
                        title = ">"+s.split("_")[0]+"\n"
                        outfile.write(title)
                        outfile.write(seq[s][0])
                        outfile.write("\n")
    print "Number of files modified =", count/len(seq)

if __name__ == '__main__':
    main()
