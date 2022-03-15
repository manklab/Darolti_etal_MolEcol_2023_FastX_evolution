#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-

''' Removes Ns in alignment
This script removes masked sites (Ns) in the phylip alignment of each orthogroup. It removes alignments
where the length after Ns removal fails the cutoff threshold specified.'''
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
parser.add_argument("invalid_folder", type=str,
                    help="An output folder for invalid orthogroups")
parser.add_argument("species", type=str,
                    help="List of species names eg Poeciliareticulata,Poeciliaformosa,Xiphophorusmaculatus,Oryziaslatipes")
parser.add_argument("-cutoff", type=str,
                    help="bp") #use cut off from previous script
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
    ''' Read Ensembl phylip file and outputs a dictionary 
    key [seq name] values [sequence] [header]
    '''
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

#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():
    ammend = 0
    invalid_orthogroups = []
    for orthogroup in list_folder(args.infolder):
        for aln in list_files(orthogroup):
            ammend = ammend+1
            final_remove = {}
            phylip = read_phylip(aln, args.species)
            for s in phylip:
                for i, p in enumerate(phylip[s][0]):
                    if p == "N":
                        final_remove[i] = 0
            original_length = float(len(phylip[s][0]))
            final_remove_length  = float(len(final_remove))
            if original_length - final_remove_length < int(args.cutoff):
                invalid_orthogroups.append(orthogroup)
            else:
                outfilename = aln[:-4]+".Nrm.fa"
                with open(outfilename,"wb") as outfile:
                    for s in phylip:
                        seq = []
                        outfile.write(">"+s+"\n")
                        for i, fs in enumerate(phylip[s][0]):
                            if i not in final_remove:
                                seq.append(fs)
                                outfile.write(fs)
                        outfile.write("\n")

    print "%i orthogroups are invalid and were removed!" % len(set(invalid_orthogroups))
    print "Total number of orthogroups processed =", ammend
    if not os.path.isdir(args.invalid_folder): os.mkdir(args.invalid_folder)
    for o in list(set(invalid_orthogroups)):
        import shutil
        shutil.move(o, args.invalid_folder)
    
if __name__ == '__main__':
    main()
