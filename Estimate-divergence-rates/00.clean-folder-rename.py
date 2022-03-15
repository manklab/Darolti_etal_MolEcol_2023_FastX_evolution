#!/usr/bin/python2.6=7
# -*- coding: utf-8 -*-

''' Clean folder. Renames files in a folder that end in string specified.'''
#==============================================================================
import argparse
import sys
import os
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder containing folders of orthogroups")
parser.add_argument("uniqrename", type=str,
                    help="Unique identifier which file to rename ends with")
parser.add_argument("uniqnewname", type=str,
                    help="New ending of file to convert")
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
            f = os.path.join(path, name)
            file_list.append(f)
    return file_list
#==============================================================================
#Main==========================================================================
#==============================================================================m
def main():
    print "Renaming ", args.uniqrename, " to ", args.uniqnewname

    count = 0 
    for orthogroup in list_folder(args.infolder):
        count += 1
        for f in list_files(orthogroup):
            if f.endswith(args.uniqrename):
                name = f
                newname = f.split(args.uniqrename)[0]+args.uniqnewname
                os.rename(name, newname)
    print "Number of folders modified =", count

if __name__ == '__main__':
    main()