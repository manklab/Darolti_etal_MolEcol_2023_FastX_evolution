#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
''' Runs Samtools Mpileup
Takes a folder containing folders of sorted bam files and an indexed genome to
run samtools mpileup'''
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
                    help="Infolder of files on which to run Samtools")
parser.add_argument("reference", type=str,
                    help="Index reference file")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder)]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if name.endswith("_sorted"):
                f = os.path.join(path, name)
                file_list.append(f)
    return file_list

def exec_in_row(cmds):
    ''' Execute commands one after the other until finished.'''
    if not cmds:
        return  # empty list
    def done(p):
        return p.poll() is not None
    def success(p):
        return p.returncode == 0
    def fail():
        sys.exit(1)
    for task in cmds:
        print task
        p = Popen(task, shell=True)
        p.wait()
    if done(p):
            if success(p):
                print "done!"
            else:
                fail()
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    infiles = []
    infolders = list_folder(args.infolder)
    for infolder in infolders:
        files = list_files(infolder)
        if len(files) == 1:
            bamfile = files[0]
            infiles.append(bamfile)
    print "Number of infiles =", len(infiles)

    directory = os.path.dirname(args.infolder) + "/mpileup/"
    if os.path.exists(directory) == False:
        os.mkdir(directory)

    #Make list of bam files
    bamlist = infiles[0]
    for file in infiles[1:]:
        bamlist = bamlist+" "+file
    print bamlist

    outputfile = directory+"output.pileup"
    commandstorun = []
    command = ["samtools mpileup -Q 20 -Bd 10000000 -f "+args.reference+" "+bamlist+" > "+ outputfile]
    commandstorun.append(command)

    print "Number of commands to run =", len(commandstorun)
    print commandstorun
    print "\n"
    exec_in_row(commandstorun)
    print "\n"

if __name__ == "__main__":
    main()
