#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

''' Run PRANK
This script automatizes a prank alignment between reciprocal orthologs.'''
#==============================================================================
import argparse
import sys
import os
from subprocess import Popen
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="Infolder containing folders of orthogroups")
parser.add_argument("tree", type=str,
                    help="A file specifying a rooted tree for Prank "
                    "eg File contains (((Poeciliareticulata,Poeciliaformosa),Xiphophorusmaculatus),Oryziaslatipes); "
                    "Species names must match names in fasta file")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def cpu_count():
    ''' Returns the number of CPUs in the system '''
    num = 1
    if sys.platform == 'win32':
        try:
            num = int(os.environ['NUMBER_OF_PROCESSORS'])
        except (ValueError, KeyError):
            pass
    elif sys.platform == 'darwin':
        try:
            num = int(os.popen('sysctl -n hw.ncpu').read())
        except ValueError:
            pass
    else:
        try:
            num = os.sysconf('SC_NPROCESSORS_ONLN')
        except (ValueError, OSError, AttributeError):
            pass

    return num

def exec_in_row(cmds):
    ''' Exec commands one after the other until finished. This is helpful
    if a program is already parallelized, but we have to submit many jobs
    after each other.'''
    if not cmds:
        return  # empty list

    def done(p):
        return p.poll() is not None

    def success(p):
        return p.returncode == 0

    def fail():
        sys.exit(1)

    for task in cmds:
        p = Popen(task)
        p.wait()

    if done(p):
            if success(p):
                print "done!"
            else:
                fail()

def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if not f.endswith(".DS_Store")]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if name.endswith("cdna.fa"):
                f = os.path.join(path, name)
                file_list.append(f)
    return file_list
#==============================================================================
#Main==========================================================================
#==============================================================================

def main():
    '''Main runs only when the file itself is executed.
    If a part of the script is imported, main will not be executed.'''
    makeprankrun = []
    infiles = list_folder(args.infolder)
    # Creates the prank commands for each file supplied.
    for infile in infiles:
        for inf in list_files(infile):
            outfile = inf[:-8]+".prank"
            makeprankrun.append(["prank", "-d="+inf, 
                                               "-f=fasta",
                                               "-t="+args.tree,
                                               "-o="+outfile,
                                               "-DNA",
                                               "-codon",
                                               "-once"
                                               ])
    # Passes the makeprankrun commands (a list of lists) to the exec_commands
    # function and executes them in parallel.
    # exec_commands(makeprankrun)
    # only use in row for Prank as can cause interference if run in parallel
    exec_in_row(makeprankrun)

if __name__ == "__main__":
    main()
