#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
''' CONCATENATE SPLICE JUNCTIONS
Takes a folder containing folders of STAR first pass mapping. One folder for each sample.
Identifies SJ.out.tab files and concatenates them together into a new folder'''
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
                    help="Infolder of folders containing SJ.out.tab files to concatenate")
parser.add_argument("outfile", type=str,
                    help="An outfile of merged SJ.out.tab files")
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
    ''' Returns the number of CPUs in the system
    '''
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
        print task
        p = Popen(task, shell=True)
        p.wait()

    if done(p):
            if success(p):
                print "done!"
            else:
                fail()

def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder)]

def list_files(current_dir):
    file_list = []
    for path, subdirs, files in os.walk(current_dir): # Walk directory tree
        for name in files:
            if name.endswith("SJ.out.tab"):
                f = os.path.join(path, name)
                file_list.append(f)
    return file_list
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    infiles = []
    infolders = list_folder(args.infolder)
    for infolder in infolders:
        files = list_files(infolder)
        if len(files) == 1:
            splicefile = files[0]
            infiles.append(splicefile)
    print "Number of infiles (SJ.out.tab) =", len(infiles)
    print infiles

    #make concatenate commands
    makeconcatenaterun = []
    infilenames = infiles[0]
    for file in infiles[1:]:
        infilenames = infilenames+" "+file
        
    concatenatecommand = ["cat "+str(infilenames)+" > "+args.outfile]
    makeconcatenaterun.append(concatenatecommand)

    #run concatenate
    print "Number of commands to run =", len(makeconcatenaterun)
    print makeconcatenaterun
    # exec_in_row(makeconcatenaterun)

if __name__ == "__main__":
    main()
