#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
''' STAR map second pass
Takes a folder containing paired-end fastq files and finds pairs of forward and reverse reads.
Takes file of merged splice junction files generated from first pass of STAR.
Executes STAR mapping commands one after the other until finished.'''
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
                    help="Infolder of fastq files on which to run STAR")
parser.add_argument("STARindex", type=str,
                    help="Path to folder containing STAR index (same as specified in --genomeDir in index step)")
parser.add_argument("STARmergedsplicejunctions", type=str,
                    help="Path to file of merged SJ.out.tab files from STAR first pass")
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
    return [os.path.join(infolder, f) for f in os.listdir(infolder) if f.endswith(".gz")]

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
    infiles = list_folder(args.infolder)
    print "Number of infiles (fastq) =", len(infiles)

    #Creates dictionary of forward and reverse pairs for each sample
    filedictionary = defaultdict(list)
    for infile in infiles:
        name = os.path.basename(infile)
        uniqid = name.split("_")[1] + "_" + name.split("_")[2] #Gets unique identifier for each sample. Assumes name of fastq files is of the format sample_female_1_R1.fastq.gz or sample_female_2_R2.fastq.gz
        directory1 = os.path.dirname(infile) + "/secondpass/"
        directory2 = os.path.dirname(infile) + "/secondpass/" + uniqid
        if os.path.exists(directory1) == False:
            os.mkdir(directory1)
        if os.path.exists(directory2) == False:
            os.mkdir(directory2)
        filedictionary[uniqid].append(infile)
    print "Number of samples =", len(filedictionary)

    #Creates list of STAR commands to run
    commandstorun = []
    for uniqid in filedictionary:
        for file in filedictionary[uniqid]:
            initial_directory = os.getcwd()
            directory = os.path.dirname(file)+"/secondpass/"+uniqid+"/"
            if file.endswith("R1.fastq.gz"):
                forward_fastq = file
            elif file.endswith("R2.fastq.gz"):
                reverse_fastq = file
        command = ["STAR --genomeDir "+args.STARindex+" --readFilesCommand zcat --readFilesIn "+forward_fastq+" "+reverse_fastq+" --runThreadN 12 --outFilterMultimapNmax 1 --outFileNamePrefix "+directory+" --sjdbFileChrStartEnd "+args.STARmergedsplicejunctions]
        commandstorun.append(command)

    print "Number of commands to run =", len(commandstorun)
    print commandstorun
    print "\n"
    exec_in_row(commandstorun)
    print "\n"

if __name__ == "__main__":
    main()