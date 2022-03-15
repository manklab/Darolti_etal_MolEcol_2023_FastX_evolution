#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
''' Get coordinates of open reading frame within transcripts.
This script should be run on output of Evolutionary-divergence-rates/06.orf-prediction.py. 
It takes the output of blastx for each orthogroup and identifies the start and stop position of the open reading frame within the transcript.
Need blast to be run with format '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe sseq'. '''
#==============================================================================
import argparse
import sys
import os
from Bio.Seq import Seq
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infolder", type=str,
                    help="An input folder containing folders of sequences of orthogroups with open reading frames"
                    "Output folder from alielw/Evolutionary-rates-analysis/06.orf-prediction.py")
parser.add_argument("species", type=str,
                    help="Species for which divergence data will be compared to polymorphism data"
                    "eg Poeciliareticulata")
parser.add_argument("outfile", type=str,
                    help="Outfile containing coordinates of open reading frame within the transcript")
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
    ''' Read fasta file and outputs a dictionary 
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
        print "File does not exist!"

def read_blastout_bitscore(source):
    '''Reads in a blastout file, picks tophit and stores the results in a dictionary
    key [query id] values [seq id, bitscore, pidentiy].
    The blastout output format supported is: "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sseq"
    For each hit only the besthit is stored, rest is discarded. Minimum 30% identity, then pick tophit with greatest bitscore. 
    If bitscores are identical pick the tophit with greatest %identity. If equal %identity then tophit randomly picked.
    '''
    try:
        with open(source, "r") as infile:
           blasthits = {}
           sbase = os.path.basename(source).split(".")
           s_db = sbase[0].split("_")[0]
           s_query = sbase[2].split("_")[0]
           for line in infile:
                line = line.rstrip().split()
                qseqid = line[0]
                sseqid = line[1]
                pidentity = float(line[2])
                evalue = float(line[10])
                bitscore =float(line[11])
                qframe = int(line[12])
                sframe = int(line[13])
                qstart = int(line[6])
                qend = int(line[7])
                if pidentity > 30: # Check for min 30% identity
                    if qseqid in blasthits:
                        if blasthits[qseqid][1] < bitscore:
                            blasthits[qseqid] = (qseqid,
                                            bitscore,
                                            sframe,
                                            qframe,
                                            qstart,
                                            qend)
                        elif blasthits[qseqid][1] == bitscore:
                            if blasthits[qseqid][2] < pidentity:
                                blasthits[qseqid] = (qseqid,
                                            bitscore,
                                            sframe,
                                            qframe,
                                            qstart,
                                            qend)
                    else:
                        blasthits[qseqid] = (qseqid,
                                            bitscore,
                                            sframe,
                                            qframe,
                                            qstart,
                                            qend)
        return blasthits
    except IOError:
        print "!----ERROR----!"
        print "File %s does not exit!" % source
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(1)

def assign_files(infolder):
    for f in list_files(infolder):
        if f.endswith("QUERY.fasta"):
            query_seq = read_fasta(f)
        elif f.endswith(".bla"):
            blastout = read_blastout_bitscore(f)
        abs_path = os.path.dirname(f)
    return query_seq, blastout, abs_path

def get_orf_position(qstart, qend, frame, seq_len):
    ''' Takes the query start of the alignment and query end of the alignment
    and calculates the temporary open reading frame positions in the the cDNA.
    For negative frames, the reverse complement is needed and the positions
    are substracted from the length of the input sequence. For positive frames
    no reversal is needed. Orfstop in reverse complement +1 for divisible by
    3, same in orfstart non reverse complement'''
    if frame < 0:
        #account for reverse complement
        orfstart = seq_len - qstart
        orfstop  = seq_len - qend + 1 
    else:
        #account for python 0 start and slice notation
        orfstart = qstart-1
        orfstop = qend
    return orfstart, orfstop

def get_offset(orfstart):
    '''Dynamically calculates the offset for a given frame from orfstart.'''
    ts = orfstart
    while ts >3:
        ts -= 3
    return ts

def translate_orf_positions(orfstart, orfstop, offset, frame, seq):
    ''' Takes orfstart and orfstop (calculated by get_orf_position) and the
        offset (calculated by get_offset) plus the frame and sequence
        to calculate the start and stop of the BLASTX alignment in the
        translated protein and translates the cDNA in the right frame.
        At the moment we think protstop is +1, not the real one'''
    if frame < 0:
        cdna = str(Seq(seq).reverse_complement())
    else:
        cdna = seq
    if (orfstart - offset) % 3 == 0:
        protstart = (orfstart - offset) / 3
    else:
        print "orfstart - offset is not divisible by 3!"
        print orfstart, offset
    if (orfstop - offset) % 3 == 0:
        protstop = (orfstop - offset) / 3
    else:
        print "orfstop - offset is not divisible by 3!"
        print orfstop, offset
    tmp_protein = str(Seq(cdna[int(offset):]).translate())
    return protstart, protstop, tmp_protein

def extend_forward(protein, protstart):
    # reduce protstop by 1 to start at the right AA
    # offset = len(protein) - len(protein[protstart:])
    offset = len(protein) - len(protein[int(protstart):])
    for i, p in enumerate(protein[int(protstart):]):
        if p == "*":
            return i+offset
        elif i == len(protein[int(protstart):])-1:
            return i+offset

def find_character_index(seq, char):
    return [i for i, c in enumerate(seq) if c == char]

def extend_backwards(protein, protstop):
    # reduce protstop by 1 to start at the right AA
    overhead = protein[:protstop]
    stop_pos = find_character_index(overhead, "*")
    start_pos = find_character_index(overhead, "M")
    if not stop_pos and not start_pos: 
        return 0
    elif not stop_pos and start_pos:
        return start_pos[0]
    elif not start_pos and stop_pos:
        return None
    elif start_pos and stop_pos:
        max_stop = max(stop_pos)
        max_start = max(start_pos)
        dist = None
        start = None
        for i in start_pos:
            if i > max_stop:
                if dist:
                    if (i - max_stop) < dist:
                        start = i
                        dist = i - max_stop
                else:
                    start = i
                    dist = i - max_stop
        return start
    else:
        "No condition is applying..aborting!"
        sys.exit(1)

def get_extended_cdna_pos(extended_protstop, 
                          extended_protstart, 
                          frame, 
                          orfstart, 
                          orfstop,
                          seq_len, offset):  
    if frame > 0:
        cdna_start = extended_protstart*3 + offset
        cdna_stop = extended_protstop*3 + offset
    else:
        cdna_start = seq_len - (extended_protstart*3 + offset)
        cdna_stop = seq_len - (extended_protstop*3 + offset)
    return cdna_start, cdna_stop

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

    print "Species =", args.species

    #calculate open reading frame
    with open(args.outfile, "w") as coord_outfile:
        coord_outfile.write("Gene")
        coord_outfile.write("\t")
        coord_outfile.write("Orthogroup")
        coord_outfile.write("\t")
        coord_outfile.write("Sequence")
        coord_outfile.write("\t")
        coord_outfile.write("Start")
        coord_outfile.write("\t")
        coord_outfile.write("Stop")
        coord_outfile.write("\n")
        print "Number of orthogroups =", len(list_folder(args.infolder))
        for orthogroup in list_folder(args.infolder):
            query_seq, blastout, abs_path = assign_files(orthogroup)
            if len(blastout) == len(query_seq):
                ortho_name = abs_path.split("/")[-1].split("_")[1]
                translated_outfile_path = abs_path+"/"+ortho_name+".stripped.check.pep.fa"
                cdna_outfile_path = abs_path+"/"+ortho_name+".stripped.check.cdna.fa"
                species_dict = defaultdict(list)
                with open(translated_outfile_path, "w") as query_outfile, open(cdna_outfile_path, "w") as cdna_outfile:
                    for hit in blastout:
                        seq, bitscore, sframe, qframe, qstart, qend = blastout[hit]
                        seq_len = len(query_seq[seq][0])
                        orfstart, orfstop = get_orf_position(qstart, qend, qframe, seq_len)
                        offset = get_offset(orfstart)
                        protstart, protstop, prot = translate_orf_positions(orfstart, 
                                                                            orfstop, 
                                                                            offset, 
                                                                            qframe, 
                                                                            query_seq[seq][0]) 

                        extended_protstop = extend_forward(prot, protstart)
                        extended_protstart = extend_backwards(prot, extended_protstop)
                        if extended_protstart != None:
                            #python counting
                            cdna_start, cdna_stop = get_extended_cdna_pos( extended_protstop, 
                                                                          extended_protstart, 
                                                                          qframe, 
                                                                          orfstart, 
                                                                          orfstop,
                                                                          seq_len,
                                                                          offset)
                        if extended_protstart != None:
                            extended_prot = prot[extended_protstart:extended_protstop+1]
                            if extended_prot[-1] == "*":
                                    extended_prot = extended_prot[:-1]
                                    if qframe < 0:
                                        first = cdna_stop  
                                        last = cdna_start
                                        extended_cdna = str(Seq(query_seq[seq][0][first:last]).reverse_complement()) 
                                        species_dict[hit].append("reversecomplement")
                                        species_dict[hit].append(first+1) #python counting so need to add 1 for first. not for last due to splice
                                        species_dict[hit].append(last)
                                    else:
                                        first = cdna_start
                                        last = cdna_stop
                                        extended_cdna = query_seq[seq][0][first:last]
                                        species_dict[hit].append("complement")
                                        species_dict[hit].append(first+1) #python counting so need to add 1 for first. not for last due to splice
                                        species_dict[hit].append(last)
                            else:
                                extended_prot = prot[extended_protstart:extended_protstop+1]
                                if qframe < 0:
                                    first = cdna_stop-3  #account for python counting of amino acids
                                    last = cdna_start
                                    extended_cdna = str(Seq(query_seq[seq][0][first:last]).reverse_complement())
                                    species_dict[hit].append("reversecomplement")
                                    species_dict[hit].append(first+1) #python counting so need to add 1 for first. not for last due to splice
                                    species_dict[hit].append(last)
                                else:
                                    first = cdna_start
                                    last = cdna_stop+3  #account for python counting of amino acids
                                    extended_cdna = query_seq[seq][0][first:last]
                                    species_dict[hit].append("complement")
                                    species_dict[hit].append(first+1) #python counting so need to add 1 for first. not for last due to splice
                                    species_dict[hit].append(last)
                            
                            tmp_protein = str(Seq(extended_cdna).translate())
                            if tmp_protein != extended_prot:
                                print "ERROR: Translated cDNA does not match extended Protein!"
                                sys.exit(0)
                            query_outfile.write(">%s %s \n" % (seq, qframe))
                            query_outfile.write(extended_prot+"\n")
                            cdna_outfile.write(">%s %s \n" % (seq, qframe))
                            cdna_outfile.write(extended_cdna+"\n")
                        else:
                            print "ERROR - invalid open reading frame"

                    for species in species_dict:
                        if species.split("_")[0] == args.species:
                            coord_outfile.write(species)
                            coord_outfile.write("\t")
                            coord_outfile.write(os.path.basename(orthogroup))
                            coord_outfile.write("\t")
                            coord_outfile.write(str(species_dict[species][0]))
                            coord_outfile.write("\t")
                            coord_outfile.write(str(species_dict[species][1]))
                            coord_outfile.write("\t")
                            coord_outfile.write(str(species_dict[species][2]))
                            coord_outfile.write("\n")

if __name__ == '__main__':
    main()
