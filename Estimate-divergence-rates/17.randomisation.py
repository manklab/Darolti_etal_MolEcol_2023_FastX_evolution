#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

''' Calculate confidence intervals for dN, dS and dNdS estimates'''
#==============================================================================
import argparse
import sys
import os
import random
import numpy as np
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("dSSdSfilter", type=str,
                    help="An infile containing genes with dSSdSfilter values")
parser.add_argument("out_dN_bootstrap", type=str,
                    help="Bootstrap values for dN estimates")
parser.add_argument("out_dS_bootstrap", type=str,
                    help="Bootstrap values for dS estimates")
parser.add_argument("out_dN_dS_bootstrap", type=str,
                    help="Bootstrap values for dNdS estimates")
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

    bootstrap = 0
    length = 0

    list_lines = []
    dN_dS_list = []
    dN_list = []
    dS_list = []

    count = 0
    with open(args.dSSdSfilter, 'r') as infile:
        for line in infile:
            if count >0:
                length += 1
                list_lines.append(line)
            count += 1
    with open(args.out_dN_bootstrap, 'w') as out_dN:
        with open(args.out_dS_bootstrap, 'w') as out_dS:
            with open(args.out_dN_dS_bootstrap, 'w') as out_dN_dS:
                while bootstrap < 1000:
                    bootstrap += 1
                    bootstrap_list = []
                    N_list = []
                    S_list = []
                    N_dN_list = []
                    S_dS_list = []
                    sum_N = 0
                    sum_S = 0
                    sum_N_dN = 0
                    sum_S_dS = 0
                    mean_dN = 0
                    mean_dS = 0
                    dN_dS = 0

                    for i in range(length):
                        bootstrap_list.append(random.choice(list_lines))
                    if len(bootstrap_list) != length:
                        print "ERROR - incorrect number sampled"
                    for x in bootstrap_list:
                        N_value = x.split(",")[2]
                        N_list.append(N_value)
                        N_dN_value = x.split(",")[7]
                        N_dN_list.append(N_dN_value)
                        S_value = x.split(",")[3]
                        S_list.append(S_value)
                        S_dS_value = x.split(",")[8]
                        S_dS_list.append(S_dS_value)

                    for i in N_list:
                        i = float(i)
                        sum_N = sum_N + i
                    for i in N_dN_list:
                        i = float(i)
                        sum_N_dN = sum_N_dN + i
                    mean_dN = float(sum_N_dN / sum_N)
                    out_dN.write(str(mean_dN))
                    out_dN.write("\n")
                    dN_list.append(mean_dN)

                    for i in S_list:
                        i = float(i)
                        sum_S = sum_S + i
                    for i in S_dS_list:
                        i = float(i)
                        sum_S_dS = sum_S_dS + i
                    mean_dS = float(sum_S_dS / sum_S)
                    out_dS.write(str(mean_dS))
                    out_dS.write("\n")
                    dS_list.append(mean_dS)

                    dN_dS = float(mean_dN / mean_dS)
                    out_dN_dS.write(str(dN_dS))
                    out_dN_dS.write("\n")
                    dN_dS_list.append(dN_dS)

    lower_quart_dN = np.percentile(dN_list,2.5)
    print ("dN_list lower quartile = "), lower_quart_dN
    upper_quart_dN = np.percentile(dN_list,97.5)
    print ("dN_list upper quartile = "), upper_quart_dN
    median_dN = np.median(dN_list)
    print ("median dN_list = "), median_dN
    print ("\n")

    lower_quart_dS = np.percentile(dS_list,2.5)
    print ("dS_list lower quartile = "), lower_quart_dS
    upper_quart_dS = np.percentile(dS_list,97.5)
    print ("dS_list upper quartile = "), upper_quart_dS
    median_dS = np.median(dS_list)
    print ("median dS_list = "), median_dS
    print ("\n")

    lower_quart_dN_dS = np.percentile(dN_dS_list,2.5)
    print ("dN_dS_list lower quartile = "), lower_quart_dN_dS
    upper_quart_dN_dS = np.percentile(dN_dS_list,97.5)
    print ("dN_dS_list upper quartile = "), upper_quart_dN_dS
    median_dN_dS = np.median(dN_dS_list)
    print ("median dN_dS_list = "), median_dN_dS

if __name__ == '__main__':
    main()