#!/usr/bin/python2.6
# -*- coding: utf-8 -*-
''' Calculate confidence intervals for pN, pS and pNpS estimates'''
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
                    help="An infile containing genes with dSfilter values")
parser.add_argument("out_pN_bootstrap", type=str,
                    help="An infile containing genes with dSfilter values")
parser.add_argument("out_pS_bootstrap", type=str,
                    help="An infile containing genes with dSfilter values")
parser.add_argument("out_pN_pS_bootstrap", type=str,
                    help="An infile containing genes with dSfilter values")
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
    for path, subdirs, files in os.walk(current_dir):
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
    pN_pS_list = []
    pN_list = []
    pS_list = []

    count = 0
    with open(args.dSSdSfilter, 'r') as infile:
        for line in infile:
            if count >0:
                length += 1
                list_lines.append(line)
            count += 1
    with open(args.out_pN_bootstrap, 'w') as out_pN:
        with open(args.out_pS_bootstrap, 'w') as out_pS:
            with open(args.out_pN_pS_bootstrap, 'w') as out_pN_pS:
                while bootstrap < 1000:
                    bootstrap += 1
                    bootstrap_list = []
                    N_list = []
                    S_list = []
                    N_pN_list = []
                    S_pS_list = []
                    sum_N = 0
                    sum_S = 0
                    sum_N_pN = 0
                    sum_S_pS = 0
                    mean_pN = 0
                    mean_pS = 0
                    pN_pS = 0

                    for i in range(length):
                        bootstrap_list.append(random.choice(list_lines))
                    if len(bootstrap_list) != length:
                        print "ERROR - incorrect number sampled"
                    for x in bootstrap_list:
                        N_value = x.split(",")[2]
                        N_list.append(N_value)
                        N_pN_value = x.split(",")[9]
                        N_pN_list.append(N_pN_value)
                        S_value = x.split(",")[3]
                        S_list.append(S_value)
                        S_pS_value = x.split(",")[10].rstrip()
                        S_pS_list.append(S_pS_value)

                    for i in N_list:
                        i = float(i)
                        sum_N = sum_N + i
                    for i in N_pN_list:
                        i = float(i)
                        sum_N_pN = sum_N_pN + i
                    mean_pN = float(sum_N_pN / sum_N)
                    out_pN.write(str(mean_pN))
                    out_pN.write("\n")
                    pN_list.append(mean_pN)

                    for i in S_list:
                        i = float(i)
                        sum_S = sum_S + i
                    for i in S_pS_list:
                        i = float(i)
                        sum_S_pS = sum_S_pS + i
                    mean_pS = float(sum_S_pS / sum_S)
                    out_pS.write(str(mean_pS))
                    out_pS.write("\n")
                    pS_list.append(mean_pS)

                    pN_pS = float(mean_pN / mean_pS)
                    out_pN_pS.write(str(pN_pS))
                    out_pN_pS.write("\n")
                    pN_pS_list.append(pN_pS)

    lower_quart_pN = np.percentile(pN_list,2.5)
    print ("pN_list lower quartile = "), lower_quart_pN
    upper_quart_pN = np.percentile(pN_list,97.5)
    print ("pN_list upper quartile = "), upper_quart_pN
    median_pN = np.median(pN_list)
    print ("median pN_list = "), median_pN
    print ("\n")

    lower_quart_pS = np.percentile(pS_list,2.5)
    print ("pS_list lower quartile = "), lower_quart_pS
    upper_quart_pS = np.percentile(pS_list,97.5)
    print ("pS_list upper quartile = "), upper_quart_pS
    median_pS = np.median(pS_list)
    print ("median pS_list = "), median_pS
    print ("\n")

    lower_quart_pN_pS = np.percentile(pN_pS_list,2.5)
    print ("pN_pS_list lower quartile = "), lower_quart_pN_pS
    upper_quart_pN_pS = np.percentile(pN_pS_list,97.5)
    print ("pN_pS_list upper quartile = "), upper_quart_pN_pS
    median_pN_pS = np.median(pN_pS_list)
    print ("median pN_pS_list = "), median_pN_pS

if __name__ == '__main__':
    main()