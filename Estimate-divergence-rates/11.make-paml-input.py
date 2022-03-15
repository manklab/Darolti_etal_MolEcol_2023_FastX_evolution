#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

''' Prepare PAML
This script prepares sequences for PAML.
Takes each folder of orthologous sequences and removes all files apart from phylip file.'''
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
parser.add_argument("tree", type=str,
                    help="path to tree for paml")
parser.add_argument("-cleandata", type=str,
                    help="")
parser.add_argument("-model", type=str,
                    help="")
parser.add_argument("-nssites", type=str,
                    help="")
parser.add_argument("-fixomega", type=str,
                    help="")
parser.add_argument("-omega", type=str,
                    help="")
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
    count = 0 
    for orthogroup in list_folder(args.infolder):
        for f in list_files(orthogroup):
            if not f.endswith(".phy"):
                os.remove(f)
            else:
                count +=1                
                control_file = f[:-3]+"m"+args.model+"ns"+args.nssites+"w"+args.fixomega+args.omega+".ctl"
                seq = f
                out_file = f[:-3]+"m"+args.model+"ns"+args.nssites+"w"+args.fixomega+args.omega+".txt"

                try:
                    with open(control_file,"wb") as outfile:
                        outfile.write("      seqfile = "+seq+" * sequence data filename"+"\n")
                        outfile.write("     treefile = "+args.tree+"      * tree structure file name"+"\n")
                        outfile.write("      outfile = "+out_file+"           * main result file name"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("        noisy = 9  * 0,1,2,3,9: how much rubbish on the screen"+"\n")
                        outfile.write("      verbose = 1  * 0: concise; 1: detailed, 2: too much"+"\n")
                        outfile.write("      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic"+"\n")
                        outfile.write("                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs"+"\n")
                        outfile.write("    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("*        ndata = 10"+"\n")
                        outfile.write("        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis"+"\n")
                        outfile.write("       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a"+"\n")
                        outfile.write("   aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)"+"\n")
                        outfile.write("                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("        model = "+args.model+"\n")
                        outfile.write("                   * models for codons:"+"\n")
                        outfile.write("                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches"+"\n")
                        outfile.write("                   * models for AAs or codon-translated AAs:"+"\n")
                        outfile.write("                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F"+"\n")
                        outfile.write("                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("      NSsites = "+args.nssites+"  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;"+"\n")
                        outfile.write("                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;"+"\n")
                        outfile.write("                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;"+"\n")
                        outfile.write("                   * 13:3normal>0"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below"+"\n")
                        outfile.write("        Mgene = 0"+"\n")
                        outfile.write("                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff"+"\n")
                        outfile.write("                   * AA: 0:rates, 1:separate"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated"+"\n")
                        outfile.write("        kappa = 2  * initial or fixed kappa"+"\n")
                        outfile.write("    fix_omega = "+args.fixomega+"  * 1: omega or omega_1 fixed, 0: estimate "+"\n")
                        outfile.write("        omega = "+args.omega+" * initial or fixed omega, for codons or codon-based AAs"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha"+"\n")
                        outfile.write("        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)"+"\n")
                        outfile.write("       Malpha = 0  * different alphas for genes"+"\n")
                        outfile.write("        ncatG = 8  * # of categories in dG of NSsites models"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates"+"\n")
                        outfile.write(" RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("   Small_Diff = .5e-6"+"\n")
                        # outfile.write("    cleandata = 0  * remove sites with ambiguity data (1:yes, 0:no)?"+"\n")
                        outfile.write("    cleandata = "+args.cleandata+"  * remove sites with ambiguity data (1:yes, 0:no)?"+"\n")
                        outfile.write("*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed"+"\n")
                        outfile.write("       method = 0  * Optimization method 0: simultaneous; 1: one branch a time"+"\n")
                        outfile.write(""+"\n")
                        outfile.write("* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,"+"\n")
                        outfile.write("* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., "+"\n")
                        outfile.write("* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., "+"\n") 
                        outfile.write("* 10: blepharisma nu."+"\n")
                        outfile.write("* These codes correspond to transl_table 1 to 11 of GENEBANK."+"\n")
                        outfile.write(""+"\n")

                except IOError:
                    print "File does not exit!"

    print "Number of folders modified =", count

if __name__ == '__main__':
    main()