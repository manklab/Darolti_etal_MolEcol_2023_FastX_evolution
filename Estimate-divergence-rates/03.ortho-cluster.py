#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

''' Get ORTHOGROUPS
This script takes the pickle output of reciprocal orthologs and identifies clusters
of reciprocal orthologs.'''
#==============================================================================
import argparse
import sys
import cPickle as pickle
from collections import defaultdict
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
                    help="A python pickle file containing a list of lists with reciprocal best-hits")
parser.add_argument("outfile", type=str,
                    help="The prefix of the outfile [script adds .pkl and .txt]")
parser.add_argument("species", type=str,
                    help="List of species names eg 'Poeciliareticulata,Poeciliaformosa,Xiphophorusmaculatus,Oryziaslatipes' ")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Functions=====================================================================
#==============================================================================
def read_pkl(pkl_file):
    '''Reads in a pickled file and returns the stored object'''
    with open(pkl_file, "rb") as infile:
        return pickle.load(infile)

def get_ortho_dict(recip_besthit_all, SPECIES):
    ''''For every gene in each species find the reciprocal best hit in every other species.
    Produces a dictionary of dictionaries. 1 dictionary per species, within that the number of genes.
    '''
    multi_species_rbh = {}
    SPECIES = SPECIES.split(",")
    for s in SPECIES:
        print s
        species_rbh = defaultdict(list)
        for rbh in recip_besthit_all:
            for g1, g2 in rbh:
                if g1.startswith(s):
                    species_rbh[g1].append(g2)
                if g2.startswith(s):
                    species_rbh[g2].append(g1)
        multi_species_rbh[s] = species_rbh
    return multi_species_rbh


def multi_species_ortholog(ms_rbh):
    '''This function determines if genes in all species are reciprocal besthits
    to each other. It uses the ortho dict and finds the intersection of all 
    reciprocal best hits. If the intersection is equal to the number of used
    species a one-to-one ortholog group has been defined.'''
    SPECIES = ms_rbh.keys()
    target = ms_rbh[SPECIES[0]]
    rest = {}
    ortho_groups = []
    for s in SPECIES[1:]:
        rest.update(ms_rbh[s])
    for gene in target:
        cluster = []
        rbh = target[gene]
        cadd = rbh+[gene]
        cluster.append(set(cadd))
        for r in rbh:
            radd = rest[r] + [r]
            cluster.append(set(radd))
        if len(set.intersection(*cluster)) == len(SPECIES):
            ortho_groups.append(list(set.intersection(*cluster)))
    return ortho_groups

#==============================================================================
#Main==========================================================================
#==============================================================================
def main():
    # Read in a pickle file that contains a list of lists. Each sublist
    # refers to the reciprocal besthit between two species. All possible
    # combinations for all species should be included.
    recip_besthit_all = read_pkl(args.infile)
    print "Reciprocal besthits for each species combination:"
    for l in recip_besthit_all:
        g1, g2 = l[0]
        print g1, g2
        print g1.split("_")[0], g2.split("_")[0], len(l)
    multi_species_rbh = get_ortho_dict(recip_besthit_all, args.species)
    print multi_species_rbh.keys()
    ortho_groups = multi_species_ortholog(multi_species_rbh)
    print "%i one-to-one ortholog groups were detected!" % len(ortho_groups)
    
    print "Starting to pickle..."
    with open(args.outfile+".pkl", "w") as outfile:
        pickle.dump(ortho_groups, outfile, -1)

    print "Writing to file ..."
    count = 0
    with open(args.outfile+".txt", "w") as outtext:
        for gene in ortho_groups:
            count += 1
            outtext.write(str(gene))
            outtext.write("\n")
    print "Number of orthogroups =", count


if __name__ == '__main__':
    main()