#!/usr/bin/python2.7
# -*- coding: utf-8 -*-

''' Test for differences in diveregnce estimates between autosomal and X-linked loci by using subsampling without replacement statistical tests '''
#==============================================================================
import argparse
import sys
import random
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("larger_dataset", type=str,
                    help="Larger file - in this case, the autosomal genes file")
parser.add_argument("smaller_dataset", type=str,
                    help="Samller file - in tis case, the X-linked genes file")
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================

#==============================================================================
def main():

    with open(args.larger_dataset, 'r') as larger:
        with open(args.smaller_dataset, 'r') as smaller:
            #get names of genes
            genelist_larger = []
            count = 0
            for line in larger:
                if count > 0:
                    genelist_larger.append(line)
                count += 1
            genelist_smaller = []
            count = 0
            for line in smaller:
                if count > 0:
                    genelist_smaller.append(line)
                count += 1

            print len(genelist_larger), len(genelist_smaller)

            # Get mean pN, pS and pNpS for the larger and smaller files
            list_N_larger = []
            list_S_larger = []
            list_NdN_larger = []
            list_SdS_larger = []
            for line in genelist_larger:
                list_N_larger.append(float(line.split(",")[2]))
                list_S_larger.append(float(line.split(",")[3]))
                list_SdS_larger.append(float(line.split(",")[10]))
                list_NdN_larger.append(float(line.split(",")[9]))
            final_sum_N_larger = sum(list_N_larger)
            final_sum_S_larger = sum(list_S_larger)
            final_sum_NdN_larger = sum(list_NdN_larger)
            final_sum_SdS_larger = sum(list_SdS_larger)
            dS_larger = float(final_sum_SdS_larger / final_sum_S_larger) + 0.000001
            dN_larger = float(final_sum_NdN_larger / final_sum_N_larger) + 0.000001
            dNdS_larger = float(dN_larger/dS_larger)
            print "check pN_larger, pS_larger, pNpS_larger = ", dN_larger, dS_larger, dNdS_larger 

            list_N_smaller = []
            list_S_smaller = []
            list_NdN_smaller = []
            list_SdS_smaller = []
            for line in genelist_smaller:
                list_N_smaller.append(float(line.split(",")[2]))
                list_S_smaller.append(float(line.split(",")[3]))
                list_SdS_smaller.append(float(line.split(",")[10]))
                list_NdN_smaller.append(float(line.split(",")[9]))
            final_sum_N_smaller = sum(list_N_smaller)
            final_sum_S_smaller = sum(list_S_smaller)
            final_sum_NdN_smaller = sum(list_NdN_smaller)
            final_sum_SdS_smaller = sum(list_SdS_smaller)
            dS_smaller = float(final_sum_SdS_smaller / final_sum_S_smaller) + 0.000001
            dN_smaller = float(final_sum_NdN_smaller / final_sum_N_smaller) + 0.000001
            dNdS_smaller = float(dN_smaller/dS_smaller)
            print "check pN_smaller, pS_smaller, pNpS_smaller = ", dN_smaller, dS_smaller, dNdS_smaller 

            # Get 1000 replicate random subsets of length(smaller) from bigger. Calculate dN, dS and dNdS on each replicate 
            # Make list of 1000 dN, dS and dNdS replicates
            count=0
            count_a=0
            dS_perm = []
            dN_perm = []
            dNdS_perm = []
            i = 0
            while i < 1000:
                smaller_list_rand=[]
                smaller_list_rand = random.sample(genelist_larger, len(genelist_smaller))
                i = i +1
                length =  len(set.intersection(set(smaller_list_rand)))
                if length == len(genelist_smaller):
                    count = count+1
                else:
                    count_a = count_a +1
                list_N=[]
                list_S=[]
                list_NdN=[]
                list_SdS=[]
                for line in smaller_list_rand:
                    list_N.append(float(line.split(",")[2]))
                    list_S.append(float(line.split(",")[3]))
                    list_SdS.append(float(line.split(",")[10]))
                    list_NdN.append(float(line.split(",")[9]))
                final_sum_N= sum(list_N)
                final_sum_S= sum(list_S)
                final_sum_NdN= sum(list_NdN)
                final_sum_SdS= sum(list_SdS)
                dS = float(final_sum_SdS / final_sum_S) + 0.000001
                dN = float(final_sum_NdN / final_sum_N) + 0.000001
                dNdS= float(dN/dS)
                dS_perm.append(dS)
                dN_perm.append(dN)
                dNdS_perm.append(dNdS)
            print "Number of times no replacement =", count
            print "Number of times with replacement =", count_a

            # Get p values by comparing the mean dN, dS and dNdS of bigger and smaller to list of 100 replicates dN, dS and dNdS 
            # 1) dS_p
            p=0
            count_p=0
            paml_sorted = sorted(dS_perm,reverse=True)
            middle = paml_sorted[500]
            first = paml_sorted[0]
            if dS_smaller > middle:
                for paml in paml_sorted:
                    if float(dS_smaller) <= paml:
                        count_p = count_p +1
            else:
                for paml in paml_sorted:
                    if float(dS_smaller) >= paml:
                        count_p = count_p +1
            p = (float(count_p)/1000)*2#ONE TAILED TO DO TWO TAILED *2
            print "pS_p = ", p

            # # 2) dN_p
            p=0
            count_p=0
            paml_sorted = sorted(dN_perm,reverse=True)
            middle = paml_sorted[500]
            first = paml_sorted[0]
            if dN_smaller > middle:
                for paml in paml_sorted:
                    if float(dN_smaller) <= paml:
                        count_p = count_p +1
            else:
                for paml in paml_sorted:
                    if float(dN_smaller) >= paml:
                        count_p = count_p +1
            p = (float(count_p)/1000)*2#ONE TAILED TO DO TWO TAILED *2
            print "pN_p = ", p

            # # 3) dNdS_p
            p=0
            count_p=0
            paml_sorted = sorted(dNdS_perm,reverse=True)
            middle = paml_sorted[500]
            first = paml_sorted[0]
            if dNdS_smaller > middle:
                for paml in paml_sorted:
                    if float(dNdS_smaller) <= paml:
                        count_p = count_p +1
            else:
                for paml in paml_sorted:
                    if float(dNdS_smaller) >= paml:
                        count_p = count_p +1
            p = (float(count_p)/1000)*2#ONE TAILED TO DO TWO TAILED *2
            print "pNpS_p = ", p

if __name__ == "__main__":
    main()