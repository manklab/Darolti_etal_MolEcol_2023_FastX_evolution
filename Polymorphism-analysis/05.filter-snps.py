#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
''' Filter VCF file
Filters a VCF file to identify sites that are valid SNPs. Valid SNPs defined as sites
where minimum no. of individuals have DP >= 20 and the alternative allele frequency is >= 0.1.
Excludes missing data. Removes SNPs where MAF=1, REF=N and triallelic sites.'''
#==============================================================================
import argparse
import sys
import os
from collections import defaultdict
from collections import Counter
import vcf
import time
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("infile", type=str,
					help="A vcf file")
parser.add_argument("coverage_threshold", type=str,
					help="Minimum no. of individuals with DP >= 20 in order to call a SNP at a given site")
parser.add_argument("outfile", type=str,
					help="An outfile file of SNPs")
# This checks if the user supplied any arguments. If not, help is printed.
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(1)
# Shorten the access to command line arguments.
args = parser.parse_args()
#==============================================================================
#Main==========================================================================
#==============================================================================
def main():

	#VCF FORMAT
	#DP= Quality Read Depth of bases with Phred score >= 20

	print "Minimum no. of individuals with DP >= 20 in order to call a SNP at a given site =", args.coverage_threshold

	startTime = time.time()
	count_valid_SNPs = 0
	count_triallelic = 0
	count_N_SNPs = 0
	count_MAF1_SNPs = 0
	vcf_lines = 0
	bp = ["A","G","T","C"]
	
	vcf_reader = vcf.Reader(open(args.infile,"r"))
	with open(args.outfile, "w") as outfile:
		# Print header
		outfile.write("CHROM")
		outfile.write("\t")
		outfile.write('POS')
		outfile.write("\t")
		outfile.write('REF')
		outfile.write("\t")
		outfile.write('ALT')
		outfile.write("\n")

		for record in vcf_reader:
			sys.stdout.write('%d\r' % (vcf_lines))
			sys.stdout.flush()
			vcf_lines += 1
			sample_DP = 0
			genotypes = []
			# Remove INDEL info, only consider SNPs
			if record.is_snp: 
				# Check INDEL filter has worked. REF should be 1bp
				if len(record.REF) == 1:
					# Check REF is not N
					if record.REF not in bp:
						count_N_SNPs += 1
					else:
						# Check not triallelic site
						if len(record.ALT) != 1:
							count_triallelic += 1
						else:
							for sample in record.samples:
		 						# Ignore missing data (./.) when calculating MAF
								if sample['GT'] != "./.":
		 							genotypes.append(sample['GT'])
		 							# Count number of individuals with DP >= 20
		 							if sample['DP'] != "None" and float(sample['DP']) >= 20:
		 								sample_DP += 1
							
			# Check if pass minimum coverage
			if sample_DP >= float(args.coverage_threshold):
				# Calculate MAF
				allele_list = []
				for genotype in genotypes:
					alleles = genotype.split("/")
					for allele in alleles:
						allele_list.append(allele)
				chromosomes = len(allele_list)
				allele_count = Counter(allele_list)
				#minor allele is always 1 as triallelic sites have been removed
				#minor allele frequency refers to alternative allele frequency
				MAF = float(allele_count["1"])/chromosomes
				# Check if pass MAF frequency
				if MAF >= 0.1:
					# Remove if MAF 1 - ambiguity in the reference?
					if MAF == 1:
						count_MAF1_SNPs += 1
					else:
						count_valid_SNPs += 1
						# Print valid SNPs
						outfile.write(str(record.CHROM))
						outfile.write("\t")
						outfile.write(str(record.POS))
						outfile.write("\t")
						outfile.write(str(record.REF))
						outfile.write("\t")
						outfile.write(str(record.ALT[0]))
						outfile.write("\n")

	print "Number of lines in vcf file =", vcf_lines
	print "Number of valid SNPs =",count_valid_SNPs
	print "Number of triallelic+ SNPs =", count_triallelic
	print "Number of SNPS where REF=N =", count_N_SNPs
	print "Number of SNPS where MAF=1 =", count_MAF1_SNPs

	print "The script took {0} second !".format(time.time() - startTime)

if __name__ == '__main__':
	main()
