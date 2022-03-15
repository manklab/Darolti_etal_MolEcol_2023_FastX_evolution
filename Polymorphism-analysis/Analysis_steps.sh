# Analysis to identify SNPs, and obtain polymorphism and divergence data for a set of genes. The polymorphism and divergence pipelines impose different
# filtering criteria to identify sites used in SNP calling or calculating rates of evolutoion. Thus, this pipeline identifed a set of codons where all
# sites pass the filtering criteria for both pipelines.

#Pipeline to:
#	1. Map RNA-seq data to reference genome
#	2. Conduct SNP calling
#	3. Identify codons passing the following criteria:
#		a. all sites pass the minimum coverage threshold
#		b. no alignment gaps
#		c. no ambiguity data (Ns)
#	4. Filter divergence data
#	5. Filter polymorphism data
#	6. Identify nonsynonymous and synonymous SNPs
#	7. Tests of selection

# We have previously trimmed RNA-seq FASTQ files for our target species (P. reticulata, P. wingei, P. picta, P. parae, P. latipinna, G. holbrooki) using Trimmomatic v0.36 (see methods in Wright et al. 2017, Darolti et al. 2021, Almeida et al. 2022, Sandkam et al. 2022)
# We have also previously assembled female-specifc genome assemblies for each target species (see methods in Wright et al. 2017, Darolti et al. 2021, Almeida et al. 2022, Sandkam et al. 2022)

## 1. Map RNA-seq data to reference genome ##

target_species=(Poeciliareticulata Poeciliawingei Poeciliapicta Poeciliaparae Poecilialatipinna Gambusiaholbrooki)

# Create STAR genome indexes
for sp in ${target_species[@]}; do
	STAR --runMode genomeGenerate --genomeDir ../genomes/${sp}/ --genomeFastaFiles ${sp}_genome_assembly.fa --runThreadN 12 --limitGenomeGenerateRAM 100000000000
done

# STAR map first pass
for sp in ${target_species[@]}; do
	python 01.STAR-map-firstpass.py ../${sp}/ ../genomes/${sp}/
done

# Concatenate splice junctions
for sp in ${target_species[@]}; do
	python 02.concatenate-splicejunctions.py ../${sp}/firstpass ../${sp}/cat/SJ.out.tab.merged
done

# STAR map second pass
for sp in ${target_species[@]}; do
	python 03.STAR-map-secondpass.py ../${sp}/ ../genomes/${sp}/ ../${sp}/cat/SJ.out.tab.merged
done


## 2. Conduct SNP calling ##

# Sort and index alignment files
for sp in ${target_species[@]}; do
	for sample in ../${sp}/secondpass/*/; do
		[ -L "${sample%/}" ] && continue
		for file in "${sample}"*Aligned.out.sam; do
			samtools view -bh ${file} | samtools sort > ${file%.*}_sorted
			samtools index ${file%.*}_sorted
		done
	done
done

# Create samtools genome indexes
for sp in ${target_species[@]}; do
	samtools faidx ../genomes/${sp}/${sp}_genome_assembly.fa
done

# Samtools mpileup
for sp in ${target_species[@]}; do
	python 04.mpileup.py ../${sp}/secondpass ../genomes/${sp}/${sp}_genome_assembly.fa
done

# Run VarScan mpileup2snp. To include sample names in the output VCF header, need to first prepare a file with one sample name per line (sample_list.txt).
for sp in ${target_species[@]}; do
	java -jar VarScan.v2.3.9.jar mpileup2snp ../${sp}/mpileup/output.pileup --min-coverage 2 --min-avg-qual 20 --p-value 0.05 --strand-filter 1 --min-freq-for-hom 0.85 --min-var-freq 1e-1 --output-vcf 1 --vcf-sample-list ../${sp}/mpileup/sample_list.txt > ../${sp}/mpileup/mpileup2snp.vcf
done

# Filter SNPs
for sp in ${target_species[@]}; do
	python 05.filter-snps.py ../${sp}/mpileup/mpileup2snp.vcf 2 ../${sp}/mpileup/mpileup2snp_filtered.vcf
done


## Identify codons passing the following criteria: 
#		a. all sites pass the minimum coverage threshold
#		b. no alignment gaps
#		c. no ambiguity data (Ns)

# Get coordinates of ORFs for every gene
for sp in ${target_species[@]}; do
	python 06.orf-prediction.py ../blastx_${sp}/ ${sp} ../${sp}/ORF/${sp}_orf_coordinates_transcript.txt
done

# Get ORF position in scaffold. Need as an input a file with positional information of genes (e.g., which scaffold they are located on)
for sp in ${target_species[@]}; do
	python 07.orf-position-within-scaffold.py ../${sp}/ORF/${sp}_orf_coordinates_transcript.txt ../${sp}/ORF/${sp}_gene_position_info.txt ../gtfs/${sp}_gtf ../${sp}/ORF/${sp}_orf_coordinates_scaffold.txt
done

# Extract sites that pass the minimum coverage threshold
for sp in ${target_species[@]}; do
	python 08.extract-sites-min-coverage.py ../${sp}/mpileup/output.pileup 20 2 ../${sp}/ORF/${sp}_gene_position_info.txt ../${sp}/codons_passig_criteria/${sp}-extract-sites-min-cov.pkl
done

# Extract codons where all sites pass the minimum coverage threshold
for sp in ${target_species[@]}; do
	python 09.extract-codons-min-coverage.py ../${sp}/codons_passig_criteria/${sp}-extract-sites-min-cov.pkl ../${sp}/ORF/${sp}_gene_position_info.txt ../${sp}/ORF/${sp}_orf_coordinates_scaffold.txt ../${sp}/codons_passig_criteria/${sp}-extract-codons-min-cov.txt
done

# Extract codons without gaps. ! Note - input folder contains folders of orthogroup sequences aligned by prank BEFORE gap removal. i.e. output of Estimate-divergence-rates/08.run-prank.py
for sp in ${target_species[@]}; do
	python 10.extract-codons-no-gaps.py ../blastx_${sp}/ ${sp} ../${sp}/codons_passig_criteria/${sp}_out_
done

# Extract codons without Ns. ! Note - input folder contains folders of orthogroup sequences aligned by prank AFTER gap removal and masking with swamp but BEFORE N removal.
for sp in ${target_species[@]}; do
	python 11.extract-codons-no-Ns.py ../blastx_${sp}/ ${sp} ${sp},Poeciliaformosa,Xiphophorusmaculatus,Oryziaslatipes ../${sp}/codons_passig_criteria/${sp}_out_codons_with_gaps.txt ../${sp}/codons_passig_criteria/${sp}_out_codons_without_Ns.txt
done

# Check filtered codons. Prefiltering folder should have the data BEFORE gap removal. Postfiltering folder should have the data AFTER gap removal, masking and N removal
for sp in ${target_species[@]}; do
	python 12.check-filtered-codons.py ../prefiltering_blastx_${sp}/ ../postfiltering_blastx_${sp}/ ${sp} ../${sp}/codons_passig_criteria/${sp}_out_codons_without_gaps.txt ../${sp}/codons_passig_criteria/${sp}_out_codons_without_Ns.txt
done


## 4. Filter divergence data ##

# Filter codons passing gap, ambiguity and coverage filters
for sp in ${target_species[@]}; do
	python 13.filter-aln-gaps-Ns-coverage.py ../prefiltering_blastx_${sp}/ ../postfiltering_blastx_${sp}/ ${sp} ../${sp}/codons_passig_criteria/${sp}_out_codons_without_gaps.txt ../${sp}/codons_passig_criteria/${sp}_out_codons_without_Ns.txt ../${sp}/codons_passig_criteria/${sp}-extract-codons-min-cov.txt ../${sp}/ORF/${sp}_orf_coordinates_scaffold.txt ../divergence_estimates/${sp}/${sp}_branchmodeltest_dSSdSfilter.txt ../${sp}/filter-aln-gaps-Ns-coverage -cutoff 300
done

# Convert alignment fasta file into phylip format
for sp in ${target_species[@]}; do
	python 14.convert-fasta-phylip.py ../${sp}/filter-aln-gaps-Ns-coverage prank.gaps.Ns.coverage.filtered.fa
done

# Prepare sequences for PAML. Clean data option is turned on. Need to first make unrooted tree file with specified branch of interest: ((Species1#1,Species2),Species3,Species4); where Species1 is in turn each of the target species, Species2 is "Poeciliaformosa", Species3 is "Xiphophorusmaculatus" and Species4 is "Oryziaslatipes"
for sp in ${target_species[@]}; do
	python ../Estimate-divergence-rates/11.make-paml-input.py ../${sp}/filter-aln-gaps-Ns-coverage ../trees/unrooted_tree_branchofinterest_${sp}.txt -cleandata 1 -model 2 -nssites 0 -fixomega 0 -omega .4
done

# Run PAML
for sp in ${target_species[@]}; do
	python ../Estimate-divergence-rates/12.run-paml.py ../${sp}/filter-aln-gaps-Ns-coverage
done

# Extract branch lengths and filter orthologs with branch-specific dS >= 2 or if SdS for that branch is <= 1. Change branch of interest accordingly.
for sp in ${target_species[@]}; do
	python ../Estimate-divergence-rates/14.paml-extract-branch-lengths-dSSdSfilter.py ../${sp}/filter-aln-gaps-Ns-coverage 6..4 ../divergence_estimates_filtered/${sp}/${sp}_branchmodeltest_dSSdSfilter.txt
done


## 5. Filter polymorphism data ##

# Filter SNPs to only keep those that are positioned in codons without gaps, without Ns and that pass the coverage threshold
for sp in ${target_species[@]}; do
	python 15.filter-polymorphism-data.py ../${sp}/mpileup/mpileup2snp_filtered.vcf ../${sp}/ORF/${sp}_orf_coordinates_transcript.txt ../${sp}/ORF/${sp}_gene_position_info.txt ../${sp}/codons_passig_criteria/${sp}_out_codons_without_gaps.txt ../${sp}/codons_passig_criteria/${sp}_out_codons_without_Ns.txt ../${sp}/codons_passig_criteria/${sp}-extract-codons-min-cov.txt ../${sp}/ORF/${sp}_orf_coordinates_scaffold.txt ../${sp}/mpileup/filtered_snps.txt
done


## 6. Identify nonsynonymous and synonymous SNPs ##
for sp in ${target_species[@]}; do
	python 16.identify-NpN-SpS.py ../${sp}/mpileup/filtered_snps.txt ../${sp}/ORF/${sp}_orf_coordinates_transcript.txt ../${sp}/ORF/${sp}_gene_position_info.txt ../${sp}/ORF/${sp}_orf_coordinates_scaffold.txt ../genomes/${sp}/${sp}_genome_assembly.fa ../${sp}/NpN_SpS/${sp}_NpN_SpS.txt
done

# Find genes that have both divergence and polymorphism data 
for sp in ${target_species[@]}; do
	python 17.filter-divergence-and-polymorphism-genes.py ../divergence_estimates_filtered/${sp}/${sp}_branchmodeltest_dSSdSfilter.txt ../${sp}/NpN_SpS/${sp}_NpN_SpS.txt ../${sp}/ORF/${sp}_orf_coordinates_transcript.txt ../${sp}/NpN_SpS/${sp}_genes_both_divergence_and_polymorphism.txt
done

# Obtain gene position on chromosomes.
for sp in ${target_species[@]}; do
	python 18.genome-position.py ../gtfs/${sp}_gtf ../${sp}/NpN_SpS/${sp}_genes_both_divergence_and_polymorphism.txt ../${sp}/NpN_SpS/${sp}_genes_both_divergence_and_polymorphism_posinfo.txt
done

# Our P. reticulata de novo genome assembly was anchored to the P. reticulata reference genome, which has a large inversion on the X chromosome. This inversion is present in the reference genome but absent from all our study populations. We thus need to correct for this inversion.
python ../Estimate-divergence-rates/15.reposition-inversion.py ../${sp}/NpN_SpS/${sp}_genes_both_divergence_and_polymorphism_posinfo.txt ../${sp}/NpN_SpS/${sp}_genes_both_divergence_and_polymorphism_posinfo_inversion.txt

# Split genes into autosomal and sex-linked categories based on their chromosomal location. The X-linked category excluded genes in the previously identified pseudoautosomal regions (0-5Mb and >26Mb of P. reticulata chromosome 12, and >20Mb of P. wingei, P. picta, P. parae, P. latipinna and G. holbrooki chromosomes that are syntenic to the guppy chromosome 12 (Darolti et al. 2019; Sandkam et al. 2021)). 
python ../Estimate-divergence-rates/16.split-categories-pret.py ../Poeciliareticulata/NpN_SpS/Poeciliareticulata_genes_both_divergence_and_polymorphism_posinfo_inversion.txt ../Poeciliareticulata/NpN_SpS/Poeciliareticulata_autosome.txt ../Poeciliareticulata/NpN_SpS/Poeciliareticulata_sexchromosome_all.txt ../Poeciliareticulata/NpN_SpS/Poeciliareticulata_sexchromosome_noPAR.txt ../Poeciliareticulata/NpN_SpS/Poeciliareticulata_PAR.txt ../Poeciliareticulata/NpN_SpS/Poeciliareticulata_S1.txt ../Poeciliareticulata/NpN_SpS/Poeciliareticulata_S2.txt

python ../Estimate-divergence-rates/16.split-categories-pwin.py ../Poeciliawingei/NpN_SpS/Poeciliawingei_genes_both_divergence_and_polymorphism_posinfo_inversion.txt ../Poeciliawingei/NpN_SpS/Poeciliawingei_autosome.txt ../Poeciliawingei/NpN_SpS/Poeciliawingei_sexchromosome_all.txt ../Poeciliawingei/NpN_SpS/Poeciliawingei_sexchromosome_noPAR.txt ../Poeciliawingei/NpN_SpS/Poeciliawingei_PAR.txt ../Poeciliawingei/NpN_SpS/Poeciliawingei_S1.txt ../Poeciliawingei/NpN_SpS/Poeciliawingei_S2.txt

other_species=(Poeciliapicta Poeciliaparae Poecilialatipinna Gambusiaholbrooki)
for sp in ${other_species[@]}; do
	python ../Estimate-divergence-rates/16.split-categories-ppic-ppar-plat-ghol.py ../${sp}/NpN_SpS/${sp}_genes_both_divergence_and_polymorphism_posinfo.txt ../${sp}/NpN_SpS/${sp}_autosome.txt ../${sp}/NpN_SpS/${sp}_sexchromosome_all.txt ../${sp}/NpN_SpS/${sp}__sexchromosome_noPAR.txt ../${sp}/NpN_SpS/${sp}_PAR.txt
done

# Estimate confidence intervals
for sp in ${target_species[@]}; do
	python 19.randomisation.py ../${sp}/NpN_SpS/${sp}_autosome.txt ../${sp}/NpN_SpS/${sp}_autosome_dN_boot.txt ../${sp}/NpN_SpS/${sp}_autosome_dS_boot.txt ../${sp}/NpN_SpS/${sp}_autosome_dNdS_boot.txt
	python 19.randomisation.py ../${sp}/NpN_SpS/${sp}__sexchromosome_noPAR.txt ../${sp}/NpN_SpS/${sp}__sexchromosome_noPAR_dN_boot.txt ../${sp}/NpN_SpS/${sp}__sexchromosome_noPAR_dS_boot.txt ../${sp}/NpN_SpS/${sp}__sexchromosome_noPAR_dNdS_boot.txt
done

# Test for significant differences in rates of divergence between autosomal and X-linked genes
for sp in ${target_species[@]}; do
	python 20.subsampling-test.py ../${sp}/NpN_SpS/${sp}_autosome.txt ../${sp}/NpN_SpS/${sp}_sexchromosome_noPAR.txt
done


## 7. Tests of selection ##

# Get matrix file
for sp in ${target_species[@]}; do
	python 21.filter-divergence-and-polymorphism-genes-matrix.py ../divergence_estimates_filtered/${sp}/${sp}_branchmodeltest_dSSdSfilter.txt ../${sp}/NpN_SpS/${sp}_NpN_SpS.txt ../${sp}/ORF/${sp}_orf_coordinates_transcript.txt ../${sp}/NpN_SpS/${sp}_genes_both_divergence_and_polymorphism_matrix.txt
done

# Run R script separately for each target species
22.tests-of-selection.R


