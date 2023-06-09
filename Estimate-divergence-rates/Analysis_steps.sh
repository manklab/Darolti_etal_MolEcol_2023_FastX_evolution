#Analysis to estimate divergence rates
#Pipeline to:
#	1. Identify reciprocal orthologs
#	2. Identify open reading frames
#	3. Align orthologs
#	4. Mask alignments for poorly aligned regions
#	5. Calculate evolutionary rates


## 1. Identify reciprocal orthologs ##

# We have previously obtained de novo transcripts for our target species (P. reticulata, P. wingei, P. picta, P. parae, P. latipinna, G. holbrooki) by mapping RNA-seq reads to species-specific genome assemblies (Wright et al. 2017, Darolti et al. 2021, Sandkam et al. 2022)

# Obtain longest isoforms from Ensembl cds fasta files. This script processes an Ensembl fasta file and picks the longest isoform for each gene. Outputs a new fasta file with ending in _longest.fasta
python 00.get-longest-isoform.py ../1.Identify_reciprocal_orthologs/Outgroup_species/Poecilia_formosa.PoeFor_5.1.2.cds.all.fa
python 00.get-longest-isoform.py ../1.Identify_reciprocal_orthologs/Outgroup_species/Xiphophorus_maculatus.Xipmac4.4.2.cds.all.fa
python 00.get-longest-isoform.py ../1.Identify_reciprocal_orthologs/Outgroup_species/Oryzias_latipes.MEDAKA1.cds.all.fa
python 00.get-longest-isoform.py ../1.Identify_reciprocal_orthologs/Outgroup_species/Poecilia_latipinna.P_latipinna-1.0.cds.all.fa
python 00.get-longest-isoform.py ../1.Identify_reciprocal_orthologs/Outgroup_species/Poecilia_reticulata.Guppy_female_1.0_MT.cds.all.fa
python 00.get-longest-isoform.py ../1.Identify_reciprocal_orthologs/Outgroup_species/Gambusia_affinis.ASM309773v1.cds.all.fa

# Rename P. reticulata and P. latipinna files from Ensembl to not confuse them with the de novo transcript fasta files for these species
mv Poeciliareticulata_longest.fasta PoeciliareticulataEnsembl_longest.fasta
mv Poecilialatipinna_longest.fasta PoecilialatipinnaEnsembl_longest.fasta

target_species=(Poeciliareticulata Poeciliawingei Poeciliapicta Poeciliaparae Poecilialatipinna Gambusiaholbrooki PoeciliareticulataEnsembl)

# Reciprocal best-hit blast. A list of input genomes is blasted against each other and the output is written in provided format.
for sp in ${target_species[@]}; do
	python 01.run-blastall.py ${sp}_longest.fasta Poeciliaformosa_longest.fasta Xiphophorusmaculatus_longest.fasta Oryziaslatipes_longest.fasta -o ../run-blastall/$sp -e 10e-10 -f "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" -p 1 -id 30 -b blastn
done

# Identify top blast hits
for sp in ${target_species[@]}; do
	python 02.top-blasthit.py ../run-blastall/$sp ../run-blastall/$sp/${sp}_reciprocal_orthologs
done

# Identify clusters of reciprocal orthologs (orthogroups). 
for sp in ${target_species[@]}; do
	python 03.ortho-cluster.py ../run-blastall/$sp/${sp}_reciprocal_orthologs.pkl ../run-blastall/$sp/${sp}_orthocluster ${sp},Poeciliaformosa,Xiphophorusmaculatus,Oryziaslatipes
done


## 2. Identify open reading frames ##

# Prepare sequences for blastx analysis. This script takes orthogroups and creates a folder for each. In each folder, a query fasta file is made with each nucleotide sequence. A db fasta file is made with the protein sequence of one of the species (ie the best annotated genome).
# Make folders containing fasta files of all species for each analysis
for sp in ${target_species[@]}; do
	python 04.make-blastx-input.py Oryzias_latipes.MEDAKA1.pep.all.fa Oryziaslatipes ../run-blastall/$sp/${sp}_orthocluster.pkl ../all_fastafiles_${sp}/ ../blastx_${sp}/
done

# Run blastx analysis. 
for sp in ${target_species[@]}; do
	python 05.run-blastx.py ../blastx_${sp}/ -e 10e-10 -f "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe sframe sseq" -b blastx
done

# Identify open reading frames
for sp in ${target_species[@]}; do
	python 06.orf-prediction.py ../blastx_${sp}/ ../invalid_ORF_orthogroups_${sp}/
done


## 3. Align orthologs ##

# Prepare orthologous sequences for prank alignment
for sp in ${target_species[@]}; do
	python 07.make-prank-input.py ../blastx_${sp}/ stripped.cdna.fa
done

# Run prank. Need to first make rooted tree files as input of the format: (((Species1,Species2),Species3),Species4); where Species1 is in turn each of the target species, Species2 is "Poeciliaformosa", Species3 is "Xiphophorusmaculatus" and Species4 is "Oryziaslatipes"
for sp in ${target_species[@]}; do
	python 08.run-prank.py ../blastx_${sp}/ ../trees/rooted_tree_${sp}.txt
done

# Remove gaps in alignments and remove alignments that are shorter than set length threshold after gap removal.
for sp in ${target_species[@]}; do
	python 09.remove-gaps.py ../blastx_${sp}/ ../invalid_length_orthogroups_${sp}/ -cutoff 300
done


## 4. Mask alignments for poorly aligned regions ##

# Convert alignment fasta file into phylip format
for sp in ${target_species[@]}; do
	python 10.convert-fasta-phylip.py ../blastx_${sp}/ prank.gapsrm.fa
done

# Prepare sequences for PAML. Clean data option is turned off as it is required for running SWAMP. Need to first make unrooted tree files as input of the format: ((Species1,Species2),Species3,Species4); where Species1 is in turn each of the target species, Species2 is "Poeciliaformosa", Species3 is "Xiphophorusmaculatus" and Species4 is "Oryziaslatipes"
for sp in ${target_species[@]}; do
	python 11.make-paml-input.py ../blastx_${sp}/ ../trees/unrooted_tree_${sp}.txt -cleandata 0 -model 0 -nssites 0 -fixomega 0 -omega .4
done

# Run PAML. (Note - PAML does not like long file names and paths so might need to rename each input folder to something shorter)
for sp in ${target_species[@]}; do
	python 12.run-paml.py ../blastx_${sp}/
done

# Run SWAMP. See manual (https://github.com/peterwharrison/SWAMP) for running SWAMP (http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4251194/). Often requires trying multiple windows and thresholds. Run branch-site model before and after masking to identify genes under positive selection. See if these genes are identified because the alignments are poor quality. Use this to guide the choice of masking criteria. Need to first make branchcodes files for each analysis (one file with all branchcodes from the PAML run and one file only with the branch code of interest)
for sp in ${target_species[@]}; do
	python SWAMP.py -i ../blastx_${sp}/ -b ../branchcodes/all_branchcodes_${sp}.txt -m 100 -t 6 -w 15 > ../swamp_stats/info_w15t6m100_${sp}.txt
	python 00.clean-folder-rename.py ../blastx_${sp}/ _masked.phy .phy
	python SWAMP.py -i ../blastx_${sp}/ -b ../branchcodes/only_targetsp_branchcode_${sp}.txt -m 100 -t 2 -w 5 > ../swamp_stats/info_w15t6_w5t2m100_${sp}.txt
done


## 5. Calculate evolutionary rates ##

# Remove masked sites (Ns) in the phylip alignment files and remove alignments that are shorter than set length threshold after Ns removal.
for sp in ${target_species[@]}; do
	python 13.remove-N.py ../blastx_${sp}/ ../invalid_Ns_${sp}/ ${sp},Poeciliaformosa,Xiphophorusmaculatus,Oryziaslatipes -cutoff 300
done

# Convert alignment fasta file into phylip format
for sp in ${target_species[@]}; do
	python 10.convert-fasta-phylip.py ../blastx_${sp}/ _masked.Nrm.fa
done

# Prepare sequences for PAML. Clean data option is turned on. Need to first make unrooted tree file with specified branch of interest: ((Species1#1,Species2),Species3,Species4); where Species1 is in turn each of the target species, Species2 is "Poeciliaformosa", Species3 is "Xiphophorusmaculatus" and Species4 is "Oryziaslatipes"
for sp in ${target_species[@]}; do
	python 11.make-paml-input.py ../blastx_${sp}/ ../trees/unrooted_tree_branchofinterest_${sp}.txt -cleandata 1 -model 2 -nssites 0 -fixomega 0 -omega .4
done

# Run PAML
for sp in ${target_species[@]}; do
	python 12.run-paml.py ../blastx_${sp}/
done

# Extract branch lengths and filter orthologs with branch-specific dS >= 2 or if SdS for that branch is <= 1. Change branch of interest accordingly.
for sp in ${target_species[@]}; do
	python 14.paml-extract-branch-lengths-dSSdSfilter.py ../blastx_${sp}/ 6..3 ../divergence_estimates/${sp}/${sp}_branchmodeltest_dSSdSfilter.txt
done

# Obtain gene position on chromosomes
for sp in ${target_species[@]}; do
	python 15.genome-position.py ../run-blastall/$sp/${sp}_orthocluster.txt Oryziaslatipes $sp ../gtfs/${sp}_gtf ../divergence_estimates/${sp}/${sp}_branchmodeltest_dSSdSfilter.txt ../divergence_estimates/${sp}/${sp}_branchmodeltest_dSSdSfilter_posinfo.txt
done

# Our P. reticulata de novo genome assembly was anchored to the P. reticulata reference genome, which has a large inversion on the X chromosome. This inversion is present in the reference genome but absent from all our study populations. We thus need to correct for this inversion.
reticulata_species=(Poeciliareticulata, PoeciliareticulataEnsembl)
for sp in ${reticulata_species[@]}; do
	python 15.reposition-inversion.py ../divergence_estimates/${sp}/${sp}_branchmodeltest_dSSdSfilter_posinfo.txt ../divergence_estimates/${sp}/${sp}_branchmodeltest_dSSdSfilter_posinfo_inversion.txt
done

# Split genes into autosomal and sex-linked categories based on their chromosomal location. The X-linked category excluded genes in the previously identified pseudoautosomal regions (0-5Mb and >26Mb of P. reticulata chromosome 12, and >20Mb of P. wingei, P. picta, P. parae, P. latipinna and G. holbrooki chromosomes that are syntenic to the guppy chromosome 12 (Darolti et al. 2019; Sandkam et al. 2021)). 
for sp in ${reticulata_species[@]}; do
	python 16.split-categories-pret.py ../divergence_estimates/${sp}/${sp}_branchmodeltest_dSSdSfilter_posinfo_inversion.txt ../divergence_estimates/${sp}/${sp}_autosome.txt ../divergence_estimates/${sp}/${sp}_sexchromosome_all.txt ../divergence_estimates/${sp}/${sp}_sexchromosome_noPAR.txt ../divergence_estimates/${sp}/${sp}_PAR.txt ../divergence_estimates/${sp}/${sp}_S1.txt ../divergence_estimates/${sp}/${sp}_S2.txt
done

python 16.split-categories-pwin.py ../divergence_estimates/Poeciliawingei/Poeciliawingei_branchmodeltest_dSSdSfilter_posinfo.txt ../divergence_estimates/Poeciliawingei/Poeciliawingei_autosome.txt ../divergence_estimates/Poeciliawingei/Poeciliawingei_sexchromosome_all.txt ../divergence_estimates/Poeciliawingei/Poeciliawingei_sexchromosome_noPAR.txt ../divergence_estimates/Poeciliawingei/Poeciliawingei_PAR.txt ../divergence_estimates/Poeciliawingei/Poeciliawingei_S1.txt ../divergence_estimates/Poeciliawingei/Poeciliawingei_S2.txt

other_species=(Poeciliapicta Poeciliaparae Poecilialatipinna Gambusiaholbrooki)
for sp in ${other_species[@]}; do
	python 16.split-categories-ppic-ppar-plat-ghol.py ../divergence_estimates/${sp}/${sp}_branchmodeltest_dSSdSfilter.txt ../divergence_estimates/${sp}/${sp}_autosome.txt ../divergence_estimates/${sp}/${sp}_sexchromosome_all.txt ../divergence_estimates/${sp}/${sp}_sexchromosome_noPAR.txt ../divergence_estimates/${sp}/${sp}_PAR.txt
done

# Estimate confidence intervals
for sp in ${target_species[@]}; do
	python 17.randomisation.py ../divergence_estimates/${sp}/${sp}_autosome.txt ../divergence_estimates/${sp}/${sp}_autosome_dN_boot.txt ../divergence_estimates/${sp}/${sp}_autosome_dS_boot.txt ../divergence_estimates/${sp}/${sp}_autosome_dNdS_boot.txt
	python 17.randomisation.py ../divergence_estimates/${sp}/${sp}_sexchromosome_noPAR.txt ../divergence_estimates/${sp}/${sp}_sexchromosome_noPAR_dN_boot.txt ../divergence_estimates/${sp}/${sp}_sexchromosome_noPAR_dS_boot.txt ../divergence_estimates/${sp}/${sp}_sexchromosome_noPAR_dNdS_boot.txt
done

# Test for significant differences in rates of divergence between autosomal and X-linked genes
for sp in ${target_species[@]}; do
	python 18.subsampling-test.py ../divergence_estimates/${sp}_autosome.txt ../divergence_estimates/${sp}_sexchromosome_noPAR.txt
done
