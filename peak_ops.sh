#PBS -l nodes=1:ppn=4
#PBS -l walltime=5:00:00 
 
#Global variables
PBAM="/fs/project/PAS1143/Michael/public_bam"
PPEAK="/fs/project/PAS1143/Michael/public_peaks"
MBAM="/fs/project/PAS1143/Michael/aln"
MPEAK="/fs/project/PAS1143/Michael/library_peaks"
MUT="/fs/project/PAS1143/Michael/mutation"
GENES="/fs/project/PAS1143/Michael/TSS_hg38/gencode.v27.annotation.gtf"
TSS="/fs/project/PAS1143/common/TSS_annotations/hg38/TSS_hg38_first3cols.bed"
 
#Remove track line of mutation files to use in R script.
for mutfile in "Immortal_AG591" "Immortal_AG695" "Immortal_AH425" "Post_stasis_AI263" "Post_stasis_AI958" "Tumorigenic_AI027";
	do
		 awk '{if (NR!= 1) {print}}' $MUT/${mutfile}_rgb_sorted.bed > $MUT/$mutfile.bed 
	done

#Public data
for pfile in "H3K9Me3_HMEC_1" "H3K9Me3_HMEC_2" "DNase_HAEC_1" "DNase_HAEC_2";
	do
		#Call peaks.
		#macs2 callpeak --broad -t $PBAM/${pfile}_filtered.bam --outdir $PPEAK/$pfile
		#Remove unneeded columns.
		awk '{printf("%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5)}' $PPEAK/$pfile/NA_peaks.broadPeak > $PPEAK/${pfile}_unsorted.bed
		#Sort the file.
		sort -k1,1V -k2,2n $PPEAK/${pfile}_unsorted.bed > $PPEAK/${pfile}.bed
		#Find distance to the closest TSS.
		bedtools closest -d -t first -a $PPEAK/$pfile.bed -b $TSS > $PPEAK/${pfile}_closest.bed
		#Filter out peaks which are proximal to TSS regions (promoters) and those which are distal.
		awk '{if ($10 <= 1500) {printf("%s\t%s\t%s\tproximal\t%s\n", $5, $6, $7, $8)}}' ${pfile}_closest.bed > ${pfile}_proximal.bed
		awk '{if ($10 > 1500) {printf("%s\t%s\t%s\tdistal\t%s\n", $5, $6, $7, $8)}}' ${pfile}_closest.bed > ${pfile}_distal.bed
	done
	
#Michael's libraries
for mfile in "Library_3" "Library_4" "Library_5";
	do
		#Call peaks.
		#macs2 callpeak --broad -t $MBAM/${mfile}_filtered.bam --outdir $MPEAK/$mfile.bed
		#Remove unneeded columns.
		#awk '{printf("%s\t%s\t%s\t%s\t%s\n", $1, $2, $3, $4, $5)}' $MPEAK/$mfile/$mfile.broadPeak > $MPEAK/$mfile.bed
		#Sort the file.
		sort -k1,1V -k2,2n $PPEAK/${pfile}_unsorted.bed > $PPEAK/${pfile}.bed
		#Find distance to the closest TSS.
		bedtools closest -d -t first -a $MPEAK/$mfile.bed -b $TSS > $MPEAK/${mfile}_closest.bed
		#Filter out peaks which are proximal to TSS regions (promoters) and those which are distal.
		awk '{if ($10 <= 1500) {printf("%s\t%s\t%s\t%s\t%s\tproximal\n", $5, $6, $7, $8, $9)}}' ${mfile}_closest.bed > ${mfile}_proximal.bed
		awk '{if ($10 > 1500) {printf("%s\t%s\t%s\t%s\t%s\tdistal\n", $5, $6, $7, $8, $9)}}' ${mfile}_closest.bed > ${mfile}_distal.bed
	done
 
#Concatenate replicates.
cat $PPEAK/DNase_HAEC_1_proximal.bed $PPEAK/DNase_HAEC_2_proximal.bed > $PPEAK/DNase_HAEC_unsorted_proximal.bed
cat $PPEAK/DNase_HAEC_1_distal.bed $PPEAK/DNase_HAEC_2_distal.bed > $PPEAK/DNase_HAEC_unsorted_distal.bed
cat $PPEAK/H3K9Me3_HMEC_1_proximal.bed $PPEAK/H3K9Me3_HMEC_2_proximal.bed > $PPEAK/H3K9Me3_HMEC_unsorted_proximal.bed
cat $PPEAK/H3K9Me3_HMEC_1_distal.bed $PPEAK/H3K9Me3_HMEC_2_distal.bed > $PPEAK/H3K9Me3_HMEC_unsorted_distal.bed
  
#Sort bed files. Concatenate proximal and distal and sort again.
for catfile in "H3K9Me3_HMEC" "DNase_HAEC";
	do
		bedtools sort -i $PPEAK/${catfile}_unsorted_proximal.bed > $PPEAK/${catfile}_proximal.bed
		bedtools sort -i $PPEAK/${catfile}_unsorted_distal.bed > $PPEAK/${catfile}_distal.bed
		cat $PPEAK/${catfile}_proximal.bed $PPEAK/${catfile}_distal.bed > $PPEAK/${catfile}_all_unsorted.bed
		bedtools sort -i $PPEAK/${catfile}_all_unsorted.bed > $PPEAK/${catfile}_all.bed
	done
for catfile2 in "Library_3" "Library_4" "Library_5";
	do
		bedtools sort -i $MPEAK/${catfile2}_unsorted_proximal.bed > $MPEAK/${catfile2}_proximal.bed
		bedtools sort -i $MPEAK/${catfile2}_unsorted_distal.bed > $MPEAK/${catfile2}_distal.bed
		cat $MPEAK/${catfile2}_proximal.bed $MPEAK/${catfile2}_distal.bed > $MPEAK/${catfile2}_all_unsorted.bed
		bedtools sort -i $MPEAK/${catfile2}_all_unsorted.bed > $MPEAK/${catfile2}_all.bed
	done
	
#Add track lines.
sed -i '1 i\track name=H3K9Me3_HMEC_proximal description="proximal H3K9Me3 peaks for HMEC cells (public data)" visibility=2' $PPEAK/H3K9Me3_HMEC_proximal.bed > $PPEAK/H3K9Me3_HMEC_proximal_head.bed
sed -i '1 i\track name=H3K9Me3_HMEC_distal description="distal H3K9Me3 peaks for HMEC cells (public data)" visibility=2' $PPEAK/H3K9Me3_HMEC_distal.bed > $PPEAK/H3K9Me3_HMEC_distal_head.bed
sed -i '1 i\track name=H3K9Me3_HMEC description="H3K9Me3 peaks for HMEC cells (public data)" visibility=2' $PPEAK/H3K9Me3_HMEC_all.bed > $PPEAK/H3K9Me3_head.bed
sed -i '1 i\track name=DNase_HAEC_proximal description="proximal DNase-seq peaks for HAEC cells (public data)" visibility=2' $PPEAK/DNase_HAEC_proximal.bed > $PPEAK/DNase_HAEC_proximal_head.bed
sed -i '1 i\track name=DNase_HAEC_distal description="distal DNase-seq peaks for HAEC cells (public data)" visibility=2' $PPEAK/DNase_HAEC_distal.bed > $PPEAK/DNase_HAEC_distal_head.bed
sed -i '1 i\track name=DNase_HAEC description="DNase-seq peaks for HAEC cells (public data)" visibility=2' $PPEAK/DNase_HAEC.bed > $PPEAK/DNase_HAEC_head.bed
sed -i '1 i\track name=HMEC_primary_proximal description="proximal ATAC-seq peaks for primary HMEC cells" visibility=2' $MPEAK/Library_3_proximal.bed > $MPEAK/Library_3_proximal_head.bed
sed -i '1 i\track name=HMEC_primary_distal description="distal ATAC-seq peaks for primary HMEC cells" visibility=2' $MPEAK/Library_3_distal.bed > $MPEAK/Library_3_distal_head.bed
sed -i '1 i\track name=HMEC_primary description="ATAC-seq peaks for primary HMEC cells" visibility=2' $MPEAK/Library_3_all.bed > $MPEAK/Library_3_head.bed
sed -i '1 i\track name=HMEC_AI958_proximal description="proximal ATAC-seq peaks for AI958 cells" visibility=2' $MPEAK/Library_4_proximal.bed > $MPEAK/Library_4_proximal_head.bed
sed -i '1 i\track name=HMEC_AI958_distal description="distal ATAC-seq peaks for AI958 cells" visibility=2' $MPEAK/Library_4_distal.bed > $MPEAK/Library_4_distal_head.bed
sed -i '1 i\track name=HMEC_AI958 description="ATAC-seq peaks for AI958 cells" visibility=2' $MPEAK/Library_4_all.bed > $MPEAK/Library_4_head.bed
sed -i '1 i\track name=HMEC_AI263_proximal description="proximal ATAC-seq peaks for AI263 cells" visibility=2' $MPEAK/Library_5_proximal.bed > $MPEAK/Library_5_proximal_head.bed
sed -i '1 i\track name=HMEC_AI263_distal description="distal ATAC-seq peaks for AI263 cells" visibility=2' $MPEAK/Library_5_distal.bed > $MPEAK/Library_5_distal_head.bed
sed -i '1 i\track name=HMEC_AI263 description="ATAC-seq peaks for AI263 cells" visibility=2' $MPEAK/Library_5_all.bed > $MPEAK/Library_5_head.bed

#Run the script to generate plots.
Rscript peak_mutation_density_plots.R