#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00 

#Download
	cd ~/atac_mutation
	wget https://www.encodeproject.org/files/ENCFF111NYT/@@download/ENCFF111NYT.bam
#Keep unique reads
    samtools view -bq 1 "ENCFF111NYT.bam" > "H3K9Me3_HMEC_2_mapq1.bam"    
#Keep unique reads(Index)
   samtools sort -o "H3K9Me3_HMEC_2_sorted.bam" -T "H3K9Me3_HMEC_2_tmp" "H3K9Me3_HMEC_2_mapq1.bam"
	samtools index "H3K9Me3_HMEC_2_sorted.bam"
#Removing mitochondrial reads
   samtools view -b "H3K9Me3_HMEC_2_sorted.bam" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > "H3K9Me3_HMEC_2_nochrM.bam"
#Removing mitochondrial reads (Index)
    samtools index "H3K9Me3_HMEC_2_nochrM.bam"
#Marking duplicates
	java -jar -Xmx8g $HOME/picard.jar MarkDuplicates INPUT="H3K9Me3_HMEC_2_nochrM.bam" OUTPUT="H3K9Me3_HMEC_2_prefiltered.bam" METRICS_FILE="H3K9Me3_HMEC_2_metrics.txt" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
#Filtering reads
    samtools view -F 1540 -b "H3K9Me3_HMEC_2_prefiltered.bam" > "H3K9Me3_HMEC_2_filtered.bam"
#Filtering reads (Index)
    samtools index "H3K9Me3_HMEC_2_filtered.bam"
#Generate WIG files.
	module load python/2.7
	gosr binbam -f 0 -n 1000 -t H3K9Me3_HMEC_2 H3K9Me3_HMEC_2_filtered.bam 50 H3K9Me3_HMEC_2 > H3K9Me3_HMEC_2.wig
#Convert to BigWig.
	wigToBigWig H3K9Me3_HMEC_2.wig $HOME/hg38.chrom.sizes H3K9Me3_HMEC_2.bw
