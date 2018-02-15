#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00 

#Download
	cd ~/atac_mutation
	wget https://www.encodeproject.org/files/ENCFF887QAW/@@download/ENCFF887QAW.bam
#Keep unique reads
    samtools view -bq 1 "ENCFF887QAW.bam" > "DNase_HAEC_2_mapq1.bam"    
#Keep unique reads(Index)
   samtools sort -o "DNase_HAEC_2_sorted.bam" -T "DNase_HAEC_2_tmp" "DNase_HAEC_2_mapq1.bam"
	samtools index "DNase_HAEC_2_sorted.bam"
#Removing mitochondrial reads
   samtools view -b "DNase_HAEC_2_sorted.bam" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > "DNase_HAEC_2_nochrM.bam"
#Removing mitochondrial reads (Index)
    samtools index "DNase_HAEC_2_nochrM.bam"
#Marking duplicates
	java -jar -Xmx8g $HOME/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT="DNase_HAEC_2_nochrM.bam" OUTPUT="DNase_HAEC_2_prefiltered.bam" METRICS_FILE="DNase_HAEC_2_metrics.txt" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
#Filtering reads
    samtools view -F 1540 -b "DNase_HAEC_2_prefiltered.bam" > "DNase_HAEC_2_filtered.bam"
#Filtering reads (Index)
    samtools index "DNase_HAEC_2_filtered.bam"
#Generate WIG files.
	module load python/2.7
	gosr binbam -f 0 -n 1000 -t DNase_HAEC_2 DNase_HAEC_2_filtered.bam 50 DNase_HAEC_2 > DNase_HAEC_2.wig
#Convert to BigWig.
	wigToBigWig DNase_HAEC_2.wig $HOME/hg38.chrom.sizes DNase_HAEC_2.bw
