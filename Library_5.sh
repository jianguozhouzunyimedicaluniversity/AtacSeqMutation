#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00 

#Download
	cd ~/atac_mutation/
	#wget https://genome.med.nyu.edu/results/external/IARC/2017-11-06/fastq/Library_5_S5_L001_R1_001.fastq.gz
	#wget https://genome.med.nyu.edu/results/external/IARC/2017-11-06/fastq/Library_5_S5_L001_R2_001.fastq.gz
#Unzip file.
	#gunzip Library_5_S5_L001_R1_001.fastq.gz
	#gunzip Library_5_S5_L001_R2_001.fastq.gz
#Remove adapters. Specify that reads must be longer than 2 bp.
	#cutadapt -a GCTACGCT -A GCTACGCT -o "Library_5_R1_trimmed.fastq" -p "Library_5_R2_trimmed.fastq" -m 2 "Library_5_S5_L001_R1_001.fastq" "Library_5_S5_L001_R2_001.fastq"
	cutadapt -a AGCGTAGC -A AGCGTAGC -o "Library_5_R1_trimmed2.fastq" -p "Library_5_R2_trimmed2.fastq" -m 2 "Library_5_R1_trimmed.fastq" "Library_5_R2_trimmed.fastq"
	
 bowtie2 -x $HOME/hg38/hg38 -1 "Library_5_R1_trimmed2.fastq" -2 "Library_5_R2_trimmed2.fastq" -S "Library_5.bam" -k 12 \
	--end-to-end --sensitive -I 0 -X 1000 --no-discordant --fr --no-unal -p 16
#Keep unique reads
    samtools view -bq 1 "Library_5.bam" > "Library_5_mapq1.bam"
#Sorting
    samtools sort -o "Library_5_sorted.bam" -T "Library_5_tmp" "Library_5_mapq1.bam" 
	samtools flagstat "Library_5_sorted.bam"   
#Keep unique reads(Index)
   samtools index "Library_5_sorted.bam"
#Removing mitochondrial reads
   samtools view -b "Library_5_sorted.bam" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > "Library_5_nochrM.bam"
#Removing mitochondrial reads (Index)
    samtools index "Library_5_nochrM.bam"
#Marking duplicates
	java -jar -Xmx8g $HOME/picard.jar MarkDuplicates INPUT="Library_5_nochrM.bam" OUTPUT="Library_5_prefiltered.bam" METRICS_FILE="Library_5_metrics.txt" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
#Filtering reads
    samtools view -F 1540 -b "Library_5_prefiltered.bam" > "Library_5_filtered.bam"
#Filtering reads (Index)
    samtools index "Library_5_filtered.bam"
#Generate WIG files.
	module load python/2.7
	gosr binbam -f 0 -n 1000 -t Library_5 Library_5_filtered.bam 50 Library_5 > Library_5.wig
#Convert to BigWig.
	wigToBigWig Library_5.wig $HOME/hg38.chrom.sizes Library_5.bw
