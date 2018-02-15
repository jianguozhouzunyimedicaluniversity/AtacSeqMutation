#PBS -l nodes=1:ppn=4
#PBS -l walltime=10:00:00 
root_dir = /fs/project/PAS1143/Michael/ATACRun_2018Jan/
chrom_sizes_dir = $HOME
picard_dir = $HOME
#Download
	module load python/3.6
	wget https://genome.med.nyu.edu/results/external/IARC/2017-11-06/fastq/Library_2_S2_L001_R1_001.fastq.gz
	wget https://genome.med.nyu.edu/results/external/IARC/2017-11-06/fastq/Library_2_S2_L001_R2_001.fastq.gz
#Remove adapters. Specify that reads must be longer than 2 bp.
	cutadapt -a TAAGGCGA -a TCGCCTTA -A TAAGGCGA -A TCGCCTTA -n 2 -o Library_2_R1_trimmed.fastq.gz -p Library_2_R2_trimmed.fastq.gz -m 2 Library_2_S2_L001_R1_001.fastq.gz Library_2_S2_L001_R2_001.fastq.gz > /fs/project/PAS1143/Michael/ATACRun_2018Jan/logfiles/Library_2_cutadapt.log
#Align.
	bowtie2 -x $chrom_sizes_dir/mm10/mm10 -1 <(gunzip -c Library_2_R1_trimmed.fastq.gz) -2 <(gunzip -c Library_2_R2_trimmed.fastq.gz)  -S $root_dir/aln/Library_2.bam -k 12\
		--end-to-end --sensitive -I 0 -X 1000 --no-discordant --fr --no-unal -p 16 &> $root_dir/logfiles/Library_2_alignment.log
	cd $root_dir/aln/
	samtools sort -o "Library_2_sorted.bam" -T "Library_2_tmp" "Library_2.bam" 
	samtools index "Library_2_sorted.bam"
	samtools idxstats "Library_2_sorted.bam" > $root_dir/logfiles/"Library_2_sorted_idxstats.log"
#Keep unique reads
    samtools view -bq 1 Library_2_sorted.bam > Library_2_mapq1.bam
	samtools sort -o "Library_2_mapq1_sorted.bam" -T "Library_2_tmp" "Library_2_mapq1.bam" 
	samtools index "Library_2_mapq1_sorted.bam"
	samtools idxstats "Library_2_mapq1_sorted.bam" > $root_dir/logfiles/"Library_2_mapq1_idxstats.log"
#Removing mitochondrial reads
	samtools view -b "Library_2_mapq1_sorted.bam" chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > "Library_2_nochrM.bam"
	samtools sort -o "Library_2_nochrM_sorted.bam" -T "Library_2_tmp" "Library_2_nochrM.bam" 
	samtools index "Library_2_nochrM_sorted.bam"
	samtools idxstats "Library_2_nochrM_sorted.bam" > $root_dir/logfiles/"Library_2_nochrM_idxstats.log"
#Marking duplicates
	java -jar -Xmx8g $picard_dir/picard.jar MarkDuplicates INPUT=Library_2_nochrM_sorted.bam OUTPUT=Library_2_prefiltered.bam METRICS_FILE=Library_2_metrics.txt\
	 USE_JDK_DEFLATER=true USE_JDK_INFLATER=true &> $root_dir/logfiles/Library_2_markduplicates.log
	samtools sort -o "Library_2_prefiltered_sorted.bam" -T "Library_2_tmp" "Library_2_prefiltered.bam" 
	samtools index "Library_2_prefiltered_sorted.bam"
	samtools idxstats "Library_2_prefiltered_sorted.bam" > $root_dir/logfiles/"Library_2_prefiltered_idxstats.log"
#Filtering reads
    samtools view -F 1540 -b "Library_2_prefiltered.bam" > "Library_2_filtered.bam"
	samtools sort -o "Library_2_filtered_sorted.bam" -T "Library_2_tmp" "Library_2_filtered.bam" 
	samtools index "Library_2_filtered_sorted.bam"
	samtools idxstats "Library_2_filtered_sorted.bam" > $root_dir/logfiles/"Library_2_filtered_idxstats.log"
#Generate WIG files.
	module load python/2.7
	gosr binbam -f 0 -n 1000 -t Library_2 Library_2_filtered.bam 50 Library_2 1> $root_dir/tracks/Library_2.wig 2> $root_dir/logfiles/Library_2_gosr.log
#Convert to BigWig.
	wigToBigWig $root_dir/tracks/Library_2.wig $chrom_sizes_dir/mm10.chrom.sizes $root_dir/tracks/Library_2.bw
