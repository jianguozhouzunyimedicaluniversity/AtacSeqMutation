#commands to create csv from gtf, run perl script, and sort the output of the perl script.
#We also remove all non-nuclear chromosomes, remove the header, and retain only the first three columns for TSS.

#sed -e 's/;/,/g' -e 's/ //g' -e 's/\t/,/g' gencode.v27.annotation.gtf > gencode.v27.annotation.csv
#./gtftobed_human.pl gencode.v27.annotation.csv TSS_hg38_unsorted.bed
#sort -k1,1V -k2,2n TSS_hg38_unsorted.bed > TSS_hg38_noheader.bed
#sed '1i\'"chr\tstart\tend\texon\texonnumber\tgeneid\ttranscriptid\tgenename\ttranscriptname\tstrand" TSS_hg38_noheader.bed > TSS_hg38.bed
awk '$1~/chr1|chr2|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr20|chr21|chr22|chrX|chrY/' > TSS_hg38_nuclear.bed
awk '{printf("%s\t%s\t%s\n", $1, $2, $3)}' TSS_hg38_noheader.bed > TSS_hg38_first3cols.bed