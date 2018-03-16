setwd("/fs/project/PAS1143/Michael")

#Mutation sets
Immortal_AG591=read.table("mutation/Immortal_AG591.bed", sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Immortal_AG591)=c("chr", "start", "stop", "mutation", "5", "6", "7", "8", "9")
Immortal_AG695=read.table("mutation/Immortal_AG695.bed", sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Immortal_AG695)=c("chr", "start", "stop", "mutation", "5", "6", "7", "8", "9")
Immortal_AH425=read.table("mutation/Immortal_AH425.bed", sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Immortal_AH425)=c("chr", "start", "stop", "mutation", "5", "6", "7", "8", "9")
Post_Stasis_AI263=read.table("mutation/Post_stasis_AI263.bed", sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Post_Stasis_AI263)=c("chr", "start", "stop", "mutation", "5", "6", "7", "8", "9")
Post_Stasis_AI958=read.table("mutation/Post_stasis_AI958.bed", sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Post_Stasis_AI958)=c("chr", "start", "stop", "mutation", "5", "6", "7", "8", "9")
Tumorigenic_AI027=read.table("mutation/Tumorigenic_AI027.bed", sep="\t", header = FALSE, stringsAsFactors = FALSE)
colnames(Tumorigenic_AI027)=c("chr", "start", "stop", "mutation", "5", "6", "7", "8", "9")

#Peak sets
dnase_dist_peaks=read.table("public_peaks/DNase_HAEC_distal.bed")
dnase_dist_peaks=dnase_dist_peaks[,c(1:3)]
colnames(dnase_dist_peaks)=c("chr", "start", "stop")

dnase_prox_peaks=read.table("public_peaks/DNase_HAEC_proximal.bed")
dnase_prox_peaks=dnase_prox_peaks[,c(1:3)]
colnames(dnase_prox_peaks)=c("chr", "start", "stop")

h3k9me3_dist_peaks=read.table("public_peaks/H3K9Me3_HMEC_distal.bed")
h3k9me3_dist_peaks=h3k9me3_dist_peaks[,c(1:3)]
colnames(h3k9me3_dist_peaks)=c("chr", "start", "stop")

h3k9me3_prox_peaks=read.table("public_peaks/H3K9Me3_HMEC_proximal.bed")
h3k9me3_prox_peaks=h3k9me3_prox_peaks[,c(1:3)]
colnames(h3k9me3_prox_peaks)=c("chr", "start", "stop")

hmec_prox_peaks=read.table("library_peaks/Library_3_proximal.bed")
hmec_prox_peaks=hmec_prox_peaks[,c(1:3)]
colnames(hmec_prox_peaks)=c("chr", "start", "stop")

hmec_dist_peaks=read.table("library_peaks/Library_3_distal.bed")
hmec_dist_peaks=hmec_dist_peaks[,c(1:3)]
colnames(hmec_dist_peaks)=c("chr", "start", "stop")

#Chromosome data frame
chroms=c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

#create a matrix that identifies the start and stop of every coordinate 
hotspots_hmec_prox = lapply(chroms, function(x){return(hmec_prox_peaks[hmec_prox_peaks$chr==x,])})
hotspots_hmec_dist=lapply(chroms, function(x){return(hmec_dist_peaks[hmec_dist_peaks$chr==x,])})
hotspots_dnase_prox=lapply(chroms, function(x){return(dnase_prox_peaks[dnase_prox_peaks$chr==x,])})
hotspots_dnase_dist=lapply(chroms, function(x){return(dnase_dist_peaks[dnase_dist_peaks$chr==x,])})
hotspots_h3k9me3_prox=lapply(chroms, function(x){return(h3k9me3_prox_peaks[h3k9me3_prox_peaks$chr==x,])})
hotspots_h3k9me3_dist=lapply(chroms, function(x){return(h3k9me3_dist_peaks[h3k9me3_dist_peaks$chr==x,])})
ag591=lapply(chroms, function(x){return(Immortal_AG591[Immortal_AG591$chr==x,])})
ag695=lapply(chroms, function(x){return(Immortal_AG695[Immortal_AG695$chr==x,])})
ah425=lapply(chroms, function(x){return(Immortal_AH425[Immortal_AH425$chr==x,])})
ai263=lapply(chroms, function(x){return(Post_Stasis_AI263[Post_Stasis_AI263$chr==x,])})
ai958=lapply(chroms, function(x){return(Post_Stasis_AI958[Post_Stasis_AI958$chr==x,])})
ai027=lapply(chroms, function(x){return(Tumorigenic_AI027[Tumorigenic_AI027$chr==x,])})
allstarts=lapply(1:length(chroms), function(x){
	return(c(hotspots_hmec_prox[[x]]$start, hotspots_hmec_dist[[x]]$start, hotspots_dnase_prox[[x]]$start, hotspots_dnase_dist[[x]]$start, hotspots_h3k9me3_prox[[x]]$start, hotspots_h3k9me3_dist[[x]]$start, ag591[[x]]$start, ag695[[x]]$start, ah425[[x]]$start, ai263[[x]]$start, ai958[[x]]$start, ai027[[x]]$start))})
allends=lapply(1:length(chroms), function(x){
	return(c(hotspots_hmec_prox[[x]]$stop, hotspots_hmec_dist[[x]]$stop, hotspots_dnase_prox[[x]]$stop, hotspots_dnase_dist[[x]]$stop, hotspots_h3k9me3_prox[[x]]$stop, hotspots_h3k9me3_dist[[x]]$stop, ag591[[x]]$stop, ag695[[x]]$stop, ah425[[x]]$stop, ai263[[x]]$stop, ai958[[x]]$stop, ai027[[x]]$stop))})
minstarts=lapply(1:length(allstarts), function(x){return(min(allstarts[[x]]))})
maxends=lapply(1:length(allstarts), function(x){return(max(allends[[x]]))})
chroms_df = as.data.frame(cbind(chroms, minstarts, maxends))
colnames(chroms_df)=c("chr", "start", "stop")
subsetallchromosomes=lapply(chroms, function(x){return(chroms_df[chroms_df$chr==x,])})

#find the maximum and minimum value of each coordinate 
minval=lapply(subsetallchromosomes, function(x){return(x[1,2])})
maxval=lapply(subsetallchromosomes, function(x){return(x[1,3])})
       
#divide in 1Mb bins 
divison=lapply(1:length(minval), function(x){return(seq(as.numeric(minval[[x]])-2000000,as.numeric(maxval[[x]])+2000000,1000000))})

#create a dataframe for the 1Mb bins
MB_bins_list=lapply(1:length(chroms), function(x){return(cbind(chroms[x], divison[[x]][1:length(divison[[x]])-1], divison[[x]][2:length(divison[[x]])]-1))})

#This actually creates the list of bins

findmutation<-function(data){  
    spont=lapply(chroms, function(x){return(data[data$chr==x,])})
    allmuts=lapply(spont, function(x){return(x$start)})
    listofmutsbychr=lapply(1:length(chroms), function(x){
				muts=matrix(data=0, nrow=nrow(MB_bins_list[[x]]), ncol=1)
				num = lapply(1:length(allmuts[[x]]), function(y){
					return(which(allmuts[[x]][y]>=as.numeric(as.character(MB_bins_list[[x]][,2])) & allmuts[[x]][y] <= as.numeric(as.character(MB_bins_list[[x]][,3]))))})
				for(i in num){
					muts[i,1] =  as.numeric(muts[i,1]+1)
				}
			return(as.data.frame(muts))})
	return(listofmutsbychr)
}

##place the mutations in the correct bin
Immortal_AG591_mut=findmutation(Immortal_AG591)
Immortal_AG695_mut=findmutation(Immortal_AG695) 
Immortal_AH425_mut=findmutation(Immortal_AH425) 
Post_Stasis_AI263_mut=findmutation(Post_Stasis_AI263) 
Post_Stasis_AI958_mut=findmutation(Post_Stasis_AI958)
Tumorigenic_AI027_mut=findmutation(Tumorigenic_AI027) 

#Convert chromosome-specific lists to matrices
MB_bins_matrix = do.call(rbind, MB_bins_list)
colnames(MB_bins_matrix)=c("chr", "start", "stop")
Immortal_AG591_mut_matrix=do.call(rbind, Immortal_AG591_mut)
Immortal_AG695_mut_matrix=do.call(rbind, Immortal_AG695_mut)
Immortal_AH425_mut_matrix=do.call(rbind, Immortal_AH425_mut)
Post_Stasis_AI263_mut_matrix=do.call(rbind, Post_Stasis_AI263_mut)
Post_Stasis_AI958_mut_matrix=do.call(rbind, Post_Stasis_AI958_mut)
Tumorigenic_AI027_mut_matrix=do.call(rbind, Tumorigenic_AI027_mut) 

#####DHS by number of overlapping peaks###########################################
if(!require(GenomicRanges)){library("GenomicRanges")}
MB_bins_df = as.data.frame(MB_bins_matrix)
MB_bins_df$start = as.numeric(as.character(MB_bins_df$start))
MB_bins_df$stop = as.numeric(as.character(MB_bins_df$stop))
MB_bins_granges <- with(MB_bins_df, GRanges(chr, IRanges(start, stop)))

#Proximal HMEC
prox_hmec_granges <- with(hmec_prox_peaks, GRanges(chr, IRanges(start, stop)))
overlaps_hmec_granges=findOverlaps(MB_bins_granges, prox_hmec_granges)
MB_bins_df=cbind(MB_bins_df, matrix(data=0, nrow=nrow(MB_bins_df), ncol=1))
MB_bins_df[,4] = as.numeric(as.character(unlist(lapply(1:nrow(MB_bins_df), function(x){return(length(which(queryHits(overlaps_hmec_granges)==x)))}))))
hmec_prox_overlapping_bins=as.matrix(MB_bins_df[,4])

#Distal HMEC
dist_hmec_granges <- with(hmec_dist_peaks, GRanges(chr, IRanges(start, stop)))
overlaps_hmec_granges=findOverlaps(MB_bins_granges, dist_hmec_granges)
MB_bins_df=cbind(MB_bins_df, matrix(data=0, nrow=nrow(MB_bins_df), ncol=1))
MB_bins_df[,5] = as.numeric(as.character(unlist(lapply(1:nrow(MB_bins_df), function(x){return(length(which(queryHits(overlaps_hmec_granges)==x)))}))))
hmec_dist_overlapping_bins=as.matrix(MB_bins_df[,5])

#Proximal dnase
prox_dnase_granges <- with(dnase_prox_peaks, GRanges(chr, IRanges(start, stop)))
overlaps_dnase_granges=findOverlaps(MB_bins_granges, prox_dnase_granges)
MB_bins_df=cbind(MB_bins_df, matrix(data=0, nrow=nrow(MB_bins_df), ncol=1))
MB_bins_df[,6] = as.numeric(as.character(unlist(lapply(1:nrow(MB_bins_df), function(x){return(length(which(queryHits(overlaps_dnase_granges)==x)))}))))
dnase_prox_overlapping_bins=as.matrix(MB_bins_df[,6])

#Distal dnase
dist_dnase_granges <- with(dnase_dist_peaks, GRanges(chr, IRanges(start, stop)))
overlaps_dnase_granges=findOverlaps(MB_bins_granges, dist_dnase_granges)
MB_bins_df=cbind(MB_bins_df, matrix(data=0, nrow=nrow(MB_bins_df), ncol=1))
MB_bins_df[,7] = as.numeric(as.character(unlist(lapply(1:nrow(MB_bins_df), function(x){return(length(which(queryHits(overlaps_dnase_granges)==x)))}))))
dnase_dist_overlapping_bins=as.matrix(MB_bins_df[,7])

#Proximal h3k9me3
prox_h3k9me3_granges <- with(h3k9me3_prox_peaks, GRanges(chr, IRanges(start, stop)))
overlaps_h3k9me3_granges=findOverlaps(MB_bins_granges, prox_h3k9me3_granges)
MB_bins_df=cbind(MB_bins_df, matrix(data=0, nrow=nrow(MB_bins_df), ncol=1))
MB_bins_df[,8] = as.numeric(as.character(unlist(lapply(1:nrow(MB_bins_df), function(x){return(length(which(queryHits(overlaps_h3k9me3_granges)==x)))}))))
h3k9me3_prox_overlapping_bins=as.matrix(MB_bins_df[,8])

#Distal h3k9me3
dist_h3k9me3_granges <- with(h3k9me3_dist_peaks, GRanges(chr, IRanges(start, stop)))
overlaps_h3k9me3_granges=findOverlaps(MB_bins_granges, dist_h3k9me3_granges)
MB_bins_df=cbind(MB_bins_df, matrix(data=0, nrow=nrow(MB_bins_df), ncol=1))
MB_bins_df[,9] = as.numeric(as.character(unlist(lapply(1:nrow(MB_bins_df), function(x){return(length(which(queryHits(overlaps_h3k9me3_granges)==x)))}))))
h3k9me3_dist_overlapping_bins=as.matrix(MB_bins_df[,9])

####################################################################################

#Create full matrix.
fullmatrix=cbind(MB_bins_df[,c(1:3)], hmec_prox_overlapping_bins, hmec_dist_overlapping_bins, dnase_prox_overlapping_bins, dnase_dist_overlapping_bins, h3k9me3_prox_overlapping_bins, h3k9me3_dist_overlapping_bins, Immortal_AG591_mut_matrix, Immortal_AG695_mut_matrix, Immortal_AH425_mut_matrix, Post_Stasis_AI263_mut_matrix, Post_Stasis_AI958_mut_matrix, Tumorigenic_AI027_mut_matrix)
colnames(fullmatrix)=c("chr", "start", "stop", "overlapping_hmec_prox", "overlapping_hmec_dist","overlapping_dnase_prox", "overlapping_dnase_dist", "overlapping_h3k9me3_prox", "overlapping_h3k9me3_dist","Immortal_AG591_mut_matrix", "Immortal_AG695_mut_matrix", "Immortal_AH425_mut_matrix", "Post_Stasis_AI263_mut_matrix", "Post_Stasis_AI958_mut_matrix", "Tumorigenic_AI027_mut_matrix")
#Save the plot for each mutation-cell line pair.
#Immortal_AG591
ratio = 3/4
nonzero = which(fullmatrix$overlapping_hmec_prox != 0 & fullmatrix$Immortal_AG591_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_prox[nonzero], fullmatrix$Immortal_AG591_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_prox_ag591.png")
plot(fullmatrix$overlapping_hmec_prox, fullmatrix$Immortal_AG591_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG591 Mutations in HMEC Primary Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_prox), y = ratio * max(fullmatrix$Immortal_AG591_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_hmec_dist != 0 & fullmatrix$Immortal_AG591_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_dist[nonzero], fullmatrix$Immortal_AG591_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_dist_ag591.png")
plot(fullmatrix$overlapping_hmec_dist, fullmatrix$Immortal_AG591_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG591 Mutations in HMEC Primary Distal Peaks", col="blue")
text(labels = txt,x = ratio * max(fullmatrix$overlapping_hmec_dist), y = ratio * max(fullmatrix$Immortal_AG591_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_prox != 0 & fullmatrix$Immortal_AG591_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_prox[nonzero], fullmatrix$Immortal_AG591_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_prox_ag591.png")
plot(fullmatrix$overlapping_dnase_prox, fullmatrix$Immortal_AG591_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG591 Mutations in HAEC DNase Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_prox), y = ratio * max(fullmatrix$Immortal_AG591_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_dist != 0 & fullmatrix$Immortal_AG591_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_dist[nonzero], fullmatrix$Immortal_AG591_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_dist_ag591.png")
plot(fullmatrix$overlapping_dnase_dist, fullmatrix$Immortal_AG591_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG591 Mutations in HAEC DNase Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_dist), y = ratio * max(fullmatrix$Immortal_AG591_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_prox != 0 & fullmatrix$Immortal_AG591_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_prox[nonzero], fullmatrix$Immortal_AG591_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_prox_ag591.png")
plot(fullmatrix$overlapping_h3k9me3_prox, fullmatrix$Immortal_AG591_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG591 Mutations in HAEC H3K9Me3 Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_prox), y = ratio * max(fullmatrix$Immortal_AG591_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_dist != 0 & fullmatrix$Immortal_AG591_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_dist[nonzero], fullmatrix$Immortal_AG591_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_dist_ag591.png")
plot(fullmatrix$overlapping_h3k9me3_dist, fullmatrix$Immortal_AG591_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG591 Mutations in HAEC H3K9Me3 Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_dist), y = ratio * max(fullmatrix$Immortal_AG591_mut_matrix))
dev.off()

#Immortal AG695
nonzero = which(fullmatrix$overlapping_hmec_prox != 0 & fullmatrix$Immortal_AG695_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_prox[nonzero], fullmatrix$Immortal_AG695_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_prox_ag695.png")
plot(fullmatrix$overlapping_hmec_prox, fullmatrix$Immortal_AG695_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG695 Mutations in HMEC Primary Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_prox), y = ratio * max(fullmatrix$Immortal_AG695_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_hmec_dist != 0 & fullmatrix$Immortal_AG695_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_dist[nonzero], fullmatrix$Immortal_AG695_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_dist_ag695.png")
plot(fullmatrix$overlapping_hmec_dist, fullmatrix$Immortal_AG695_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG695 Mutations in HMEC Primary Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_dist), y = ratio * max(fullmatrix$Immortal_AG695_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_prox != 0 & fullmatrix$Immortal_AG695_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_prox[nonzero], fullmatrix$Immortal_AG695_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_prox_ag695.png")
plot(fullmatrix$overlapping_dnase_prox, fullmatrix$Immortal_AG695_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG695 Mutations in HAEC DNase Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_prox), y = ratio * max(fullmatrix$Immortal_AG695_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_dist != 0 & fullmatrix$Immortal_AG695_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_dist[nonzero], fullmatrix$Immortal_AG695_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_dist_ag695.png")
plot(fullmatrix$overlapping_dnase_dist, fullmatrix$Immortal_AG695_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG695 Mutations in HAEC DNase Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_dist), y = ratio * max(fullmatrix$Immortal_AG695_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_prox != 0 & fullmatrix$Immortal_AG695_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_prox[nonzero], fullmatrix$Immortal_AG695_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_prox_ag695.png")
plot(fullmatrix$overlapping_h3k9me3_prox, fullmatrix$Immortal_AG695_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG695 Mutations in HAEC H3K9Me3 Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_prox), y = ratio * max(fullmatrix$Immortal_AG695_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_dist != 0 & fullmatrix$Immortal_AG695_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_dist[nonzero], fullmatrix$Immortal_AG695_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_dist_ag695.png")
plot(fullmatrix$overlapping_h3k9me3_dist, fullmatrix$Immortal_AG695_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AG695 Mutations in HAEC H3K9Me3 Distal Peaks", col="blue")
text(labels = txt,  x = ratio * max(fullmatrix$overlapping_h3k9me3_dist), y = ratio * max(fullmatrix$Immortal_AG695_mut_matrix))
dev.off()

#Immortal AH425
nonzero = which(fullmatrix$overlapping_hmec_prox != 0 & fullmatrix$Immortal_AH425_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_prox[nonzero], fullmatrix$Immortal_AH425_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_prox_ah425.png")
plot(fullmatrix$overlapping_hmec_prox, fullmatrix$Immortal_AH425_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AH425 Mutations in HMEC Primary Proximal Peaks", col="blue")
text(labels = txt,  x = ratio * max(fullmatrix$overlapping_hmec_prox), y = ratio * max(fullmatrix$Immortal_AH425_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_hmec_dist != 0 & fullmatrix$Immortal_AH425_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_dist[nonzero], fullmatrix$Immortal_AH425_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_dist_ah425.png")
plot(fullmatrix$overlapping_hmec_dist, fullmatrix$Immortal_AH425_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AH425 Mutations in HMEC Primary Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_dist), y = ratio * max(fullmatrix$Immortal_AH425_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_prox != 0 & fullmatrix$Immortal_AH425_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_prox[nonzero], fullmatrix$Immortal_AH425_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_prox_ah425.png")
plot(fullmatrix$overlapping_dnase_prox, fullmatrix$Immortal_AH425_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AH425 Mutations in HAEC DNase Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_prox), y = ratio * max(fullmatrix$Immortal_AH425_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_dist != 0 & fullmatrix$Immortal_AH425_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_dist[nonzero], fullmatrix$Immortal_AH425_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_dist_ah425.png")
plot(fullmatrix$overlapping_dnase_dist, fullmatrix$Immortal_AH425_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AH425 Mutations in HAEC DNase Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_dist), y = ratio * max(fullmatrix$Immortal_AH425_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_prox != 0 & fullmatrix$Immortal_AH425_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_prox[nonzero], fullmatrix$Immortal_AH425_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_prox_ah425.png")
plot(fullmatrix$overlapping_h3k9me3_prox, fullmatrix$Immortal_AH425_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AH425 Mutations in HAEC H3K9Me3 Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_prox), y = ratio * max(fullmatrix$Immortal_AH425_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_dist != 0 & fullmatrix$Immortal_AH425_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_dist[nonzero], fullmatrix$Immortal_AH425_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_dist_ah425.png")
plot(fullmatrix$overlapping_h3k9me3_dist, fullmatrix$Immortal_AH425_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Immortal AH425 Mutations in HAEC H3K9Me3 Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_dist), y = ratio * max(fullmatrix$Immortal_AH425_mut_matrix))
dev.off()

#Post-Stasis AI263
nonzero = which(fullmatrix$overlapping_hmec_prox != 0 & fullmatrix$Post_Stasis_AI263_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_prox[nonzero], fullmatrix$Post_Stasis_AI263_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_prox_ai263.png")
plot(fullmatrix$overlapping_hmec_prox, fullmatrix$Post_Stasis_AI263_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI263 Mutations in HMEC Primary Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_prox), y = ratio * max(fullmatrix$Post_Stasis_AI263_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_hmec_dist != 0 & fullmatrix$Post_Stasis_AI263_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_dist[nonzero], fullmatrix$Post_Stasis_AI263_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_dist_ai263.png")
plot(fullmatrix$overlapping_hmec_dist, fullmatrix$Post_Stasis_AI263_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI263 Mutations in HMEC Primary Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_dist), y = ratio * max(fullmatrix$Post_Stasis_AI263_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_prox != 0 & fullmatrix$Post_Stasis_AI263_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_prox[nonzero], fullmatrix$Post_Stasis_AI263_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_prox_ai263.png")
plot(fullmatrix$overlapping_dnase_prox, fullmatrix$Post_Stasis_AI263_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI263 Mutations in HAEC DNase Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_prox), y = ratio * max(fullmatrix$Post_Stasis_AI263_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_dist != 0 & fullmatrix$Post_Stasis_AI263_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_dist[nonzero], fullmatrix$Post_Stasis_AI263_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_dist_ai263.png")
plot(fullmatrix$overlapping_dnase_dist, fullmatrix$Post_Stasis_AI263_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI263 Mutations in HAEC DNase Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_dist), y = ratio * max(fullmatrix$Post_Stasis_AI263_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_prox != 0 & fullmatrix$Post_Stasis_AI263_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_prox[nonzero], fullmatrix$Post_Stasis_AI263_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_prox_ai263.png")
plot(fullmatrix$overlapping_h3k9me3_prox, fullmatrix$Post_Stasis_AI263_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI263 Mutations in HAEC H3K9Me3 Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_prox), y = ratio * max(fullmatrix$Post_Stasis_AI263_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_dist != 0 & fullmatrix$Post_Stasis_AI263_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_dist[nonzero], fullmatrix$Post_Stasis_AI263_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_dist_ai263.png")
plot(fullmatrix$overlapping_h3k9me3_dist, fullmatrix$Post_Stasis_AI263_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI263 Mutations in HAEC H3K9Me3 Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_dist), y = ratio * max(fullmatrix$Post_Stasis_AI263_mut_matrix))
dev.off()

#Post-Stasis AI958
nonzero = which(fullmatrix$overlapping_hmec_prox != 0 & fullmatrix$Post_Stasis_AI958_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_prox[nonzero], fullmatrix$Post_Stasis_AI958_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_prox_ai958.png")
plot(fullmatrix$overlapping_hmec_prox, fullmatrix$Post_Stasis_AI958_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI958 Mutations in HMEC Primary Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_prox), y = ratio * max(fullmatrix$Post_Stasis_AI958_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_hmec_dist != 0 & fullmatrix$Post_Stasis_AI958_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_dist[nonzero], fullmatrix$Post_Stasis_AI958_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_dist_ai958.png")
plot(fullmatrix$overlapping_hmec_dist, fullmatrix$Post_Stasis_AI958_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI958 Mutations in HMEC Primary Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_dist), y = ratio * max(fullmatrix$Post_Stasis_AI958_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_prox != 0 & fullmatrix$Post_Stasis_AI958_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_prox[nonzero], fullmatrix$Post_Stasis_AI958_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_prox_ai958.png")
plot(fullmatrix$overlapping_dnase_prox, fullmatrix$Post_Stasis_AI958_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI958 Mutations in HAEC DNase Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_prox), y = ratio * max(fullmatrix$Post_Stasis_AI958_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_dist != 0 & fullmatrix$Post_Stasis_AI958_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_dist[nonzero], fullmatrix$Post_Stasis_AI958_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_dist_ai958.png")
plot(fullmatrix$overlapping_dnase_dist, fullmatrix$Post_Stasis_AI958_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI958 Mutations in HAEC DNase Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_dist), y = ratio * max(fullmatrix$Post_Stasis_AI958_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_prox != 0 & fullmatrix$Post_Stasis_AI958_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_prox[nonzero], fullmatrix$Post_Stasis_AI958_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_prox_ai958.png")
plot(fullmatrix$overlapping_h3k9me3_prox, fullmatrix$Post_Stasis_AI958_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI958 Mutations in HAEC H3K9Me3 Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_prox), y = ratio * max(fullmatrix$Post_Stasis_AI958_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_dist != 0 & fullmatrix$Post_Stasis_AI958_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_dist[nonzero], fullmatrix$Post_Stasis_AI958_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_dist_ai958.png")
plot(fullmatrix$overlapping_h3k9me3_dist, fullmatrix$Post_Stasis_AI958_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Post-Stasis AI958 Mutations in HAEC H3K9Me3 Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_dist), y = ratio * max(fullmatrix$Post_Stasis_AI958_mut_matrix))
dev.off()

#Tumorigenic AI027
nonzero = which(fullmatrix$overlapping_hmec_prox != 0 & fullmatrix$Tumorigenic_AI027_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_prox[nonzero], fullmatrix$Tumorigenic_AI027_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_dist_ai027.png")
plot(fullmatrix$overlapping_hmec_prox, fullmatrix$Tumorigenic_AI027_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Tumorigenic AI027 Mutations in HMEC Primary Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_dist), y = ratio * max(fullmatrix$Tumorigenic_AI027_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_hmec_dist != 0 & fullmatrix$Tumorigenic_AI027_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_hmec_dist[nonzero], fullmatrix$Tumorigenic_AI027_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/hmec_dist_ai027.png")
plot(fullmatrix$overlapping_hmec_dist, fullmatrix$Tumorigenic_AI027_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Tumorigenic AI027 Mutations in HMEC Primary Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_hmec_dist), y = ratio * max(fullmatrix$Tumorigenic_AI027_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_prox != 0 & fullmatrix$Tumorigenic_AI027_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_prox[nonzero], fullmatrix$Tumorigenic_AI027_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_prox_ai027.png")
plot(fullmatrix$overlapping_dnase_prox, fullmatrix$Tumorigenic_AI027_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Tumorigenic AI027 Mutations in HAEC DNase Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_prox), y = ratio * max(fullmatrix$Tumorigenic_AI027_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_dnase_dist != 0 & fullmatrix$Tumorigenic_AI027_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_dnase_dist[nonzero], fullmatrix$Tumorigenic_AI027_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/dnase_dist_ai027.png")
plot(fullmatrix$overlapping_dnase_dist, fullmatrix$Tumorigenic_AI027_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Tumorigenic AI027 Mutations in HAEC DNase Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_dnase_dist), y = ratio * max(fullmatrix$Tumorigenic_AI027_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_prox != 0 & fullmatrix$Tumorigenic_AI027_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_prox[nonzero], fullmatrix$Tumorigenic_AI027_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_prox_ai027.png")
plot(fullmatrix$overlapping_h3k9me3_prox, fullmatrix$Tumorigenic_AI027_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Tumorigenic AI027 Mutations in HAEC H3K9Me3 Proximal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_prox), y = ratio * max(fullmatrix$Tumorigenic_AI027_mut_matrix))
dev.off()

nonzero = which(fullmatrix$overlapping_h3k9me3_dist != 0 & fullmatrix$Tumorigenic_AI027_mut_matrix != 0)
test = cor.test(fullmatrix$overlapping_h3k9me3_dist[nonzero], fullmatrix$Tumorigenic_AI027_mut_matrix[nonzero], alternative = "greater", method = "spearman")
rho = test$estimate
txt = paste("rho = ", rho)
png("/fs/project/PAS1143/Michael/figs/h3k9me3_dist_ai027.png")
plot(fullmatrix$overlapping_h3k9me3_dist, fullmatrix$Tumorigenic_AI027_mut_matrix, xlab="Read Counts", ylab="Number of Mutations", pch=20, main="Tumorigenic AI027 Mutations in HAEC H3K9Me3 Distal Peaks", col="blue")
text(labels = txt, x = ratio * max(fullmatrix$overlapping_h3k9me3_dist), y = ratio * max(fullmatrix$Tumorigenic_AI027_mut_matrix))
dev.off()
