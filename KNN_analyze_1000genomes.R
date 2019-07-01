library(PopGenome)
library(readxl)
source("kNNCallElkiALL.R")
source("kNN_tau_window.R")
source("kNN_calc_best_k.R")

# define the samples 
samples <- as.data.frame(read_xlsx("20130606_sample_info.xlsx"))
CEU     <- samples[which(samples[,3]=="CEU"),1]
CHB     <- samples[which(samples[,3]=="CHB"),1]
YRI     <- samples[which(samples[,3]=="YRI"),1]

genome = readVCF("ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", 1000, "2", 1, 240000000,samplenames=c(CEU,CHB,YRI), include.unknown=TRUE)

#genome = readVCF("ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", 100, "6", 1, 170000000,samplenames=c(CEU,CHB,YRI), include.unknown=TRUE)

#genome = readVCF("ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", 100, "6", 1, 100000,samplenames=c(CEU,CHB,YRI), include.unknown=TRUE)

# set the populations 
genome = set.populations(genome,list(CEU,CHB,YRI),diploid=TRUE)

# Split the data into 50kb windows 
slide = sliding.window.transform(genome,jump=100000,width=100000, type=2)

# Get the genomic positions
genome.pos = sapply(slide@region.names,function(x){
  split = strsplit(x," ")[[1]][c(1,3)]
  val = mean(as.numeric(split))
  return(val)
})

# pairwise FST 
slide <- F_ST.stats(slide, mode="nucleotide")
# Get pairwise FST 
pairFST  <- t(slide@nuc.F_ST.pairwise)

#### KNN-based methods ################################
pairFST <- round(pairFST, digits=5)
pairFST[pairFST<0] <- 0
pairFST[is.na(pairFST)] <- 0

#Write output 
write.table(pairFST, file="ELKI-IN", col.names=FALSE, row.names=FALSE)

# perform kNN 
kNNCallElkiALL("ELKI-IN", "LOF",n.regions=length(slide@region.names))
kNNCallElkiALL("ELKI-IN", "KNN",n.regions=length(slide@region.names))
kNNCallElkiALL("ELKI-IN", "LOOP",n.regions=length(slide@region.names))
kNNCallElkiALL("ELKI-IN", "INFLO",n.regions=length(slide@region.names))
kNNCallElkiALL("ELKI-IN", "ODIN",n.regions=length(slide@region.names))
kNNCallElkiALL("ELKI-IN", "KNNW",n.regions=length(slide@region.names))
kNNCallElkiALL("ELKI-IN", "SIMPLIFIEDLOF",n.regions=length(slide@region.names))
kNNCallElkiALL("ELKI-IN", "LDF", n.regions=length(slide@region.names))

# calculate tau windows 
tau_lof   <- kNN_tau_window("LOFELKI-IN")
tau_knn   <- kNN_tau_window("KNNELKI-IN")
tau_loop  <- kNN_tau_window("LOOPELKI-IN")
tau_inflo <- kNN_tau_window("INFLOELKI-IN")
tau_odin  <- kNN_tau_window("ODINELKI-IN")
tau_knnw  <- kNN_tau_window("KNNWELKI-IN")
tau_simplifiedlof  <- kNN_tau_window("SIMPLIFIEDLOFELKI-IN")
tau_ldf   <- kNN_tau_window("LDFELKI-IN")

# plot diagnostic plot 
plot(tau_lof, ylim=c(0,1), col="black", xaxt="n", ylab="tau", pch=19, type="b", xlab="k")
points(tau_knn, col="grey", xaxt="n", pch=19, type="b")
points(tau_loop, col="green", xaxt="n", pch=19, type="b")
points(tau_inflo, col="orange", xaxt="n", pch=19, type="b")
points(tau_odin, col="purple",  xaxt="n", pch=19, type="b")
points(tau_knnw, col="yellow",  xaxt="n", pch=19, type="b")
points(tau_simplifiedlof, col="dark green",  xaxt="n", pch=19, type="b")
points(tau_ldf, col="red", xaxt="n", pch=19, type="b")
	legend("bottomleft",c("lof","knn","loop","inflo","odin","knnw","simplifiedlof","ldf"),
	fill=c("black","grey","green","orange","purple","yellow","dark green","red"))
axis(1,1:length(tau_knn),names(tau_knn), las=2)

# calculate the best k
k_lof   <- kNN_calc_best_k(tau_lof)
k_knn   <- kNN_calc_best_k(tau_knn)
k_loop  <- kNN_calc_best_k(tau_loop)
k_inflo <- kNN_calc_best_k(tau_inflo)
k_odin  <- kNN_calc_best_k(tau_odin)
k_knnw  <- kNN_calc_best_k(tau_knnw)
k_simplifiedlof <- kNN_calc_best_k(tau_simplifiedlof)
k_ldf   <- kNN_calc_best_k(tau_ldf)

# call the kNN methods with best k
kNNCallElkiSINGLE("ELKI-IN", "LOF",k_lof)
kNNCallElkiSINGLE("ELKI-IN", "KNN",k_knn)
kNNCallElkiSINGLE("ELKI-IN", "LOOP",k_loop)
kNNCallElkiSINGLE("ELKI-IN", "INFLO",k_inflo)
kNNCallElkiSINGLE("ELKI-IN", "ODIN",k_odin)
kNNCallElkiSINGLE("ELKI-IN", "KNNW",k_knnw)
kNNCallElkiSINGLE("ELKI-IN", "SIMPLIFIEDLOF",k_simplifiedlof)
kNNCallElkiSINGLE("ELKI-IN", "LDF", k_ldf)

# read out the kNN scores
lof_scores <- readElki(paste("LOFELKI-IN/Elki_k",k_lof,sep=""))
knn_scores <- readElki(paste("KNNELKI-IN/Elki_k",k_knn,sep=""))
loop_scores <- readElki(paste("LOOPELKI-IN/Elki_k",k_loop,sep=""))
inflo_scores <- readElki(paste("INFLOELKI-IN/Elki_k",k_inflo,sep=""))
odin_scores <- readElki(paste("ODINELKI-IN/Elki_k",k_odin,sep=""))
knnw_scores <- readElki(paste("KNNWELKI-IN/Elki_k",k_knnw,sep=""))
simplifiedlof_scores <- readElki(paste("SIMPLIFIEDLOFELKI-IN/Elki_k",k_simplifiedlof,sep=""))
ldf_scores <- readElki(paste("LDFELKI-IN/Elki_k",k_ldf,sep=""))

#plot the scores
par(mfrow=c(3,3))
plot(genome.pos, lof_scores, pch=19)
plot(genome.pos, knn_scores, pch=19)
plot(genome.pos, loop_scores, pch=19)
plot(genome.pos, inflo_scores, pch=19)
plot(genome.pos, odin_scores, pch=19)
plot(genome.pos, knnw_scores, pch=19)
plot(genome.pos, simplifiedlof_scores, pch=19)
plot(genome.pos, ldf_scores, pch=19)



























