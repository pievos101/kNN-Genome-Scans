library(PopGenome)
#library(DDoutlier)
#library(Rlof)
library(pcadapt) 
library(BlockFeST)
library(pROC)
library(PRROC)

source("~/GitHub/kNN-Genome-Scans/simulations/convertTolfmm.R")
source("~/GitHub/kNN-Genome-Scans/simulations/R-Elki-Read-Out.R")
source("~/GitHub/kNN-Genome-Scans/simulations/R-Elki-Call.R")
source("~/GitHub/kNN-Genome-Scans/kNN_calc_best_k.R")
source("~/GitHub/kNN-Genome-Scans/kNN_tau_window.R")
source("~/GitHub/kNN-Genome-Scans/kNNCallElkiALL.R")

neut <- readMS("KNN_neut")
sel  <- readMS("KNN_sel")

data <- concatenate.classes(list(neut,sel))
data <- set.populations(data, list(1:100,101:200,201:300))
data <- F_ST.stats(data, mode="nucleotide")
#FST
FST  <- as.vector(data@nucleotide.F_ST)
#FST[FST<=0] <- 0
#AUC of FST
pred <- numeric(1000)
pred[951:1000] <- 1

nas <- which(sapply(data@region.data@biallelic.sites,length)==0)

if(length(nas)!=0){
pred <- pred[-nas]
FST  <- FST[-nas]
}

auc_FST <- auc(pred~FST)
auc_FST

pr_FST  <- pr.curve(scores.class0 = FST[pred==1], scores.class1 = FST[pred==0], curve = T)$auc.integral

pr_FST

#AUC of pcadapt 
gc()
convert_to_lfmm2(data)
# new pcadapt works 
filename <- read.pcadapt("LFMM-sel", type = "lfmm")
res <- pcadapt(filename, K=2, ploidy=1, min.maf = 0)
#create groups
GROUP <- sort(rep(1:1000,50))
#GROUP <- getBayes(data, snps=TRUE)$FUNC
P <- tapply(res$pvalue,GROUP,function(x){val<-log(x);sum(val[is.finite(val)])})

P <- as.numeric(P)

auc_PCADAPT <- auc(pred~P)

auc_PCADAPT

pr_PCADAPT  <- pr.curve(scores.class0 = abs(P)[pred==1], scores.class1 = abs(P)[pred==0], curve = T)$auc.integral

pr_PCADAPT

# FLK
convert_to_ped(data)
# call hapflk
# ./hapflk --file /home/bastian/kNN-project/REVISIONS/SINGLE-SNP_SIM/SNPS-0-7/FLK -K 1
# read hapflk results 
FLK <- read.table("hapflk.flk",stringsAsFactors=FALSE)[,6]
FLK <- as.numeric(FLK)
FLK <- FLK[-1]
FLK <- tapply(FLK,GROUP,function(x){val<-log(x);sum(val[is.finite(val)])})
auc_FLK <- auc(pred~FLK);auc_FLK

FLK <- as.numeric(FLK)

pr_FLK  <- pr.curve(scores.class0 = abs(FLK)[pred==1], scores.class1 =  abs(FLK)[pred==0], curve = T)$auc.integral

pr_FLK

## BlockFeST
#snps <- getBayes(data, snps=TRUE)
#BB   <- BlockFeST(snps, GROUP = snps$FUNC)


#### KNN-based methods ################################
# Get pairwise FST 
pairFSTx <- t(data@nuc.F_ST.pairwise)

pairFST  <- pairFSTx

if(length(nas)!=0){
pairFST  <- pairFSTx[-nas,]
}

pairFST  <- round(pairFST, digits=5)
pairFST[pairFST<0] <- 0
pairFST[is.na(pairFST)] <- 0

#Write output 
#1: 12
#2: 13
#3: 23
write.table(pairFST, file="ELKI-IN", col.names=FALSE, row.names=FALSE)

#methods <- c("LOF","KNN","LOOP","FASTABOD","LDOF","INFLO","COF","ODIN",
#"KNNW","SIMPLIFIEDLOF","LDF")
#}

CallElkiALL("ELKI-IN", "LOF")
CallElkiALL("ELKI-IN", "KNN")
CallElkiALL("ELKI-IN", "LOOP")
CallElkiALL("ELKI-IN", "INFLO")
CallElkiALL("ELKI-IN", "ODIN")
CallElkiALL("ELKI-IN", "KNNW")
CallElkiALL("ELKI-IN", "SIMPLIFIEDLOF")
CallElkiALL("ELKI-IN", "LDF")

lof   <- readElkiOUT("LOFELKI-IN", pred)
knn   <- readElkiOUT("KNNELKI-IN", pred)
loop  <- readElkiOUT("LOOPELKI-IN", pred)
inflo <- readElkiOUT("INFLOELKI-IN", pred)
odin  <- readElkiOUT("ODINELKI-IN", pred)
knnw  <- readElkiOUT("KNNWELKI-IN", pred)
simplifiedlof <- readElkiOUT("SIMPLIFIEDLOFELKI-IN", pred)
ldf   <- readElkiOUT("LDFELKI-IN", pred)

##plot the results
plot(lof,col="black", type="b", pch=1, ylim=c(0,1),xaxt="n", xlab="k", ylab="AUC")
lines(knn, col="grey", type="b", pch=2)
lines(loop, col="green", type="b", pch=3)
#lines(fastabod, col="blue", type="b", pch=4)
lines(inflo, col="orange", type="b", pch=5)
lines(odin, col="purple", type="b", pch=6)
lines(knnw, col="yellow", type="b", pch=7)
lines(simplifiedlof, col="dark green", type="b", pch=8)
lines(ldf, col="red", type="b", pch=9)
abline(h=auc_FST)
abline(h=auc_PCADAPT)
abline(h=auc_FLK)
breaks <- c(1,10,20,30,40,50,60,70,80,90,99)
axis(1,breaks,as.character(breaks*10))
#legend
legend("bottomleft",c("lof","knn","loop","inflo","odin","knnw","simplifiedlof","ldf"),
pch = c(1,2,3,5,6,7,8,9), 
col=c("black","grey","green","orange","purple","yellow","dark green","red"), 
lwd=c(1,1,1,1,1,1,1,1))

## Calculate the AUCs for kNN

# calculate tau windows 
tau_lof   <- kNN_tau_window("LOFELKI-IN")
tau_knn   <- kNN_tau_window("KNNELKI-IN")
tau_loop  <- kNN_tau_window("LOOPELKI-IN")
tau_inflo <- kNN_tau_window("INFLOELKI-IN")
tau_odin  <- kNN_tau_window("ODINELKI-IN")
tau_knnw  <- kNN_tau_window("KNNWELKI-IN")
tau_simplifiedlof  <- kNN_tau_window("SIMPLIFIEDLOFELKI-IN")
tau_ldf   <- kNN_tau_window("LDFELKI-IN")

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

lof_auc   <- auc(pred~lof_scores)
knn_auc   <- auc(pred~knn_scores)
loop_auc  <- auc(pred~loop_scores)
inflo_auc <- auc(pred~inflo_scores)
odin_auc  <- auc(pred~odin_scores)
knnw_auc  <- auc(pred~knnw_scores)
simplifiedlof_auc <- auc(pred~simplifiedlof_scores)
ldf_auc   <- auc(pred~ldf_scores)

RES        <- c(lof_auc,knn_auc,loop_auc,inflo_auc,odin_auc,knnw_auc,ldf_auc, auc_FST, auc_PCADAPT, auc_FLK)
names(RES) <- c("lof","knn","loop","inflo","odin","knnw","ldf","FST","PCADAPT","FLK")

barplot(RES, las=2, main="AUC")

# calculate PR values 

lof_pr <- pr.curve(scores.class0 = lof_scores[pred==1], scores.class1 =  lof_scores[pred==0], curve = T)$auc.integral
knn_pr <- pr.curve(scores.class0 = knn_scores[pred==1], scores.class1 =  knn_scores[pred==0], curve = T)$auc.integral
loop_pr <- pr.curve(scores.class0 = loop_scores[pred==1], scores.class1 =  loop_scores[pred==0], curve = T)$auc.integral
inflo_pr <- pr.curve(scores.class0 = inflo_scores[pred==1], scores.class1 =  inflo_scores[pred==0], curve = T)$auc.integral
odin_pr <- pr.curve(scores.class0 = -odin_scores[pred==1], scores.class1 =  -odin_scores[pred==0], curve = T)$auc.integral
knnw_pr <- pr.curve(scores.class0 = knnw_scores[pred==1], scores.class1 =  knnw_scores[pred==0], curve = T)$auc.integral
simplifiedlof_pr <- pr.curve(scores.class0 = simplifiedlof_scores[pred==1], scores.class1 =  simplifiedlof_scores[pred==0], curve = T)$auc.integral
ldf_pr <- pr.curve(scores.class0 = ldf_scores[pred==1], scores.class1 =  ldf_scores[pred==0], curve = T)$auc.integral


RES        <- c(lof_pr,knn_pr,loop_pr,inflo_pr,odin_pr,knnw_pr,ldf_pr, pr_FST, pr_PCADAPT, pr_FLK)
names(RES) <- c("lof","knn","loop","inflo","odin","knnw","ldf","FST","PCADAPT","FLK")

barplot(RES, las=2, main="PRAUC")













