calcSIL <- function(x11,x22){

res <- numeric(length(x11)*length(x22))
count <- 0
for(xx in 1:length(x11)){ 
  for(yy in 1:length(x22)){
  x1 <- x11[xx]
  x2 <- x22[yy]
  count <- count + 1
  res[count] <- (x1 - x2) ^ 2 #sqrt(sum((x1 - x2) ^ 2))
  }
}

B    <- mean(res)

A1   <- numeric(choose(length(x11),2))
count <- 0
for(xx in 1:(length(x11)-1)){ 
  for(yy in xx:(length(x11))){
  x1 <- x11[xx]
  x2 <- x11[yy]
  count <- count + 1
  A1[count] <- (x1 - x2) ^ 2 #sqrt(sum((x1 - x2) ^ 2))
  }
}

x22 <- c(x11,x22)
A2   <- numeric(choose(length(x22),2))
count <- 0
for(xx in 1:(length(x22)-1)){ 
  for(yy in xx:(length(x22))){
  x1 <- x22[xx]
  x2 <- x22[yy]
  count <- count + 1
  A2[count] <- (x1 - x2) ^ 2 #sqrt(sum((x1 - x2) ^ 2))
  }
}

A  <- mean(c(A1,A2))
#A <- mean(A1)

return((B-A)/(max(A,B)))
return(A)

}
#########################################################################

euc.dist  <- function(x11,x22){

#vvv <- c(x11,x22)
#x11 <- vvv
#x22 <- vvv

#return((var(x11)+var(x22)))
res <- numeric(length(x11)*length(x22))
count <- 0
for(xx in 1:length(x11)){ 
  for(yy in 1:length(x22)){
  x1 <- x11[xx]
  x2 <- x22[yy]
  count <- count + 1
  res[count] <- (x1 - x2) ^ 2 #sqrt(sum((x1 - x2) ^ 2))
  }
}
return(mean(res))#/(var(x11)+var(x22)))
}

######################################################################

calcK <- function(folder, fst){

sel.ids <- c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,
40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,
76,78,80,82,84,86,88,90,92,94,96,98,100)

#euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

k <- seq(10,990,by=10)
res <- numeric(length(k))
count <- 1

for (xx in k){

tt <- as.character(read.table(paste(folder,"/Elki_k",xx,"/cluster.txt",sep=""))[,5]) #5#3
tt <- strsplit(tt,'=')
tt <- sapply(tt,function(x){return(x[2])})
tt <- as.numeric(tt)

aa <- sort(tt, decreasing=TRUE)[51:1000]
bb <- sort(tt, decreasing=TRUE)[1:50]
#aa  <- tt[-sel.ids]
#bb  <- tt[sel.ids]
#P  <- table(tt)/sum(table(tt))
#V  <- as.numeric(names(table(tt)))
#E  <- sum(log(V)+log(P))
res[count] <- euc.dist(aa,bb)
#sum(sapply(split(c(aa,bb),1:k),function(x){euc.dist(x,x)}))
#euc.dist(aa,bb) 
#sum((aa-bb)^2)/sd(tt) 
#euc.dist(aa,bb)
#euc.dist(aa,bb) 
#euc.dist(aa,bb)/sqrt(var(aa)+var(bb)) 
#euc.dist(aa,bb) 
#sqrt(euc.dist(tt,tt)/sd(tt))
#calcSIL(aa,bb)
#euc.dist(aa,bb)
#sum((aa-mean(aa))^2) + sum((bb-mean(bb))^2) 
#euc.dist(mean(tt[-sel.ids]), mean(sort(tt, decreasing=TRUE)[1:50])) 
#euc.dist(max(tt[-sel.ids]), max(sort(tt, decreasing=TRUE)[1:50])) 
#euc.dist(min(tt), min(sort(tt, decreasing=TRUE)[1:50])) 
#euc.dist(min(tt), max(sort(tt, decreasing=TRUE)[1:50]))
#euc.dist(aa,bb)
#(median(aa) - mean(c(bb,aa)))^2
count <- count + 1
}

return(res)

}

readElkiOUT <- function(folder, pred){

require(pROC)

k <- seq(10,length(pred)-10,by=10)

#selection
#sel.ids <- 950:1000
#introgression
#sel.ids <- c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99)
# introgression 
#sel.ids <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64, 66, 68, 70, 72, 74, 76, 78, 80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100)

#group <- group <- sort(rep(1:1000,50))


ROCvalues <- rep(NaN,length(k))
#vec       <- numeric(1000)
#vec[sel.ids] <- 1
vec <- pred

count <- 1
for (xx in k){

#tt <- as.character(read.table(paste(folder,"/Elki_k",xx,"/cluster.txt",sep=""))[,5])  
tt <- as.character(read.table(paste(folder,"/Elki_k",xx,"/cluster.txt",sep=""))[,5])  
#5#3#8
#print(tt)
tt <- strsplit(tt,'=')
tt <- sapply(tt,function(x){return(x[2])})
tt <- as.numeric(tt)
#tt <- tapply(tt,group,sum)
#print(tt)
ROCvalues[count] <- auc(vec,tt)[1]
count <- count + 1
}

return(ROCvalues)
}

readElkiOUTALL <- function(ElkiIN, methods){

null  <- rep(NaN,99)

LOF <- null
KNN <- null
LOOP <- null
FASTABOD <- null
LDOF <- null
INFLO <- null
COF <- null
ODIN <- null
KNNW <- null
SIMPLIFIEDLOF <- null
LDF <- null

if(is.element("ALL", methods)){
methods <- c("LOF","KNN","LOOP","FASTABOD","LDOF","INFLO","COF","ODIN","KNNW","SIMPLIFIEDLOF","LDF")
}

#LOF
print("LOF")
if(is.element("LOF", methods)){
folder <- paste("LOF",ElkiIN, sep="")
LOF <- readElkiOUT(folder)
}
#KNN
print("KNN")
if(is.element("KNN", methods)){
folder <- paste("KNN",ElkiIN, sep="")
KNN <- readElkiOUT(folder)
}
#LoOP
print("LOOP")
if(is.element("LOOP", methods)){
folder <- paste("LOOP",ElkiIN, sep="")
LOOP <- readElkiOUT(folder)
}
#FastABOD
print("FASTABOD")
if(is.element("FASTABOD", methods)){
folder <- paste("FASTABOD",ElkiIN, sep="")
FASTABOD <- readElkiOUT(folder)
}
#LDOF
print("LDOF")
if(is.element("LDOF", methods)){
folder <- paste("LDOF",ElkiIN, sep="")
LDOF <- readElkiOUT(folder)
}
#INFLO
print("INFLO")
if(is.element("INFLO", methods)){
folder <- paste("INFLO",ElkiIN, sep="")
INFLO <- readElkiOUT(folder)
}
#COF
print("COF")
if(is.element("COF", methods)){
folder <- paste("COF",ElkiIN, sep="")
COF <- readElkiOUT(folder)
}
#ODIN
print("ODIN")
if(is.element("ODIN", methods)){
folder <- paste("ODIN",ElkiIN, sep="")
ODIN <- readElkiOUT(folder)
}
#KNNW
print("KNNW")
if(is.element("KNNW", methods)){
folder <- paste("KNNW",ElkiIN, sep="")
KNNW <- readElkiOUT(folder)
}
#SimplifiedLOF
print("SIMPLIFIEDLOF")
if(is.element("SIMPLIFIEDLOF", methods)){
folder <- paste("SIMPLIFIEDLOF",ElkiIN, sep="")
SIMPLIFIEDLOF <- readElkiOUT(folder)
}
#LDF
print("LDF")
if(is.element("LDF", methods)){
folder <- paste("LDF",ElkiIN, sep="")
LDF <- readElkiOUT(folder)
}

k <- rep(seq(10,990,by=10),11)
methods <- c(rep("LOF",99),rep("KNN",99),rep("LOOP",99),rep("FASTABOD",99),rep("LDOF",99),rep("INFLO",99),rep("COF",99),rep("ODIN",99),rep("KNNW",99),rep("SIMPLIFIEDLOF",99),rep("LDF",99))
AUC <- c(LOF,KNN,LOOP,FASTABOD,LDOF,INFLO,COF,ODIN,KNNW,SIMPLIFIEDLOF,LDF)

DATA <- cbind(methods,k,AUC)
DATA <- as.data.frame(DATA)   
DATA[,2] <- as.numeric(k)
DATA[,3] <- as.numeric(AUC)

# max mean sqrt(n)=30
MAX  <- c(max(LOF),max(KNN),max(LOOP),max(FASTABOD),max(LDOF),max(INFLO),max(COF),max(ODIN),max(KNNW),max(SIMPLIFIEDLOF),max(LDF))
KMAX <-  c(k[which.max(LOF)],k[which.max(KNN)],k[which.max(LOOP)],k[which.max(FASTABOD)],k[which.max(LDOF)],k[which.max(INFLO)],k[which.max(COF)],k[which.max(ODIN)],k[which.max(KNNW)],k[which.max(SIMPLIFIEDLOF)],k[which.max(LDF)])
MEAN <- c(mean(LOF),mean(KNN),mean(LOOP),mean(FASTABOD),mean(LDOF),mean(INFLO),mean(COF),mean(ODIN),mean(KNNW),mean(SIMPLIFIEDLOF),mean(LDF))
MIN  <- c(min(LOF),min(KNN),min(LOOP),min(FASTABOD),min(LDOF),min(INFLO),min(COF),min(ODIN),min(KNNW),min(SIMPLIFIEDLOF),min(LDF))
KMIN <-  c(k[which.min(LOF)],k[which.min(KNN)],k[which.min(LOOP)],k[which.min(FASTABOD)],k[which.min(LDOF)],k[which.min(INFLO)],k[which.min(COF)],k[which.min(ODIN)],k[which.min(KNNW)],k[which.min(SIMPLIFIEDLOF)],k[which.min(LDF)])

MAXKMEAN <- cbind(MAX,KMAX,MIN,KMIN,MEAN)
colnames(MAXKMEAN) <- c("max","kmax","min","kmin","mean")
rownames(MAXKMEAN) <- c("LOF","KNN","LOOP","FASTABOD","LDOF","INFLO","COF","ODIN","KNNW","SIMPLIFIEDLOF","LDF")

return(list(DATA=DATA,MAXKMEAN=MAXKMEAN))

}
