kNN_calc_deltaFST <- function(pairFST, kNNscores, region.names, quantile=0.995){

XpairFST <- pairFST

ids <- which(knnw_scores > quantile(knnw_scores, 0.995))

#print(ids)

XpairFST <- pairFST[-ids,]

#calculate the medoid (without the ouliers) 

DIST     <- as.matrix(dist(XpairFST))

mean_dist <- apply(DIST,1,mean)

min_val <- min(mean_dist[mean_dist!=0])

id_min  <- which(mean_dist==min_val)[1]

MEDOID  <- XpairFST[id_min,]

#MEDOID 
deltaFST <- matrix(NaN, length(ids),3)

for(xx in 1:length(ids)){
deltaFST[xx,] <- pairFST[ids[xx],] - MEDOID
}

rownames(deltaFST) <- region.names[ids]
colnames(deltaFST) <- colnames(pairFST)

return(list(deltaFST=deltaFST, MEDOID=MEDOID))

}
