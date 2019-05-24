##Subfunctions #################################################################
kNN_tau_window <- function(folder){

	KKK <- list.files(folder)
        KKK <- as.character(sort(as.numeric(sapply(strsplit(KKK,"k"),function(x){x[3]}))))


	MAT <- NULL
	TAU <- numeric(length(KKK)-2)
	names(TAU) <- KKK[2:(length(KKK)-1)]
	count <- 1
	for(xx in 2:(length(KKK)-1)){

	    m1 <- readElki(paste(folder,"/Elki_k",KKK[xx-1],sep=""))
	    m2 <- readElki(paste(folder,"/Elki_k",KKK[xx],sep=""))
	    m3 <- readElki(paste(folder,"/Elki_k",KKK[xx+1],sep=""))
	   
        cat(length(m1),length(m2),length(m3), "-",KKK[xx],"\n")

	MAT <- cbind(m1,m2,m3)
	#calculate pairwise Tau
	pairs <- combn(3,2)
	res   <- apply(pairs,2,function(x){cor(MAT[,x[1]],MAT[,x[2]], method="kendall")})
        #res   <- apply(pairs,2,function(x){cor(MAT[,x[1]],MAT[,x[2]], method="spearman")})
	TAU[count] <- median(res)
	MAT <- NULL
	count <- count + 1
	}
	return(TAU)
}

readElki <- function(folder){

tt <- as.character(read.table(paste(folder,"/cluster.txt",sep=""))[,5]) #5#3
tt <- strsplit(tt,'=')
tt <- sapply(tt,function(x){return(x[2])})
tt <- as.numeric(tt)
return(tt)

}

