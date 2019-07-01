
kNN_calc_best_k <- function(tau, thres=0.90){

RUNS <- vector("list", length(tau))
r    <- 1

	for(xx in 1:length(tau)){
	   
		if(is.na(tau[xx])){
		r <- r+1
		next
		}

		if(tau[xx]>thres){
		RUNS[[r]] <- c(RUNS[[r]], names(tau[xx]))
		}else{
		r <- r + 1
		}
		
	}
n  <- sapply(RUNS, length)
id <- which.max(n)
bestk <- median(as.numeric(RUNS[[id]]))
return(ceiling(bestk))
}
