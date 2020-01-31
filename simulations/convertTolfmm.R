convert_to_lfmm2 <- function(genome.class, n=1000){

lfmm <- NULL
for(xx in 1:n){
	snp  <- genome.class@region.data@biallelic.matrix[[xx]][,]
        if(length(snp)==0){
        #snp <- numeric(300)
        #lfmm <- cbind(lfmm, snp)
        next
        }
	lfmm <- cbind(lfmm, snp)
        close(genome.class@region.data@biallelic.matrix[[xx]])
}

write.table(lfmm, file="LFMM-sel", col.names=FALSE, row.names=FALSE, sep="\t")

#return(lfmm)
}

convert_to_lfmm <- function(genome.class, n=1000){

lfmm <- NULL
for(xx in 1:n){
	snp  <- genome.class@region.data@biallelic.matrix[[xx]][,]
	lfmm <- cbind(lfmm, snp)
        close(genome.class@region.data@biallelic.matrix[[xx]])
}

write.table(t(lfmm), file="LFMM-sel", col.names=FALSE, row.names=FALSE, sep="\t")

#return(lfmm)
}

convert_to_ped <- function(genome.class, n=1000){

ped <- NULL
count <- 0
for(xx in 1:n){
	snp  <- genome.class@region.data@biallelic.matrix[[xx]][,]
        if(length(snp)==0){
        next
        }
        snp2 <- snp
        snp2[snp==1] <- "T T"
        snp2[snp==0] <- "A A"
	ped  <- cbind(ped, snp2)
        close(genome.class@region.data@biallelic.matrix[[xx]])
count <- count + dim(snp)[2] 
}

ped <- as.data.frame(ped)


# create population info 
pop_info <- data@region.data@populations[[1]]
m <- sum(sapply(pop_info,length))

pop_vec  <- numeric(m)
pop_vec[pop_info[[1]]] <- 1 
pop_vec[pop_info[[2]]] <- 2
pop_vec[pop_info[[3]]] <- 3
 
#pop_vec <- paste("P", pop_vec, sep="")

ped <- cbind(pop_vec, 1:m, 0, 0, 0, 0, ped)
ped <- as.data.frame(ped)

# store the PED file 
write.table(ped, file="FLK.ped", col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)

# Create the Map file 
write.table(cbind("chr",paste("chr",1:count,sep=""),0,1:count), file="FLK.map", col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)


#return(lfmm)
}

#####################################################################
convert_to_dxy <- function(genome.class,name){

dxy <- NULL
for(xx in 1:n){
	dxy  <- genome.class@nuc.diversity.between
	#dxy <- cbind(dxy, snp)
}

write.table(t(dxy), file=name, col.names=FALSE, row.names=FALSE, sep="\t")

#return(lfmm)
}
convert_to_fst <- function(genome.class, name){

fst <- NULL
for(xx in 1:n){
	fst <- genome.class@nuc.F_ST.pairwise
	#fst <- c(fst, snp)
}

write.table(t(fst), file=name, col.names=FALSE, row.names=FALSE, sep="\t")

#return(lfmm)
}
