kNNCallElkiALL <- function(ElkiIN, methods, samples.k=100, n.regions=1000){

if(is.element("ALL", methods)){
methods <- c("LOF","KNN","LOOP","FASTABOD","LDOF","INFLO","COF","ODIN","KNNW","SIMPLIFIEDLOF","LDF")
}

reg    <- 1:n.regions
chunks <- split(reg, ceiling(seq_along(reg)/(length(reg)/samples.k)))
k      <- sapply(chunks, median)
k      <- floor(k)
names(k) <- NULL

#LOF
LOFtime <-system.time(
if(is.element("LOF", methods)){
folder <- paste("LOF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.LOF    -lof.k",kc,"-dbc.in ", ElkiIN," -out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#KNN
KNNtime <- system.time(
if(is.element("KNN", methods)){
folder <- paste("KNN",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.KNNOutlier -knno.k",kc,"-dbc.in ", ElkiIN," -out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#LOOP
LOOPtime <- system.time(
if(is.element("LOOP", methods)){
folder <- paste("LOOP",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.LoOP    -loop.kcomp",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#FASTABOD
FASTABODtime <- system.time(
if(is.element("FASTABOD", methods)){
folder <- paste("FASTABOD",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.FastABOD    -fastabod.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#LDOF
LDOFtime <- system.time(
if(is.element("LDOF", methods)){
folder <- paste("LDOF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.LDOF    -ldof.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#INFLO
INFLOtime <- system.time(
if(is.element("INFLO", methods)){
folder <- paste("INFLO",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.INFLO    -inflo.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#COF
COFtime <- system.time(
if(is.element("COF", methods)){
folder <- paste("COF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.COF    -cof.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#ODIN
ODINtime <- system.time(
if(is.element("ODIN", methods)){
folder <- paste("ODIN",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.ODIN    -odin.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#KNNW
KNNWtime <- system.time(
if(is.element("KNNW", methods)){
folder <- paste("KNNW",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.KNNWeightOutlier   -knnwod.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#SIMPLIFIEDLOF
SIMPLIFIEDLOFtime <- system.time(
if(is.element("SIMPLIFIEDLOF", methods)){
folder <- paste("SIMPLIFIEDLOF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.SimplifiedLOF   -lof.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#LDF
LDFtime <- system.time(
if(is.element("LDF", methods)){
folder <- paste("LDF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.LDF   -ldf.k",kc,"-ldf.h 1 -dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
TIME <- c(LOFtime, KNNtime, LOOPtime, FASTABODtime, LDOFtime, INFLOtime, COFtime, ODINtime, KNNWtime, SIMPLIFIEDLOFtime, LDFtime)
names(TIME) <- c("LOF","KNN","LOOP","FASTABOD","LDOF","INFLO","COF","ODIN","KNNW","SIMPLLOF","LDF")
print(TIME)
return(TIME)
}# End of function 

kNNCallElkiSINGLE <- function(ElkiIN, methods, k=NA){

if(is.element("ALL", methods)){
methods <- c("LOF","KNN","LOOP","FASTABOD","LDOF","INFLO","COF","ODIN","KNNW","SIMPLIFIEDLOF","LDF")
}

#LOF
LOFtime <-system.time(
if(is.element("LOF", methods)){
folder <- paste("LOF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.LOF    -lof.k",kc,"-dbc.in ", ElkiIN," -out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#KNN
KNNtime <- system.time(
if(is.element("KNN", methods)){
folder <- paste("KNN",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.KNNOutlier -knno.k",kc,"-dbc.in ", ElkiIN," -out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#LOOP
LOOPtime <- system.time(
if(is.element("LOOP", methods)){
folder <- paste("LOOP",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.LoOP    -loop.kcomp",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#FASTABOD
FASTABODtime <- system.time(
if(is.element("FASTABOD", methods)){
folder <- paste("FASTABOD",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.FastABOD    -fastabod.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#LDOF
LDOFtime <- system.time(
if(is.element("LDOF", methods)){
folder <- paste("LDOF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.LDOF    -ldof.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#INFLO
INFLOtime <- system.time(
if(is.element("INFLO", methods)){
folder <- paste("INFLO",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.INFLO    -inflo.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#COF
COFtime <- system.time(
if(is.element("COF", methods)){
folder <- paste("COF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.COF    -cof.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#ODIN
ODINtime <- system.time(
if(is.element("ODIN", methods)){
folder <- paste("ODIN",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.ODIN    -odin.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#KNNW
KNNWtime <- system.time(
if(is.element("KNNW", methods)){
folder <- paste("KNNW",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.KNNWeightOutlier   -knnwod.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#SIMPLIFIEDLOF
SIMPLIFIEDLOFtime <- system.time(
if(is.element("SIMPLIFIEDLOF", methods)){
folder <- paste("SIMPLIFIEDLOF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.SimplifiedLOF   -lof.k",kc,"-dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
#LDF
LDFtime <- system.time(
if(is.element("LDF", methods)){
folder <- paste("LDF",ElkiIN, sep="")
for (xx in 1:length(k)){

kc  <- k[xx]
string <- paste("java -cp elki.jar de.lmu.ifi.dbs.elki.application.KDDCLIApplication -algorithm  outlier.lof.LDF   -ldf.k",kc,"-ldf.h 1 -dbc.in ", ElkiIN,"-out", paste(folder,"/Elki_k",kc, sep=""))
cat(string,"\n")
system(string)
}
}
)[3]
TIME <- c(LOFtime, KNNtime, LOOPtime, FASTABODtime, LDOFtime, INFLOtime, COFtime, ODINtime, KNNWtime, SIMPLIFIEDLOFtime, LDFtime)
names(TIME) <- c("LOF","KNN","LOOP","FASTABOD","LDOF","INFLO","COF","ODIN","KNNW","SIMPLLOF","LDF")
print(TIME)
return(TIME)
}# End of function 

