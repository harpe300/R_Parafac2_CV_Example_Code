
# Split Half Analysis to validate the CV-MSE chosen model using Tucker congruence coefficients as a metric of similarity 
# If the same (or a similar approximation) solution is replicated across halves, this is validation that the model reflects genuine factors that parsimoniously explain the data across multiple dimensions.

library(dplyr)
library(multiway)
library(abind)
library(lattice)
library(latticeExtra)

setwd('/LOCAL/')

# A 'just for fun' function to generate random seed values using alphanumeric inputs (from https://stackoverflow.com/a/10913336)
convert2seed <- function(x) {
  #require("digest")
  hexval <- paste0("0x", digest::digest(x,"crc32"))
  intval <- type.convert(hexval) %% .Machine$integer.max
  return(intval)
}


# Misc variables and functions
library(eegkit)
biosemi <- c("FP1", "FPZ", "FP2", "AF8", "AF4", "AFZ", "AF3", "AF7", "F7", "F5", "F3", "F1", "FZ", "F2", "F4", "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCZ", "FC1", "FC3", "FC5", "FT7", "T7", "C5", "C3", "C1", "CZ", "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2", "CPZ", "CP1", "CP3", "CP5", "TP7", "P7", "P5", "P3", "P1", "PZ", "P2", "P4", "P6", "P8", "PO8", "PO4", "POZ", "PO3", "PO7", "O1", "OZ", "O2")
data(eegcoord)
cidx <- match(biosemi, rownames(eegcoord))

ncol <- 256
colseq <- colorRampPalette(rev(c("darkred","red","#EF8A62","white","#67A9CF","blue","darkblue")))(ncol)

bin2msec <- function(x, tbin=129, samplerate=64) {
  res <- 1000*(x-tbin)/samplerate
  res
}
logspace <- function(from, to, length.out) {
  # logarithmic spaced sequence
  exp(seq(log(from), log(to), length.out = length.out))
                                           }
xvals <- bin2msec(129:(129+64))
yvals <- round(logspace(2, 30, 40),2)[7:40]



load('./data_for_analysis.RData')
ls()
rm(esmf_go)
gc()

## Make Four Random Splits ##
# This randomly allocates individuals to a fold while keeping the clustered nature of the data intact, thus ensuring full independence across splits. 
# This avoids spuriously inflating the congruence coefficients by calculating congruence between two models from non-independent splits of data.

nfolds <- 4

my_ids <- ids

my_fams <- unique(floor(my_ids/100))

idyr <- rle(floor(my_ids/100));
twinpair <- NULL # vector of length(id) where 1 = unmatched twin and 2 = complete pair twin
for(i in seq(along <- my_ids)) {
  idx = match(floor(my_ids[i]/100),idyr$values);
  twinpair[i] = idyr$lengths[idx];
}

my_ids_completepair = my_ids[twinpair %in% 2] # individual-level id for full cluster data
my_fams_completepair = unique(floor(my_ids[twinpair %in% 2]/100)) # cluster-level id

my_ids_singletons   = my_ids[twinpair %in% 1] # individual-level id for incomplete cluster data

seed <- convert2seed('Pompeii')
seed
set.seed(seed) 

comp_pair_folds = sample(rep(1:nfolds, length.out = length(my_fams_completepair))) # rand allocation of full cluster data to a fold
singleton_folds = sample(rep(1:nfolds, length.out = length(my_ids_singletons))) # rand allocation of incomplete cluster data to a fold

foldids <- rep(NA, length(my_ids))
for ( ff in 1:nfolds ) foldids[which(floor(my_ids/100) %in% my_fams_completepair[which(comp_pair_folds %in% ff)])] <- ff
for ( ff in 1:nfolds ) foldids[which(my_ids %in% my_ids_singletons[which(singleton_folds %in% ff)])] <- ff
foldids

set.seed(NULL)


#  Split Half analysis for Parafac2 model (Harshman, 1994)

# Prepare data for AB / CD splits and run Parafac2 model on AB and CD halves 
Xs_AB <- nscale(esmf_nogo[,7:40,,(foldids %in% c(1,2))], mode=2) # scale freq mode
Xs_CD <- nscale(esmf_nogo[,7:40,,(foldids %in% c(3,4))], mode=2) # scale freq mode

gc()

R <- 5

cl <- makeCluster(10)
ce <- clusterEvalQ(cl, library(multiway))
clusterSetRNGStream(cl, 42)
resAB <- parafac2(Xs_AB, nfac=R, const=c("uncons","unimod","uncons","nonneg"),nstart=10,parallel=TRUE,cl=cl)
stopCluster(cl)
  
cl <- makeCluster(10)
ce <- clusterEvalQ(cl, library(multiway))
clusterSetRNGStream(cl, 42)
resCD <- parafac2(Xs_CD, nfac=R, const=c("uncons","unimod","uncons","nonneg"),nstart=10,parallel=TRUE,cl=cl)
stopCluster(cl)

rm(list = c('Xs_AB','Xs_CD'))
gc()

save.image(file='./ABCD_Nogo_splithalf_111121.RData')

rm(list=ls())
gc()


load('./ABCD_Nogo_splithalf_111121.RData')


# Visually validate similarity in factor structure across split halves

# Reorder or resign if necessary (ordering and sign of factors is arbitrary, i.e., does not affect model fit)
resAB <- reorder(resAB, neworder = c(1,5,2,4,3))
resCD <- reorder(resCD, neworder = c(1,3,4,2,5))
resAB <- resign( resAB, newsign = c(1,1,1,1,1), mode = 'C', absorb = 'A')
resCD <- resign( resCD, newsign = c(1,1,1,1,1), mode = 'C', absorb = 'A')
  
# For split AB, calculate the mean outer product of the time course weights (i.e., AkDk) and frequency weights (Mode B) to visualize modes A and B together
pcs_AB = sapply(paste("pc",1:R,sep="."), function(x) list())
for (jj in 1:R) {
temp <- matrix(0, nrow=dim(resAB$A[[1]]), ncol=dim(resAB$B)[1])

nrecs <- dim(resAB$D)[1]

for ( l in 1:dim(resAB$D)[1] ) {
	 temp <- temp + (tcrossprod(resAB$A[[l]][,jj] * resAB$D[l,jj], resAB$B[,jj])/dim(resAB$D)[1])
}
pcs_AB[[jj]] = levelplot(temp, row.values=xvals, column.values= yvals, scales=list(y=list(at=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30), log=TRUE), x=list(at=c(0, 250, 500, 750, 1000),labels=c(0, 250, 500, 750, 1000)), tck = c(.1,0), alternating=1,draw=FALSE,cex=.0), asp=I(3/5), col.regions= colseq, xlab=list(label="Time (ms)",cex=.0), ylab=list(label="Frequency (Hz)",cex=.00) ,ylim = range(yvals),contour=F, cuts=100, colorkey=FALSE); 
}

# For split CD, calculate the mean outer product of the time course weights (i.e., AkDk) and frequency weights (Mode B) to visualize modes A and B together
pcs_CD = sapply(paste("pc",1:R,sep="."), function(x) list())
for (jj in 1:R) {
temp <- matrix(0, nrow=dim(resCD$A[[1]]), ncol=dim(resCD$B)[1])

nrecs <- dim(resCD$D)[1]
for ( l in 1:dim(resCD$D)[1] ) {
	 temp <- temp + (tcrossprod(resCD$A[[l]][,jj] * resCD$D[l,jj], resCD$B[,jj])/dim(resCD$D)[1])
}
pcs_CD[[jj]] = levelplot(temp, row.values=xvals, column.values= yvals, scales=list(y=list(at=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30), log=TRUE), x=list(at=c(0, 250, 500, 750, 1000),labels=c(0, 250, 500, 750, 1000)), tck = c(.1,0), alternating=1,draw= FALSE,cex=.00), asp=I(3/5), col.regions= colseq, xlab=list(label="Time (ms)",cex=.00), ylab=list(label="Frequency (Hz)",cex=.00) ,ylim = range(yvals),contour=F, cuts=100, colorkey=FALSE); 
} 

gridExtra::grid.arrange(grobs = c(pcs_AB,pcs_CD), ncol = R)

# Visualize Mode C weights
dev.new();par(mfrow=c(2,R))
for (ii in 1:R) { eegspace(eegcoord[cidx, 4:5], resAB$C[,ii], mycolors = colseq, colorbar = FALSE) }
for (ii in 1:R) { eegspace(eegcoord[cidx, 4:5], resCD$C[,ii], mycolors = colseq, colorbar = FALSE) }


# Use TCCs to empirically test split half similarity

# Calc avg time courses (AkDk) for AB (necessary for TCC)
modeA_AB <- matrix(0, nrow=dim(resAB$A[[1]])[1], ncol=R)
for (jj in 1:R) {
  for ( l in 1:length(resAB$A) ) {
     modeA_AB[,jj] <- modeA_AB[,jj] + ((resAB$A[[l]][,jj] * resAB$D[l,jj]) / dim(resAB$D)[1])
  }
}

# Calc avg time courses (AkDk) for CD (necessary for TCC)
modeA_CD <- matrix(0, nrow=dim(resCD$A[[1]])[1], ncol=R)
for (jj in 1:R) {
  for ( l in 1:length(resCD$A) ) {
     modeA_CD[,jj] <- modeA_CD[,jj] + ((resCD$A[[l]][,jj] * resCD$D[l,jj]) / dim(resCD$D)[1])
  }
}

# Compute TCC for each mode
AB_CD_modeA_tcc = congru(modeA_AB,modeA_CD)  
AB_CD_modeB_tcc = congru(resAB$B,resCD$B)
AB_CD_modeC_tcc = congru(resAB$C,resCD$C)

# Multiply the TCCs of each mode as in Smilde, Bro, Geladi (2008, pg 268) to get an overall TCC similarity metric to validate the chosen model
AB_CD_modeA_tcc * AB_CD_modeB_tcc * AB_CD_modeC_tcc

