library(dplyr)
library(multiway)
library(abind)
library(lattice)
library(latticeExtra)


setwd('/Users/harpe300/Google Drive/statistics/e11f23x_c11f34x/multiway/data/')


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
  # blatantly stolen from library("emdbook"), because need only this
  exp(seq(log(from), log(to), length.out = length.out))
                                           }
timems <- bin2msec(129:(129+64))
freqhz <- round(logspace(2, 30, 40),2)


load('./data_for_analysis.RData')
ls()

rm(list = c('esmf_go', 'esmf_nogo'))
gc()

dataset <- NULL

# Load Parafac2 model and visualize using levelplot
load('./Nogo_parafac2_5c_50start.RData')

res_No_5fac_50starts <- resign(res_No_5fac_50starts, newsign = c(1,1,1,1,1), mode = 'C', absorb = 'A')

res <- res_No_5fac_50starts

R <- dim(res$D)[2]

nrecs <- dim(res$D)[1]
pcs_nogo = sapply(paste("pc",1:R,sep="."), function(x) list())

for (jj in 1:R) {

temp <- matrix(0, nrow=dim(res$A[[1]])[1], ncol=dim(res$B)[1])

# Calculate the mean outer product of the time course weights (i.e., AkDk) and frequency weights (Mode B) to visualize modes A and B together
for ( l in 1:nrecs ) {
   temp <- temp + (tcrossprod(res$A[[l]][,jj] * res$D[l,jj] , res$B[,jj])/dim(res$D)[1])
}

ff <- levelplot(temp, row.values=timems, column.values=freqhz[7:40], scales=list(y=list(at=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30), log=TRUE), x=list(at=c(0, 250, 500, 750, 1000),labels=c(0, 250, 500, 750, 1000)), tck = c(.1,0), alternating=1,draw=TRUE,cex=.75), asp=I(3/5), col.regions= colseq, xlab=list(label="Time (ms)",cex=.75), ylab=list(label="Frequency (Hz)",cex=.75) ,ylim = range(freqhz[7:40]),contour=F, cuts=100, colorkey=list( col.regions= colseq, width = 0.5), main=list(label=paste("N Factor ", jj, ' %VarExp: ', 100*round(res$Rsq * colSums(res$D^2) / sum(res$D^2),3)[jj]))); 

pcs_nogo[[jj]] = ff

}


# Load Parafac2 model and visualize using levelplot
load('./Go_parafac2_4c_50start.RData')

res_Go_4fac_50starts <- resign(res_Go_4fac_50starts, newsign = c(1,1,1,1), mode = 'C', absorb = 'A')

res <- res_Go_4fac_50starts

R <- dim(res$D)[2]

nrecs <- dim(res$D)[1]
pcs_go = sapply(paste("pc",1:R,sep="."), function(x) list())

for (jj in 1:R) {

temp <- matrix(0, nrow=dim(res$A[[1]])[1], ncol=dim(res$B)[1])

# Calculate the mean outer product of the time course weights (i.e., AkDk) and frequency weights (Mode B) to visualize modes A and B together
for ( l in 1:nrecs ) {
   temp <- temp + (tcrossprod(res$A[[l]][,jj] * res$D[l,jj], res$B[,jj])/dim(res$D)[1])
}

ff <- levelplot(temp, row.values=timems, column.values=freqhz[7:40], scales=list(y=list(at=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30), log=TRUE), x=list(at=c(0, 250, 500, 750, 1000),labels=c(0, 250, 500, 750, 1000)), tck = c(.1,0), alternating=1,draw=TRUE,cex=.75), asp=I(3/5), col.regions= colseq, xlab=list(label="Time (ms)",cex=.75), ylab=list(label="Frequency (Hz)",cex=.75) ,ylim = range(freqhz[7:40]),contour=F, cuts=100, colorkey=list( col.regions= colseq, width = 0.5), main=list(label=paste("G Factor ", jj, ' %VarExp: ', 100*round(res$Rsq * colSums(res$D^2) / sum(res$D^2),3)[jj]))); 

pcs_go[[jj]] = ff

}


lay = rbind(c(1, 6),
            c(2, 7),
            c(3, 8),
            c(4, 9),
            c(5,NA))
            
            
gridExtra::grid.arrange(grobs = c(pcs_nogo,pcs_go), layout_matrix = lay)


# Vizualize Mode C factor weights according to spatial location
dev.new();par(mfcol=c(5,2))

for (ii in 1:5) { eegspace(eegcoord[cidx, 4:5], res_No_5fac_50starts$C[,ii], mycolors = colseq, ncolor = 100, colorlab = paste("Component: ", ii, ' Mode C weights')) }
for (ii in 1:4) { eegspace(eegcoord[cidx, 4:5], res_Go_4fac_50starts$C[,ii], mycolors = colseq, ncolor = 100, colorlab = paste("Component: ", ii, ' Mode C weights')) }


