library(dplyr)
library(multiway)
library(abind)
library(lattice)
library(latticeExtra)

# Initialize data frame
my_list <- c("EnrS_g", "EnrS_n")
ds <- sapply(my_list, function(x) data.frame())

fpath <- "/LOCAL/"

nelec <- 61

# Read first five lines, get class property for all variables
fname <- file.path(fpath, paste0("DATA.dat"))
first5 <- read.table(fname, header=TRUE, sep="\t", nrows=5, comment.char="", stringsAsFactors=FALSE)
my_classes <- sapply(first5, class)

cmd <- paste("wc -l", fname)
nrecs <- as.numeric(sub("(\\d+).+", "\\1", system(cmd, intern=TRUE)))

gng <- read.table(fname, sep="\t", header=TRUE, na="-999", colClasses=my_classes, comment.char="", nrows=nrecs, stringsAsFactors=FALSE)
for ( cond in c("G","N") ) {
  dset <- grep(paste0("EnrS_", tolower(cond)), my_list, value=TRUE)
  temp.dat <- subset(gng, catcodes==cond)
  complete_data <- rle(temp.dat$id)$values[rle(temp.dat$id)$lengths==nelec]
  temp.dat <- subset(temp.dat, id %in% complete_data)
  ds[[dset]] <- temp.dat
}
rm(list=c('gng','temp.dat', 'complete_data', 'dset'))
gc()

fq_bins <- as.numeric(sub(".+t\\d+f(\\d+).dat", "\\1", fname))
srate <- 256
npts_per_sec <- as.numeric(sub(".+t(\\d+)f\\d+.dat", "\\1", fname)) 
time_pts <- 4*npts_per_sec
tbin <- unique(ds[['EnrS_g']]$tbin)


# misc functions
vec2mat <- function(x, nr, nc) {
  my_mat <- matrix(x, nrow=nr, ncol=nc, byrow=TRUE)
  my_mat
}
mat2vec <- function(x) {
  my_vec <- matrix(x, nrow=1, ncol=prod(dim(x)))
  my_vec
}
bin2msec <- function(x, tbin=129, samplerate=64) {
  res <- 1000*(x-tbin)/samplerate
  res
}
msec2bin <- function(x, tbin=129, samplerate=64) {
  sampling_interval <- 1000/samplerate
  res <- round(x/sampling_interval + tbin)
  res
}
dbt <- function(x1, x2) {
  res <- 10*log10(x1/x2)
  res
}
logspace <- function(from, to, length.out) {
  exp(seq(log(from), log(to), length.out = length.out))
                                           }


ids <- unique(ds[['EnrS_g']]['id'])[,1] # ID vector
elecs <- trimws(unique(ds[['EnrS_g']]['elecname'])[,1]) # elec vector

# reduce martices to just the data to reduce memory load
ds[['EnrS_g']] <- subset(ds[['EnrS_g']], id %in% ids, select=grep("^F\\d+", names(ds[['EnrS_g']])))
ds[['EnrS_n']] <- subset(ds[['EnrS_n']], id %in% ids, select=grep("^F\\d+", names(ds[['EnrS_n']])))

# transpose to facilitate reshaping below
dat_g <- t(ds[['EnrS_g']])
dat_n <- t(ds[['EnrS_n']])

rm(ds)
gc()

nrecs <- ncol(dat_n)
nsubs <- length(ids)

# structure into three dimensional array and transpose array by it permuting to the proper ordering
mats_g <- array(dat_g, dim=c(time_pts, fq_bins, nsubs*nelec))
mats_n <- array(dat_n, dim=c(time_pts, fq_bins, nsubs*nelec))
mats_g <- aperm(mats_g, c(3, 2, 1))
mats_n <- aperm(mats_n, c(3, 2, 1))

# Preprocess time-series data
mats_n_baseline = apply(mats_n[,,msec2bin(-450):msec2bin(-250)], 1:2, mean)
mats_n_db = mats_n

for (jj in 1:dim(mats_n_db)[2]) {
  mats_n_db[,jj,] <- (mats_n_db[,jj,]-mats_n_baseline[,jj])/mats_n_baseline[,jj]
}

mats_g_baseline = apply(mats_g[,,msec2bin(-450):msec2bin(-250)], 1:2, mean)
mats_g_db = mats_g

for (jj in 1:dim(mats_g_db)[2]) {
  mats_g_db[,jj,] <- (mats_g_db[,jj,]-mats_g_baseline[,jj])/mats_g_baseline[,jj]
}

# check dimensions
dim(mats_n_db)
ncond <- 1 
nelec <- 61
ncond*nelec*nsubs == dim(mats_n_db)[1] # verify that data is the expected size

# arrange from 3-way array (records, fq, time) to 4-way array (elec, subs, freq, time)
go2 <- array(mats_g_db, dim=c(nelec, ncond*nsubs, dim(mats_g_db)[2], dim(mats_g_db)[3]))
nogo2 <- array(mats_n_db, dim=c(nelec, ncond*nsubs, dim(mats_n_db)[2], dim(mats_n_db)[3]))

# verify that the data was restructured properly (should be identical == TRUE)
identical(go2[22,1,,], mats_g_db[22,,])
identical(nogo2[22,1,,], mats_n_db[22,,])

# rearrange 4-way array (elec, subs, freq, time) to 4-way array (time by freq by elecs by subs)
esmf_go <- aperm(go2, c(4, 3, 1, 2))
esmf_go <- esmf_go[msec2bin(0):(tbin+npts_per_sec),,,] # reduce to desired time range

esmf_nogo <- aperm(nogo2, c(4, 3, 1, 2))
esmf_nogo <- esmf_nogo[msec2bin(0):(tbin+npts_per_sec),,,] # reduce to desired time range

rm(list=c('ds',"dat_g","dat_n","first5","mats_g","mats_g_baseline","mats_n","mats_n_baseline","mats_n_db","mats_g_db","my_classes","my_list"))
gc()

timems <- bin2msec(msec2bin(0):(tbin+npts_per_sec))
freqhz <- round(logspace(2, 30, 40),2)

save(ids, elecs, esmf_go, esmf_nogo, timems, freqhz, file = '/LOCAL/data_for_analysis.RData')
