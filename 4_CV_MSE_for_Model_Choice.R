
setwd('/LOCAL/')

#source("../scripts/cv.4way_edit_parallel.R") # code provided at the end of this script

# A 'just for fun' function to generate random seed values using alphanumeric inputs (from https://stackoverflow.com/a/10913336)
convert2seed <- function(x) {
  #require("digest")
  hexval <- paste0("0x", digest::digest(x,"crc32"))
  intval <- type.convert(hexval) %% .Machine$integer.max
  return(intval)
}


load('./data_for_analysis.RData')


# Use a simple 2-fold for this example. Randomly allocate individuals to a fold while keeping the clustered nature of the data intact. ...
# This avoids the spurious attentuation of the mean squared error (MSE) that would likely occur by training and testing the model on ...
# non-independent sets of data. Repeat the 2-fold random allocation procedure three times to get three separate sets of random folds. ...
# Perform the CV procedure once for each of the three sets,  then average the MSE across folds.

# Construct folds for CV
nfolds <- 2
foldid = sapply(paste("foldids_",1:3,sep=""), function(x) data.frame())

my_ids <- ids

my_fams <- unique(floor(my_ids/100)) #

idyr <- rle(floor(my_ids/100));
twinpair <- NULL # vector of length(id) where 1 = unmatched twin and 2 = complete pair twin
for(i in seq(along <- my_ids)) {
  idx = match(floor(my_ids[i]/100),idyr$values);
  twinpair[i] = idyr$lengths[idx];
}

seed <- convert2seed('Miami')
seed
set.seed(seed) 

pair_folds = sample(rep(1:nfolds, length.out = length(my_fams))) # rand allocation of families to each fold
foldids <- rep(NA, length(my_ids))
for ( ff in 1:nfolds ) foldids[which(floor(my_ids/100) %in% my_fams[which(pair_folds %in% ff)])] <- ff
foldid[[1]] = foldids

seed <- convert2seed('Daylight Matters')
seed
set.seed(seed) 

pair_folds = sample(rep(1:nfolds, length.out = length(my_fams))) # rand allocation of families to each fold
foldids <- rep(NA, length(my_ids))
for ( ff in 1:nfolds ) foldids[which(floor(my_ids/100) %in% my_fams[which(pair_folds %in% ff)])] <- ff
foldid[[2]] = foldids

seed <- convert2seed('Home to You')
seed
set.seed(seed)

pair_folds = sample(rep(1:nfolds, length.out = length(my_fams))) # rand allocation of families to each fold
foldids <- rep(NA, length(my_ids))
for ( ff in 1:nfolds ) foldids[which(floor(my_ids/100) %in% my_fams[which(pair_folds %in% ff)])] <- ff
foldid[[3]] = foldids


set.seed(NULL)

# Prepare data (65*34*61*1422 = 191,699,820 entries)
Xs <- nscale(es_nogo[,7:40,,], mode=2)

rm('es_nogo');gc();

# Create empty data frame for CV
es_n_cv.2folds3sets = sapply(paste("cv_",1:3,sep=""), function(x) data.frame())


# Fit Harshman's Parafac2 model to the 4-way array from the training fold, then fit a model to the test fold applying the training model weights, ...
# and calculate the mean-square error (MSE) for the test data relative to the fit. Do this process for factor size 1:10 to use the CV-MSE as a method of choosing the appropriate model
for (aa in 1:length(foldid)) {
  es_n_cv.2folds3sets[[aa]] <- cv.4way_parallel(Xs, 1:10, NA, model="parafac2", nfolds=2, foldid=foldid[[aa]], nstarts = 10, constraints = c("uncons","unimod","uncons","nonneg"), numclust = 10)
  gc();
  set.seed(NULL);
}


# Average the CV-MSE across the three sets
es_n_cv.2folds3sets_avg = NULL
es_n_cv.2folds3sets_avg$nfac = es_n_cv.2folds3sets[[1]][1]
es_n_cv.2folds3sets_avg$mse = rowMeans(data.frame(lapply(es_n_cv.2folds3sets, function(i) unlist(i[2]))))
es_n_cv.2folds3sets_avg$ccd = rowMeans(data.frame(lapply(es_n_cv.2folds3sets, function(i) unlist(i[3]))))


# Evaluate change in the averaged CV-MSE as a function of the number of factors to estimate the most appropriate model dimensionality. Plots of the CV-MSE resemble and can be interpreted as scree plots.
mse <- es_n_cv.2folds3sets_avg$mse

dev.new();
plot(1:length(mse),mse,main="DATA",xlab = 'Factors',ylab = 'mse')
dev.new();
plot(2:length(mse),rev(diff(rev(mse))),type='p',main="DATA",xlab = 'Factors',ylab = 'Δ mse')
dev.new();
barplot(rev(diff(rev(mse))),names.arg = c(2:length(mse)),main="DATA",xlab = 'Factors',ylab = 'Δ mse')
dev.new();
txtplot::txtplot(2:length(mse),rev(diff(rev(mse))),,xlab = 'Factors',ylab = 'Δ mse')

# Alternative algorithm that uses change in the CV-MSE to decide the optimal dimensionality (equation from Helwig and Snodgress 2019 (MSEt−1 − MSEt)/σMSE < δ)
print(paste('Chosen comp:', max(which(rev(diff(rev(mse)))>sd(mse) * .2)+1)))



cv.4way_parallel <-
  function(X, nfac, df = NA,
           model = c("parafac", "parafac2"),
           nfolds = 2, foldid = NULL, 
           verbose = TRUE, nstarts = 10, constraints = NULL, numclust, ...){
    # k-fold cross-validation for 4-way model
    # Nathaniel E. Helwig
    # Edited by JH on Sept 15 2021 to implement parallel computation and more flexibility for model fitting
    
    
    ######***######   INITIAL CHECKS   ######***######
    
    if((identical(grep('package:parall',search(),fixed=TRUE), integer(0))) == TRUE ){
    	require(parallel)
    } 
    
    ### check 'X' input
    xdim <- dim(X)
    lxdim <- length(xdim)
    if(lxdim != 4L) stop("Input 'X' must be 4-way array")
    if(any(is.nan(X)) | any(is.infinite(X))) stop("Input 'X' cannot contain NaN or Inf values") 
    if(any(is.na(X))){
      missingdata <- TRUE
      naid <- which(is.na(X))
    } else {
      missingdata <- FALSE
      xcx <- sumsq(X)
    }
    
    ### check 'nfac' input
    if(missing(nfac)){
      nfac <- 1:8
    } else {
      nfac <- sort(as.integer(nfac))
      if(nfac[1] < 1L) stop("Input 'nfac' must be a vector of positive integers.")
    }
    
    ### check 'model' input
    mode <- 4L
    model <- as.character(model[1])
    if(!any(model == c("parafac", "parafac2"))) stop("Input 'model' must be 'parafac' or 'parafac2'.")
    
    ### check 'nfolds' and 'foldid'
    if(is.null(foldid)){
      nfolds <- as.integer(nfolds[1])
      foldid <- sample(rep(1:nfolds, length.out = xdim[mode]))
    } else {
      if(length(foldid) != xdim[mode]) stop("Input 'foldid' must be of length dim(X)[mode]")
      foldid <- as.integer(foldid)
      nfolds <- length(unique(foldid))
    }
    if(any(tapply(foldid,foldid,length)==1L)){
      warning("Fold assignments resulted in one or more folds with only one observation")
    }
    
     ### check 'numclust' input for parellel processing
    if(missing(numclust)){
      numclust <- detectCores()
    }    
        
    ### check if df was provided
    if(is.na(df[1])){
      
      ### tuning only the number of factors...
      
      ### initialize data frame to hold results
      numfac <- length(nfac)
      result <- data.frame(nfac = nfac,
                           mse = rep(0, numfac),
                           ccd = rep(0, numfac))
      rownames(result) <- nfac
      
      ### k-fold tuning...
      for(kk in 1:nfolds){
        
        ## print progress
        if(verbose) cat("Fitting for fold ", kk, " (out of ",nfolds,")\n", sep = "")
        
        ## get training and testing data
        indx <- which(foldid == kk)
        Xtrain <- X[,,,-indx]
        Xtest <- X[,,,indx]
        nxtest <- prod(dim(X[,,,indx]))
        
        ## loop through nfac factors
        for(rr in 1:numfac){
          
          # print progress
          if(verbose) cat("...number of factors ", rr, " (out of ",numfac,")\n", sep = "")
          
          # parafac or parafac2
          if(model == "parafac"){
            
            cl <- makeCluster(numclust) #cl <- makeCluster(detectCores())
            ce <- clusterEvalQ(cl, library(multiway))
            # fit to training dataset
            mod <- parafac(X = Xtrain, nfac = nfac[rr], 
                           verbose = FALSE, nstart = nstarts,  
                           const = constraints, parallel = TRUE, cl = cl, ...)
            stopCluster(cl)

            cl <- makeCluster(numclust) #cl <- makeCluster(detectCores())
            ce <- clusterEvalQ(cl, library(multiway))
            # predictions for testing data
            fit <- parafac(X = Xtest, nfac = nfac[rr], 
                           Afixed = mod$A, Bfixed = mod$B,
                           Cfixed = mod$C, verbose = FALSE,   
                           nstart = nstarts, const = constraints,   
                           parallel = TRUE, cl = cl, ...)
            stopCluster(cl)
                        
          } else {
            
            cl <- makeCluster(numclust) #cl <- makeCluster(detectCores())
            ce <- clusterEvalQ(cl, library(multiway))
            # fit to training dataset
            mod <- parafac2(X = Xtrain, nfac = nfac[rr], 
                            verbose = FALSE, nstart = nstarts,  
                            const = constraints, parallel = TRUE, cl = cl, ...)
            stopCluster(cl)

            cl <- makeCluster(numclust) #cl <- makeCluster(detectCores())
            ce <- clusterEvalQ(cl, library(multiway))            
            # predictions for testing data
            fit <- parafac2(X = Xtest, nfac = nfac[rr], 
                            Gfixed = smpower(mod$Phi), 
                            Bfixed = mod$B, Cfixed = mod$C,
                            verbose = FALSE, nstart = nstarts,   
                            const = constraints, parallel = TRUE, cl = cl, ...)
            stopCluster(cl)
          }
          
          # get sse.new and ccd.new
          mse.new <- mean((Xtest - fitted(fit))^2)
          ccd.new <- corcondia(Xtest, fit)
          
          # add to previous results
          result$mse[rr] <- result$mse[rr] + mse.new
          result$ccd[rr] <- result$ccd[rr] + ccd.new
          
        } # end for(rr in 1:numfac)
        
      } # end for(kk in 1:nfolds)
      
      ## normalize by number of folds
      result[,2:3] <- result[,2:3] / nfolds
      return(result)
      
    } 
