library(dplyr)
library(multiway)
library(abind)
library(lattice)
library(latticeExtra)

setwd('/LOCAL/')

load('./data_for_analysis.RData')
ls()

# Slab-scale within levels of second mode (reduces the influence of the power law on magnitude)
Xs <- nscale(esmf_nogo[,7:40,,], mode=2)

sqrt(rowMeans(aperm(Xs, c(2,1,3,4))^2)) # verify that sum of squares are 1, as expected

rm(esmf_nogo,esmf_go)
gc()

# Initialize empty data frame for output
pf2 <- sapply(paste("pf2",1:10,sep="."), function(x) data.frame())

<code># Fit a four-mode extension of Harshman's Parafac2 model (Helwig, Hong, & Bokhari, 2012) to the 4-way array to parsimoniously reduce dimensionality by explaining the maximum amount of variance across all four dimensions simultaneously. Contrast this with PCA/ICA/etc. reduction methods that only work on 2-way data and require the researcher to either 'throw out' portions of the data or unfold a multiway array into a two-dimensional matrix, which causes significant interpretational difficulties. Unlike PCA/ICA/etc., Parafac2 allows for variation across a mode, meaning that the model can flexibly capture and represent important between-variable differences in the component loadings (e.g., regional variation in sales, between-participant latency shifts in a time-series vector, etc.) </code>
 
# Use 10 random starts for each model, choosing the one with the lowest SSE as best-fitting. Optimizes the chances that the chosen solution has within-sample uniqueness and the optimality of an obtained solution, rather than a solution that was caught in highly suboptimal local minima during optimization.

# Fit models with 1 to 10 factors 
for ( R in 1:10 ) {
  print(paste("N of components:", R), quote=FALSE)
  
  # parallelize across 6 cores to speed up processing
  cl <- makeCluster(6)
  ce <- clusterEvalQ(cl, library(multiway))
  clusterSetRNGStream(cl, 42)
  res <- parafac2(Xs, nfac=R, const=c("uncons","unimod","uncons","nonneg"),nstart=10,parallel=TRUE,cl=cl)
  stopCluster(cl)
  
  pf2[[R]] <- res
}

save(pf2, file = './DATA_parafac2_1to10c_10start.RData')

rm(Xs)
gc()


# Evaluate change in sum of squared error as a function of factor size to estimate the most appropriate model size. Plots of the SSE approximate a Scree plot
pf2_nogo <- pf2
rm(pf2)

sse <- sapply(pf2_nogo, function(x) x$SSE)
r2 <- sapply(pf2_nogo, function(x) x$Rsq)
normerr <- 1-r2  # Equivalent to sse/sumsq

dev.new();
plot(1:length(sse),sse,main="DATA",xlab = 'Factors',ylab = 'SSE')
dev.new();
plot(2:length(sse),rev(diff(rev(sse))),type='p',main="DATA",xlab = 'Factors',ylab = 'Δ SSE')
dev.new();
barplot(rev(diff(rev(sse))),names.arg = c(2,3,4,5,6,7,8,9,10),main="DATA",xlab = 'Factors',ylab = 'Δ SSE')
dev.new();
barplot(rev(diff(rev(normerr))),names.arg = c(2,3,4,5,6,7,8,9,10),main="DATA",xlab = 'Factors',ylab = 'Δ error')
txtplot::txtplot(2:10,rev(diff(rev(sse))),,xlab = 'Factors',ylab = 'Δ SSE')




