
# Use levelplot and false color 2D plots to visualize the raw time-frequency transformed signal data

library(dplyr)
library(multiway)
library(abind)
library(lattice)
library(latticeExtra)

setwd('/LOCAL/data/')


load('./data_for_analysis.RData')
ls()
esmf_nogo <- esmf_nogo[,7:40,,]

###  Plotting TF surfaces  ###

ncol <- 256
colseq <- colorRampPalette(rev(c("darkred","red","#EF8A62","white","#67A9CF","blue","darkblue")))(ncol)

zlim_MF <- max(abs(range(apply(esmf_nogo[,,which(elecs %in% c('CZ','FCZ')),], 1:2, mean))))

nogo_MF <- levelplot(apply(esmf_nogo[,,which(elecs %in% c('CZ','FCZ')),], 1:2, mean), 
  row.values=timems, column.values=freqhz[7:40], 
  scales=list(y=list(at=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30), log=TRUE),x=list(at=c(0, 250, 500, 750, 1000),labels=c(0, 250, 500, 750, 1000)),tck = c(.1,0),alternating=1,draw=TRUE,cex=1), 
  asp=I(3/5), 
  col.regions=colseq, 
  at=c(-Inf, seq(-1* zlim_MF, zlim_MF, length=100),Inf), 
  xlab = "Time (ms)", ylab="Frequency (Hz)",ylim = range(freqhz[7:40]),
  contour=F, cuts=100, colorkey=list(width = 0.8), 
  main=list(label=c('Nogo - Medial Frontal'),cex=1,vjust=.8)); 

# Alternatively, plot surfaces and contours

nogo_MF_contour <- levelplot(apply(esmf_nogo[,,which(elecs %in% c('CZ','FCZ')),], 1:2, mean), 
  row.values=timems, column.values=freqhz[7:40], 
  scales=list(y=list(at=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30), log=TRUE),
x=list(at=c(0, 250, 500, 750, 1000),labels=c(0, 250, 500, 750, 1000)),
tck = c(.1,0),alternating=1,draw=TRUE,cex=1), 
  asp=I(3/5), 
  col.regions=colseq, 
  at=c(-Inf, seq(-1* zlim_MF, zlim_MF, length=100),Inf), 
  xlab = "Time (ms)", ylab="Frequency (Hz)",ylim = range(freqhz[7:40]),
  contour=F, cuts=100, colorkey=list(width = 0.8), 
  main=list(label=c('Nogo - Medial Frontal'),cex=1,vjust=.8)) + 
contourplot(apply(esmf_nogo[,,which(elecs %in% c('CZ','FCZ')),], 1:2, mean), 
  row.values=timems, column.values=freqhz[7:40], 
  scales=list(y=list(at=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30), log=TRUE),x=list(at=c(0, 250, 500, 750, 1000),labels=c(0, 250, 500, 750, 1000)),tck = c(.1,0),alternating=1,draw=TRUE,cex=1), 
  asp=I(3/5), 
  col.regions=colseq, 
  at=c(-Inf, seq(-1* zlim_MF, zlim_MF, length=10),Inf), 
  xlab = "Time (ms)", ylab="Frequency (Hz)",
  ylim = range(freqhz[7:40]),
  contour=T, 
  colorkey=list(width = 0.8), 
  main=list(label=c('Nogo - Medial Frontal'),cex=1,vjust=.8),labels = FALSE)

# Print a blank plot to use a legend in Affinity Designer figure creation
blank_plot <- levelplot(matrix(0,65,34), 
  row.values=timems, column.values=freqhz[7:40], 
  scales=list(y=list(at=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30), log=TRUE),tck = c(.5,0),alternating=1,draw=TRUE,x=list(at=seq(0,1000,250))), 
  asp=I(3/5), 
  col.regions=colseq, 
  at=seq(-1000, 1000, length=100), 
  xlab = "Time (ms)", ylab="Frequency (Hz)",ylim = range(freqhz[7:40]),
  contour=T, cuts=100, 
  colorkey=list(title='Relative Δ from baseline', at=seq(-1000,1000, length=100), col.regions=colseq, width = 0.5, labels=list(at=c(-1000, 0, 1000), labels=c("-", "0.0", "+"))),
  main=list(label='TF Legend',cex=.75))

gridExtra::grid.arrange(grobs = list(nogo_MF, nogo_MF_contour, blank_plot), ncol = 1)


# Alternate way to visualize the data according to spatial location

library(ggplot2)

cols <-pals::parula(9)
ncol <- 100
colseq <- colorRampPalette(cols)(ncol)

# Initialize list to store the ggplot objects
pcs = sapply(paste("elec",1:length(elecs),sep="."), function(x) list())

for (jj in 1:length(elecs)) {

zlim <- 1.5

dat <- expand.grid(time = timems, freq = freqhz[7:40])
z <- expand.grid(apply(esmf_nogo[,,jj,], 1:2, mean))

dat$pow <- z[,1]

dat <- tibble(dat)

ff <- ggplot(dat) + 
  aes(x = time, y = freq, z = pow, fill = pow) + 
  ggtitle(elecs[jj]) +
  geom_tile(aes(fill = pow)) +
  scale_fill_gradientn( colours = colseq, limits = c(zlim*-1,zlim)) + 
  theme_void() + 
  scale_y_log10(breaks=c(3.03 ,5 ,8 ,12 ,17 ,23 ,30),labels=c(3 ,5 ,8 ,12 ,17 ,23 ,30)) +
  theme(aspect.ratio = 3/5, plot.title = element_text(hjust = 0.5,vjust=-2,size=10), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  guides(fill= 'none')

pcs[[jj]] = ff

}

lay = rbind(c(NA,NA, 1,NA, 2,NA, 3,NA,NA),
            c(NA, 8, 7,NA, 6,NA, 5, 4,NA),
            c( 9,10,11,12,13,14,15,16,17),
            c(26,25,24,23,22,21,20,19,18),
            c(27,28,29,30,31,32,33,34,35),
            c(44,43,42,41,40,39,38,37,36),
            c(45,46,47,48,49,50,51,52,53),
            c(NA,58,57,NA,56,NA,55,54,NA),
            c(NA,NA,59,NA,60,NA,61,NA,NA))

gridExtra::grid.arrange(grobs = pcs, layout_matrix = lay)

