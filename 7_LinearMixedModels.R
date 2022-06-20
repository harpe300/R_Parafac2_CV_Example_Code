library(dplyr)
library(multiway)
library(abind)
library(lattice)
library(latticeExtra)
library(car)
library(lme4)
library(lmerTest)

#source("/LOCAL/lm_beta_lmer.R");
# function to calculate standardized beta from lme4 lmer objects (which only return unstandardized beta)
lm.beta.lmer <- function(mod) {
   b <- fixef(mod)[-1]
   sd.x <- apply(cbind(getME(mod,"X")[,-1]),2,sd) # added cbind to account for 1-IV models 02/27/18
  #sd.x <- apply(getME(mod,"X")[,-1],2,sd) # 
   sd.y <- sd(getME(mod,"y"))
   b*sd.x/sd.y
}

setwd('/LOCAL/')

library(eegkit)
biosemi <- c("FP1", "FPZ", "FP2", "AF8", "AF4", "AFZ", "AF3", "AF7", "F7", "F5", "F3", "F1", "FZ", "F2", "F4", "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCZ", "FC1", "FC3", "FC5", "FT7", "T7", "C5", "C3", "C1", "CZ", "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2", "CPZ", "CP1", "CP3", "CP5", "TP7", "P7", "P5", "P3", "P1", "PZ", "P2", "P4", "P6", "P8", "PO8", "PO4", "POZ", "PO3", "PO7", "O1", "OZ", "O2")
data(eegcoord)
cidx <- match(biosemi, rownames(eegcoord))

# cmap
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
yvals <- round(logspace(2, 30, 40),2)


load('./data_for_analysis.RData')
ls()
rm(list = c('esmf_go', 'esmf_nogo'))
gc()

dataset <- NULL

#Nogo
load('./Nogo_parafac2_5c_50start.RData')

res_No_5fac_50starts <- resign(res_No_5fac_50starts, newsign = c(1,1,1,1,1), mode = 'C', absorb = 'A')

dataset$ID             <- ids
dataset$N_alpha        <- res_No_5fac_50starts$D[,1]
dataset$N_mftheta      <- res_No_5fac_50starts$D[,2]
dataset$N_highalphbeta <- res_No_5fac_50starts$D[,3]
dataset$N_earlytheta   <- res_No_5fac_50starts$D[,4]
dataset$N_frontalbeta  <- res_No_5fac_50starts$D[,5]

# Go
load('./Go_parafac2_4c_50start.RData')

res_Go_4comp <- resign(res_Go_4fac_50starts, newsign = c(1,1,1,1), mode = 'C', absorb = 'A')

dataset$G_alpha        <- res_Go_4fac_50starts$D[,1]
dataset$G_earlytheta   <- res_Go_4fac_50starts$D[,2]
dataset$G_highalphbeta <- res_Go_4fac_50starts$D[,3]
dataset$G_smabeta      <- res_Go_4fac_50starts$D[,4]

dataset <- as.data.frame(dataset)

gc()

# merge with additional datasets and rename some variables
agesx   <- read.table("/LOCAL/agesx_10132019.txt", header=TRUE, sep=",", na="NA");
agesx   <- subset(agesx, select = c("ID","fam_clust","IDSC","IDES","AGE_FU4","AGE_FU6","IDSEX","ZYGOSITY"));
colnames(agesx) <- c("ID","fam_clust","IDSC","IDES","AGE_FU4",'AGE_FU6',"sex","zyg")
dataset <- merge(dataset,agesx,by="ID",all.x=TRUE);
dataset$AGE <- apply(X=cbind(dataset$AGE_FU4,dataset$AGE_FU6), MARGIN=1, FUN=max,na.rm=T)

# Make factors/etc.
dataset$sex[dataset$sex %in% 1] <- -1
dataset$sex[dataset$sex %in% 2] <-  1

dataset$zyg[dataset$zyg %in% 1] <-  1
dataset$zyg[dataset$zyg %in% 2] <- -1

dataset$cohort = NULL
dataset$cohort[dataset$IDES == 9] = 0
dataset$cohort[dataset$IDES != 9] = 1

dataset$sex <- as.factor(dataset$sex)
dataset$cohort <- as.factor(dataset$cohort)

dataset$AGE_center <- scale(dataset$AGE, center=TRUE, scale=FALSE)

dataset$sex_wec <- as.factor(dataset$sex)
contrasts(dataset$sex_wec) <- wec::contr.wec(dataset$sex_wec, omitted="-1")


# Add behavioral data
ESbeh.data <- read.table("./es_nogo_behdata_acc_rt.txt", header=TRUE, sep="\t", na=c(-999, "NaN"));
MFbeh.data <- read.table("./mf_nogo_behdata_acc_rt.txt", header=TRUE, sep="\t", na=c(-999, "NaN"));
beh.data <- rbind(ESbeh.data,MFbeh.data)
beh.data$ID <- as.numeric(paste(substr(beh.data$subname,1,5),substr(beh.data$subname,7,8),sep=''))

dataset <- merge(dataset, beh.data,by="ID",all.x=TRUE);
rm(list = c('ESbeh.data', 'MFbeh.data', 'beh.data'))

# PID-5
PID_ESMF.data <- foreign::read.spss(file = './PID5_ESMF.sav', to.data.frame = TRUE)
myvars <- c("ID", "Anh", "Anh_MISS", "Anx", "Anx_MISS", "Att", "Att_MISS", "Cal", "Cal_MISS", "Dec", "Dec_MISS", "Dep", "Dep_MISS", "Dis", "Dis_MISS", "Ecc", "Ecc_MISS", "Emo", "Emo_MISS", "Hos", "Hos_MISS", "Imp", "Imp_MISS", "Int", "Int_MISS", "Irr", "Irr_MISS", "Man", "Man_MISS", "PeD", "PeD_MISS", "Per", "Per_MISS", "Res", "Res_MISS", "Rig", "Rig_MISS", "Ris", "Ris_MISS", "Sep", "Sep_MISS", "Sub", "Sub_MISS", "Sus", "Sus_MISS", "Unu", "Unu_MISS", "Wit", "Wit_MISS", "NegAff", "Detach", "Antago", "Disinh", "Psycho", "INV1", "INV2", "INVALID");
PID_ESMF.data <- PID_ESMF.data[myvars]

dataset <- merge(dataset, PID_ESMF.data,by="ID", all.x = TRUE);
rm(PID_ESMF.data)

# MPQ PBYA
MPQ_ES.data <- read.table("/LOCAL/ESFU4_MPQ_PBYA.dat", header=TRUE, sep="\t", na=c("-99","NA"));
myvars <- c("ID", "AG", "NMISS", "AL", "CON", "HA", "RT", "SR", "WB", "MBPD");
MPQ_ES.data <- MPQ_ES.data[myvars]

MPQ_MF.data <- read.table("/LOCAL/MFFU6_MPQ_PBYA.dat", header=TRUE, sep="\t", na=c("-99","NA"));
MPQ_MF.data <- MPQ_MF.data[myvars]

identical(names(MPQ_MF.data),names(MPQ_ES.data))
MPQ.data <- rbind(MPQ_ES.data, MPQ_MF.data)
dataset <- merge(dataset, MPQ.data, by="ID", all.x = TRUE);

rm(list = c('MPQ_ES.data', 'MPQ_MF.data', 'MPQ.data'))
gc()

# Function to calculate nonparametric residual bootstrap SEs and bias-corrected and accelerated CIs (Carpenter, Goldstein, & Rasbash, 2003; Leeden, Meijer, & Busing, 2007) from lme4 lmer objects
boots <- function(x) {
  set.seed(42)
  b_par <- bootmlm::bootstrap_mer(x,  FUN = fixef, nsim = 5000L, type = "residual_cgr")
  b_sum <- summary(b_par)
  zs <- b_sum$original/b_sum$bootSE
  ps <- 2*pnorm(abs(zs), lower.tail = FALSE)

  lbound <- rep(0,dim(b_sum)[1]); ubound = rep(0,dim(b_sum)[1]);

  Ls <- bootmlm::empinf_merm(x, fixef)

  for(i in 1:dim(b_sum)[1]) {
  	tempbca <- boot::boot.ci(b_par, conf = 0.95, index = i, type = "bca", L = Ls[,i])
    lbound[i] <- tempbca$bca[4]
    ubound[i] <- tempbca$bca[5]
                               }
                               
  lmerSE <- summary(x)$coefficients[,2]  # note that these are not KR SE (which are a little larger)
  print(cbind(round(b_sum,3)[,c(2,4)],round(lmerSE,3),round(zs,3),round(ps,4),round(lbound,3),round(ubound,3)))
  return(b_par)
}


# Log transform these three variables to adjust large right skew
dataset$N_highalphbeta_log <- log1p(dataset$N_highalphbeta)
dataset$N_frontalbeta_log  <- log1p(dataset$N_frontalbeta)
dataset$G_highalphbeta_log <- log1p(dataset$G_highalphbeta)

phenos = names(dataset)[c(2:3,113,5,114,7:8,115,10)]

# Linear mixed model predicting d_prime with all phenos factor scores, adjusting for sex and cohort effects, with a random intercept to account for clustering
formula_LMM = reformulate(paste(c(phenos, ' sex + cohort + (1 | fam_clust)')),response = 'd_prime');
print(formula_LMM);
LMM.model = lmer(formula_LMM, data = dataset)
print(summary(LMM.model, ddf = "Kenward-Roger"));
print(lm.beta.lmer(LMM.model));
print(boots(LMM.model))


# Individual linear mixed model predicting d_prime with one phenotype at a time, adjusting for sex and cohort effects, with a random intercept to account for clustering
for(i in seq_along(phenos)){
writeLines('\n');
writeLines('\n');

writeLines(c('Current variable: ', phenos[i], '\n'));
    formula_LMM = reformulate(paste(c(phenos[i], '+ sex + cohort + (1 | fam_clust)')),response = 'd_prime');
    print(formula_LMM);
    LMM.model = lmer(formula_LMM, data = dataset)
    print(summary(LMM.model, ddf = "Kenward-Roger"));
    print(lm.beta.lmer(LMM.model));
    print(boots(LMM.model))
}



# PID-5 Disinh
dataset$Disinh_sqrt <- sqrt(dataset$Disinh)
dataset_validpid_compdisinh <- dataset[dataset$INVALID %in% 0 & complete.cases(dataset$Disinh),]


# Linear mixed model with both Go-Beta and Nogo-MFTheta predicting Disinh_sqrt, adjusting for sex and cohort effects, and with a random intercept to account for clustering
LMMdis <- lmer(Disinh_sqrt ~ N_mftheta + G_smabeta + sex + cohort + (1 | fam_clust),data = dataset_validpid_compdisinh)
summary(LMMdis, ddf = "Kenward-Roger")
lm.beta.lmer(LMMdis)
boots(LMMdis)

# Visualize the linear mixed model regression results using a partial residual plot
p1 = visreg::visreg(LMMdis,"N_mftheta",overlay=TRUE,line=list(col=c("#000000")),fill=list(col=c("#000000")),points=list(col=c("#000000"), pch=c(1),size=1.5), xlab="Nogo Medial Frontal \n Theta Component", ylab="Disinhibition (sqrt)",gg=TRUE) + ggplot2::theme_classic()
p2  = visreg::visreg(LMMdis,"G_smabeta",overlay=TRUE,line=list(col=c("#000000")),fill=list(col=c("#000000")),points=list(col=c("#000000"), pch=c(1),size=1.5), xlab="Go Sensorimotor \n Beta Component", ylab="Disinhibition (sqrt)",gg=TRUE) + ggplot2::theme_classic()
dev.new();gridExtra:: grid.arrange(p1,p2,nrow=1)



# Estimate reliability of the measures using intraclass correlations

ICC_lme4 <-
function(mod){
	var_B = mod$vcov[1];
	var_W = mod$vcov[2];
	icc_val = round(var_B/(var_B+var_W),digits=3);
	list(icc_val=icc_val, var_B=var_B, var_W=var_W);
					}

fisherz <-
function(r,n,rho0=0,sig=.95){
  z=log((1+r)/(1-r))/2;
  z0=log((1+rho0)/(1-rho0))/2;
  zstar=(z-z0)*sqrt(n-3);
  pval=2*(1-pnorm(abs(zstar)));
  zcrit = cbind(-1,1)*(qnorm((1+sig)/2)/sqrt(n-3));
  zCI = z+zcrit;
  zCI = (exp(2*zCI) - 1)/(exp(2*zCI) + 1);
  list(r=r,n=n,z=z,zstar=zstar,pval=pval,zCI=zCI);
                              }


n_pairs = length(unique(dataset$fam_clust));
n_pairs_MZ = length(unique(dataset$fam_clust[dataset$zyg== 1]));
n_pairs_DZ = length(unique(dataset$fam_clust[dataset$zyg==-1]));
 
for(i in seq_along(phenos)){

writeLines(c('Current variable: ', phenos[i], '\n'));

print(paste("Mean: ", round(mean(dataset[,phenos[i]],na.rm=T),4)))
print(paste("SD: ",   round(sd(dataset[,phenos[i]],na.rm=T),4)))
print(paste("Range: ",round(range(dataset[,phenos[i]],na.rm=T),4)))
print(psych::describe(dataset[,phenos[i]],na.rm=T),digits=3)

writeLines('\n');

# Use the between- and within-cluster variance from the linear mixed model to calculate the intraclass correlations
writeLines('ICC (lme4::lmer, adjusted for sex and cohort effects): ' );

formula_ICC_cov = reformulate(paste(' + sex + cohort + (1 | fam_clust)'),response = phenos[i]);
print(formula_ICC_cov);

writeLines('\n');
LMM_MZ_cov.model = lmer(formula_ICC_cov, data = dataset[dataset$zyg== 1,]);
MZ_cov = ICC_lme4(as.data.frame(VarCorr(LMM_MZ_cov.model)));
mz_cov_icc = fisherz(MZ_cov$icc_val, n_pairs_MZ, 0);

LMM_DZ_cov.model = lmer(formula_ICC_cov, data = dataset[dataset$zyg==-1,]);
DZ_cov = ICC_lme4(as.data.frame(VarCorr(LMM_DZ_cov.model)));
dz_cov_icc = fisherz(DZ_cov$icc_val, n_pairs_DZ, 0);

print(as.data.frame(cbind(mz_cov_icc, dz_cov_icc)));

writeLines('### \n');
}


# Restructure data from long to wide format (one row for each fam_clust) by idsc (individual within a fam_clust)
theData = dataset;
names(theData) <- tolower(names(theData))
theData$zygosity  = theData$zyg
theData$idsex     = theData$sex
theData <- reshape(theData, direction="wide", timevar="idsc", idvar=c("fam_clust", "idsex", "zygosity"), sep="_")
names(theData)<-gsub("(.+)\\.([cishon]_[01])$", "\\1_\\2", names(theData))
