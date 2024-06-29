################################################################################
##################  1| Fit Vital Rates x Species x Grid    #####################
################################################################################
#--| Housekeeping |-------------------------------------------------------->----
rm(list=ls(all=TRUE)) ## Clear the workspace

library(tidyr)
library(treemendous)
library(tidyverse)
library(reshape)
library(splitstackshape)
library(MASS)
library(beepr)
library(pscl)  # for computation of psuedo R2 values
library(lme4)
library(stringr)
library(fmsb) ## for getting R square of survival model (NagelkerkeR2,m1)
library(MuMIn)
setwd("~/[name of wd]")
#___________________________________________________________________________----
# For model fitting | Skip to line 320 and load quality controlled data.
#--| Data entry |---------------------------------------------------------->----

################################################################################
##################  1| Fit Vital Rates x Species x Grid    #####################
################################################################################
#--| Housekeeping |-------------------------------------------------------->----
rm(list=ls(all=TRUE)) ## Clear the workspace

library(tidyr)
library(treemendous)
library(tidyverse)
library(reshape)
library(splitstackshape)
library(MASS)
library(beepr)
library(pscl)  # for computation of psuedo R2 values
library(lme4)
library(stringr)
library(fmsb) ## for getting R square of survival model (use funciton NagelkerkeR2(m1))
library(MuMIn)
setwd("~/[fill in name of wd]")
#___________________________________________________________________________----
# For model fitting | Skip to line 256 and load quality controlled data.
#--| Data entry |---------------------------------------------------------->----
load('./data/[subset of open access data].Rdata')
grid_info <- read.csv("./data/250k.csv") ## database length - 7809829
names(grid_info )[names(grid_info ) == "tmt_plot_id"] <- "tmt.plot.id" 
dat <- TreeMort_TreeData[,c("database.code","tmt.plot.id","tmt.tree.id", 
                            "species.cor","family.cor","d", "d2", "tree.status2", 
                            "year","year2", "interval","mode.death2")]
stacked_dat <- merge(dat, grid_info, by = 'tmt.plot.id') 
length(unique(grid_info$tmt.plot.id[!grid_info %in% stacked_dat$tmt.plot.id]))  

stacked_dat[is.na(stacked_dat$year2),'database.code']
unique(stacked_dat[is.na(stacked_dat$year2),'database.code'])

length(unique(stacked_dat$hexgrid)) ## 129 unique grids
tapply(stacked_dat$database.code,stacked_dat$hexgrid,length)
length(unique(stacked_dat$species.cor)) ## 6354 species
(length(unique(TreeMort_TreeData$tmt.plot.id)))-(length(unique(stacked_dat$tmt.plot.id))) ## missing 2 plots

# FIAN | Forest Inventory of America, North 
# FIANE | Forest Inventory of America, Northeast 
# FIANW | Forest Inventory of America, Northwest
# FIARM | Forest Inventory of America, Rocky Mountains
# FIAS | Forest Inventory of America, South
# NAL | National Inventory, Alberta
# NBC | National Inventory, British Columbia
# NQU | National Inventory, Quebec
# NSA | National Inventory, Saskatchewan
# FGE | Forest Geo (Central and South America only). FGE.nj = nanjenshan taiwan, FGE.fu = Fushan taiwan, 
# FGE.co - cocoli panama, FGE.sh =sinharaja sri lanka, FGE.wa = wisconsin, FGE.wy = Wytham Woods England
# FPN | Forest Plots Network (Central and south America)
#___________________________________________________________________________----
#  removing databases with third party data restrictions. For sharing with MS 
#--| subset to open access networks. For data sharing with MS |------------>----
stacked_dat <- droplevels(subset(stacked_dat, database.code=="FIAN"
                                 | database.code=="FIANE" | database.code=="FIANW"
                                 | database.code=="FIARM" | database.code=="FIAS"
))

stacked_dat <- stacked_dat[,c("tmt.plot.id","tmt.tree.id", 
                                     "species.cor","family.cor","d", "d2", "tree.status2", 
                                     "year","year2", "interval","mode.death2")] 
stacked_dat <- subset(stacked_dat, interval > 0 & interval < 11 )  ## remove intervals greated than 11 yrs

## sub-sampled the North American to balance the dataset. 
# North <- droplevels(subset(stacked_dat, database.code=="FIAN"
#                                  | database.code=="FIANE" | database.code=="FIANW"
#                                  | database.code=="FIARM" | database.code=="FIAS"
#                                  | database.code=="SYN")); length(North) ## 2,404,588
# North <- North[sample(nrow(North), 1145588), ]
# stack <- droplevels(subset(stacked_dat, !database.code=="FIAN"
#                            & !database.code=="FIANE" & !database.code=="FIANW"
#                            & !database.code=="FIARM" & !database.code=="FIAS"
#                            & !database.code=="SYN")) ## 3,393,098
# stacked_dat <- rbind(North,stack) # balanced dataset. Note: random sampling will influence sp # and traits. 

## remove unnatural mortality due to human dist.
## remove all harvested trees (i.e., unnatural death)
## 1= natural death (unknown if standing or failing, i.e. lying dead) (1s = natural standing, 1f = natural failing, (wind related)
## 1h = natural loss followed by salvage logging - mainly applies for Swiss NFI, 2 = harvested, 3= stem not present, vanished, unknown whether death was natural or tree was harvested
## 4 = no information at all about the tree death, NA= tree is alive
levels(as.factor(stacked_dat$mode.death2))
stacked_dat$mode.death2 <- as.factor(stacked_dat$mode.death2)
stacked_dat$mode.death2  = factor(stacked_dat$mode.death2, 
                                  levels=c(levels(stacked_dat$mode.death2), 
                                           "alive"))
stacked_dat$mode.death2[is.na(stacked_dat$mode.death2)] = "alive" 
stacked_dat <- droplevels(subset(stacked_dat, mode.death2!= "2"))


##correct species names. This will need to be modified for subset in random sample
#--Species Name check----------------------------------------------------------|
levels(stacked_dat$species.cor) <- c(levels(stacked_dat$species.cor), "Populus canadensis")
stacked_dat$species.cor[stacked_dat$species.cor == 'Populus x canadensis'] <- 'Populus canadensis'
stacked_dat <- droplevels(within(stacked_dat, species.cor[species.cor == "Taxus baccata "] <- 'Taxus baccata')) 

sp <- levels(as.factor(stacked_dat$species.cor));length(sp) ## 5947
input <- sp %>%
  tibble::as_tibble_col(column_name = 'binomial') %>%
  tidyr::separate(col = 'binomial', into = c('Genus', 'Species'))%>%
  dplyr::select(Genus, Species) %>%# select columns
  dplyr::distinct(Genus, Species) %>%# remove duplicate binomials
  dplyr::mutate(Genus = stringr::str_to_title(Genus)) %>%# capitalize Genus
  dplyr::mutate(Species = stringr::str_remove(Species, ".*?\\s")) %>%# remove everything before first space
  tidyr::drop_na(c('Genus', 'Species')) %>%# remove rows with NA's
  dplyr::arrange(Genus, Species) # sort names

result <- input %>% 
  matching() %>% 
  enforce_matching('WCVP')
  ## match using sequential backbones. WCVP is the backbone used by phylo.tree 
result %>% summarize_output() 

sp.cor <- as.data.frame(cbind(result$Orig.Species,result$matched,
                              paste(result$Orig.Genus,result$Orig.Species),
                              paste(result$Matched.Genus,result$Matched.Species))) ## 

colnames(sp.cor)[colnames(sp.cor) == "V1"] <- "old.sp" ## old species names, used to remove indet and spec species from the dataset below
colnames(sp.cor)[colnames(sp.cor) == "V2"] <- "match" ## matched in Treemendous
colnames(sp.cor)[colnames(sp.cor) == "V3"] <- "species.cor" ## old name of species
colnames(sp.cor)[colnames(sp.cor) == "V4"] <- "species.name" ## corrected species name, based on Treemendous
stacked_dat <- merge(stacked_dat,sp.cor,by = 'species.cor')
unknown <- c("indet", "spec", "x") ## need to remove unknown species
stacked_dat <- subset(stacked_dat,!old.sp%in%unknown & match=="TRUE") ## drops from 5340 total sp


## remove indet species, which are unidentified species
# stacked_dat <- droplevels(subset(stacked_dat, family.cor !="Arecaceae")) ## remove palms from dataset
# stacked_dat <- droplevels(subset(stacked_dat, family.cor !="Musaceae")) ## remove banana and other herb. muscease
# stacked_dat <- droplevels(subset(stacked_dat, family.cor !="Cyatheaceae"))## remove tree ferns
# length(levels(as.factor(stacked_dat$species.cor))) #reduced 5,499 species 
# length(levels(as.factor(stacked_dat$genus.cor))) # reduced to 959 genuses
# 

## exclude species with low sample sizes
spp.n <- tapply(stacked_dat$d,stacked_dat$species.name,length)
spp <- levels(as.factor(stacked_dat$species.name)); length(spp) 
tree.n <- matrix(NA,nrow=length(spp))
spp.x.g.list <- as.matrix(stacked_dat$species.name)
tmt.tree.id.list <- as.matrix(stacked_dat$tmt.tree.id)
str(tmt.tree.id.list)
for (i in 1:length(spp)) {
  tree.n[i] <- length(unique(tmt.tree.id.list[spp.x.g.list==spp[i]]))
}

cbind(tree.n,spp.n)[tree.n<=10,]
nrow(cbind(tree.n,spp.n)[spp.n<100,]) 
nrow(cbind(tree.n,spp.n)[tree.n>=100,]) 

minNtree=20; minNpertree=5; # Set these acceptance thresholds
minNtree*minNpertree # min Nsamples
in.names <- spp[which(as.matrix(spp.n)>=minNtree*minNpertree & tree.n>=minNtree)]
length(in.names) ## 1491
stacked_dat <- subset(stacked_dat,stacked_dat$species.name%in%in.names)
tapply(stacked_dat$d,stacked_dat$species.name,length)
range(tapply(stacked_dat$d,stacked_dat$species.name,length))
length(unique(stacked_dat$species.name)) ## 1491 species

spp <- levels(as.factor(stacked_dat$species.name)); length(spp) 
# Check for immortal trees in grids
minNdeaths <- 5 # Set the minimum number of deaths to include spp.x.grid
weak.links <- c()
for (i in 1:length(spp)) {
  dat <- subset(stacked_dat,stacked_dat$species.name==spp[i])
  survTplus1 <- abs(dat$tree.status2-1)
  if(sum(dat$tree.status2)==0) print('Potential issue')
  if(sum(dat$tree.status2)<=minNdeaths) {
    print(paste('Only',sum(dat$tree.status2),'deaths observed'))
    weak.links <- c(weak.links,dat$species.name[1])
  }
}
length(weak.links)/length(spp) ## 22% of the species 
sum(stacked_dat$species.name%in%weak.links)/sum(!stacked_dat$species.name%in%weak.links) ## 2.7% of the obs. 
tapply(stacked_dat$species.name,stacked_dat$database.code,length)

# subset the weak links.
stacked_dat <- subset(stacked_dat,!species.name%in%weak.links)
tapply(stacked_dat$d,stacked_dat$database.code,length)
length(unique(stacked_dat$species.name))  
spp <- levels(as.factor(stacked_dat$species.name)); length(spp)

## FGE   FIAN  FIANE  FIANW  FIARM   FIAS    FPN    NAL    NBC    NQU    NSA    SYN 
## 180002 303264 170424  96804  63917 374777 633406 313091 702018 252176  84207  24070 

stacked_dat$Resolve_Biome <- round(stacked_dat$Resolve_Biome,digit=0)
tapply(stacked_dat$tmt.tree.id,stacked_dat$Resolve_Biome,length)
length(unique(stacked_dat$species.name)) ## 1153 species
spp <- levels(as.factor(stacked_dat$species.name)); length(spp) 
stacked_dat$gridID = as.factor(stacked_dat$hexgrid)# Maximum number of grids observed for any species
freq.tab = tapply(!is.na(stacked_dat$d),list(stacked_dat$species.name,
                                             stacked_dat$gridID),sum)
freq.tab[is.na(freq.tab)] = 0
n.grid <- rowSums(freq.tab>0)

# Variable interval length | One complicating aspect of this database is that 
# trees were censused at different time intervals. What we ultimately want are 
# growth and survival metrics over a standard time interval, whether this be 
# one, five, ten, etc. years. To do this we need to standardize the time step. 
# Here we explore a couple ways of doing this.  
# 1 | Convert size to annual time step using the annual relative growth rate
# [(s_t/s_0)^1/t] # s_0 is initial size and s_t is size at future time t. 
# t is therefore the time interval
# Example
#s0 = 10 # initial size of 10cm
#s5 = 20 # size in five years is 20cm
#rgr = (s5/s0)^(1/5); rgr # Rel. Growth Rate is 1.148, i.e., a 14.8% annual growth rate
#s0*rgr^5 # At this rgr, we get 20cm in five years of growth. s0*rgr^t recapitulates St
#s0*rgr # one year of growth is 11.48, ** This ignores the intercept **
rgr = ((stacked_dat$d2)/(stacked_dat$d))^(1/stacked_dat$interval)
cbind(stacked_dat$d2,stacked_dat$d*rgr^stacked_dat$interval)[100:115,] 
# This annual rgr recapitulates future size from current size
stacked_dat$d1<- stacked_dat$d * rgr # Calculate a new size at only one year.
#plot(stacked_dat$d,stacked_dat$d1)
#abline(1,1)

# Exclude extremely large growth or shrinkage
# hist(rgr)
quantile(rgr,probs=c(0.0001,0.001,0.01,0.025,0.5,0.975,0.99,0.999,0.9999),na.rm=TRUE)
sum(!is.na(rgr) & rgr==1.000)
summary(rgr)
sum(!is.na(rgr) & rgr>=1.5)
sum(!is.na(rgr) & rgr<=0.95)
unique(stacked_dat$interval[!is.na(rgr) & rgr>=1.5])
unique(stacked_dat$database.code[!is.na(rgr) & rgr>=1.15])
unique(stacked_dat$database.code[!is.na(rgr) & rgr<=0.95])
rgr.thresh = quantile(rgr,probs=c(0.001,0.999),na.rm=TRUE)
stacked_dat$d1[!is.na(rgr) & rgr<rgr.thresh[1] & rgr>rgr.thresh[2]] <- NA
summary(stacked_dat$d1)
# save(stacked_dat,spp,file='./data/TreeData_QC.Rdata') 
#___________________________________________________________________________----
#--| Growth |-------------------------------------------------------------->----
# rm(list=ls(all=TRUE))
# load('./data/TreeData_QC.Rdata') # Load Quality controlled data set, stacked_dat and spp
# 2 | Create a nonlinear function that directly computes the annual slope and 
# intercept of the growth function that allows for variable time intervals. 
# Method 1 above implicitly assumes a zero intercept in the growth function. 
# I.e., future size is simply current size multiplied by the rate of growth. 
# However, growth functions typically have a nonzero intercept. So that the 
# annual intercept value is included in the growth multiplication in the 
# subsequent year. We worked through a growing series algebraically and found  
# the following function, which we can fit directly using the non-linear least 
# squares function, nls. 
nls(log(d2)~ b0*( (1-b1^(interval)) / (1-b1) ) + log(d)*b1^interval,data=dat, 
    start=list(b0=1,b1=0.5),control=list(maxiter = 1e3, tol = 1e-06, 
                                         minFactor = 1e-6,warnOnly = TRUE))

# plotting parameters
cols <- c('black','gray40','magenta','gray80','cyan4','cyan1'); 
cexlabs <- 1.7; cextits <- 1.1; pchs <- 21; cexs <- 1.5; lwds <- 2; 
cexaxs <- 1.1; ltys <- 2; cexcoef <- 0.9
w=8.5; h=11; nrows=6; ncols=4;
newsheet <- seq(1,length(spp),ncols*nrows)

# graphics.off(); x11(width=w,height=h) #
# pdf(file='Figures/GrowthFunction_GraphicalSummary_Book.pdf',
#     width=w, height=h, onefile=TRUE)
# Growth loop
fitg_B <- fitg_B1 <- fitg_B2 <- fitg_Bp <- matrix(NA,nrow=length(spp),ncol=2)
fitg_Bg <- matrix(NA,nrow=length(spp),ncol=max(n.grid))
fitg_sigma <- fitg_r2 <- matrix(NA,nrow=length(spp),ncol=1)
fitg_F <- matrix(NA,nrow=length(spp),ncol=4)
min.sz <- max.sz <- matrix(NA,nrow=length(spp),ncol=1)
newSize <- matrix(NA,nrow=length(spp),ncol=50)
for (i in 1:length(spp)) {
  dat <- subset(stacked_dat,stacked_dat$species.name==spp[i])
  gridLevs = levels(droplevels(dat$gridID))
  dat <- droplevels(subset(dat,!is.na(dat$d1))) ## for 1-yr time step
  if(i%in%newsheet) {
   graphics.off(); quartz(width=w,height=h); 
   #graphics.off(); x11(width=w,height=h) #
    par(mfrow=c(nrows,ncols),oma=c(3,3.5,1,1),mar=c(1.75,1.75,0.5,0.5))
  }
  # Fit Growth Model [with random effects]
  if(n.grid[i]>1 & length(levels(dat$gridID))>1) {
    fitg <- lmer(log(d1)~log(d)+(1|gridID),na.action=na.omit,data=dat) 
    fitg_B[i,1] <- fixef(fitg)[1]
    fitg_B[i,2] <- fixef(fitg)[2]
    fitg_r2[i,1] <- r.squaredGLMM(fitg)[1] ## added this line April 2023 for R2
    r.int = as.matrix(ranef(fitg)$gridID+fixef(fitg)[1])
    if(length(r.int)!=n.grid[i]) { # This accounts for spp that grids with no surviving plants. These have survival estimates, but no growth parameters.
      notthese = which(!gridLevs %in% rownames(r.int))
      fitg_Bg[i,(1:n.grid[i])[-notthese]] <- r.int; rm(r.int)
    } else fitg_Bg[i,1:n.grid[i]] <- r.int; rm(r.int)
  } else {
    fitg <- lm(log(dat$d1)~log(dat$d),na.action=na.omit) 
    # growth on a one year timestep
    # regression parameters for growth model
    fitg_B[i,1] <- fitg$coeff[1]
    fitg_B[i,2] <- fitg$coeff[2]
    fitg_Bg[i,1] <- fitg$coeff[1] 
    fitg_r2[i,1] <- summary(fitg)$r.squared ### added this line April 2023 for R2
  }
  fitg_sigma[i] <- summary(fitg)$sigma
  min.sz[i] <- log(100)
  max.sz[i] <- max(c(log(dat$d),log(dat$d1)),na.rm=TRUE) # for one year time step
  newSize[i,] <- seq(min.sz[i],max.sz[i],length.out=50)
  plot(newSize[i,],newSize[i,],type='n',xlab='',ylab='',cex=cexs,
       cex.lab=cexlabs,cex.axis=cexaxs,pch=pchs,bg=cols[2])
  points(log(dat$d),log(dat$d1),pch=pchs,bg=cols[2]) ## for one year time-step
  abline(0,1,col=cols[6],lwd=lwds)
  title(main=spp[i],line=-1,cex.main=cextits)
  lines(newSize[i,],fitg_B[i,1]+fitg_B[i,2]*newSize[i,],col=cols[3], lwd=lwds)
  if(i%in%(newsheet[-1]-1)) {
    title(xlab=expression(ln(Size[t])),ylab=expression(ln(Size[t+1])),
          outer=TRUE,line=1.3,cex.lab=cexlabs)
    thisone = which((newsheet[-1]-1)==i)
    fileout = paste0('./figures/GrowthFunction/Gyx_DataFitPlot_',thisone,'.pdf')
    dev.copy(pdf,file=fileout,w,h);
    dev.off();
    rm(thisone,fileout)
  }
}; beepr::beep('fanfare')
title(xlab=expression(ln(Size[t])),ylab=expression(ln(Size[t+1])),
      outer=TRUE,line=1.3,cex.lab=cexlabs)
#fileout = paste0('figures/GrowthFunction/Gyx_DataFitPlot_',length(newsheet),'.pdf')
#dev.copy(pdf,file=fileout,w,h); dev.off();rm(fileout)

which(n.grid!=rowSums(!is.na(fitg_Bg))) # Check sum = 0

# Quick Fix | Some grids had no trees survive. We skipped them above, but we 
# would like to add in the species level mean estimates.
for (i in 1:length(spp)) {
  if(n.grid[i]!=sum(!is.na(fitg_Bg[i,]))) {
    dat <- subset(stacked_dat,stacked_dat$species.name==spp[i])
    gridLevs = levels(droplevels(as.factor(dat$gridID)))
    dat <- droplevels(subset(dat,!is.na(dat$d1)))
    gridSub = levels(droplevels(as.factor(dat$gridID)))
    these = which(!gridLevs %in% gridSub)
    fitg_Bg[i,these] <- fitg_B[i,1]
    rm(these,gridLevs,gridSub)
  } 
}; beepr::beep('fanfare')
which(n.grid!=rowSums(!is.na(fitg_Bg))) # Check sum = 0

save(fitg_B,fitg_Bg,n.grid,file='./data/TreeMort_Growth_1yr_woGrid.Rdata')
save(fitg_sigma,file='./data/TreeMort_sigmag_1yr_woGrid.Rdata') 
#___________________________________________________________________________----
#--| Survival |------------------------------------------------------------>----
# Tree status in this database is 0 for alive and 1 for death. The log-log link function we use below is a way of modeling mortality, so we will analyze deaths here, but below we model survival [survival = 1 - mortality].
summary(as.factor(stacked_dat$tree.status2))

# Polynomial regression
# poly function works same as brute force polynomial set up when raw = TRUE
# Complementary log log link function used to account for time interval differences. Time is included as offset.
norder = 2 # select order of polynomial, e.g., 2 = quadratic
fits_B <- fits_Bp <- matrix(NA,nrow=length(spp),ncol=(1+norder))
fits_Bg <- matrix(NA,nrow=length(spp),ncol=max(n.grid))
fits_df <- matrix(NA,nrow=length(spp),ncol=2)
mu.size <- matrix(NA,nrow=length(spp),ncol=1)
# save standard deviation (sigma) of the residuals
fits_sigma <- fits_r2 <- fits_X2p <- matrix(NA,nrow=length(spp),ncol=1)

# Polynomial regression model fit line function
poly.fit = function(x,B,nyrs) {
  nc = length(B)
  X = matrix(NA,nrow=length(x),ncol=nc)
  for(k in 1:nc) X[,k] = x^(k-1)
  p = X%*%B + log(nyrs)
  return(exp(-exp(p)))
}

# plotting parameters
cols <- c('black','gray40','magenta','gray80','cyan4','cyan1'); 
cexlabs <- 1.7; cextits <- 1.1; pchs <- 21; cexs <- 1.5; lwds <- 2; 
cexaxs <- 1.1; ltys <- 2; cexcoef <- 0.9
w=8.5; h=11; nrows=6; ncols=4;
newsheet <- seq(1,length(spp),ncols*nrows)
min.sz <- max.sz <- matrix(NA,nrow=length(spp),ncol=1)
newSize <- matrix(NA,nrow=length(spp),ncol=50)

# Loop over all spp and fit survival function
for (i in 1:length(spp)) {
  dat <- droplevels(subset(stacked_dat,stacked_dat$species.name==spp[i]))
  if(i%in%newsheet) {
    graphics.off(); quartz(width=w,height=h); 
    #graphics.off(); x11(width=w,height=h);
    par(mfrow=c(nrows,ncols),oma=c(3,3.5,1,1),mar=c(1.75,1.75,0.5,0.5))
  }
  # Fit Survival model [Polynomial with comp. log-log link and time offset][with random effects]
  mortTplus1 <- dat$tree.status2 ## for 1-yr timestep
  # Raw-dogging not a good idea here.
  mu.size[i] = mean(log(dat$d))
  center.d = log(dat$d) - mu.size[i]
  d.sq = center.d^2
  if(n.grid[i]>1) {
    m1 <- glm(mortTplus1~center.d+d.sq+offset(log(dat$interval)), 
              family=binomial(link="cloglog"))
    fits_r2 <- NagelkerkeR2(m1)[2]
    fits <- glmer(mortTplus1~center.d+d.sq+(1|dat$gridID)+
                    offset(log(dat$interval)),family=binomial(link="cloglog"),
                  start=list(fixef=coef(m1)),
                  control=glmerControl(nAGQ0initStep=FALSE)) ; rm(m1)
    # regression parameters and model statistics for survival model 
    fits_B[i,] <- fixef(fits)
    fits_r2[i] <- fits_r2
    fits_Bg[i,1:n.grid[i]] <- as.matrix(ranef(fits)$`dat$gridID`+fixef(fits)[1])
  } else {
    mortTplus1 <- dat$tree.status2 ## for 1-yr timestep
    fits <- glm(mortTplus1~center.d+d.sq+offset(log(dat$interval)), 
                family=binomial(link="cloglog"))
    fits_r2 <- NagelkerkeR2(fits)[2]
    # regression parameters and model statistics for survival model 
    # for(k in (0:norder)+1) {
    fits_B[i,] <- fits$coeff
    fits_Bg[i,1] <- fits$coeff[1] 
    fits_r2[i] <- fits_r2
  }
  # Plotting
  min.sz[i] <- log(100)
  max.sz[i] <- max(c(log(dat$d),log(dat$d1)),na.rm=TRUE) 
  newSize[i,] <- seq(min.sz[i],max.sz[i],length.out=50)
  os<-order(log(dat$d)); os.sizeT<-(log(dat$d))[os]; 
  os.survTplus1<-(abs(dat$tree.status2-1))[os];
  psz<-tapply(os.sizeT,as.numeric(cut(os.sizeT,50)),mean,na.rm=T)
  ps<-tapply(os.survTplus1,as.numeric(cut(os.sizeT,50)),mean,na.rm=T)
  plot(psz,ps,xlab='',ylab='',cex=cexs,cex.lab=cexlabs,cex.axis=cexaxs,
       pch=pchs,bg=cols[2],ylim=c(-0.12,1),xlim=c(min.sz[i],max.sz[i]))
  title(sub=spp[i],line=-1.25,cex.sub=cextits)
  mtime = median(dat$interval,na.rm=TRUE)
  lines(newSize[i,],poly.fit(newSize[i,]-mu.size[i],fits_B[i,],mtime),
        col=cols[3], lwd=lwds) # cloglog is l=log(-log(1-p)) instead of logit log(p/(1-p))
  if(i%in%(newsheet[-1]-1)) {
    title(xlab=expression(ln(Size[t])),ylab='Pr(Survival)',
          outer=TRUE,line=1.3,cex.lab=cexlabs)
    thisone = which((newsheet[-1]-1)==i)
    fileout = paste0('./figures/SurvivalFunction/Sx_DataFitPlot_',thisone,'.pdf')
    dev.copy(pdf,file=fileout,w,h);
    dev.off();
    rm(thisone,fileout)
  }
  rm(dat,center.d,d.sq,mortTplus1)
}; beepr::beep('fanfare'); print('You were right, Lala!')
title(xlab=expression(ln(Size[t])),ylab='Pr(Survival)',outer=TRUE,line=1.3,
      cex.lab=cexlabs)
fileout = paste0('./figures/SurvivalFunction/Sx_DataFitPlot_',length(newsheet),'.pdf')
dev.copy(pdf,file=fileout,w,h); dev.off();rm(fileout)

# Brief check
which(n.grid!=rowSums(!is.na(fits_Bg)))

save(fits_B,fits_Bg,mu.size,file='./data/TreeMort_Survival_1yr_woGrid.Rdata')
#___________________________________________________________________________----
#--| Species-level info needed for IPM |----------------------------------->----
## Subset data by spp and save in list of datasets
# load('./data/TreeData_QC.Rdata') # Load Quality controlled data set, stacked_dat and spp
stacked_dat$gridID = as.factor(stacked_dat$hexgrid)
freq.tab = tapply(!is.na(stacked_dat$d),list(stacked_dat$species.name,
                                             stacked_dat$gridID),sum)
freq.tab[is.na(freq.tab)] = 0
n.grid <- rowSums(freq.tab>0)

# Data wrangling 
min.sz <- max.sz <- matrix(NA,nrow=length(spp),ncol=1)
min.sz[] <- log(100)
info.spp = data.frame("species.name"=spp,"species.name"=NA, "database.code"=NA,
                      "old.sp"=NA,"family.cor"=NA)
thesecol = c("species.name","database.code","old.sp","family.cor")

GridMat = BiomeMat = matrix(NA,nrow=length(spp),ncol=max(n.grid))
# Set size thresholds for passage times
yf.spp = hf.spp = matrix(NA,nrow=length(spp),ncol=2)
dimnames(yf.spp) <- dimnames(hf.spp) <- list(spp,c('70th','90th'))

for (i in 1:length(spp)) {
  dat <- subset(stacked_dat,stacked_dat$species.name==spp[i])
  info.spp[i,thesecol] <- dat[1,c("species.name","database.code",
                                          "old.sp","family.cor")]
  GridMat[i,1:n.grid[i]] = as.numeric(levels(as.factor(dat$hexgrid)))
  for(g in 1:n.grid[i]) {
    BiomeMat[i,g] = dat$Resolve_Biome[dat$hexgrid==GridMat[i,g]][1]
  }
  for(j in 1:2) {
    hf.spp[i,j] = quantile(dat$height,prob=c(0.7,0.9)[j],na.rm=TRUE)
    yf.spp[i,j] = quantile(dat$d,prob=c(0.7,0.9)[j],na.rm=TRUE)
  }
  dat <- subset(dat,!is.na(dat$d1))
  max.sz[i] <- max(c(log(dat$d),log(dat$d1)),na.rm=TRUE)
  rm(dat)
  if(i%in%floor(seq(115,length(spp),115)))
    print(paste('Loop is',round(i/length(spp),2)*100,'% complete!'))
}

save(min.sz,max.sz,info.spp,GridMat,BiomeMat,yf.spp,hf.spp,spp,
      file='Min_Max_Size.Rdata')
#___________________________________________________________________________----
