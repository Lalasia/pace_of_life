################################################################################
################  5| Analyse environmental effect on LHTs   ####################
################################################################################

rm(list=ls(all=TRUE)) ## Clear the workspace

# library(devtools)
# devtools::install_github("jinyizju/V.PhyloMaker2")
# Installation of ggtree: https://bioconductor.org/install/
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18")
# BiocManager::install(c("ggtree"))

library(stringr) # to split species names
library(V.PhyloMaker2) # to build phylogenetic trees
library(MCMCglmm) # for bayesian models
library(phangorn)  # for nnls.tree to make ultrameric
library(ggtree) # for adding color to phylogeny plot
library(ggplot2) # needed for some plots
library(corrplot) # for correlation matrix matrix

setwd("~/[set wd]")
trait_env <- read.csv(file='./data/LHT_x_Env_dataset.csv')
################################################################################
## Add phylogeny
################################################################################
trait_env$genus <- str_split_fixed(trait_env$species.name, " ", 2)[,1]

dat = data.frame(species=trait_env$species.name,genus=trait_env$genus,
                 family=trait_env$Family)
# phyloTree <- phylo.maker(sp.list=dat,scenarios= c("S3"))
# Error on initial run: "Taxonomic classification not consistent between sp.list and tree."
# Leguminosae   Fabaceae
# Compositae    Asteraceae --> Update these to the current family names
trait_env[dat$family=='Leguminosae','Family'] = 'Fabaceae'
trait_env[dat$family=='Compositae','Family'] = 'Asteraceae'
# Fix these genus specific family classifications
# Genus         Old Family        Change to this family
# Poraqueiba    Icacinaceae       Metteniusaceae
# Nyssa         Cornaceae         Nyssaceae
# Gallesia      Phytolaccaceae    Petiveriaceae
# Dendrobangia  Cardiopteridaceae Metteniusaceae
# Chaetocarpus  Peraceae          Euphorbiaceae
# Calophyllum   Clusiaceae        Calophyllaceae
trait_env[dat$genus=='Poraqueiba' & 
            dat$family=='Icacinaceae','Family'] = 'Metteniusaceae'
trait_env[dat$genus=='Nyssa' & 
            dat$family=='Cornaceae','Family'] = 'Nyssaceae'
trait_env[dat$genus=='Gallesia' & 
            dat$family=='Phytolaccaceae','Family'] = 'Petiveriaceae'
trait_env[dat$genus=='Dendrobangia' & 
            dat$family=='Cardiopteridaceae','Family'] = 'Metteniusaceae'
trait_env[dat$genus=='Chaetocarpus' & 
            dat$family=='Peraceae','Family'] = 'Euphorbiaceae'
trait_env[dat$genus=='Calophyllum' & 
            dat$family=='Clusiaceae','Family'] = 'Calophyllaceae'
# Error on initial run: "Note: 8 taxa fail to be binded to the tree,"
# [1] "Sclerolobium_paniculatum"   "Sclerolobium_chrysophyllum"
# [3] "Sclerolobium_eriopetalum"   "Sclerolobium_melinonii"    
# [5] "Sclerolobium_bracteosum"    "Lecointea_amazonica"       
# [7] "Clathrotropis_macrocarpa"   "Marmaroxylon_basijugum" 
trait_env[dat$species=='Sclerolobium paniculatum','species.name'] = 'Tachigali paniculata'
trait_env[dat$species=='Sclerolobium chrysophyllum','species.name'] = 'Tachigali chrysophylla'
trait_env[dat$species=='Sclerolobium eriopetalum','species.name'] = 'Tachigali eriopetala'
trait_env[dat$species=='Sclerolobium melinonii','species.name'] = 'Tachigali melinonii'
trait_env[dat$species=='Sclerolobium bracteosum','species.name'] = 'Tachigali bracteosa'
trait_env[dat$species=='Marmaroxylon basijugum','species.name'] = 'Zygia basijuga'

dat = data.frame(species=trait_env$species.name,genus=trait_env$genus,
                 family=trait_env$Family)
spp.dat <- data.frame(dat[!duplicated(dat$species),1:3],
                  'species.relative'=NA,'genus.relative'=NA)
prspp = c("Clathrotropis macrocarpa","Lecointea amazonica",
          "Tachigali eriopetala","Tachigali melinonii","Zygia basijuga")
spp.dat$genus.relative[spp.dat$species=="Clathrotropis macrocarpa"] = "Ormosia"
spp.dat$genus.relative[spp.dat$species=="Lecointea amazonica"] = "Uribea"
spp.dat$genus.relative[spp.dat$species=="Tachigali eriopetala"] = "Tachigali"
spp.dat$genus.relative[spp.dat$species=="Tachigali melinonii"] = "Tachigali"
spp.dat$genus.relative[spp.dat$species=="Zygia basijuga"] = "Zygia"
# Placed with sister genera (Cardoso et al. 2012, Fabaceae phylogeny)
spp.dat[spp.dat$species%in%prspp,]

# Uncomment and run for first run through ##
# rel <- bind.relative(spp.dat)
# phyloTree <- phylo.maker(sp.list=rel$species.list, tree=rel$phylo,
#                          nodes=rel$nodes.info, scenarios= c("S3"))
# phyloTree$species.list[phyloTree$species.list$status=='fail to bind',]
# phyloTree$species.list[phyloTree$species.list$species%in%prspp,]
# write.tree(phyloTree$scenario.3, "./data/phylogeny.tre")
phy = read.tree("./data/phylogeny.tre")
phy$tip.label

# Plot phylogeny 
w = 7; h = 7
pdf("./figures/phyloTree.pdf",width=w,height=h) 
plot.phylo(phy, cex = 0.2,type = "fan") 

circ <- ggtree(phy, layout = "circular")
df <- data.frame(trait_env$species,trait_env$LE_10)
colnames(df) <- c("sp","LE.10")
df2 <- aggregate(list(LE.10=df$LE.10), by=list(df$sp),FUN=mean,na.rm=TRUE);
colnames(df2) <- c("sp","LE.10")
df2 <-df2[order(match(df2$sp,phy$tip.label)),]
rownames(df2) <- phy$tip.label
df2 <- data.frame(df2[,2]);colnames(df2) <- c("LE.10")
rownames(df2) <- phy$tip.label

pdf("./figures/phyloTree2.pdf",width=w,height=h) 

gheatmap(circ, df2, offset=0, width=.2) +
  scale_fill_viridis_c(option = "plasma",ggtitle("Life exp from 10cm"))
dev.off(); dev.off();

# post process phylogeny
trait_env$animal = gsub(' ','_',trait_env$species.name)
# Add animal column becuase MCMCglmm requires the species column to be called 
# 'animal' to identify it was a phylo model 
phy$tip.label
length(phy$tip.label) # 1123
length(unique(phy$tip.label)) # 1123
length(unique(trait_env$animal)) # 1123
# Check node labels. All must be unique for inverseA function
length(unique(phy$node.label)) # 43
length(phy$node.label) # 858
phy$node.label = 1:length(phy$node.label)

which(!phy$tip.label%in%trait_env$animal) # Check sum = 0
which(!trait_env$animal%in%phy$tip.label) # Check sum = 0
phy$tip.label[phy$edge.length==0]
mcmctree <- nnls.tree(cophenetic(phy),phy,method='ultrametric')

sum(is.na(mcmctree$tip.label))
# Ainv<-inverseA(mcmctree,nodes='ALL')$Ainv
INphylo = inverseA(mcmctree,nodes='ALL')
INphylo$pedigree 
is.rooted(mcmctree)

################################################################################
## Phylogenetic mixed effects model
################################################################################
# Scale variables for mixed effects modeling.
names(trait_env)
head(trait_env)
trait_env$sLE_10 <- scale(log(trait_env$LE_10))                  
trait_env$stau_20 <- scale(trait_env$tau_20)                         
trait_env$sLE.f_20 <- scale(log(trait_env$LE.f_20))                         
trait_env$stau_70 <- scale(trait_env$tau_70)                   
trait_env$sLifespan.age <- scale(log(trait_env$Lifespan.age))
trait_env$sLifespan.size <- scale(trait_env$Lifespan.size) # already log transformed
ntraits = 4
priorNull <- list(R = list(V = diag(ntraits)*0.02, nu = 7), 
                   G = list(G1 = list(V = diag(ntraits), nu = 0.002)))
priors<-list(R=list(V=diag(ntraits)*1, nu=0.002), 
             G=list(G1=list(V=diag(ntraits)*1, nu=0.002)))

## units is the response variable value, and trait is the response variable name, 
# which corresponds to the categories. By specifying rcov = ~us(trait):units, 
# you are allowing the residual variance to be heterogeneous across "traits" 
# (response categories) so that all elements of the residual variance-covariance 
# matrix will be estimated.
mcmcdata = trait_env
names(mcmcdata)
# Check sum = 0
sum(is.na(mcmcdata[,c('stau_20','sLE.f_20',"stau_70",
                      'sLifespan.age','sLifespan.size',
                      'NPP','Growing_season','SG_SOC_Content_015cm', 
                      'CHELSA_BIO_Annual_Precipitation',
                      'CHELSA_BIO_Temperature_Seasonality')]))

# Set MCMC parameters
burns <- 50; thins <- 1; nitts <- 1.0*10*thins + burns # for 1000 MCMC samples
burns <- 5000; thins <- 10; nitts <- 1.0*1000*thins + burns # for 1000 MCMC samples
nsam <- (nitts - burns) / thins; nsam; nitts
mod_full <- MCMCglmm(cbind(sLE_10,stau_20,stau_70,sLifespan.size)
                  ~ trait + 
                    trait:soil_pca +  ## I removed NPP here: was ~ scale(NPP) + scale(Growing_season) + 
                    trait:precip_pca + 
                    trait:temp_pca -1,
                  random= ~ us(trait):animal, rcov = ~ us(trait):units,
                  family = c("gaussian","gaussian", "gaussian",
                             "gaussian"), 
                  pedigree = INphylo$pedigree, #mcmctree, # phylogeny
                  # ginv = list(animal = Ainv),
                  data = mcmcdata,scale = TRUE,nitt = nitts,burnin = burns,
                  thin = thins,prior = priorNull)

summary(mod_full)
priorNull2 <- list(R = list(V = diag(ntraits)*0.02, nu = 7))

mod_full_noPhy <- MCMCglmm(cbind(sLE_10, stau_20,stau_70, sLifespan.size)
                        ~ trait +
                          trait:soil_pca + 
                          trait:precip_pca + 
                          trait:temp_pca -1,
                         rcov = ~ us(trait):units,
                         family = c("gaussian", "gaussian",
                                    "gaussian", "gaussian"),
                         data = mcmcdata,
                         scale = TRUE,
                         nitt = nitts,
                         burnin = burns,
                         thin = thins, 
                         prior = priorNull2)
mod_full$DIC # -8863.892 | 6 Feb 2024
mod_full_noPhy$DIC # 67968.33  | 6 Feb 2024

mod_full$DIC-mod_full_noPhy$DIC # -76832.22  | 6 Feb 2024
summary(mod_full)
summary(mod_full_noPhy)

# Function ------->
# Plot parameter estimates | plots means as points with error bars:
# Receives output from MCMCglmm object
plot.estimates <- function(x) {
  if (class(x) != "summary.mcmc")
    x <- summary(x)
  n <- dim(x$solutions)[1]
  par(mar=c(2, 20, 4, 1))
  plot(x$solutions[,1], n:1,
       yaxt="n", ylab="",
       xlim=range(x$solutions[,2:3])*1.2,
       pch=19,
       main="Posterior means and 95% credible intervals")
  grid()
  axis(2, at=n:1, rownames(x$solutions), las=2)
  arrows(x$solutions[,2], n:1, x$solutions[,3], n:1, code=0)
  abline(v=0, lty=2)
}

quartz(7,5)
x11(width=7, height=5)
plot.estimates(mod_full) # function defined above.
# Convergence diagnostics (copied from Titus von der Malsburg, MCMCglmm-intro on GitHub)
quartz(8,16);par(mfrow=c(5,2), mar=c(2,2,1,0))
x11(width=8,height=16);par(mfrow=c(5,2), mar=c(2,2,1,0))
plot(mod_full$Sol, auto.layout=F)
modfullest <- posterior.mode(mod_full$Sol)
HPDinterval(mod_full$VCV)

autocorr.plot(mod_full$Sol) # great drop off.
autocorr.plot(mod_full$VCV) # great drop off
summary(effectiveSize(mod_full$Sol)) # 1000 minimum, mostly = 1,000, up to 1,117
effectiveSize(mod_full$VCV) # minimum 845.
# Phylogenetic heritablity
# Collect the posterior modes for each variance
dim(mod_full$VCV)
dimnames(mod_full$VCV)
Gi = 1:(ntraits*ntraits) # indexes G matrix estimates
Ri = ntraits*ntraits+1:(ntraits*ntraits) # indexes residual matrix estimates
GvT <- matrix(posterior.mode(mod_full$VCV[,Gi]),
              nrow=ntraits,ncol=ntraits)
traitnames <- c('LifeEx_10','tau_20','tau_70','LifeSp_size')
dimnames(GvT) = list(traitnames,traitnames)
# Residual variance-covariance
RvT <- matrix(posterior.mode(mod_full$VCV[,Ri]),
              nrow=ntraits,ncol=ntraits)
dimnames(RvT) = list(traitnames,traitnames)
# Total Phenotypic variance-covariance
PvT <- GvT + RvT;
herit <-
  mod_full[["VCV"]][ ,Gi] /
  (mod_full[["VCV"]][ ,Gi] + mod_full[["VCV"]][ ,Ri])
effectiveSize(herit)
lambda = cbind(posterior.mode(herit)[1+seq(0,ntraits*ntraits-1,ntraits+1)],
              HPDinterval(herit)[1+seq(0,ntraits*ntraits-1,ntraits+1),]); lambda
effectiveSize(mod_full[["VCV"]][ , 1:25])
posterior.mode(mod_full[["VCV"]][ , 1:25])

# Compute correlations
GcT <- GvT/sqrt(diag(ntraits)%*%diag(GvT)%*%diag(GvT))
RcT <- RvT/sqrt(diag(ntraits)%*%diag(RvT)%*%diag(RvT))
PvT <- PvT/sqrt(diag(ntraits)%*%diag(PvT)%*%diag(PvT))

# Trait values
ZT <- posterior.mode(mod_full$Sol)[1:ntraits] 
GvT
RvT
PvT
GcT
RcT
PvT/sqrt(diag(ntraits)%*%diag(PvT)%*%diag(PvT))
ZT
# save(list=ls(all=TRUE),file='./data/Env_Phylo_GLMM_6Feb2024.Rdata')
load('./data/Env_Phylo_GLMM_6Feb2024.Rdata')

#--Individual trait models------------------------------------------------------|
plot.add.est <- function(x,ran,cols) {
  n <- dim(x)[1]
  par(mar=c(1, 1, 0, 1))
  plot(x[,1], n:1,yaxt="n", ylab="",bty='n',bg='gray80',
       xlim=ran*1.0, xaxt='n',pch=19,main="",col=cols,cex=cexs)
  grid()
  arrows(x[,2], n:1, x[,3], n:1, code=0,col=cols,lwd=lwds)
  abline(v=0, lty=2)
}

## graphing fixed effects
w = 5; h = 6; cexs=1.5; lwds=2
cols <- c("dodgerblue","dodgerblue3","dodgerblue4")
env.names = c('L.E. from 10cm','Yrs to 20cm','Yrs to 70th%','Size max. age')
x11(width=w,height=h); par(mfrow=c(ntraits,1),oma=c(4,4,2,2))
quartz(width=w,height=h); par(mfrow=c(ntraits,1),oma=c(4,4,2,2))
# Order 2, 3, 4, 1
x = summary(mod_full)$solutions[(ntraits+1):(ntraits*4),1:3]
ran = range(x)
plot.add.est(x[seq(2,(ntraits*3),ntraits),],ran,cols) # function defined above.
mtext(env.names[2],side=2,line=1)
for(n in 3:1)
  text(ran[1]+0.01,n,c('Temp PCA','Precip PCA','Soil PCA')[n],col=cols[3:1][n])
for(i in c(3,4,1)) {
  plot.add.est(x[seq(i,(ntraits*3),ntraits),],ran,cols) # function defined above.
  mtext(env.names[i],side=2,line=1)
}
axis(1,pretty(ran))
mtext('Model coefficients',side=1,line=3)
title(main='Phylogenetic model',outer=TRUE)
dev.copy(pdf,'./figures/CoEff_plot_PCAs_PhyCor.pdf', width=w, height=h)
dev.off()
# non phylogenetic model
x = summary(mod_full_noPhy)$solutions[(ntraits+1):(ntraits*4),1:3]
ran = range(x)
x11(width=w,height=h); par(mfrow=c(ntraits,1),oma=c(4,4,2,2))
quartz(width=w,height=h); par(mfrow=c(ntraits,1),oma=c(4,4,2,2))
# Order 2, 3, 4, 1
plot.add.est(x[seq(2,(ntraits*3),ntraits),],ran,cols) # function defined above.
mtext(env.names[2],side=2,line=1)
for(n in 3:1)
  text(ran[1]+0.01,n,c('Temp PCA','Precip PCA','Soil PCA')[n],col=cols[3:1][n])
for(i in c(3,4,1)) {
  plot.add.est(x[seq(i,(ntraits*3),ntraits),],ran,cols) # function defined above.
  mtext(env.names[i],side=2,line=1)
}
axis(1,pretty(ran))
mtext('Model coefficients',side=1,line=3)
title(main='Base model',outer=TRUE)
dev.copy(pdf,'./figures/CoEff_plot_PCAs_noPhy.pdf', width=w, height=h)
dev.off()

summary_fullMod <-summary(mod_full)$solutions
#write.csv(summary_fullMod, "./output/summary_fullMod.csv",row.names=TRUE)
#---------------------------------------------------------------------
phylo_ests <- data.frame(lambda)
colnames(phylo_ests) <- c("post.mean", "Lower_95_CI","Upper_95_CI")

phylo_ests$Trait <- c("Life exp from 10 cm",
                      "Yrs to 20 cm",
                      "Yrs to 70th%",
                      "Size at maximal age")

phylo_ests$Trait <- ordered(phylo_ests$Trait, 
                            levels = c("Size at maximal age",
                                       "Life exp from 10 cm",
                                       "Yrs to 70th%",
                                       "Yrs to 20 cm"))
w=4; h=3
quartz(width=w,height=h);
x11(width=w,height=h);
par(mfrow=c(1,1),oma=c(2,2,1,1),mar=c(2,2,.2,.2))

zp1 <- ggplot(phylo_ests)
zp1 <- zp1 + geom_hline(yintercept = .50, colour = gray(2/3), lty = 2)
zp1 <- zp1 + geom_pointrange(aes(x = Trait, y = post.mean, 
                                 ymin = Lower_95_CI,
                                 ymax = Upper_95_CI),
                             lwd = 0.75, position = position_dodge(width = 1/2),
                             shape = 21, fill = "grey40", colour = "grey40")
zp1 <- zp1 + coord_flip() + theme_bw()
zp1 <- zp1 + theme(axis.title.y=element_blank())
zp1 <- zp1 + theme(axis.text.x = element_text(colour = "black"))
zp1 <- zp1 + theme(axis.text.y = element_text(colour = "black"))
print(zp1)
dev.copy(pdf,'./figures/PhyH2_PascalsLambda_plot.pdf',width=w,height=h)
dev.off()

###### plotting for figure 3
# names(trait_env)
# trait_env[1:30,c('Resolve_Biome','Biome_num')]
biocol <- adjustcolor(c("#3B99B1BF","#45ADA3BF","#7CBA96BF","#EACB2BBF","#E9AE1DBF","#E78F0ABF","#EA6900BF"),alpha.f=0.75)
trait_env$Resolve_Biome[trait_env$Resolve_Biome == 'Tropical grassland'] <- 'Tropical savanna'
trait_env$Resolve_Biome[trait_env$Resolve_Biome == 'Temperate grassland'] <- 'Temperate savanna'

region = c("Tropical moist","Tropical dry","Tropical savanna",
           "Temperate savanna","Temperate broadleaf","Temperate conifer",
           "Boreal")  
col_line='grey29' ## color for dashed line
trait_env <- trait_env[order(trait_env$species.name),]
trait_env <- trait_env[order(trait_env$Biome_num),]
lwd = 3; tks = -.02; tks.wd=.8; cex.leg =.8; sizeaxis =1
graphics.off(); w = 8; h = 7.5
quartz(width=w,height=h);
x11(width=w,height=h);
par(mfrow=c(2,2),oma=c(.4,.4,0.1,.1),mar=c(4,4,2,2))
rownames(PvT) <- c("LE", "Gr 20", "Gr 70th", "Size")
colnames(PvT) <- c("LE", "Gr 20", "Gr 70th", "Size")
corrplot(PvT, type="lower", tl.col="black",tl.cex = sizeaxis,addCoef.col = 1,diag = FALSE)
mtext('A: Phenotypic correlation',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
plot(trait_env$tau_20~(trait_env$temp_pca),col=biocol[trait_env$Biome_num], frame = FALSE,
     pch=16, ylab = "", xlab = "",log='y',yaxt="n",xaxt="n",ylim=c(2,200))
axis(1, at = c(-6,-4,-2,0,2,4),labels=TRUE)
axis(2, at = c(1,10,100),labels=TRUE)
axis(2, at = c(1:5,seq(10,200,10)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
mtext('Yrs to grow from 10-20cm dbh',side=2,line=2.5,outer=F,cex=sizeaxis)
mtext('B: Temp. effect on Gr 20',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
legend('bottomleft',region,col=biocol,pch=19,bty='n',cex=cex.leg)
mtext('Temperature index, PC Axis 1',side=1,line=2.5,cex=sizeaxis)
trait_env$Lifespan.size_cm <- exp(trait_env$Lifespan.size)*.1 ## need to first convert to natural scale, and then change from mm to cm
plot(trait_env$Lifespan.size_cm~trait_env$temp_pca,col=biocol[trait_env$Biome_num], frame = FALSE,
     pch=16, ylab = "", yaxt="n",xaxt="n",xlab = "",log='y')
axis(1, at = c(-6,-4,-2,0,2,4),labels=TRUE)
axis(2, at = c(1,10,100,1000),labels=TRUE)
axis(2, at = c(1:5,seq(10, 100,10),seq(100,1000,100)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
mtext('Maximal size (cm dbh)',side=2,line=2.5,outer=F,cex=sizeaxis)
mtext('C: Temp. effect on size',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
mtext('Temperature index, PC Axis 1',side=1,line=2.5,cex=sizeaxis)
plot(trait_env$LE_10~trait_env$temp_pca,col=biocol[trait_env$Biome_num], frame = FALSE,
     pch=16, ylab = "", xlab = "",log='y',yaxt="n",xaxt="n")
axis(1, at = c(-6,-4,-2,0,2,4),labels=TRUE)
axis(2, at = c(1,10,100,1000),labels=TRUE)
axis(2, at = c(1:5,seq(10, 100,10),seq(100,1000,100),seq(1000, 3000,1000)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
mtext('Life exp from 10cm dbh (yrs)',side=2,line=2.5,outer=F,cex=sizeaxis)
mtext('D: Temp. effects on LE',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
mtext('Temperature index, PC Axis 1',side=1,line=2.5,cex=sizeaxis)
dev.copy(pdf,'./figures/Fig3.pdf',width=w,height=h)
dev.off()

