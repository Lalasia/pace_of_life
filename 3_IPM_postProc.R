################################################################################
############  3| Compile all IPM runs and Life History output  #################
################################################################################
#--| Housekeeping |-------------------------------------------------------->----
rm(list=ls(all=TRUE)) ## Clear the workspace
setwd("~/[add wd here]")
# load min max size with species info.
load('data/Min_Max_Size.Rdata')
#___________________________________________________________________________----
#--| Compile set 10 |------------------------------------------------------>----
# Load and compile each block of 10 spp
# 4 Jan block 10 [900-1000] hit a snag. i = 962 has an immortal survival-growth 
# matrix. I reran code in blocks of 10, and for run # 7 I skipped 962. 
load('data/temp/TreeMort_SurvGrow_Kernels10_1.Rdata')
load('data/temp/TreeMort_LifeHistory_Output10_1.Rdata')

bP_spp=P_spp; bN_spp=N_spp
btau=tau; bLE.f=LE.f; bLE.0=LE.0; blifespan=lifespan
btau.g=tau.g; bLE.0.g=LE.0.g; bLE.f.g=LE.f.g; blifespan.g=lifespan.g
i.end = seq(910,1000,10); i.end # for iterations 901-1000
i.start = c(901,i.end[-10]+1); i.start # for iterations 901-1000
i.list = c(961,963:970) # for iterations 961,963-970

for(i in 2:10) {
  load(paste0('data/temp/TreeMort_SurvGrow_Kernels10_',i,'.Rdata'))
  load(paste0('data/temp/TreeMort_LifeHistory_Output10_',i,'.Rdata'))
  these = i.start[i]:i.end[i]
  bP_spp[these,,]=P_spp[these,,];   bN_spp[these,,]=N_spp[these,,]
  btau[these,]=tau[these,];         bLE.f[these,]=LE.f[these,];
  bLE.0[these,]=LE.0[these,];       blifespan[these,]=lifespan[these,]
  btau.g[these,,]=tau.g[these,,];   bLE.0.g[these,,]=LE.0.g[these,,];
  bLE.f.g[these,,]=LE.f.g[these,,]; blifespan.g[these,,]=lifespan.g[these,,]
}
# revert to original object names
P_spp=bP_spp; N_spp=bN_spp
tau=btau; LE.f=bLE.f; LE.0=bLE.0; lifespan=blifespan
tau.g=btau.g; LE.0.g=bLE.0.g; LE.f.g=bLE.f.g; lifespan.g=blifespan.g

# save(P_spp,N_spp,file='data/temp/TreeMort_SurvGrow_Kernels10.Rdata')
# save(tau,LE.f,LE.0,lifespan,tau.g,LE.0.g,LE.f.g,lifespan.g,
#      file='data/temp/TreeMort_LifeHistory_Output10.Rdata')
#___________________________________________________________________________----
#--| Compile full set |---------------------------------------------------->----
# Load and compile each block of 100 spp
# 4 Feb 2024: two species did not run, i = c(30, 962). These were skipped and 
# should be removed
load('data/temp/TreeMort_SurvGrow_Kernels1.Rdata')
load('data/temp/TreeMort_LifeHistory_Output1.Rdata')

bP_spp=P_spp; bN_spp=N_spp
btau=tau; bLE.f=LE.f; bLE.0=LE.0; blifespan=lifespan
btau.g=tau.g; bLE.0.g=LE.0.g; bLE.f.g=LE.f.g; blifespan.g=lifespan.g

i.end = c(seq(100,1100,100),length(spp)); i.end
i.start = c(1,i.end[-12]+1); i.start

for(i in 2:12) {
  load(paste0('data/temp/TreeMort_SurvGrow_Kernels',i,'.Rdata'))
  load(paste0('data/temp/TreeMort_LifeHistory_Output',i,'.Rdata'))
  these = i.start[i]:i.end[i]
  bP_spp[these,,]=P_spp[these,,];   bN_spp[these,,]=N_spp[these,,]
  btau[these,]=tau[these,];         bLE.f[these,]=LE.f[these,];
  bLE.0[these,]=LE.0[these,];       blifespan[these,]=lifespan[these,]
  btau.g[these,,]=tau.g[these,,];   bLE.0.g[these,,]=LE.0.g[these,,];
  bLE.f.g[these,,]=LE.f.g[these,,]; blifespan.g[these,,]=lifespan.g[these,,]
}
# revert to original object names
P_spp=bP_spp; N_spp=bN_spp
tau=btau; LE.f=bLE.f; LE.0=bLE.0; lifespan=blifespan
tau.g=btau.g; LE.0.g=bLE.0.g; LE.f.g=bLE.f.g; lifespan.g=blifespan.g

# save(P_spp,N_spp,file='data/TreeMort_SurvGrow_Kernels.Rdata')
# save(tau,LE.f,LE.0,lifespan,tau.g,LE.0.g,LE.f.g,lifespan.g,
#      file='data/TreeMort_LifeHistory_Output.Rdata')
#___________________________________________________________________________----
#--| Clean new data objects |---------------------------------------------->----
sum(is.na(tau[,1]))
which(is.na(tau[,1])) # here are the problem matrices. The P matrix has columns 
# that sum to one, i.e., immortality.
# Alibertia verrucosa Senefeldera macrophylla 
# 30                     962 
not.i = c(-30,-962) # we need to exclude these.

##---------------------------------------Need to remove two immortals.
not.i = c(-30,-962)

dim(tau)
dim(LE.f)
dim(LE.0)
dim(tau.g)
dim(LE.0.g)
dim(LE.f.g)
dim(lifespan)
dim(lifespan.g)

dim(min.sz)
dim(max.sz)
dim(info.spp)
dim(GridMat)
dim(BiomeMat)
dim(yf.spp)
dim(hf.spp)
length(spp)

tau = tau[not.i,]
LE.f = LE.f[not.i,]
LE.0 = LE.0[not.i,]
tau.g = tau.g[not.i,,]
LE.0.g = LE.0.g[not.i,,]
LE.f.g = LE.f.g[not.i,,]
lifespan = lifespan[not.i,]
lifespan.g = lifespan.g[not.i,,]

min.sz = min.sz[not.i]
max.sz = max.sz[not.i]
info.spp = info.spp[not.i,]
GridMat = GridMat[not.i,]
BiomeMat = BiomeMat[not.i,]
yf.spp = yf.spp[not.i,]
hf.spp = hf.spp[not.i,]
spp = spp[not.i]

# save(tau,LE.f,LE.0,lifespan,tau.g,LE.0.g,LE.f.g,lifespan.g,
#      file='data/TreeMort_LifeHistory_Output.Rdata')
# save(min.sz,max.sz,info.spp,GridMat,BiomeMat,yf.spp,hf.spp,spp,
#      file='data/Min_Max_Size.Rdata')

# Extra checks. Running these lines is not critical----------------------------| 
x11(); 
hist(LE.0[,1])
plot(lifespan[,1],LE.0[,1],ylim=c(0,8e3))
max(colSums(P_spp[1,,]))
tau.g[1:30,,]
LE.0.g[1,,]
LE.f.g[1,,]

for (i in (1:length(spp))[not.i]) {
  P = P_spp[i,,]
  if(max(colSums(P))>1.0) print(paste("Problem matrix. Species #",i,
                                      'has max column sums =',max(colSums(P))))
} # If all is well, nothing should print.

# Check biome coding: biome codes should be integer values
unique(BiomeMat[,1])
unique(BiomeMat[!is.na(BiomeMat)]) # No problems here.
#___________________________________________________________________________----
#--| Clean and merge data |------------------------------------------------>----
# library(dplyr)
# library(tidyr)
# If not already loaded, load the following.
load('data/TreeMort_LifeHistory_Output.Rdata')
load('data/Min_Max_Size.Rdata')
#-- Data Merge 1 | Species mean output ----------------------------------------|
GridMat_mean <- GridMat[,1]
BiomeMat_mean <- BiomeMat[,1]
dimnames(LE.0)
dimnames(yf.spp)
dat_mean = data.frame('species.name'=info.spp$species.name,
                      'grid.num'=c(GridMat_mean),
                      'biome'=c(BiomeMat_mean),
                      'DBH_70'=yf.spp[,1],
                      'DBH_90'=yf.spp[,2],
                      # 'Ht_70'=hf.spp[,1],
                      # 'Ht_90'=hf.spp[,2],
                      'tau_20'=tau[,1],
                      'tau_70'=tau[,2],
                      'tau_90'=tau[,3],
                      'LE_10'=LE.0[,1],
                      'LE.f_20'=LE.f[,1],
                      'LE.f_70'=LE.f[,2],
                      'LE.f_90'=LE.f[,3],
                      'Lifespan.age'=lifespan[,1],
                      'Lifespan.size'=lifespan[,2],
                      'Family'=info.spp$family.cor)
unique(dat_mean$biome) 
str(dat_mean)
tapply(dat_mean$biome,dat_mean$biome,length)
# 1, Tropical Moist Forests 
# 2, Tropical Dry Forest 
# -- not 3, Tropical Coniferous Forests 
# 4, Temperate Broadleaf Forests 
# 5, Temperate Conifer Forests 
# 6, Boreal Forests 
# 7, Tropical Grasslands 
# 8, Temperate Grasslands 
# 9, Flooded Grasslands 
# -- not 10,Montane Grasslands 
# -- not 11, Tundra 
# -- not 12, Mediterranean Forests 
# 13, Deserts 
# -- not 14, Mangroves

#-- Data Merge 2 | Species x grid output --------------------------------------|
dat_grid = data.frame('species.name'=rep(info.spp$species.name,ncol(BiomeMat)), ## use grid.num=="746407" to check. something off here.
                      'grid.num'=c(GridMat),
                      'biome'=c(BiomeMat),
                      'DBH_70'=rep(yf.spp[,1],ncol(BiomeMat)),
                      'DBH_90'=rep(yf.spp[,2],ncol(BiomeMat)),
                      # 'Ht_70'=rep(hf.spp[,1],ncol(BiomeMat)),
                      # 'Ht_90'=rep(hf.spp[,2],ncol(BiomeMat)),
                      'tau_20'=c(tau.g[,1,]),
                      'tau_70'=c(tau.g[,2,]),
                      'tau_90'=c(tau.g[,3,]),
                      'LE_10'=c(LE.0.g[,1,]),
                      'LE.f_20'=c(LE.f.g[,1,]),
                      'LE.f_70'=c(LE.f.g[,2,]),
                      'LE.f_90'=c(LE.f.g[,3,]),
                      'Lifespan.age'=c(lifespan.g[,1,]),
                      'Lifespan.size'= c(lifespan.g[,2,]),
                      'Family'=info.spp$family.cor)
unique(dat_grid$biome) 
tapply(!is.na(dat_grid$biome),dat_grid$biome,sum)   
# only 13 flooded grasslands, 14 tundra, and 1 Mediterranean.

#-- Data Merge 3 | Clean Biome indicator --------------------------------------|
dat_grid <- droplevels(subset(dat_grid, biome==1|biome==2| biome==4 |biome==5 | biome==6 | 
                                biome==7 |biome==8 |biome==12|biome==13))
dat_grid$biome <- as.factor(dat_grid$biome)
# 1, Tropical Moist Forests 
# 2, Tropical Dry Forest 
# 3, Tropical Coniferous Forests 
# 4, Temperate Broadleaf Forests 
# 5, Temperate Conifer Forests 
# 6, Boreal Forests 
# 7, Tropical Savanna 
# 8, Temperate Sanvanna 
# 9, Flooded Grasslands 
# 10, Montane Grasslands 
# 11, Tundra 
# 12, Mediterranean Forests 
# 13, Deserts 
# 14, Mangroves
BiomeLevs = c("Tropical moist","Tropical dry","Tropical coniferous",
              "Temperate broadleaf","Temperate conifer","Boreal",
              "Tropical savanna","Temperate savanna","Flooded grassland",
              "Montane grassland","Tundra","Mediterranean","Desert","Mangrove")
biome_order=c(1:2,7:8,4:5,6,9,11:14)
BiomeLevs[biome_order]
unique(BiomeLevs[dat_grid$biome])
# There should not be any NAs here. 

dat_grid$Resolve_Biome <- dat_mean$Resolve_Biome <- NA 
for(i in 1:14) {
  dat_grid$Resolve_Biome[dat_grid$biome==i] <- BiomeLevs[i] 
  dat_mean$Resolve_Biome[dat_mean$biome==i] <- BiomeLevs[i] 
}
dat_grid$Resolve_Biome <- factor(dat_grid$Resolve_Biome,levels=BiomeLevs[biome_order])
dat_mean$Resolve_Biome <- factor(dat_mean$Resolve_Biome,levels=BiomeLevs[biome_order])

region = BiomeLevs[biome_order][1:7]
sum(dat_grid$Resolve_Biome%in%region)

dat_grid <- droplevels(subset(dat_grid,Resolve_Biome %in% region))
dat_mean <- droplevels(subset(dat_mean,Resolve_Biome %in% region))

dat_grid$Biome_num = as.numeric(dat_grid$Resolve_Biome)
dat_mean$Biome_num = as.numeric(dat_mean$Resolve_Biome)

# Subset out the quasi-immortals
sum(is.infinite(dat_grid$Lifespan.age)) 
# 127 species-grid-cells measured more than 10K yrs old, after which we ceased 
# simulations --> Replace Inf with NA and exclude these.
unique(dat_grid$species.name[is.infinite(dat_grid$Lifespan.age)])
# 38 species
dat_grid$Lifespan.age[is.infinite(dat_grid$Lifespan.age)] = NA
dat_grid <- dat_grid[rowSums(!is.na(dat_grid))==ncol(dat_grid),]

write.csv(dat_mean,file='./data/Mean_Species_LHT_dataset.csv',row.names=FALSE)
write.csv(dat_grid,file='./data/Species_x_grid_LHT_dataset.csv',row.names=FALSE)

#-- Data Merge 4 | Environmental data -----------------------------------------|
# The following file contains the average environmental data across all land 
# mass within the equal-area hexagon, i.e., grid cell. 
env_dat <- read.csv('./data/20240126_treemort_hex_250k.csv')
head(env_dat)
unique(env_dat$Resolve_Biome)
unique(as.factor(env_dat$Resolve_Biome))
levels(as.factor(env_dat$Resolve_Biome))
env_dat$Resolve_Biome <- round(env_dat$Resolve_Biome ,digit=0)
names(env_dat)[names(env_dat) == "hexgrid"] <- "grid.num"

grid_reduce <- env_dat %>% distinct(grid.num, .keep_all = TRUE) ## remove all the duplicate info, created by having multiple plots per grid cell
grid_reduce <- subset(grid_reduce,!is.na(grid.num))
nrow(grid_reduce)
# write.csv(grid_reduce,file='./data/Treemort_hex_250k.csv',row.names=FALSE)

# Soil [Soil Grid]
soil <- as.data.frame(scale(grid_reduce[,c(44:54)])); names(grid_reduce)[c(44:54)]
soil_cor = cor(soil, use="complete.obs", method = "spearman");round(soil_cor,2)
soil_pca <- princomp(soil) ## Soil organic carbon and bulk density heavy weighted in PC1
soil_pca$loadings
summary(soil_pca) # pc1 only captures 39% of the variance. pcs 1-3 capture 78%
soil_pca <- soil_pca$scores[,1]
# Precipitation [Chelsa +]
precip <- scale(grid_reduce[,c(2,5,12:16,28,42)]); names(grid_reduce)[c(2,5,12:16,28,42)]
precip_cor = cor(precip, use="complete.obs", method = "spearman");round(precip_cor,2)
precip_pca <- princomp(precip)
precip_pca$loadings  ## Annual_Precipitation and Precipitation_of_Warmest_Quarter heavy weighted in PC1, followed by AridityIndex
summary(precip_pca) # pc1 captures 49% of the variance. pcs 1-3 capture 80%
precip_pca <- precip_pca$scores[,1]
# Temperature [Chelsa]
temp <- scale(grid_reduce[,c(4,6:11,18,27,29)]); names(grid_reduce)[c(4,6:11,18,27,29)]
temp_cor = cor(temp, use="complete.obs", method = "spearman");round(temp_cor,2)
temp_pca <- princomp(temp) 
temp_pca$loadings ## CGIAR_Aridity_Index  and Temperature_Seasonality
summary(temp_pca) # pc1 captures 67% of the variance. pcs 1-3 capture 91%
temp_pca <- temp_pca$scores[,1]

#-- Data Merge 5 | Merge Life history with Environmental data -----------------|
# Extract treemort plot id, grid num, and three strongest variables ini soil, precip, temp
env_pcas <- grid_reduce[,c(60,59,5,53,18)]; names(grid_reduce)[c(60,59,5,53,18)] 
NPP <- (grid_reduce[,c(30)]); names(grid_reduce)[30]
Growing_season <- (grid_reduce[,c(27)]); names(grid_reduce)[27]
# all other growing season measures were super correlated, so they were dropped 
# from consideration
env_pcas <- as.data.frame(cbind(env_pcas,NPP,Growing_season,temp_pca,precip_pca,soil_pca))

# Merge life history info with environmental data here, then save.
unique(dat_grid$grid.num[!dat_grid$grid.num %in% env_dat$grid.num]) 
# Check output = numeric(0) indicates a full match, i.e., no grids missing
length(unique(env_dat$grid.num)) ## a total of 134 
dat_grid$grid.num <- as.factor(dat_grid$grid.num)
env_pcas$grid.num <- as.factor(env_pcas$grid.num)

trait_env <- merge(dat_grid, env_pcas, by = c('grid.num'), all.x=TRUE)  ## ok this works. boom
summary(trait_env)
nrow(trait_env);nrow(env_pcas);nrow(dat_grid)
write.csv(trait_env,file='./data/LHT_x_Env_dataset.csv',row.names=FALSE)
#___________________________________________________________________________----