################################################################################
################  4| Compute life history trait clusters   #####################
################################################################################

rm(list=ls(all=TRUE)) ## Clear the workspace
library(factoextra) ## used for cluster analysis
library(car) # For data ellipses
setwd("~/[add wd here]")
dat_grid <- read.csv(file='./Species_x_grid_LHT_dataset.csv')

#_______________________________________________________________________________
#--| loading the life hisotry trait data from IPM script ---------------------->
#_______________________________________________________________________________
dat_grid <- dat_grid[rowSums(!is.na(dat_grid))==ncol(dat_grid),]
dat_grid$Lifespan.size <- (exp(dat_grid$Lifespan.size))*.1 ## convert back from log and change to cm
dat_grid$DBH_70 <- dat_grid$DBH_70*.1
sp <- levels(as.factor(dat_grid$species.name));length(sp)

#_______________________________________________________________________________
#--| PCA of traits and identifying # of trait clusters ------------------------>
#_______________________________________________________________________________
PCA_traits <- as.data.frame(scale(dat_grid[,c(4:14)])) ## doing a PCA for all the traits
LH_corr <- cor(PCA_traits, use="complete.obs", method = "spearman");round(LH_corr,2)
LH_cor_plot <- ggcorrplot(LH_corr, 
                          hc.order = TRUE, 
                          type = "lower",
                          lab = TRUE)

## reducing some of the correlated variables, which capture redundant trait information
## and conducting univariate trait correlations.
PCA_traits2 <- as.data.frame(scale(dat_grid[,c(4,6:7,9,13:14)]))
colnames(PCA_traits2)[1] = "Size at 70%"
colnames(PCA_traits2)[2] = "Passage time,  yrs | 10-20 cm dbh"
colnames(PCA_traits2)[3] = "Passage time, yrs to 70% quantile size"
colnames(PCA_traits2)[4] = "Life exp. from 10 cm dbh"
colnames(PCA_traits2)[5] = "Maximal lifespan age"
colnames(PCA_traits2)[6] = "Maximal size"


## removing highly correlated traits before doing the PCA clustering 
PCA_t <- scale(PCA_traits[,c(3,4,11,6)]) ## PCA of traits
princomp(PCA_t)$loadings 
pcomp = princomp(PCA_t)
data=pcomp$scores

k.max <- 30 ## identifying the clusters of traits
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=25,iter.max = 10)$tot.withinss})

fviz_nbclust(data, kmeans, method = "silhouette")
fviz_nbclust(data, kmeans, method = "wss")

sil <- fviz_nbclust(data, kmeans, method = "silhouette",linecolor = "darkblue") + ## indicates 2 or 4 clusters
  geom_vline(xintercept = 3, linetype = 2,col = "darkblue")+
  labs(title= "Silhouette")
elbow <- fviz_nbclust(data, kmeans, method = "wss",linecolor = "darkblue") +
  geom_vline(xintercept = 4, linetype = 2,col = "darkblue")+
  labs(title = "Elbow method")

## setting the functional trait clusters at 4,
clust = kmeans(data, 4, nstart=25,iter.max = 10) 
clust$cluster <- as.factor(clust$cluster)

#_______________________________________________________________________________
#--| Plotting for Fig 2 ------------------------------------------------------->
#_______________________________________________________________________________
region = c("Tropical moist","Tropical dry","Tropical savanna",
           "Temperate savanna","Temperate broadleaf","Temperate conifer",
           "Boreal")
biome <- c("Tr. moist", "Tr. dry", "Tr savanna","Temp savanna", "Temp broadleaf", "Temp conifer", "Boreal")
biocol <- adjustcolor(hcl.colors(8,palette='Zissou 1'),alpha.f=0.75) ## for biomes
biocol2 <-adjustcolor(c("#000204","grey"),alpha.f=0.6) ## for density plot
biocol3  <- adjustcolor(c("#004A4A", "#DFAF26","#D45E62","#5EA777"),alpha.f=0.8) ## for clusters

sizeaxis = 1.5
cexs= .9
col_line = "grey44"
lwd = 3
cex.leg =.8
sizeaxis =1
tks = -.02
tks.wd=.8

quartz(width=6,height=3);
#x11(width=7,height=9);
par(mfrow=c(1,2),oma=c(.4,.4,0.4,.4),mar=c(4,4,2,2))
grid.arrange(elbow,sil, ncol=2)
#figure_path <- paste0('~/[name path]', "Clusters",".pdf",sep = "")
#dev.copy(pdf,figure_path, width=6,height=3)
#dev.off()

graphics.off(); 
quartz(width=8,height=6.5);
#x11(width=7,height=9);
par(mfrow=c(2,2),oma=c(.4,.4,0.4,.4),mar=c(4,4,2,2))
trop <- as.data.frame(subset(dat_grid, dat_grid$Resolve_Biome %in% c("Tropical moist","Tropical dry","Tropical savanna")))
e_trop <- as.data.frame(subset(dat_grid, dat_grid$Resolve_Biome %in% c("Boreal","Temperate broadleaf","Temperate conifer","Temperate savanna")))
trop$LE_10_log <- log(trop$LE_10)
e_trop$LE_10_log <- log(e_trop$LE_10)
dens_t <- density(log(trop$LE_10))
dens_e <- density(log(e_trop$LE_10))

# plot density, Fig 2 panel A
plot(dens_e, frame = F,xaxt="n",ylim=c(0, .4),ylab= "",xlab= "", col=biocol2[1],main="")
polygon(dens_e, col = biocol2[1])
lines(density(log(e_trop$LE_10)), lwd = 3, col = "grey")
lines(dens_t,xaxt="n",xlim=c(0, 2.5), ylim=c(0, 1.5), ylab= "",xlab= "", col=biocol2[2],main="")
polygon(dens_t, col = biocol2[2])
lines(density(log(trop$LE_10)), lwd = 3, col = "grey")

v1 <- summary(trop$LE_10_log)[3]
abline(v=v1,col=biocol2[2],lty=3, lwd=3)
exp(summary(trop$LE_10_log)[4])
text(x = v1-1, y = .05, label = "60 yrs",col="black")

v2 <- summary(e_trop$LE_10_log)[4]
abline(v=v2,col=biocol2[1],lty=3, lwd=3)
exp(summary(e_trop$LE_10_log)[4])
text(x = v2+1, y = .05, label = "95 yrs",col="black")

legend("topleft", c("Tropics","Extratropics"),
       lty=3, col = c(biocol2[2],biocol2[1]),bty='n',cex=cex.leg,lwd=3)

axis(1, at = c(log(1),log(10),log(100),log(1000),log(2000)),labels=c("1","10","100","1000","2000"))
axis(1, at = c(log(1),log(2),log(3),log(4),log(5),
               log(10),log(20),log(30),log(40),log(50),log(60),
               log(70),log(80),log(90),
               log(100),log(200),log(300),log(400),log(500),
               log(600),log(700),log(800),log(900),
               log(1000),log(2000)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
mtext('Density',side=2,line=2.5,outer=FALSE,cex=cexs)
mtext('Life exp from 10cm dbh (yrs)',side=1,line=2.5,outer=FALSE,cex=cexs)
mtext('A',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)

## Fig 2 panel B
plot(dat_grid$LE_10~dat_grid$tau_20,col=biocol[dat_grid$Biome_num], frame=F,
     pch=19, ylab = "", xlab = "",log='yx',yaxt="n",xaxt="n",ylim=c(1,6000),xlim=c(2,170))
axis(2, at = c(1,10,100,1000),labels=TRUE)
axis(2, at = c(0,1),xaxt = "n",labels=TRUE)
axis(2, at = c(1:5,seq(10, 100,10),seq(100,1000,100),seq(1000, 3000,1000)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
axis(1, at = c(10, 100),labels=TRUE)
axis(1, at = c(0,1),xaxt = "n",labels=TRUE)
axis(1, at = c(1:5,seq(10,200,10)),labels=FALSE,lwd.ticks=tks.wd,tck=tks)
mtext('Years to grow from 10-20cm dbh',side=1,line=2.5,outer=FALSE,cex=cexs)
mtext('Life exp from 10cm dbh (yrs)',side=2,line=2.5,outer=F,cex=cexs)
legend('topleft',biome,col=biocol,pch=19,bty='n',cex=cex.leg)
mtext('B',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)

## Fig 2 panel C
plot(dat_pca$Comp.2~dat_pca$Comp.1,col=biocol3[dat_pca$cluster],pch=16, ylab= "", yaxt="n",xaxt="n",
     xlim = c(-3, 6),ylim = c(-8, 4),xlab="",frame=F)#,yaxt="n",xaxt="n") ## plotting 1st pca
axis(1, at = c(-2,2,6),labels=TRUE)
axis(1, at = c(-7,-2,2),yaxt = "n",labels=TRUE)
axis(2, at = c(-6,-2,2),labels=TRUE)
axis(2, at = c(-8,8),xaxt = "n",labels=FALSE)

points(clust$centers[,1],clust$centers[,2],col=biocol3,pch=18,cex=2)
Trait <- c("Cluster 1", "Cluster 2", "Cluster 3","Cluster 4")
group_col <- c("#D45E62","#004A4A","#DFAF26","#5EA777")
box(bty='l')
legend('bottomleft',Trait,col=group_col,pch=19,bty='n',cex=cex.leg)
abline(v=0,col='grey44',lty=3, lwd=1.5)
abline(h=0,col='grey44',lty=3, lwd=1.5)
summary(pcomp) ## variation explained by each PCA axis
mtext('PC axis 2 (28%), max. size (small-large) ',side=2,line=2.5,outer=FALSE,cex=cexs)
mtext('PC axis 1 (46%), growth strategy (fast-slow)',side=1,line=2.5,outer=FALSE,cex=cexs)
mtext('C',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
for(i in 1:1067) lines(c(dat_pca[i,1],clust$centers[dat_pca$cluster[i],1]),
                       c(dat_pca[i,2],clust$centers[dat_pca$cluster[i],2]),
                       col=biocol3[dat_pca$cluster[i]])
## Fig 2 panel D
plot(dat_pca$Comp.3~dat_pca$Comp.2,col=biocol3[dat_pca$cluster],pch=16, yaxt="n",xaxt="n",
     xlim = c(-8, 3),ylim = c(-4, 12),ylab= "", xlab="",frame=F)#,yaxt="n",xaxt="n") ## plotting 1st pca
points(clust$centers[,2],clust$centers[,3],col=biocol3,pch=18,cex=2)
axis(1, at = c(-6,-2,2),labels=TRUE)
axis(2, at = c(-3,4.5,12),labels=TRUE)
axis(2, at = c(-5,12),xaxt = "n",labels=FALSE)
axis(1, at = c(-8,4),yaxt = "n",labels=FALSE)
mtext('PC axis 3 (21%), life exp (low-high)',side=2,line=2.5,outer=FALSE,cex=cexs)
mtext('PC axis 2 (28%), max. size (small-large)',side=1,line=2.5,outer=FALSE,cex=cexs)
for(i in 1:1067) lines(c(dat_pca[i,2],clust$centers[dat_pca$cluster[i],2]),
                       c(dat_pca[i,3],clust$centers[dat_pca$cluster[i],3]),
                       col=biocol3[dat_pca$cluster[i]])
mtext('D',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
abline(v=0,col='grey44',lty=3, lwd=1.5)
abline(h=0,col='grey44',lty=3, lwd=1.5)
box(bty='l')

#figure_path <- paste0('~/[name path]/', "Figure 2",".pdf",sep = "")
#dev.copy(pdf,figure_path,width=8,height=6.5)
#dev.off()

#_______________________________________________________________________________
#--| plot for Fig S8 |--------------------------------------------------------->
#_______________________________________________________________________________
graphics.off(); 
quartz(width=8,height=4);
#x11(width=7,height=9);
par(mfrow=c(1,2),oma=c(.4,.4,0.4,.4),mar=c(4,4,2,2))
ran.size <- range((dat_grid$Lifespan.size))*c(0.9,1.1)
ran.tau <- range((dat_grid$tau_20))*c(0.90,1.1)
plot(dat_grid$Lifespan.size~dat_grid$tau_20,col=biocol[dat_grid$Biome_num], pch=19, ylab = "", xlab = "",log='yx',ylim=ran.size,xlim=ran.tau,yaxt="n",xaxt="n")
axis(2, at = c(1,10,100,1000),labels=TRUE)
axis(2, at = c(1:5,seq(10, 100,10),seq(100,1000,100),seq(1000, 3000,1000)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
axis(1, at = c(1,10,100,1000),labels=TRUE)
axis(1, at = c(1:5,seq(10, 100,10),seq(100,1000,100),seq(1000, 3000,1000)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
mtext('Years to 20 cm dbh',side=1,line=2.5,outer=FALSE,cex=cexs)
mtext('Max size (cm dbh)',side=2,line=2.5,outer=F,cex=cexs)
mtext('A',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
legend('topleft',biome,col=biocol,pch=19,bty='n',cex=cex.leg)

ran.le <- range((dat_grid$LE_10))*c(0.90,1.1)
plot(dat_grid$LE_10~dat_grid$Lifespan.size,col=biocol[dat_grid$Biome_num], pch=19, ylab = "", xlab = "",log='yx',xlim=ran.size,ylim=ran.le,yaxt="n",xaxt="n")
axis(2, at = c(1,10,100,1000),labels=TRUE)
axis(2, at = c(1:5,seq(10, 100,10),seq(100,1000,100),seq(1000, 3000,1000)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
axis(1, at = c(1,10,100,1000),labels=TRUE)
axis(1, at = c(1:5,seq(10, 100,10),seq(100,1000,100),seq(1000, 3000,1000)),labels=FALSE, lwd.ticks=tks.wd,tck=tks)
mtext('Max size (cm dbh)',side=1,line=2.5,outer=FALSE,cex=cexs)
mtext('Life exp from 10 cm dbh (yrs)',side=2,line=2.5,outer=F,cex=cexs)
mtext('B',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)

#figure_path <- paste0('~/[name path]/', "S1_trait_relationships",".pdf",sep = "")
#dev.copy(pdf,figure_path,width=8,height=4)
#dev.off()

