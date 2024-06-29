################################################################################
#################  6| Calculate and analyse convex hulls   #####################
################################################################################
rm(list=ls(all=TRUE)) ## Clear the workspace
setwd("~/[name]")
library(car) # For data ellipses
library(tidyverse)
library(geometry)
library(corrplot)
library(margins)
library(broom)

dat_grid <- read.csv(file='./data/Species_x_grid_LHT_dataset.csv')
dat_grid <- dat_grid[rowSums(!is.na(dat_grid))==ncol(dat_grid),]
dat_grid$Lifespan.size <- (exp(dat_grid$Lifespan.size))*.1 ## convert back from log and change to cm
dat_grid$DBH_70 <- dat_grid$DBH_70*.1

PCA_t <- scale(dat_grid[,c(6:7,9,14)])
pcomp = princomp(PCA_t)
var.ex = pcomp$sdev^2/sum(pcomp$sdev^2)
data=as.data.frame(pcomp$scores)
data$Comp.1 <- data$Comp.1*var.ex[1] ## scale by proportion of variance explained by PCA axis
data$Comp.2 <- data$Comp.2*var.ex[2] 
data$Comp.3 <- data$Comp.3*var.ex[3] 

#___________________________________________________________________________----
#--| calculating the demographic trait diversity, using PCA 1, 2, and 3 |-->----
hull <- as.data.frame(cbind(dat_grid[1],dat_grid[2],dat_grid[16],data)) ## First add biome and grid number
hull$grid.num <- as.factor(hull$grid.num)
# the convex hull calculation requires one more data point than the number of dimensions
# If we do a 3-d hull, then we need four or more data per grid cell.
these.gs = levels(hull$grid.num)[which(tapply(hull$grid.num,hull$grid.num,length)<=4)]
hull <- droplevels(subset(hull,!grid.num %in% these.gs))  
summary(hull$grid.num,options(max.print=100000))

##this code calculates the volume of the convex hull for each gridcell based on the Comp.1, Comp.2, and Comp.3 coordinates.
grids = levels(hull$grid.num)
volume = matrix(NA,nr=length(grids),nc=1)
for (i in 1:length(grids)) {
  volume[i] = convhulln(data.matrix(hull[hull$grid.num==grids[i] , 
                                        c("Comp.1","Comp.2",'Comp.3')]), option="FA")$vol
}
grid.num <- levels(as.factor(hull$grid.num))

##sp_richness, species richness at the grid cell level
sp_richness <- hull %>%
  group_by(grid.num) %>%
  summarize(species.name = n_distinct(species.name),.groups = 'drop')

dat <- as.data.frame(cbind(grid.num,volume,sp_richness[2]))
colnames(dat)[3] <- "sp.richness"

grid_reduce <- hull %>% distinct(grid.num, .keep_all = TRUE) ## reduces dataset for one row per grid.num
dat <- merge(dat, grid_reduce, by = c('grid.num'), all.x=TRUE)  ## 
dat$volume <- as.numeric(dat$volume)

env_data = read.csv(file='./data/hex_250k.csv')
grid_ID <- levels(as.factor(grid_reduce$grid.num));length(grid_ID) 
dat2 <- merge(dat, env_data , by = c('grid.num'))
names(dat2)[names(dat2) == "Resolve_Biome.x"] <- "biome"
names(dat2)[names(dat2) == "CHELSA_exBIO_NPP"] <- "NPP" ## grams carbon per meter^2  *10
names(dat2)[names(dat2) == "CHELSA_BIO_Annual_Mean_Temperature"] <- "temp"

#___________________________________________________________________________----
#--| Demographic diversity patterns, using trait PCA 1, 2, and 3 |--------->----

cor(dat2$sp.richness,dat2$volume) ## species diversity and species richness are correlated, with diversity saturating at higher levels of diversity. 
cor(log(dat2$NPP),log(dat2$volume)) ## corr of .70
richness_mod <- lm(log(dat2$volume) ~ log(dat2$sp.richness) + I(log(dat2$sp.richness)^2)); summary(richness_mod)
#write.csv( glance(richness_mod) , "./output/spFig4.1_anova.csv" )
#write.csv( tidy(richness_mod) , "./output/spFig4.1_coeff.csv" )

cor(dat2$temp,dat2$NPP) ## super correlated variables (corr = 0.94). It doesn't make sense to include both in a glm model
temp_mod <- lm(formula =(log(volume))~(temp),data=dat2); summary(temp_mod) ## diversity varied predictable across a temp gradient
summary(temp_mod)$adj.r.squared ##R-square = 0.40
#write.csv( glance(temp_mod) , "./output/tempFig4.2_anova.csv" )
#write.csv( tidy(temp_mod) , "./output/tempFig4.2_coeff.csv" )

npp_mod <- lm(formula =scale(log(volume))~scale(log(NPP))-1,data=dat2); summary(npp_mod)##R-square = 0.49
#write.csv( glance(npp_mod) , "./output/NPPFig4.3_anova.csv" )
#write.csv( tidy(npp_mod) , "./output/NPPFig4.3_coeff.csv" )

volume_mod <- lm(formula = scale(log(NPP))~scale(log(volume))+scale((temp))-1,data=dat2);summary(volume_mod)
summary(margins(volume_mod, variables = 'volume'))
summary(margins(volume_mod, variables = 'temp')) ## both richness and diversity are predictive of productivity. 
cplot(volume_mod, "volume")
summary(volume_mod)$adj.r.squared ##0.90
#write.csv( glance(volume_mod) , "./output/DDFig4.4_anova.csv" )
#write.csv( tidy(volume_mod) , "./output/DDFig4.4_coeff.csv" )

#___________________________________________________________________________----
#--| Plotting results, Fig. 4 in main text |------------------------------->----
sizeaxis = 1
cexs= .9
col_line = "grey44"
lwd = 4
cex.leg =.9
lty = 4
al =2.5 ## x and y labels
pch=19
cex.point=1.3
ran.temp <- range(dat2$temp)*c(1.07,1.07)
ran.vol <- range(log(dat2$volume))*c(1.07,1.07)
ran.npp <- range(scale(dat2$NPP))*c(0.93,1.07)

region = c("Tropical moist","Tropical dry","Tropical savanna",
           "Temperate savanna","Temperate broadleaf","Temperate conifer",
           "Boreal")
dat2$biome <- ordered(dat2$biome, levels = c("Tropical moist","Tropical dry","Tropical savanna",
                                                       "Temperate savanna","Temperate broadleaf","Temperate conifer",
                                                       "Boreal"))
dat2$biome <- as.factor(dat2$biome)
## plots of the main figures presented in Fig 4
quartz(width=8,height=6.5);
#x11(width=7,height=9);
biocol <- hcl.colors(8,palette='Zissou 1')
par(mfrow=c(2,2),oma=c(.4,.4,0.1,.1),mar=c(4,4,2,2))
ran.vol <- range(log(dat2$volume))*c(1.07,1.07)
ran.sp <- range(log(dat2$sp.richness))*c(0.93,1.07)
plot(log(dat2$volume)~log(dat2$sp.richness),col=biocol[dat2$biome],pch=pch,frame=F,
     ylab = "", xlab = "",las=2, xaxt="n",yaxt="n",ylim=ran.vol,xlim=ran.sp,cex=cex.point)
box(bty='l')
axis(1, at = c(1,2,4,6),labels=TRUE)
#axis(1, at = c(1,2,4,6),labels=FALSE)
axis(2,pretty(ran.vol,n=4))
mtext('A',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
mtext('Demographic trait diversity (log)',side=2,line=al,outer=F,cex=sizeaxis)
mtext('Species richness (log)',side=1,line=al,outer=F,cex=sizeaxis)
legend('bottomright',region,col=biocol,pch=19,bty='n',cex=cex.leg)
summary(richness_mod)
#use model to get predicted values
pred <- predict(richness_mod)
ix <- sort(dat2$sp.richness, index.return=T)$ix
#add polynomial curve to plot
lines(log(dat2$sp.richness)[ix], pred[ix], col='darkgrey', lty=lty,lwd=lwd)
mtext(expression("p < 0.001, adj R" ^ "2"*" = .65"),side=3,line=-1.5,outer=FALSE,cex=.8,adj=0.05)


plot(log(dat2$volume)~dat2$temp,col=biocol[dat2$biome],pch=pch, ylab = "", xlab = "", 
     frame=F,ylim=ran.vol,xlim=ran.temp,yaxt="n",xaxt="n",cex=cex.point)
mtext('Demographic trait diversity (log)',side=2,line=al,outer=F,cex=sizeaxis)
mtext(expression(paste("Mean annual temp. (",degree~C,")")),side=1,line=al,outer=F,cex=sizeaxis)
mtext('B',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
summary(temp_mod)
mtext(expression("p < 0.001, adj R" ^ "2"*" = .40"),side=3,line=-1.5,outer=FALSE,cex=.8,adj=0.05)
abline(lm(log(volume)~temp, data = dat2),lty=lty,lwd=lwd,col = "darkgrey")
box(bty='l')
axis(2,pretty(ran.vol,n=4))
axis(1,pretty(ran.temp,n=4))

ran.npp <- range(scale(dat2$NPP))*c(1.5,1.07)
summary(npp_mod)
plot(log(dat2$volume)~scale(log(dat2$NPP)),col=biocol[dat2$biome], pch=pch, ylab = "", xlab = "", 
     frame=F,xlim=ran.npp,ylim=ran.vol,yaxt="n",xaxt="n",cex=cex.point)
axis(2,pretty(ran.vol,n=4))
axis(1,pretty(ran.npp,n=4))
mtext(expression(paste("NPP, g m"^"-2", " yr"^"-1", " (log)")),side=1,line=al,outer=F,cex=sizeaxis)
mtext('Demographic trait diversity (log)',side=2,line=al,outer=F,cex=sizeaxis)
mtext('C',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
mtext(expression("p < 0.001, adj R" ^ "2"*" = .49"),side=3,line=-1.5,outer=FALSE,cex=.8,adj=0.05)
abline(lm(log(volume)~scale(NPP), data = dat2),lty=lty,lwd=lwd,col = "darkgrey")
box(bty='l')

summary(volume_mod)
summary(margins(volume_mod, variables = 'volume'))
summary(margins(volume_mod, variables = 'temp')) ## both richness and diversity are predictive of productivity. 
plot(scale(log(dat2$NPP))~(log(dat2$volume)),col=biocol[dat2$biome],pch=pch, ylab = "", xlab = "", 
     frame=F,xlim=ran.vol,ylim=ran.npp,cex=cex.point,yaxt="n",xaxt="n")
axis(2,pretty(ran.npp,n=4))
axis(1,pretty(ran.vol,n=4))
mtext('D',side=3,line=0.4,outer=FALSE,cex=sizeaxis,adj=0)
mtext('Demographic trait diversity (log)',side=1,line=al,outer=F,cex=sizeaxis)
mtext(expression(paste("NPP, g m"^"-2", " yr"^"-1"," (log)")),side=2,line=al,outer=F,cex=sizeaxis)
mtext(expression('Demo. diversity, AME = 1.43, p < 0.01'),side=3,line=-1.5,outer=FALSE,cex=.8,adj=0.02)
abline(lm(scale(log(dat2$NPP))~(log(dat2$volume)), data = dat2), lty=lty,lwd=lwd,col = "darkgrey")
box(bty='l')

#figure_path <- paste0('~/[path]/', "Hull volume",".pdf",sep = "")
#dev.copy(pdf,figure_path, width=8,height=6.5)
#dev.off()

# ggsave("plot.pdf", width = 15,  height = 10,  units = "cm",dpi = 300)
# png("test.png",width=3.25,height=3.25,units="in",res=1200)
# par(mar=c(5,5,2,2),xaxs = "i",yaxs = "i",cex.axis=1.3,cex.lab=1.4)
# plot(perf,avg="vertical",spread.estimate="stddev",col="black",lty=3, lwd=3)
# dev.off()
## science format guidelines: 300 dpi, 5.7 cm (2.24 inches or 1 column) or 12.1 cm (4.76 inches or 2 columns), or 18.4 cm (7.24 inches or 3 columns)

#___________________________________________________________________________----
#--| Tropical versus extra-tropical trends |------------------------------->----
trop <- as.data.frame(subset(dat2, dat2$biome %in% c("Tropical moist","Tropical dry","Tropical savanna")))
e_trop <- as.data.frame(subset(dat2, dat2$biome %in% c("Boreal","Temperate broadleaf","Temperate conifer","Temperate savanna")))

temp_mod_trop = lm(scale(log(NPP))~scale((temp))+scale(log(volume))-1,data=trop)
temp_mod_trop2 = lm(scale(log(NPP))~scale((temp))-1,data=trop)
anova(temp_mod_trop,temp_mod_trop2)
summary(temp_mod_trop)
write.csv( glance(mt) , "./output/temp_mod_trop_anova.csv" )
write.csv( tidy(mt) , "./output/temp_mod_trop_coeff.csv" )

temp_mod_etrop = lm(scale(log(NPP))~scale((temp))+scale(log(volume))-1,data=e_trop)
temp_mod_etrop2 = lm(scale(log(NPP))~scale((temp))-1,data=e_trop)
anova(temp_mod_etrop,temp_mod_etrop2)
write.csv( glance(me) , "./output/NPP_etropics_anova.csv" )
write.csv( tidy(me) , "./output/NPP_etropics_coeff.csv" )

## tropics and extra-tropics analysis, NPP predicting volume
npp_trop = lm(log(volume)~scale(log(NPP)),data=trop); summary(npp_trop)
npp_etrop = lm(log(volume)~scale(log(NPP)),data=e_trop); summary(npp_etrop)


