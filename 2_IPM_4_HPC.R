################################################################################
##########  2| IPM x Species x Grid | Developed to run on HPC  #################
################################################################################
#--| Housekeeping |-------------------------------------------------------->----
rm(list=ls(all=TRUE)) ## Clear the workspace
# setwd("C:/Users/mceld/Dropbox/Tree_Longevity_paper/Dryad_files") # RMc computer
# The following packages will be called. Ensure that the slurm job script calls 
# these pacakages on the selected nodel.
# library(statmod) --> locally loaded
# library(dplyr) --> available
# library(data.table) --> available
library("statmod",lib.loc="./R/x86_64-pc-linux-gnu-library/4.1/") # load library from local folder

# load data and the model parameters.
load('./Rdata/TreeMort_Survival_1yr_woGrid.Rdata')
load('./Rdata/TreeMort_Growth_1yr_woGrid.Rdata')
load('./Rdata/TreeMort_sigmag_1yr_woGrid.Rdata')
load('./Rdata/Min_Max_Size.Rdata')

i.end = c(seq(100,1100,100),length(spp))
i.start = c(1,i.end[-12]+1)
index = as.numeric(commandArgs(trailingOnly = TRUE)) # SLURM array interation saved here
print(paste0('Index # ', index,', i = ',i.start[index],':',i.end[index]))

print(paste('Check sum zero =',sum(n.grid!=rowSums(!is.na(fitg_Bg))))) # Check sum = 0
print(paste('Check sum zero =',sum(n.grid!=rowSums(!is.na(fits_Bg))))) # Check sum = 0
print(paste('Check sum zero =',sum(is.na(fitg_B)))) # Check sum = 0
print(paste('Check sum zero =',sum(is.na(fits_B)))) # Check sum = 0
print(paste('Check sum zero =',sum(is.na(fitg_sigma)))) # Check sum = 0
#___________________________________________________________________________----
#--| A. sx and gyx FUNCTIONs |--------------------------------------------->----
sx<-function(x,B_s,mu.x){
  nc = length(B_s)
  cen.x = x - mu.x
  X = matrix(NA,nrow=length(x),ncol=nc)
  for(k in 1:nc) X[,k] = cen.x^(k-1)
  p = X%*%B_s
  s<-exp(-exp(p) ) ## for a 1-yr timestep. log(-log(1-p)) instead of logit log(p/(1-p))
  #s<-exp(xbeta)/( 1+exp(xbeta) ) ## prob of survival
  return(s);
}

## B. GROWTH FUNCTION g(y,x)
gyx<-function(y,x,B_g,sig_g){
  mux<-B_g[1] + B_g[2]*x ## mean
  fac1<-sqrt(2*pi)*sig_g;        #use sigma from growth model
  fac2<-((y-mux)^2)/(2*sig_g^2); 
  # return(exp(-fac2)/fac1);
  return(dnorm(y,mean=mux,sd=sig_g))
}
g_mean<-function(x,B_g){
  mux<-B_g[1] + B_g[2]*x ## mean
  return(mux)
}

## C. SURVIVAL-GROWTH p(y,x) = s(x)*g(y,x)
pyx<-function(y,x,B_s,mu.x,B_g,sig_g){
  p<-sx(x,B_s,mu.x)*gyx(y,x,B_g,sig_g)
  return(p)
}

# Functions to compute the survival-growth kernel (P) and the fundamental kernel (N)
require(statmod); 
# Gauss-Legendre quadrature on interval (L,U) 
gaussQuadInt <- function(L,U,order) {
  # nodes and weights on [-1,1]
  out <- gauss.quad(order); #GL is the default 
  w <- out$weights; x <- out$nodes;  
  weights=0.5*(U-L)*w; 
  nodes=0.5*(U+L) + 0.5*(U-L)*x; 
  return(list(weights=weights,nodes=nodes)); 
}

Pmat<- function(m,L,U,B_g,sig_g,B_s,mu.x,level = 420,diag_tresh= 50, 
                correction = "ceiling") {
  h <- (U - L) / m
  meshpts <- L + ((1:m) - 1/2) * h
  df <- expand.grid(idx1=1:m, idx=1:m)
  if(!correction %in% c("ceiling", "sizeExtremes", "constant", "none"))
    stop("correction should be in ceiling, sizeExtremes, constant, or none")
  require(dplyr)
  df <- df %>% mutate(id = 1:n(),
                      G = 0,
                      x1 = meshpts[idx1],
                      x = meshpts[idx],
                      pos_sel = x1>=x,
                      diag_dist = sqrt(((x+x1)/2-x)^2+(x1-(x+x1)/2)^2),
                      diag_sel = diag_dist < diag_tresh & pos_sel,
                      up_diag_sel_mb = pos_sel & diag_dist >= diag_tresh)
  require(data.table)
  out1 <- gaussQuadInt(-h/2, h/2, 3) # This is where we choose different number of nodes
  out2 <- gaussQuadInt(-h/2, h/2, floor(level/3))
  # create data for mean growth along x for the integration on dim 3
  quadx <- expand.grid(idx=1:m,
                       map1=seq.int(length(out1$weights)))
  quadx <- quadx %>% mutate(x = meshpts[idx] + out1$nodes[map1],
                            idx_3 = paste0(idx, c("A", "B", "C")[map1]))
  quadx <- as.data.table(quadx)
  quadx <- quadx[ , g_m:= g_mean(x, B_g)]
  g_mt <- quadx$g_m
  names(g_mt) <- quadx$idx_3
  # create data for integration on 3 x floor(level)
  quad <- expand.grid(id=df$id[df$diag_sel],
                      map1=seq.int(length(out1$weights)),
                      map2=seq.int(length(out2$weights)))
  # geat g_mean for all x values
  quad <- quad %>% mutate(x = df$x[id] + out1$nodes[map1],
                          x1 = df$x1[id] + out2$nodes[map2],
                          idx = df$idx[id],
                          idx1 = df$idx1[id],
                          idx_3 = paste0(idx, c("A", "B", "C")[map1]),
                          weights1 = out1$weights[map1],
                          weights2 = out2$weights[map2],
                          g_m= g_mt[idx_3])
  quad <- as.data.table(quad)
  # compute growth kernel
  quad <- quad[ , fvals := gyx(x1,x,B_g,sig_g)]
  quad <- quad[ , G:= fvals *weights1*weights2]
  quad <- quad[, Gs:=sum(G), by = id]
  uniquad <- subset(unique(quad,by = "id"))
  rm(quad)
  res_GLe<- uniquad$Gs/h
  rm(uniquad)
  # mid-bin integration for the rest
  df2 <- df %>% filter(up_diag_sel_mb)
  df2 <- as.data.table(df2)
  df2 <- df2[ , G := gyx(x1,x,B_g,sig_g)*h]
  df$G[df$diag_sel] <- res_GLe
  df$G[df$up_diag_sel_mb] <- df2$G
  rm(df2, res_GLe)
  G <- df$G
  dim(G) <- c(m,m)
  rm(df)
  # eviction correction as in IPMpack
  if (correction == "constant"){# Based on IPMpack
    nvals <- colSums(G)
    P <- t((t(G)/nvals) * sx(x=meshpts,B_s,mu.x))
  }
  if (correction == "sizeExtremes"){# Based on IPMpack
    select_size_t <- meshpts > (U- 2*diag_tresh)
    DiffNvals <- pmax(1- colSums(G), 0)
    G[m, select_size_t] <- G[m, select_size_t] + DiffNvals[select_size_t]
    P <- G * sx(x=meshpts,B_s,mu.x)
  }
  if (correction == "ceiling"){# Based on Williams et al. 2012 Ecology
    DiffNvals <- pmax(1- colSums(G), 0)
    s_x <- sx(x=c(meshpts,U),B_s,mu.x)
    G <- cbind(G, rep(0, length.out = length(meshpts)))
    G <- rbind(G, c(DiffNvals,1))
    # P<- t(t(G) * s_x) # This seems to no longer work
    P<- G
    for(k in (0:m)+1) P[k,]<- G[k,] * s_x
    meshpts = c(meshpts,U)
  }
  if (correction == "none"){
    P <- t(t(G) * sx(x=meshpts,B_s,mu.x))
  }
  return(list(meshpts = meshpts, P = P))
}

Nmat<-function(P,y0,yf,meshpts) { # y0 (initial size) and yf (target size)
  dimP = nrow(P)
  I = diag(1,dimP)
  y=meshpts;
  yf.i = which(abs(y - yf)==min(abs(y - yf)))
  y0.i = which(abs(y - y0)==min(abs(y - y0)))
  if(length(yf.i)>1) yf.i = max(yf.i)
  if(length(y0.i)>1) y0.i = max(y0.i)
  N = solve(I-P)
  P_prime = P
  P_prime[,yf.i] <- 0
  N_prime <- solve(I-P_prime)
  tau <- (N_prime%*%N_prime)[yf.i,y0.i] / N_prime[yf.i,y0.i]
  LE.0 <- sum(N[,y0.i])
  LE.f <- tau + sum(N[,yf.i])
  return(list('N'=N,'tau'=tau,'LE.0'=LE.0,'LE.f'=LE.f))
}

# Lifespan [age at which cumulative survivorship is 95%]
# Age from stage, quasi stationary dist, and age-dep mortality following 
# Horvitz and Tuljapurkar 2008. Start with a cohort of newborns (age=0) at t=0. 
# We are interested in the proportion of the initial cohort surviving to age x.
# We normalize the initial cohort to sum to one

Lspan = function(P,survivorship,meshpts) {
  ntime = 1e2 # Time horizon to capture 95% cohort mortality.
  inc = 1e2 # incremental increase until mortality horizon reached.
  m = nrow(P)
  n <- matrix(0,nrow=m,ncol=ntime)
  n[1,1] <- 1
  # Use the recursion equation  n(t+1) = P*n(t)
  # Time here equals age
  for (t in 2:ntime) {
    n[,t] <- P%*%n[,t-1]
  }
  while (sum(n[,ntime])>survivorship) { # Set limit of 10,000 yrs
    if(ntime>1e4) {print("run time error: time horizon too big"); break }
    ntime = ntime + inc # add 100 additional years if needed
    n <- cbind(n,matrix(0,nrow=m,ncol=inc))
    # Continue the recursion equation  n(t+1) = P*n(t)
    for (t in (ntime-inc):ntime) {
      n[,t] <- P%*%n[,t-1]
    }
  }
  # Survivorship to Age x (lx) = sum(Ni(x))
  lx <- colSums(n)
  age = min(which(lx<=survivorship))
  dbh = sum(meshpts*n[,age]/sum(n[,age]))
  return(list('age'=age,'dbh'=dbh))
}
#___________________________________________________________________________----
#--| Calculate Life History |---------------------------------------------->----
# Compute the survival-growth kernel (P) and the fundamental kernel (N)
dimP <- 600 # Ceiling eviction correction adds an extra row and column
P_spp <- N_spp <- array(NA,dim=c(length(spp),dimP+1,dimP+1))
tau <- LE.0 <- LE.f <- matrix(NA,nrow=length(spp),ncol=3)
dimnames(tau) <- list(spp,c('20cm','70th','90th'))
dimnames(LE.f) <- dimnames(LE.0) <- dimnames(tau)
tau.g <- LE.0.g <- LE.f.g <- array(NA,dim=c(length(spp),3,max(n.grid)))
dimnames(tau.g) <- list(spp,c('20cm','70th','90th'),paste0('g',1:max(n.grid)))
dimnames(LE.f.g) <- dimnames(LE.0.g) <- dimnames(tau.g)
lifespan <- matrix(NA,nrow=length(spp),ncol=2)
dimnames(lifespan) <- list(spp,c('age','dbh'))
lifespan.g <- array(NA,dim=c(length(spp),2,max(n.grid)))
dimnames(lifespan.g) <- list(spp,c('age','dbh'),paste0('g',1:max(n.grid)))
# Set initial and target sizes for passage time and life expectancy.
y0 = log(100) # initial size is currently fixed at 10 cm
yf = log(200) # This is our target size
#--| Begin Loop |--------------------------------------------------------------|
start<-as.POSIXlt(Sys.time()); 
print(paste('Begin loop at:',start))
# for (i in 1:length(spp)) {
for (i in i.start[index]:i.end[index]) { 
  # Calculate prediction surfaces
  Pout = Pmat(m=dimP,L=min.sz[i],U=max.sz[i],fitg_B[i,],fitg_sigma[i],
              fits_B[i,],mu.size[i],level=420,diag_tresh=50, correction="ceiling")
  if(sum(is.na(Pout$P))>=1) 
    print(paste('Problem with',spp[i],', index',i,', mean matrix'))
  P_spp[i,,] <- Pout$P
  Nlist <- Nmat(Pout$P,y0,yf,Pout$meshpts)
  N_spp[i,,] <- Nlist$N
  tau[i,1] <- Nlist$tau
  LE.0[i,1] <- Nlist$LE.0
  LE.f[i,1] <- Nlist$LE.f
  lx = Lspan(Pout$P,0.05,Pout$meshpts) # Set at 0.05 survivorship
  lifespan[i,1] = lx$age
  lifespan[i,2] = lx$dbh
  rm(Nlist,lx); gc()
  for(j in 1:2) {
    # Calculate prediction surfaces
    Nlist <- Nmat(Pout$P,y0,log(yf.spp[i,j]),Pout$meshpts)
    tau[i,j+1] <- Nlist$tau
    LE.0[i,j+1] <- Nlist$LE.0
    LE.f[i,j+1] <- Nlist$LE.f
    rm(Nlist)
  }
  rm(Pout); gc()
  if(n.grid[i]>1) { # run same sequence for each grid cell
    for(g in 1:n.grid[i]) {
      # Calculate prediction surfaces
      Pout = Pmat(m=dimP,L=min.sz[i],U=max.sz[i],c(fitg_Bg[i,g],fitg_B[i,2]),
                  fitg_sigma[i],c(fits_Bg[i,g],fits_B[i,2]),mu.size[i],level=420,
                  diag_tresh=50, correction="ceiling")
      if(sum(is.na(Pout$P))>=1) 
        print(paste('Problem with',spp[i],', index',i,', grid number',g))
      Nlist <- Nmat(Pout$P,y0,yf,Pout$meshpts)
      tau.g[i,1,g] <- Nlist$tau
      LE.0.g[i,1,g] <- Nlist$LE.0
      LE.f.g[i,1,g] <- Nlist$LE.f
      lx = Lspan(Pout$P,0.05,Pout$meshpts) # Set at 0.05 survivorship
      lifespan.g[i,1,g] = lx$age
      lifespan.g[i,2,g] = lx$dbh
      rm(Nlist); gc()
      for(j in 1:2) {
        # Calculate prediction surfaces
        Nlist <- Nmat(Pout$P,y0,log(yf.spp[i,j]),Pout$meshpts)
        tau.g[i,j+1,g] <- Nlist$tau
        LE.0.g[i,j+1,g] <- Nlist$LE.0
        LE.f.g[i,j+1,g] <- Nlist$LE.f
        rm(Nlist); gc()
      }
      rm(Pout); gc()
    }
  }
  if(n.grid[i]==1) { # run once for the lone grid cell :(
    # Fill in with "mean" for the n of 1's 
    tau.g[i,1:3,1] <- tau[i,1:3]
    LE.0.g[i,1:3,1] <- LE.0[i,1:3]
    LE.f.g[i,1:3,1] <- LE.f[i,1:3]
    lifespan.g[i,1,1] = lifespan[i,1]
    lifespan.g[i,2,1] = lifespan[i,2]
  }
  # Run Time printout
  if (i %% 10 == 0) { 
    print(paste('Index',i,',',i-i.start[index]+1,'% Complete'))
    elapsedtime <- as.POSIXlt(Sys.time())-start;
    print(paste("Elapsed Time =",round(elapsedtime,2),
                attr(elapsedtime,which='units')))
  }
}
#--| End Loop |----------------------------------------------------------------|
elapsedtime <- as.POSIXlt(Sys.time())-start;
print(paste('Array #',index,'100% Complete'))
print(paste("Total Run Time =",round(elapsedtime,2),
            attr(elapsedtime,which='units')))
save(P_spp,N_spp,file=paste0('./Rdata/TreeMort_SurvGrow_Kernels',index,'.Rdata'))
print(paste0('saved ./Rdata/TreeMort_SurvGrow_Kernels',index,'.Rdata'))
save(tau,LE.f,LE.0,lifespan,tau.g,LE.0.g,LE.f.g,lifespan.g,
     file=paste0('./Rdata/TreeMort_LifeHistory_Output',index,'.Rdata'))
print(paste0('saved ./Rdata/TreeMort_LifeHistory_Output',index,'.Rdata'))
#___________________________________________________________________________----