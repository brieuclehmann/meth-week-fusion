#initialise
library(WSGeometry)

########################################
#### Worked out example with 4 dots ####
########################################

pos.list<-vector("list",4)
weights.list<-vector("list",4)

#in each list element are 4 points --> samples from our subposterior
pos.list[[1]]<-matrix(c(0,0,1,1,1,0,0,1),nrow=4,ncol=2)/10
pos.list[[2]]<-matrix(c(9,9,10,10,10,9,9,10),nrow=4,ncol=2)/10
pos.list[[3]]<-matrix(c(9,9,10,10,1,0,0,1),nrow=4,ncol=2)/10
pos.list[[4]]<-matrix(c(0,0,1,1,10,9,9,10),nrow=4,ncol=2)/10

#plot points
plot(0, 0, xlab = "", ylab = "", type = "n", xlim = c(0, 1), ylim = c(0, 1))
for(i in 1:4) points(pos.list[[i]][,1], pos.list[[i]][,2], col = i)

#build weights list (no default), calc and plot barycentre
#weights.list[[1]]<-c(.1,.2,.3,.4)
#weights.list[[2]]<-c(.4,.3,.2,.1)
#weights.list[[3]]<-c(.2,.1,.3,.4)
#weights.list[[4]]<-c(.3,.4,.1,.2)
for(i in 1:4) weights.list[[i]]<-rep(1/4,4)
bary<-barycenter_lp(pos.list,weights.list)
points(bary$positions[,1],bary$positions[,2], col = "orange", pch = 13)





############################################################
#### Wasserstein barycentre and distance for MC samples #### ONLY IN d>2 !!! ####
############################################################

N=20
subposterior_samples=list(); subposterior_hist=list(); weights_list=list()
for(i in 1:4) {
  subposterior_samples[[i]]=matrix(rnorm(2*N),N,2)
  weights_list[[i]]=rep(1,N)/N
  #subposterior_hist[[i]]=hist(subposterior_samples[[i]], freq=F, plot = F)
  #support=range(unlist(subposterior_samples))
  #plot( subposterior_hist[[i]], col=rgb(runif(1),runif(1),runif(1),1/4), xlim=support, ylim = c(0,10), add=ifelse(i==1,F,T))
  }

#plot common background with subposterior samples in different colora
plot(0, 0, xlab = "", ylab = "", type = "n", xlim = c(-5, 5), ylim = c(-5, 5))
for(i in 1:4) points(subposterior_samples[[i]][,1], subposterior_samples[[i]][,2], col = i)

#calculate and plot wasserstein barycentre *sample*
bary<-barycenter_lp(subposterior_samples, weights_list)
points(bary$positions[,1],bary$positions[,2], col = "orange", pch = 13)


#calculate wasserstein distance between two MC samples
P1=wpp(subposterior_samples[[1]], rep(1/N,N))
P2=wpp(subposterior_samples[[2]], rep(1/N,N))
ws_dist(P1,P2)


####