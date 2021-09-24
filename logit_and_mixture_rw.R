##############################
#### Logit data and model ####
##############################

# fake data logit reg d=80 ##
N=128
d=2
X=matrix(rnorm(N*d),N,d) #continuous covariates
#X=matrix( sample(0:4, N*d, replace=T) ,N,d) #discrete covariates
X[,1]=rep(1,N)
beta=seq(-2,2,length.out = d)
y = rbinom(N,1, 1/(1+exp(-X %*% beta)) )
#write.table(cbind(y,X),'../Code_cpp/tests/logreg-100-10-categ.txt', row.names=FALSE, col.names = FALSE) #use matlab to turn it into binary

table(y)
GLM<-glm(y ~ . -1 ,family=binomial, data=as.data.frame(X) )
#summary(GLM)
cat(GLM$coef,sep = ", ")



#to get reference posterior
rwLogit=function(y,X,xi0,sigma,run,A){
  ptm=proc.time()
  xi=xi0
  d=length(xi0)
  U=function(q){as.numeric(sum(log(1+exp(X %*% q)))- t(y) %*% (X %*% q) + (t(q) %*% q)/2)}
  tAch=t(chol(A))
  out=matrix(0,run,d)
  count=0
  for(n in 1:run){
    xinew <- xi + sigma * tAch %*% rnorm(d)
    if(runif(1)<exp(U(xi)-U(xinew))){
      xi=xinew
      count=count+1
    }
    if(n%%(run/10)==0) cat("Iteration n.",n,"\n") #track progress1
    out[n,]=xi
  }
  time=proc.time()-ptm
  cat('Run Time',round(time[1]/60,digits=2),'Min','\n')
  cat("Acc Rate",count/run,"\n")
  return(list(out,time[1]))
}

AA=diag(d)
resRW=rwLogit(y, X, GLM$coefficients, 0.8, 100000,  AA)
AA=cov(resRW[[1]]); #AA=diag(3)
resRW=rwLogit(y, X, GLM$coefficients, 2.3, 100000,  AA)
plot.ts(resRW[[1]][1:10000,])
hist(resRW[[1]][,1],breaks=200)



###########################
######### mixture ######### proportions of 2 bivariate normals fixed at .5
###########################

N=100
d=4
X=matrix(NA,N,d/2)
X[1:(N/2),]=matrix(rnorm(N,1.5,1),N/2,d/2)
X[(N/2+1):N,]=matrix(rnorm(N,-1.5,1),N/2,d/2)
plot(X)
hist(X[,1])
hist(X[,2])



#to get reference posterior
rwMixture=function(X,xi0,sigma,run,A){
  ptm=proc.time()
  xi=xi0
  d=length(xi0)
  U=function(q){
    q1=q[1:(d/2)]
    q2=q[(d/2+1):d]
    as.numeric( norm(q1,"2")^2/2 + norm(q2,"2")^2/2 - sum(
      log(.5*exp(-apply(X-q1, 1, function(x) norm(x,"2")^2)/2) + .5*exp(-apply(X-q2, 1, function(x) norm(x,"2")^2)/2)  ) ) ) 
    }
  tAch=t(chol(A))
  out=matrix(0,run,d)
  count=0
  for(n in 1:run){
    xinew <- xi + sigma * tAch %*% rnorm(d)
    if(runif(1)<exp(U(xi)-U(xinew))){
      xi=xinew
      count=count+1
    }
    if(n%%(run/10)==0) cat("Iteration n.",n,"\n") #track progress1
    out[n,]=xi
  }
  time=proc.time()-ptm
  cat('Run Time',round(time[1]/60,digits=2),'Min','\n')
  cat("Acc Rate",count/run,"\n")
  return(list(out,time[1]))
}

AA=diag(d)
resRW=rwMixture(X, c(1,1,1,1), 0.2, 1000,  AA)
AA=cov(resRW[[1]]); #AA=diag(3)
resRW=rwMixture(X, colMeans(resRW[[1]]), .8, 1000,  AA)
plot.ts(resRW[[1]][1:10000,])
hist(resRW[[1]][,1],breaks=200)

####