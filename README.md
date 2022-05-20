# AR-Cstep
An algorithm for penalized maximum trimmed likelihood estimation in mislabeled Omics data. The AR-Cstep algorithm combining the accept-reject algorithm and C-step algorithm solves the problem that C-step algorithm does not converge due to the change of penalty parameters. In order to reach the maximum of truncated likelihood function in algorithm iteration and avoid falling into local optimal,Combines a probabilistic accept-reject algorithm of the Metropolis type.
# Installation
The required R package needs to be installed before analysis,You can copy the following code.
install.packages(c("mvtnorm","glmnet","glmpath","glmnetUtils","MASS","pacman"))
library(pacman) 
p_load(mvtnorm,glmnet,glmpath,glmnetUtils,MASS)
install_github("hwsun2000/AR-Cstep")
# Usage
lts.logistic(x,y,gamma=0.25,ncsteps=2,formersubset=500,latersubset=10)
# Example
set.seed(101) 
p=400  
n=200  
epsilon=0.05  
epsilon1=epsilon*(1/3) 
epsilon2=epsilon*(2/3) 

beta1=rep(1,30)
beta0=rep(0,p-length(beta1))
beta=c(beta1,beta0)

ssigma =0.9  
sigma =ssigma^t(sapply(1:p, function(i, j) abs(i-j), 1:p))
x=rmvnorm(n,sigma=sigma) 
feta = x%*%beta; 
fprob = exp(feta)/(1+exp(feta))
y=rbinom(n,1,fprob)       
index1=which(y==1) 
index0=which(y==0)
n1=length(index1)
n0=length(index0)

k0=sample(index0,floor(n*epsilon1))  
k1=sample(index1,floor(n*epsilon2))   

y[k0]=1 
y[k1]=0 
outlier=c(k0,k1)
     
fit <-lts.logistic(x,y)
beta=as.matrix(fit$rwt.betas)
outliers=(1:n)[-fit$rwt.indices]  
# Development
This R package are develpoed by Hongwei Sun and wenjing Zhang.
