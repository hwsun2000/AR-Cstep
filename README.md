# AR-Cstep
An algorithm for penalized maximum trimmed likelihood estimation in mislabeled Omics data. The AR-Cstep algorithm combining the accept-reject algorithm and C-step algorithm solves the problem that C-step algorithm does not converge due to the change of penalty parameters. In order to reach the maximum of truncated likelihood function in algorithm iteration and avoid falling into local optimal,Combines a probabilistic accept-reject algorithm of the Metropolis type.
Installation
The required R package needs to be installed before analysis,You can copy the following code. install.packages(c("mvtnorm","glmnet","glmpath","glmnetUtils","pacman")) library(pacman) p_load(mvtnorm,glmnet,glmpath,glmnetUtils)
Usage
lts.logistic(x,y,gamma=0.25,ncsteps=2,formersubset=500,latersubset=10)
Development
This R code are develpoed by Hongwei Sun and wenjing Zhang.
