\name{AR-Cstep}
\alias{lts.logistic}
\title{AR-Cstep: An algorithm for penalized maximum trimmed likelihood estimation in mislabeled Omics data}
\usage{
lts.logistic(x,y,gamma=0.25,ncsteps=2,formersubset=500,latersubset=10)
}
\arguments{
\item{x}{a numeric matrix containing the predictor variables.}
\item{y}{a numeric vector containing the binary response variable.}
\item{gamma}{a numeric value giving the percentage of the residuals for which the elastic net-type  penalized likelihood should be maximized(the default is 0.75).}
\item{ncsteps}{a positive integer giving the number of AR-Csteps to perform on all subsamples in the first phase of the algorithm (the default is to perform two C-steps).}
\item{formersubset}{a numeric value giving the number of subsamples in the first phase of  the algorithm. The default is to first perform nstep C-steps on 500 intitial subsamples,  and then to keep the 10 subsamples with the largest likelihood function.}
\item{latersubset}{a numeric value giving the number of subsamples in the second phase of  the algorithm.The default is 10 which means keep 10 subsamples with lowest value of the  objective function after performing nstep C-steps on 500 initial samples.}
}
\description{
The AR-Cstep algorithm combining the accept-reject algorithm and C-step algorithm solves the problem that C-step algorithm does not converge due to the change of penalty parameters. In order to reach the maximum of truncated likelihood function in algorithm iteration and avoid falling into local optimal,Combines a probabilistic accept-reject algorithm of the Metropolis type.
}
