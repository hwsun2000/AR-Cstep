#Title:AR-Cstep: An algorithm for penalized maximum trimmed likelihood estimation in mislabeled Omics data

#Date:2020-6-17
#Author:Hongwei Sun


#variable list
#x             a numeric matrix containing the predictor variables
#y             a numeric vector containing the binary response variable.
#h             a numeric value giving the percentage of the residuals for which the elastic net-type  penalized likelihood should be maximized(the default is 0.75).
#ncsteps       a positive integer giving the number of AR-Csteps to perform on all subsamples in the first phase of the algorithm (the default is to perform two C-steps)
#formersubset  a numeric value giving the number of subsamples in the first phase of  the algorithm. The default is to first perform nstep C-steps on 500 intitial subsamples,  and then to keep the 10 subsamples with the largest likelihood function.
#latersubset   a numeric value giving the number of subsamples in the second phase of  the algorithm.The default is 10 which means keep 10 subsamples with lowest value of the  objective function after performing nstep C-steps on 500 initial samples.


#' AR-Cstep: An algorithm for penalized maximum trimmed likelihood estimation in mislabeled Omics data
#'
#' @param x a numeric matrix containing the predictor variables.
#' @param y a numeric vector containing the binary response variable.
#' @param gamma a numeric value giving the percentage of the residuals for which the elastic net-type  penalized likelihood should be maximized(the default is 0.75).
#' @param ncsteps a positive integer giving the number of AR-Csteps to perform on all subsamples in the first phase of the algorithm (the default is to perform two C-steps).
#' @param formersubset a numeric value giving the number of subsamples in the first phase of  the algorithm. The default is to first perform nstep C-steps on 500 intitial subsamples,  and then to keep the 10 subsamples with the largest likelihood function.
#' @param latersubset a numeric value giving the number of subsamples in the second phase of  the algorithm.The default is 10 which means keep 10 subsamples with lowest value of the  objective function after performing nstep C-steps on 500 initial samples.
#'
#' @return rwt.betas，lpredict，rwt.indices


library(glmpath)
library(glmnetUtils)
library(MASS)

lts.logistic<-function(x,y,gamma=0.25,ncsteps=2,formersubset=500,latersubset=10 ){


  n<-length(y)  # sample size
  p<-dim(x)[2]  # number of variables,

   center=apply(x,2,mean)
   x.c = sweep(x, 2, center)
   sd.c=apply(x,2,sd)
   xsd=sweep(x.c,2,sd.c,"/")

  hn=floor((n+1)*(1-gamma))
  n1=length(which(y==1))
  n0=length(which(y==0))
  index1<-which(y==1)
  index0<-which(y==0)

  sn1=floor((n1+1)*(1-gamma))
  sn0=hn-sn1

  kmax=1000
  rseqmax=3
  D=0.05
  max.lik=-10000

  #robustfied point-biserial correlation or Pearson correlation

  rob_r=rep(0,p)
  for (j in 1:p)
   {
  if (mad(x[,j])!=0)
   {rob_r[j]=sqrt(n0*n1/(n*(n-1)))*(median(x[index1,j])-median(x[index0,j]))/mad(x[,j])
   }   else
   {  rob_r[j]=cor(x[,j],y) }
   }
   maxcc=max(abs(rob_r))   # the largest value of lambda
   frac=seq(0.05,1,by=0.05)*maxcc  # lambda sequence for AR-Cstep.

   ##two-step AR-Cstep

   cstep<-function(indices){
   for (i in 1:ncsteps){

    center=apply(x[indices,],2,mean)
    x.c = sweep(x, 2, center)
    sd.c=apply(x[indices,],2,sd)
    xsd=sweep(x.c,2,sd.c,"/")

     EN=cv.glmnet(xsd[indices,],y[indices],standardize=FALSE,lambda=frac,alpha=0.5,family="binomial")
     lpredict= predict(EN,xsd[1:n,],type="response")  #predictive probablility
     ww=coef(EN)       #estimated coefficients
     mybetas=as.matrix(ww)
     res=y*log(lpredict)+(1-y)*log((1-lpredict))   #likelihood function

     o1<-index1[order(res[index1],decreasing=TRUE)[1:sn1]]
     o0<-index0[order(res[index0],decreasing=TRUE)[1:sn0]]
     indices<-sort(c(o1,o0))
   }
     crit=sum(res[indices])
     result<-list(indices=indices,crit=crit,mybetas=mybetas)
     return(result)
    }


   subset2=matrix(0,6,formersubset)

   for (j in 1:formersubset)
  {
   k1=sample(index1,3)
   k0=sample(index0,3)
   subset2[,j]=t(c(k1,k0))
   }

   coef=matrix(0,hn,ncol(subset2))
   output=rep(0,ncol(subset2))
   xishu=matrix(0,p+1,ncol(subset2))

   ##The first phase: perform nstep AR-Csteps on 500 intitial subsamples

   for (k in 1:ncol(subset2))
   {
    jieguo<-cstep(subset2[,k])
    output[k]<-jieguo$crit    #likelihood function
    coef[,k]<-jieguo$indices  #subscript of subsamples
    xishu[,k]<-jieguo$mybetas  # coefficient
    }

   ##keep latersubset subsets with largest value of the likelihood function
     betterindices=order(output,decreasing=TRUE)[1:latersubset]
     a=coef[,betterindices]

     ##Cstep2 function: AR-Cstep until convergence

     cstep2<-function(indices){

     ##performed penalized regression on the first subset, and retain the hn samples with largest likelihood as optimal subset at current step.
     poscaso=indices      ###the optimal subset

     #normalize by trimmed mean and standard

      center.trim=apply(x[poscaso,],2,mean)    ##standardization
      sd.trim=apply(x[poscaso,],2,sd)
      sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
      center.tc<-ifelse(sd.trim==0,center,center.trim)
      x.c = sweep(x, 2, center.tc)
      xsd=sweep(x.c,2,sd.tc,"/")

      xsd.sub=xsd[poscaso,]
      EN=cv.glmnet(x=xsd.sub,y=y[poscaso],family='binomial',alpha=0.5,standardize=FALSE)
      beta=as.matrix(coef(EN,s="lambda.1se"))
      lpredict=as.matrix(predict(EN,newx=xsd[1:n,],type="response"))
      lpredict[lpredict<=0.0000001]=0.0000001
      lpredict[lpredict>=0.9999999]=0.9999999

      res=y*log(lpredict)+(1-y)*log((1-lpredict))  #likelihood
      lik=sum(res[poscaso])

      likold=-10000
       k=1
      bestlik=lik
      b2=beta
      p2=poscaso
      rseq=0     #the iteration number of the subset was not replaced

      while(k<kmax && rseq<rseqmax) {

      Tk=log(k+1,2)/D
      k=k+1

      if(lik-likold!=0) {rseq=0}
      if(lik-likold==0) {rseq=rseq+1}

      likold=lik
      o1<-index1[order(res[index1],decreasing=TRUE)[1:sn1]]
      o0<-index0[order(res[index0],decreasing=TRUE)[1:sn0]]

      pos2<-sort(c(o1,o0))

     ## pos2 is composed of samples with the hn largest likelihood. Below solve the likelihood function corresponding to pos2.
     #normalize by trimmed mean and standard

      center.trim=apply(x[pos2,],2,mean)
      sd.trim=apply(x[pos2,],2,sd)
      sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
      center.tc<-ifelse(sd.trim==0,center,center.trim)
       x.c = sweep(x, 2, center.tc)
       xsd=sweep(x.c,2,sd.tc,"/")

      xsd.iter=xsd[pos2,]
      y.iter=y[pos2]

      EN=cv.glmnet(x=xsd.iter,y=y.iter,family='binomial',alpha=0.5,standardize=FALSE)
      betaff=as.matrix(coef(EN,s="lambda.1se"))
      lpredict=as.matrix(predict(EN,newx=xsd[1:n,],type="response"))
      lpredict[lpredict<=0.0000001]=0.0000001
      lpredict[lpredict>=0.9999999]=0.9999999
       res=y*log(lpredict)+(1-y)*log((1-lpredict))

      lik.cand=sum(res[pos2])     ##the likelihood function corresponding to pos2

      pii=min(exp(Tk*(lik.cand-lik)),1)

      if(rbinom(1,1,pii)==1)  {poscaso=pos2       #The subscript of subset pos2 is assigned to poscaso, so that the likelihhod can be calculated after iteration
                         likold=lik              #the original likelihood lik is assigned to likold
                         lik=lik.cand            #Assign likelihood function lik.cand of new subset to lik
                         beta=betaff            ## Assign the coefficients of the new subset to beta
                         res=res                #Record the new likelihood function res of each individual

    if(lik>bestlik) {
    bestlik=lik
    b2=beta
    p2=pos2
               }

     }   #If the likelihood function is better than bestlik, record it.

    }     ##Complete the iteration process

lik=bestlik
beta=b2
pos2=p2

if(lik>max.lik)     ####Compared with the likelihood function corresponding to the initial subset.
{max.lik=lik
max.pos2=pos2
max.beta=beta
        }

result<-list(max.pos2=max.pos2, max.lik=max.lik,max.beta=max.beta)
      return(result)

       }



      output3=rep(0,latersubset)
     rawindices=matrix(0,hn,latersubset)
      mybetas3=matrix(0,p+1,latersubset)


   ######### for latersubset subsets, perform cstep2 function

    for (k in 1:latersubset)
    {
   jieguo2<-cstep2(a[,k])
  output3[k]<-jieguo2$max.lik     #likelihood function
   rawindices[,k]<-jieguo2$max.pos2  #subscript of subsamples
   mybetas3[,k]<-jieguo2$max.beta # coefficient
    }


    #keep the optimal subset with largest likelihood
    lastindices=which.max(output3)
    raw_indices=rawindices[,lastindices]
    lastcrit=output3[lastindices]

    # perform EN on the optimal subset

      center.trim=apply(x[raw_indices,],2,mean)
      sd.trim=apply(x[raw_indices,],2,sd)
      sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
      center.tc<-ifelse(sd.trim==0,center,center.trim)
       x.c = sweep(x, 2, center.tc)
       xsd=sweep(x.c,2,sd.tc,"/")

       xsd_pre<-cbind(rep(1,dim(xsd)[1]),xsd)

       my.alpha <- seq(0.1,1,0.1)
      cvm=rep(0,length(my.alpha))
       var.selected.EN <- matrix(0,p+1,length(my.alpha))
       link<- matrix(0,n,length(my.alpha))

     fit.EN.cv=cva.glmnet(x=xsd[raw_indices,],y=y[raw_indices],alpha=my.alpha,family="binomial",standardize=FALSE)
      for (j in 1:length(my.alpha)){
           cvm[j]=fit.EN.cv$modlist[[j]]$cvm[which(fit.EN.cv$modlist[[j]]$cvm==min(fit.EN.cv$modlist[[j]]$cvm))]
          var.selected.EN[,j]=as.matrix(coef(fit.EN.cv$modlist[[j]],s="lambda.min"))
          link[,j]=xsd_pre%*%var.selected.EN[,j]
       }

    raw.betas=var.selected.EN[,which(cvm==min(cvm))]
    link_raw=link[,which(cvm==min(cvm))]
    rawpredict<-exp(link_raw)/(1+exp(link_raw))

     # reweighted step
     perres=(y-rawpredict)/sqrt(rawpredict*(1-rawpredict))
     q <- qnorm(0.9875)
     ok<- as.integer(abs(perres) <= q)
     rwt.indices<-sort(which(ok>0))

     # perform EN on the reweighted subset
     if (sum(y[rwt.indices])>1&sum(y[rwt.indices])<length(y[rwt.indices])-1)
      {

      center.trim=apply(x[rwt.indices,],2,mean)
      sd.trim=apply(x[rwt.indices,],2,sd)
      sd.tc<-ifelse(sd.trim==0,sd.c,sd.trim)
      center.tc<-ifelse(sd.trim==0,center,center.trim)
       x.c = sweep(x, 2, center.tc)
       xsd=sweep(x.c,2,sd.tc,"/")
       xsd_pre<-cbind(rep(1,dim(xsd)[1]),xsd)

     my.alpha <- seq(0.1,1,0.1)
      cvm=rep(0,length(my.alpha))
      var.selected.EN <- matrix(0,p+1,length(my.alpha))
       link<- matrix(0,n,length(my.alpha))

     fit.EN.cv=cva.glmnet(x=xsd[rwt.indices,],y=y[rwt.indices],alpha=my.alpha,family="binomial",standardize=FALSE)

      for (j in 1:length(my.alpha)){

           cvm[j]=fit.EN.cv$modlist[[j]]$cvm[which(fit.EN.cv$modlist[[j]]$cvm==min(fit.EN.cv$modlist[[j]]$cvm))]

          var.selected.EN[,j]=as.matrix(coef(fit.EN.cv$modlist[[j]],s="lambda.min"))
          link[,j]=xsd_pre%*%var.selected.EN[,j]

       }

    rwt.betas=var.selected.EN[,which(cvm==min(cvm))]
    link_rwt=link[,which(cvm==min(cvm))]
   lpredict<-exp(link_rwt)/(1+exp(link_rwt))

   }  else
    {  rwt.betas=raw.betas
       lpredict=rawpredict
    }


   # residual diagnostic
   perres2=(y-lpredict)/sqrt(lpredict*(1-lpredict))
   q2<- qnorm(0.9875)
   ok2<- as.integer(abs(perres2)<=q2)
   rwt.indices2<-sort(which(ok2>0))


   result<-list(lastindices=lastindices,lastcrit=lastcrit,raw_indices=raw_indices,perres=perres,raw.betas=raw.betas,rawpredict=rawpredict,ok=ok,rwt.indices=rwt.indices,rwt.betas=rwt.betas,lpredict=lpredict,perres2=perres2,rwt.indices2=rwt.indices2)
   return(result)

    }




