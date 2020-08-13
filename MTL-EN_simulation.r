


install.packages("mvtnorm")
install.packages("glmnet")


library(mvtnorm)
library(glmnet)


n=500
p=1000
nn=5

beta1=rep(1,30)

beta0=rep(0,p-length(beta1))
beta=c(beta1,beta0)

beta1 <- which(beta!=0)##index variable
beta0 <-which(beta==0)  


     lt2Sp=rep(0,nn)     
     lt2Sn=rep(0,nn)    
      lt2num=rep(0,nn)  

     lt2error=rep(0,nn)  
     lt2rmspe=rep(0,nn)


     onum=rep(0,nn)    
      power=rep(0,nn)  
       Ierror=rep(0,nn)  
    Yorden=rep(0,nn)  

     onum2=rep(0,nn)   
     power2=rep(0,nn) 
     Ierror2=rep(0,nn)  
     Yorden2=rep(0,nn)  

     data=list()
     data1=list()

     file=list()
     file1=list()

     txtfile=paste(1:nn,"train.txt",sep="")
     txtfile1=paste(1:nn,"test.txt",sep="")   


       for(m in 1:nn)
       {

       file[[m]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile[m],sep="")

       data[[m]]=read.table(file=file[[m]],header=TRUE)

      y=as.matrix(data[[m]][1])

      outlier_indi<-as.matrix(data[[m]][2])
     outlier_indi=as.vector(outlier_indi)

       routliers<-which(outlier_indi==1)
       rinliers<-which(outlier_indi==0)
    
      x=as.matrix(data[[m]][-c(1,2)])

        file1[[m]]=paste("D:/2019/simulation/n500p1000eps5%/",txtfile1[m],sep="")

       data1[[m]]=read.table(file=file1[[m]],header=TRUE)

       y1=as.matrix(data1[[m]][1])
       x1=as.matrix(data1[[m]][-1])


     center=apply(x,2,mean)
      x.c = sweep(x, 2, center)
      sd.c=apply(x,2,sd)
      xsd=sweep(x.c,2,sd.c,"/")   


     center=apply(x1,2,mean)
      x.c = sweep(x1, 2, center)
      sd.c=apply(x1,2,sd)
      xsd1=sweep(x.c,2,sd.c,"/")   

    xsd1_pre<-cbind(rep(1,dim(xsd1)[1]),xsd1)


    fit <-lts.logistic(x,y)
                
    l=as.matrix(fit$rwt.betas)
    ll2=as.matrix(l)[-1] 


     lt2predict=exp(xsd1_pre%*%l)/(1+exp(xsd1_pre%*%l))  

        lt2predict2=as.matrix(fit$lpredict)   

       lt2num[m]=length(which(ll2!=0))   

    
       TP=length(which(ll2[beta1]!=0))  
       FN=length(which(ll2[beta1]==0))  
       TN=length(which(ll2[beta0]==0))   
       FP=length(which(ll2[beta0]!=0))   
       lt2Sn[m]=TP/(TP+FN)
       if((TP+FP)==0)
        {lt2Sp[m]=0} else
      {lt2Sp[m]=TP/(TP+FP)}
       
       lt2rmspe[m]=sqrt(sum((y1-lt2predict)^2)/n)
 
    lt2yhat=rep(0,n)
    lt2yhat[which(lt2predict>=0.5)]=1
    lt2error[m]= (n-length(which(y1==lt2yhat)))/n  

      
       outliers=(1:n)[-fit$rwt.indices]           
       onum[m]=length(outliers)   

      outlier_ident=rep(0,n)
      outlier_ident[outliers]=1      

      odent=outliers
      ident=(1:n)[-odent]

      TPresult<- length(intersect(odent,routliers)) 
	FPresult<- length(intersect(odent,rinliers)) 
      TNresult<- length(intersect(ident,rinliers)) 
	FNresult<- length(intersect(ident,routliers)) 

     power[m]=TPresult/length(routliers)
     Ierror[m]=FPresult/length(rinliers)  
     Yorden[m]=power[m]-Ierror[m]
   
      
      yhat<-rep(0,dim(x)[1])       
      yhat[which(lt2predict2>=0.5)]=1           
     outliers2<-which(y!=yhat)     



    onum2[m]=length(outliers2)   


 }

 lt2SnSp=sqrt(lt2Sn*lt2Sp)


paste("mean(onum) is",mean(onum))
paste("mean(power) is",mean(power))
paste("mean(Ierror) is",mean(Ierror))
paste("mean(Yorden) is",mean(Yorden))
paste("mean(num) is",mean(lt2num))
paste("mean(fdr) is",mean(1-lt2Sp))
paste("mean(Sn) is",mean(lt2Sn))
paste("mean(SnSp) is",mean(lt2SnSp))
paste("mean(MSE) is",mean(lt2rmspe))
paste("mean(Miscl) is",mean(lt2error))



paste("sd(onum) is",sd(onum))
paste("sd(power) is",sd(power))
paste("sd(Ierror) is",sd(Ierror))
paste("sd(Yorden) is",sd(Yorden))
paste("sd(num) is",sd(lt2num))
paste("sd(fdr) is",sd(1-lt2Sp))
paste("sd(Sn) is",sd(lt2Sn))
paste("sd(SnSp) is",sd(lt2SnSp))
paste("sd(MSE) is",sd(lt2rmspe))
paste("sd(Miscl) is",sd(lt2error))





