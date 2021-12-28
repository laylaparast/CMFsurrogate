
Kern.FUN <- function(zz,zi,bw) 
{ 
  out = (VTM(zz,length(zi))- zi)/bw
  dnorm(out)/bw
  
}


Kern.FUN.M <- function(zz,zi,bw) 
{ 
  out=NULL
  pro=1
  for (j in 1:ncol(zz)){
  out[[j]] = (VTM(zz[,j],nrow(zi))- zi[,j])/bw[j] 
  pro=pro*(dnorm(out[[j]])/bw[j]) 
  }
  pro
  
}


Kern.FUN.d <- function(zz,zi,bw) 
{ 
  out=NULL
  pro=1
  for (j in 1:ncol(zz)){
    out[[j]] = (VTM(zz[,j],nrow(zi))- zi[,j])/bw[j] 
    pro=pro*(dnorm(out[[j]])/bw[j])*(-out[[j]])/bw[j]
  }
  pro
  
}


VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}



gen.bootstrap.weights=function( n, num.perturb=500){
  sapply(1:num.perturb,function(x) sample(1:n,n,replace=T))

}


resam<- function(index,yob,sob,aob,n){
  yob=yob[index]
  sob=sob[index,]
  aob=aob[index]
  ######## parametric all the way
  sobsob=cbind(sob[, 1],sob[, 2],sob[, 3],sob[, 4],sob[, 1]*sob[, 2],sob[, 1]*sob[, 3],sob[, 1]*sob[, 4],
               sob[, 2]*sob[, 3],sob[, 2]*sob[, 4],sob[, 3]*sob[, 4],sob[,1]^2,sob[,2]^2,sob[,3]^2,sob[,4]^2)#
  
  y=yob[aob==0]
  # x=sob[aob==0,]
  # temp=lm(y~x)
  # m0.ob=cbind(rep(1,nrow(sob)),sob)%*%temp$coefficients
  xx=sobsob[aob==0,]
  temp=lm(y~xx)
  m0.ob=cbind(rep(1,nrow(sob)),sobsob)%*%temp$coefficients
  
  y=yob[aob==1]
  # x=sob[aob==1,]
  # temp=lm(y~x)
  # m1.ob=cbind(rep(1,nrow(sob)),sob)%*%temp$coefficients
  xx=sobsob[aob==1,]
  temp=lm(y~xx)
  m1.ob=cbind(rep(1,nrow(sob)),sobsob)%*%temp$coefficients
  
  temp=glm(aob~sobsob,family=binomial)
  p1=temp$fitted.values
  p0=1-p1
  
  ms.p=m1.ob*p1+m0.ob*p0
  c.hat=mean(yob*(1-aob))/mean(1-aob)-mean(ms.p*(1-aob))/mean(1-aob)
  temp=mean(p0*(1-aob))/mean((1-aob))
  gs.p=ms.p+p0*c.hat/temp
  
  # causal=mean(yob*aob)/mean(aob)-mean(yob*(1-aob))/mean(1-aob)
  # causals=mean(gs.p*aob)/mean(aob)-mean(gs.p*(1-aob))/mean((1-aob))
  # pte.pa=causals/causal
  
  ################ additive linear 
  sobs=cbind(ns(sob[,1],df=4),ns(sob[,2],df=4),ns(sob[,3],df=4),ns(sob[,4],df=4))
  
  temp=t(sobs)%*%(sobs*aob)/sum(aob)+t(sobs)%*%(sobs*(1-aob))/sum(1-aob)-
    matrix(apply(sobs*aob,2,sum)/sum(aob),ncol=1)%*%apply(sobs*(1-aob),2,sum)/sum(1-aob)-
    matrix(apply(sobs*(1-aob),2,sum)/sum(1-aob),ncol=1)%*%apply(sobs*aob,2,sum)/sum(aob)
  temp2=apply(sobs*yob*aob,2,sum)/sum(aob)+apply(sobs*yob*(1-aob),2,sum)/sum(1-aob)-
    sum(yob*aob)/sum(aob)*apply(sobs*(1-aob),2,sum)/sum(1-aob)-
    sum(yob*(1-aob))/sum(1-aob)*apply(sobs*aob,2,sum)/sum(aob)
  beta=ginv(temp)%*%temp2
  gs.l=sobs%*%beta
  
  # causal=mean(yob*aob)/mean(aob)-mean(yob*(1-aob))/mean(1-aob)
  # causals=mean(gs.l*aob)/mean(aob)-mean(gs.l*(1-aob))/mean((1-aob))
  # pte.l=causals/causal
  
  ######## convex combination
  model.new <- function(temp) {
    sob.r=as.numeric(temp*gs.p+(1-temp)*gs.l)
    bw = 1.06*sd(sob.r)*n^(-1/5)/(n^0.2)
    kern = Kern.FUN(zz=sob.r,zi=sob.r,bw)
    ms.r=apply(yob*kern,2,sum)/apply(kern,2,sum)
    c.hat=mean(yob*(1-aob))/mean(1-aob)-mean(ms.r*(1-aob))/mean(1-aob)
    p0s.hat=apply((1-aob)*kern,2,sum)/apply(kern,2,sum)
    temp=mean(p0s.hat*(1-aob))/mean((1-aob))
    gs.r=ms.r+p0s.hat*c.hat/temp

    causal=mean(yob*aob)/mean(aob)-mean(yob*(1-aob))/mean(1-aob)
    causals=mean(gs.r*aob)/mean(aob)-mean(gs.r*(1-aob))/mean((1-aob))
    -causals/causal
  }
  # opt=optimize(model.new,lower=0, upper=1)
  # pte.convex=-opt$objective
  opt=nlm(model.new, 0, hessian = FALSE, iterlim = 100)
  pte.convex=-opt$minimum
  pte.pa=-model.new(1)
  pte.l=-model.new(0)
  
  out=c( pte.pa,pte.l,pte.convex)
  
}


pte.estimate.multiple = function(sob, yob, aob, var = TRUE ,  rep=500) {
  n=length(yob)
    
  ################ parametric all the way
  sobsob=cbind(sob[, 1],sob[, 2],sob[, 3],sob[, 4],sob[, 1]*sob[, 2],sob[, 1]*sob[, 3],sob[, 1]*sob[, 4],
               sob[, 2]*sob[, 3],sob[, 2]*sob[, 4],sob[, 3]*sob[, 4],sob[,1]^2,sob[,2]^2,sob[,3]^2,sob[,4]^2)#
  
  y=yob[aob==0]
  xx=sobsob[aob==0,]
  temp=lm(y~xx)
  m0.ob=cbind(rep(1,nrow(sob)),sobsob)%*%temp$coefficients
  
  y=yob[aob==1]
   xx=sobsob[aob==1,]
  temp=lm(y~xx)
  m1.ob=cbind(rep(1,nrow(sob)),sobsob)%*%temp$coefficients
  
  temp=glm(aob~sobsob,family=binomial)
  p1=temp$fitted.values
  p0=1-p1
  
  ms.p=m1.ob*p1+m0.ob*p0
  c.hat=mean(yob*(1-aob))/mean(1-aob)-mean(ms.p*(1-aob))/mean(1-aob)
  temp=mean(p0*(1-aob))/mean((1-aob))
  gs.p=ms.p+p0*c.hat/temp
  
  causal=mean(yob*aob)/mean(aob)-mean(yob*(1-aob))/mean(1-aob)
  causals=mean(gs.p*aob)/mean(aob)-mean(gs.p*(1-aob))/mean((1-aob))
  pte.pa=causals/causal
  
  ################ additive linear
  sobs=cbind(ns(sob[,1],df=4),ns(sob[,2],df=4),ns(sob[,3],df=4),ns(sob[,4],df=4))
  
  temp=t(sobs)%*%(sobs*aob)/sum(aob)+t(sobs)%*%(sobs*(1-aob))/sum(1-aob)-
    matrix(apply(sobs*aob,2,sum)/sum(aob),ncol=1)%*%apply(sobs*(1-aob),2,sum)/sum(1-aob)-
    matrix(apply(sobs*(1-aob),2,sum)/sum(1-aob),ncol=1)%*%apply(sobs*aob,2,sum)/sum(aob)
  temp2=apply(sobs*yob*aob,2,sum)/sum(aob)+apply(sobs*yob*(1-aob),2,sum)/sum(1-aob)-
    sum(yob*aob)/sum(aob)*apply(sobs*(1-aob),2,sum)/sum(1-aob)-
    sum(yob*(1-aob))/sum(1-aob)*apply(sobs*aob,2,sum)/sum(aob)
  beta=ginv(temp)%*%temp2
  gs.l=sobs%*%beta
  
  causal=mean(yob*aob)/mean(aob)-mean(yob*(1-aob))/mean(1-aob)
  causals=mean(gs.l*aob)/mean(aob)-mean(gs.l*(1-aob))/mean((1-aob))
  pte.l=causals/causal
  
  ######## convex combination
  model.new <- function(temp) {
    sob.r=as.numeric(temp*gs.p+(1-temp)*gs.l)
    bw = 1.06*sd(sob.r)*n^(-1/5)/(n^0.2)
    kern = Kern.FUN(zz=sob.r,zi=sob.r,bw)
    ms.r=apply(yob*kern,2,sum)/apply(kern,2,sum)
    c.hat=mean(yob*(1-aob))/mean(1-aob)-mean(ms.r*(1-aob))/mean(1-aob)
    p0s.hat=apply((1-aob)*kern,2,sum)/apply(kern,2,sum)
    temp=mean(p0s.hat*(1-aob))/mean((1-aob))
    gs.r=ms.r+p0s.hat*c.hat/temp
    
    causal=mean(yob*aob)/mean(aob)-mean(yob*(1-aob))/mean(1-aob)
    causals=mean(gs.r*aob)/mean(aob)-mean(gs.r*(1-aob))/mean((1-aob))
    -causals/causal
  }
   opt=nlm(model.new, 0, hessian = FALSE, iterlim = 100)
  pte.convex=-opt$minimum
  index.max=opt$estimate
  pte.pa=-model.new(1)
  pte.l=-model.new(0)
  
   #### variance
  if (var==T){
    re=rep
    index=gen.bootstrap.weights( n, num.perturb=re)
    temp=apply(index,2,resam,yob,sob,aob,n)
    pte.se=apply(temp,1,sd)[3]
    
    return(list('pte.es'=pte.convex,'pte.se'=pte.se))
  }else{
    return(list('pte.es'=pte.convex))
  }
  
  }
  

 

 

 

 

