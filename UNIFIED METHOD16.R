rm(list=ls());
##install.packages("plot3D");
##library("plot3D");
library(combinat);
library(gtools);
## kernels
##-----1 Uniform Kernel---------------------
UniKer=function(u)
{
out=1/2+0*u;
return(out);
}
##---2 Triangular Kernel-------------------
TriKer=function(u)
{
out=1-abs(u);
return(out);
}

##----3 Biweight (Quartic) Kernel---------------------------
Biwt=function(u)   ##Biwt:Biweight (or quartic)
{
out=15/16*(1-u^2)^2;
return(out);
}

##--4 Epanechnikov (parabolic) Kernel------------
EpaKer=function(u)
{
out=3/4*(1-u^2);
return(out);
}
##----- Triweight Kernel--------------------------
##Triwt=function(u)
##{
##out=35/32*(1-u^2)^3;
##return(out);
##}
##----5 Tricube Kernel ---------------------------
Tricube=function(u)
{
out=70/81*(1-abs(u)^3)^3;
return(out);
}

##----6 Cosine Kernel-----------------------------
CosineKer=function(u)
{
out=pi/4*cos(pi/2*u);
return(out);
}





###  paned indexes within a winodw
#####################################################################
PanedVecs=function(globalInputVec,globalOutputVec,centre,radius)
{
l_window=centre-radius;
r_window=centre+radius;
trueFalseVec=globalInputVec>=l_window & globalInputVec<=r_window;
trueVecInd=which(trueFalseVec==TRUE);   #trueVecInd: the indexes of the true values
localInputVec=globalInputVec[trueVecInd];
localOutputVec=globalOutputVec[trueVecInd];
sLocalInputVec=(localInputVec-centre)/radius;   ##sLocalInputVec:scaled paned-X vector
return(list(sLocalInputVec,localOutputVec));
}

####bounds for a given set of data
Bounds=function(panedVecs,centre,radius,alpha,tau,rho)
{
xVec=panedVecs[[1]];
yVec=panedVecs[[2]];
len_xVec=length(xVec);
uVec=vector(,len_xVec+2);
uVec[1]=-1; uVec[len_xVec+2]=1;
uVec[2:(len_xVec+1)]=xVec;
bounds=matrix(nrow=2,ncol=len_xVec);
   for(i in 1:len_xVec)
   {
   lowerbound=uVec[i+1]-alpha*tau(uVec[i],uVec[i+1]);
   upperbound=uVec[i+1]+alpha*tau(uVec[i+1],uVec[i+2]);
   bounds[1,i]=lowerbound;
   bounds[2,i]=upperbound;
   }
return(bounds);
}

#### weights
##CASE 1
Weight=function(panedVecs,Kernel,bounds,alpha)   ##bounds is a 2-by-n matrix, where n is the size of the samples
{
xVec=panedVecs[[1]];
yVec=panedVecs[[2]];
len_xVec=length(xVec);
uVec=vector(,len_xVec+2);
uVec[1]=-1; uVec[len_xVec+2]=1;
uVec[2:(len_xVec+1)]=xVec;

n=dim(bounds)[2];
weight=vector(,n);
  for(i in 1:n)
  {
  area=integrate(Kernel,bounds[1,i],bounds[2,i]);
  area=as.numeric(area[1]);
  weight[i]=1/alpha*1/(uVec[i+2]-uVec[i])*area;
  }
return(unlist(weight));
}
##CASE 2
Weight=function(panedVecs,Kernel,bounds,alpha)   ##bounds is a 2-by-n matrix, where n is the size of the samples
{
xVec=panedVecs[[1]];
yVec=panedVecs[[2]];
len_xVec=length(xVec);
uVec=vector(,len_xVec+2);
uVec[1]=-1; uVec[len_xVec+2]=1;
uVec[2:(len_xVec+1)]=xVec;

n=dim(bounds)[2];
weight=vector(,n);
  for(i in 1:n)
  {
  area=integrate(Kernel,bounds[1,i],bounds[2,i]);
  area=as.numeric(area[1]);
  weight[i]=1/alpha*area;
  }
return(unlist(weight));
}

###unifed weights
UniWeight=function(weight)
{
return(weight/sum(weight))
}

## CASE 3
Weight=function(panedVecs,Kernel,bounds,alpha)   ##bounds is a 2-by-n matrix, where n is the size of the samples
{
xVec=panedVecs[[1]];
yVec=panedVecs[[2]];
len_xVec=length(xVec);
uVec=vector(,len_xVec+2);
uVec[1]=-1; uVec[len_xVec+2]=1;
uVec[2:(len_xVec+1)]=xVec;

n=dim(bounds)[2];
weight=vector(,n);
  for(i in 1:n)
  {
  area=integrate(Kernel,bounds[1,i],bounds[2,i]);
  area=as.numeric(area[1]);
  weight[i]=1/alpha*1/(uVec[i+1]-uVec[i])*area;
  }
return(unlist(weight));
}


## CASE 4
Weight=function(panedVecs,Kernel,bounds,alpha)   ##bounds is a 2-by-n matrix, where n is the size of the samples
{
xVec=panedVecs[[1]];
yVec=panedVecs[[2]];
len_xVec=length(xVec);
uVec=vector(,len_xVec+2);
uVec[1]=-1; uVec[len_xVec+2]=1;
uVec[2:(len_xVec+1)]=xVec;

n=dim(bounds)[2];
weight=vector(,n);
  for(i in 1:n)
  {
  area=integrate(Kernel,bounds[1,i],bounds[2,i]);
  area=as.numeric(area[1]);
  weight[i]=1/alpha*1/(3*abs(uVec[i+2]-uVec[i+1])+2*abs(uVec[i+1]-uVec[i]))*area;
  }
return(unlist(weight));
}

## CASE 5
Weight=function(panedVecs,Kernel,bounds,alpha)   ##bounds is a 2-by-n matrix, where n is the size of the samples
{
xVec=panedVecs[[1]];
yVec=panedVecs[[2]];
len_xVec=length(xVec);
uVec=vector(,len_xVec+2);
uVec[1]=-1; uVec[len_xVec+2]=1;
uVec[2:(len_xVec+1)]=xVec;

n=dim(bounds)[2];
weight=vector(,n);
  for(i in 1:n)
  {
  area=integrate(Kernel,bounds[1,i],bounds[2,i]);
  area=as.numeric(area[1]);
  weight[i]=1/alpha*1/(3*abs(uVec[i+2]-uVec[i+1])+2*abs(uVec[i+1]-uVec[i]))*area;
  }
return(unlist(weight));
}

###unifed weights
UniWeight=function(weight)
{
return(weight/sum(weight))
}
### discretise the continuous fitted function to form vector of centres
VecCentre=function(globalInputVec,numberIntervals)
{
numInt=numberIntervals;                        ##number of intervals
min=min(globalInputVec);
max=max(globalInputVec);
lenInt=(max-min)/numInt;                       ##length of each interval
vecCentre=seq(min,max,by=lenInt);
return(vecCentre);
}


## calculate the y values corresponding to vecCentre;
Yvalues=function(globalInputVec,globalOutputVec,Kernel,radius,alpha,tau,rho,numberIntervals)
{
vecCentre=VecCentre(globalInputVec,numberIntervals);
lenVecCentre=length(vecCentre);
yvalues=vector(,lenVecCentre);
  for(k in 1:lenVecCentre)
  {
  centre=vecCentre[k];
  panedVecs=PanedVecs(globalInputVec,globalOutputVec,centre,radius);
  bounds=Bounds(panedVecs,centre,radius,alpha,tau,rho);
  weight=Weight(panedVecs,Kernel,bounds,alpha);
  uniWeight=UniWeight(weight);
  yvec=unlist(panedVecs[2]);
  yvalues[k]=uniWeight%*%yvec;
  }
return(yvalues);
}


###mean of sum of squared errors
Mse=function(globalInputVec,globalOutputVec,Kernel,radius,alpha,tau,rho,numberIntervals)   ##mean of sum of squared errors
{
vecCentre=VecCentre(globalInputVec,numberIntervals);
matc=match(vecCentre,globalInputVec);
ind=which(matc!="NA");
yvalues=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha,tau,rho,numberIntervals);
fittedY=yvalues[ind];
sse=(globalOutputVec-fittedY)%*%(globalOutputVec-fittedY);
mse=sse/length(globalInputVec);
return(mse);
}



##############################################################################
#########################      data analysis    ##############################
##############################################################################
##simulated data
X=runif(26,30,100);
Y=runif(26,20,32);
X=round(X,digits=2);
Y=round(Y,digits=2);
lX=length(X);
DATA=matrix(nrow=2,ncol=lX);
DATA[1,]=X;
DATA[2,]=Y;
sortDATA=DATA[ ,order(DATA[1,],decreasing = FALSE)];
globalInputVec=sortDATA[1,];
globalOutputVec=sortDATA[2,];

##real data
IMPORT=read.csv(file.choose());
X=IMPORT[[7]][1:26];
Y=IMPORT[[19]][1:26];
XY=matrix(0,2,26);XY[1,]=X;XY[2,]=Y;XY=XY[,order(XY[1,], decreasing = F)];XY;
lX=length(X);
DATA=matrix(nrow=2,ncol=lX);
DATA[1,]=X;
DATA[2,]=Y;
sortDATA=DATA[ ,order(DATA[1,],decreasing = FALSE)];
globalInputVec=sortDATA[1,];
globalOutputVec=sortDATA[2,];

##setting r and numInt
radius=max(diff(globalInputVec));
numberIntervals=10^4;

####setting \tau and \rho
##CASE 1
tau1=function(x,y)     
{
return(y-x);
}

rho1=function(x,y,z)
{
v=tau1(x,y)+tau1(y,z);
return(v);
}


##CASE 2
tau2=function(x,y)     
{
return(y-x);
}

rho2=function(x,y,z)
{
v=1+0*x*y*z;
return(v);
}

##CASE 3
tau3=function(x,y)     
{
return(y-x);
}

rho3=function(x,y,z)
{
v=(tau3(x,y)+tau3(y,z))/(y-x);
return(v);
}


##CASE 4
tau4=function(x,y)     
{
return(abs(y-x)^(1/2));
}

rho4=function(x,y,z)
{
v=(tau4(x,y)+tau4(y,z))/(3*abs(z-y)+2*abs(y-x));
return(v);
}


##CASE 5
##tau=function(x,y)     
##{
##return(abs(y^2-x^2));
##}

##rho=function(x,y,z)
##{
##v=(tau(x,y)+tau(y,z))/(3*abs(z-y)+2*abs(y-x));
##return(v);
##}
###Kernel=UniKer;
Kernel=UniKer;
alphaVec=seq(0.01,0.99,by=0.01);
na=length(alphaVec);
mseVec1=vector(,na);
for(k in 1:na)
{
mse=Mse(globalInputVec,globalOutputVec,Kernel,radius,alphaVec[k],tau,rho,numberIntervals);
mseVec1[k]=mse;  
}
which(mseVec1==max(mseVec1));
which(mseVec1==min(mseVec1));
max(mseVec1);
min(mseVec1);

###Kernel=TriKer;
Kernel=TriKer;
alphaVec=seq(0.01,0.99,by=0.01);
na=length(alphaVec);
mseVec2=vector(,na);
for(k in 1:na)
{
mse=Mse(globalInputVec,globalOutputVec,Kernel,radius,alphaVec[k],tau,rho,numberIntervals);
mseVec2[k]=mse;  
}
which(mseVec2==max(mseVec2));
which(mseVec2==min(mseVec2));
max(mseVec2);
min(mseVec2);

##Kernel=Biwt;
Kernel=Biwt;
alphaVec=seq(0.01,0.99,by=0.01);
na=length(alphaVec);
mseVec3=vector(,na);
for(k in 1:na)
{
mse=Mse(globalInputVec,globalOutputVec,Kernel,radius,alphaVec[k],tau,rho,numberIntervals);
mseVec3[k]=mse;  
}
which(mseVec3==max(mseVec3));
which(mseVec3==min(mseVec3));
max(mseVec3);
min(mseVec3);


##Kernel=EpaKer;
Kernel=EpaKer;
alphaVec=seq(0.01,0.99,by=0.01);
na=length(alphaVec);
mseVec4=vector(,na);
for(k in 1:na)
{
mse=Mse(globalInputVec,globalOutputVec,Kernel,radius,alphaVec[k],tau,rho,numberIntervals);
mseVec4[k]=mse;  
}
which(mseVec4==max(mseVec4));
which(mseVec4==min(mseVec4));
max(mseVec4);
min(mseVec4);

##Kernel=Tricube;
Kernel=Tricube;
alphaVec=seq(0.01,0.99,by=0.01);
na=length(alphaVec);
mseVec5=vector(,na);
for(k in 1:na)
{
mse=Mse(globalInputVec,globalOutputVec,Kernel,radius,alphaVec[k],tau,rho,numberIntervals);
mseVec5[k]=mse;  
}
which(mseVec5==max(mseVec5));
which(mseVec5==min(mseVec5));
max(mseVec5);
min(mseVec5);

##Kernel=CosineKer;
Kernel=CosineKer;
alphaVec=seq(0.01,0.99,by=0.01);
na=length(alphaVec);
mseVec6=vector(,na);
for(k in 1:na)
{
mse=Mse(globalInputVec,globalOutputVec,Kernel,radius,alphaVec[k],tau,rho,numberIntervals);
mseVec6[k]=mse;  
}
which(mseVec6==max(mseVec6));
which(mseVec6==min(mseVec6));
max(mseVec6);
min(mseVec6);

##plotting
par(mfrow=c(3,2));
plot(alphaVec,mseVec1,ylim=c(min(mseVec1),max(mseVec1)),type="l",xlab="alpha",ylab="MSE",main="Uniform");
plot(alphaVec,mseVec2,ylim=c(min(mseVec2),max(mseVec2)),type="l",xlab="alpha",ylab="MSE",main="Triangular");
plot(alphaVec,mseVec3,ylim=c(min(mseVec3),max(mseVec3)),type="l",xlab="alpha",ylab="MSE",main="Biweight");
plot(alphaVec,mseVec4,ylim=c(min(mseVec4),max(mseVec4)),type="l",xlab="alpha",ylab="MSE",main="Epanechnikov");
plot(alphaVec,mseVec5,ylim=c(min(mseVec5),max(mseVec5)),type="l",xlab="alpha",ylab="MSE",main="Tricube");
plot(alphaVec,mseVec6,ylim=c(min(mseVec6),max(mseVec6)),type="l",xlab="alpha",ylab="MSE",main="Cosine");



plot(alphaVec,mseVec1,ylim=c(0,10000),type="l",xlab="alpha",ylab="MSE",main="CASE 2");
lines(alphaVec,mseVec2);
lines(alphaVec,mseVec3);
lines(alphaVec,mseVec4);
lines(alphaVec,mseVec5);
lines(alphaVec,mseVec6);

plot(alphaVec,mseVec1,ylim=c(10,11000),type="l",xlab="alpha",ylab="MSE",main="CASE 3");
lines(alphaVec,mseVec2);
lines(alphaVec,mseVec3);
lines(alphaVec,mseVec4);
lines(alphaVec,mseVec5);
lines(alphaVec,mseVec6);

plot(alphaVec,mseVec1,ylim=c(34,38.5),type="l",xlab="alpha",ylab="MSE",main="CASE 4");
lines(alphaVec,mseVec2);
lines(alphaVec,mseVec3);
lines(alphaVec,mseVec4);
lines(alphaVec,mseVec5);
lines(alphaVec,mseVec6);


text(0.8,34.6,"Uniform");
text(0.05,37.9,"Triangular");
text(0.05,38.3,"Biweight");
text(0.15,36.9,"Epanechnikov");
text(0.15,38.1,"Tricube");
text(0.3,37,"Cosine");




##presentations for latex 1
#real data


alpha1=0.01;
alpah2=0.25;
alpah3=0.50;
alpha4=0.75;
alpha5=0.99;

Kernel=UniKer;
yvalues_a1=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha1,tau,rho,numberIntervals);
yvalues_a2=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah2,tau,rho,numberIntervals);
yvalues_a3=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah3,tau,rho,numberIntervals);
yvalues_a4=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha4,tau,rho,numberIntervals);
yvalues_a5=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha5,tau,rho,numberIntervals);
plot(VecCentre(globalInputVec,numberIntervals),yvalues_a1,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)),ylab="outdoor death rate",xlab="index for CO",main="Uniform Kernel under Case 4");
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a2,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a3,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a4,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a5,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
points(globalInputVec,globalOutputVec);


Kernel=TriKer;
yvalues_a1=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha1,tau,rho,numberIntervals);
yvalues_a2=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah2,tau,rho,numberIntervals);
yvalues_a3=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah3,tau,rho,numberIntervals);
yvalues_a4=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha4,tau,rho,numberIntervals);
yvalues_a5=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha5,tau,rho,numberIntervals);
plot(VecCentre(globalInputVec,numberIntervals),yvalues_a1,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)),ylab="outdoor death rate",xlab="index for CO",main="Triangular Kernel under Case 4");
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a2,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a3,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a4,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a5,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
points(globalInputVec,globalOutputVec);

Kernel=Biwt;
yvalues_a1=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha1,tau,rho,numberIntervals);
yvalues_a2=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah2,tau,rho,numberIntervals);
yvalues_a3=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah3,tau,rho,numberIntervals);
yvalues_a4=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha4,tau,rho,numberIntervals);
yvalues_a5=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha5,tau,rho,numberIntervals);
plot(VecCentre(globalInputVec,numberIntervals),yvalues_a1,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)),ylab="outdoor death rate",xlab="index for CO",main="Biweight Kernel under Case 4");
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a2,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a3,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a4,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a5,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
points(globalInputVec,globalOutputVec);


Kernel=EpaKer;
yvalues_a1=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha1,tau,rho,numberIntervals);
yvalues_a2=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah2,tau,rho,numberIntervals);
yvalues_a3=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah3,tau,rho,numberIntervals);
yvalues_a4=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha4,tau,rho,numberIntervals);
yvalues_a5=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha5,tau,rho,numberIntervals);
plot(VecCentre(globalInputVec,numberIntervals),yvalues_a1,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)),ylab="outdoor death rate",xlab="index for CO",main="Epanechnikov Kernel under Case 4");
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a2,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a3,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a4,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a5,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
points(globalInputVec,globalOutputVec);

Kernel=Tricube;
yvalues_a1=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha1,tau,rho,numberIntervals);
yvalues_a2=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah2,tau,rho,numberIntervals);
yvalues_a3=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah3,tau,rho,numberIntervals);
yvalues_a4=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha4,tau,rho,numberIntervals);
yvalues_a5=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha5,tau,rho,numberIntervals);
plot(VecCentre(globalInputVec,numberIntervals),yvalues_a1,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)),ylab="outdoor death rate",xlab="index for CO",main="Tricube Kernel under Case 4");
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a2,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a3,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a4,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a5,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
points(globalInputVec,globalOutputVec);

Kernel=CosineKer;
yvalues_a1=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha1,tau,rho,numberIntervals);
yvalues_a2=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah2,tau,rho,numberIntervals);
yvalues_a3=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpah3,tau,rho,numberIntervals);
yvalues_a4=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha4,tau,rho,numberIntervals);
yvalues_a5=Yvalues(globalInputVec,globalOutputVec,Kernel,radius,alpha5,tau,rho,numberIntervals);
plot(VecCentre(globalInputVec,numberIntervals),yvalues_a1,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)),ylab="outdoor death rate",xlab="index for CO",main="Cosine Kernel under Case 4");
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a2,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a3,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a4,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
lines(VecCentre(globalInputVec,numberIntervals),yvalues_a5,type="l",lty=1,ylim=c(min(globalOutputVec),max(globalOutputVec)));
points(globalInputVec,globalOutputVec);


###### presentations for latex 2
##CASE 1
plot(alphaVec,mseVec1,xlim=c(0,1),ylim=c(14,25),type="l",xlab="alpha",ylab="MSE",main="CASE 1");
lines(alphaVec,mseVec2,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
lines(alphaVec,mseVec3,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
lines(alphaVec,mseVec4,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
lines(alphaVec,mseVec5,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
lines(alphaVec,mseVec6,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
text(0.9,15,"Uniform");
text(0.8,20.5,"Triangular");
text(0.2,24,"Biweight");
text(0.2,20.4,"Epanechnikov");
text(0.06,23,"Tricube");
text(0.4,20.8,"Cosine");

##CASE 2
plot(alphaVec,mseVec1,xlim=c(0,1),ylim=c(14,25),type="l",xlab="alpha",ylab="MSE",main="CASE 2");
lines(alphaVec,mseVec2,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
lines(alphaVec,mseVec3,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
lines(alphaVec,mseVec4,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
lines(alphaVec,mseVec5,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
lines(alphaVec,mseVec6,xlim=c(0,1),ylim=c(0,25),xlab="alpha",ylab="MSE");
text(0.9,15,"Uniform");
text(0.8,20.5,"Triangular");
text(0.2,24,"Biweight");
text(0.2,20.4,"Epanechnikov");
text(0.06,23,"Tricube");
text(0.4,20.8,"Cosine");

####============cross validation method

IndTracker=function(X,X_partial)
{
indtrackerVec=(X_partial-min(X))/0.01+1;
return(indtrackerVec);
}


m=650; ##number of sample points
n=4; ##number of testing points 
TT=combn(26,n);
SampleInd=sample(1:dim(TT)[2],m,replace=F);
TT=TT[,SampleInd];
sizeSamplesSpace=dim(TT)[2];sizeSamplesSpace;
ERR_UniKer=ERR_TriKer=ERR_Biwt=ERR_EpaKer=ERR_Tricube=ERR_CosineKer=
rep(0,sizeSamplesSpace);



for(s in 1:sizeSamplesSpace)
{
Input_testing_ind=TT[,s];
Input_training_ind=setdiff(1:26,TT[,s]);

globalInputVec=XY[1,][Input_training_ind];
testingInputVec=XY[1,][Input_testing_ind];

globalOutputVec=XY[2,][Input_training_ind];
testingOutputVec=XY[2,][Input_testing_ind];

radius=max(diff(globalInputVec))-min(diff(globalInputVec));
alpha=0.50;
numberIntervals=10000;

Kernel_UniKer=UniKer;
Kernel_TriKer=TriKer;
Kernel_Biwt=Biwt;
Kernel_EpaKer=EpaKer;
Kernel_Tricube=Tricube;
Kernel_CosineKer=CosineKer;

y_values_UniKer=Yvalues(globalInputVec,globalOutputVec,Kernel_UniKer,radius,alpha,tau,rho,numberIntervals)
y_values_TriKer=Yvalues(globalInputVec,globalOutputVec,Kernel_TriKer,radius,alpha,tau,rho,numberIntervals)
y_values_Biwt=Yvalues(globalInputVec,globalOutputVec,Kernel_Biwt,radius,alpha,tau,rho,numberIntervals)
y_values_EpaKer=Yvalues(globalInputVec,globalOutputVec,Kernel_EpaKer,radius,alpha,tau,rho,numberIntervals)
y_values_Tricube=Yvalues(globalInputVec,globalOutputVec,Kernel_Tricube,radius,alpha,tau,rho,numberIntervals)
y_values_CosineKer=Yvalues(globalInputVec,globalOutputVec,Kernel_CosineKer,radius,alpha,tau,rho,numberIntervals)

y_hat_testing_UniKer=y_values_UniKer[Input_testing_ind];
y_hat_testing_TriKer=y_values_TriKer[Input_testing_ind];
y_hat_testing_Biwt=y_values_Biwt[Input_testing_ind];
y_hat_testing_EpaKer=y_values_EpaKer[Input_testing_ind];
y_hat_testing_Tricube=y_values_Tricube[Input_testing_ind];
y_hat_testing_CosineKer=y_values_CosineKer[Input_testing_ind];

error_UniKer=(y_hat_testing_UniKer-testingOutputVec)%*%(y_hat_testing_UniKer-testingOutputVec);
error_TriKer=(y_hat_testing_TriKer-testingOutputVec)%*%(y_hat_testing_TriKer-testingOutputVec);
error_Biwt=(y_hat_testing_Biwt-testingOutputVec)%*%(y_hat_testing_Biwt-testingOutputVec);
error_EpaKer=(y_hat_testing_EpaKer-testingOutputVec)%*%(y_hat_testing_EpaKer-testingOutputVec);
error_Tricube=(y_hat_testing_Tricube-testingOutputVec)%*%(y_hat_testing_Tricube-testingOutputVec);
error_CosineKer=(y_hat_testing_CosineKer-testingOutputVec)%*%(y_hat_testing_CosineKer-testingOutputVec);

ERR_UniKer[s]=error_UniKer;
ERR_TriKer[s]=error_TriKer;
ERR_Biwt[s]=error_Biwt;
ERR_EpaKer[s]=error_EpaKer;
ERR_Tricube[s]=error_Tricube;
ERR_CosineKer[s]=error_CosineKer;
}

plot(sqrt(ERR_UniKer),xlab="k'th sampled testing set",ylab="squared errors",main="UniKer");
plot(sqrt(ERR_TriKer),xlab="k'th sampled testing set",ylab="squared errors",main="TriKer");
plot(sqrt(ERR_Biwt),xlab="k'th sampled testing set",ylab="squared errors",main="Biwt");
plot(sqrt(ERR_EpaKer),xlab="k'th sampled testing set",ylab="squared errors",main="EpaKer");
plot(sqrt(ERR_Tricube),xlab="k'th sampled testing set",ylab="squared errors",main="Tricube");
plot(sqrt(ERR_CosineKer),xlab="k'th sampled testing set",ylab="squared errors",main="CosineKer");


mean(sqrt(ERR_UniKer));sqrt(var(sqrt(ERR_UniKer)));
mean(sqrt(ERR_TriKer));sqrt(var(sqrt(ERR_TriKer)));
mean(sqrt(ERR_Biwt));sqrt(var(sqrt(ERR_Biwt)));
mean(sqrt(ERR_EpaKer));sqrt(var(sqrt(ERR_EpaKer)));
mean(sqrt(ERR_Tricube));sqrt(var(sqrt(ERR_Tricube)));
mean(sqrt(ERR_CosineKer));sqrt(var(sqrt(ERR_CosineKer)));



==============stratified sampling (size 50)====
straSamp=rep(0,50);
for(i in 1:50)
{
space=((i-1)*299+1):(i*299);
straSamp[i]=sample(space,1);
}



TT=combn(26,4);
SampleInd=straSamp;
TT=TT[,SampleInd];
sizeSamplesSpace=dim(TT)[2];sizeSamplesSpace;
ERR_UniKer=ERR_TriKer=ERR_Biwt=ERR_EpaKer=ERR_Tricube=ERR_CosineKer=
rep(0,sizeSamplesSpace);
tau=tau4;
rho=rho4;

for(s in 1:sizeSamplesSpace)
{
Input_testing_ind=TT[,s];
Input_training_ind=setdiff(1:26,TT[,s]);

globalInputVec=XY[1,][Input_training_ind];
testingInputVec=XY[1,][Input_testing_ind];

globalOutputVec=XY[2,][Input_training_ind];
testingOutputVec=XY[2,][Input_testing_ind];

radius=max(diff(globalInputVec))-min(diff(globalInputVec));
alpha=0.01;
numberIntervals=10000;

Kernel_UniKer=UniKer;
Kernel_TriKer=TriKer;
Kernel_Biwt=Biwt;
Kernel_EpaKer=EpaKer;
Kernel_Tricube=Tricube;
Kernel_CosineKer=CosineKer;

y_values_UniKer=Yvalues(globalInputVec,globalOutputVec,Kernel_UniKer,radius,alpha,tau,rho,numberIntervals)
y_values_TriKer=Yvalues(globalInputVec,globalOutputVec,Kernel_TriKer,radius,alpha,tau,rho,numberIntervals)
y_values_Biwt=Yvalues(globalInputVec,globalOutputVec,Kernel_Biwt,radius,alpha,tau,rho,numberIntervals)
y_values_EpaKer=Yvalues(globalInputVec,globalOutputVec,Kernel_EpaKer,radius,alpha,tau,rho,numberIntervals)
y_values_Tricube=Yvalues(globalInputVec,globalOutputVec,Kernel_Tricube,radius,alpha,tau,rho,numberIntervals)
y_values_CosineKer=Yvalues(globalInputVec,globalOutputVec,Kernel_CosineKer,radius,alpha,tau,rho,numberIntervals)


y_hat_testing_UniKer=y_values_UniKer[Input_testing_ind];
y_hat_testing_TriKer=y_values_TriKer[Input_testing_ind];
y_hat_testing_Biwt=y_values_Biwt[Input_testing_ind];
y_hat_testing_EpaKer=y_values_EpaKer[Input_testing_ind];
y_hat_testing_Tricube=y_values_Tricube[Input_testing_ind];
y_hat_testing_CosineKer=y_values_CosineKer[Input_testing_ind];

error_UniKer=(y_hat_testing_UniKer-testingOutputVec)%*%(y_hat_testing_UniKer-testingOutputVec);
error_TriKer=(y_hat_testing_TriKer-testingOutputVec)%*%(y_hat_testing_TriKer-testingOutputVec);
error_Biwt=(y_hat_testing_Biwt-testingOutputVec)%*%(y_hat_testing_Biwt-testingOutputVec);
error_EpaKer=(y_hat_testing_EpaKer-testingOutputVec)%*%(y_hat_testing_EpaKer-testingOutputVec);
error_Tricube=(y_hat_testing_Tricube-testingOutputVec)%*%(y_hat_testing_Tricube-testingOutputVec);
error_CosineKer=(y_hat_testing_CosineKer-testingOutputVec)%*%(y_hat_testing_CosineKer-testingOutputVec);

ERR_UniKer[s]=error_UniKer;
ERR_TriKer[s]=error_TriKer;
ERR_Biwt[s]=error_Biwt;
ERR_EpaKer[s]=error_EpaKer;
ERR_Tricube[s]=error_Tricube;
ERR_CosineKer[s]=error_CosineKer;
}

mean(sqrt(ERR_UniKer));sqrt(var(sqrt(ERR_UniKer)));
mean(sqrt(ERR_TriKer));sqrt(var(sqrt(ERR_TriKer)));
mean(sqrt(ERR_Biwt));sqrt(var(sqrt(ERR_Biwt)));
mean(sqrt(ERR_EpaKer));sqrt(var(sqrt(ERR_EpaKer)));
mean(sqrt(ERR_Tricube));sqrt(var(sqrt(ERR_Tricube)));
mean(sqrt(ERR_CosineKer));sqrt(var(sqrt(ERR_CosineKer)));

###=============multivariate input scenario=============
IMPORT=read.csv(file.choose());
CO=IMPORT[[7]][1:26];
NO=IMPORT[[13]][1:26];
DR=IMPORT[[19]][1:26];  

CODR=matrix(0,2,26);CODR[1,]=CO;CODR[2,]=DR;CODR=CODR[,order(CODR[1,], decreasing = F)];CODR;
NODR=matrix(0,2,26);NODR[1,]=NO;NODR[2,]=DR;NODR=NODR[,order(NODR[1,], decreasing = F)];NODR;

m=100; ##number of sample points
n=4; ##number of testing points 
TT=combn(26,n);
SampleInd=sample(1:dim(TT)[2],m,replace=F);
TT=TT[,SampleInd];
sizeSamplesSpace=dim(TT)[2];sizeSamplesSpace;
ERR_UniKer=ERR_TriKer=ERR_Biwt=ERR_EpaKer=ERR_Tricube=ERR_CosineKer=
rep(0,sizeSamplesSpace);

y_hat_testing_UniKer_CO=
y_hat_testing_TriKer_CO=
y_hat_testing_Biwt_CO=
y_hat_testing_EpaKer_CO=
y_hat_testing_Tricube_CO=
y_hat_testing_CosineKer_CO=
matrix(0,m,n);

y_hat_testing_UniKer_NO=
y_hat_testing_TriKer_NO=
y_hat_testing_Biwt_NO=
y_hat_testing_EpaKer_NO=
y_hat_testing_Tricube_NO=
y_hat_testing_CosineKer_NO=
matrix(0,m,n);


for(s in 1:sizeSamplesSpace)
{
Input_testing_ind=TT[,s];
Input_training_ind=setdiff(1:26,TT[,s]);

globalInputVec_CO=CODR[1,][Input_training_ind];
testingInputVec_CO=CODR[1,][Input_testing_ind];


globalInputVec_NO=NODR[1,][Input_training_ind];
testingInputVec_NO=NODR[1,][Input_testing_ind];

globalOutputVec_CO=CODR[2,][Input_training_ind];
globalOutputVec_NO=NODR[2,][Input_training_ind];

radius_CO=max(diff(globalInputVec_CO))-min(diff(globalInputVec_CO));
radius_NO=max(diff(globalInputVec_NO))-min(diff(globalInputVec_NO));

tau=tau4;
rho=rho4;
alpha=0.50;
numberIntervals=10000;

Kernel_UniKer=UniKer;
Kernel_TriKer=TriKer;
Kernel_Biwt=Biwt;
Kernel_EpaKer=EpaKer;
Kernel_Tricube=Tricube;
Kernel_CosineKer=CosineKer;

y_values_UniKer_CO=Yvalues(globalInputVec_CO,globalOutputVec_CO,Kernel_UniKer,radius_CO,alpha,tau,rho,numberIntervals)
y_values_TriKer_CO=Yvalues(globalInputVec_CO,globalOutputVec_CO,Kernel_TriKer,radius_CO,alpha,tau,rho,numberIntervals)
y_values_Biwt_CO=Yvalues(globalInputVec_CO,globalOutputVec_CO,Kernel_Biwt,radius_CO,alpha,tau,rho,numberIntervals)
y_values_EpaKer_CO=Yvalues(globalInputVec_CO,globalOutputVec_CO,Kernel_EpaKer,radius_CO,alpha,tau,rho,numberIntervals)
y_values_Tricube_CO=Yvalues(globalInputVec_CO,globalOutputVec_CO,Kernel_Tricube,radius_CO,alpha,tau,rho,numberIntervals)
y_values_CosineKer_CO=Yvalues(globalInputVec_CO,globalOutputVec_CO,Kernel_CosineKer,radius_CO,alpha,tau,rho,numberIntervals)

y_values_UniKer_NO=Yvalues(globalInputVec_NO,globalOutputVec_NO,Kernel_UniKer,radius_NO,alpha,tau,rho,numberIntervals)
y_values_TriKer_NO=Yvalues(globalInputVec_NO,globalOutputVec_NO,Kernel_TriKer,radius_NO,alpha,tau,rho,numberIntervals)
y_values_Biwt_NO=Yvalues(globalInputVec_NO,globalOutputVec_NO,Kernel_Biwt,radius_NO,alpha,tau,rho,numberIntervals)
y_values_EpaKer_NO=Yvalues(globalInputVec_NO,globalOutputVec_NO,Kernel_EpaKer,radius_NO,alpha,tau,rho,numberIntervals)
y_values_Tricube_NO=Yvalues(globalInputVec_NO,globalOutputVec_NO,Kernel_Tricube,radius_NO,alpha,tau,rho,numberIntervals)
y_values_CosineKer_NO=Yvalues(globalInputVec_NO,globalOutputVec_NO,Kernel_CosineKer,radius_NO,alpha,tau,rho,numberIntervals)

indTracker_CO=IndTracker(sort(CO),testingInputVec_CO);
indTracker_NO=IndTracker(sort(NO),testingInputVec_NO);

y_hat_testing_UniKer_CO[s,]=y_values_UniKer_CO[indTracker_CO];
y_hat_testing_TriKer_CO[s,]=y_values_TriKer_CO[indTracker_CO];
y_hat_testing_Biwt_CO[s,]=y_values_Biwt_CO[indTracker_CO];
y_hat_testing_EpaKer_CO[s,]=y_values_EpaKer_CO[indTracker_CO];
y_hat_testing_Tricube_CO[s,]=y_values_Tricube_CO[indTracker_CO];
y_hat_testing_CosineKer_CO[s,]=y_values_CosineKer_CO[indTracker_CO];

y_hat_testing_UniKer_NO[s,]=y_values_UniKer_NO[indTracker_NO];
y_hat_testing_TriKer_NO[s,]=y_values_TriKer_NO[indTracker_NO];
y_hat_testing_Biwt_NO[s,]=y_values_Biwt_NO[indTracker_NO];
y_hat_testing_EpaKer_NO[s,]=y_values_EpaKer_NO[indTracker_NO];
y_hat_testing_Tricube_NO[s,]=y_values_Tricube_NO[indTracker_NO];
y_hat_testing_CosineKer_NO[s,]=y_values_CosineKer_NO[indTracker_NO];
}


Y_testing=matrix(0,m,n);
for(s in 1:m)
{
Y_testing[s,]=sort(DR)[c(unlist(TT[,s]))];
}


Error_ML=function(theta,Y_testing,Y_hat_testing_CO,Y_hat_testing_NO)
{
Y_combined=theta*Y_hat_testing_CO+(1-theta)*Y_hat_testing_NO;
error_theta=sum((Y_testing-Y_combined)*(Y_testing-Y_combined));
return(error_theta)
}

Theta=seq(0.001,0.999,by=0.001);theta;
error_ML_Uniker=error_ML_Triker=error_ML_Biwt=error_ML_Epaker=
error_ML_Tricube=error_ML_Cosineker=rep(0,length(Theta));


for(i in 1:length(Theta))
{
error_ML_Uniker[i]=Error_ML(Theta[i],Y_testing,y_hat_testing_UniKer_CO,y_hat_testing_UniKer_NO);
error_ML_Triker[i]=Error_ML(Theta[i],Y_testing,y_hat_testing_TriKer_CO,y_hat_testing_TriKer_NO);
error_ML_Biwt[i]=Error_ML(Theta[i],Y_testing,y_hat_testing_Biwt_CO,y_hat_testing_Biwt_NO);
error_ML_Epaker[i]=Error_ML(Theta[i],Y_testing,y_hat_testing_EpaKer_CO,y_hat_testing_EpaKer_NO);
error_ML_Tricube[i]=Error_ML(Theta[i],Y_testing,y_hat_testing_Tricube_CO,y_hat_testing_Tricube_NO);
error_ML_Cosineker[i]=Error_ML(Theta[i],Y_testing,y_hat_testing_CosineKer_CO,y_hat_testing_CosineKer_NO);
}

plot(error_ML_Uniker,type="l",xlab="theta (unit 0.001)", ylab="sum of errors", main="Sum of Errors under alpha=0.5, CASE4");
lines(error_ML_Tricube,col="blue");
lines(error_ML_Epaker,col="red");
lines(error_ML_Biwt,col="green");
lines(error_ML_Triker,col="yellow");
lines(error_ML_Uniker,col="purple");

error_UniKer=(y_hat_testing_UniKer-testingOutputVec)%*%(y_hat_testing_UniKer-testingOutputVec);
error_TriKer=(y_hat_testing_TriKer-testingOutputVec)%*%(y_hat_testing_TriKer-testingOutputVec);
error_Biwt=(y_hat_testing_Biwt-testingOutputVec)%*%(y_hat_testing_Biwt-testingOutputVec);
error_EpaKer=(y_hat_testing_EpaKer-testingOutputVec)%*%(y_hat_testing_EpaKer-testingOutputVec);
error_Tricube=(y_hat_testing_Tricube-testingOutputVec)%*%(y_hat_testing_Tricube-testingOutputVec);
error_CosineKer=(y_hat_testing_CosineKer-testingOutputVec)%*%(y_hat_testing_CosineKer-testingOutputVec);

ERR_UniKer[s]=error_UniKer;
ERR_TriKer[s]=error_TriKer;
ERR_Biwt[s]=error_Biwt;
ERR_EpaKer[s]=error_EpaKer;
ERR_Tricube[s]=error_Tricube;
ERR_CosineKer[s]=error_CosineKer;




plot(sqrt(ERR_UniKer),xlab="k'th sampled testing set",ylab="squared errors",main="UniKer");
plot(sqrt(ERR_TriKer),xlab="k'th sampled testing set",ylab="squared errors",main="TriKer");
plot(sqrt(ERR_Biwt),xlab="k'th sampled testing set",ylab="squared errors",main="Biwt");
plot(sqrt(ERR_EpaKer),xlab="k'th sampled testing set",ylab="squared errors",main="EpaKer");
plot(sqrt(ERR_Tricube),xlab="k'th sampled testing set",ylab="squared errors",main="Tricube");
plot(sqrt(ERR_CosineKer),xlab="k'th sampled testing set",ylab="squared errors",main="CosineKer");


mean(sqrt(ERR_UniKer));sqrt(var(sqrt(ERR_UniKer)));
mean(sqrt(ERR_TriKer));sqrt(var(sqrt(ERR_TriKer)));
mean(sqrt(ERR_Biwt));sqrt(var(sqrt(ERR_Biwt)));
mean(sqrt(ERR_EpaKer));sqrt(var(sqrt(ERR_EpaKer)));
mean(sqrt(ERR_Tricube));sqrt(var(sqrt(ERR_Tricube)));
mean(sqrt(ERR_CosineKer));sqrt(var(sqrt(ERR_CosineKer)));





