"REA.Gibbs" <-
function(X, X0, lambda0, Y, N= 15, save.every=1500, burn.in= 100*
save.every)
 {
   #vanilla version from Tebaldi et al. 2005
# set up some indices 
GIBBS<- save.every* N + burn.in  #total number of iterations of MCMC

# coerce input data that are vectors vectors to 1 row matrices
# this would correspond to the case where just one region is being analyzed 
if( !is.matrix(X)) X<- matrix(X, nrow=1)
if( !is.matrix(X0)) X0<- matrix(X0, nrow=1)
if( !is.matrix(lambda0)) lambda0<- matrix(lambda0, nrow=1)
if( !is.matrix(Y)) Y<- matrix(Y, nrow=1)

R<-nrow( X)    #number of regions, j= 1,....R
M<-ncol( X)   #number of AOGCMs, i=1,....M

# for debuggin uncomment the next code line to  
# fix the seed so that results are the same for checking 
#DEBUG    set.seed(123)

# parameters for priors 
a<-0.01    #parameters for gamma priors
b<-0.01

beta.x<-rep(0,R)   #initializing bias of future temperatures  
                   #   (part that is a function of control temperatures) 
mu<-rep(mean(X0),R)     #initial guess for present temperature
nu<-rep(mean(X0),R)     #initial guess for future temperature
theta<-rep(1,R)    #region specific precision factor (future)
lambda<-matrix(rep(1,R),R,M)   #model specific precision 



#
#
# output matrices and arrays,  preserving the column (region) names
reg.names<- dimnames( X0)[[1]]

mu.out<- nu.out<- theta.out<- beta.out<-matrix( NA, nrow=N,ncol=R, dimnames=list( NULL, reg.names))

lambda.out <- array( NA, c( R,M,N)) 


# center temperatures around initial parameters
Xprime<-X-matrix(mu,ncol=M,nrow=R)
Yprime<-Y-matrix(nu,ncol=M,nrow=R)

Ncount<- 1
cat("Samples completed:", fill=T)
for(itera in 1:GIBBS){
#This is the Gibbs interation loop sample each kind of parameter 
#group successively conditional on the values of the others

#theta
  ss<-as.vector((lambda*((Yprime-beta.x*Xprime)^2))%*%rep(1,M))
   
  theta<-rgamma(n=R, shape=rep(a+M/2,R), rate = a+ss/2) #R long vector

#lambda
  ss<-as.vector((((Yprime-beta.x*Xprime)^2)*matrix(theta,ncol=M,nrow=R))+
                (Xprime^2))
  lambda<-matrix(rgamma(n=M*R, shape=rep(b+1,M*R), rate = b+ss/2),nrow=R,ncol=M)  #R by M matrix

  
#beta.x
  sigis<-1/(as.vector((lambda*(Xprime^2))%*%rep(1,M))*
            theta)
  minnie<-as.vector((lambda*Xprime*Yprime)%*%rep(1,M))*
    theta*sigis
  beta.x<-rnorm(n=R,mean=minnie,sd=sqrt(sigis))  #R long vector

#mu
  sigis<-1/(theta*beta.x^2*as.vector(lambda%*%rep(1,M))+
            as.vector(lambda%*%rep(1,M))+lambda0)
  minnie<-(as.vector((lambda*X)%*%rep(1,M))- beta.x*theta*
           as.vector((lambda*(Yprime-matrix(beta.x,R,M)*X))%*%rep(1,M))
             +lambda0*as.vector(X0))*sigis
  mu<-rnorm(n=R,mean=minnie,sd=sqrt(sigis))  #R long vector

#update Xprime
  Xprime<-X-matrix(mu,ncol=M,nrow=R)
  
#nu
  sigis<-1/(theta*as.vector(lambda%*%rep(1,M)))
  minnie<-theta*
    as.vector((lambda*(Y-matrix(beta.x,R,M)*Xprime))%*%rep(1,M))*sigis
  nu<-rnorm(n=R,mean=minnie,sd=sqrt(sigis))   #R long vector

# update Yprime  
  Yprime<-Y-matrix(nu,ncol=M,nrow=R)

# save value at steps of save.every and if burn in period is done
  if((itera%%save.every==0)& ( itera > burn.in)){
if( Ncount%%10==0){
    cat( Ncount, " ")}
    mu.out[Ncount,] <- mu
    nu.out[Ncount,] <- nu
    theta.out[Ncount,] <- theta
    beta.out[Ncount,] <- beta.x
    lambda.out[,,Ncount] <- lambda

Ncount<- Ncount+1
  }

}

cat( " ",fill=T)
cat( "All done!", fill=T)

#
# return all sets of parameters together as a list along with the data.
# 

return( list(mu=mu.out, nu=nu.out, theta=theta.out, lambda=lambda.out,
beta.x= beta.out,X=X, X0=X0, lambda0=lambda0,Y=Y))
}













