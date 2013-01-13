\name{FastRCS}
\alias{FastRCS}
\title{Computes the FastRCS multivariate outlyingness index.}
\usage{
  FastRCS(x,y,nsamp,alpha=0.5,seed=NULL)
}
\arguments{
  \item{x}{A numeric n (n>5*p) by p (p>1) matrix or data frame. Should 
	not contain an intercept.}
  \item{y}{A numeric nvector.}
   \item{nsamp}{a positive integer giving the number of resamples required;
          \code{"nsamp"} may not be reached if too many of the p-subsamples,
          chosen out of the observed vectors, are in a hyperplane.  If \code{"nsamp"} is
          omitted, it is calculated so that the probability of getting at least one
      uncontaminated starting point is always at least 99 percent when there are n/2 outliers.}
 \item{alpha}{numeric parameter controlling the size of the active subsets,
	 i.e., \code{"h=quanf(alpha,n,p)"}.  Allowed
          values are between 0.5 and 1 and the default is 0.5.}
\item{seed}{starting value for random generator. Default is seed=NULL}
}
\value{
\item{nsamp}{The value of nsamp used.}
\item{alpha}{The value of alpha used.}
\item{obj}{The value of the FastRCS objective function of the optimal h subset.}
\item{raw.fit}{A list with the following components:}
\itemize{
    \item best: The index of the h observation with smallest outlyingness indexes.
    \item ....: The results of an lm fit with weights set to 1 for the members of best and 0 otherwise.
}
\item{rew.fit}{A list with the following components:}
\itemize{
    \item best: The index of the J observation with outlyingness smaller than the rejection threshold.
    \item ....: The results of an lm fit with weights set to 1 for the members of best and 0 otherwise.
}
}
\description{
   Computes a fast and robust regression outlyingness index for a n by p matrix of 
	multivariate continuous regressors and a single dependent variable.
}
\examples{
## generate data
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
y0<-rnorm(n)
z<-c(rep(0,30),rep(1,70))
x0[1:30,]<-matrix(rnorm(30*p,5,1/100),nc=p)
y0[1:30]<-rnorm(30,10,1/100)
ns<-NumStarts(p=p,eps=0.4);
results<-FastRCS(x=x0,y=y0,alpha=0.5,nsamp=ns)
z[results$raw.fit$best]

#testing exact fit
set.seed(123)
n<-100
p<-3
x0<-matrix(rnorm(n*p),nc=p)
y0<-rep(1,n)
z<-c(rep(0,30),rep(1,70))
x0[1:30,]<-matrix(rnorm(30*p,5,1/100),nc=p)
y0[1:30]<-rnorm(30,10,1/100)
ns<-NumStarts(p=p,eps=0.4);
results<-FastRCS(x=x0,y=y0,alpha=0.5,nsamp=ns,seed=1)
z[results$raw.fit$best]
results$obj

#testing regression equivariance
n<-100
p<-3
x0<-matrix(rnorm(n*(p-1)),nc=p-1)
y0<-rnorm(n)
ns<-NumStarts(p=p,eps=0.4);
y1<-y0+cbind(1,x0)\%*\%rep(-1,p)
results1<-FastRCS(y=y0,x=x0,nsamp=ns,seed=1)$raw$coef
results2<-FastRCS(y=y1,x=x0,nsamp=ns,seed=1)$raw$coef
results1+rep(-1,p)
#should be the same:
results2
}
\references{
	Vakili, K. and Schmitt, E. (2014).
	Finding Regression Outliers With FastRCS.
	(http://arxiv.org/abs/1307.4834)
}
\author{
	Kaveh Vakili
}
\keyword{multivariate}
\keyword{robust}

