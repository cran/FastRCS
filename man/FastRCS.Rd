\name{FastRCS}
\alias{FastRCS}
\title{Computes the FastRCS multivariate outlyingness index.}
\usage{
  FastRCS(x,y,nsamp,alpha=0.5,seed=1)
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
\item{seed}{starting value for random generator. A positive integer. Default is seed = 1}
}
\value{
\item{nsamp}{The value of nsamp used.}
\item{alpha}{The value of alpha used.}
\item{obj}{The value of the FastRCS objective function of the optimal h subset.}
\item{raw}{A list with the following components:}
\itemize{
    \item best: The index of the h observation with smallest outlyingness indexes.
    \item coefficients: The vector of coefficients of the hyperplane fitted to the members of \code{$raw$best}.
    \item fitted.values: the fitted mean values: \code{cbind(1,x)\%*\%raw$coefficients}.
	\item residuals: the residuals, that is response minus fitted values.
    	\item rank: the numeric rank of the fitted linear model.
	\item  weights: (only for weighted fits) the specified weights.
	\item df.residual: the residual degrees of freedom.
}
\item{rew}{A list with the following components:}
\itemize{
    \item best: The index of the J observation with outlyingness smaller than the rejection threshold.
    \item coefficients: The vector of coefficients of the hyperplane fitted to the members of \code{$rew$best}.
	\item fitted.values: the fitted mean values: \code{cbind(1,x)\%*\%rew$coefficients}.
	\item residuals: the residuals, that is response minus fitted values.
    	\item rank: the numeric rank of the fitted linear model.
	\item  weights: (only for weighted fits) the specified weights.
	\item df.residual: the residual degrees of freedom.
}
\item{rawC}{Same structure as \code{$raw}, but after carrying  CSteps off of the members of \code{$raw} until convergence. See details.}
\item{rewC}{Same structure as \code{$rew}, but after carrying  CSteps off of the members of \code{$raw} until convergence and only then
 performing a re-weighting step. See details.}
}
\details{
	The slots \code{rew} and \code{raw} are provided for backward compatibility only. 
	FastRCS is regression and affine equivariant and thus consistent at the elliptical
	 model (Grubel and Rock (1990)). 
	
	Nonetheless, in experiments, we found that carrying CSteps (Rousseeuw and van Driessen (2006)) 
	off of the members of \code{$raw} improves the speed of convergence without increasing the bias
	 of the final estimates.
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
z[results$raw$best]

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
z[results$raw$best]
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
	Grubel, R. and Rocke, D. M. (1990).
	On the cumulants of affine equivariant estimators in elliptical families.
	Journal of Multivariate Analysis, Vol. 35, p. 203--222. Journal of Multivariate Analysis

	Rousseeuw, P. J., and van Driessen, K. (2006). Computing lts regression for large data sets. 
	Data mining and	Knowledge Discovery, 12, 29--45

	Vakili, K. and Schmitt, E. (2014).
	Finding Regression Outliers With FastRCS.
	(http://arxiv.org/abs/1307.4834)
}
\author{
	Kaveh Vakili
}
\keyword{multivariate}
\keyword{robust}

