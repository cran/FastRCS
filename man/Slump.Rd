\name{Slump}
\alias{Slump}
\title{Concrete Slump Test}
\description{
Concrete is a highly complex material. The slump flow of concrete is not 
only determined by the water content, but that is also influenced by other concrete 
ingredients. The observations for which  the values
 of either the slag and fly ash variable is exactly zero have 
been excluded.
}
\usage{
Slump
}
\format{
  \describe{
	\item{\code{28_compressive_strenght}}{}
	\item{\code{Cement}}{quantity of cement (kg/m3)}
	\item{\code{Slag}}{blast furnace slag (kg/m3)}	
	\item{\code{Fly.ash}}{fly ash (kg/m3)}
	\item{\code{Water}}{water (kg/m3)}
	\item{\code{SP}}{super-plasticizer (kg/m3)}
	\item{\code{Coarse.Aggr.}}{coarse aggregate (kg/m3)}
	\item{\code{Fine.Aggr.}}{fine aggregate (kg/m3)}
  }
}
\source{
Yeh, I. (2007). 
Modeling slump flow of concrete using second-order regressions and artificial neural networks. 
Cement and Concrete Composites, Vol.29, No. 6, 474--480.
}
\examples{
data(Slump)
alpha<-0.5
p<-ncol(Slump)
ns<-NumStarts(p=p,eps=(1-alpha)*4/5)
Fit<-FastRCS(x=Slump[,-1],y=Slump[,1],nsamp=ns,seed=1)
yname<-"Standardized robust residuals"
plot(abs(Fit$rew$resid)/summary(Fit$rew)$sigma,pch=16,ylab=yname)
}
\keyword{datasets}
