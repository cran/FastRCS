\name{DontGetKicked}
\alias{DontGetKicked}
\title{Sales Data for the Chrysler Town \& Country}
\description{
  A data frame with the subset of the 'Don't Get Kicked!'
 dataset corresponding to sales of the Chrysler Town \& Country. 
}
\usage{
DontGetKicked
}
\format{
  \describe{
	\item{\code{VehBCost}}{Acquisition cost paid for the vehicle at time of purchase.}
	\item{\code{MMRAcquisitionAuctionAveragePrice}}{Acquisition price for this vehicle in average condition at time of purchase.}
	\item{\code{MMRAcquisitonRetailCleanPrice}}{Acquisition price for this vehicle in the above Average condition at time of purchase.}	
	\item{\code{MMRAcquisitionRetailAveragePrice}}{Acquisition price for this vehicle in the retail market in average condition at time of purchase.}
	\item{\code{MMRAcquisitonRetailCleanPrice}}{Acquisition price for this vehicle in the retail market in above average condition at time of purchase.}
	\item{\code{MMRCurrentAuctionAveragePrice}}{Acquisition price for this vehicle in average condition as of current day.}
	\item{\code{MMRCurrentAuctionCleanPrice}}{Acquisition price for this vehicle in above condition as of current day.}
	\item{\code{MMRCurrentRetailAveragePrice}}{Acquisition price for this vehicle on the retail market in average condition as of current day.}
	\item{\code{MMRCurrentRetailCleanPrice}}{Acquisition price for this vehicle on the retail market in above average condition as of current day.}
	\item{\code{WarrantyCost}}{Warranty price (term=36month  and millage=36K).}
	\item{\code{VehOdo}}{The vehicle's odometer reading.}
  }
}
\source{
Kaggle Data Challenge (2012). 
Don't Get Kicked!. 
\url{http://www.kaggle.com/c/DontGetKicked}.
}
\examples{
data(DontGetKicked)
alpha<-0.5
p<-ncol(DontGetKicked)
ns<-numStarts(p=p,eps=(1-alpha)*4/5)
Fit<-FastRCS(x=DontGetKicked[,-1],y=DontGetKicked[,1],nSamp=ns,seed=1)
plot(Fit)
}
\keyword{datasets}
