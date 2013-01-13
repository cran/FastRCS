NumStarts<-function(p,gamma=0.99,eps=0.5){
	if(p>25)	stop("p too large.")
	if(gamma>=1)	stop("gamma should be smaller than 1.")
	if(gamma<=0)	stop("gamma should be larger than 0.")
	if(eps>0.5)	stop("eps should be smaller than 1/2.")
	if(eps<=0)	stop("eps should be larger than 0.")	
	ns0<-ceiling(log(1-gamma)/log(1-(1-(eps))^(p+1)))
	ord<-10^floor(log10(ns0))
	ceiling(ns0/ord)*ord
}
FastRCS<-function(x,y,nsamp=NULL,alpha=0.5,seed=NULL){#x<-x0;y<-y0;nsamp<-ns;alpha<-0.5
	k1<-25;k0<-25;J<-3;
	if(is.null(seed))	seed<-floor(runif(1,-2^31,2^31))
	seed<-as.integer(seed)+1
	x<-data.matrix(x)
	y<-data.matrix(y)
	na.x<-complete.cases(cbind(x,y))
	if(!is.numeric(alpha))	stop("alpha should be numeric")
	if(alpha<0.5 | alpha>=1)stop("alpha should be in (0.5,1(.")
	if(sum(na.x)!=nrow(x))  stop("Your data contains NA.")
	if(nrow(x)<(5*ncol(x))) stop("n<5p. You need more observations")
	cx<-cbind(1,x)
	n<-nrow(cx)
	if(nrow(unique(cbind(x,y)))<n)	stop("Your dataset contains duplicated rows. Please remove them.") 
	p<-ncol(cx)
	if(p<2)			stop("Univariate RCS is not implemented.")
	if(p>25)		stop("FastRCS only works for dimensions <=25.")
	if(is.null(nsamp)) 	nsamp<-NumStarts(p,eps=(1-alpha)) 
	h<-quanf(n=n,p=p,alpha=alpha)
	h0<-quanf(n=n,p=p,alpha=0.5);
	Dp<-rep(1.00,n);
	k0<-max(k0,p+2);
	k1<-max(k1,p+2);
	objfunC<-1e3;
	icandid<-1:n-1
	ni<-length(icandid)
	fitd<-.C("fastrcs",as.integer(nrow(cx)),as.integer(ncol(cx)),as.integer(k0),as.single(cx),as.single(y),as.integer(k1),as.single(Dp),as.integer(nsamp),as.integer(J),as.single(objfunC),as.integer(seed),as.integer(icandid),as.integer(ni),PACKAGE="FastRCS")
	outd<-as.numeric(fitd[[7]])
	if(is.nan(outd)[1])	stop("too many singular subsets encoutered!")
	best<-which(outd<=median(outd))
	weit<-as.numeric((1:n)%in%best)
	rawF<-lm(y~x,weights=weit)
	best<-which(abs(rawF$resid)<=quantile(abs(rawF$resid),h/n))
	weit<-as.numeric((1:n)%in%best)
	rawF<-lm(y~x,weights=weit);rawF$best<-best
	weit<-as.numeric(abs(rawF$resid)/median(abs(rawF$resid))/qnorm(1-alpha/2)<=sqrt(qchisq(0.975,df=1)))
       	rewF<-lm(y~x,weights=weit);rewF$best<-which(weit==1)	
	list(alpha=alpha,nsamp=nsamp,obj=as.numeric(fitd[[10]]),raw.fit=rawF,rew.fit=rewF)
}
quanf<-function(n,p,alpha)	return(floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alpha))
