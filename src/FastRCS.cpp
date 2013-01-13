#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <fstream>
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <vector>
#include <random>

#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/Cholesky>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::VectorXi;
std::mt19937 mt;

struct IdLess {					//internal function.
    template <typename T>
    IdLess(T iter) : values(&*iter) {}
    bool operator()(int left,int right){
        return values[left]<values[right];
    }
    float const* values;
};
double GetUniform(){
    static std::uniform_real_distribution<double> Dist(0,1);
    return Dist(mt);
}
void GetSmallest(const VectorXf& r,int& h,const MatrixXf& x,VectorXf& y,MatrixXf& xSub,VectorXf& ySub,VectorXi& RIndex){
	const int n=x.rows();
	VectorXi SIndx2(n);
	SIndx2.setLinSpaced(n,0,n-1);
	std::nth_element(SIndx2.data(),SIndx2.data()+h,SIndx2.data()+SIndx2.size(),IdLess(r.data()));
	for (int i=0;i<h;i++){
	 	xSub.row(i)=x.row(SIndx2(i));
		ySub(i)=y(SIndx2(i));
	}
	RIndex.head(h)=SIndx2.head(h);	
}
VectorXi SampleR(const int m,const int p){
	int i,j,nn=m;
	VectorXi ind(nn);
	VectorXi y(p);
	ind.setLinSpaced(nn,0,nn-1);
    	for(i=0;i<p;i++){
		j=GetUniform()*nn;
		y(i)=ind(j);
		ind(j)=ind(--nn);
    	}
	return y;		
}
VectorXf FindLine(const MatrixXf& xSub,const VectorXf& ySub,const int h){
	const int p=xSub.cols();
	VectorXi  QIndexp(p);
	VectorXf  bt=VectorXf::Ones(p);
	QIndexp=SampleR(h,p);
	MatrixXf  A(p,p);
	for(int i=0;i<p;i++){
		A.row(i)=xSub.row(QIndexp(i));
		bt(i)=ySub(QIndexp(i));
	}
	return(A.lu().solve(bt));
}
VectorXf OneProj(const MatrixXf& x,const VectorXf& y,const MatrixXf& xSub,VectorXf& ySub,const int h,const VectorXi& RIndex,const int h_m){
	const int n=x.rows();
	VectorXf praj(n);
	praj=((x*FindLine(xSub,ySub,h)).array()-y.array()).array().abs();
	VectorXf prej(h);
	for(int i=0;i<h;i++)	prej(i)=praj(RIndex(i));
	float prem=prej.head(h).mean(),tol=1e-7;
	if(prem<tol){	
		const int n=praj.size();
		VectorXf d_resd=VectorXf::Zero(n);
		d_resd=(praj.array()<tol).select(1.0,d_resd);
		if((d_resd.sum())>=h_m){
			prem=1.0;
		} else {
			float maxin=praj.maxCoeff();
			d_resd=(praj.array()<tol).select(maxin,praj);
			prem=d_resd.minCoeff();
		}
	}
	return praj/=prem;
}
float SubsetRankFun(const MatrixXf& x,const VectorXf& y,const MatrixXf& xSub,const VectorXf& ySub,const int h,const VectorXi& RIndex){
	const int n=x.rows();
	VectorXf praj(n);
	VectorXf prej(h);
	praj=((x*FindLine(xSub,ySub,h)).array()-y.array()).array().abs();
	for(int i=0;i<h;i++)	prej(i)=praj(RIndex(i));
	nth_element(praj.data(),praj.data()+h,praj.data()+praj.size());	
	float prem=praj.head(h).mean(),fin=(prem>1e-7)?(prej.head(h).mean()/prem):(1.0);
	return fin;
}
float Main(MatrixXf& x,VectorXf& y,const int k0,const int J,const int k1,VectorXf& dP,const int h_m,VectorXi& samset){
	int p=x.cols(),n=x.rows(),h=p+1,ni=samset.size();
	MatrixXf xSub(h_m,p);
	VectorXf ySub(h_m);
	VectorXi RIndex(n);
	VectorXf fin(k1);
	VectorXi hl(J+1);
	RIndex.head(h)=SampleR(ni,h);			//draws random p-subset
	for(int i=0;i<h;i++){
		xSub.row(i)=x.row(samset(RIndex(i)));			
		ySub(i)=y(samset(RIndex(i)));
	}
	hl.setLinSpaced(J+1,h,h_m);
	h=hl(0);
	for(int j=0;j<J;j++){					//growing step
		dP=VectorXf::Zero(n);
		for(int i=0;i<k0;i++) dP+=OneProj(x,y,xSub,ySub,h,RIndex,h_m);
		h=hl(j+1);
		GetSmallest(dP,h,x,y,xSub,ySub,RIndex);
	}
	for(int i=0;i<k1;i++) fin(i)=SubsetRankFun(x,y,xSub,ySub,h,RIndex);
	return fin.array().log().mean(); 
}
VectorXi CStep(VectorXf& dP,const MatrixXf& x,const VectorXf& y,const int h){
	const int n=x.rows(),p=x.cols();
	float w1=0,w0=0;
	int w2=1,i,j=0;
	MatrixXf xSub(h,p);
	VectorXf ySub(h);
	VectorXf m_cf(p);
	MatrixXf b=MatrixXf::Identity(p,p);
	VectorXi dI(n);
	MatrixXf Sig(p,p);

	dI.setLinSpaced(n,0,n-1);
	std::nth_element(dI.data(),dI.data()+h,dI.data()+dI.size(),IdLess(dP.data()));
	for(i=0;i<h;i++) 	xSub.row(i)=x.row(dI(i));
	for(i=0;i<h;i++) 	ySub(i)=y(dI(i));
	ColPivHouseholderQR<MatrixXf> QR(xSub);
	m_cf=QR.solve(ySub);
	dP=((x*m_cf).array()-y.array()).abs();
	w1=log(((xSub*m_cf).array()-ySub.array()).abs().sum()/(float)(h-1));
	while(w2){	
		dI.setLinSpaced(n,0,n-1);
		std::nth_element(dI.data(),dI.data()+h,dI.data()+n,IdLess(dP.data()));
		for(i=0;i<h;i++) 	xSub.row(i)=x.row(dI(i));
		for(i=0;i<h;i++) 	ySub(i)=y(dI(i));
		HouseholderQR<MatrixXf> QR(xSub);
		m_cf=QR.solve(ySub);
		dP=((x*m_cf).array()-y.array()).abs();
		w0=w1;
		j++;
		w1=log(((xSub*m_cf).array()-ySub.array()).abs().sum()/(float)(h-1));
		(w0-w1<1e-3 | j>10)?(w2=0):(w2=1);
	}
	return(dI.head(h).array()+1);
} 
extern "C"{
	void fastrcs(int* n,int* p,int* k0,float* xi,float* yi,int* k1,float* DpC,int* nsamp,int* J,float* objfunC,int* seed,int* ck,int* ni,int* n1,int* n2,int* h){
		const int ik0=*k0,iJ=*J,ik1=*k1,h_m=*h;
		float objfunA,objfunB=*objfunC;
		mt.seed(*seed);
		MatrixXf x=Map<MatrixXf>(xi,*n,*p);	
		VectorXi icK=Map<VectorXi>(ck,*ni);
		VectorXf y=Map<VectorXf>(yi,*n);	
		VectorXf DpA=VectorXf::Zero(*n);
		VectorXf DpB=VectorXf::Zero(*n);
		VectorXi dpH(*n);

		for(int i=0;i<*nsamp;i++){			//for i=0 to i<#p-subsets.
			objfunA=Main(x,y,ik0,iJ,ik1,DpA,h_m,icK);
			if(objfunA<objfunB){
				objfunB=objfunA;
				DpB=DpA;
			}
		}
		Map<VectorXf>(DpC,*n)=DpB.array();
		dpH.setLinSpaced(*n,0,*n-1);
 		std::nth_element(dpH.data(),dpH.data()+h_m,dpH.data()+*n,IdLess(DpB.data()));
		Map<VectorXi>(n1,h_m)=dpH.head(h_m).array()+1;		
		Map<VectorXi>(n2,h_m)=CStep(DpB,x,y,h_m);
		*objfunC=objfunB;
	}
}
