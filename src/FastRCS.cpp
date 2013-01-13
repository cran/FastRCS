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

#include <Eigen/Dense>
#include <Eigen/LU>

using namespace std;
using namespace Eigen;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::VectorXi;

struct IdLess {					//internal function.
    template <typename T>
    IdLess(T iter) : values(&*iter) {}
    bool operator()(int left,int right){
        return values[left]<values[right];
    }
    float const* values;
};
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
		j=rand()%nn;
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
	praj=((x*FindLine(xSub,ySub,h)).array()-y.array()).array().abs2();
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
	praj=((x*FindLine(xSub,ySub,h)).array()-y.array()).array().abs2();
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
extern "C"{
	void fastrcs(int* n,int* p,int* k0,float* xi,float* yi,int* k1,float* DpC,int* nsamp,int* J,float* objfunC,int* seed,int* ck,int* ni){
		const int ik0=*k0,iJ=*J,ik1=*k1,h_m=(*n+*p+1)/2-1;
		float objfunA,objfunB=*objfunC;
		srand(*seed);
		MatrixXf x=Map<MatrixXf>(xi,*n,*p);	
		VectorXi icK=Map<VectorXi>(ck,*ni);
		VectorXf y=Map<VectorXf>(yi,*n);	
		VectorXf DpA=VectorXf::Zero(*n);
		VectorXf DpB=VectorXf::Zero(*n);

		for(int i=0;i<*nsamp;i++){			//for i=0 to i<#p-subsets.
			objfunA=Main(x,y,ik0,iJ,ik1,DpA,h_m,icK);
			if(objfunA<objfunB){
				objfunB=objfunA;
				DpB=DpA;
			}
		}

 		Map<VectorXf>(DpC,*n)=DpB.array();
		*objfunC=objfunB;
	}
}
