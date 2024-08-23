/*
 *  dummybart: Dummy Bayesian Additive Regression Trees
 *  Copyright (C) 2019-2020 Xuetao Lu and Robert McCulloch
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#ifndef GUARD_dummybart_h
#define GUARD_dummybart_h

#include <ctime>
#include "dummyutilities.h"

//#include "dummyutilities.h"

class dummybart {
public:
//--------------------------------------------------
//help functions

	void printArray(double *x, int n, std::string name );

   //------------------------------
   //constructor/destructor
   dummybart();
   dummybart(size_t m);
   dummybart(const dummybart&);
   ~dummybart();
   //------------------------------
   //operators
   dummybart& operator=(const dummybart&);
   //------------------------------
   //get,set
   size_t getm() {return m;}
   void setm(size_t m);
   void setdata(size_t p, size_t n, double *x, double *y, size_t nc=100);
   void setdata(size_t p, size_t n, double *x, double *y, int* nc);
   
   void setpi(pinfo& pi) {this->pi = pi;}
   
   void setprior(double alpha, double beta, double tau)
      {pi.alpha=alpha; pi.mybeta = beta; pi.tau=tau;}
  
   void setdart(double _a, double _b, double _rho, bool _aug, bool _dart, 
		double _theta=0., double _omega=1.) {
     this->a=_a; this->b=_b; this->rho=_rho; this->aug=_aug; 
     this->dart=_dart; this->omega=_omega; 
     if(_theta==0.){
       this->const_theta=false;
       this->theta=1.;
     }
     else{
       this->const_theta=true;
       this->theta=_theta;
     }
}

   void startdart() {this->dartOn=!(this->dartOn);}

   void setTeA(Eigen::MatrixXd& SigmaInv);

   void settau(double tau) {pi.tau=tau;}

   tree& gettree(size_t i ) { return t[i];}

   xinfo& getxinfo() {return xi;}
   // Add for test
   dinfo& getdinfo() {return di;}

   void setxinfo(xinfo& _xi);
   std::vector<size_t>& getnv() {return nv;}
   // Add for test
   void setnv(std::vector<size_t>& _nv) {nv=_nv;}

   std::vector<double>& getpv() {return pv;}
   double gettheta() {return theta;}
   //------------------------------
   //public methods
   void birth(size_t i, size_t nid,size_t v, size_t c, double ml, double mr)
         {t[i].birth(nid,v,c,ml,mr);}
   void death(size_t i,size_t nid, double mu)
         {t[i].death(nid,mu);}
   void pr();
   void tonull() {for(size_t i=0;i!=t.size();i++) t[i].tonull();}
   void predict(size_t p, size_t n, double *x, double *fp);
   void draw(const Eigen::MatrixXd& SigmaInv, rn& gen, bool sigmaNoChange);
//   void draw_s(rn& gen);
   double f(size_t i) {return allfit[i];}
protected:
   size_t m;  //number of trees
   std::vector<tree> t; //the trees
   pinfo pi; //prior and mcmc info
   //data
   size_t p,n; //x has dim p, n obserations
   double *x,*y;  //x is column stack, pxn
   xinfo xi; //cutpoint info
   //working
   double *allfit; //if the data is set, should be f(x)
   double *r;
   double *ftemp;
   dinfo di;
   bool dart,dartOn,aug,const_theta;
   double a,b,rho,theta,omega;
   std::vector<size_t> nv;
   std::vector<double> pv, lpv;
};

#endif
