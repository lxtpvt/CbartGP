/*
 *  dummybart: Dummy Bayesian Additive Regression Trees
 *  Copyright (C) 2019-2020 Xuetao Lu and Robert McCulloch
 *
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

#include "dummybart.h"

// There is only one function is different from the old bart.cpp.
// draw(const Eigen::MatrixXd& SigmaInv, rn& gen)

//--------------------------------------------------
//help functions

void dummybart::printArray(double *x, int n, std::string name )
{
  
  std::cout << name << ":";
  for(size_t k=0;k<n;k++) {
    std::cout << x[k] << " ";
  }
  std::cout << std::endl;
}

//--------------------------------------------------
//constructor
// same to the old BART 
dummybart::dummybart():m(200),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
dummybart::dummybart(size_t im):m(im),t(m),pi(),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false) {}
dummybart::dummybart(const dummybart& ib):m(ib.m),t(m),pi(ib.pi),p(0),n(0),x(0),y(0),xi(),allfit(0),r(0),ftemp(0),di(),dartOn(false)
{
   this->t = ib.t;
}
dummybart::~dummybart()
{
   if(allfit) delete[] allfit;
   if(r) delete[] r;
   if(ftemp) delete[] ftemp;
}

//--------------------------------------------------
//operators
//same to the old BART
dummybart& dummybart::operator=(const dummybart& rhs)
{
   if(&rhs != this) {

      this->t = rhs.t;
      this->m = t.size();

      this->pi = rhs.pi;

      p=0;n=0;x=0;y=0;
      xi.clear();

      if(allfit) {delete[] allfit; allfit=0;}
      if(r) {delete[] r; r=0;}
      if(ftemp) {delete[] ftemp; ftemp=0;}

   }
   return *this;
}
//--------------------------------------------------
//get,set
void dummybart::setm(size_t m)
{
   t.resize(m);
   this->m = t.size();

   if(allfit && (xi.size()==p)) predict(p,n,x,allfit);
}

//--------------------------------------------------
void dummybart::setxinfo(xinfo& _xi)
{
   size_t p=_xi.size();
   xi.resize(p);
   for(size_t i=0;i<p;i++) {
     size_t nc=_xi[i].size();
      xi[i].resize(nc);
      for(size_t j=0;j<nc;j++) xi[i][j] = _xi[i][j];
   }
}
//--------------------------------------------------
void dummybart::setdata(size_t p, size_t n, double *x, double *y, size_t numcut)
{
  int* nc = new int[p];
  for(size_t i=0; i<p; ++i) nc[i]=numcut;
  this->setdata(p, n, x, y, nc);
  delete [] nc;
}

void dummybart::setdata(size_t p, size_t n, double *x, double *y, int *nc)
{
   this->p=p; this->n=n; this->x=x; this->y=y;
   if(xi.size()==0) makexinfo(p,n,&x[0],xi,nc);

   if(allfit) delete[] allfit;
   allfit = new double[n];
   predict(p,n,x,allfit);

   if(r) delete[] r;
   r = new double[n];

   if(ftemp) delete[] ftemp;
   ftemp = new double[n];

   di.n=n; di.p=p; di.x = &x[0]; di.y=r;
   for(size_t j=0;j<p;j++){
     nv.push_back(0);
     pv.push_back(1/(double)p);
   }
}

void dummybart::setTeA(Eigen::MatrixXd& SigmaInv)
{
  for (int i = 0; i < this->m; ++i)
   {
        this->t[i].inserteA(1,1,SigmaInv.sum()); // initialize A map with (0,0) and SigmaInv.sum()
   } 
}

//--------------------------------------------------
// same to the old BART
void dummybart::predict(size_t p, size_t n, double *x, double *fp)
//uses: m,t,xi
{
   double *fptemp = new double[n];

   for(size_t j=0;j<n;j++) fp[j]=0.0;
   for(size_t j=0;j<m;j++) {
      // here, we use the fit function in treefuns, because we need prediction rather than reordering
      fit(t[j],xi,p,n,x,fptemp); 
      for(size_t k=0;k<n;k++) fp[k] += fptemp[k];
   }

   delete[] fptemp;
}
//--------------------------------------------------
// 
void dummybart::draw(const Eigen::MatrixXd& SigmaInv, rn& gen, bool sigmaNoChange)
{
  dummyutilities::index_v xid_nid;
   for(size_t j=0;j<m;j++) {
      xid_nid.clear();
      // here, we use the fit function in dummyutilities, because we need to reorder the observations
      dummyutilities::fit(t[j],xi,p,n,x,ftemp,xid_nid);
		//std::cout << "in dummybart::draw():\n" << std::endl;
		//dummyutilities::printVec(xid_nid);
      for(size_t k=0;k<n;k++) { // initialize r
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      // new birth and death function work for dependent data
      bool bd_flag;
      bd_flag = dummyutilities::bd(t[j],xi,di,pi,r,SigmaInv,nv,pv,aug,gen,xid_nid,sigmaNoChange);
      // we have to get new xid_nid after, if tree birth or death.
      if (bd_flag)
      {
      	dummyutilities::fit(t[j],xi,p,n,x,ftemp,xid_nid);
      }
      // new drmu function draws corelated bottom nodes
      dummyutilities::drmu(t[j],pi,r,SigmaInv,xid_nid,sigmaNoChange); 
      // using the fit function in bartfuns.h. fit the new or old(in new step) tree
      fit(t[j],xi,p,n,x,ftemp); 
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k]; // update allfit with new draw
   }
    //std::cout << "in dummybart::draw():" << std::endl;
    //printArray(allfit, n, "allfit:");
   if(dartOn) {
     draw_s(nv,lpv,theta,gen);
     draw_theta0(const_theta,theta,lpv,a,b,rho,gen);
     for(size_t j=0;j<p;j++) pv[j]=::exp(lpv[j]);
   }
}
//--------------------------------------------------
//public functions
void dummybart::pr() //print to screen
{
   cout << "*****dummybart object:\n";
   cout << "m: " << m << std::endl;
   cout << "t[0]:\n " << t[0] << std::endl;
   cout << "t[m-1]:\n " << t[m-1] << std::endl;
   cout << "prior and mcmc info:\n";
   pi.pr();
   if(dart){
     cout << "*****dart prior (On):\n";
     cout << "a: " << a << std::endl;
     cout << "b: " << b << std::endl;
     cout << "rho: " << rho << std::endl;
     cout << "augmentation: " << aug << std::endl;
   }
   else cout << "*****dart prior (Off):\n";
   if(p) cout << "data set: n,p: " << n << ", " << p << std::endl;
   else cout << "data not set\n";
}
