/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017 Robert McCulloch and Rodney Sparapani
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

#include "cdbartSigmaFix.h"

// #define TRDRAW(a, b) trdraw[a][b]
// #define TEDRAW(a, b) tedraw[a][b]

void cdbartSigmaFix(
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   double* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   size_t m,		//number of trees
   int* numcut,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   double mybeta,
   double alpha,
   double tau,
   bool dart,
   double theta,
   double omega,
   double a,
   double b,
   double rho,
   bool aug,
   size_t nkeeptrain,
   size_t nkeeptest,
   size_t nkeeptestme,
   size_t nkeeptreedraws,
   size_t printevery,
   double* trmean,
   double* temean,
   Eigen::MatrixXd& SigmaInv,
   bool showMessage
)
{

   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;
   // std::vector<double*> trdraw(nkeeptrain);
   // std::vector<double*> tedraw(nkeeptest);
   // 
   // for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   // for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];
   
   //random number generation
   //arn gen(n1,n2); 
   arn gen;
   //heterbart bm(m);
   dummybart dbm(m);

   for(size_t i=0;i<n;i++) trmean[i]=0.0;
   for(size_t i=0;i<np;i++) temean[i]=0.0;
  

   //-----------------------------------------------------------

   size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1; // skip train
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;  // skip test
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   //--------------------------------------------------
   if(showMessage){
     printf("*****Into main of cbart\n");
     //print args
     printf("*****Data:\n");
     printf("data:n,p,np: %zu, %zu, %zu\n",n,p,np);
     printf("y1,yn: %lf, %lf\n",iy[0],iy[n-1]);
     printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
     if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
     printf("*****Number of Trees: %zu\n",m);
     printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
     printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
     printf("*****Prior:beta,alpha,tau: %lf,%lf,%lf\n", mybeta,alpha,tau);
     cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
          << dart << ',' << theta << ',' << omega << ',' << a << ',' 
          << b << ',' << rho << ',' << aug << endl;
     printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
            nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
     printf("*****printevery: %zu\n",printevery);
     printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);
     printf("\nMCMC\n");
   }

   //--------------------------------------------------
   dbm.setprior(alpha,mybeta,tau);
   dbm.setdata(p,n,ix,iy,numcut);
   dbm.setdart(a,b,rho,aug,dart,theta,omega);
   // very important for sigma fixed
   dbm.setTeA(SigmaInv);
   //--------------------------------------------------
   std::stringstream treess;  //string stream to write trees to  
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;
   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);

   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; //posterior mean for prediction
   if(np) { fhattest = new double[np]; }

   //--------------------------------------------------
   //mcmc

   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   size_t temecnt=0; //count test draws into posterior mean
   size_t treedrawscnt=0; //count kept bart draws
   bool keeptest,keeptestme,keeptreedraw;

   time_t tp;
   int time1 = time(&tp);
   xinfo& xi = dbm.getxinfo();

   for(size_t i=0;i<(nd+burn);i++) {
     if(showMessage){
       if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
     }
      if(i==(burn/2)&&dart) dbm.startdart(); // start dart when half burn
      //draw bart
      // (1) ===============================

      dbm.draw(SigmaInv,gen,true); // draw bart with Sigma inverse and random number generator
      // Save data
      if(i>=burn) {
         for(size_t k=0;k<n;k++) trmean[k]+=dbm.f(k);
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
           // for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=dbm.f(k);
           trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         keeptestme = nkeeptestme && (((i-burn+1) % skipteme) ==0) && np;

         if(keeptest || keeptestme) dbm.predict(p,np,ixp,fhattest);
         if(keeptest) {
           // for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k];
           tecnt+=1;
         }
         if(keeptestme) {
            for(size_t k=0;k<np;k++) temean[k]+=fhattest[k];
            temecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {

           for(size_t j=0;j<m;j++) {
             treess << dbm.gettree(j);
           }

           varcnt.push_back(dbm.getnv());
           varprb.push_back(dbm.getpv());
           
           treedrawscnt +=1;
        }
      }
   }
   int time2 = time(&tp);
   for(size_t k=0;k<n;k++) trmean[k]/=nd;
   for(size_t k=0;k<np;k++) temean[k]/=temecnt;
   if(showMessage){
     printf("time: %ds\n",time2-time1);
     printf("check counts\n");
     printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
   }
   //--------------------------------------------------
   if(fhattest) delete[] fhattest;
}
