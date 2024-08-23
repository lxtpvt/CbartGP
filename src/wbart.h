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

#include "tree.h"
#include "treefuns.h"
#include "info.h"
#include "bartfuns.h"
#include "bd.h"
#include "bart.h"
#include "heterbart.h"


Rcpp::List cwbart(
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
    double nu,
    double lambda,
    double sigma,
    double* iw,
    bool dart,
    double theta,
    double omega,
    int *grp,
    double a,
    double b,
    double rho,
    bool aug,
    size_t nkeeptrain,
    size_t nkeeptest,
    size_t nkeeptestme,
    size_t nkeeptreedraws,
    size_t printevery,
    Rcpp::NumericMatrix xifo,
    bool showMessage
);