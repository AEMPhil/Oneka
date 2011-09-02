//=============================================================================
// test_gaussian.cpp
//
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil Engineering
//    University of Minnesota
//
// version:
//    18 July 2011
//=============================================================================

//=============================================================================//
// Copyright 2011, Randal J. Barnes. 
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
//
//   1. Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright 
//      notice, this list of conditions and the following disclaimer in the 
//      documentation and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY RANDAL J BARNES ``AS IS'' AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
// EVENT SHALL RANDAL J BARNES OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF 
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//=============================================================================
#include "test_gaussian.h"

#include <cassert>
#include <cmath>
#include <iomanip>

#include "..\Engine\gaussian.h"
#include "utility.h"

namespace oneka{

//-----------------------------------------------------------------------------
// TestGaussianCDF
//-----------------------------------------------------------------------------
bool TestGaussianCDF()
{
   bool flag = true;

   const int N = 9;
   const double x[] = {-4, -3, -2, -1, 0, 1, 2, 3, 4};

   // These test values we computed using MATLAB's normcdf.
   const double y[] = {
      3.167124183312e-005,
      0.0013498980316301,
      0.0227501319481792,
      0.158655253931457,
      0.5,
      0.841344746068543,
      0.977249868051821,
      0.99865010196837,
      0.999968328758167 };

   for (int i=0; i<N; ++i)
   {
      flag &= ApproxEqual(GaussianCDF(x[i]),y[i],1e-9);
   }

   return flag;
}

//-----------------------------------------------------------------------------
// TestGaussianRNG
//-----------------------------------------------------------------------------
bool TestGaussianRNG()
{
   InitializeRNG();

   const int M = 14;                   // # of bins for Chi-square test
   const int N = 100000;

   // Compute the observed counts.
   Matrix Oi( M, 1, 0.0 );             // -inf, -3, -2.5, ... 2, 2.5, 3, inf

   for (int i=0; i<N; ++i)
   {
      double z = GaussianRNG();

      if (z<-3)
      {
         Oi(0,0) = Oi(0,0) + 1;
      }
      else if (z>3)
      {
         Oi(M-1,0) = Oi(M-1,0) + 1;
      }
      else
      {
         int j = static_cast<int>( ceil(2*(z+3)) );
         Oi(j,0) = Oi(j,0) + 1;
      }
   }

   // Compute the expected counts.
   double Ei_data[] = {
      N * 0.001349898,
      N * 0.004859767,
      N * 0.016540466,
      N * 0.044057069,
      N * 0.091848052,
      N * 0.149882284,
      N * 0.191462461,
      N * 0.191462461,
      N * 0.149882284,
      N * 0.091848052,
      N * 0.044057069,
      N * 0.016540466,
      N * 0.004859767,
      N * 0.001349898 };
   Matrix Ei(M,1,Ei_data);

   double ChiSquare = 0;
   for (int i=0; i<M; ++i)
   {
      ChiSquare += pow(Oi(i,0)-Ei(i,0), 2) / Ei(i,0);
   }

   const double X2crit = 34.528;    // X2inv(p=0.999, v=13)
   bool flag = ChiSquare <= X2crit;

   return flag;
}

//-----------------------------------------------------------------------------
// TestMVNormalRNG
//-----------------------------------------------------------------------------
bool TestMVNormalRNG()
{
   InitializeRNG();

   const int M = 100000;
   const int N = 3;

   Matrix Mu("1,2,3");
   Matrix Sigma("4,1,-1; 1,3,0; -1,0,2");
   Matrix Zero("0,0,0");

   // Generate a whole bunch of random vectors.
   Matrix X;
   MVNormalRNG( M, Mu, Sigma, X );
   
   // Check the means.
   Matrix Xbar(1,3);
   ColumnSum(X,Xbar);
   Multiply_aM( 1.0/M, Xbar, Xbar );

   Matrix zscore(1,3);
   for (int j=0; j<N; ++j)
   {
      zscore(0,j) = (Xbar(0,j)-Mu(0,j)) / sqrt(Sigma(j,j)/M);
   }

   bool flag = true;

   const double zcrit = 3.09;    // norminv(p=0.999);
   flag &= ApproxEqual(zscore,Zero,zcrit);

   // Check the covariances.
   Matrix B(X);
   for (int i=0; i<M; ++i)
   {
      for (int j=0; j<N; ++j)
      {
         B(i,j) -= Xbar(0,j);
      }
   }

   Matrix C, Cov;
   Multiply_MtM( B, B, C );
   Multiply_aM( 1.0/M, C, Cov );

   const double mcrit = 0.0595;   // computed by Monte Carlo simulation (M=100000,percentile=99.9).
   flag &= ApproxEqual(Sigma,Cov,mcrit);
   return flag;
}


} // namespace oneka