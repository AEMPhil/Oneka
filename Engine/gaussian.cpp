//=============================================================================
// gaussian.cpp
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
#include "gaussian.h"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>

#include "linear_systems.h"
#include "matrix.h"

namespace oneka{

//=============================================================================
// GaussianCDF
//
//    Returns the value of the Standard Normal cumulative distribution
//    function at the argument.
//
// notes:
// o  This routine is relatively slow, but the associated error is less than
//    1e-15 for all x.
//
// references:
// o  Marsaglia, George, 2004, Evaluating the Normal Distribution, Journal of
//    Statistical Software, v. 11, n. 4, June 2004.  Available on-line
//    <http://www.jstatsoft.org/v11/a05/paper>.
//=============================================================================
double GaussianCDF(double x)
{
   if (x < -8.0)
      return 0.0;
   else if (x > 8.0)
      return 1.0;
   else
   {
      long double s=x, t=0, b=x, q=x*x, i=1;
      while (s!=t)
         s = (t=s) + (b *= q/(i+=2));
      return 0.5 + s*exp(-0.5*q - 0.91893853320467274178L);
   }
}

//=============================================================================
// InitializeRNG
//
//    Initialize the random number generator. 
//=============================================================================

//-----------------------------------------------------------------------------
void InitializeRNG( int seed )
{
   srand( seed );
}

//-----------------------------------------------------------------------------
void InitializeRNG()
{
   srand( unsigned( time(NULL) ) );
}

//=============================================================================
// GaussianRNG
//
//    Generate a Gaussian (standard Normal) pseudo-random number.
//
// Return:
//    A single pseudo-random deviate from a standard Normal distribution 
//    (mean = 0, standard deviation = 1).
//
// Notes:
//
// o   This code is based upon Flannery et al. (1986, p. 203). 
//
// o   This routine computes two Gaussian pseudo-random number simultaneously,
//     using the Box-Mueller transformation.
//
// References:
//
// o   Press, W. B. Flannery, S. Teukolsky, and W. Vetterling, 1986, "Numerical
//     Recipes - The Art of Scientific Computing", Cambridge University Press, 
//     Cambridge, 818 pp., ISBN 0-521-30811-9.
//=============================================================================
double GaussianRNG()
{
   static bool Saved = false;
   static double S;
   
   if (Saved)
   {
      Saved = false;
      return S;
   }
   else
   {
      while (true)
      {
         double U1 = 2.0*rand()/RAND_MAX - 1.0;
         double U2 = 2.0*rand()/RAND_MAX - 1.0;

         double R = U1*U1 + U2*U2;
         if (R<1)
         {
            double P = sqrt( -2*log(R)/R );
            S = P*U1;
            Saved = true;
            return P*U2;
         }
      }
   }
}

//=============================================================================
// GaussianRNG
//
//    Generate a column matrix of uncorrelated Gaussian (standard Normal) 
//    pseudo-random numbers.
//
// Arguments:
//    M     # of random vectors to generate.
//
//    N     # of uncorrelated components in each random vector.
//
//    Z     on exit, an (MxN) matrix of pseudo-random standard Normal deviates.
//=============================================================================
bool GaussianRNG( int M, int N, Matrix& Z )
{
   assert( M >= 1 );
   assert( N >= 1 );

   Z.Resize(M,N);

   for (int i=0; i<M; ++i)
   {
      for (int j=0; j<N; ++j)
      {
         Z(i,j) = GaussianRNG();
      }
   }

   return true;
}

//=============================================================================
// MVNormalRNG
//
//    Generates a set of pseudo-random vectors from a multivariate normal 
//    distribution with mean vector Mu, and variance-covariance matrix Sigma.
//
// Arguments:  
//
//    M     # of random vectors to generate.
//
//    Mu    (1xN) vector of means.
//
//    Sigma (NxN) variance/covariance matrix.
//
//    X     (MxN) matrix of random deviates, on exit.
//
// Notes:
//
// o  The rows of X are independent random vectors.  Within each row, the 
//    columns of each are are correlated Normal random variates.
//
// o  Sigma is the covariance matrix and NOT the correlation matrix.  
//
// o  Sigma must be symmetric and positive definite -- degenerate cases are 
//    NOT allowed.
//=============================================================================
bool MVNormalRNG( int M, const Matrix& Mu, const Matrix& Sigma, Matrix& X )
{
   assert( Mu.nRows() == 1 );
   assert( Mu.nCols() >= 1 );
   assert( Mu.nCols() == Sigma.nCols() );
   assert( Sigma.nRows() == Sigma.nCols() );

   Matrix L, U;
   if( !CholeskyDecomposition( Sigma, L ) ) return false;
   Transpose(L,U);

   GaussianRNG(M, Mu.nCols(), X);
   AffineTransformation(X,U,Mu,X);

   return true;
}

} // namespace oneka