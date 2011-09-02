//=============================================================================
// oneka_engine.cpp
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
#include "oneka_engine.h"

#include <cmath>

#include "gaussian.h"
#include "linear_systems.h"
#include "matrix.h"
#include "now.h"
#include "version.h"

namespace oneka{

//-----------------------------------------------------------------------------
// Engine
//
// Arguments:
//    k     hydraulic conductivity [L/T].
//    H     aquifer thickness [L].
//    Base  elevation of the aquifer base [L].
//
//    W     number of discharge specified wells [#].
//    Xw    (W x 1) array of well x-coordinates [L].
//    Yw    (W x 1) array of well y-coordinates [L].
//    Qw    (W x 1) array of well discharges [L^3/T].
//
//    P     number of piezometrers [#].
//    Xw    (P x 1) array of piezometer x-coordinates [L].
//    Yw    (P x 1) array of piezometer y-coordinates [L].
//    Ew    (P x 1) array of Expected Values of heads [L].
//    Sw    (P x 1) array of Standard Deviations of heads [L].
//
//    Xo    x-coordinate of model origin [L].
//    Yo    y-coordinate of model origin [L].
//
//    nSims number of realizations to generate [#].
//
// Returns:
//
//    struct EngineReturn
//    {
//       double Mu[6];        // conditional mean vector of the coefficients.
//       double Cov[6][6];    // conditional covariance matrix of the coefficients.
//       int nSims;           // number of simulations.
//       double** a;          // matrix of simulated coefficient vectors.
//    };
//
// Notes:
// o  The 2D return array, "a", is (nSims x 6).  Each row of "a" is an
//    equi-probable simulated realization of the six Oneka coefficents: 
//    [A,B,C,D,E,F].
//
// o  Since the 2D return array, "a", is being allocated in this function, 
//    there is risk of a memory leak.  We should consider changing this by 
//    either returning a managed class, or having the allocation be carried 
//    out by the calling routine.
//-----------------------------------------------------------------------------
EngineReturn Engine( 
   double k, 
   double H, 
   double Base,
   int W, 
   double* Xw, 
   double* Yw, 
   double* Qw, 
   int P, 
   double* Xp, 
   double* Yp, 
   double* Ep, 
   double* Sp, 
   double Xo,
   double Yo,
   int nSims )
{
	const double FOUR_PI	= 12.56637061435917295385057;

   // Initialize.
   Matrix A(P,6);
   Matrix b(P,1);

   // Setup the system of Oneka equations.
   for( int p = 0; p < P; ++p)
   {
      // Compute the mean and variance of Phi at piezometer p.
      double Avg, Std;

      double head = Ep[p] - Base;
      if( head < H )
      {
         Avg  = 0.5*k*(head*head + Sp[p]*Sp[p]);
         Std = k*head*Sp[p];
      }
      else
      {
         Avg  = k*H*(head - 0.5*H);
         Std = k*H*Sp[p];
      }

      // Compute the combined well potential at piezometer p.
      double Phiw = 0;

      for( int w = 0; w < W; ++w )
      {
         double dX = Xp[p] - Xw[w];
         double dY = Yp[p] - Yw[w];
         Phiw += Qw[w]/FOUR_PI * log( dX*dX + dY*dY );
      }

      // Fill in the p'th row of X, b, and V.
      double dX = Xp[p] - Xo;
      double dY = Yp[p] - Yo;

      A(p,0) = dX*dX / Std;
      A(p,1) = dY*dY / Std;
      A(p,2) = dX*dY / Std;
      A(p,3) = dX    / Std;
      A(p,4) = dY    / Std;
      A(p,5) = 1     / Std;

      b(p,0) = (Avg - Phiw)/Std;
   }

   // Compute the statistics.
   Matrix Mu;
   Matrix Cov(6,6);

   Multiply_MtM( A, A, Cov );
   if( !RSPDInv( Cov, Cov ) ) throw oneka::Exception_SingularSystem();

   // Compute the least squares fit.
   if( !LeastSquaresSolve( A, b, Mu ) ) throw oneka::Exception_SingularSystem();

   // Generate the realizations.
   Matrix X( nSims, 6 );
   Matrix Mut;
   Transpose(Mu,Mut);                           // The RNG requires a row not a column.
   MVNormalRNG( nSims, Mut, Cov, X );

   // Fill the return structure and be done.
   EngineReturn S;

   S.Version = EngineVersion();
   S.RunTime = Now();

   for (int i=0; i<6; ++i)
   {
      S.Mu[i] = Mu(i,0);
      for (int j=0; j<6; ++j)
      {
         S.Cov[i][j] = Cov(i,j);
      }
   }

   S.nSims = nSims;

   S.a = new double*[nSims];
   for (int i=0; i<nSims; ++i)
   {
      S.a[i] = new double[6];

      for (int j=0; j<6; ++j)
      {
         S.a[i][j] = X(i,j);
      }
   }

   return S;
}

} // namespace oneka