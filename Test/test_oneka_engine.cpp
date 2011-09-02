//=============================================================================
// test_oneka_engine.cpp
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
#include "test_oneka_engine.h"

#include <cassert>
#include <cmath>

#include "..\Engine\oneka_engine.h"
#include "utility.h"

namespace oneka{

//-----------------------------------------------------------------------------
bool TestEngine()
{
   bool flag = true;

   // First test.
   {
      double k;
      double H;
      double Base;
      int W;
      double* Xw;
      double* Yw; 
      double* Qw; 
      int P;
      double* Xp;
      double* Yp; 
      double* Ep; 
      double* Sp; 
      double Xo;
      double Yo;
      int nSims;

      // Setup the test case.
      k = 1;
      H = 50;
      Base = 0;

      W = 1;
      Xw = new double[W];
      Yw = new double[W];
      Qw = new double[W];

      Xw[0] = 0.0;
      Yw[0] = 0.0;
      Qw[0] = 30;

      P = 8;
      Xp = new double[P];
      Yp = new double[P];
      Ep = new double[P];
      Sp = new double[P];

      Xp[0] = 100;
      Yp[0] = 0;
      Ep[0] = 45.2103543000137;
      Sp[0] = 1;

      Xp[1] = 100;
      Yp[1] = 100;
      Ep[1] = 45.4674132751695;
      Sp[1] = 1;

      Xp[2] = 0;
      Yp[2] = 100;
      Ep[2] = 51.4397613593277;
      Sp[2] = 1;

      Xp[3] = -100;
      Yp[3] = 100;
      Ep[3] = 53.2728566993506;
      Sp[3] = 1;

      Xp[4] = -100;
      Yp[4] = 0;
      Ep[4] = 53.4397613593277;
      Sp[4] = 1;

      Xp[5] = -100;
      Yp[5] = -100;
      Ep[5] = 49.6717794118054;
      Sp[5] = 1;

      Xp[6] = 0;
      Yp[6] = -100;
      Ep[6] = 47.3706252432113;
      Sp[6] = 1;

      Xp[7] = 100;
      Yp[7] = -100;
      Ep[7] = 40.3396290257491;
      Sp[7] = 1;

      Xo = 0;
      Yo = 0;

      nSims = 1;

      // FROM A FULL ONEKA RUN
      // version: gcc 08.25.10
      // run on:  Thu Aug 26 15:10:40 2010
      //
      // Fitted Model Parameters
      // ---------------------------------------------------------------------------
      //         Average      Std Dev                    Correlations
      // ---------------------------------------------------------------------------
      // A:  -0.9989E-02   0.4145E-02       1.00
      // B:  -0.9989E-02   0.4067E-02       0.31   1.00
      // C:   0.1013E-02   0.2318E-02       0.03   0.03   1.00
      // D:  -0.1998E+01   0.1914E+00      -0.07  -0.02   0.05   1.00
      // E:   0.9984E+00   0.1927E+00       0.00   0.03  -0.12   0.03   1.00
      // F:   0.1300E+04   0.5325E+02      -0.78  -0.76  -0.03   0.02  -0.00   1.00
      // ---------------------------------------------------------------------------

      double Mu_ans[]  = { -0.9989E-02, -0.9989E-02, 0.1013E-02, -0.1998E+01, 0.9984E+00, 0.1300E+04 };
      double Std_ans[] = {  0.4145E-02,  0.4067E-02, 0.2318E-02,  0.1914E+00, 0.1927E+00, 0.5325E+02 };

      // Run the test case.
      EngineReturn S;
      S = Engine( k, H, Base, W, Xw, Yw, Qw, P, Xp, Yp, Ep, Sp, Xo, Yo, nSims );

      if( 
         fabs( S.Mu[0] - Mu_ans[0] ) > 0.0001E-02 ||
         fabs( S.Mu[1] - Mu_ans[1] ) > 0.0001E-02 ||
         fabs( S.Mu[2] - Mu_ans[2] ) > 0.0001E-02 ||
         fabs( S.Mu[3] - Mu_ans[3] ) > 0.0001E+01 ||
         fabs( S.Mu[4] - Mu_ans[4] ) > 0.0001E+00 ||
         fabs( S.Mu[5] - Mu_ans[5] ) > 0.0001E+04 ||
         fabs( sqrt(S.Cov[0][0]) - Std_ans[0] ) > 0.0001E-02 ||
         fabs( sqrt(S.Cov[1][1]) - Std_ans[1] ) > 0.0001E-02 ||
         fabs( sqrt(S.Cov[2][2]) - Std_ans[2] ) > 0.0001E-02 ||
         fabs( sqrt(S.Cov[3][3]) - Std_ans[3] ) > 0.0001E+00 ||
         fabs( sqrt(S.Cov[4][4]) - Std_ans[4] ) > 0.0001E+00 ||
         fabs( sqrt(S.Cov[5][5]) - Std_ans[5] ) > 0.0001E+02
         )
      {
         flag = false;
      }
   }

   return flag;
}

} // namespace oneka
