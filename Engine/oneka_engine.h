//=============================================================================
// oneka_engine.h
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
#pragma once
#ifndef ONEKA_ENGINE_H
#define ONEKA_ENGINE_H

#include <iostream>
#include <string>

namespace oneka{

//--------------------------------------------------------------------------
// oneka engine
//
// TODO:
// o  We MUST eliminate the "double**" return in the EngineReturn 
//    strucure.  As coded, this is a memory leak just waiting to happen.
//
//--------------------------------------------------------------------------
struct EngineReturn
{
   std::string Version;    // OnekaLite version.
   std::string RunTime;    // Run date and time.

   double Mu[6];           // conditional mean vector of the coefficients.
   double Cov[6][6];       // conditional covariance matrix of the coefficients.
   int nSims;              // number of simulations.
   double** a;             // 2d array of simulated coefficient vectors.
};

EngineReturn Engine( 
   double k, double H, double Base,
   int W, double* Xw, double* Yw, double* Qw, 
   int P, double* Xp, double* Yp, double* Ep, double* Sp, 
   double Xo, double Yo,
   int nSims );


//--------------------------------------------------------------------------
// Exception classes.
//--------------------------------------------------------------------------
class Exception_SingularSystem{};


} // namespace oneka

//=============================================================================
#endif  // ONEKA_ENGINE_H