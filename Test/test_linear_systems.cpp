//=============================================================================
// test_linear_systems.cpp
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
#include "test_linear_systems.h"

#include <cassert>
#include <cmath>

#include "..\Engine\linear_systems.h"
#include "utility.h"

namespace{
   const double TOLERANCE = 1e-9;
}

namespace oneka{

//-----------------------------------------------------------------------------
// TestCholeskyDecomposition
//-----------------------------------------------------------------------------
bool TestCholeskyDecomposition()
{
   Matrix A("4,6,4,4; 6,10,9,7; 4,9,17,11; 4,7,11,18");
   Matrix L;
   CholeskyDecomposition(A,L);
   Matrix B("2,0,0,0; 3,1,0,0; 2,3,2,0; 2,1,2,3");

   return ApproxEqual(L,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestRSPDInv
//-----------------------------------------------------------------------------
bool TestRSPDInv()
{
   Matrix A("4,6,4,4; 6,10,9,7; 4,9,17,11; 4,7,11,18");
   Matrix B("945,-690,174,-48; -690,532,-140,32; 174,-140,52,-16; -48,32,-16,16");
   Matrix Ainv;
   Multiply_aM(1.0/144.0, B, Ainv);
   RSPDInv(A,B);

   return ApproxEqual(Ainv,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestLeastSquaresSolve
//-----------------------------------------------------------------------------
bool TestLeastSquaresSolve()
{
   Matrix A("5,2,8,1; 4,6,5,5; 7,1,1,3; 2,6,1,1; 4,6,7,4; 8,6,4,2; 5,8,7,1; 7,8,2,2; 6,7,5,2; 5,5,6,2");
   Matrix B("1,7,1; 6,7,2; 3,3,2; 5,2,5; 6,5,5; 4,6,1; 5,4,8; 4,2,6; 1,8,6; 4,1,1");
   Matrix X;
   LeastSquaresSolve(A,B,X);
   Matrix C("-0.122286918422277,0.266063484829536,-0.0575443373772838; 0.464217553042304,-0.0279214573318259,0.846505417553293; -0.00883317831785533,0.470311201138176,-0.027798955351842; 0.836316520297104,0.470195843209534,-0.259472798611811");

   return ApproxEqual(X,C,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestAffineTransformation
//-----------------------------------------------------------------------------
bool TestAffineTransformation()
{
   Matrix A("7,8,6; 6,3,7; 6,1,6; 2,1,4; 1,8,8; 8,2,6; 5,5,6; 6,6,2");
   Matrix B("7,2,4; 5,1,2; 5,7,7");
   Matrix C("6,2,8");
   Matrix D;
   AffineTransformation(A,B,C,D);
   Matrix DD("125,66,94; 98,66,87; 83,57,76; 45,35,46; 93,68,84; 102,62,86; 96,59,80; 88,34,58");
   
   return ApproxEqual(D,DD,TOLERANCE);
}


} // namespace oneka