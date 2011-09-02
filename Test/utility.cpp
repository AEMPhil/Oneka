//=============================================================================
// utility.cpp
//    
//    Utility functions to streamline the unit testing.
//
// author:
//    Dr. Randal J. Barnes
//    Department of Civil Engineering
//    University of Minnesota
//
// version:
//    18 July 2011
//=============================================================================

//=============================================================================
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
#include "utility.h"

#include <cassert>
#include <cmath>

namespace oneka{

//-----------------------------------------------------------------------------
bool ApproxEqual( const double x, const double y, double tol )
{
   return (fabs(x-y) <= tol);
}

//-----------------------------------------------------------------------------
bool ApproxEqual( const Matrix& A, const Matrix& B, double tol )
{
   // Compare the sizes first.
   if (A.nRows() != B.nRows() || A.nCols() != B.nCols())
      return false;

   // Compare the contents.
   Matrix C;

   Subtract_MM( A, B, C );
   return (MaxAbs(C) <= tol);
}

//-----------------------------------------------------------------------------
bool RelativeEqual( const double x, const double y, double tol )
{
   return (fabs(x-y) <= tol*fabs(y));
}


} // namespace oneka
