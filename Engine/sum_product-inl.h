//=============================================================================
// sum_product-inl.h
//
//    A simple implementation of a core linear algebra computational component.
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
#pragma once
#ifndef SUM_PRODUCT_H
#define SUM_PRODUCT_H

namespace oneka{

//-----------------------------------------------------------------------------
// This routine computes a dot product between two vectors.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the first vector.
//    y     pointer to the first element of the second vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, const double* y )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
      Sum += (*x++) * (*y++);

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between two vectors.  Both vectors
// allow for a non-unit stride.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the first vector.
//    dx    stride between subsequent elements in the first vector.
//    y     pointer to the first element of the second vector.
//    dy    stride between subsequent elements in the second vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, int dx, const double* y, int dy )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x) * (*y);
      x += dx;
      y += dy;
   }

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between two vectors.  The y vector
// allows for a non-unit stride.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the first vector.
//    y     pointer to the first element of the second vector.
//    dy    stride between subsequent elements in the second vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, const double* y, int dy )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x++) * (*y);
      y += dy;
   }

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between two vectors.  The x vector
// allows for a non-unit stride.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the first vector.
//    dx    stride between subsequent elements in the first vector.
//    y     pointer to the first element of the second vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, int dx, const double* y )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x) * (*y++);
      x += dx;
   }

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between a vector and itself.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x) * (*x);
      ++x;
   }

   return Sum;
}

//-----------------------------------------------------------------------------
// This routine computes a dot product between a vector and itself.  The 
// vector allows for a non-unit stride.
//
// Arguments:
//
//    n     total number of elements in each vector.
//    x     pointer to the first element of the vector.
//    dx    stride between subsequent elements in the vector.
//-----------------------------------------------------------------------------
inline double SumProduct( int n, const double* x, int dx )
{
   double Sum = 0.0;

   for (int i=0; i<n; ++i)
   {
      Sum += (*x) * (*x);
      x += dx;
   }

   return Sum;
}


} // namespace oneka

//=============================================================================
#endif  // SUM_PRODUCT_H