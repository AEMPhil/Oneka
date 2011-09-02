//=============================================================================
// matrix.h
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
#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>

namespace oneka{

//=============================================================================
// Matrix
//=============================================================================
class Matrix
{
public:
   // Life cycle
   Matrix();                                          // null constructor
   Matrix( const Matrix& A );                         // copy constructor
   Matrix( int nrows, int ncols );                    // dimensioned constructor
   Matrix( int nrows, int ncols, double a );          // constructor w/ scalar fill
   Matrix( int nrows, int ncols, const double* a );   // constructor w/ array fill
   Matrix( const std::string& str );

   ~Matrix();                                         // destructor
   void Resize( int nrows, int ncols );               // destructive resize.

   // Operators
   Matrix& operator=( const Matrix& A );              // assignment operator
   Matrix& operator=( double a );                     // scalar assignment

   double& operator()( int row, int col );            // access
   double  operator()( int row, int col ) const;      // const access

   // Inquiry.
   int nRows() const;                                 // return the row size
   int nCols() const;                                 // return the column size

   // Access to the raw storage.
   const double* Base() const;                        // r/o access
   const double* Base( int row, int col ) const;      // r/o access

   double* Base();                                    // r/w access 
   double* Base( int row, int col );                  // r/w access with offset

private:
   int     m_nRows;                                   // allocated # of rows
   int     m_nCols;                                   // allocated # of columns
   double* m_Data;                                    // allocated memory
};


//=============================================================================
// IO Stream
//=============================================================================
std::ostream& operator << ( std::ostream& ostr, const Matrix& A );


//=============================================================================
// Matrix sums, measures and norms.
//=============================================================================
void ColumnSum( const Matrix& A, Matrix& x );
void RowSum( const Matrix& A, Matrix& x );
double Trace( const Matrix& A );

double MaxAbs( const Matrix& A );      // maximum absolute value 
double L1Norm( const Matrix& A );      // max column sum of abs      
double LInfNorm( const Matrix& A );    // max row sum of abs
double FNorm( const Matrix& A );       // sqrt of sum of sqrs


//=============================================================================
// Unary Matrix operations.
//=============================================================================
void Transpose( const Matrix& A, Matrix& C );                        // C = A'
void Negative(  const Matrix& A, Matrix& C );                        // C = -A
void Identity( Matrix& A, int n );                                   // A = I(n)

//=============================================================================
// scalar/Matrix arithmetic routines.
//=============================================================================
void Add_aM( double a, const Matrix& A, Matrix& C );                  // C = a+A
void Multiply_aM( double a, const Matrix& A, Matrix& C );            // C = a*A


//=============================================================================
// Matrix/Matrix addition and subtraction.
//=============================================================================
void Add_MM     ( const Matrix& A, const Matrix& B, Matrix& C );     // C = A + B
void Subtract_MM( const Matrix& A, const Matrix& B, Matrix& C );     // C = A - B


//=============================================================================
// Matrix/Matrix multiplication
//=============================================================================
void Multiply_MM  ( const Matrix& A, const Matrix& B, Matrix& C );   // C = AB
void Multiply_MtM ( const Matrix& A, const Matrix& B, Matrix& C );   // C = A'B
void Multiply_MMt ( const Matrix& A, const Matrix& B, Matrix& C );   // C = AB'
void Multiply_MtMt( const Matrix& A, const Matrix& B, Matrix& C );   // C = A'B'

//=============================================================================
// Quadratic forms
//=============================================================================
double QuadraticForm_MtMM( const Matrix& a, const Matrix& B, const Matrix& c );  // a'Bc
double QuadraticForm_MMM ( const Matrix& a, const Matrix& B, const Matrix& c );  // aBc


} // namespace onkea_bayes

//=============================================================================
#endif  // MATRIX_H
