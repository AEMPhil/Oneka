//=============================================================================
// test_matrix.cpp
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
#include "test_matrix.h"

#include <cassert>
#include <cmath>

#include "..\Engine\matrix.h"
#include "utility.h"

namespace{
   const double TOLERANCE = 1e-9;
}

namespace oneka{

//=============================================================================
// Test Matrix Class
//
//    Test each of the class methods.  These tests are not exhaustive, but 
//    every class method is exercised and the results compared to a known
//    correct result.
//=============================================================================

//-----------------------------------------------------------------------------
// TestMatrixNullConstructor
//-----------------------------------------------------------------------------
bool TestMatrixNullConstructor()
{
   Matrix A;

   return (A.nRows() == 0 && A.nCols() == 0);
}

//-----------------------------------------------------------------------------
// TestMatrixCopyConstructor
//-----------------------------------------------------------------------------
bool TestMatrixCopyConstructor()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B( A );

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixDimensionedConstructor
//-----------------------------------------------------------------------------
bool TestMatrixDimensionedConstructor()
{
   Matrix A(2,3);
   Matrix B("0,0,0;0,0,0");

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixConstructorWithScalarFill
//-----------------------------------------------------------------------------
bool TestMatrixConstructorWithScalarFill()
{
   Matrix A(2,3, 1.2);
   Matrix B("1.2,1.2,1.2;1.2,1.2,1.2");

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixConstructorWithArrayFill
//-----------------------------------------------------------------------------
bool TestMatrixConstructorWithArrayFill()
{
   double A_data[] = {1.0,2.0,3.0,4.0,5.0,6.0};
   Matrix A(2,3, A_data);
   Matrix B("1,2,3;4,5,6");

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixConstructorWithStringFill
//-----------------------------------------------------------------------------
bool TestMatrixConstructorWithStringFill()
{
   Matrix A("1,,;4,5,");
   Matrix B("1,0,0;4,5,0");

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixDestructiveResize
//-----------------------------------------------------------------------------
bool TestMatrixDestructiveResize()
{
   Matrix A("1,2,3;4,5,6");
   A.Resize(2,2);
   Matrix B("0,0;0,0");

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixAssignmentOperator
//-----------------------------------------------------------------------------
bool TestMatrixAssignmentOperator()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B("0,1,1,0");
   B = A;

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixScalarAssignment
//-----------------------------------------------------------------------------
bool TestMatrixScalarAssignment()
{
   Matrix A("1,2,3;4,5,6");
   A = 0.0;
   Matrix B("0,0,0;0,0,0");

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixAccess
//-----------------------------------------------------------------------------
bool TestMatrixAccess()
{
   Matrix A(2,3);
   Matrix B("1,2,3;4,5,6");

   for (int i=0; i<2; ++i)
      for (int j=0; j<3; ++j)
         A(i,j) = B(i,j);

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixRowAndColumnSize
//-----------------------------------------------------------------------------
bool TestMatrixRowAndColumnSize()
{
   Matrix A("1,2,3;4,5,6");

   return (A.nRows() == 2 && A.nCols() == 3);
}

//-----------------------------------------------------------------------------
// TestMatrixAccessToRawStorage
//-----------------------------------------------------------------------------
bool TestMatrixAccessToRawStorage()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B(2,3);

   const double* pA = A.Base();
   double* pB = B.Base();

   for (int i=0; i<6; ++i)
      *pB++ = *pA++;

   return ApproxEqual(A,B,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixAccessToRawStorageWithOffset
//-----------------------------------------------------------------------------
bool TestMatrixAccessToRawStorageWithOffset()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B(2,3);

   for (int i=0; i<2; ++i)
      for (int j=0; j<3; ++j)
         *B.Base(i,j) = *A.Base(i,j);

   return ApproxEqual(A,B,TOLERANCE);
}


//=============================================================================
// TestMatrixAssociatedFunctions
//
//    Test each of the functions associated with the Matrix class.  These 
//    tests are not exhaustive, but every function is exercised and the 
//    results compared to a known correct result.
//=============================================================================

//-----------------------------------------------------------------------------
// TestMatrixColumnSum
//-----------------------------------------------------------------------------
bool TestMatrixColumnSum()
{
   Matrix A("1,2,3;4,5,6;7,8,9");
   Matrix x;
   ColumnSum( A, x );
   Matrix col_sum("12,15,18");

   return ApproxEqual(x,col_sum,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixRowSum
//-----------------------------------------------------------------------------
bool TestMatrixRowSum()
{
   Matrix A("1,2,3;4,5,6;7,8,9");
   Matrix x;
   RowSum( A, x );
   Matrix row_sum("6;15;24");

   return ApproxEqual(x,row_sum,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixTrace
//-----------------------------------------------------------------------------
bool TestMatrixTrace()
{
   Matrix A("1,2,3;4,5,6;7,8,9");

   return ApproxEqual( Trace(A), 15.0, TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixMaxAbs
//-----------------------------------------------------------------------------
bool TestMatrixMaxAbs()
{
   Matrix A("-1,2,-3;4,-5,6;-7,8,-9");

   return ApproxEqual( MaxAbs(A), 9.0, TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixL1Norm
//-----------------------------------------------------------------------------
bool TestMatrixL1Norm()
{
   Matrix A("1,2,3;4,5,6;7,8,9");

   return ApproxEqual( L1Norm(A), 18.0, TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixLInfNorm
//-----------------------------------------------------------------------------
bool TestMatrixLInfNorm()
{
   Matrix A("1,2,3;4,5,6;7,8,9");

   return ApproxEqual( LInfNorm(A), 24.0, TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixFNorm
//-----------------------------------------------------------------------------
bool TestMatrixFNorm()
{
   Matrix A("1,2,3;4,5,6;7,8,9");

   return ApproxEqual( FNorm(A), 16.8819430161341, TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixTranspose
//-----------------------------------------------------------------------------
bool TestMatrixTranspose()
{
   Matrix A("1,2,3;4,5,6;7,8,9");
   Matrix C;
   Transpose( A, C );
   Matrix At("1,4,7; 2,5,8; 3,6,9");

   return ApproxEqual(C,At,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixNegative
//-----------------------------------------------------------------------------
bool TestMatrixNegative()
{
   Matrix A("1,2,3;4,5,6;7,8,9");
   Matrix C;
   Negative( A, C );
   Matrix B("-1,-2,-3;-4,-5,-6;-7,-8,-9");

   return ApproxEqual(B,C,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixIdentity
//-----------------------------------------------------------------------------
bool TestMatrixIdentity()
{
   Matrix A(3,2,1.0);
   Identity( A, 4 );
   Matrix I("1,0,0,0; 0,1,0,0; 0,0,1,0; 0,0,0,1");

   return ApproxEqual(A,I,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixAdd_aM
//-----------------------------------------------------------------------------
bool TestMatrixAdd_aM()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B;
   Add_aM(2,A,B);
   Matrix Ap2("3,4,5;6,7,8");

   return ApproxEqual(B,Ap2,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixMultiply_aM
//-----------------------------------------------------------------------------
bool TestMatrixMultiply_aM()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B;
   Multiply_aM(2,A,B);
   Matrix Ax2("2,4,6;8,10,12");

   return ApproxEqual(B,Ax2,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixAdd_MM
//-----------------------------------------------------------------------------
bool TestMatrixAdd_MM()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B("1,0,1;0,0,1");
   Matrix C;
   Add_MM(A,B,C);
   Matrix ApB("2,2,4;4,5,7");

   return ApproxEqual(C,ApB,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixSubtract_MM
//-----------------------------------------------------------------------------
bool TestMatrixSubtract_MM()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B("1,0,1;0,0,1");
   Matrix C;
   Subtract_MM(A,B,C);
   Matrix AmB("0,2,2;4,5,5");

   return ApproxEqual(C,AmB,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixMultiply_MM
//-----------------------------------------------------------------------------
bool TestMatrixMultiply_MM()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B("1,2;3,4;5,6");
   Matrix C;
   Multiply_MM(A,B,C);
   Matrix AxB("22,28; 49,64");

   return ApproxEqual(C,AxB,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixMultiply_MtM
//-----------------------------------------------------------------------------
bool TestMatrixMultiply_MtM()
{
   Matrix A("1,4;2,5;3,6");
   Matrix B("1,2;3,4;5,6");
   Matrix C;
   Multiply_MtM(A,B,C);
   Matrix AtxB("22,28; 49,64");

   return ApproxEqual(C,AtxB,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixMultiply_MMt
//-----------------------------------------------------------------------------
bool TestMatrixMultiply_MMt()
{
   Matrix A("1,2,3;4,5,6");
   Matrix B("1,3,5;2,4,6");
   Matrix C;
   Multiply_MMt(A,B,C);
   Matrix AxBt("22,28; 49,64");

   return ApproxEqual(C,AxBt,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixMultiply_MtMt
//-----------------------------------------------------------------------------
bool TestMatrixMultiply_MtMt()
{
   Matrix A("1,4;2,5;3,6");
   Matrix B("1,3,5;2,4,6");
   Matrix C;
   Multiply_MtMt(A,B,C);
   Matrix AtxBt("22,28; 49,64");

   return ApproxEqual(C,AtxBt,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixQuadraticForm_MtMM
//-----------------------------------------------------------------------------
bool TestMatrixQuadraticForm_MtMM()
{
   Matrix a("1;2;3");
   Matrix B("1,2,3;4,5,6;7,8,9");
   Matrix c("4;5;6");

   double q = QuadraticForm_MtMM(a,B,c);

   return ApproxEqual(q,552.0,TOLERANCE);
}

//-----------------------------------------------------------------------------
// TestMatrixQuadraticForm_MMM
//-----------------------------------------------------------------------------
bool TestMatrixQuadraticForm_MMM()
{
   Matrix a("1,2,3");
   Matrix B("1,2,3;4,5,6;7,8,9");
   Matrix c("4;5;6");

   double q = QuadraticForm_MMM(a,B,c);

   return ApproxEqual(q,552.0,TOLERANCE);
}


} // namespace oneka
