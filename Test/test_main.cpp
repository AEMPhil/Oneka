//=============================================================================
// TestEngine.cpp
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
#include "stdafx.h"
#include <assert.h>
#include <iostream>

#include "test_gaussian.h"
#include "test_matrix.h"
#include "test_linear_systems.h"
#include "test_oneka_engine.h"

#include "..\Engine\now.h"
#include "..\Engine\version.h"

#define RUN_TEST( F ) (F ? true : (std::cerr << "FAILED: " << #F << std::endl).good() && false )


//-----------------------------------------------------------------------------
int _tmain(int argc, _TCHAR* argv[])
{
   using namespace oneka;
   bool flag = true;

   // Write out a head.
   std::cerr << std::endl << "Engine: ";
   std::cerr << oneka::EngineVersion() << std::endl;
   std::cerr << oneka::Now() << std::endl;

   // Test oneka::matrix class.
   flag &= RUN_TEST( TestMatrixNullConstructor() );
   flag &= RUN_TEST( TestMatrixCopyConstructor() );
   flag &= RUN_TEST( TestMatrixDimensionedConstructor() );
   flag &= RUN_TEST( TestMatrixConstructorWithScalarFill() );
   flag &= RUN_TEST( TestMatrixConstructorWithArrayFill() );
   flag &= RUN_TEST( TestMatrixConstructorWithStringFill() );
   flag &= RUN_TEST( TestMatrixDestructiveResize() );
   flag &= RUN_TEST( TestMatrixAssignmentOperator() );
   flag &= RUN_TEST( TestMatrixScalarAssignment() );
   flag &= RUN_TEST( TestMatrixAccess() );
   flag &= RUN_TEST( TestMatrixRowAndColumnSize() );
   flag &= RUN_TEST( TestMatrixAccessToRawStorage() );
   flag &= RUN_TEST( TestMatrixAccessToRawStorageWithOffset() );

   // Test oneka::Matrix associated functions.
   flag &= RUN_TEST( TestMatrixColumnSum() );
   flag &= RUN_TEST( TestMatrixRowSum() );
   flag &= RUN_TEST( TestMatrixTrace() );
   flag &= RUN_TEST( TestMatrixMaxAbs() );
   flag &= RUN_TEST( TestMatrixL1Norm() );
   flag &= RUN_TEST( TestMatrixLInfNorm() );
   flag &= RUN_TEST( TestMatrixFNorm() );
   flag &= RUN_TEST( TestMatrixTranspose() );
   flag &= RUN_TEST( TestMatrixNegative() );
   flag &= RUN_TEST( TestMatrixIdentity() );
   flag &= RUN_TEST( TestMatrixAdd_aM() );
   flag &= RUN_TEST( TestMatrixMultiply_aM() );
   flag &= RUN_TEST( TestMatrixAdd_MM() );
   flag &= RUN_TEST( TestMatrixSubtract_MM() );
   flag &= RUN_TEST( TestMatrixMultiply_MM() );
   flag &= RUN_TEST( TestMatrixMultiply_MtM() );
   flag &= RUN_TEST( TestMatrixMultiply_MMt() );
   flag &= RUN_TEST( TestMatrixMultiply_MtMt() );

   // Test oneka::linear_systems
   flag &= RUN_TEST( TestCholeskyDecomposition() );
   flag &= RUN_TEST( TestRSPDInv() );
   flag &= RUN_TEST( TestLeastSquaresSolve() );
   flag &= RUN_TEST( TestAffineTransformation() );

   // Test oneka::gaussian
   flag &= RUN_TEST( TestGaussianCDF() );
   flag &= RUN_TEST( TestGaussianRNG() );
   flag &= RUN_TEST( TestMVNormalRNG() );

   // Test oneka::oneka_engine
   flag &= RUN_TEST( TestEngine() );

   // A happy message...
   if (flag)
   {
      std::cerr << "ALL TESTS PASSED\n\n";
   }
   else
   {
      std::cerr << "\n>>> AT LEAST ONE TEST FAILED <<< \n\n";
   }

   return flag;
}