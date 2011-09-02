//=============================================================================
// test_matrix.h
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
#ifndef TEST_MATRIX_H
#define TEST_MATRIX_H

namespace oneka{

// Test oneka::matrix class.
bool TestMatrixNullConstructor();
bool TestMatrixCopyConstructor();
bool TestMatrixDimensionedConstructor();

bool TestMatrixConstructorWithScalarFill();
bool TestMatrixConstructorWithArrayFill();
bool TestMatrixConstructorWithStringFill();

bool TestMatrixDestructiveResize();

bool TestMatrixAssignmentOperator();
bool TestMatrixScalarAssignment();

bool TestMatrixAccess();
bool TestMatrixRowAndColumnSize();

bool TestMatrixAccessToRawStorage();
bool TestMatrixAccessToRawStorageWithOffset();

// Test oneka::matrix associated functions.
bool TestMatrixColumnSum();
bool TestMatrixRowSum();
bool TestMatrixTrace();

bool TestMatrixMaxAbs();
bool TestMatrixL1Norm();
bool TestMatrixLInfNorm();
bool TestMatrixFNorm();

bool TestMatrixTranspose();
bool TestMatrixNegative();
bool TestMatrixIdentity();

bool TestMatrixAdd_aM();
bool TestMatrixMultiply_aM();

bool TestMatrixAdd_MM();
bool TestMatrixSubtract_MM();

bool TestMatrixMultiply_MM();
bool TestMatrixMultiply_MtM();
bool TestMatrixMultiply_MMt();
bool TestMatrixMultiply_MtMt();

bool TestMatrixQuadraticForm_MtMM();
bool TestMatrixQuadraticForm_MMM();

} // namespace onkea_bayes

//=============================================================================
#endif  // TEST_MATRIX_H
