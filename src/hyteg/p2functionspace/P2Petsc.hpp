/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <hyteg/edgedofspace/EdgeDoFPetsc.hpp>
#include <hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFPetsc.hpp>
#include <hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFPetsc.hpp>
#include <hyteg/p1functionspace/P1Petsc.hpp>
#include <hyteg/p2functionspace/P2Function.hpp>
#include <hyteg/sparseassembly/SparseMatrixProxy.hpp>
#include <hyteg/sparseassembly/VectorProxy.hpp>

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {
namespace petsc {

inline void createVectorFromFunction( const P2Function< PetscScalar >&      function,
                                      const P2Function< PetscInt >&         numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   createVectorFromFunction( function.getVertexDoFFunction(), numerator.getVertexDoFFunction(), vec, level, flag );
   edgedof::createVectorFromFunction( function.getEdgeDoFFunction(), numerator.getEdgeDoFFunction(), vec, level, flag );
}

inline void createFunctionFromVector( const P2Function< PetscScalar >&      function,
                                      const P2Function< PetscInt >&         numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   createFunctionFromVector( function.getVertexDoFFunction(), numerator.getVertexDoFFunction(), vec, level, flag );
   edgedof::createFunctionFromVector( function.getEdgeDoFFunction(), numerator.getEdgeDoFFunction(), vec, level, flag );
}

inline void applyDirichletBC( const P2Function< PetscInt >& numerator, std::vector< PetscInt >& mat, uint_t level )
{
   applyDirichletBC( numerator.getVertexDoFFunction(), mat, level );
   edgedof::applyDirichletBC( numerator.getEdgeDoFFunction(), mat, level );
}

template < class OperatorType >
inline void createMatrix( const OperatorType&                         opr,
                          const P2Function< PetscInt >&               src,
                          const P2Function< PetscInt >&               dst,
                          const std::shared_ptr< SparseMatrixProxy >& mat,
                          uint_t                                      level,
                          DoFType                                     flag )
{
   createMatrix( opr.getVertexToVertexOpr(), src.getVertexDoFFunction(), dst.getVertexDoFFunction(), mat, level, flag );
   EdgeDoFToVertexDoF::createMatrix(
       opr.getEdgeToVertexOpr(), src.getEdgeDoFFunction(), dst.getVertexDoFFunction(), mat, level, flag );
   VertexDoFToEdgeDoF::createMatrix(
       opr.getVertexToEdgeOpr(), src.getVertexDoFFunction(), dst.getEdgeDoFFunction(), mat, level, flag );
   edgedof::createMatrix( opr.getEdgeToEdgeOpr(), src.getEdgeDoFFunction(), dst.getEdgeDoFFunction(), mat, level, flag );
}

} // namespace petsc
} // namespace hyteg

#endif