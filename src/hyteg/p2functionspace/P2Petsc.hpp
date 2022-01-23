/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl
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

#include "hyteg/edgedofspace/EdgeDoFPetsc.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFPetsc.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFPetsc.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

// ============
//  P2Function
// ============
inline void createVectorFromFunction( const P2Function< real_t >&           function,
                                      const P2Function< idx_t >&            numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   createVectorFromFunction( function.getVertexDoFFunction(), numerator.getVertexDoFFunction(), vec, level, flag );
   createVectorFromFunction( function.getEdgeDoFFunction(), numerator.getEdgeDoFFunction(), vec, level, flag );
}

inline void createFunctionFromVector( const P2Function< real_t >&           function,
                                      const P2Function< idx_t >&            numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   createFunctionFromVector( function.getVertexDoFFunction(), numerator.getVertexDoFFunction(), vec, level, flag );
   createFunctionFromVector( function.getEdgeDoFFunction(), numerator.getEdgeDoFFunction(), vec, level, flag );
}

inline void applyDirichletBC( const P2Function< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   applyDirichletBC( numerator.getVertexDoFFunction(), mat, level );
   applyDirichletBC( numerator.getEdgeDoFFunction(), mat, level );
}

// ==================
//  P2VectorFunction
// ==================
inline void createVectorFromFunction( const P2VectorFunction< real_t >&     function,
                                      const P2VectorFunction< idx_t >&      numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   for ( uint_t k = 0; k < function.getDimension(); k++ )
   {
      createVectorFromFunction( function[k], numerator[k], vec, level, flag );
   }
}

inline void createFunctionFromVector( const P2VectorFunction< real_t >&     function,
                                      const P2VectorFunction< idx_t >&      numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   for ( uint_t k = 0; k < function.getDimension(); k++ )
   {
      createFunctionFromVector( function[k], numerator[k], vec, level, flag );
   }
}

inline void applyDirichletBC( const P2VectorFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   for ( uint_t k = 0; k < numerator.getDimension(); k++ )
   {
      applyDirichletBC( numerator[k], mat, level );
   }
}

} // namespace hyteg
