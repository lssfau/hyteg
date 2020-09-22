/*
 * Copyright (c) 2020 Marcus Mohr.
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

#define INSTANTIATE_EXTERNAL_DIAGONAL_ASSEMBLY
#include "DiagonalNonConstantOperator.hpp"
#undef INSTANTIATE_EXTERNAL_DIAGONAL_ASSEMBLY

// As long as we cannot use FunctionIterator< P2Function > we specialise here
namespace hyteg {
namespace workaround {

template <>
void externalDiagonalAssembly< P1Function< real_t > >( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                       const P1Function< real_t >&                 diagVals,
                                                       const P1Function< PetscInt >&               numerator,
                                                       uint_t                                      level,
                                                       DoFType                                     flag )
{
   if ( flag != All )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Input flag ignored in externalDiagonalAssembly(); using flag = All" );
   }

   FunctionIterator< P1Function< PetscInt > > idxIter( numerator, level );
   for ( auto valIter : FunctionIterator< P1Function< real_t > >( diagVals, level ) )
   {
      WALBERLA_ASSERT( valIter.isVertexDoF() );
      uint_t diagIdx = uint_c( ( *idxIter ).value() );
      mat->addValue( diagIdx, diagIdx, valIter.value() );
      idxIter++;
   }
}

template <>
void externalDiagonalAssembly< P2Function< real_t > >( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                       const P2Function< real_t >&                 diagVals,
                                                       const P2Function< PetscInt >&               numerator,
                                                       uint_t                                      level,
                                                       DoFType                                     flag )
{
   if ( flag != All )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Input flag ignored in externalDiagonalAssembly(); using flag = All" );
   }

   using vertexdof::VertexDoFFunction;

   FunctionIterator< VertexDoFFunction< PetscInt > > idxIterV( numerator.getVertexDoFFunction(), level );
   for ( auto valIter : FunctionIterator< VertexDoFFunction< real_t > >( diagVals.getVertexDoFFunction(), level ) )
   {
      WALBERLA_ASSERT( valIter.isVertexDoF() );
      uint_t diagIdx = uint_c( ( *idxIterV ).value() );
      mat->addValue( diagIdx, diagIdx, valIter.value() );
      idxIterV++;
   }

   FunctionIterator< EdgeDoFFunction< PetscInt > > idxIterE( numerator.getEdgeDoFFunction(), level );
   for ( auto valIter : FunctionIterator< EdgeDoFFunction< real_t > >( diagVals.getEdgeDoFFunction(), level ) )
   {
      WALBERLA_ASSERT( valIter.isEdgeDoF() );
      uint_t diagIdx = uint_c( ( *idxIterE ).value() );
      mat->addValue( diagIdx, diagIdx, valIter.value() );
      idxIterE++;
   }
}

} // namespace workaround
} // namespace hyteg
