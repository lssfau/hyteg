/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/composites/P1StokesBlockPreconditioner.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

inline void createVectorFromFunction( const P1StokesFunction< real_t >&     function,
                                      const P1StokesFunction< idx_t >&      numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   createVectorFromFunction( function.uvw[0], numerator.uvw[0], vec, level, flag );
   createVectorFromFunction( function.uvw[1], numerator.uvw[1], vec, level, flag );
   if ( function.uvw[0].getStorage()->hasGlobalCells() )
   {
      createVectorFromFunction( function.uvw[2], numerator.uvw[2], vec, level, flag );
   }
   createVectorFromFunction( function.p, numerator.p, vec, level, flag );
}

inline void createFunctionFromVector( const P1StokesFunction< real_t >&     function,
                                      const P1StokesFunction< idx_t >&      numerator,
                                      const std::shared_ptr< VectorProxy >& vec,
                                      uint_t                                level,
                                      DoFType                               flag )
{
   createFunctionFromVector( function.uvw[0], numerator.uvw[0], vec, level, flag );
   createFunctionFromVector( function.uvw[1], numerator.uvw[1], vec, level, flag );
   if ( function.uvw[0].getStorage()->hasGlobalCells() )
   {
      createFunctionFromVector( function.uvw[2], numerator.uvw[2], vec, level, flag );
   }
   createFunctionFromVector( function.p, numerator.p, vec, level, flag );
}

inline void applyDirichletBC( const P1StokesFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   applyDirichletBC( numerator.uvw[0], mat, level );
   applyDirichletBC( numerator.uvw[1], mat, level );
   if ( numerator.uvw[0].getStorage()->hasGlobalCells() )
   {
      applyDirichletBC( numerator.uvw[2], mat, level );
   }
   //  applyDirichletBC(numerator.p, mat, level);
}

} // namespace hyteg
