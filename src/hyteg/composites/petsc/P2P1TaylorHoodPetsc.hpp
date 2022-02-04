/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesBlockPreconditioner.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesBlockPreconditioner.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p2functionspace/P2Petsc.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

inline void createVectorFromFunction( const P2P1TaylorHoodFunction< real_t >& function,
                                      const P2P1TaylorHoodFunction< idx_t >&  numerator,
                                      const std::shared_ptr< VectorProxy >&   vec,
                                      uint_t                                  level,
                                      DoFType                                 flag )
{
  function.toVector( numerator, vec, level, flag );
}

inline void createFunctionFromVector( const P2P1TaylorHoodFunction< real_t >& function,
                                      const P2P1TaylorHoodFunction< idx_t >&  numerator,
                                      const std::shared_ptr< VectorProxy >&   vec,
                                      uint_t                                  level,
                                      DoFType                                 flag )
{
  function.fromVector( numerator, vec, level, flag );
}

inline void applyDirichletBC( const P2P1TaylorHoodFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
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
