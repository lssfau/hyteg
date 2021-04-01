/*
 * Copyright (c) 2020 Daniel Drzisga
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

#include "hyteg/Operator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/composites/petsc/P2P1TaylorHoodPetsc.hpp"
#include "hyteg/edgedofspace/EdgeDoFProjectNormalOperator.hpp"
#include "hyteg/p1functionspace/P1ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

using walberla::real_t;

class P2ProjectNormalOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2ProjectNormalOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                            size_t                                                   minLevel,
                            size_t                                                   maxLevel,
                            const std::function< void( const Point3D&, Point3D& ) >& normal_function );

   ~P2ProjectNormalOperator() override = default;

   void project( const P2Function< real_t >& dst_u,
                 const P2Function< real_t >& dst_v,
                 const P2Function< real_t >& dst_w,
                 size_t                      level,
                 DoFType                     flag ) const;

   void project( const P2P1TaylorHoodFunction< real_t >& dst, size_t level, DoFType flag ) const;

#ifdef HYTEG_BUILD_WITH_PETSC
   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param numU  P2Function for determining row indices
   /// \param numV  P2Function for determining row indices
   /// \param numW  P2Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  determines on which primitives this operator is assembled
   ///
   void assembleLocalMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                             const P2Function< PetscInt >&               numU,
                             const P2Function< PetscInt >&               numV,
                             const P2Function< PetscInt >&               numW,
                             uint_t                                      level,
                             DoFType                                     flag );

   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param num   P2VectorFunction for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  determines on which primitives this operator is assembled
   ///
   void assembleLocalMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                             const P2VectorFunction< PetscInt >&         num,
                             uint_t                                      level,
                             DoFType                                     flag );
#endif

 private:
   P1ProjectNormalOperator      p1Operator;
   EdgeDoFProjectNormalOperator edgeDoFOperator;
};


} // namespace hyteg
