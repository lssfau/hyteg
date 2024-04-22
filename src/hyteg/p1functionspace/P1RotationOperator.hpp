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

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/composites//P1StokesFunction.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

using walberla::real_t;

class P1RotationOperator : public Operator< P1VectorFunction< real_t >, P1VectorFunction< real_t > >
{
 public:
   P1RotationOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                       size_t                                                   minLevel,
                       size_t                                                   maxLevel,
                       const std::function< void( const Point3D&, Point3D& ) >& normal_function );

   ~P1RotationOperator() override = default;

   void rotate( const P1Function< real_t >& dst_u,
                const P1Function< real_t >& dst_v,
                const P1Function< real_t >& dst_w,
                size_t                      level,
                DoFType                     flag,
                bool                        transpose = false ) const;

   void rotate( const P1VectorFunction< real_t >& dst, size_t level, DoFType flag, bool transpose = false ) const;

   void rotate( const P1StokesFunction< real_t >& dst, size_t level, DoFType flag, bool transpose = false ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< idx_t >&                  numU,
                  const P1Function< idx_t >&                  numV,
                  const P1Function< idx_t >&                  numW,
                  uint_t                                      level,
                  DoFType                                     flag ) const;

   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param numU  P1Function for determining row indices
   /// \param numV  P1Function for determining row indices
   /// \param numW  P1Function for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  determines on which primitives this operator is assembled
   ///
   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1VectorFunction< idx_t >&            src,
                  const P1VectorFunction< idx_t >&            dst,
                  uint_t                                      level,
                  DoFType                                     flag ) const override
   {
      WALBERLA_UNUSED( dst );
      toMatrix( mat, src[0], src[1], src[2], level, flag );
   }

 private:
   const std::function< void( const Point3D&, Point3D& ) > normal_function_;
};

} // namespace hyteg
