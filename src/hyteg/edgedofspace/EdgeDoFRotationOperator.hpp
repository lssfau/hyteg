/*
 * Copyright (c) 2023 Ponsuganth Ilangovan P
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

#include "hyteg/Algorithms.hpp"
#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFOperatorTypeDefs.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

class EdgeDoFRotationOperator final : public Operator< hyteg::EdgeDoFFunction< real_t >, hyteg::EdgeDoFFunction< real_t > >
{
 public:
   EdgeDoFRotationOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                            size_t                                                   minLevel,
                            size_t                                                   maxLevel,
                            const std::function< void( const Point3D&, Point3D& ) >& normal_function );
   ~EdgeDoFRotationOperator() = default;

   void rotate( const EdgeDoFFunction< real_t >& dst_u,
                const EdgeDoFFunction< real_t >& dst_v,
                const EdgeDoFFunction< real_t >& dst_w,
                uint_t                           level,
                DoFType                          flag,
                bool                             transpose = false ) const;

   /// Assemble operator as sparse matrix
   ///
   /// \param mat   a sparse matrix proxy
   /// \param numU  EdgeDoFFunction for determining row indices
   /// \param numV  EdgeDoFFunction for determining row indices
   /// \param numW  EdgeDoFFunction for determining row indices
   /// \param level level in mesh hierarchy for which local operator is to be assembled
   /// \param flag  determines on which primitives this operator is assembled
   ///
   void assembleLocalMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                             const EdgeDoFFunction< idx_t >&             numU,
                             const EdgeDoFFunction< idx_t >&             numV,
                             const EdgeDoFFunction< idx_t >&             numW,
                             uint_t                                      level,
                             DoFType                                     flag ) const;

 private:
   const std::function< void( const Point3D&, Point3D& ) > normal_function_;
};

} // namespace hyteg
