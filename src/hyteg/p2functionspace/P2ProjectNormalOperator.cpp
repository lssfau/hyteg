/*
 * Copyright (c) 2020 Daniel Drzisga.
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

#include "P2ProjectNormalOperator.hpp"

namespace hyteg {

P2ProjectNormalOperator::P2ProjectNormalOperator( const std::shared_ptr< PrimitiveStorage >&               storage,
                                                  size_t                                                   minLevel,
                                                  size_t                                                   maxLevel,
                                                  const std::function< void( const Point3D&, Point3D& ) >& normal_function )
: Operator( storage, minLevel, maxLevel )
, p1Operator( storage, minLevel, maxLevel, normal_function )
, edgeDoFOperator( storage, minLevel, maxLevel, normal_function )
{}

void P2ProjectNormalOperator::apply( const P2Function< real_t >& dst_u,
                                     const P2Function< real_t >& dst_v,
                                     const P2Function< real_t >& dst_w,
                                     size_t                      level,
                                     DoFType                     flag ) const
{
   p1Operator.apply( dst_u.getVertexDoFFunction(), dst_v.getVertexDoFFunction(), dst_w.getVertexDoFFunction(), level, flag );
   edgeDoFOperator.apply( dst_u.getEdgeDoFFunction(), dst_v.getEdgeDoFFunction(), dst_w.getEdgeDoFFunction(), level, flag );
}

void P2ProjectNormalOperator::apply( const P2P1TaylorHoodFunction< real_t >& dst, size_t level, DoFType flag ) const
{
   apply( dst.uvw.u, dst.uvw.v, dst.uvw.w, level, flag );
}

#ifdef HYTEG_BUILD_WITH_PETSC

void P2ProjectNormalOperator::assembleLocalMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                   const P2Function< PetscInt >&               numU,
                                                   const P2Function< PetscInt >&               numV,
                                                   const P2Function< PetscInt >&               numW,
                                                   uint_t                                      level,
                                                   DoFType                                     flag )
{
   p1Operator.assembleLocalMatrix(
       mat, numU.getVertexDoFFunction(), numV.getVertexDoFFunction(), numW.getVertexDoFFunction(), level, flag );
   edgeDoFOperator.assembleLocalMatrix(
       mat, numU.getEdgeDoFFunction(), numV.getEdgeDoFFunction(), numW.getEdgeDoFFunction(), level, flag );
}
#endif

} // namespace hyteg
