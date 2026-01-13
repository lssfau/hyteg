/*
 * Copyright (c) 2020-2025 Daniel Drzisga, Marcus Mohr.
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

void P2ProjectNormalOperator::project( const P2Function< real_t >& dst_u,
                                       const P2Function< real_t >& dst_v,
                                       const P2Function< real_t >& dst_w,
                                       size_t                      level,
                                       DoFType                     flag ) const
{
   p1Operator.project( dst_u.getVertexDoFFunction(), dst_v.getVertexDoFFunction(), dst_w.getVertexDoFFunction(), level, flag );
   edgeDoFOperator.project( dst_u.getEdgeDoFFunction(), dst_v.getEdgeDoFFunction(), dst_w.getEdgeDoFFunction(), level, flag );
}

void P2ProjectNormalOperator::project( const P2Function< real_t >& dst_u,
                                       const P2Function< real_t >& dst_v,
                                       size_t                      level,
                                       DoFType                     flag ) const
{
   p1Operator.project( dst_u.getVertexDoFFunction(), dst_v.getVertexDoFFunction(), level, flag );
   edgeDoFOperator.project( dst_u.getEdgeDoFFunction(), dst_v.getEdgeDoFFunction(), level, flag );
}

void P2ProjectNormalOperator::project( const P2VectorFunction< real_t >& dst, size_t level, DoFType flag ) const
{
   if ( dst.getDimension() == 3 )
   {
      project( dst[0], dst[1], dst[2], level, flag );
   }
   else if ( dst.getDimension() == 2 )
   {
      project( dst[0], dst[1], level, flag );
   }
   else
   {
      WALBERLA_ABORT( "Cannot deal with P2VectorFunction of dimension " << dst.getDimension() );
   }
}

void P2ProjectNormalOperator::project( const P2P1TaylorHoodFunction< real_t >& dst, size_t level, DoFType flag ) const
{
   if ( dst.uvw().getDimension() == 3 )
   {
      project( dst.uvw()[0], dst.uvw()[1], dst.uvw()[2], level, flag );
   }
   else if ( dst.uvw().getDimension() == 2 )
   {
      project( dst.uvw()[0], dst.uvw()[1], level, flag );
   }
   else
   {
      WALBERLA_ABORT( "Cannot deal with P2P1TaylorHoodFunction of dimension " << dst.uvw().getDimension() );
   }
}

void P2ProjectNormalOperator::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                        const P2Function< idx_t >&                  numU,
                                        const P2Function< idx_t >&                  numV,
                                        const P2Function< idx_t >&                  numW,
                                        uint_t                                      level,
                                        DoFType                                     flag ) const
{
   p1Operator.toMatrix( mat, numU.getVertexDoFFunction(), numV.getVertexDoFFunction(), numW.getVertexDoFFunction(), level, flag );
   edgeDoFOperator.assembleLocalMatrix(
       mat, numU.getEdgeDoFFunction(), numV.getEdgeDoFFunction(), numW.getEdgeDoFFunction(), level, flag );
}

void P2ProjectNormalOperator::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                        const P2VectorFunction< idx_t >&            num,
                                        uint_t                                      level,
                                        DoFType                                     flag,
                                        bool                                        transpose ) const
{
   if ( num.getDimension() == 3 )
   {
      toMatrix( mat, num[0], num[1], num[2], level, flag );
   }
   else if ( num.getDimension() == 2 )
   {
      p1Operator.toMatrix( mat, num[0].getVertexDoFFunction(), num[1].getVertexDoFFunction(), level, flag );
      edgeDoFOperator.assembleLocalMatrix( mat, num[0].getEdgeDoFFunction(), num[1].getEdgeDoFFunction(), level, flag );
   }
   else
   {
      WALBERLA_ABORT( "Encountered unsupported field dimension " << num.getDimension() );
   }
}

} // namespace hyteg
