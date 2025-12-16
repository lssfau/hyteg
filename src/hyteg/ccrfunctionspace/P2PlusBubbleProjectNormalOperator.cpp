/*
 * Copyright (c) 2025 Marcus Mohr.
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

#include "P2PlusBubbleProjectNormalOperator.hpp"

namespace hyteg {

P2PlusBubbleProjectNormalOperator::P2PlusBubbleProjectNormalOperator(
    const std::shared_ptr< PrimitiveStorage >&               storage,
    size_t                                                   minLevel,
    size_t                                                   maxLevel,
    const std::function< void( const Point3D&, Point3D& ) >& normalFunction )
: Operator( storage, minLevel, maxLevel )
, p1Operator_( storage, minLevel, maxLevel, normalFunction )
, edgeDoFOperator_( storage, minLevel, maxLevel, normalFunction )
{}

void P2PlusBubbleProjectNormalOperator::project( const P2PlusBubbleVectorFunction< real_t >& dst,
                                                 size_t                                      level,
                                                 DoFType                                     flag ) const
{
   // we only need to project the P2 component, but not the bubble part

   if ( dst.getDimension() == 3 )
   {
      WALBERLA_ABORT( "P2PlusBubbleVectorFunction of dimension " << dst.getDimension() << ". Does this already work?" );
   }
   else if ( dst.getDimension() == 2 )
   {
      p1Operator_.project( dst[0].getVertexDoFFunction(), dst[1].getVertexDoFFunction(), level, flag );
      edgeDoFOperator_.project( dst[0].getEdgeDoFFunction(), dst[1].getEdgeDoFFunction(), level, flag );
   }
   else
   {
      WALBERLA_ABORT( "Cannot deal with P2PlusBubbleVectorFunction of dimension " << dst.getDimension() );
   }
}

void P2PlusBubbleProjectNormalOperator::project( const CCRStokesFunction< real_t >& dst, size_t level, DoFType flag ) const
{
   project( dst.uvw(), level, flag );
}

void P2PlusBubbleProjectNormalOperator::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                  const P2PlusBubbleVectorFunction< idx_t >&  numSrc,
                                                  const P2PlusBubbleVectorFunction< idx_t >&  numDst,
                                                  uint_t                                      level,
                                                  DoFType                                     flag ) const
{
   WALBERLA_UNUSED( numDst );

   if ( numSrc.getDimension() == 3 )
   {
      WALBERLA_ABORT( "P2PlusBubbleVectorFunction of dimension " << numSrc.getDimension() << ". Does this already work?" );
   }
   else if ( numSrc.getDimension() == 2 )
   {
      p1Operator_.toMatrix( mat, numSrc[0].getVertexDoFFunction(), numSrc[1].getVertexDoFFunction(), level, flag );
      edgeDoFOperator_.assembleLocalMatrix( mat, numSrc[0].getEdgeDoFFunction(), numSrc[1].getEdgeDoFFunction(), level, flag );

      // The bubble component only lives on faces. IT not affected by the projection, thus, we would need to add an
      // identity matrix for this part. However, we implement this projection operator for use in treating free-slip
      // boundary conditions. In 2D this will not involve faces. So we do currently _not_ implement identity matrix
      // assembly.
      for ( const auto& it : numSrc.getStorage()->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = numSrc.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            WALBERLA_ABORT( "Did not expect to have to deal with faces in assembly. Implement assembly of an identity matrix"
                            << " for the bubble component, if this really is what you need." )
         }
      }
   }
   else
   {
      WALBERLA_ABORT( "Cannot deal with P2PlusBubbleVectorFunction of dimension " << numSrc.getDimension() );
   }
}

} // namespace hyteg
