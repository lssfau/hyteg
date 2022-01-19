/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/dgfunctionspace/DGOperator.hpp"

#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"

namespace hyteg {
namespace dg {

DGOperator::DGOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                        uint_t                                     minLevel,
                        uint_t                                     maxLevel,
                        const std::shared_ptr< DGForm >&           form )
: Operator< DGFunction< real_t >, DGFunction< real_t > >( storage, minLevel, maxLevel )
, form_( form )
{}

void DGOperator::apply( const DGFunction< real_t >& src,
                        const DGFunction< real_t >& dst,
                        size_t                      level,
                        DoFType                     flag,
                        UpdateType                  updateType ) const
{
   // TODO: communicate

   applyInner( src, dst, level, flag, updateType );

   // TODO: apply inner facets
   // TODO: apply boundary facets
};

void DGOperator::toMatrix( const std::shared_ptr< SparseMatrixProxy >&             mat,
                           const typename srcType::template FunctionType< idx_t >& src,
                           const typename dstType::template FunctionType< idx_t >& dst,
                           size_t                                                  level,
                           DoFType                                                 flag ) const
{
   WALBERLA_ABORT( "toMatrix() not implemented in DGOperator!" );
}

void DGOperator::smooth_jac( const DGFunction< real_t >& dst,
                             const DGFunction< real_t >& rhs,
                             const DGFunction< real_t >& tmp,
                             real_t                      relax,
                             size_t                      level,
                             DoFType                     flag ) const
{
   WALBERLA_ABORT( "DGOperator: weighted Jacobi not implemented." );
}

std::shared_ptr< DGFunction< real_t > > DGOperator::getInverseDiagonalValues() const
{
   WALBERLA_ABORT( "DGOperator: inverse diagonal not implemented." );
}

void DGOperator::applyInner( const DGFunction< real_t >& src,
                             const DGFunction< real_t >& dst,
                             size_t                      level,
                             DoFType                     flag,
                             UpdateType                  updateType ) const
{
   if ( this->getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "DG apply not implemented in 3D." );
   }
   else
   {
      for ( auto faceIt : this->getStorage()->getFaces() )
      {
         const auto faceId = faceIt.first;
         const auto face   = *faceIt.second;

         const auto srcPolyDegree = src.polynomialDegree( faceId );
         const auto dstPolyDegree = dst.polynomialDegree( faceId );

         const auto numSrcDofs = src.basis()->numDoFsPerElement( srcPolyDegree );
         const auto numDstDofs = dst.basis()->numDoFsPerElement( dstPolyDegree );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( auto elementIdx : facedof::macroface::Iterator( level, faceType ) )
            {
               // First we assemble the local element matrix.
               // For that we need the (affine) coordinates of the local element.

               const auto microVertexIndices = facedof::macroface::getMicroVerticesFromMicroFace( elementIdx, faceType );

               std::array< Eigen::Matrix< real_t, 2, 1 >, 3 > vertexCoords;
               for ( uint_t i = 0; i < 3; i++ )
               {
                  const auto coord     = vertexdof::macroface::coordinateFromIndex( level, face, microVertexIndices[i] );
                  vertexCoords[i]( 0 ) = coord[0];
                  vertexCoords[i]( 1 ) = coord[1];
               }

               Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > mat;
               form_->integrateVolume( vertexCoords, *src.basis(), *dst.basis(), srcPolyDegree, dstPolyDegree, mat );

               // Now we need to perform the multiplication with the source DoFs.

               Eigen::Matrix< real_t, Eigen::Dynamic, 1 > srcDofs;
               Eigen::Matrix< real_t, Eigen::Dynamic, 1 > dstDofs;
               srcDofs.resize( numSrcDofs, Eigen::NoChange_t::NoChange );
               dstDofs.resize( numDstDofs, Eigen::NoChange_t::NoChange );

               for ( int srcDofIdx = 0; srcDofIdx < numSrcDofs; srcDofIdx++ )
               {
                  // TODO: this is not optimal -> better to gather the pointer outside of this loop and then only access
                  srcDofs( srcDofIdx ) = src.volumeDoFFunction()->dof( faceId, elementIdx, srcDofIdx, faceType, level );
               }

               dstDofs = mat * srcDofs;

               for ( int dstDofIdx = 0; dstDofIdx < numDstDofs; dstDofIdx++ )
               {
                  // TODO: this is not optimal -> better to gather the pointer outside of this loop and then only access
                  dst.volumeDoFFunction()->dof( faceId, elementIdx, dstDofIdx, faceType, level ) = dstDofs( dstDofIdx );
               }
            }
         }
      }
   }
}

} // namespace dg
} // namespace hyteg
