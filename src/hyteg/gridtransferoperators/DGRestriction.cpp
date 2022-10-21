/*
 * Copyright (c) 2022 Andreas Wagner.
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

#include "DGRestriction.hpp"

#include "hyteg/indexing/MacroEdgeIndexing.hpp"

namespace hyteg {

void RestrictionFormDG1::integrate2D( const std::vector< Point >&                             dst,
                                      const std::vector< Point >&                             src,
                                      Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const
{
   using Point2 = Eigen::Matrix< real_t, 2, 1 >;

   const Point2 b  = src[0].block( 0, 0, 2, 1 );
   const Point2 a1 = ( src[1] - src[0] ).block( 0, 0, 2, 1 );
   const Point2 a2 = ( src[2] - src[0] ).block( 0, 0, 2, 1 );

   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > A( 2, 2 );
   A.col( 0 ) = a1;
   A.col( 1 ) = a2;

   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > Ainv = A.inverse();

   Point2 p0 = Ainv * ( dst[0].block( 0, 0, 2, 1 ) - b );
   Point2 p1 = Ainv * ( dst[1].block( 0, 0, 2, 1 ) - b );
   Point2 p2 = Ainv * ( dst[2].block( 0, 0, 2, 1 ) - b );

   auto phi0 = []( auto x ) { return 1 - x[0] - x[1]; };
   auto phi1 = []( auto x ) { return x[0]; };
   auto phi2 = []( auto x ) { return x[1]; };

   localMat.resize( 3, 3 );
   localMat( 0, 0 ) = phi0( p0 );
   localMat( 0, 1 ) = phi0( p1 );
   localMat( 0, 2 ) = phi0( p2 );
   localMat( 1, 0 ) = phi1( p0 );
   localMat( 1, 1 ) = phi1( p1 );
   localMat( 1, 2 ) = phi1( p2 );
   localMat( 2, 0 ) = phi2( p0 );
   localMat( 2, 1 ) = phi2( p1 );
   localMat( 2, 2 ) = phi2( p2 );
}

void DGRestriction::restrict( const dg::DGFunction< real_t >& function,
                              const walberla::uint_t&         fineLevel,
                              const DoFType&                  flag ) const
{
   WALBERLA_UNUSED( flag );

   using volumedofspace::indexing::ElementNeighborInfo;

   const uint_t coarseLevel = fineLevel - 1;

   const uint_t dim     = 2;
   const auto   storage = function.getStorage();

   if ( storage->hasGlobalCells() )
      WALBERLA_ABORT( "Prolongation currently only supports 2D." )

   std::vector< PrimitiveID > pids = storage->getFaceIDs();

   const uint_t numMicroVolTypes = ( storage->hasGlobalCells() ? 6 : 2 );

   for ( const auto& pid : pids )
   {
      const auto polyDegree = function.polynomialDegree( pid );
      const auto numDofs    = function.basis()->numDoFsPerElement( dim, polyDegree );

      const auto coarseDofMemory = function.volumeDoFFunction()->dofMemory( pid, coarseLevel );
      const auto fineDofMemory   = function.volumeDoFFunction()->dofMemory( pid, fineLevel );

      const auto memLayout = function.volumeDoFFunction()->memoryLayout();

      for ( uint_t coarseMicroVolType = 0; coarseMicroVolType < numMicroVolTypes; coarseMicroVolType++ )
      {
         auto coarseFaceType = facedof::allFaceTypes[coarseMicroVolType];
         auto coarseCellType = celldof::allCellTypes[coarseMicroVolType];

         auto itFace = facedof::macroface::Iterator( coarseLevel, coarseFaceType ).begin();
         auto itCell = celldof::macrocell::Iterator( coarseLevel, coarseCellType ).begin();

         while ( ( dim == 2 && itFace != itFace.end() ) || ( dim == 3 && itCell != itCell.end() ) )
         {
            indexing::Index coarseElementIdx;

            if ( dim == 2 )
            {
               coarseElementIdx = *itFace;
               itFace++;
            }
            else
            {
               coarseElementIdx = *itCell;
               itCell++;
            }

            ElementNeighborInfo coarseNeighborInfo;
            if ( dim == 2 )
            {
               coarseNeighborInfo = ElementNeighborInfo(
                   coarseElementIdx, coarseFaceType, coarseLevel, function.getBoundaryCondition(), pid, storage );
            }
            else
            {
               WALBERLA_LOG_INFO_ON_ROOT( "not implemented" );
            }

            Eigen::Matrix< real_t, Eigen::Dynamic, 1 > coarseDofs;
            coarseDofs.resize( numDofs, Eigen::NoChange_t::NoChange );
            coarseDofs.setZero();

            std::vector< hyteg::indexing::Index > fineElementIndices;
            std::vector< facedof::FaceType >      fineFaceTypes;

            volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                coarseElementIdx, coarseFaceType, fineElementIndices, fineFaceTypes );

            for ( uint_t fineIdx = 0; fineIdx < fineElementIndices.size(); fineIdx += 1 )
            {
               auto fineElementIdx = fineElementIndices[fineIdx];
               auto fineFaceType   = fineFaceTypes[fineIdx];

               volumedofspace::indexing::ElementNeighborInfo fineNeighborInfo;
               if ( dim == 2 )
               {
                  fineNeighborInfo = ElementNeighborInfo(
                      fineElementIdx, fineFaceType, fineLevel, function.getBoundaryCondition(), pid, storage );
               }
               else
               {
                  WALBERLA_LOG_INFO_ON_ROOT( "not implemented" );
               }

               Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic > localMat;
               localMat.resize( numDofs, numDofs );

               // assemble local matrix
               if ( dim == 2 )
               {
                  form_->integrate2D(
                      fineNeighborInfo.elementVertexCoords(), coarseNeighborInfo.elementVertexCoords(), localMat );
               }
               else
               {
                  WALBERLA_LOG_INFO_ON_ROOT( "not implemented" );
               }

               Eigen::Matrix< real_t, Eigen::Dynamic, 1 > fineDofs;
               fineDofs.resize( numDofs, Eigen::NoChange_t::NoChange );
               for ( uint_t fineDofIdx = 0; fineDofIdx < numDofs; fineDofIdx++ )
               {
                  if ( dim == 2 )
                  {
                     fineDofs( fineDofIdx ) = fineDofMemory[volumedofspace::indexing::index(
                         fineElementIdx.x(), fineElementIdx.y(), fineFaceType, fineDofIdx, numDofs, fineLevel, memLayout )];
                  }
                  else
                  {
                     WALBERLA_LOG_INFO_ON_ROOT( "not implemented" );
                  }
               }

               coarseDofs += localMat * fineDofs;
            }

            for ( uint_t coarseDofIdx = 0; coarseDofIdx < numDofs; coarseDofIdx++ )
            {
               if ( dim == 2 )
               {
                  coarseDofMemory[volumedofspace::indexing::index( coarseElementIdx.x(),
                                                                   coarseElementIdx.y(),
                                                                   coarseFaceType,
                                                                   coarseDofIdx,
                                                                   numDofs,
                                                                   coarseLevel,
                                                                   memLayout )] = coarseDofs( coarseDofIdx );
               }
               else
               {
                  WALBERLA_LOG_INFO_ON_ROOT( "not implemented" );
               }
            }
         }
      }
   }
}

} // namespace hyteg