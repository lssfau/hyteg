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
#pragma once

#include "core/Abort.h"
#include "core/DataTypes.h"

#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

class ProlongationForm
{
 public:
   using Point = Eigen::Matrix< real_t, 3, 1 >;

   virtual void integrate2D( const std::vector< Point >&                              dst,
                             const std::vector< Point >&                              src,
                             Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const
   {
      WALBERLA_LOG_INFO_ON_ROOT( "not implemented" )
   }
};

class ProlongationFormDG1 : public ProlongationForm
{
 public:
   using Point3 = Eigen::Matrix< real_t, 3, 1 >;
   using Point2 = Eigen::Matrix< real_t, 2, 1 >;

   void integrate2D( const std::vector< Point3 >&                             dst,
                     const std::vector< Point3 >&                             src,
                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const override
   {
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
      localMat( 0, 1 ) = phi1( p0 );
      localMat( 0, 2 ) = phi2( p0 );
      localMat( 1, 0 ) = phi0( p1 );
      localMat( 1, 1 ) = phi1( p1 );
      localMat( 1, 2 ) = phi2( p1 );
      localMat( 2, 0 ) = phi0( p2 );
      localMat( 2, 1 ) = phi1( p2 );
      localMat( 2, 2 ) = phi2( p2 );
   }
};

class ProlongationFormDG0 : public ProlongationForm
{
 public:
   using Point3 = Eigen::Matrix< real_t, 3, 1 >;

   void integrate2D( const std::vector< Point3 >&                             dst,
                     const std::vector< Point3 >&                             src,
                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& localMat ) const override
   {
      localMat.resize( 1, 1 );
      localMat( 0, 0 ) = 1;
   }
};

class DGProlongation : public ProlongationOperator< dg::DGFunction< real_t > >
{
 public:
   explicit DGProlongation( std::shared_ptr< ProlongationForm > form )
   : form_( form )
   {}

   void prolongate( const dg::DGFunction< real_t >& function,
                    const walberla::uint_t&         srcLevel,
                    const DoFType&                  flag ) const override
   {
      using volumedofspace::indexing::ElementNeighborInfo;

      const uint_t dstLevel = srcLevel + 1;

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

         const auto srcDofMemory = function.volumeDoFFunction()->dofMemory( pid, srcLevel );
         const auto dstDofMemory = function.volumeDoFFunction()->dofMemory( pid, dstLevel );

         const auto memLayout = function.volumeDoFFunction()->memoryLayout();

         for ( uint_t srcMicroVolType = 0; srcMicroVolType < numMicroVolTypes; srcMicroVolType++ )
         {
            auto srcFaceType = facedof::allFaceTypes[srcMicroVolType];
            auto srcCellType = celldof::allCellTypes[srcMicroVolType];

            auto itFace = facedof::macroface::Iterator( srcLevel, srcFaceType ).begin();
            auto itCell = celldof::macrocell::Iterator( srcLevel, srcCellType ).begin();

            while ( ( dim == 2 && itFace != itFace.end() ) || ( dim == 3 && itCell != itCell.end() ) )
            {
               indexing::Index srcElementIdx;

               if ( dim == 2 )
               {
                  srcElementIdx = *itFace;
                  itFace++;
               }
               else
               {
                  srcElementIdx = *itCell;
                  itCell++;
               }

               volumedofspace::indexing::ElementNeighborInfo srcNeighborInfo;
               if ( dim == 2 )
               {
                  srcNeighborInfo =
                      ElementNeighborInfo( srcElementIdx, srcFaceType, srcLevel, function.getBoundaryCondition(), pid, storage );
               }
               else
               {
                  WALBERLA_LOG_INFO_ON_ROOT( "not implemented" );
               }

               Eigen::Matrix< real_t, Eigen::Dynamic, 1 > srcDofs;
               srcDofs.resize( numDofs, Eigen::NoChange_t::NoChange );

               for ( uint_t srcDofIdx = 0; srcDofIdx < numDofs; srcDofIdx++ )
               {
                  if ( dim == 2 )
                  {
                     srcDofs( srcDofIdx ) = srcDofMemory[volumedofspace::indexing::index(
                         srcElementIdx.x(), srcElementIdx.y(), srcFaceType, srcDofIdx, numDofs, srcLevel, memLayout )];
                  }
                  else
                  {
                     WALBERLA_LOG_INFO_ON_ROOT( "not implemented" );
                  }
               }

               std::vector< hyteg::indexing::Index > dstElementIndices;
               std::vector< facedof::FaceType >      dstFaceTypes;

               volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                   srcElementIdx, srcFaceType, dstElementIndices, dstFaceTypes );

               for ( uint_t dstIdx = 0; dstIdx < dstElementIndices.size(); dstIdx += 1 )
               {
                  auto dstElementIdx = dstElementIndices[dstIdx];
                  auto dstFaceType   = dstFaceTypes[dstIdx];

                  volumedofspace::indexing::ElementNeighborInfo dstNeighborInfo;
                  if ( dim == 2 )
                  {
                     dstNeighborInfo = ElementNeighborInfo(
                         dstElementIdx, dstFaceType, dstLevel, function.getBoundaryCondition(), pid, storage );
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
                     form_->integrate2D( dstNeighborInfo.elementVertexCoords(), srcNeighborInfo.elementVertexCoords(), localMat );
                  }
                  else
                  {
                     WALBERLA_LOG_INFO_ON_ROOT( "not implemented" );
                  }

                  Eigen::Matrix< real_t, Eigen::Dynamic, 1 > dstDofs;
                  dstDofs.resize( numDofs, Eigen::NoChange_t::NoChange );
                  dstDofs.setZero();

                  dstDofs += localMat * srcDofs;

                  for ( uint_t dstDofIdx = 0; dstDofIdx < numDofs; dstDofIdx++ )
                  {
                     if ( dim == 2 )
                     {
                        dstDofMemory[volumedofspace::indexing::index(
                            dstElementIdx.x(), dstElementIdx.y(), dstFaceType, dstDofIdx, numDofs, dstLevel, memLayout )] =
                            dstDofs( dstDofIdx );
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
   };

 protected:
   std::shared_ptr< ProlongationForm > form_;
};

class DG1Prolongation : public DGProlongation
{
 public:
   DG1Prolongation()
   : DGProlongation( std::make_shared< ProlongationFormDG1 >() )
   {}
};

class DG0Prolongation : public DGProlongation
{
 public:
   DG0Prolongation()
       : DGProlongation( std::make_shared< ProlongationFormDG0 >() )
   {}
};

} // namespace hyteg
