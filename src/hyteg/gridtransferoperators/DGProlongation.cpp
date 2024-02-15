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

#include "DGProlongation.hpp"

#include "hyteg/indexing/MacroEdgeIndexing.hpp"

namespace hyteg {

inline real_t phi0( const Point2D& x )
{
   return 1 - x( 0 ) - x( 1 );
}

inline real_t phi1( const Point2D& x )
{
   return x( 0 );
}

inline real_t phi2( const Point2D& x )
{
   return x( 1 );
}

inline real_t phi0( const Point3D& x )
{
   return 1 - x( 0 ) - x( 1 ) - x( 2 );
}

inline real_t phi1( const Point3D& x )
{
   return x( 0 );
}

inline real_t phi2( const Point3D& x )
{
   return x( 1 );
}

inline real_t phi3( const Point3D& x )
{
   return x( 2 );
}

void ProlongationFormDG1::integrate2D( const std::vector< Point3D >& dst,
                                       const std::vector< Point3D >& src,
                                       MatrixXr&                     localMat ) const
{
   Point2D b;
   b << src[0]( 0 ), src[0]( 1 );
   const Point3D s10 = src[1] - src[0];
   const Point3D s20 = src[2] - src[0];

   Point2D a1;
   Point2D a2;
   a1 << s10( 0 ), s10( 1 );
   a2 << s20( 0 ), s20( 1 );

   Matrix2r A( 2, 2 );
   A << a1( 0 ), a2( 0 ), a1( 1 ), a2( 1 );

   Matrix2r Ainv{ A.inverse() };

   Point2D d0;
   Point2D d1;
   Point2D d2;
   d0 << dst[0]( 0 ), dst[0]( 1 );
   d1 << dst[1]( 0 ), dst[1]( 1 );
   d2 << dst[2]( 0 ), dst[2]( 1 );

   d0 -= b;
   d1 -= b;
   d2 -= b;

   const Point2D p0{ Ainv * d0 };
   const Point2D p1{ Ainv * d1 };
   const Point2D p2{ Ainv * d2 };

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

void ProlongationFormDG1::integrate3D( const std::vector< Point3D >& dst,
                                       const std::vector< Point3D >& src,
                                       MatrixXr&                     localMat ) const
{
   const Point3D b{ src[0] };

   const Point3D s10{ src[1] - src[0] };
   const Point3D s20{ src[2] - src[0] };
   const Point3D s30{ src[3] - src[0] };

   Matrix3r A( 3, 3 );
   A << s10, s20, s30;

   Matrix3r Ainv{ A.inverse() };

   Point3D d0 = dst[0];
   Point3D d1 = dst[1];
   Point3D d2 = dst[2];
   Point3D d3 = dst[3];

   d0 -= b;
   d1 -= b;
   d2 -= b;
   d3 -= b;

   const Point3D p0{ Ainv * d0 };
   const Point3D p1{ Ainv * d1 };
   const Point3D p2{ Ainv * d2 };
   const Point3D p3{ Ainv * d3 };

   localMat.resize( 4, 4 );
   localMat( 0, 0 ) = phi0( p0 );
   localMat( 0, 1 ) = phi1( p0 );
   localMat( 0, 2 ) = phi2( p0 );
   localMat( 0, 3 ) = phi3( p0 );
   localMat( 1, 0 ) = phi0( p1 );
   localMat( 1, 1 ) = phi1( p1 );
   localMat( 1, 2 ) = phi2( p1 );
   localMat( 1, 3 ) = phi3( p1 );
   localMat( 2, 0 ) = phi0( p2 );
   localMat( 2, 1 ) = phi1( p2 );
   localMat( 2, 2 ) = phi2( p2 );
   localMat( 2, 3 ) = phi3( p2 );
   localMat( 3, 0 ) = phi0( p3 );
   localMat( 3, 1 ) = phi1( p3 );
   localMat( 3, 2 ) = phi2( p3 );
   localMat( 3, 3 ) = phi3( p3 );
}

void DGProlongation::prolongate( const dg::DGFunction< real_t >& function,
                                 const walberla::uint_t&         coarseLevel,
                                 const DoFType&                  flag ) const
{
   WALBERLA_UNUSED( flag );

   using volumedofspace::indexing::ElementNeighborInfo;

   const uint_t fineLevel = coarseLevel + 1;

   const auto   storage = function.getStorage();
   const uint_t dim     = storage->hasGlobalCells() ? 3 : 2;

   std::vector< PrimitiveID > pids;
   if ( dim == 2 )
   {
      pids = storage->getFaceIDs();
   }
   else
   {
      pids = storage->getCellIDs();
   }

   const uint_t numMicroVolTypes = ( storage->hasGlobalCells() ? 6 : 2 );

   for ( const auto& pid : pids )
   {
      const auto polyDegree = function.polynomialDegree( pid );
      const int  numDofs = static_cast< int >( function.basis()->numDoFsPerElement( dim, static_cast< uint_t >( polyDegree ) ) );

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
               coarseNeighborInfo = ElementNeighborInfo(
                   coarseElementIdx, coarseCellType, coarseLevel, function.getBoundaryCondition(), pid, storage );
            }

            VectorXr coarseDofs;
            coarseDofs.resize( numDofs, Eigen::NoChange_t::NoChange );

            for ( int coarseDofIdx = 0; coarseDofIdx < numDofs; coarseDofIdx++ )
            {
               if ( dim == 2 )
               {
                  coarseDofs( coarseDofIdx ) = coarseDofMemory[volumedofspace::indexing::index( coarseElementIdx.x(),
                                                                                                coarseElementIdx.y(),
                                                                                                coarseFaceType,
                                                                                                uint_c( coarseDofIdx ),
                                                                                                uint_c( numDofs ),
                                                                                                coarseLevel,
                                                                                                memLayout )];
               }
               else
               {
                  coarseDofs( coarseDofIdx ) = coarseDofMemory[volumedofspace::indexing::index( coarseElementIdx.x(),
                                                                                                coarseElementIdx.y(),
                                                                                                coarseElementIdx.z(),
                                                                                                coarseCellType,
                                                                                                uint_c( coarseDofIdx ),
                                                                                                uint_c( numDofs ),
                                                                                                coarseLevel,
                                                                                                memLayout )];
               }
            }

            std::vector< hyteg::indexing::Index > fineElementIndices;
            std::vector< facedof::FaceType >      fineFaceTypes;
            std::vector< celldof::CellType >      fineCellTypes;

            if ( dim == 2 )
            {
               volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                   coarseElementIdx, coarseFaceType, fineElementIndices, fineFaceTypes );
            }
            else
            {
               volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                   coarseElementIdx, coarseCellType, fineElementIndices, fineCellTypes );
            }

            for ( uint_t fineIdx = 0; fineIdx < fineElementIndices.size(); fineIdx += 1 )
            {
               auto fineElementIdx = fineElementIndices[fineIdx];

               facedof::FaceType fineFaceType{ facedof::FaceType::GRAY };
               celldof::CellType fineCellType{ celldof::CellType::BLUE_DOWN };
               if ( dim == 2 )
               {
                  fineFaceType = fineFaceTypes[fineIdx];
               }
               else
               {
                  fineCellType = fineCellTypes[fineIdx];
               }

               volumedofspace::indexing::ElementNeighborInfo fineNeighborInfo;
               if ( dim == 2 )
               {
                  fineNeighborInfo = ElementNeighborInfo(
                      fineElementIdx, fineFaceType, fineLevel, function.getBoundaryCondition(), pid, storage );
               }
               else
               {
                  fineNeighborInfo = ElementNeighborInfo(
                      fineElementIdx, fineCellType, fineLevel, function.getBoundaryCondition(), pid, storage );
               }

               MatrixXr localMat;
               localMat.resize( numDofs, numDofs );

               // assemble local matrix
               if ( dim == 2 )
               {
                  form_->integrate2D(
                      fineNeighborInfo.elementVertexCoords(), coarseNeighborInfo.elementVertexCoords(), localMat );
               }
               else
               {
                  form_->integrate3D(
                      fineNeighborInfo.elementVertexCoords(), coarseNeighborInfo.elementVertexCoords(), localMat );
               }

               VectorXr fineDofs;
               fineDofs.resize( numDofs, Eigen::NoChange_t::NoChange );

               fineDofs = localMat * coarseDofs;

               for ( int fineDofIdx = 0; fineDofIdx < numDofs; fineDofIdx++ )
               {
                  if ( dim == 2 )
                  {
                     fineDofMemory[volumedofspace::indexing::index( fineElementIdx.x(),
                                                                    fineElementIdx.y(),
                                                                    fineFaceType,
                                                                    uint_c( fineDofIdx ),
                                                                    uint_c( numDofs ),
                                                                    fineLevel,
                                                                    memLayout )] = fineDofs( fineDofIdx );
                  }
                  else
                  {
                     fineDofMemory[volumedofspace::indexing::index( fineElementIdx.x(),
                                                                    fineElementIdx.y(),
                                                                    fineElementIdx.z(),
                                                                    fineCellType,
                                                                    uint_c( fineDofIdx ),
                                                                    uint_c( numDofs ),
                                                                    fineLevel,
                                                                    memLayout )] = fineDofs( fineDofIdx );
                  }
               }
            }
         }
      }
   }
}

} // namespace hyteg