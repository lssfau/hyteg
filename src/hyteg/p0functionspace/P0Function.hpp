/*
* Copyright (c) 2017-2024 Ponsuganth Ilangovan, Nils Kohl, Marcus Mohr.
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

#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

namespace hyteg {

using namespace dg;

template < typename ValueType >
class P0Function : public Function< P0Function< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = P0Function< VType >;

   P0Function( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               BoundaryCondition                          bc )
   : Function< P0Function< ValueType > >( name, storage, minLevel, maxLevel )
   {
      auto basis  = std::make_shared< DGBasisLinearLagrange_Example >();
      dgFunction_ = std::make_shared< DGFunction< ValueType > >( name, storage, minLevel, maxLevel, basis, 0, bc );
   }

   P0Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : P0Function( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   const std::shared_ptr< DGFunction< ValueType > > getDGFunction() const { return dgFunction_; }

   void setBoundaryCondition( BoundaryCondition bc ) { dgFunction_->setBoundaryCondition( bc ); }

   BoundaryCondition getBoundaryCondition() const { return dgFunction_->getBoundaryCondition(); }

   uint_t getDimension() const override final { return dgFunction_->getDimension(); };

   void setDoNotWarnOnInterpolateFlag() { doNotWarnOnInterpolate_ = true; }

   // template < typename SenderType, typename ReceiverType >
   void communicate( const uint_t& level ) const { dgFunction_->communicate( level ); }

   void add( const ValueType scalar, uint_t level, DoFType flag = All ) const { dgFunction_->add( scalar, level, flag ); };

   void add( const std::vector< ValueType >                                                scalars,
             const std::vector< std::reference_wrapper< const P0Function< ValueType > > >& functions,
             uint_t                                                                        level,
             DoFType                                                                       flag = All ) const
   {
      std::vector< ValueType > new_scalars( scalars );
      new_scalars.push_back( 1.0 );
      std::vector< std::reference_wrapper< const P0Function< ValueType > > > new_functions( functions );
      new_functions.push_back( *this );
      assign( new_scalars, new_functions, level, flag );
   };

   void multElementwise( const std::vector< std::reference_wrapper< const P0Function< ValueType > > >& functions,
                         uint_t                                                                        level,
                         DoFType                                                                       flag = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   void interpolate( ValueType constant, uint_t level, DoFType dofType = All ) const
   {
      if ( !doNotWarnOnInterpolate_ )
         WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::interpolate() 'interpolates' values at the centroid." );
      if ( this->storage_->hasGlobalCells() )
      {
         for ( auto& it : this->getStorage()->getCells() )
         {
            const auto  cellID = it.first;
            const auto& cell   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( cellID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 3, 0 ), 1 );

            const auto memLayout = getDGFunction()->volumeDoFFunction()->memoryLayout();
            auto       dofs      = getDGFunction()->volumeDoFFunction()->dofMemory( cellID, level );

            for ( auto cellType : celldof::allCellTypes )
            {
               for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
               {
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayout )] =
                      constant;
               }
            }
         }
      }
      else
      {
         for ( auto& it : this->getStorage()->getFaces() )
         {
            const auto  faceID = it.first;
            const auto& face   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( faceID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 2, 0 ), 1 );

            const auto memLayout = getDGFunction()->volumeDoFFunction()->memoryLayout();
            auto       dofs      = getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
               {
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayout )] = constant;
               }
            }
         }
      }
   }

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType dofType = All ) const
   {
      WALBERLA_UNUSED( dofType );
      if ( !doNotWarnOnInterpolate_ )
         WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::interpolate() 'interpolates' values at the centroid." );

      if ( this->storage_->hasGlobalCells() )
      {
         for ( auto& it : this->getStorage()->getCells() )
         {
            const auto  cellID = it.first;
            const auto& cell   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( cellID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 3, 0 ), 1 );

            const auto memLayout = getDGFunction()->volumeDoFFunction()->memoryLayout();
            auto       dofs      = getDGFunction()->volumeDoFFunction()->dofMemory( cellID, level );

            for ( auto cellType : celldof::allCellTypes )
            {
               for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
               {
                  const std::array< indexing::Index, 4 > vertexIndices =
                      celldof::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );
                  std::array< Point3D, 4 > elementVertices;
                  for ( uint_t i = 0; i < 4; i++ )
                  {
                     const auto elementVertex = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[i] );
                     elementVertices[i]( 0 )  = elementVertex[0];
                     elementVertices[i]( 1 )  = elementVertex[1];
                     elementVertices[i]( 2 )  = elementVertex[2];
                  }

                  const Point3D centroid =
                      ( elementVertices[0] + elementVertices[1] + elementVertices[2] + elementVertices[3] ) / real_c( 4 );

                  const auto val = expr( Point3D( centroid( 0 ), centroid( 1 ), centroid( 2 ) ) );

                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayout )] =
                      ValueType( val );
               }
            }
         }
      }
      else
      {
         for ( auto& it : this->getStorage()->getFaces() )
         {
            const auto  faceID = it.first;
            const auto& face   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( faceID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 2, 0 ), 1 );

            const auto memLayout = getDGFunction()->volumeDoFFunction()->memoryLayout();
            auto       dofs      = getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
               {
                  const std::array< indexing::Index, 3 > vertexIndices =
                      facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
                  std::array< Point2D, 3 > elementVertices;
                  for ( uint_t i = 0; i < 3; i++ )
                  {
                     const auto elementVertex = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                     elementVertices[i]( 0 )  = elementVertex[0];
                     elementVertices[i]( 1 )  = elementVertex[1];
                  }

                  const Point2D centroid = ( elementVertices[0] + elementVertices[1] + elementVertices[2] ) / real_c( 3 );

                  const auto val = expr( Point3D( centroid( 0 ), centroid( 1 ), 0 ) );

                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayout )] =
                      ValueType( val );
               }
            }
         }
      }
   }

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expressions,
                     uint_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   };

   void interpolate( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expression,
                     const std::vector< std::reference_wrapper< const P0Function< ValueType > > >&        srcFunctions,
                     uint_t                                                                               level,
                     DoFType                                                                              flag = All ) const
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::interpolate() 'interpolates' values at the centroid." );
      if ( this->storage_->hasGlobalCells() )
      {
         for ( auto& it : this->getStorage()->getCells() )
         {
            const auto  cellID = it.first;
            const auto& cell   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( cellID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 3, 0 ), 1 );

            std::vector< ValueType* > dofs;
            dofs.reserve( srcFunctions.size() + 1 );

            std::vector< volumedofspace::indexing::VolumeDoFMemoryLayout > memLayouts;
            memLayouts.reserve( srcFunctions.size() + 1 );

            std::vector< ValueType > srcValues;
            srcValues.reserve( srcFunctions.size() );

            dofs.push_back( getDGFunction()->volumeDoFFunction()->dofMemory( cellID, level ) );
            memLayouts.push_back( getDGFunction()->volumeDoFFunction()->memoryLayout() );

            for ( const auto src : srcFunctions )
            {
               dofs.push_back( src.get().getDGFunction()->volumeDoFFunction()->dofMemory( cellID, level ) );
               memLayouts.push_back( src.get().getDGFunction()->volumeDoFFunction()->memoryLayout() );
            }

            for ( auto cellType : celldof::allCellTypes )
            {
               for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
               {
                  const std::array< indexing::Index, 4 > vertexIndices =
                      celldof::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );
                  std::array< Point3D, 4 > elementVertices;
                  for ( uint_t i = 0; i < 4; i++ )
                  {
                     const auto elementVertex = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[i] );
                     elementVertices[i]( 0 )  = elementVertex[0];
                     elementVertices[i]( 1 )  = elementVertex[1];
                     elementVertices[i]( 2 )  = elementVertex[2];
                  }

                  const Point3D centroid =
                      ( elementVertices[0] + elementVertices[1] + elementVertices[2] + elementVertices[3] ) / real_c( 4 );

                  for ( size_t k = 0; k < srcFunctions.size(); ++k )
                  {
                     srcValues[k] = dofs[k + 1][volumedofspace::indexing::index(
                         idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayouts[k + 1] )];
                  }

                  const auto val = expression( Point3D( centroid( 0 ), centroid( 1 ), centroid( 2 ) ), srcValues );

                  dofs[0]
                      [volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayouts[0] )] =
                          ValueType( val );
               }
            }
         }
      }
      else
      {
         for ( auto& it : this->getStorage()->getFaces() )
         {
            const auto  faceID = it.first;
            const auto& face   = *it.second;

            WALBERLA_CHECK_EQUAL( getDGFunction()->polynomialDegree( faceID ), 0 );
            WALBERLA_CHECK_EQUAL( getDGFunction()->basis()->numDoFsPerElement( 2, 0 ), 1 );

            std::vector< ValueType* > dofs;
            dofs.reserve( srcFunctions.size() + 1 );

            std::vector< volumedofspace::indexing::VolumeDoFMemoryLayout > memLayouts;
            memLayouts.reserve( srcFunctions.size() + 1 );

            std::vector< ValueType > srcValues;
            srcValues.reserve( srcFunctions.size() );

            dofs.push_back( getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level ) );
            memLayouts.push_back( getDGFunction()->volumeDoFFunction()->memoryLayout() );

            for ( const auto src : srcFunctions )
            {
               dofs.push_back( src.get().getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level ) );
               memLayouts.push_back( src.get().getDGFunction()->volumeDoFFunction()->memoryLayout() );
            }

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
               {
                  const std::array< indexing::Index, 3 > vertexIndices =
                      facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
                  std::array< Point2D, 3 > elementVertices;
                  for ( uint_t i = 0; i < 3; i++ )
                  {
                     const auto elementVertex = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                     elementVertices[i]( 0 )  = elementVertex[0];
                     elementVertices[i]( 1 )  = elementVertex[1];
                  }

                  const Point2D centroid = ( elementVertices[0] + elementVertices[1] + elementVertices[2] ) / 3.;

                  for ( size_t k = 0; k < srcFunctions.size(); ++k )
                  {
                     srcValues[k] =
                         dofs[k + 1]
                             [volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayouts[k + 1] )];
                  }

                  const auto newValue = expression( Point3D( centroid( 0 ), centroid( 1 ), 0 ), srcValues );

                  dofs[0][volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayouts[0] )] =
                      ValueType( newValue );
               }
            }
         }
      }
   };

   void transferToLowerLevel( uint_t level )
   {
      if( this->storage_->hasGlobalCells() )
      {
         for ( auto& it : this->storage_->getCells() )
         {
            PrimitiveID cellID = it.first;
            Cell&       cell   = *( it.second );

            uint_t coarseLevel = level - 1;
            uint_t fineLevel   = level;

            const auto coarseDofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( cellID, coarseLevel );
            const auto fineDofMemory   = this->getDGFunction()->volumeDoFFunction()->dofMemory( cellID, fineLevel );

            for ( auto coarseCellType : celldof::allCellTypes )
            {
               for ( const auto& itSrc : celldof::macrocell::Iterator( coarseLevel, coarseCellType ) )
               {
                  const indexing::Index& coarseElementIdx = itSrc;

                  std::vector< hyteg::indexing::Index > fineElementIndices;
                  std::vector< celldof::CellType >      fineCellTypes;

                  real_t fineAverage = 0.0;
                  real_t fineInvAverage = 0.0;
                  real_t fineGeoAverage = 1.0;

                  volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                     coarseElementIdx, coarseCellType, fineElementIndices, fineCellTypes );

                  WALBERLA_CHECK_EQUAL( fineElementIndices.size(), fineCellTypes.size() );

                  uint_t nFineElements = fineElementIndices.size();

                  for ( uint_t fineIdx = 0; fineIdx < fineElementIndices.size(); fineIdx += 1 )
                  {
                     auto fineElementIdx = fineElementIndices[fineIdx];
                     auto fineCellType   = fineCellTypes[fineIdx];

                     real_t fineVal =
                        fineDofMemory[volumedofspace::indexing::index( fineElementIdx.x(),
                                                                        fineElementIdx.y(),
                                                                        fineElementIdx.z(),
                                                                        fineCellType,
                                                                        0u,
                                                                        1u,
                                                                        fineLevel,
                                                                        volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )];
                     
                     fineAverage += fineVal;
                     fineInvAverage += (1.0 / fineVal);
                     fineGeoAverage *= std::pow(fineVal, 1.0 / nFineElements);
                  }

                  // WALBERLA_LOG_INFO_ON_ROOT( "fineElementIndices.size() = " << fineElementIndices.size() );

                  coarseDofMemory[volumedofspace::indexing::index( coarseElementIdx.x(),
                                                                  coarseElementIdx.y(),
                                                                  coarseElementIdx.z(),
                                                                  coarseCellType,
                                                                  0u,
                                                                  1u,
                                                                  coarseLevel,
                                                                  volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )] =
                     // fineElementIndices.size() / fineInvAverage;
                     fineAverage / fineElementIndices.size();
                     // fineGeoAverage;
               }
            }
         }
      }
      else
      {
         // WALBERLA_ABORT("Not implemented");
         for ( auto& it : this->storage_->getFaces() )
         {
            PrimitiveID faceID = it.first;
            Face&       face   = *( it.second );

            uint_t coarseLevel = level - 1;
            uint_t fineLevel   = level;

            const auto coarseDofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( faceID, coarseLevel );
            const auto fineDofMemory   = this->getDGFunction()->volumeDoFFunction()->dofMemory( faceID, fineLevel );

            for ( auto coarseFaceType : facedof::allFaceTypes )
            {
               for ( const auto& itSrc : facedof::macroface::Iterator( coarseLevel, coarseFaceType ) )
               {
                  const indexing::Index& coarseElementIdx = itSrc;

                  std::vector< hyteg::indexing::Index > fineElementIndices;
                  std::vector< facedof::FaceType >      fineFaceTypes;

                  real_t fineAverage = 0.0;
                  real_t fineInvAverage = 0.0;

                  volumedofspace::indexing::getFineMicroElementsFromCoarseMicroElement(
                     coarseElementIdx, coarseFaceType, fineElementIndices, fineFaceTypes );

                  WALBERLA_CHECK_EQUAL( fineElementIndices.size(), fineFaceTypes.size() );

                  for ( uint_t fineIdx = 0; fineIdx < fineElementIndices.size(); fineIdx += 1 )
                  {
                     auto fineElementIdx = fineElementIndices[fineIdx];
                     auto fineFaceType   = fineFaceTypes[fineIdx];

                     real_t fineVal =
                        fineDofMemory[volumedofspace::indexing::index( fineElementIdx.x(),
                                                                        fineElementIdx.y(),
                                                                        fineFaceType,
                                                                        0u,
                                                                        1u,
                                                                        fineLevel,
                                                                        volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )];

                     fineAverage += fineVal;
                     fineInvAverage += (1.0 / fineVal);
                  }

                  // WALBERLA_LOG_INFO_ON_ROOT( "fineElementIndices.size() = " << fineElementIndices.size() );

                  coarseDofMemory[volumedofspace::indexing::index( coarseElementIdx.x(),
                                                                  coarseElementIdx.y(),
                                                                  coarseFaceType,
                                                                  0u,
                                                                  1u,
                                                                  coarseLevel,
                                                                  volumedofspace::indexing::VolumeDoFMemoryLayout::SoA )] =
                     fineElementIndices.size() / fineInvAverage;
                     // fineAverage / fineElementIndices.size();
               }
            }
         }
      }
   }

   void transferToAllLowerLevels( uint_t sourceLevel )
   {
      for ( uint_t level = sourceLevel; level > this->getMinLevel(); level-- )
      {
         transferToLowerLevel( level );
      }
   }

   inline real_t evaluateSampledAverage( std::array< Point3D, 3 > microTriangles, std::array< real_t, 3 > valueTriangles )
   {
      // WALBERLA_CHECK( nStrides > 1u, "nStrides must be greater than 1u" );

      Point3D microTet0 = microTriangles[0];
      Point3D microTet1 = microTriangles[1];
      Point3D microTet2 = microTriangles[2];

      real_t valueTet0 = valueTriangles[0];
      real_t valueTet1 = valueTriangles[1];
      real_t valueTet2 = valueTriangles[2];

      Point3D coordinates = (microTet0 + microTet1 + microTet2) / 3.0;

      auto xLocal = vertexdof::macroface::transformToLocalTri( microTet0, microTet1, microTet2, coordinates );

      auto value = valueTet0 * ( real_c( 1.0 ) - xLocal[0] - xLocal[1] ) + valueTet1 * xLocal[0] +
                  valueTet2 * xLocal[1];

      // return (value + valueTet0 + valueTet1 + valueTet2) / 4.0;
      return (valueTet0 + valueTet1 + valueTet2) / 3.0;
   }

   inline real_t evaluateSampledAverage( std::array< Point3D, 4 > microTets, std::array< real_t, 4 > valueTets )
   {
      // WALBERLA_CHECK( nStrides > 1u, "nStrides must be greater than 1u" );

      Point3D microTet0 = microTets[0];
      Point3D microTet1 = microTets[1];
      Point3D microTet2 = microTets[2];
      Point3D microTet3 = microTets[3];

      real_t valueTet0 = valueTets[0];
      real_t valueTet1 = valueTets[1];
      real_t valueTet2 = valueTets[2];
      real_t valueTet3 = valueTets[3];

      Point3D coordinates = (microTet0 + microTet1 + microTet2 + microTet3) / 4.0;

      // auto xLocal = vertexdof::macrocell::detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, coordinates );

      std::function< real_t(const Point3D&) > locallyEvaluate = [&](const Point3D& x)
      {
         return valueTet0 * ( real_c( 1.0 ) - x[0] - x[1] - x[2] ) + valueTet1 * x[0] +
                  valueTet2 * x[1] + valueTet3 * x[2];
      };

      Point3D qp1 = Point3D(0.25, 0.25, 0.25);
      Point3D qp2 = Point3D(0.16666667, 0.16666667, 0.5);
      Point3D qp3 = Point3D(0.16666667, 0.5, 0.16666667);
      Point3D qp4 = Point3D(0.5, 0.16666667, 0.16666667);
      Point3D qp5 = Point3D(0.16666667, 0.16666667, 0.16666667);

      real_t val1 = locallyEvaluate(qp1);
      real_t val2 = locallyEvaluate(qp2);
      real_t val3 = locallyEvaluate(qp3);
      real_t val4 = locallyEvaluate(qp4);
      real_t val5 = locallyEvaluate(qp5);

      // auto value = valueTet0 * ( real_c( 1.0 ) - xLocal[0] - xLocal[1] - xLocal[2] ) + valueTet1 * xLocal[0] +
      //             valueTet2 * xLocal[1] + valueTet3 * xLocal[2];

      return (val1 + val2 + val3 + val4 + val5) / 5.0;
      // return std::pow(val1, 1.0 / 5.0) * std::pow(val2, 1.0 / 5.0) 
      // * std::pow(val3, 1.0 / 5.0) * std::pow(val4, 1.0 / 5.0) * std::pow(val5, 1.0 / 5.0);
      // return 5.0 / ( (1.0 / value) + (1.0 / valueTet0) + (1.0 / valueTet1) + (1.0 / valueTet2) + (1.0 / valueTet3) );
      // return 4.0 / ( (1.0 / val1) + (1.0 / val2) + (1.0 / val3) + (1.0 / val4) );
      // return 5.0 / ( (1.0 / val1) + (1.0 / val2) + (1.0 / val3) + (1.0 / val4) + (1.0 / val5) );
   }

   void averageFromP1( P1Function< real_t > src, uint_t level )
   {
      if( this->storage_->hasGlobalCells() )
      {
         for ( auto it : this->storage_->getCells() )
         {
            PrimitiveID cellId = it.first;
            Cell&       cell   = *( it.second );

            const auto p0DofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( cellId, level );

            auto       p1FuncId   = src.getCellDataID();
            const auto p1FuncData = cell.getData( p1FuncId )->getPointer( level );

            for ( auto cellType : celldof::allCellTypes )
            {
               for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
               {
                  uint_t p0DofIdx = volumedofspace::indexing::index(
                     idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::SoA );

                  const std::array< indexing::Index, 4 > vertexIndices =
                     celldof::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );

                  auto microTet0 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[0] );
                  auto microTet1 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[1] );
                  auto microTet2 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[2] );
                  auto microTet3 = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[3] );

                  auto valueTet0 = p1FuncData[vertexdof::macrocell::index(
                     level, vertexIndices[0].x(), vertexIndices[0].y(), vertexIndices[0].z() )];
                  auto valueTet1 = p1FuncData[vertexdof::macrocell::index(
                     level, vertexIndices[1].x(), vertexIndices[1].y(), vertexIndices[1].z() )];
                  auto valueTet2 = p1FuncData[vertexdof::macrocell::index(
                     level, vertexIndices[2].x(), vertexIndices[2].y(), vertexIndices[2].z() )];
                  auto valueTet3 = p1FuncData[vertexdof::macrocell::index(
                     level, vertexIndices[3].x(), vertexIndices[3].y(), vertexIndices[3].z() )];

                  real_t sampledAverage = evaluateSampledAverage(
                     { microTet0, microTet1, microTet2, microTet3 }, { valueTet0, valueTet1, valueTet2, valueTet3 } );

                  // WALBERLA_LOG_INFO_ON_ROOT( "valueTet0 = " << valueTet0 );
                  // WALBERLA_LOG_INFO_ON_ROOT( "valueTet1 = " << valueTet1 );
                  // WALBERLA_LOG_INFO_ON_ROOT( "valueTet2 = " << valueTet2 );
                  // WALBERLA_LOG_INFO_ON_ROOT( "valueTet3 = " << valueTet3 );
                  // WALBERLA_LOG_INFO_ON_ROOT( "sampledAverage = " << sampledAverage );

                  p0DofMemory[p0DofIdx] = sampledAverage;
               }
            }
         }
      }
      else
      {
         for ( auto& it : this->getStorage()->getFaces() )
         {
            PrimitiveID faceId = it.first;
            Face&       face   = *( it.second );

            const auto p0DofMemory = this->getDGFunction()->volumeDoFFunction()->dofMemory( faceId, level );

            auto       p1FuncId   = src.getFaceDataID();
            const auto p1FuncData = face.getData( p1FuncId )->getPointer( level );

            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
               {
                  uint_t p0DofIdx = volumedofspace::indexing::index(
                     idxIt.x(), idxIt.y(), faceType, 0, 1, level, volumedofspace::indexing::VolumeDoFMemoryLayout::SoA );

                  const std::array< indexing::Index, 3 > vertexIndices =
                      facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
                  std::array< Point3D, 3 > elementVertices;
                  for ( uint_t i = 0; i < 3; i++ )
                  {
                     const auto elementVertex = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                     elementVertices[i]( 0 )  = elementVertex[0];
                     elementVertices[i]( 1 )  = elementVertex[1];
                     elementVertices[i]( 2 )  = 0.0;
                  }

                  auto valueTri0 = p1FuncData[vertexdof::macroface::indexFromVertex(
                     level, vertexIndices[0].x(), vertexIndices[0].y(), stencilDirection::VERTEX_C )];
                  auto valueTri1 = p1FuncData[vertexdof::macroface::indexFromVertex(
                     level, vertexIndices[1].x(), vertexIndices[1].y(), stencilDirection::VERTEX_C )];
                  auto valueTri2 = p1FuncData[vertexdof::macroface::indexFromVertex(
                     level, vertexIndices[2].x(), vertexIndices[2].y(), stencilDirection::VERTEX_C )];

                  real_t sampledAverage = evaluateSampledAverage(
                     elementVertices, { valueTri0, valueTri1, valueTri2 } );

                  p0DofMemory[p0DofIdx] = sampledAverage;
               }
            }
         }
      }
      
   }

   /// Set all function DoFs to zero including the ones in the halos
   void setToZero( const uint_t level ) const override final { dgFunction_->setToZero( level ); };

   void swap( const P0Function< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      WALBERLA_UNUSED( flag );
      dgFunction_->swap( *other.getDGFunction(), level, flag );
   };

   void copyFrom( const P0Function< ValueType >&         other,
                  const uint_t&                          level,
                  const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   };

   bool evaluate( const Point3D& coordinates, uint_t level, ValueType& value, real_t searchToleranceRadius = 1e-05 ) const
   {
      return dgFunction_->evaluate( coordinates, level, value, searchToleranceRadius );
   }

   void assign( const std::vector< ValueType >&                                               scalars,
                const std::vector< std::reference_wrapper< const P0Function< ValueType > > >& functions,
                uint_t                                                                        level,
                DoFType                                                                       flag = All ) const
   {
      std::vector< std::reference_wrapper< const DGFunction< ValueType > > > dgFunctions;
      for ( auto f : functions )
      {
         dgFunctions.push_back( *f.get().getDGFunction() );
      }
      dgFunction_->assign( scalars, dgFunctions, level, flag );
   }

   ValueType dotGlobal( const P0Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      return dgFunction_->dotGlobal( *rhs.getDGFunction(), level );
   }

   ValueType sumGlobal( uint_t level, const DoFType& flag = All ) const { return dgFunction_->sumGlobal( level ); }

   ValueType dotLocal( const P0Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      return dgFunction_->dotLocal( *rhs.getDGFunction(), level );
   }

   void enumerate( uint_t level ) const { dgFunction_->enumerate( level ); }

   void enumerate( uint_t level, ValueType& offset ) const { dgFunction_->enumerate( level, offset ); }

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const P0Function< idx_t >&            numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      dgFunction_->toVector( *numerator.getDGFunction(), vec, level, flag );
   }

   void fromVector( const P0Function< idx_t >&            numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const
   {
      dgFunction_->fromVector( *numerator.getDGFunction(), vec, level, flag );
   }
   /// @}

   uint_t getNumberOfLocalDoFs( uint_t level ) const { return dgFunction_->getNumberOfLocalDoFs( level ); }

   uint_t getNumberOfGlobalDoFs( uint_t          level,
                                 const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                 const bool&     onRootOnly   = false ) const
   {
      return dgFunction_->getNumberOfGlobalDoFs( level, communicator, onRootOnly );
   }

   template < typename OtherValueType >
   void copyBoundaryConditionFromFunction( const P0Function< OtherValueType >& other )
   {
      dgFunction_->copyBoundaryConditionFromFunction( *other.getDGFunction() );
   }

   ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const
   {
      if ( flag != All && flag != Inner )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::getMaxMagnitude -> DoFType flag will be ignored!" );
      }
      return dgFunction_->getMaxMagnitude( level, mpiReduce );
   }

   ValueType getMaxValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const
   {
      if ( flag != All && flag != Inner )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::getMaxValue -> DoFType flag will be ignored!" );
      }
      return dgFunction_->getMax( level, mpiReduce );
   }

   ValueType getMinValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const
   {
      if ( flag != All && flag != Inner )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "P0Function::getMinValue -> DoFType flag will be ignored!" );
      }
      return dgFunction_->getMin( level, mpiReduce );
   }

 private:
   std::shared_ptr< DGFunction< ValueType > > dgFunction_;
   bool                                       doNotWarnOnInterpolate_ = true;
};

namespace dg {
/// \brief removes the mean value of func to project it onto the space of functions with mean value 0
///
/// \param func             [in] function to project
/// \param level            [in] refinement level
inline void projectMean( P0Function< real_t >& func, const uint_t& level )
{
   const uint_t numGlobalVertices = func.getNumberOfGlobalDoFs( level );
   const real_t sum               = func.sumGlobal( level, Inner );
   func.add( -sum / ( real_c( numGlobalVertices ) ), level, Inner );
}
} // namespace dg

} // namespace hyteg
