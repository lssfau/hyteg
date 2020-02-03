/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include <core/math/MatrixMxN.h>

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

enum class TimeSteppingScheme
{
   ExplicitEuler,
   RK4,
};

std::array< std::array< real_t, 1 >, 1 > RK_A_ExplicitEuler = {{{0}}};
std::array< real_t, 1 >                  RK_b_ExplicitEuler = {1};

std::array< std::array< real_t, 4 >, 4 > RK_A_RK4 = {{{0, 0, 0, 0}, {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}}};
std::array< real_t, 4 >                  RK_b_RK4 = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};

template < typename FunctionType, uint_t stages >
inline Point3D stepRK( const Point3D&                                            startPosition,
                       const real_t&                                             velocityFactor,
                       const std::array< std::array< real_t, stages >, stages >& A,
                       const std::array< real_t, stages >&                       b,
                       const FunctionType&                                       ux,
                       const FunctionType&                                       uy,
                       const FunctionType&                                       uz,
                       const real_t&                                             dt,
                       const uint_t&                                             level )
{
   std::array< real_t, stages > kx;
   std::array< real_t, stages > ky;
   std::array< real_t, stages > kz;

   kx[0] = velocityFactor * ux.evaluate( startPosition, level );
   ky[0] = velocityFactor * uy.evaluate( startPosition, level );
   kz[0] = velocityFactor * uz.evaluate( startPosition, level );

   for ( uint_t i = 1; i < stages; i++ )
   {
      Point3D evaluationPoint = startPosition;
      for ( uint_t j = 0; j < i; j++ )
      {
         evaluationPoint += dt * A[i][j] * Point3D( {kx[j], ky[j], kz[j]} );
      }
      kx[i] = velocityFactor * ux.evaluate( evaluationPoint, level );
      ky[i] = velocityFactor * uy.evaluate( evaluationPoint, level );
      kz[i] = velocityFactor * uz.evaluate( evaluationPoint, level );
   }

   real_t kxSum = 0;
   real_t kySum = 0;
   real_t kzSum = 0;
   for ( uint_t i = 0; i < stages; i++ )
   {
      kxSum += b[i] * kx[i];
      kySum += b[i] * ky[i];
      kzSum += b[i] * kz[i];
   }

   return startPosition + dt * Point3D( {kxSum, kySum, kzSum} );
}

template< typename FunctionType >
Point3D performInnerRKTimeSteps( const Point3D&      startingPosition,
                               const uint_t&       innerSteps,
                               const FunctionType& ux,
                               const FunctionType& uy,
                               const FunctionType& uz,
                               const uint_t&       level,
                               const real_t&       dt,
                               const TimeSteppingScheme & timeSteppingScheme )
{
   Point3D posPast = startingPosition;
   for ( uint_t step = 0; step < innerSteps; step++ )
   {
      if ( timeSteppingScheme == TimeSteppingScheme::ExplicitEuler )
      {
         posPast = stepRK< FunctionType, 1 >( posPast, -1.0, RK_A_ExplicitEuler, RK_b_ExplicitEuler, ux, uy, uz, dt, level );
      }
      else if ( timeSteppingScheme == TimeSteppingScheme::RK4 )
      {
         posPast = stepRK< FunctionType, 4 >( posPast, -1.0, RK_A_RK4, RK_b_RK4, ux, uy, uz, dt, level );
      }
   }
   return posPast;
}

template< typename FunctionType >
void integrateNodes( const PrimitiveStorage& ,
                     const FunctionType&     ,
                     const FunctionType&     ,
                     const FunctionType&     ,
                     const FunctionType&     ,
                     const FunctionType&     ,
                     const real_t&           ,
                     const uint_t&           ,
                     const DoFType&,
                     const uint_t&,
                     const TimeSteppingScheme & )
{
   WALBERLA_ABORT( "Integration not implemented for this function type." );
}

template<>
void integrateNodes( const PrimitiveStorage& storage,
                     const P1Function< real_t >&     c,
                     const P1Function< real_t >&     cOld,
                     const P1Function< real_t >&     ux,
                     const P1Function< real_t >&     uy,
                     const P1Function< real_t >&     uz,
                     const real_t&           dt,
                     const uint_t&           level,
                     const DoFType&,
                     const uint_t& steps,
                     const TimeSteppingScheme & timeSteppingScheme )
{
   WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1 );

   cOld.assign( {1.0}, {c}, level, All );

   communication::syncFunctionBetweenPrimitives( cOld, level );
   communication::syncFunctionBetweenPrimitives( ux, level );
   communication::syncFunctionBetweenPrimitives( uy, level );
   communication::syncFunctionBetweenPrimitives( uz, level );

   for ( const auto& oldVertexIt : storage.getVertices() )
   {
      const auto& vertex = *oldVertexIt.second;

      if ( storage.onBoundary( vertex.getID() ) )
         continue;

      auto cData = vertex.getData( c.getVertexDataID() )->getPointer( level );

      const Point3D coordinate = vertex.getCoordinates();
      const uint_t  idx        = 0;

      auto posPast = performInnerRKTimeSteps< P1Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

      // evaluate new temperature
      const auto newTemp = cOld.evaluate( posPast, level );
      cData[idx]         = newTemp;
   }

   for ( const auto& oldEdgeIt : storage.getEdges() )
   {
      const auto& edge = *oldEdgeIt.second;

      if ( storage.onBoundary( edge.getID() ) )
         continue;

      auto cData = edge.getData( c.getEdgeDataID() )->getPointer( level );

      for ( const auto& it : vertexdof::macroedge::Iterator( level, 1 ) )
      {
         const Point3D coordinate = vertexdof::macroedge::coordinateFromIndex( level, edge, it );
         const uint_t  idx        = vertexdof::macroedge::index( level, it.x() );

         auto posPast = performInnerRKTimeSteps< P1Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         const auto newTemp = cOld.evaluate( posPast, level );
         cData[idx]         = newTemp;
      }
   }

   for ( const auto& oldFaceIt : storage.getFaces() )
   {
      const auto& face = *oldFaceIt.second;

      auto cData = face.getData( c.getFaceDataID() )->getPointer( level );

      for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
      {
         const Point3D coordinate = vertexdof::macroface::coordinateFromIndex( level, face, it );
         const uint_t  idx        = vertexdof::macroface::index( level, it.x(), it.y() );

         auto posPast = performInnerRKTimeSteps< P1Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         // evaluate new temperature
         const auto newTemp = cOld.evaluate( posPast, level );
         cData[idx]         = newTemp;
      }
   }
}


template<>
void integrateNodes( const PrimitiveStorage& storage,
                     const P2Function< real_t >&     c,
                     const P2Function< real_t >&     cOld,
                     const P2Function< real_t >&     ux,
                     const P2Function< real_t >&     uy,
                     const P2Function< real_t >&     uz,
                     const real_t&           dt,
                     const uint_t&           level,
                     const DoFType&,
                     const uint_t& steps,
                     const TimeSteppingScheme & timeSteppingScheme )
{
   WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1 );

   cOld.assign( {1.0}, {c}, level, All );

   communication::syncFunctionBetweenPrimitives( cOld, level );
   communication::syncFunctionBetweenPrimitives( ux, level );
   communication::syncFunctionBetweenPrimitives( uy, level );
   communication::syncFunctionBetweenPrimitives( uz, level );

   for ( const auto& oldVertexIt : storage.getVertices() )
   {
      const auto& vertex = *oldVertexIt.second;

      if ( storage.onBoundary( vertex.getID() ) )
         continue;

      auto cData = vertex.getData( c.getVertexDoFFunction().getVertexDataID() )->getPointer( level );

      const Point3D coordinate = vertex.getCoordinates();
      const uint_t  idx        = 0;

      auto posPast = performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

      // evaluate new temperature
      const auto newTemp = cOld.evaluate( posPast, level );
      cData[idx]         = newTemp;
   }

   for ( const auto& oldEdgeIt : storage.getEdges() )
   {
      const auto& edge = *oldEdgeIt.second;

      if ( storage.onBoundary( edge.getID() ) )
         continue;

      auto cDataV = edge.getData( c.getVertexDoFFunction().getEdgeDataID() )->getPointer( level );
      auto cDataE = edge.getData( c.getEdgeDoFFunction().getEdgeDataID() )->getPointer( level );

      for ( const auto& it : vertexdof::macroedge::Iterator( level, 1 ) )
      {
         const Point3D coordinate = vertexdof::macroedge::coordinateFromIndex( level, edge, it );
         const uint_t  idx        = vertexdof::macroedge::index( level, it.x() );

         auto posPast = performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         const auto newTemp = cOld.evaluate( posPast, level );
         cDataV[idx]         = newTemp;
      }

      for ( const auto& it : edgedof::macroedge::Iterator( level, 0 ) )
      {
         const Point3D coordinate = edgedof::macroedge::coordinateFromIndex( level, edge, it );
         const uint_t  idx        = edgedof::macroedge::index( level, it.x() );

         auto posPast = performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         const auto newTemp = cOld.evaluate( posPast, level );
         cDataE[idx]         = newTemp;
      }
   }

   for ( const auto& oldFaceIt : storage.getFaces() )
   {
      const auto& face = *oldFaceIt.second;

      auto cDataV = face.getData( c.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto cDataE = face.getData( c.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

      for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
      {
         const Point3D coordinate = vertexdof::macroface::coordinateFromIndex( level, face, it );
         const uint_t  idx        = vertexdof::macroface::index( level, it.x(), it.y() );

         auto posPast = performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         // evaluate new temperature
         const auto newTemp = cOld.evaluate( posPast, level );
         cDataV[idx]         = newTemp;
      }

      for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
      {
         if ( it.row() != 0 )
         {
            const Point3D coordinate = edgedof::macroface::coordinateFromIndex( level, face, it, edgedof::EdgeDoFOrientation::X );
            const uint_t  idx        = edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::X );

            auto posPast = performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

            // evaluate new temperature
            const auto newTemp = cOld.evaluate( posPast, level );
            cDataE[idx]         = newTemp;
         }

         if ( it.col() != 0 )
         {
            const Point3D coordinate = edgedof::macroface::coordinateFromIndex( level, face, it, edgedof::EdgeDoFOrientation::Y );
            const uint_t  idx        = edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::Y );

            auto posPast = performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

            // evaluate new temperature
            const auto newTemp = cOld.evaluate( posPast, level );
            cDataE[idx]         = newTemp;
         }

         if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
         {
            const Point3D coordinate = edgedof::macroface::coordinateFromIndex( level, face, it, edgedof::EdgeDoFOrientation::XY );
            const uint_t  idx        = edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::XY );

            auto posPast = performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

            // evaluate new temperature
            const auto newTemp = cOld.evaluate( posPast, level );
            cDataE[idx]         = newTemp;
         }
      }
   }
}


template < typename FunctionType >
class MMOCTransport
{
 public:
   MMOCTransport( const std::shared_ptr< PrimitiveStorage >& storage,
                  const uint_t                               minLevel,
                  const uint_t                               maxLevel,
                  const TimeSteppingScheme&                  timeSteppingSchemeConvection )
   : storage_( storage )
   , cOld_( "cOld", storage, minLevel, maxLevel )
   , timeSteppingSchemeConvection_( timeSteppingSchemeConvection )
   {}

   void step( const FunctionType& c,
              const FunctionType& ux,
              const FunctionType& uy,
              const FunctionType& uz,
              const uint_t&       level,
              const DoFType&      flag,
              const real_t&       dt,
              const uint_t&       innerSteps )
   {
      integrateNodes< FunctionType >( *storage_, c, cOld_, ux, uy, uz, dt, level, flag, innerSteps, timeSteppingSchemeConvection_ );
   }

 private:

   const std::shared_ptr< PrimitiveStorage > storage_;
   FunctionType                              cOld_;
   TimeSteppingScheme                        timeSteppingSchemeConvection_;
};

} // namespace hyteg
