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

#include <convection_particles/data/ParticleStorage.h>
#include <convection_particles/mpi/SyncNextNeighborsNoGhosts.h>
#include <core/math/MatrixMxN.h>
#include <core/mpi/MPIWrapper.h>

#include "hyteg/FunctionIterator.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/communication/convection_particles/SyncNextNeighborsByPrimitiveID.h"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2MacroCell.hpp"
#include "hyteg/p2functionspace/P2MacroFace.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/convection_particles_coupling/SetupPrimitiveStorageConvectionParticlesInterface.hpp"

namespace hyteg {

using namespace walberla::convection_particles;

enum class TimeSteppingScheme
{
   ExplicitEuler,
   RK4,
};

const static std::vector< std::vector< real_t > > RK_A_ExplicitEuler = {{{0}}};
const static std::vector< real_t >                RK_b_ExplicitEuler = {1};

const static std::vector< std::vector< real_t > > RK_A_RK4 = {{{0, 0, 0, 0}, {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}}};
const static std::vector< real_t >                RK_b_RK4 = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};

const static std::map< TimeSteppingScheme, std::vector< std::vector< real_t > > > RK_A = {
    {TimeSteppingScheme::ExplicitEuler, RK_A_ExplicitEuler},
    {TimeSteppingScheme::RK4, RK_A_RK4}};

const static std::map< TimeSteppingScheme, std::vector< real_t > > RK_b = {
    {TimeSteppingScheme::ExplicitEuler, RK_b_ExplicitEuler},
    {TimeSteppingScheme::RK4, RK_b_RK4}};

template < typename FunctionType >
inline Point3D stepRK( const Point3D&                              startPosition,
                       const real_t&                               velocityFactor,
                       const std::vector< std::vector< real_t > >& A,
                       const std::vector< real_t >&                b,
                       const FunctionType&                         ux,
                       const FunctionType&                         uy,
                       const FunctionType&                         uz,
                       const real_t&                               dt,
                       const uint_t&                               level )
{
   const uint_t          stages = b.size();
   std::vector< real_t > kx( stages );
   std::vector< real_t > ky( stages );
   std::vector< real_t > kz( stages );

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

template < typename FunctionType >
Point3D performInnerRKTimeSteps( const Point3D&            startingPosition,
                                 const uint_t&             innerSteps,
                                 const FunctionType&       ux,
                                 const FunctionType&       uy,
                                 const FunctionType&       uz,
                                 const uint_t&             level,
                                 const real_t&             dt,
                                 const TimeSteppingScheme& timeSteppingScheme )
{
   Point3D posPast = startingPosition;
   for ( uint_t step = 0; step < innerSteps; step++ )
   {
      posPast = stepRK< FunctionType >(
          posPast, -1.0, RK_A.at( timeSteppingScheme ), RK_b.at( timeSteppingScheme ), ux, uy, uz, dt, level );
   }
   return posPast;
}

template < typename FunctionType >
void integrateNodes( const PrimitiveStorage&,
                     const FunctionType&,
                     const FunctionType&,
                     const FunctionType&,
                     const FunctionType&,
                     const FunctionType&,
                     const real_t&,
                     const uint_t&,
                     const DoFType&,
                     const uint_t&,
                     const TimeSteppingScheme&,
                     const real_t& )
{
   WALBERLA_ABORT( "Integration not implemented for this function type." );
}

template <>
void integrateNodes( const PrimitiveStorage&     storage,
                     const P1Function< real_t >& c,
                     const P1Function< real_t >& cOld,
                     const P1Function< real_t >& ux,
                     const P1Function< real_t >& uy,
                     const P1Function< real_t >& uz,
                     const real_t&               dt,
                     const uint_t&               level,
                     const DoFType&,
                     const uint_t&             steps,
                     const TimeSteppingScheme& timeSteppingScheme,
                     const real_t& )
{
   WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1 );

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

      auto posPast =
          performInnerRKTimeSteps< P1Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

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

         auto posPast =
             performInnerRKTimeSteps< P1Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

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

         auto posPast =
             performInnerRKTimeSteps< P1Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         // evaluate new temperature
         const auto newTemp = cOld.evaluate( posPast, level );
         cData[idx]         = newTemp;
      }
   }
}

template <>
void integrateNodes( const PrimitiveStorage&     storage,
                     const P2Function< real_t >& c,
                     const P2Function< real_t >& cOld,
                     const P2Function< real_t >& ux,
                     const P2Function< real_t >& uy,
                     const P2Function< real_t >& uz,
                     const real_t&               dt,
                     const uint_t&               level,
                     const DoFType&,
                     const uint_t&             steps,
                     const TimeSteppingScheme& timeSteppingScheme,
                     const real_t&             initialOffset )
{
   WALBERLA_CHECK_EQUAL( walberla::mpi::MPIManager::instance()->numProcesses(), 1 );

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

      auto posPast =
          performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

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

         auto posPast =
             performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         const auto newTemp = cOld.evaluate( posPast, level );
         cDataV[idx]        = newTemp;
      }

      for ( const auto& it : edgedof::macroedge::Iterator( level, 0 ) )
      {
         const Point3D coordinate = edgedof::macroedge::coordinateFromIndex( level, edge, it );
         const uint_t  idx        = edgedof::macroedge::index( level, it.x() );

         auto posPast =
             performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         const auto newTemp = cOld.evaluate( posPast, level );
         cDataE[idx]        = newTemp;
      }
   }

   for ( const auto& oldFaceIt : storage.getFaces() )
   {
      const auto& face = *oldFaceIt.second;

      if ( storage.onBoundary( face.getID() ) )
         continue;

      auto cDataV = face.getData( c.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto cDataE = face.getData( c.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

      auto uxDataV = face.getData( ux.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto uxDataE = face.getData( ux.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

      auto uyDataV = face.getData( uy.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto uyDataE = face.getData( uy.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

      auto uzDataV = face.getData( uz.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto uzDataE = face.getData( uz.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

      for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
      {
         Point3D      coordinate = vertexdof::macroface::coordinateFromIndex( level, face, it );
         const uint_t idx        = vertexdof::macroface::index( level, it.x(), it.y() );

         if ( initialOffset > 0 )
         {
            auto velX = uxDataV[idx];
            auto velY = uyDataV[idx];
            auto velZ = uzDataV[idx];

            coordinate += initialOffset * Point3D( {velX, velY, velZ} );
         }

         auto posPast =
             performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         // evaluate new temperature
         const auto newTemp = cOld.evaluate( posPast, level );
         cDataV[idx]        = newTemp;
      }

      for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
      {
         if ( it.row() != 0 )
         {
            Point3D      coordinate = edgedof::macroface::coordinateFromIndex( level, face, it, edgedof::EdgeDoFOrientation::X );
            const uint_t idx        = edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::X );

            if ( initialOffset > 0 )
            {
               auto velX = uxDataE[idx];
               auto velY = uyDataE[idx];
               auto velZ = uzDataE[idx];

               coordinate += initialOffset * Point3D( {velX, velY, velZ} );
            }

            auto posPast =
                performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

            // evaluate new temperature
            const auto newTemp = cOld.evaluate( posPast, level );
            cDataE[idx]        = newTemp;
         }

         if ( it.col() != 0 )
         {
            Point3D      coordinate = edgedof::macroface::coordinateFromIndex( level, face, it, edgedof::EdgeDoFOrientation::Y );
            const uint_t idx        = edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::Y );

            if ( initialOffset > 0 )
            {
               auto velX = uxDataE[idx];
               auto velY = uyDataE[idx];
               auto velZ = uzDataE[idx];

               coordinate += initialOffset * Point3D( {velX, velY, velZ} );
            }

            auto posPast =
                performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

            // evaluate new temperature
            const auto newTemp = cOld.evaluate( posPast, level );
            cDataE[idx]        = newTemp;
         }

         if ( it.col() + it.row() != ( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
         {
            Point3D      coordinate = edgedof::macroface::coordinateFromIndex( level, face, it, edgedof::EdgeDoFOrientation::XY );
            const uint_t idx        = edgedof::macroface::index( level, it.x(), it.y(), edgedof::EdgeDoFOrientation::XY );

            if ( initialOffset > 0 )
            {
               auto velX = uxDataE[idx];
               auto velY = uyDataE[idx];
               auto velZ = uzDataE[idx];

               coordinate += initialOffset * Point3D( {velX, velY, velZ} );
            }

            auto posPast =
                performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

            // evaluate new temperature
            const auto newTemp = cOld.evaluate( posPast, level );
            cDataE[idx]        = newTemp;
         }
      }
   }

   for ( const auto& oldCellIt : storage.getCells() )
   {
      const auto& cell = *oldCellIt.second;

      auto cDataV = cell.getData( c.getVertexDoFFunction().getCellDataID() )->getPointer( level );
      auto cDataE = cell.getData( c.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

      auto uxDataV = cell.getData( ux.getVertexDoFFunction().getCellDataID() )->getPointer( level );
      auto uxDataE = cell.getData( ux.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

      auto uyDataV = cell.getData( uy.getVertexDoFFunction().getCellDataID() )->getPointer( level );
      auto uyDataE = cell.getData( uy.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

      auto uzDataV = cell.getData( uz.getVertexDoFFunction().getCellDataID() )->getPointer( level );
      auto uzDataE = cell.getData( uz.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

      for ( const auto& it : vertexdof::macrocell::Iterator( level, 1 ) )
      {
         Point3D      coordinate = vertexdof::macrocell::coordinateFromIndex( level, cell, it );
         const uint_t idx        = vertexdof::macrocell::index( level, it.x(), it.y(), it.z() );

         if ( initialOffset > 0 )
         {
            auto velX = uxDataV[idx];
            auto velY = uyDataV[idx];
            auto velZ = uzDataV[idx];

            coordinate += initialOffset * Point3D( {velX, velY, velZ} );
         }

         auto posPast =
             performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         // evaluate new temperature
         const auto newTemp = cOld.evaluate( posPast, level );
         cDataV[idx]        = newTemp;
      }

      for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) )
      {
         for ( const auto& orientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
         {
            if ( edgedof::macrocell::isInnerEdgeDoF( level, it, orientation ) )
            {
               Point3D      coordinate = edgedof::macrocell::coordinateFromIndex( level, cell, it, orientation );
               const uint_t idx        = edgedof::macrocell::index( level, it.x(), it.y(), it.z(), orientation );

               if ( initialOffset > 0 )
               {
                  auto velX = uxDataE[idx];
                  auto velY = uyDataE[idx];
                  auto velZ = uzDataE[idx];

                  coordinate += initialOffset * Point3D( {velX, velY, velZ} );
               }

               auto posPast = performInnerRKTimeSteps< P2Function< real_t > >(
                   coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

               // evaluate new temperature
               const auto newTemp = cOld.evaluate( posPast, level );
               cDataE[idx]        = newTemp;
            }
         }
      }

      for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) )
      {
         Point3D      coordinate = edgedof::macrocell::coordinateFromIndex( level, cell, it, edgedof::EdgeDoFOrientation::XYZ );
         const uint_t idx        = edgedof::macrocell::index( level, it.x(), it.y(), it.z(), edgedof::EdgeDoFOrientation::XYZ );

         if ( initialOffset > 0 )
         {
            auto velX = uxDataE[idx];
            auto velY = uyDataE[idx];
            auto velZ = uzDataE[idx];

            coordinate += initialOffset * Point3D( {velX, velY, velZ} );
         }

         auto posPast =
             performInnerRKTimeSteps< P2Function< real_t > >( coordinate, steps, ux, uy, uz, level, dt, timeSteppingScheme );

         // evaluate new temperature
         const auto newTemp = cOld.evaluate( posPast, level );
         cDataE[idx]        = newTemp;
      }
   }
}

void updateParticlePosition( const SetupPrimitiveStorage&                                     setupStorage,
                             walberla::convection_particles::data::ParticleStorage::Particle& particle,
                             const walberla::math::Vector3< real_t >&                         newPosition )
{
   std::set< PrimitiveID > containingPrimitives;
   auto                    pt = newPosition;
   if ( setupStorage.getNumberOfCells() == 0 )
   {
      for ( const auto& faceIt : setupStorage.getFaces() )
      {
         auto faceID = faceIt.first;
         auto face   = faceIt.second;

         Point2D pointOfInterest( {pt[0], pt[1]} );

         if ( circleTriangleIntersection( pointOfInterest,
                                          1e-05,
                                          Point2D( {face->getCoordinates().at( 0 )[0], face->getCoordinates().at( 0 )[1]} ),
                                          Point2D( {face->getCoordinates().at( 1 )[0], face->getCoordinates().at( 1 )[1]} ),
                                          Point2D( {face->getCoordinates().at( 2 )[0], face->getCoordinates().at( 2 )[1]} ) ) )
         {
            containingPrimitives.insert( faceID );
         }
      }
   }
   else
   {
      for ( const auto& cellIt : setupStorage.getCells() )
      {
         auto cellID = cellIt.first;
         auto cell   = cellIt.second;

         Point3D pointOfInterest( {pt[0], pt[1], pt[2]} );

         if ( sphereTetrahedronIntersection( pointOfInterest,
                                             1e-05,
                                             cell->getCoordinates().at( 0 ),
                                             cell->getCoordinates().at( 1 ),
                                             cell->getCoordinates().at( 2 ),
                                             cell->getCoordinates().at( 3 ) ) )
         {
            containingPrimitives.insert( cellID );
         }
      }
   }

   if ( containingPrimitives.empty() )
   {
      // do not change macro nor position
      // the particle left the domain or neighborhood
      // TODO: define behaviour here
      return;
   }
   else
   {
      // std::set is required to be sorted by the standard.
      // ::begin() returns the smallest element.
      particle->setPosition( newPosition );
      particle->setContainingPrimitive( *containingPrimitives.begin() );
   }
}

real_t evaluateAtParticlePosition( PrimitiveStorage&                                                      storage,
                                   const P2Function< real_t >&                                            function,
                                   const walberla::convection_particles::data::ParticleStorage::Particle& particle,
                                   const uint_t&                                                          level )
{
   if ( !storage.hasGlobalCells() )
   {
      WALBERLA_CHECK( storage.faceExistsLocally( particle.getContainingPrimitive() ) );
      Face& face = *storage.getFace( particle.getContainingPrimitive() );
      return P2::macroface::evaluate( level,
                                      face,
                                      toPoint3D( particle.getPosition() ),
                                      function.getVertexDoFFunction().getFaceDataID(),
                                      function.getEdgeDoFFunction().getFaceDataID() );
   }
   else
   {
      WALBERLA_CHECK( storage.cellExistsLocally( particle.getContainingPrimitive() ) );
      Cell& cell = *storage.getCell( particle.getContainingPrimitive() );
      return P2::macrocell::evaluate( level,
                                      cell,
                                      toPoint3D( particle.getPosition() ),
                                      function.getVertexDoFFunction().getCellDataID(),
                                      function.getEdgeDoFFunction().getCellDataID() );
   }
}

uint_t initializeParticles( walberla::convection_particles::data::ParticleStorage& particleStorage,
                            const SetupPrimitiveStorage&                           setupStorage,
                            PrimitiveStorage&                                      storage,
                            const P2Function< real_t >&                            c,
                            const P2Function< real_t >&                            ux,
                            const P2Function< real_t >&                            uy,
                            const P2Function< real_t >&                            uz,
                            const uint_t&                                          level,
                            const DoFType&,
                            const TimeSteppingScheme& timeSteppingScheme,
                            const real_t&             initialOffset )
{
   communication::syncFunctionBetweenPrimitives( ux, level );
   communication::syncFunctionBetweenPrimitives( uy, level );
   communication::syncFunctionBetweenPrimitives( uz, level );

   // initialize particles locally

   // store
   // - start position (coords)
   // - start index
   // - start process
   // - start macro ID
   // - start function (edgedof or vertexdof)
   // - start edgedof orientation
   // - initial velocity (optimization especially for single step, store in k[0], no evaluate call needed)
   //
   // also store number of created particles to verify the number of processes messages later

   particleStorage.clear();

   const uint_t                                rank     = uint_c( walberla::mpi::MPIManager::instance()->rank() );
   const std::vector< std::vector< real_t > >& A        = RK_A.at( timeSteppingScheme );
   const std::vector< real_t >&                b        = RK_b.at( timeSteppingScheme );
   const uint_t                                rkStages = b.size();

   for ( auto it : FunctionIterator< vertexdof::VertexDoFFunction< real_t > >( c.getVertexDoFFunction(), level ) )
   {
      if ( storage.onBoundary( it.primitiveID(), true ) )
         continue;

      auto particleIt = particleStorage.create();
      particleIt->setOwner( (int) rank );
      particleIt->setPosition( toVec3( it.coordinates() ) );
      particleIt->setStartPosition( toVec3( it.coordinates() ) );
      particleIt->setStartDoFType( 0 ); // 0 == vertexdof
      particleIt->setStartEdgeDoFOrientation( it.edgeDoFOrientation() );
      particleIt->setStartPrimitiveID( it.primitiveID() );
      particleIt->setStartIndex( it.index() );
      particleIt->setStartProcess( rank );
      particleIt->getKRef().resize( rkStages );

      PrimitiveID containingPrimitiveID = it.primitiveID();
      if ( !storage.hasGlobalCells() && !it.isOnMacroFace() )
      {
         // 2D
         containingPrimitiveID = *std::min_element( storage.getPrimitive( it.primitiveID() )->neighborFaces().begin(),
                                                    storage.getPrimitive( it.primitiveID() )->neighborFaces().end() );
      }
      else if ( storage.hasGlobalCells() && !it.isOnMacroCell() )
      {
         // 3D
         containingPrimitiveID = *std::min_element( storage.getPrimitive( it.primitiveID() )->neighborCells().begin(),
                                                    storage.getPrimitive( it.primitiveID() )->neighborCells().end() );
      }
      particleIt->setContainingPrimitive( containingPrimitiveID );

      if ( it.isOnMacroVertex() )
      {
         particleIt->getKRef()[0][0] = storage.getVertex( it.primitiveID() )
                                           ->getData( ux.getVertexDoFFunction().getVertexDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         particleIt->getKRef()[0][1] = storage.getVertex( it.primitiveID() )
                                           ->getData( uy.getVertexDoFFunction().getVertexDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         if ( storage.hasGlobalCells() )
            particleIt->getKRef()[0][2] = storage.getVertex( it.primitiveID() )
                                              ->getData( uz.getVertexDoFFunction().getVertexDataID() )
                                              ->getPointer( level )[it.arrayIndex()];
      }
      if ( it.isOnMacroEdge() )
      {
         particleIt->getKRef()[0][0] = storage.getEdge( it.primitiveID() )
                                           ->getData( ux.getVertexDoFFunction().getEdgeDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         particleIt->getKRef()[0][1] = storage.getEdge( it.primitiveID() )
                                           ->getData( uy.getVertexDoFFunction().getEdgeDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         if ( storage.hasGlobalCells() )
            particleIt->getKRef()[0][2] = storage.getEdge( it.primitiveID() )
                                              ->getData( uz.getVertexDoFFunction().getEdgeDataID() )
                                              ->getPointer( level )[it.arrayIndex()];
      }
      if ( it.isOnMacroFace() )
      {
         particleIt->getKRef()[0][0] = storage.getFace( it.primitiveID() )
                                           ->getData( ux.getVertexDoFFunction().getFaceDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         particleIt->getKRef()[0][1] = storage.getFace( it.primitiveID() )
                                           ->getData( uy.getVertexDoFFunction().getFaceDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         if ( storage.hasGlobalCells() )
            particleIt->getKRef()[0][2] = storage.getFace( it.primitiveID() )
                                              ->getData( uz.getVertexDoFFunction().getFaceDataID() )
                                              ->getPointer( level )[it.arrayIndex()];
      }
   }

   for ( auto it : FunctionIterator< EdgeDoFFunction< real_t > >( c.getEdgeDoFFunction(), level ) )
   {
      if ( storage.onBoundary( it.primitiveID(), true ) )
         continue;

      auto particleIt = particleStorage.create();
      particleIt->setOwner( (int) rank );
      particleIt->setPosition( toVec3( it.coordinates() ) );
      particleIt->setStartPosition( toVec3( it.coordinates() ) );
      particleIt->setStartDoFType( 1 ); // 1 == edgedof
      particleIt->setStartEdgeDoFOrientation( it.edgeDoFOrientation() );
      particleIt->setStartPrimitiveID( it.primitiveID() );
      particleIt->setStartIndex( it.index() );
      particleIt->setStartProcess( rank );
      particleIt->getKRef().resize( rkStages );

      PrimitiveID containingPrimitiveID = it.primitiveID();
      if ( !storage.hasGlobalCells() && !it.isOnMacroFace() )
      {
         // 2D
         containingPrimitiveID = *std::min_element( storage.getPrimitive( it.primitiveID() )->neighborFaces().begin(),
                                                    storage.getPrimitive( it.primitiveID() )->neighborFaces().end() );
      }
      else if ( storage.hasGlobalCells() && !it.isOnMacroCell() )
      {
         // 3D
         containingPrimitiveID = *std::min_element( storage.getPrimitive( it.primitiveID() )->neighborCells().begin(),
                                                    storage.getPrimitive( it.primitiveID() )->neighborCells().end() );
      }
      particleIt->setContainingPrimitive( containingPrimitiveID );

      WALBERLA_CHECK( !it.isOnMacroVertex() );
      if ( it.isOnMacroEdge() )
      {
         particleIt->getKRef()[0][0] = storage.getEdge( it.primitiveID() )
                                           ->getData( ux.getEdgeDoFFunction().getEdgeDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         particleIt->getKRef()[0][1] = storage.getEdge( it.primitiveID() )
                                           ->getData( uy.getEdgeDoFFunction().getEdgeDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         if ( storage.hasGlobalCells() )
            particleIt->getKRef()[0][2] = storage.getEdge( it.primitiveID() )
                                              ->getData( uz.getEdgeDoFFunction().getEdgeDataID() )
                                              ->getPointer( level )[it.arrayIndex()];
      }
      if ( it.isOnMacroFace() )
      {
         particleIt->getKRef()[0][0] = storage.getFace( it.primitiveID() )
                                           ->getData( ux.getEdgeDoFFunction().getFaceDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         particleIt->getKRef()[0][1] = storage.getFace( it.primitiveID() )
                                           ->getData( uy.getEdgeDoFFunction().getFaceDataID() )
                                           ->getPointer( level )[it.arrayIndex()];
         if ( storage.hasGlobalCells() )
            particleIt->getKRef()[0][2] = storage.getFace( it.primitiveID() )
                                              ->getData( uz.getEdgeDoFFunction().getFaceDataID() )
                                              ->getPointer( level )[it.arrayIndex()];
      }
   }

   const uint_t numberOfCreatedParticles = particleStorage.size();
   // WALBERLA_LOG_INFO( "Particles after creation: " << particleStorage.size() );

   // now sync particles so that all particles are "assigned" a macro-cell (macro-face in 2D) automatically
   // walberla::convection_particles::mpi::SyncNextNeighborsNoGhosts SNN;
   walberla::convection_particles::mpi::SyncNextNeighborsByPrimitiveID SNN;

   SNN( particleStorage, setupStorage );

   // WALBERLA_LOG_INFO( "Particles after init sync: " << particleStorage.size() );
   return numberOfCreatedParticles;
}

void particleIntegration( walberla::convection_particles::data::ParticleStorage& particleStorage,
                          const SetupPrimitiveStorage&                           setupStorage,
                          PrimitiveStorage&                                      storage,
                          const P2Function< real_t >&                            ux,
                          const P2Function< real_t >&                            uy,
                          const P2Function< real_t >&                            uz,
                          const real_t&                                          dt,
                          const uint_t&                                          level,
                          const DoFType&,
                          const uint_t&             steps,
                          const TimeSteppingScheme& timeSteppingScheme )
{
   communication::syncFunctionBetweenPrimitives( ux, level );
   communication::syncFunctionBetweenPrimitives( uy, level );
   communication::syncFunctionBetweenPrimitives( uz, level );

   walberla::convection_particles::mpi::SyncNextNeighborsByPrimitiveID SNN;

   const uint_t                                rank     = uint_c( walberla::mpi::MPIManager::instance()->rank() );
   const std::vector< std::vector< real_t > >& A        = RK_A.at( timeSteppingScheme );
   const std::vector< real_t >&                b        = RK_b.at( timeSteppingScheme );
   const uint_t                                rkStages = b.size();

   SNN( particleStorage, setupStorage );

   for ( uint_t step = 0; step < steps; step++ )
   {
      // WALBERLA_LOG_INFO_ON_ROOT( "Starting inner time step " << step << " ..." )
      // TODO: sort particle storage by governing macro (optimization, tbd)

      // RK stage 0
      // skip setting to start pos (already happened)
      // skip synchronization (already happened)

      for ( auto p : particleStorage )
      {
         p->getKRef()[0][0] = -evaluateAtParticlePosition( storage, ux, p, level );
         p->getKRef()[0][1] = -evaluateAtParticlePosition( storage, uy, p, level );
         if ( storage.hasGlobalCells() )
            p->getKRef()[0][2] = -evaluateAtParticlePosition( storage, uz, p, level );
      }

      // RK stages [1, ..., stages - 1]
      for ( uint_t stage = 1; stage < rkStages; stage++ )
      {
         // determine function evaluation points for this stage
         for ( auto p : particleStorage )
         {
            auto evaluationPoint = p.getStartPosition();
            for ( uint_t j = 0; j < stage; j++ )
            {
               evaluationPoint += dt * A[stage][j] * p->getK().at( j );
            }
            updateParticlePosition( setupStorage, p, evaluationPoint );
         }

         // sync particles to be able to evaluate the velocity at that point
         SNN( particleStorage, setupStorage );

         // evaluate velocity at current particle positions and update k[stage]
         for ( auto p : particleStorage )
         {
            p->getKRef()[stage][0] = -evaluateAtParticlePosition( storage, ux, p, level );
            p->getKRef()[stage][1] = -evaluateAtParticlePosition( storage, uy, p, level );
            if ( storage.hasGlobalCells() )
               p->getKRef()[stage][2] = -evaluateAtParticlePosition( storage, uz, p, level );
         }
      }

      // all k[i] are now calculated, set final integration result for each particle and
      // assign the position accordingly
      for ( auto p : particleStorage )
      {
         auto finalPosition = p->getStartPosition();
         for ( uint_t i = 0; i < rkStages; i++ )
         {
            finalPosition += dt * b[i] * p->getKRef()[i];
         }
         updateParticlePosition( setupStorage, p, finalPosition );
         p->setStartPosition( p->getPosition() );
      }

      // sync particles as position was finally updated
      SNN( particleStorage, setupStorage );
   }
}

void evaluateTemperature( walberla::convection_particles::data::ParticleStorage& particleStorage,
                          PrimitiveStorage&                                      storage,
                          const P2Function< real_t >&                            c,
                          const P2Function< real_t >&                            cOld,
                          const uint_t&                                          level,
                          const DoFType&,
                          const uint_t& numberOfCreatedParticles )
{
   communication::syncFunctionBetweenPrimitives( cOld, level );

   // evaluate temperature at final position
   for ( auto p : particleStorage )
   {
      const auto finalTemperature = evaluateAtParticlePosition( storage, cOld, p, level );
      p->setFinalTemperature( finalTemperature );
   }

   // Communicate temperatures in two steps:
   // 1. via MPI_ANY_SOURCE, send a dummy message (number of particles,
   //    allows for check it all msgs were received) to original process (where particle was created)
   // 2. use the buffer system as now as usual, the receiver knows the sender ranks by now

   // part I

   std::map< uint_t, unsigned > numParticlesToSendToRank;
   std::map< uint_t, unsigned > numParticlesToReceiveFromRank;

   std::map< uint_t, MPI_Request > sendRequests;

   for ( const auto& p : particleStorage )
   {
      if ( numParticlesToSendToRank.count( p.getStartProcess() ) == 0 )
         numParticlesToSendToRank[p.getStartProcess()] = 0;
      numParticlesToSendToRank[p.getStartProcess()]++;
      sendRequests[p.getStartProcess()] = MPI_Request();
   }

   for ( auto it : numParticlesToSendToRank )
   {
      // WALBERLA_LOG_INFO( "Particle communcation prep: rank " << rank << " -> " << it.first << ": " << it.second )
      MPI_Isend( &it.second,
                 1,
                 MPI_UNSIGNED,
                 (int) it.first,
                 0,
                 walberla::mpi::MPIManager::instance()->comm(),
                 &sendRequests[it.first] );
   }

   for ( auto it : sendRequests )
   {
      MPI_Status status;
      MPI_Wait( &it.second, &status );
   }

   uint_t numReceivedParticleLocations = 0;
   while ( numReceivedParticleLocations < numberOfCreatedParticles )
   {
      unsigned   numParticles;
      MPI_Status status;
      MPI_Recv(
          &numParticles, 1, MPI_UNSIGNED, MPI_ANY_SOURCE, MPI_ANY_TAG, walberla::mpi::MPIManager::instance()->comm(), &status );
      // WALBERLA_LOG_INFO( "Particle communcation prep: receiving " << numParticles << " particles from rank " << status.MPI_SOURCE );
      numReceivedParticleLocations += numParticles;
      numParticlesToReceiveFromRank[uint_c( status.MPI_SOURCE )] = numParticles;
      // WALBERLA_LOG_INFO( "total received particle infos: " << numReceivedParticleLocations )
      // WALBERLA_LOG_INFO( "created particles: " << numberOfCreatedParticles )
   }

   // part II
   std::set< walberla::mpi::MPIRank > ranksToReceiveFrom;
   for ( auto r : numParticlesToReceiveFromRank )
   {
      ranksToReceiveFrom.insert( (walberla::mpi::MPIRank) r.first );
   }
   walberla::mpi::BufferSystem bufferSystem( walberla::mpi::MPIManager::instance()->comm() );
   bufferSystem.setReceiverInfo( ranksToReceiveFrom, true );

   for ( const auto& p : particleStorage )
   {
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getStartPrimitiveID();
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getStartIndex();
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getStartDoFType();
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getStartEdgeDoFOrientation();
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getFinalTemperature();
   }

   bufferSystem.sendAll();

   for ( auto i = bufferSystem.begin(); i != bufferSystem.end(); ++i )
   {
      PrimitiveID                 primitiveID;
      indexing::Index             index;
      uint_t                      dofType;
      edgedof::EdgeDoFOrientation orientation;
      real_t                      temp;

      while ( !i.buffer().isEmpty() )
      {
         i.buffer() >> primitiveID;
         i.buffer() >> index;
         i.buffer() >> dofType;
         i.buffer() >> orientation;
         i.buffer() >> temp;

         WALBERLA_CHECK( storage.primitiveExistsLocally( primitiveID ) );
         if ( storage.vertexExistsLocally( primitiveID ) )
         {
            auto vertex = storage.getVertex( primitiveID );
            WALBERLA_CHECK_EQUAL( dofType, 0 );
            vertex->getData( c.getVertexDoFFunction().getVertexDataID() )->getPointer( level )[0] = temp;
         }
         else if ( storage.edgeExistsLocally( primitiveID ) )
         {
            auto edge = storage.getEdge( primitiveID );
            if ( dofType == 0 )
               edge->getData( c.getVertexDoFFunction().getEdgeDataID() )
                   ->getPointer( level )[vertexdof::macroedge::index( level, index.x() )] = temp;
            else
               edge->getData( c.getEdgeDoFFunction().getEdgeDataID() )
                   ->getPointer( level )[edgedof::macroedge::index( level, index.x() )] = temp;
         }
         else if ( storage.faceExistsLocally( primitiveID ) )
         {
            auto face = storage.getFace( primitiveID );
            if ( dofType == 0 )
               face->getData( c.getVertexDoFFunction().getFaceDataID() )
                   ->getPointer( level )[vertexdof::macroface::index( level, index.x(), index.y() )] = temp;
            else
               face->getData( c.getEdgeDoFFunction().getFaceDataID() )
                   ->getPointer( level )[edgedof::macroface::index( level, index.x(), index.y(), orientation )] = temp;
         }
         else if ( storage.cellExistsLocally( primitiveID ) )
         {
            auto cell = storage.getCell( primitiveID );
            if ( dofType == 0 )
               cell->getData( c.getVertexDoFFunction().getCellDataID() )
                   ->getPointer( level )[vertexdof::macrocell::index( level, index.x(), index.y(), index.z() )] = temp;
            else
               cell->getData( c.getEdgeDoFFunction().getCellDataID() )
                   ->getPointer( level )[edgedof::macrocell::index( level, index.x(), index.y(), index.z(), orientation )] = temp;
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
   , cTmp_( "cTmp", storage, minLevel, maxLevel )
   , cPlus_( "cPlus", storage, minLevel, maxLevel )
   , cMinus_( "cMinus", storage, minLevel, maxLevel )
   , cAdjusted_( "cAdjusted", storage, minLevel, maxLevel )
   , timeSteppingSchemeConvection_( timeSteppingSchemeConvection )
   , numberOfCreatedParticles_( 0 )
   , particleStorage_( 10000 )
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
      cOld_.assign( {1.0}, {c}, level, All );

      integrateNodes< FunctionType >(
          *storage_, c, cOld_, ux, uy, uz, dt, level, flag, innerSteps, timeSteppingSchemeConvection_, 0 );
   }

   void step( const std::shared_ptr< SetupPrimitiveStorage >& setupStorage,
              const FunctionType&                             c,
              const FunctionType&                             ux,
              const FunctionType&                             uy,
              const FunctionType&                             uz,
              const uint_t&                                   level,
              const DoFType&                                  flag,
              const real_t&                                   dt,
              const uint_t&                                   innerSteps,
              const bool&                                     resetParticles = true )
   {
      if ( resetParticles )
      {
         cOld_.assign( {1.0}, {c}, level, All );
         numberOfCreatedParticles_ = initializeParticles(
             particleStorage_, *setupStorage, *storage_, c, ux, uy, uz, level, Inner, timeSteppingSchemeConvection_, 0 );
      }

      particleIntegration(
          particleStorage_, *setupStorage, *storage_, ux, uy, uz, dt, level, Inner, innerSteps, timeSteppingSchemeConvection_ );
      evaluateTemperature( particleStorage_, *storage_, c, cOld_, level, Inner, numberOfCreatedParticles_ );
   }

   template < typename MassOperator >
   void step( const FunctionType& c,
              const FunctionType& ux,
              const FunctionType& uy,
              const FunctionType& uz,
              const uint_t&       level,
              const DoFType&      flag,
              const real_t&       dt,
              const uint_t&       innerSteps,
              const MassOperator& massOperator,
              const real_t&       allowedRelativeMassDifference,
              const real_t&       adjustedAdvectionOffset )
   {
      cOld_.assign( {1.0}, {c}, level, All );

      // calculate old mass
      massOperator.apply( cOld_, cTmp_, level, flag );
      auto massBefore = cTmp_.sumGlobal( level, flag );

      integrateNodes< FunctionType >(
          *storage_, c, cOld_, ux, uy, uz, dt, level, flag, innerSteps, timeSteppingSchemeConvection_, 0 );

      // calculate new mass
      massOperator.apply( c, cTmp_, level, flag );
      auto massAfter = cTmp_.sumGlobal( level, flag );

      auto relativeMassDifference = std::abs( ( massAfter - massBefore ) / massBefore );

      if ( relativeMassDifference <= allowedRelativeMassDifference )
         return;

      // perform adjusted advection steps
      integrateNodes< FunctionType >( *storage_,
                                      cPlus_,
                                      cOld_,
                                      ux,
                                      uy,
                                      uz,
                                      dt,
                                      level,
                                      flag,
                                      innerSteps,
                                      timeSteppingSchemeConvection_,
                                      adjustedAdvectionOffset );
      integrateNodes< FunctionType >( *storage_,
                                      cMinus_,
                                      cOld_,
                                      ux,
                                      uy,
                                      uz,
                                      dt,
                                      level,
                                      flag,
                                      innerSteps,
                                      timeSteppingSchemeConvection_,
                                      -adjustedAdvectionOffset );

      // max/min assign functions
      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > maxAssignment =
          []( const Point3D&, const std::vector< real_t >& values ) { return std::max( values[0], values[1] ); };

      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > minAssignment =
          []( const Point3D&, const std::vector< real_t >& values ) { return std::min( values[0], values[1] ); };

      if ( massAfter <= massBefore )
      {
         cAdjusted_.interpolate( maxAssignment, {cPlus_, cMinus_}, level );
      }
      else
      {
         cAdjusted_.interpolate( minAssignment, {cPlus_, cMinus_}, level );
      }

      // calculate adjustment mass
      massOperator.apply( cAdjusted_, cTmp_, level, flag );
      auto massAdjusted = cTmp_.sumGlobal( level, flag );

      auto theta = ( massBefore - massAdjusted ) / ( massAfter - massAdjusted );

      c.assign( {theta, 1 - theta}, {c, cAdjusted_}, level, flag );
   }

 private:
   const std::shared_ptr< PrimitiveStorage >             storage_;
   FunctionType                                          cOld_;
   FunctionType                                          cTmp_;
   FunctionType                                          cPlus_;
   FunctionType                                          cMinus_;
   FunctionType                                          cAdjusted_;
   TimeSteppingScheme                                    timeSteppingSchemeConvection_;
   uint_t                                                numberOfCreatedParticles_;
   walberla::convection_particles::data::ParticleStorage particleStorage_;
};

} // namespace hyteg
