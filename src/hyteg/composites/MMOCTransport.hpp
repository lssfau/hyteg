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
using namespace walberla::mpistubs;

enum class TimeSteppingScheme
{
   ExplicitEuler, // first order
   RK3,           // third order
   Ralston,       // third order
   RK4,           // fourth order
};

const static std::vector< std::vector< real_t > > RK_A_ExplicitEuler = {{{0}}};
const static std::vector< real_t >                RK_b_ExplicitEuler = {1};

const static std::vector< std::vector< real_t > > RK_A_RK3 = {{{0, 0, 0}, {0.5, 0, 0}, {-1., 2., 0}}};
const static std::vector< real_t >                RK_b_RK3 = {1. / 6., 2. / 3., 1. / 6.};

const static std::vector< std::vector< real_t > > RK_A_Ralston = {{{0, 0, 0}, {0.5, 0, 0}, {0, 0.75, 0}}};
const static std::vector< real_t >                RK_b_Ralston = {2. / 9., 1. / 3., 4. / 9.};

const static std::vector< std::vector< real_t > > RK_A_RK4 = {{{0, 0, 0, 0}, {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}}};
const static std::vector< real_t >                RK_b_RK4 = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};

const static std::map< TimeSteppingScheme, std::vector< std::vector< real_t > > > RK_A = {
    {TimeSteppingScheme::ExplicitEuler, RK_A_ExplicitEuler},
    {TimeSteppingScheme::RK3, RK_A_RK3},
    {TimeSteppingScheme::Ralston, RK_A_Ralston},
    {TimeSteppingScheme::RK4, RK_A_RK4}};

const static std::map< TimeSteppingScheme, std::vector< real_t > > RK_b = {
    {TimeSteppingScheme::ExplicitEuler, RK_b_ExplicitEuler},
    {TimeSteppingScheme::RK3, RK_b_RK3},
    {TimeSteppingScheme::Ralston, RK_b_Ralston},
    {TimeSteppingScheme::RK4, RK_b_RK4}};

std::vector< PrimitiveID > getLocalAndNeighboringPrimitives( const SetupPrimitiveStorage & setupStorage )
{
   std::set< PrimitiveID > localAndNeighboringPrimitives;
   if ( setupStorage.getNumberOfCells() > 0 )
   {
      for ( const auto& it : setupStorage.getCells() )
      {
         auto cellID = it.first;
         auto cell   = it.second;
         if ( setupStorage.getTargetRank( cellID ) != uint_c( walberla::mpi::MPIManager::instance()->rank() ) )
            continue;
         localAndNeighboringPrimitives.insert( cellID );
         for ( const auto& vertexID : cell->neighborVertices() )
         {
            auto vertex = setupStorage.getVertex( vertexID );
            for ( const auto& neighborCellID : vertex->neighborCells() )
            {
               if ( setupStorage.getTargetRank( neighborCellID ) == uint_c( walberla::mpi::MPIManager::instance()->rank() ) )
                  continue;
               localAndNeighboringPrimitives.insert( neighborCellID );
            }
         }
      }
   }
   else
   {
      for ( const auto& it : setupStorage.getFaces() )
      {
         auto faceID = it.first;
         auto face   = it.second;
         if ( setupStorage.getTargetRank( faceID ) != uint_c( walberla::mpi::MPIManager::instance()->rank() ) )
            continue;
         localAndNeighboringPrimitives.insert( faceID );
         for ( const auto& vertexID : face->neighborVertices() )
         {
            auto vertex = setupStorage.getVertex( vertexID );
            for ( const auto& neighborFaceID : vertex->neighborFaces() )
            {
               if ( setupStorage.getTargetRank( neighborFaceID ) == uint_c( walberla::mpi::MPIManager::instance()->rank() ) )
                  continue;
               localAndNeighboringPrimitives.insert( neighborFaceID );
            }
         }
      }
   }
   return std::vector< PrimitiveID >( localAndNeighboringPrimitives.begin(), localAndNeighboringPrimitives.end() );
}

void updateParticlePosition( const SetupPrimitiveStorage&                           setupStorage,
                             const std::vector< PrimitiveID > & localAndNeighboringPrimitives,
                             walberla::convection_particles::data::ParticleStorage& particleStorage )
{
   std::set< PrimitiveID > containingPrimitives;
   if ( setupStorage.getNumberOfCells() == 0 )
   {
      for ( const auto& faceID : localAndNeighboringPrimitives )
      {
         auto face   = setupStorage.getFace( faceID );

         for ( auto p : particleStorage )
         {
            Point2D pointOfInterest( {p->getPosition()[0], p->getPosition()[1]} );

            if ( isPointInTriangle( pointOfInterest,
                                    Point2D( {face->getCoordinates().at( 0 )[0], face->getCoordinates().at( 0 )[1]} ),
                                    Point2D( {face->getCoordinates().at( 1 )[0], face->getCoordinates().at( 1 )[1]} ),
                                    Point2D( {face->getCoordinates().at( 2 )[0], face->getCoordinates().at( 2 )[1]} ) ) )
            {
               p->setContainingPrimitive( faceID );
            }
         }
      }
   }
   else
   {
      for ( const auto& cellID : localAndNeighboringPrimitives )
      {
         auto cell   = setupStorage.getCell( cellID );

         for ( auto p : particleStorage )
         {
            auto pointOfInterest = toPoint3D( p->getPosition() );

            if ( isPointInTetrahedron( pointOfInterest,
                                       cell->getCoordinates().at( 0 ),
                                       cell->getCoordinates().at( 1 ),
                                       cell->getCoordinates().at( 2 ),
                                       cell->getCoordinates().at( 3 ) ) )
            {
               p->setContainingPrimitive( cellID );
            }
         }
      }
   }
}

real_t evaluateAtParticlePosition( PrimitiveStorage&                                                      storage,
                                   const P2Function< real_t >&                                            function,
                                   const walberla::convection_particles::data::ParticleStorage::Particle& particle,
                                   const uint_t&                                                          level )
{
   real_t result;
   if ( !storage.hasGlobalCells() )
   {
      WALBERLA_CHECK( storage.faceExistsLocally( particle.getContainingPrimitive() ) );
      Face& face = *storage.getFace( particle.getContainingPrimitive() );
      result     = P2::macroface::evaluate( level,
                                        face,
                                        toPoint3D( particle.getPosition() ),
                                        function.getVertexDoFFunction().getFaceDataID(),
                                        function.getEdgeDoFFunction().getFaceDataID() );
   }
   else
   {
      WALBERLA_CHECK( storage.cellExistsLocally( particle.getContainingPrimitive() ) );
      Cell& cell = *storage.getCell( particle.getContainingPrimitive() );
      result     = P2::macrocell::evaluate( level,
                                        cell,
                                        toPoint3D( particle.getPosition() ),
                                        function.getVertexDoFFunction().getCellDataID(),
                                        function.getEdgeDoFFunction().getCellDataID() );
   }
   return result;
}

void evaluateAtParticlePosition( PrimitiveStorage&                                                      storage,
                                 const std::vector< P2Function< real_t > >&                             functions,
                                 const walberla::convection_particles::data::ParticleStorage::Particle& particle,
                                 const uint_t&                                                          level,
                                 std::vector< real_t >&                                                 results )
{
   if ( !storage.hasGlobalCells() )
   {
      WALBERLA_CHECK( storage.faceExistsLocally( particle.getContainingPrimitive() ) );
      Face& face = *storage.getFace( particle.getContainingPrimitive() );

      for ( uint_t i = 0; i < functions.size(); i++ )
      {
         results[i] = P2::macroface::evaluate( level,
                                               face,
                                               toPoint3D( particle.getPosition() ),
                                               functions[i].getVertexDoFFunction().getFaceDataID(),
                                               functions[i].getEdgeDoFFunction().getFaceDataID() );
      }
   }
   else
   {
      std::vector< PrimitiveDataID< FunctionMemory< real_t >, Cell > > vertexDataIDs;
      std::vector< PrimitiveDataID< FunctionMemory< real_t >, Cell > > edgeDataIDs;

      for ( uint_t i = 0; i < functions.size(); i++ )
      {
         vertexDataIDs.push_back( functions[i].getVertexDoFFunction().getCellDataID() );
         edgeDataIDs.push_back( functions[i].getEdgeDoFFunction().getCellDataID() );
      }

      WALBERLA_CHECK( storage.cellExistsLocally( particle.getContainingPrimitive() ) );
      Cell& cell = *storage.getCell( particle.getContainingPrimitive() );

      P2::macrocell::evaluate( level, cell, toPoint3D( particle.getPosition() ), vertexDataIDs, edgeDataIDs, results );
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
                          const std::vector< PrimitiveID > &                     localAndNeighboringPrimitives,
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

   storage.getTimingTree()->start( "Sync particles" );
   SNN( particleStorage, setupStorage );
   storage.getTimingTree()->stop( "Sync particles" );

   for ( uint_t step = 0; step < steps; step++ )
   {
      // WALBERLA_LOG_INFO_ON_ROOT( "Starting inner time step " << step << " ..." )
      // TODO: sort particle storage by governing macro (optimization, tbd)

      // RK stage 0
      // skip setting to start pos (already happened)
      // skip synchronization (already happened)

      storage.getTimingTree()->start( "Evaluate at particle position" );

      std::vector< real_t > results( {0, 0} );
      std::vector< P2Function< real_t > > functions = { ux, uy };
      if ( storage.hasGlobalCells() )
      {
         results.push_back( 0 );
         functions.push_back( uz );
      }

      for ( auto p : particleStorage )
      {
         evaluateAtParticlePosition( storage, functions, p, level, results );
         p->getKRef()[0][0] = - results[0];
         p->getKRef()[0][1] = - results[1];
         if ( storage.hasGlobalCells() )
            p->getKRef()[0][2] = - results[2];
      }
      storage.getTimingTree()->stop( "Evaluate at particle position" );

      // RK stages [1, ..., stages - 1]
      for ( uint_t stage = 1; stage < rkStages; stage++ )
      {
         // determine function evaluation points for this stage
         storage.getTimingTree()->start( "Update particle position" );
         for ( auto p : particleStorage )
         {
            auto evaluationPoint = p.getStartPosition();
            for ( uint_t j = 0; j < stage; j++ )
            {
               evaluationPoint += dt * A[stage][j] * p->getK().at( j );
            }
            p->setPosition( evaluationPoint );
         }
         updateParticlePosition( setupStorage, localAndNeighboringPrimitives, particleStorage );
         storage.getTimingTree()->stop( "Update particle position" );

         // sync particles to be able to evaluate the velocity at that point
         storage.getTimingTree()->start( "Sync particles" );
         SNN( particleStorage, setupStorage );
         storage.getTimingTree()->stop( "Sync particles" );

         // evaluate velocity at current particle positions and update k[stage]
         storage.getTimingTree()->start( "Evaluate at particle position" );
         for ( auto p : particleStorage )
         {
            evaluateAtParticlePosition( storage, functions, p, level, results );
            p->getKRef()[stage][0] = - results[0];
            p->getKRef()[stage][1] = - results[1];
            if ( storage.hasGlobalCells() )
               p->getKRef()[stage][2] = - results[2];
         }
         storage.getTimingTree()->stop( "Evaluate at particle position" );
      }

      // all k[i] are now calculated, set final integration result for each particle and
      // assign the position accordingly
      storage.getTimingTree()->start( "Update particle position" );
      for ( auto p : particleStorage )
      {
         auto finalPosition = p->getStartPosition();
         for ( uint_t i = 0; i < rkStages; i++ )
         {
            finalPosition += dt * b[i] * p->getKRef()[i];
         }
         p->setPosition( finalPosition );
         p->setStartPosition( p->getPosition() );
      }
      updateParticlePosition( setupStorage, localAndNeighboringPrimitives, particleStorage );
      storage.getTimingTree()->stop( "Update particle position" );

      // sync particles as position was finally updated
      storage.getTimingTree()->start( "Sync particles" );
      SNN( particleStorage, setupStorage );
      storage.getTimingTree()->stop( "Sync particles" );
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

   std::map< uint_t, int > numParticlesToSendToRank;
   std::map< uint_t, int > numParticlesToReceiveFromRank;

   std::map< uint_t, MPI_Request > sendRequests;

   std::set< walberla::mpi::MPIRank > ranksToReceiveFrom;

   WALBERLA_MPI_SECTION()
   {
      const int TAG = 98234;

      for ( const auto& p : particleStorage )
      {
         if ( numParticlesToSendToRank.count( p.getStartProcess() ) == 0 )
            numParticlesToSendToRank[p.getStartProcess()] = 0;
         numParticlesToSendToRank[p.getStartProcess()]++;
         sendRequests[p.getStartProcess()] = MPI_Request();
      }

      for ( auto& it : numParticlesToSendToRank )
      {
         //         WALBERLA_LOG_INFO( "Particle communcation prep: rank " << walberla::mpi::MPIManager::instance()->rank() << " -> "
         //                                                                << it.first << ": " << it.second )
         MPI_Isend( &it.second,
                    1,
                    MPI_INT,
                    (int) it.first,
                    TAG,
                    walberla::mpi::MPIManager::instance()->comm(),
                    &sendRequests[it.first] );
      }

      for ( auto& it : sendRequests )
      {
         MPI_Status status;
         MPI_Wait( &it.second, &status );
      }

      int numReceivedParticleLocations = 0;
      while ( numReceivedParticleLocations < numberOfCreatedParticles )
      {
         MPI_Status status;

         int numParticlesSum = 0;

         MPI_Recv( &numParticlesSum, 1, MPI_INT, MPI_ANY_SOURCE, TAG, walberla::mpi::MPIManager::instance()->comm(), &status );

         //         WALBERLA_LOG_INFO( "Particle communcation prep: receiving " << numParticlesSum << " particles from rank "
         //                                                                     << status.MPI_SOURCE );
         numReceivedParticleLocations += numParticlesSum;
         numParticlesToReceiveFromRank[uint_c( status.MPI_SOURCE )] = numParticlesSum;
         //         WALBERLA_LOG_INFO( "total received particle infos: " << numReceivedParticleLocations )
         //         WALBERLA_LOG_INFO( "created particles: " << numberOfCreatedParticles )
      }

      // part II
      for ( auto r : numParticlesToReceiveFromRank )
      {
         ranksToReceiveFrom.insert( (walberla::mpi::MPIRank) r.first );
      }
   }

   WALBERLA_NON_MPI_SECTION()
   {
      numParticlesToSendToRank[0]      = (int) particleStorage.size();
      numParticlesToReceiveFromRank[0] = (int) particleStorage.size();
      ranksToReceiveFrom.insert( 0 );
   }

   walberla::mpi::BufferSystem bufferSystem( walberla::mpi::MPIManager::instance()->comm(), 654654555 );

   for ( const auto& p : particleStorage )
   {
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getStartPrimitiveID();
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getStartIndex();
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getStartDoFType();
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getStartEdgeDoFOrientation();
      bufferSystem.sendBuffer( p.getStartProcess() ) << p.getFinalTemperature();
   }

   bufferSystem.setReceiverInfo( ranksToReceiveFrom, true );

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
                  const std::shared_ptr< SetupPrimitiveStorage >& setupStorage,
                  const uint_t                               minLevel,
                  const uint_t                               maxLevel,
                  const TimeSteppingScheme&                  timeSteppingSchemeConvection )
   : storage_( storage )
   , setupStorage_( setupStorage )
   , cOld_( "cOld", storage, minLevel, maxLevel )
   , cTmp_( "cTmp", storage, minLevel, maxLevel )
   , cPlus_( "cPlus", storage, minLevel, maxLevel )
   , cMinus_( "cMinus", storage, minLevel, maxLevel )
   , cAdjusted_( "cAdjusted", storage, minLevel, maxLevel )
   , timeSteppingSchemeConvection_( timeSteppingSchemeConvection )
   , numberOfCreatedParticles_( 0 )
   , particleStorage_( 10000 )
   , localAndNeighboringPrimitives_( getLocalAndNeighboringPrimitives( *setupStorage ) )
   {}


   void step( const FunctionType&                             c,
              const FunctionType&                             ux,
              const FunctionType&                             uy,
              const FunctionType&                             uz,
              const uint_t&                                   level,
              const DoFType&                                  flag,
              const real_t&                                   dt,
              const uint_t&                                   innerSteps,
              const bool&                                     resetParticles = true )
   {
      storage_->getTimingTree()->start( "MMOCTransport" );
      if ( resetParticles )
      {
         cOld_.assign( {1.0}, {c}, level, All );
         storage_->getTimingTree()->start( "Particle initialization" );
         numberOfCreatedParticles_ = initializeParticles(
             particleStorage_, *setupStorage_, *storage_, c, ux, uy, uz, level, Inner, timeSteppingSchemeConvection_, 0 );
         storage_->getTimingTree()->stop( "Particle initialization" );
      }

      storage_->getTimingTree()->start( "Particle integration" );
      particleIntegration(
          particleStorage_, *setupStorage_, *storage_, localAndNeighboringPrimitives_, ux, uy, uz, dt, level, Inner, innerSteps, timeSteppingSchemeConvection_ );
      storage_->getTimingTree()->stop( "Particle integration" );

      storage_->getTimingTree()->start( "Temperature evaluation" );
      evaluateTemperature( particleStorage_, *storage_, c, cOld_, level, Inner, numberOfCreatedParticles_ );
      storage_->getTimingTree()->stop( "Temperature evaluation" );

      storage_->getTimingTree()->stop( "MMOCTransport" );
   }

//   template < typename MassOperator >
//   void step( const FunctionType& c,
//              const FunctionType& ux,
//              const FunctionType& uy,
//              const FunctionType& uz,
//              const uint_t&       level,
//              const DoFType&      flag,
//              const real_t&       dt,
//              const uint_t&       innerSteps,
//              const MassOperator& massOperator,
//              const real_t&       allowedRelativeMassDifference,
//              const real_t&       adjustedAdvectionOffset )
//   {
//      cOld_.assign( {1.0}, {c}, level, All );
//
//      // calculate old mass
//      massOperator.apply( cOld_, cTmp_, level, flag );
//      auto massBefore = cTmp_.sumGlobal( level, flag );
//
//      integrateNodes< FunctionType >(
//          *storage_, c, cOld_, ux, uy, uz, dt, level, flag, innerSteps, timeSteppingSchemeConvection_, 0 );
//
//      // calculate new mass
//      massOperator.apply( c, cTmp_, level, flag );
//      auto massAfter = cTmp_.sumGlobal( level, flag );
//
//      auto relativeMassDifference = std::abs( ( massAfter - massBefore ) / massBefore );
//
//      if ( relativeMassDifference <= allowedRelativeMassDifference )
//         return;
//
//      // perform adjusted advection steps
//      integrateNodes< FunctionType >( *storage_,
//                                      cPlus_,
//                                      cOld_,
//                                      ux,
//                                      uy,
//                                      uz,
//                                      dt,
//                                      level,
//                                      flag,
//                                      innerSteps,
//                                      timeSteppingSchemeConvection_,
//                                      adjustedAdvectionOffset );
//      integrateNodes< FunctionType >( *storage_,
//                                      cMinus_,
//                                      cOld_,
//                                      ux,
//                                      uy,
//                                      uz,
//                                      dt,
//                                      level,
//                                      flag,
//                                      innerSteps,
//                                      timeSteppingSchemeConvection_,
//                                      -adjustedAdvectionOffset );
//
//      // max/min assign functions
//      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > maxAssignment =
//          []( const Point3D&, const std::vector< real_t >& values ) { return std::max( values[0], values[1] ); };
//
//      std::function< real_t( const Point3D&, const std::vector< real_t >& ) > minAssignment =
//          []( const Point3D&, const std::vector< real_t >& values ) { return std::min( values[0], values[1] ); };
//
//      if ( massAfter <= massBefore )
//      {
//         cAdjusted_.interpolate( maxAssignment, {cPlus_, cMinus_}, level );
//      }
//      else
//      {
//         cAdjusted_.interpolate( minAssignment, {cPlus_, cMinus_}, level );
//      }
//
//      // calculate adjustment mass
//      massOperator.apply( cAdjusted_, cTmp_, level, flag );
//      auto massAdjusted = cTmp_.sumGlobal( level, flag );
//
//      auto theta = ( massBefore - massAdjusted ) / ( massAfter - massAdjusted );
//
//      c.assign( {theta, 1 - theta}, {c, cAdjusted_}, level, flag );
//   }

 private:
   const std::shared_ptr< PrimitiveStorage >             storage_;
   const std::shared_ptr< SetupPrimitiveStorage >        setupStorage_;
   FunctionType                                          cOld_;
   FunctionType                                          cTmp_;
   FunctionType                                          cPlus_;
   FunctionType                                          cMinus_;
   FunctionType                                          cAdjusted_;
   TimeSteppingScheme                                    timeSteppingSchemeConvection_;
   uint_t                                                numberOfCreatedParticles_;
   walberla::convection_particles::data::ParticleStorage particleStorage_;
   std::vector< PrimitiveID >                            localAndNeighboringPrimitives_;
};

} // namespace hyteg
