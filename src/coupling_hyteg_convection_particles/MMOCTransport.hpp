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

#include "core/math/MatrixMxN.h"
#include "core/mpi/MPIWrapper.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/functions/FunctionIterator.hpp"
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

#include "convection_particles/data/ParticleStorage.h"
#include "convection_particles/mpi/SyncNextNeighborsNoGhosts.h"
#include "coupling_hyteg_convection_particles/communication/SyncNextNeighborsByPrimitiveID.h"
#include "coupling_hyteg_convection_particles/primitivestorage/PrimitiveStorageConvectionParticlesInterface.hpp"

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
const static std::vector< real_t >                RK_c_ExplicitEuler = {0};

const static std::vector< std::vector< real_t > > RK_A_RK3 = {{{0, 0, 0}, {0.5, 0, 0}, {-1., 2., 0}}};
const static std::vector< real_t >                RK_b_RK3 = {1. / 6., 2. / 3., 1. / 6.};
const static std::vector< real_t >                RK_c_RK3 = {0, 0.5, 1.};

const static std::vector< std::vector< real_t > > RK_A_Ralston = {{{0, 0, 0}, {0.5, 0, 0}, {0, 0.75, 0}}};
const static std::vector< real_t >                RK_b_Ralston = {2. / 9., 1. / 3., 4. / 9.};
const static std::vector< real_t >                RK_c_Ralston = {0, 0.5, 0.75};

const static std::vector< std::vector< real_t > > RK_A_RK4 = {{{0, 0, 0, 0}, {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}}};
const static std::vector< real_t >                RK_b_RK4 = {1. / 6., 1. / 3., 1. / 3., 1. / 6.};
const static std::vector< real_t >                RK_c_RK4 = {0, 0.5, 0.5, 1.};

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

const static std::map< TimeSteppingScheme, std::vector< real_t > > RK_c = {
    {TimeSteppingScheme::ExplicitEuler, RK_c_ExplicitEuler},
    {TimeSteppingScheme::RK3, RK_c_RK3},
    {TimeSteppingScheme::Ralston, RK_c_Ralston},
    {TimeSteppingScheme::RK4, RK_c_RK4}};

inline std::vector< PrimitiveID > getNeighboringPrimitives( const PrimitiveID& primitiveID, const PrimitiveStorage& storage )
{
   std::set< PrimitiveID > neighboringPrimitives;

   WALBERLA_CHECK( storage.primitiveExistsLocally( primitiveID ) );

   if ( storage.hasGlobalCells() )
   {
      auto cellID = primitiveID;
      auto cell   = storage.getCell( cellID );

      for ( const auto& vertexID : cell->neighborVertices() )
      {
         auto vertex = storage.getVertex( vertexID );
         for ( const auto& neighborCellID : vertex->neighborCells() )
         {
            neighboringPrimitives.insert( neighborCellID );
         }
      }
   }
   else
   {
      auto faceID = primitiveID;
      auto face   = storage.getFace( faceID );

      for ( const auto& vertexID : face->neighborVertices() )
      {
         auto vertex = storage.getVertex( vertexID );
         for ( const auto& neighborFaceID : vertex->neighborFaces() )
         {
            neighboringPrimitives.insert( neighborFaceID );
         }
      }
   }
   return std::vector< PrimitiveID >( neighboringPrimitives.begin(), neighboringPrimitives.end() );
}

inline void updateParticlePosition( const PrimitiveStorage&                                storage,
                                    walberla::convection_particles::data::ParticleStorage& particleStorage,
                                    const real_t&                                          particleLocationRadius )
{
   if ( !storage.hasGlobalCells() )
   {
      for ( auto p : particleStorage )
      {
         p.setOutsideDomain( 0 );
         bool foundByPointLocation = false;

         // check for current cell (probability is high that we find the particle here...)
         const auto faceID = p->getContainingPrimitive();
         WALBERLA_ASSERT( storage.faceExistsLocally( faceID ) || storage.faceExistsInNeighborhood( faceID ) );
         const auto face = storage.getFace( faceID );

         Point3D computationalLocation;
         face->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocation );
         Point2D computationalLocation2D( {computationalLocation[0], computationalLocation[1]} );

         if ( isPointInTriangle( computationalLocation2D,
                                 Point2D( {face->getCoordinates().at( 0 )[0], face->getCoordinates().at( 0 )[1]} ),
                                 Point2D( {face->getCoordinates().at( 1 )[0], face->getCoordinates().at( 1 )[1]} ),
                                 Point2D( {face->getCoordinates().at( 2 )[0], face->getCoordinates().at( 2 )[1]} ) ) )
         {
            p->setContainingPrimitive( faceID );
            foundByPointLocation = true;
         }
         else
         {
            // check for neighbor cells if we did not find it in its previous cell
            const auto neighboringFaces = getNeighboringPrimitives( faceID, storage );
            for ( const auto& neighborFaceID : neighboringFaces )
            {
               WALBERLA_ASSERT( storage.faceExistsLocally( neighborFaceID ) ||
                                storage.faceExistsInNeighborhood( neighborFaceID ) );
               const auto neighborFace = storage.getFace( neighborFaceID );
               Point3D    computationalLocationNeighbor;
               neighborFace->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocationNeighbor );
               Point2D computationalLocationNeighbor2D( {computationalLocationNeighbor[0], computationalLocationNeighbor[1]} );

               if ( isPointInTriangle(
                        computationalLocationNeighbor2D,
                        Point2D( {neighborFace->getCoordinates().at( 0 )[0], neighborFace->getCoordinates().at( 0 )[1]} ),
                        Point2D( {neighborFace->getCoordinates().at( 1 )[0], neighborFace->getCoordinates().at( 1 )[1]} ),
                        Point2D( {neighborFace->getCoordinates().at( 2 )[0], neighborFace->getCoordinates().at( 2 )[1]} ) ) )
               {
                  // set it to the first neighbor we found to contain the particle
                  p->setContainingPrimitive( neighborFaceID );
                  foundByPointLocation = true;
                  break;
               }
            }
         }

         if ( !foundByPointLocation )
         {
            // At this point there are still three possible scenarios regarding the location of the particle:
            // 1. The particle is outside the neighborhood -> timestep too large, we do not care and crash.
            // 2. The particle is outside of the entire domain -> we set the outsideDomain flag.
            // 3. The particle is in the neighborhood patch, but floating-point errors made all point location
            //    calculations return false. We therefore check with a larger radius.
            if ( sphereTriangleIntersection( computationalLocation,
                                             particleLocationRadius,
                                             face->getCoordinates().at( 0 ),
                                             face->getCoordinates().at( 1 ),
                                             face->getCoordinates().at( 2 ) ) )
            {
               p->setContainingPrimitive( faceID );
               foundByPointLocation = true;
            }
            else
            {
               const auto neighboringFaces = getNeighboringPrimitives( faceID, storage );
               for ( const auto& neighborFaceID : neighboringFaces )
               {
                  WALBERLA_ASSERT( storage.faceExistsLocally( neighborFaceID ) ||
                                   storage.faceExistsInNeighborhood( neighborFaceID ) );
                  const auto neighborFace = storage.getFace( neighborFaceID );
                  Point3D    computationalLocationNeighbor;
                  neighborFace->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocationNeighbor );

                  if ( sphereTriangleIntersection( computationalLocationNeighbor,
                                                   particleLocationRadius,
                                                   neighborFace->getCoordinates().at( 0 ),
                                                   neighborFace->getCoordinates().at( 1 ),
                                                   neighborFace->getCoordinates().at( 2 ) ) )
                  {
                     p->setContainingPrimitive( neighborFaceID );
                     foundByPointLocation = true;
                     break;
                  }
               }
            }
         }

         if ( !foundByPointLocation )
         {
            p->setOutsideDomain( 1 );
         }
      }
   }
   else
   {
      for ( auto p : particleStorage )
      {
         p.setOutsideDomain( 0 );
         bool foundByPointLocation = false;

         // check for current cell (probability is high that we find the particle here...)
         const auto cellID = p->getContainingPrimitive();
         WALBERLA_ASSERT( storage.cellExistsLocally( cellID ) || storage.cellExistsInNeighborhood( cellID ) );
         const auto cell = storage.getCell( cellID );

         Point3D computationalLocation;
         cell->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocation );

         if ( isPointInTetrahedron( computationalLocation,
                                    cell->getCoordinates().at( 0 ),
                                    cell->getCoordinates().at( 1 ),
                                    cell->getCoordinates().at( 2 ),
                                    cell->getCoordinates().at( 3 ),
                                    cell->getFaceInwardNormal( 0 ),
                                    cell->getFaceInwardNormal( 1 ),
                                    cell->getFaceInwardNormal( 2 ),
                                    cell->getFaceInwardNormal( 3 ) ) )
         {
            p->setContainingPrimitive( cellID );
            foundByPointLocation = true;
         }
         else
         {
            // check for neighbor cells if we did not find it in its previous cell
            const auto& neighboringCells = cell->getIndirectNeighborCellIDs();
            for ( const auto& neighborCellID : neighboringCells )
            {
               WALBERLA_ASSERT( storage.cellExistsLocally( neighborCellID ) ||
                                storage.cellExistsInNeighborhood( neighborCellID ) );
               const auto neighborCell = storage.getCell( neighborCellID );
               Point3D    computationalLocationNeighbor;
               neighborCell->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocationNeighbor );

               if ( isPointInTetrahedron( computationalLocationNeighbor,
                                          neighborCell->getCoordinates().at( 0 ),
                                          neighborCell->getCoordinates().at( 1 ),
                                          neighborCell->getCoordinates().at( 2 ),
                                          neighborCell->getCoordinates().at( 3 ),
                                          neighborCell->getFaceInwardNormal( 0 ),
                                          neighborCell->getFaceInwardNormal( 1 ),
                                          neighborCell->getFaceInwardNormal( 2 ),
                                          neighborCell->getFaceInwardNormal( 3 ) ) )
               {
                  // set it to the first neighbor we found to contain the particle
                  p->setContainingPrimitive( neighborCellID );
                  foundByPointLocation = true;
                  break;
               }
            }
         }

         if ( !foundByPointLocation )
         {
            // At this point there are still three possible scenarios regarding the location of the particle:
            // 1. The particle is outside the neighborhood -> timestep too large, we do not care and crash.
            // 2. The particle is outside of the entire domain -> we set the outsideDomain flag.
            // 3. The particle is in the neighborhood patch, but floating-point errors made all point location
            //    calculations return false. We therefore check with a larger radius.
            if ( sphereTetrahedronIntersection( computationalLocation,
                                                particleLocationRadius,
                                                cell->getCoordinates().at( 0 ),
                                                cell->getCoordinates().at( 1 ),
                                                cell->getCoordinates().at( 2 ),
                                                cell->getCoordinates().at( 3 ) ) )
            {
               p->setContainingPrimitive( cellID );
               foundByPointLocation = true;
            }
            else
            {
               const auto& neighboringCells = cell->getIndirectNeighborCellIDs();
               for ( const auto& neighborCellID : neighboringCells )
               {
                  WALBERLA_ASSERT( storage.cellExistsLocally( neighborCellID ) ||
                                   storage.cellExistsInNeighborhood( neighborCellID ) );
                  const auto neighborCell = storage.getCell( neighborCellID );
                  Point3D    computationalLocationNeighbor;
                  neighborCell->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocationNeighbor );

                  if ( sphereTetrahedronIntersection( computationalLocationNeighbor,
                                                      particleLocationRadius,
                                                      neighborCell->getCoordinates().at( 0 ),
                                                      neighborCell->getCoordinates().at( 1 ),
                                                      neighborCell->getCoordinates().at( 2 ),
                                                      neighborCell->getCoordinates().at( 3 ) ) )
                  {
                     p->setContainingPrimitive( neighborCellID );
                     foundByPointLocation = true;
                     break;
                  }
               }
            }
         }

         if ( !foundByPointLocation )
         {
            p->setOutsideDomain( 1 );
         }
      }
   }
}

template < typename FunctionType >
inline real_t evaluateAtParticlePosition( PrimitiveStorage&                                                      storage,
                                          const FunctionType&                                                    function,
                                          const walberla::convection_particles::data::ParticleStorage::Particle& particle,
                                          const uint_t&                                                          level,
                                          const bool& setParticlesOutsideDomainToZero )
{
   if ( setParticlesOutsideDomainToZero && particle.getOutsideDomain() == 1 )
   {
      return real_c( 0 );
   }

   real_t result;
   if ( !storage.hasGlobalCells() )
   {
      WALBERLA_CHECK( storage.faceExistsLocally( particle.getContainingPrimitive() ) );
      Face&   face = *storage.getFace( particle.getContainingPrimitive() );
      Point3D computationalLocation;
      face.getGeometryMap()->evalFinv( toPoint3D( particle.getPosition() ), computationalLocation );

      if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
      {
         result = vertexdof::macroface::evaluate( level, face, computationalLocation, function.getFaceDataID() );
      }
      else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
      {
         result = P2::macroface::evaluate( level,
                                           face,
                                           computationalLocation,
                                           function.getVertexDoFFunction().getFaceDataID(),
                                           function.getEdgeDoFFunction().getFaceDataID() );
      }
      else
      {
         WALBERLA_ABORT( "Not implemented for this discretization." )
      }
   }
   else
   {
      WALBERLA_CHECK( storage.cellExistsLocally( particle.getContainingPrimitive() ) );
      Cell&   cell = *storage.getCell( particle.getContainingPrimitive() );
      Point3D computationalLocation;
      cell.getGeometryMap()->evalFinv( toPoint3D( particle.getPosition() ), computationalLocation );

      if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
      {
         result = vertexdof::macrocell::evaluate( level, cell, computationalLocation, function.getCellDataID() );
      }
      else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
      {
         result = P2::macrocell::evaluate( level,
                                           cell,
                                           computationalLocation,
                                           function.getVertexDoFFunction().getCellDataID(),
                                           function.getEdgeDoFFunction().getCellDataID() );
      }
      else
      {
         WALBERLA_ABORT( "Not implemented for this discretization." )
      }
   }
   return result;
}

template < typename FunctionType >
inline void evaluateAtParticlePosition( PrimitiveStorage&                                                      storage,
                                        const std::vector< FunctionType >&                                     functions,
                                        const walberla::convection_particles::data::ParticleStorage::Particle& particle,
                                        const uint_t&                                                          level,
                                        std::vector< real_t >&                                                 results,
                                        const bool& setParticlesOutsideDomainToZero )
{
   if ( setParticlesOutsideDomainToZero && particle.getOutsideDomain() == 1 )
   {
      for ( uint_t i = 0; i < functions.size(); i++ )
      {
         results[i] = real_c( 0 );
      }
      return;
   }

   if ( !storage.hasGlobalCells() )
   {
      WALBERLA_CHECK( storage.faceExistsLocally( particle.getContainingPrimitive() ) );
      Face& face = *storage.getFace( particle.getContainingPrimitive() );

      Point3D computationalLocation;
      face.getGeometryMap()->evalFinv( toPoint3D( particle.getPosition() ), computationalLocation );

      for ( uint_t i = 0; i < functions.size(); i++ )
      {
         if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
         {
            results[i] = vertexdof::macroface::evaluate( level, face, computationalLocation, functions[i].getFaceDataID() );
         }
         else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
         {
            results[i] = P2::macroface::evaluate( level,
                                                  face,
                                                  computationalLocation,
                                                  functions[i].getVertexDoFFunction().getFaceDataID(),
                                                  functions[i].getEdgeDoFFunction().getFaceDataID() );
         }
         else
         {
            WALBERLA_ABORT( "Not implemented for this discretization." )
         }
      }
   }
   else
   {
      std::vector< PrimitiveDataID< FunctionMemory< real_t >, Cell > > vertexDataIDs;
      std::vector< PrimitiveDataID< FunctionMemory< real_t >, Cell > > edgeDataIDs;

      for ( uint_t i = 0; i < functions.size(); i++ )
      {
         if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
         {
            vertexDataIDs.push_back( functions[i].getCellDataID() );
         }
         else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
         {
            vertexDataIDs.push_back( functions[i].getVertexDoFFunction().getCellDataID() );
            edgeDataIDs.push_back( functions[i].getEdgeDoFFunction().getCellDataID() );
         }
         else
         {
            WALBERLA_ABORT( "Not implemented for this discretization." )
         }
      }

      WALBERLA_CHECK( storage.cellExistsLocally( particle.getContainingPrimitive() ) );
      Cell&   cell = *storage.getCell( particle.getContainingPrimitive() );
      Point3D computationalLocation;
      cell.getGeometryMap()->evalFinv( toPoint3D( particle.getPosition() ), computationalLocation );
      if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
      {
         for ( uint_t i = 0; i < functions.size(); i++ )
         {
            results[i] = vertexdof::macrocell::evaluate( level, cell, computationalLocation, functions[i].getCellDataID() );
         }
      }
      else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
      {
         P2::macrocell::evaluate( level, cell, computationalLocation, vertexDataIDs, edgeDataIDs, results );
      }
      else
      {
         WALBERLA_ABORT( "Not implemented for this discretization." )
      }
   }
}

template < typename FunctionType >
inline uint_t initializeParticles( walberla::convection_particles::data::ParticleStorage& particleStorage,
                                   PrimitiveStorage&                                      storage,
                                   const FunctionType&                                    c,
                                   const FunctionType&                                    ux,
                                   const FunctionType&                                    uy,
                                   const FunctionType&                                    uz,
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
   // const std::vector< std::vector< real_t > >& A        = RK_A.at( timeSteppingScheme ); //this seems to be unused?
   const std::vector< real_t >&                b        = RK_b.at( timeSteppingScheme );
   const uint_t                                rkStages = b.size();

   if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
   {
      for ( auto it : FunctionIterator< vertexdof::VertexDoFFunction< real_t > >( c, level ) )
      {
         // if ( storage.onBoundary( it.primitiveID(), true ) )
         //   continue;

         Point3D physicalLocation;
         auto    primitive = storage.getPrimitive( it.primitiveID() );
         primitive->getGeometryMap()->evalF( it.coordinates(), physicalLocation );

         auto particleIt = particleStorage.create();
         particleIt->setOwner( (int) rank );
         particleIt->setPosition( toVec3( physicalLocation ) );
         particleIt->setStartPosition( toVec3( physicalLocation ) );
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
            particleIt->getKRef()[0][0] =
                storage.getVertex( it.primitiveID() )->getData( ux.getVertexDataID() )->getPointer( level )[it.arrayIndex()];
            particleIt->getKRef()[0][1] =
                storage.getVertex( it.primitiveID() )->getData( uy.getVertexDataID() )->getPointer( level )[it.arrayIndex()];
            if ( storage.hasGlobalCells() )
               particleIt->getKRef()[0][2] =
                   storage.getVertex( it.primitiveID() )->getData( uz.getVertexDataID() )->getPointer( level )[it.arrayIndex()];
         }
         if ( it.isOnMacroEdge() )
         {
            particleIt->getKRef()[0][0] =
                storage.getEdge( it.primitiveID() )->getData( ux.getEdgeDataID() )->getPointer( level )[it.arrayIndex()];
            particleIt->getKRef()[0][1] =
                storage.getEdge( it.primitiveID() )->getData( uy.getEdgeDataID() )->getPointer( level )[it.arrayIndex()];
            if ( storage.hasGlobalCells() )
               particleIt->getKRef()[0][2] =
                   storage.getEdge( it.primitiveID() )->getData( uz.getEdgeDataID() )->getPointer( level )[it.arrayIndex()];
         }
         if ( it.isOnMacroFace() )
         {
            particleIt->getKRef()[0][0] =
                storage.getFace( it.primitiveID() )->getData( ux.getFaceDataID() )->getPointer( level )[it.arrayIndex()];
            particleIt->getKRef()[0][1] =
                storage.getFace( it.primitiveID() )->getData( uy.getFaceDataID() )->getPointer( level )[it.arrayIndex()];
            if ( storage.hasGlobalCells() )
               particleIt->getKRef()[0][2] =
                   storage.getFace( it.primitiveID() )->getData( uz.getFaceDataID() )->getPointer( level )[it.arrayIndex()];
         }
      }
   }
   else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
   {
      for ( auto it : FunctionIterator< vertexdof::VertexDoFFunction< real_t > >( c.getVertexDoFFunction(), level ) )
      {
         // if ( storage.onBoundary( it.primitiveID(), true ) )
         //   continue;

         Point3D physicalLocation;
         auto    primitive = storage.getPrimitive( it.primitiveID() );
         primitive->getGeometryMap()->evalF( it.coordinates(), physicalLocation );

         auto particleIt = particleStorage.create();
         particleIt->setOwner( (int) rank );
         particleIt->setPosition( toVec3( physicalLocation ) );
         particleIt->setStartPosition( toVec3( physicalLocation ) );
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
         //      if ( storage.onBoundary( it.primitiveID(), true ) )
         //         continue;

         Point3D physicalLocation;
         auto    primitive = storage.getPrimitive( it.primitiveID() );
         primitive->getGeometryMap()->evalF( it.coordinates(), physicalLocation );

         auto particleIt = particleStorage.create();
         particleIt->setOwner( (int) rank );
         particleIt->setPosition( toVec3( physicalLocation ) );
         particleIt->setStartPosition( toVec3( physicalLocation ) );
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
   }
   else
   {
      WALBERLA_ABORT( "Not implemented for this discretization." )
   }

   const uint_t numberOfCreatedParticles = particleStorage.size();
   // WALBERLA_LOG_INFO( "Particles after creation: " << particleStorage.size() );

   // now sync particles so that all particles are "assigned" a macro-cell (macro-face in 2D) automatically
   // walberla::convection_particles::mpi::SyncNextNeighborsNoGhosts SNN;
   walberla::convection_particles::mpi::SyncNextNeighborsByPrimitiveID SNN;

   SNN( particleStorage, storage );

   // WALBERLA_LOG_INFO( "Particles after init sync: " << particleStorage.size() );
   return numberOfCreatedParticles;
}

template < typename FunctionType >
inline void particleIntegration( walberla::convection_particles::data::ParticleStorage& particleStorage,
                                 PrimitiveStorage&                                      storage,
                                 const FunctionType&                                    ux,
                                 const FunctionType&                                    uy,
                                 const FunctionType&                                    uz,
                                 const FunctionType&                                    uxLastTimeStep,
                                 const FunctionType&                                    uyLastTimeStep,
                                 const FunctionType&                                    uzLastTimeStep,
                                 const real_t&                                          dt,
                                 const uint_t&                                          level,
                                 const DoFType&,
                                 const uint_t&             steps,
                                 const TimeSteppingScheme& timeSteppingScheme,
                                 const real_t&             particleLocationRadius,
                                 const bool&               setParticlesOutsideDomainToZero )
{
   communication::syncFunctionBetweenPrimitives( ux, level );
   communication::syncFunctionBetweenPrimitives( uy, level );
   communication::syncFunctionBetweenPrimitives( uz, level );
   communication::syncFunctionBetweenPrimitives( uxLastTimeStep, level );
   communication::syncFunctionBetweenPrimitives( uyLastTimeStep, level );
   communication::syncFunctionBetweenPrimitives( uzLastTimeStep, level );

   walberla::convection_particles::mpi::SyncNextNeighborsByPrimitiveID SNN;

   const uint_t                                rank     = uint_c( walberla::mpi::MPIManager::instance()->rank() );
   const std::vector< std::vector< real_t > >& A        = RK_A.at( timeSteppingScheme );
   const std::vector< real_t >&                b        = RK_b.at( timeSteppingScheme );
   const std::vector< real_t >&                c        = RK_c.at( timeSteppingScheme );
   const uint_t                                rkStages = b.size();

   storage.getTimingTree()->start( "Sync particles" );
   SNN( particleStorage, storage, true );
   storage.getTimingTree()->stop( "Sync particles" );

   for ( uint_t step = 0; step < steps; step++ )
   {
      // WALBERLA_LOG_INFO_ON_ROOT( "Starting inner time step " << step << " ..." )
      // TODO: sort particle storage by governing macro (optimization, tbd)

      // RK stage 0
      // skip setting to start pos (already happened)
      // skip synchronization (already happened)
      // assuming (true for explicit RK schemes) that c[0] == 0

      storage.getTimingTree()->start( "Evaluate at particle position" );

      std::vector< real_t >       results( {0, 0} );
      std::vector< real_t >       resultsLastTimeStep( {0, 0} );
      std::vector< FunctionType > functions             = {ux, uy};
      std::vector< FunctionType > functionsLastTimeStep = {uxLastTimeStep, uyLastTimeStep};
      if ( storage.hasGlobalCells() )
      {
         results.push_back( 0 );
         resultsLastTimeStep.push_back( 0 );
         functions.push_back( uz );
         functionsLastTimeStep.push_back( uzLastTimeStep );
      }

      for ( auto p : particleStorage )
      {
         evaluateAtParticlePosition( storage, functions, p, level, results, setParticlesOutsideDomainToZero );
         p->getKRef()[0][0] = -results[0];
         p->getKRef()[0][1] = -results[1];
         if ( storage.hasGlobalCells() )
            p->getKRef()[0][2] = -results[2];
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
         updateParticlePosition( storage, particleStorage, particleLocationRadius );
         storage.getTimingTree()->stop( "Update particle position" );

         // sync particles to be able to evaluate the velocity at that point
         storage.getTimingTree()->start( "Sync particles" );
         SNN( particleStorage, storage, true );
         storage.getTimingTree()->stop( "Sync particles" );

         // evaluate velocity at current particle positions and update k[stage]
         // we perform a linear interpolation here using the c-weights of the RK method
         // and the "current" and "last" velocity functions
         storage.getTimingTree()->start( "Evaluate at particle position" );
         for ( auto p : particleStorage )
         {
            evaluateAtParticlePosition( storage, functions, p, level, results, setParticlesOutsideDomainToZero );
            evaluateAtParticlePosition(
                storage, functionsLastTimeStep, p, level, resultsLastTimeStep, setParticlesOutsideDomainToZero );
            p->getKRef()[stage][0] = -( ( 1.0 - c[stage] ) * results[0] + c[stage] * resultsLastTimeStep[0] );
            p->getKRef()[stage][1] = -( ( 1.0 - c[stage] ) * results[1] + c[stage] * resultsLastTimeStep[1] );
            if ( storage.hasGlobalCells() )
               p->getKRef()[stage][2] = -( ( 1.0 - c[stage] ) * results[2] + c[stage] * resultsLastTimeStep[2] );
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
      updateParticlePosition( storage, particleStorage, particleLocationRadius );
      storage.getTimingTree()->stop( "Update particle position" );

      // sync particles as position was finally updated
      storage.getTimingTree()->start( "Sync particles" );
      SNN( particleStorage, storage, true );
      storage.getTimingTree()->stop( "Sync particles" );
   }
}

template < typename FunctionType >
inline void evaluateTemperature( walberla::convection_particles::data::ParticleStorage& particleStorage,
                                 PrimitiveStorage&                                      storage,
                                 const FunctionType&                                    c,
                                 const FunctionType&                                    cOld,
                                 const uint_t&                                          level,
                                 const DoFType&,
                                 const uint_t& numberOfCreatedParticles,
                                 const bool    globalMaxLimiter                = true,
                                 const bool    setParticlesOutsideDomainToZero = false )
{
   communication::syncFunctionBetweenPrimitives( cOld, level );

   // limiting temperature evaluation
   real_t minTempCOld = 0;
   real_t maxTempCOld = 0;

   if ( globalMaxLimiter )
   {
      minTempCOld = cOld.getMinValue( level );
      maxTempCOld = cOld.getMaxValue( level );
   }

   // evaluate temperature at final position
   for ( auto p : particleStorage )
   {
      auto finalTemperature = evaluateAtParticlePosition( storage, cOld, p, level, setParticlesOutsideDomainToZero );
      if ( globalMaxLimiter )
      {
         finalTemperature = std::max( finalTemperature, minTempCOld );
         finalTemperature = std::min( finalTemperature, maxTempCOld );
      }
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
#ifdef _MSC_VER
      //need a first receive to avoid blocking on windows while communicating with itself
      int selfCommMessage = 0;
      MPI_Request selfCommRequest;
      MPI_Irecv(&selfCommMessage, 1, MPI_INT, walberla::mpi::MPIManager::instance()->rank(), TAG, walberla::mpi::MPIManager::instance()->comm(), &selfCommRequest);
#endif
      for ( const auto& p : particleStorage )
      {
         if ( numParticlesToSendToRank.count( p.getStartProcess() ) == 0 ){
            numParticlesToSendToRank[p.getStartProcess()] = 0;
            sendRequests[p.getStartProcess()] = MPI_Request();
         }
         numParticlesToSendToRank[p.getStartProcess()]++;
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
      #ifdef _MSC_VER
      //get self communication
      {
         int ready;
         MPI_Status status;
         MPI_Test(&selfCommRequest, &ready, &status);
         if(ready){
            numReceivedParticleLocations = selfCommMessage;
            numParticlesToReceiveFromRank[uint_c( status.MPI_SOURCE )] = numReceivedParticleLocations;
         }else{
            WALBERLA_LOG_INFO("somethings very wrong here");
         }

      }
      #endif
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
            if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
            {
               vertex->getData( c.getVertexDataID() )->getPointer( level )[0] = temp;
            }
            else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
            {
               vertex->getData( c.getVertexDoFFunction().getVertexDataID() )->getPointer( level )[0] = temp;
            }
            else
            {
               WALBERLA_ABORT( "Not implemented for this discretization." );
            }
         }
         else if ( storage.edgeExistsLocally( primitiveID ) )
         {
            auto edge = storage.getEdge( primitiveID );
            if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
            {
               WALBERLA_CHECK_EQUAL( dofType, 0 );
               edge->getData( c.getEdgeDataID() )->getPointer( level )[vertexdof::macroedge::index( level, index.x() )] = temp;
            }
            else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
            {
               if ( dofType == 0 )
               {
                  edge->getData( c.getVertexDoFFunction().getEdgeDataID() )
                      ->getPointer( level )[vertexdof::macroedge::index( level, index.x() )] = temp;
               }
               else
               {
                  edge->getData( c.getEdgeDoFFunction().getEdgeDataID() )
                      ->getPointer( level )[edgedof::macroedge::index( level, index.x() )] = temp;
               }
            }
            else
            {
               WALBERLA_ABORT( "Not implemented for this discretization." );
            }
         }
         else if ( storage.faceExistsLocally( primitiveID ) )
         {
            auto face = storage.getFace( primitiveID );
            if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
            {
               WALBERLA_CHECK_EQUAL( dofType, 0 );
               face->getData( c.getFaceDataID() )
                   ->getPointer( level )[vertexdof::macroface::index( level, index.x(), index.y() )] = temp;
            }
            else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
            {
               if ( dofType == 0 )
               {
                  face->getData( c.getVertexDoFFunction().getFaceDataID() )
                      ->getPointer( level )[vertexdof::macroface::index( level, index.x(), index.y() )] = temp;
               }
               else
               {
                  face->getData( c.getEdgeDoFFunction().getFaceDataID() )
                      ->getPointer( level )[edgedof::macroface::index( level, index.x(), index.y(), orientation )] = temp;
               }
            }
            else
            {
               WALBERLA_ABORT( "Not implemented for this discretization." );
            }
         }
         else if ( storage.cellExistsLocally( primitiveID ) )
         {
            auto cell = storage.getCell( primitiveID );
            if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
            {
               WALBERLA_CHECK_EQUAL( dofType, 0 );
               cell->getData( c.getCellDataID() )
                   ->getPointer( level )[vertexdof::macrocell::index( level, index.x(), index.y(), index.z() )] = temp;
            }
            else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
            {
               if ( dofType == 0 )
               {
                  cell->getData( c.getVertexDoFFunction().getCellDataID() )
                      ->getPointer( level )[vertexdof::macrocell::index( level, index.x(), index.y(), index.z() )] = temp;
               }
               else
               {
                  cell->getData( c.getEdgeDoFFunction().getCellDataID() )
                      ->getPointer( level )[edgedof::macrocell::index( level, index.x(), index.y(), index.z(), orientation )] =
                      temp;
               }
            }
            else
            {
               WALBERLA_ABORT( "Not implemented for this discretization." );
            }
         }
      }
   }
}

template < typename FunctionType >
class MMOCTransport
{
 public:
   typedef typename FunctionTrait< FunctionType >::AssocVectorFunctionType vecfun_t;

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
   {
      WALBERLA_CHECK_GREATER_EQUAL( storage->getAdditionalHaloDepth(),
                                    1,
                                    "For the particle transport implementation, "
                                    "the additional halo depth of the PrimitiveStorage "
                                    "must at least be set to 1." )
      particleLocationRadius_ = 0.1 * MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   }

   void step( const FunctionType&                                                    c,
              const vecfun_t& u,
              const vecfun_t& uLastTimeStep,
              const uint_t&                                                          level,
              const DoFType&                                                         flag,
              const real_t&                                                          dt,
              const uint_t&                                                          innerSteps,
              const bool&                                                            resetParticles                  = true,
              const bool&                                                            globalMaxLimiter                = true,
              const bool&                                                            setParticlesOutsideDomainToZero = false )
   {
      uint_t aux = u.getDimension() == 3 ? 2 : 0;
      step( c,
            u[0],
            u[1],
            u[aux],
            uLastTimeStep[0],
            uLastTimeStep[1],
            uLastTimeStep[aux],
            level,
            flag,
            dt,
            innerSteps,
            resetParticles,
            globalMaxLimiter,
            setParticlesOutsideDomainToZero );
   }

   template < typename MassOperator >
   void step( const FunctionType&       c,
              const vecfun_t& u,
              const vecfun_t& uLastTimeStep,
              const uint_t&             level,
              const DoFType&            flag,
              const real_t&             dt,
              const uint_t&             innerSteps,
              const MassOperator&       massOperator,
              const real_t&             allowedRelativeMassDifference,
              const real_t&             dtPertubationAdjustedAdvection,
              const bool&               globalMaxLimiter                = true,
              const bool&               setParticlesOutsideDomainToZero = false )
   {
      uint_t aux = u.getDimension() == 3 ? 2 : 0;
      step( c,
            u[0],
            u[1],
            u[aux],
            uLastTimeStep[0],
            uLastTimeStep[1],
            uLastTimeStep[aux],
            level,
            flag,
            dt,
            innerSteps,
            massOperator,
            allowedRelativeMassDifference,
            dtPertubationAdjustedAdvection,
            globalMaxLimiter,
            setParticlesOutsideDomainToZero );
   }

 private:
   void step( const FunctionType& c,
              const FunctionType& ux,
              const FunctionType& uy,
              const FunctionType& uz,
              const FunctionType& uxLastTimeStep,
              const FunctionType& uyLastTimeStep,
              const FunctionType& uzLastTimeStep,
              const uint_t&       level,
              const DoFType&      flag,
              const real_t&       dt,
              const uint_t&       innerSteps,
              const bool&         resetParticles                  = true,
              const bool&         globalMaxLimiter                = true,
              const bool&         setParticlesOutsideDomainToZero = false )
   {
      storage_->getTimingTree()->start( "MMOCTransport" );

      cOld_.copyBoundaryConditionFromFunction( c );
      cTmp_.copyBoundaryConditionFromFunction( c );
      cPlus_.copyBoundaryConditionFromFunction( c );
      cMinus_.copyBoundaryConditionFromFunction( c );
      cAdjusted_.copyBoundaryConditionFromFunction( c );

      if ( resetParticles )
      {
         cOld_.assign( {1.0}, {c}, level, All );
         storage_->getTimingTree()->start( "Particle initialization" );
         numberOfCreatedParticles_ =
             initializeParticles( particleStorage_, *storage_, c, ux, uy, uz, level, Inner, timeSteppingSchemeConvection_, 0 );
         storage_->getTimingTree()->stop( "Particle initialization" );
      }

      storage_->getTimingTree()->start( "Particle integration" );
      particleIntegration( particleStorage_,
                           *storage_,
                           ux,
                           uy,
                           uz,
                           uxLastTimeStep,
                           uyLastTimeStep,
                           uzLastTimeStep,
                           dt,
                           level,
                           Inner,
                           innerSteps,
                           timeSteppingSchemeConvection_,
                           particleLocationRadius_,
                           setParticlesOutsideDomainToZero );
      storage_->getTimingTree()->stop( "Particle integration" );

      storage_->getTimingTree()->start( "Temperature evaluation" );
      evaluateTemperature( particleStorage_,
                           *storage_,
                           c,
                           cOld_,
                           level,
                           Inner,
                           numberOfCreatedParticles_,
                           globalMaxLimiter,
                           setParticlesOutsideDomainToZero );
      storage_->getTimingTree()->stop( "Temperature evaluation" );

      storage_->getTimingTree()->stop( "MMOCTransport" );
   }

   template < typename MassOperator >
   void step( const FunctionType& c,
              const FunctionType& ux,
              const FunctionType& uy,
              const FunctionType& uz,
              const FunctionType& uxLastTimeStep,
              const FunctionType& uyLastTimeStep,
              const FunctionType& uzLastTimeStep,
              const uint_t&       level,
              const DoFType&      flag,
              const real_t&       dt,
              const uint_t&       innerSteps,
              const MassOperator& massOperator,
              const real_t&       allowedRelativeMassDifference,
              const real_t&       dtPertubationAdjustedAdvection,
              const bool&         globalMaxLimiter                = true,
              const bool&         setParticlesOutsideDomainToZero = false )
   {
      cOld_.copyBoundaryConditionFromFunction( c );
      cTmp_.copyBoundaryConditionFromFunction( c );
      cPlus_.copyBoundaryConditionFromFunction( c );
      cMinus_.copyBoundaryConditionFromFunction( c );
      cAdjusted_.copyBoundaryConditionFromFunction( c );

      cPlus_.assign( {1.0}, {c}, level, All );
      cMinus_.assign( {1.0}, {c}, level, All );

      // calculate old mass
      massOperator.apply( c, cTmp_, level, flag );
      auto massBefore = cTmp_.sumGlobal( level, flag );

      step( c,
            ux,
            uy,
            uz,
            uxLastTimeStep,
            uyLastTimeStep,
            uzLastTimeStep,
            level,
            flag,
            dt,
            innerSteps,
            true,
            globalMaxLimiter,
            setParticlesOutsideDomainToZero );

      // calculate new mass
      massOperator.apply( c, cTmp_, level, flag );
      auto massAfter = cTmp_.sumGlobal( level, flag );

      auto relativeMassDifference = std::abs( ( massAfter - massBefore ) / massBefore );

      if ( relativeMassDifference <= allowedRelativeMassDifference )
         return;

      // perform adjusted advection
      storage_->getTimingTree()->start( "MMOC with adjusted advection" );
      step( cPlus_,
            ux,
            uy,
            uz,
            uxLastTimeStep,
            uyLastTimeStep,
            uzLastTimeStep,
            level,
            flag,
            dt + dtPertubationAdjustedAdvection,
            innerSteps,
            true,
            globalMaxLimiter,
            setParticlesOutsideDomainToZero );
      step( cMinus_,
            ux,
            uy,
            uz,
            uxLastTimeStep,
            uyLastTimeStep,
            uzLastTimeStep,
            level,
            flag,
            dt - dtPertubationAdjustedAdvection,
            innerSteps,
            true,
            globalMaxLimiter,
            setParticlesOutsideDomainToZero );
      storage_->getTimingTree()->stop( "MMOC with adjusted advection" );

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

   const std::shared_ptr< PrimitiveStorage >             storage_;
   FunctionType                                          cOld_;
   FunctionType                                          cTmp_;
   FunctionType                                          cPlus_;
   FunctionType                                          cMinus_;
   FunctionType                                          cAdjusted_;
   TimeSteppingScheme                                    timeSteppingSchemeConvection_;
   uint_t                                                numberOfCreatedParticles_;
   walberla::convection_particles::data::ParticleStorage particleStorage_;
   real_t                                                particleLocationRadius_;
};

} // namespace hyteg
