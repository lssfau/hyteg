/*
 * Copyright (c) 2025 Ponsuganth Ilangovan
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

namespace hyteg {
template < typename FunctionType >
class MeshTransferWithParticles
{
 public:
   MeshTransferWithParticles()
   : particleStorage_( 10000 )
   {}

   void transfer( const FunctionType& src, const FunctionType& dst, const uint_t levelSrc, const uint_t levelDst )
   {
      const std::shared_ptr< PrimitiveStorage >& storageSrc = src.getStorage();
      const std::shared_ptr< PrimitiveStorage >& storageDst = dst.getStorage();

      particleStorage_.clear();

      const uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

      initialiseParticles( dst, levelDst, storageSrc );

      const uint_t numberOfCreatedParticles = particleStorage_.size();

      walberla::convection_particles::mpi::SyncNextNeighborsByPrimitiveID SNN;

      findSrcMacroWithParticle( *storageSrc, particleStorage_, 0.1 * MeshQuality::getMinimalEdgeLength( storageSrc, levelSrc ) );

      SNN( particleStorage_, *storageSrc );

      evaluateParticles( src, levelSrc );

      communicateParticles( dst, levelDst, numberOfCreatedParticles );
   }

 private:
   walberla::convection_particles::data::ParticleStorage particleStorage_;

   inline void evaluateParticles( const FunctionType& src, const uint_t levelSrc )
   {
      const std::shared_ptr< PrimitiveStorage >& storageSrc = src.getStorage();

      for ( auto p : particleStorage_ )
      {
         if ( storageSrc->hasGlobalCells() )
         {
            WALBERLA_CHECK( storageSrc->cellExistsLocally( p.getContainingPrimitive() ) );
            Cell&   cell = *( storageSrc->getCell( p.getContainingPrimitive() ) );
            Point3D computationalLocation;
            Point3D positionToEvaluate = toPoint3D( p.getPosition() );
            cell.getGeometryMap()->evalFinv( positionToEvaluate, computationalLocation );

            real_t finalValue = 0.0;

            if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
            {
               finalValue = vertexdof::macrocell::evaluate( levelSrc, cell, computationalLocation, src.getCellDataID() );
            }
            else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
            {
               finalValue = P2::macrocell::evaluate( levelSrc,
                                                     cell,
                                                     computationalLocation,
                                                     src.getVertexDoFFunction().getCellDataID(),
                                                     src.getEdgeDoFFunction().getCellDataID() );
            }
            else
            {
               WALBERLA_ABORT( "Unsupported function type" );
            }

            p->setFinalTemperature( finalValue );
         }
         else
         {
            WALBERLA_CHECK( storageSrc->faceExistsLocally( p.getContainingPrimitive() ) );
            Face&   face = *( storageSrc->getFace( p.getContainingPrimitive() ) );
            Point3D computationalLocation;
            Point3D positionToEvaluate = toPoint3D( p.getPosition() );
            face.getGeometryMap()->evalFinv( positionToEvaluate, computationalLocation );

            real_t finalValue = 0.0;

            if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
            {
               finalValue = vertexdof::macroface::evaluate( levelSrc, face, computationalLocation, src.getFaceDataID() );
            }
            else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
            {
               finalValue = P2::macroface::evaluate( levelSrc,
                                                     face,
                                                     computationalLocation,
                                                     src.getVertexDoFFunction().getFaceDataID(),
                                                     src.getEdgeDoFFunction().getFaceDataID() );
            }
            else
            {
               WALBERLA_ABORT( "Not implemented for this discretization." )
            }

            p->setFinalTemperature( finalValue );
         }
      }
   }

   inline void initialiseParticles( const FunctionType&                        dst,
                                    const uint_t                               levelDst,
                                    const std::shared_ptr< PrimitiveStorage >& storageSrc )
   {
      const std::shared_ptr< PrimitiveStorage >& storageDst = dst.getStorage();

      const uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

      if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
      {
         for ( auto it : FunctionIterator< vertexdof::VertexDoFFunction< real_t > >( dst.getVertexDoFFunction(), levelDst ) )
         {
            const PrimitiveID somePid = storageSrc->getPrimitiveIDs()[0];

            Point3D physicalLocation;
            auto    primitive = storageDst->getPrimitive( it.primitiveID() );
            primitive->getGeometryMap()->evalF( it.coordinates(), physicalLocation );

            auto particleIt = particleStorage_.create();
            particleIt->setOwner( (int) rank );
            particleIt->setPosition( toVec3( physicalLocation ) );
            particleIt->setStartPosition( toVec3( physicalLocation ) );
            particleIt->setStartDoFType( 0 ); // 0 == vertexdof
            particleIt->setStartEdgeDoFOrientation( it.edgeDoFOrientation() );
            particleIt->setStartPrimitiveID( it.primitiveID() );
            particleIt->setStartIndex( it.index() );
            particleIt->setStartProcess( rank );

            particleIt->setContainingPrimitive( somePid );

            // Update position
            storageDst->getPrimitive( it.primitiveID() )->getGeometryMap()->evalF( it.coordinates(), physicalLocation );
            particleIt->setPosition( toVec3( physicalLocation ) );
         }

         for ( auto it : FunctionIterator< EdgeDoFFunction< real_t > >( dst.getEdgeDoFFunction(), levelDst ) )
         {
            const PrimitiveID somePid = storageSrc->getPrimitiveIDs()[0];

            Point3D physicalLocation;
            auto    primitive = storageDst->getPrimitive( it.primitiveID() );
            primitive->getGeometryMap()->evalF( it.coordinates(), physicalLocation );

            auto particleIt = particleStorage_.create();
            particleIt->setOwner( (int) rank );
            particleIt->setPosition( toVec3( physicalLocation ) );
            particleIt->setStartPosition( toVec3( physicalLocation ) );
            particleIt->setStartDoFType( 1 ); // 1 == edgedof
            particleIt->setStartEdgeDoFOrientation( it.edgeDoFOrientation() );
            particleIt->setStartPrimitiveID( it.primitiveID() );
            particleIt->setStartIndex( it.index() );
            particleIt->setStartProcess( rank );

            particleIt->setContainingPrimitive( somePid );

            // Update position
            storageDst->getPrimitive( it.primitiveID() )->getGeometryMap()->evalF( it.coordinates(), physicalLocation );
            particleIt->setPosition( toVec3( physicalLocation ) );
         }
      }
      else if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
      {
         for ( auto it : FunctionIterator< vertexdof::VertexDoFFunction< real_t > >( dst, levelDst ) )
         {
            const PrimitiveID somePid = storageSrc->getPrimitiveIDs()[0];

            Point3D physicalLocation;
            auto    primitive = storageDst->getPrimitive( it.primitiveID() );
            primitive->getGeometryMap()->evalF( it.coordinates(), physicalLocation );

            auto particleIt = particleStorage_.create();
            particleIt->setOwner( (int) rank );
            particleIt->setPosition( toVec3( physicalLocation ) );
            particleIt->setStartPosition( toVec3( physicalLocation ) );
            particleIt->setStartDoFType( 0 ); // 0 == vertexdof
            particleIt->setStartEdgeDoFOrientation( it.edgeDoFOrientation() );
            particleIt->setStartPrimitiveID( it.primitiveID() );
            particleIt->setStartIndex( it.index() );
            particleIt->setStartProcess( rank );

            particleIt->setContainingPrimitive( somePid );

            // Update position
            storageDst->getPrimitive( it.primitiveID() )->getGeometryMap()->evalF( it.coordinates(), physicalLocation );
            particleIt->setPosition( toVec3( physicalLocation ) );
         }
      }
      else
      {
         WALBERLA_ABORT( "Not implemented for this discretisation" );
      }
   }

   inline void communicateParticles( const FunctionType& dst, const uint_t levelDst, const uint_t numberOfCreatedParticles )
   {
      const std::shared_ptr< PrimitiveStorage >& storageDst = dst.getStorage();

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
         static const int TAG1 = communication::MPITagProvider::getMPITag();

         for ( const auto& p : particleStorage_ )
         {
            if ( numParticlesToSendToRank.count( p.getStartProcess() ) == 0 )
            {
               numParticlesToSendToRank[p.getStartProcess()] = 0;
               sendRequests[p.getStartProcess()]             = MPI_Request();
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
                       TAG1,
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

            MPI_Recv(
                &numParticlesSum, 1, MPI_INT, MPI_ANY_SOURCE, TAG1, walberla::mpi::MPIManager::instance()->comm(), &status );

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
         numParticlesToSendToRank[0]      = (int) particleStorage_.size();
         numParticlesToReceiveFromRank[0] = (int) particleStorage_.size();
         ranksToReceiveFrom.insert( 0 );
      }

      static const int            TAG2 = communication::MPITagProvider::getMPITag();
      walberla::mpi::BufferSystem bufferSystem( walberla::mpi::MPIManager::instance()->comm(), TAG2 );

      for ( const auto& p : particleStorage_ )
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

            WALBERLA_CHECK( storageDst->primitiveExistsLocally( primitiveID ) );
            if ( storageDst->vertexExistsLocally( primitiveID ) )
            {
               auto vertex = storageDst->getVertex( primitiveID );
               WALBERLA_CHECK_EQUAL( dofType, 0 );
               if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
               {
                  vertex->getData( dst.getVertexDataID() )->getPointer( levelDst )[0] = temp;
               }
               else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
               {
                  vertex->getData( dst.getVertexDoFFunction().getVertexDataID() )->getPointer( levelDst )[0] = temp;
               }
               else
               {
                  WALBERLA_ABORT( "Not implemented for this discretization." );
               }
            }
            else if ( storageDst->edgeExistsLocally( primitiveID ) )
            {
               auto edge = storageDst->getEdge( primitiveID );
               if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
               {
                  WALBERLA_CHECK_EQUAL( dofType, 0 );
                  edge->getData( dst.getEdgeDataID() )
                      ->getPointer( levelDst )[vertexdof::macroedge::index( levelDst, index.x() )] = temp;
               }
               else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
               {
                  if ( dofType == 0 )
                  {
                     edge->getData( dst.getVertexDoFFunction().getEdgeDataID() )
                         ->getPointer( levelDst )[vertexdof::macroedge::index( levelDst, index.x() )] = temp;
                  }
                  else
                  {
                     edge->getData( dst.getEdgeDoFFunction().getEdgeDataID() )
                         ->getPointer( levelDst )[edgedof::macroedge::index( levelDst, index.x() )] = temp;
                  }
               }
               else
               {
                  WALBERLA_ABORT( "Not implemented for this discretization." );
               }
            }
            else if ( storageDst->faceExistsLocally( primitiveID ) )
            {
               auto face = storageDst->getFace( primitiveID );
               if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
               {
                  WALBERLA_CHECK_EQUAL( dofType, 0 );
                  face->getData( dst.getFaceDataID() )
                      ->getPointer( levelDst )[vertexdof::macroface::index( levelDst, index.x(), index.y() )] = temp;
               }
               else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
               {
                  if ( dofType == 0 )
                  {
                     face->getData( dst.getVertexDoFFunction().getFaceDataID() )
                         ->getPointer( levelDst )[vertexdof::macroface::index( levelDst, index.x(), index.y() )] = temp;
                  }
                  else
                  {
                     face->getData( dst.getEdgeDoFFunction().getFaceDataID() )
                         ->getPointer( levelDst )[edgedof::macroface::index( levelDst, index.x(), index.y(), orientation )] =
                         temp;
                  }
               }
               else
               {
                  WALBERLA_ABORT( "Not implemented for this discretization." );
               }
            }
            else if ( storageDst->cellExistsLocally( primitiveID ) )
            {
               auto cell = storageDst->getCell( primitiveID );
               if constexpr ( std::is_same< FunctionType, P1Function< real_t > >::value )
               {
                  WALBERLA_CHECK_EQUAL( dofType, 0 );
                  cell->getData( dst.getCellDataID() )
                      ->getPointer( levelDst )[vertexdof::macrocell::index( levelDst, index.x(), index.y(), index.z() )] = temp;
               }
               else if constexpr ( std::is_same< FunctionType, P2Function< real_t > >::value )
               {
                  if ( dofType == 0 )
                  {
                     cell->getData( dst.getVertexDoFFunction().getCellDataID() )
                         ->getPointer( levelDst )[vertexdof::macrocell::index( levelDst, index.x(), index.y(), index.z() )] =
                         temp;
                  }
                  else
                  {
                     cell->getData( dst.getEdgeDoFFunction().getCellDataID() )
                         ->getPointer(
                             levelDst )[edgedof::macrocell::index( levelDst, index.x(), index.y(), index.z(), orientation )] =
                         temp;
                  }
               }
               else
               {
                  WALBERLA_ABORT( "Not implemented for this discretization." );
               }
            }
            else
            {
               WALBERLA_ABORT( "Some PID not found" );
            }
         }
      }
   }

   inline void findSrcMacroWithParticle( const PrimitiveStorage&                                storage,
                                         walberla::convection_particles::data::ParticleStorage& particleStorage,
                                         const real_t                                           particleLocationRadius )
   {
      if ( storage.hasGlobalCells() )
      {
         for ( auto p : particleStorage )
         {
            p.setOutsideDomain( 0 );
            bool foundByPointLocation = false;

            std::vector< PrimitiveID > searchCellIDs;
            std::vector< PrimitiveID > neighborCellIDs;

            storage.getCellIDs( searchCellIDs );
            storage.getNeighboringCellIDs( neighborCellIDs );

            searchCellIDs.insert( searchCellIDs.end(), neighborCellIDs.begin(), neighborCellIDs.end() );

            for ( const auto& neighborCellID : searchCellIDs )
            {
               WALBERLA_ASSERT( storage.cellExistsLocally( neighborCellID ) ||
                                storage.cellExistsInNeighborhood( neighborCellID ) );
               const auto neighborCell = storage.getCell( neighborCellID );
               Point3D    computationalLocationNeighbor;
               neighborCell->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocationNeighbor );

               bool pointInNeighbourTetrahedronCheck = isPointInTetrahedron( computationalLocationNeighbor,
                                                                             neighborCell->getCoordinates().at( 0 ),
                                                                             neighborCell->getCoordinates().at( 1 ),
                                                                             neighborCell->getCoordinates().at( 2 ),
                                                                             neighborCell->getCoordinates().at( 3 ),
                                                                             neighborCell->getFaceInwardNormal( 0 ),
                                                                             neighborCell->getFaceInwardNormal( 1 ),
                                                                             neighborCell->getFaceInwardNormal( 2 ),
                                                                             neighborCell->getFaceInwardNormal( 3 ) );

               bool pointPairingNeighbourCheck = neighborCell->getGeometryMap()->verifyPointPairing(
                   computationalLocationNeighbor, toPoint3D( p->getPosition() ) );

               if ( pointInNeighbourTetrahedronCheck && pointPairingNeighbourCheck )
               {
                  // set it to the first neighbor we found to contain the particle
                  p->setContainingPrimitive( neighborCellID );
                  foundByPointLocation = true;
                  break;
               }
            }

            if ( !foundByPointLocation )
            {
               for ( const auto& neighborCellID : searchCellIDs )
               {
                  WALBERLA_ASSERT( storage.cellExistsLocally( neighborCellID ) ||
                                   storage.cellExistsInNeighborhood( neighborCellID ) );
                  const auto neighborCell = storage.getCell( neighborCellID );
                  Point3D    computationalLocationNeighbor;
                  neighborCell->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocationNeighbor );

                  bool sphereTetrahedronNeighbourIntersectionCheck =
                      sphereTetrahedronIntersection( computationalLocationNeighbor,
                                                     particleLocationRadius,
                                                     neighborCell->getCoordinates().at( 0 ),
                                                     neighborCell->getCoordinates().at( 1 ),
                                                     neighborCell->getCoordinates().at( 2 ),
                                                     neighborCell->getCoordinates().at( 3 ) );
                  //
                  // We remove the pointPairingNeighbourCheck here to reduce stringency,
                  // this could be useful for DoFs on the boundaries
                  // where some particles may not satisfy the point paring due to different blending maps
                  // but could still be within the sphere tetrahedron threshold
                  //
                  if ( sphereTetrahedronNeighbourIntersectionCheck )
                  {
                     p->setContainingPrimitive( neighborCellID );
                     foundByPointLocation = true;
                     break;
                  }
               }
            }

            if ( !foundByPointLocation )
            {
               p->setOutsideDomain( 1 );
               WALBERLA_ABORT( "Interpolated mesh transfer failed at "
                               << toPoint3D( p->getPosition() ) << " with norm = " << toPoint3D( p->getPosition() ).norm() );
            }
         }
      }
      else
      {
         for ( auto p : particleStorage )
         {
            p.setOutsideDomain( 0 );
            bool foundByPointLocation = false;

            std::vector< PrimitiveID > searchFaceIDs;
            std::vector< PrimitiveID > neighborFaceIDs;

            storage.getFaceIDs( searchFaceIDs );
            storage.getNeighboringFaceIDs( neighborFaceIDs );

            searchFaceIDs.insert( searchFaceIDs.end(), neighborFaceIDs.begin(), neighborFaceIDs.end() );

            for ( const auto& neighborFaceID : searchFaceIDs )
            {
               WALBERLA_ASSERT( storage.faceExistsLocally( neighborFaceID ) ||
                                storage.faceExistsInNeighborhood( neighborFaceID ) );
               const auto neighborFace = storage.getFace( neighborFaceID );
               Point3D    computationalLocationNeighbor;
               neighborFace->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocationNeighbor );
               Point2D computationalLocationNeighbor2D( computationalLocationNeighbor[0], computationalLocationNeighbor[1] );

               bool pointInNeighbourTriangleCheck = isPointInTriangle(
                   computationalLocationNeighbor2D,
                   Point2D( neighborFace->getCoordinates().at( 0 )[0], neighborFace->getCoordinates().at( 0 )[1] ),
                   Point2D( neighborFace->getCoordinates().at( 1 )[0], neighborFace->getCoordinates().at( 1 )[1] ),
                   Point2D( neighborFace->getCoordinates().at( 2 )[0], neighborFace->getCoordinates().at( 2 )[1] ) );

               bool pointPairingNeighbourCheck = neighborFace->getGeometryMap()->verifyPointPairing(
                   computationalLocationNeighbor, toPoint3D( p->getPosition() ) );

               if ( pointInNeighbourTriangleCheck && pointPairingNeighbourCheck )
               {
                  // set it to the first neighbor we found to contain the particle
                  p->setContainingPrimitive( neighborFaceID );
                  foundByPointLocation = true;
                  break;
               }
            }

            if ( !foundByPointLocation )
            {
               for ( const auto& neighborFaceID : searchFaceIDs )
               {
                  WALBERLA_ASSERT( storage.faceExistsLocally( neighborFaceID ) ||
                                   storage.faceExistsInNeighborhood( neighborFaceID ) );
                  const auto neighborFace = storage.getFace( neighborFaceID );
                  Point3D    computationalLocationNeighbor;
                  neighborFace->getGeometryMap()->evalFinv( toPoint3D( p->getPosition() ), computationalLocationNeighbor );

                  WALBERLA_LOG_INFO("neighborFace->getCoordinates().at( 0 ) = " << neighborFace->getCoordinates().at( 0 ));
                  WALBERLA_LOG_INFO("neighborFace->getCoordinates().at( 1 ) = " << neighborFace->getCoordinates().at( 1 ));
                  WALBERLA_LOG_INFO("neighborFace->getCoordinates().at( 2 ) = " << neighborFace->getCoordinates().at( 2 ));

                  bool sphereTriangleNeighbourIntersectionCheck =
                      sphereTriangleIntersection( computationalLocationNeighbor,
                                                  particleLocationRadius,
                                                  neighborFace->getCoordinates().at( 0 ),
                                                  neighborFace->getCoordinates().at( 1 ),
                                                  neighborFace->getCoordinates().at( 2 ) );

                  bool pointPairingNeighbourCheck = neighborFace->getGeometryMap()->verifyPointPairing(
                      computationalLocationNeighbor, toPoint3D( p->getPosition() ) );

                  if ( sphereTriangleNeighbourIntersectionCheck )
                  {
                     p->setContainingPrimitive( neighborFaceID );
                     foundByPointLocation = true;
                     break;
                  }
               }
            }

            if ( !foundByPointLocation )
            {
               p->setOutsideDomain( 1 );
               WALBERLA_ABORT( "Interpolated mesh transfer failed at "
                               << toPoint3D( p->getPosition() ) << " with norm = " << toPoint3D( p->getPosition() ).norm() );
            }
         }
      }
   }
};
} // namespace hyteg
