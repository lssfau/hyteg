/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "P2Function.hpp"

#include <algorithm>

#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/geometry/GeometryHelpers.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2MacroCell.hpp"
#include "hyteg/p2functionspace/P2MacroEdge.hpp"
#include "hyteg/p2functionspace/P2MacroFace.hpp"
#include "hyteg/p2functionspace/P2MacroVertex.hpp"
#include "hyteg/p2functionspace/P2Multigrid.hpp"
#include "hyteg/p2functionspace/P2TransferOperators.hpp"

namespace hyteg {

template < typename ValueType >
P2Function< ValueType >::P2Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage )
: Function< P2Function< ValueType > >( name, storage )
, vertexDoFFunction_( vertexdof::VertexDoFFunction< ValueType >( name + "_VertexDoF_dummy", storage ) )
, edgeDoFFunction_( EdgeDoFFunction< ValueType >( name + "__EdgeDoF_dummy", storage ) )
{}

template < typename ValueType >
P2Function< ValueType >::P2Function( const std::string&                         name,
                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                     uint_t                                     minLevel,
                                     uint_t                                     maxLevel )
: P2Function( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
{}

template < typename ValueType >
P2Function< ValueType >::P2Function( const std::string&                         name,
                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                     uint_t                                     minLevel,
                                     uint_t                                     maxLevel,
                                     BoundaryCondition                          boundaryCondition )
: Function< P2Function< ValueType > >( name, storage, minLevel, maxLevel )
, vertexDoFFunction_(
      vertexdof::VertexDoFFunction< ValueType >( name + "_VertexDoF", storage, minLevel, maxLevel, boundaryCondition ) )
, edgeDoFFunction_( EdgeDoFFunction< ValueType >( name + "_EdgeDoF", storage, minLevel, maxLevel, boundaryCondition ) )
{
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      /// one has to use the communicators of the vertexDoF and edgeDoF function to communicate
      /// TODO: find better solution
      communicators_[level] = nullptr;
   }
}

template < typename ValueType >
bool P2Function< ValueType >::evaluate( const Point3D& coordinates,
                                        uint_t         level,
                                        ValueType&     value,
                                        real_t         searchToleranceRadius ) const
{
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( coordinates );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( value );
      WALBERLA_UNUSED( searchToleranceRadius );
      WALBERLA_ABORT( "P2Function< ValueType >::evaluate not implemented for requested template parameter" );
   }
   else
   {
      if ( !this->getStorage()->hasGlobalCells() )
      {
         auto [found, faceID] = findFaceIDForPointIn2D( this->getStorage(), coordinates, searchToleranceRadius );
         if ( found )
         {
            value = P2::macroface::evaluate( level,
                                             *( this->getStorage()->getFace( faceID ) ),
                                             coordinates,
                                             vertexDoFFunction_.getFaceDataID(),
                                             edgeDoFFunction_.getFaceDataID() );
            return true;
         }
      }

      else
      {
         auto [found, cellID] = findCellIDForPointIn3D( this->getStorage(), coordinates, searchToleranceRadius );
         if ( found )
         {
            value = P2::macrocell::evaluate( level,
                                             *( this->getStorage()->getCell( cellID ) ),
                                             coordinates,
                                             vertexDoFFunction_.getCellDataID(),
                                             edgeDoFFunction_.getCellDataID() );
            return true;
         }
      }
   }

   // we place it here to make certain compilers happy; do you hear that intel?
   return false;
}

template < typename ValueType >
void P2Function< ValueType >::evaluateGradient( const Point3D& coordinates, uint_t level, Point3D& gradient ) const
{
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( coordinates );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( gradient );
      WALBERLA_ABORT( "P2Function< ValueType >::evaluateGradient not implemented for requested template parameter" );
   }
   else
   {
      // negative value would exclude this alternative feature in finding primitive ID
      real_t searchToleranceRadius = real_c( 1e-12 );

      // Check if 2D or 3D function
      if ( !this->getStorage()->hasGlobalCells() )
      {
         auto [found, faceID] = findFaceIDForPointIn2D( this->getStorage(), coordinates, searchToleranceRadius );
         if ( found )
         {
            P2::macroface::evaluateGradient( level,
                                             *( this->getStorage()->getFace( faceID ) ),
                                             coordinates,
                                             vertexDoFFunction_.getFaceDataID(),
                                             edgeDoFFunction_.getFaceDataID(),
                                             gradient );
            return;
         }
      }
      else
      {
         WALBERLA_ABORT( " P2Function< real_t >::evaluateGradient not implemented for 3D case" );
      }

      WALBERLA_ABORT( "There is no local macro element including a point at the given coordinates " << coordinates );
   }
}

template < typename ValueType >
void P2Function< ValueType >::interpolate( ValueType constant, uint_t level, DoFType flag ) const
{
   vertexDoFFunction_.interpolate( constant, level, flag );
   edgeDoFFunction_.interpolate( constant, level, flag );
}

template < typename ValueType >
void P2Function< ValueType >::interpolate( ValueType constant, uint_t level, BoundaryUID boundaryUID ) const
{
   vertexDoFFunction_.interpolate( constant, level, boundaryUID );
   edgeDoFFunction_.interpolate( constant, level, boundaryUID );
}

template < typename ValueType >
void P2Function< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                           uint_t                                              level,
                                           DoFType                                             flag ) const
{
   vertexDoFFunction_.interpolate( expr, level, flag );
   edgeDoFFunction_.interpolate( expr, level, flag );
}

template < typename ValueType >
void P2Function< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                           uint_t                                              level,
                                           BoundaryUID                                         boundaryUID ) const
{
   vertexDoFFunction_.interpolate( expr, level, boundaryUID );
   edgeDoFFunction_.interpolate( expr, level, boundaryUID );
}

template < typename ValueType >
void P2Function< ValueType >::interpolate(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< std::reference_wrapper< const P2Function< ValueType > > >&        srcFunctions,
    uint_t                                                                               level,
    DoFType                                                                              flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              edgeDoFFunctions;

   for ( const P2Function< ValueType >& function : srcFunctions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
   }

   vertexDoFFunction_.interpolate( expr, vertexDoFFunctions, level, flag );
   edgeDoFFunction_.interpolate( expr, edgeDoFFunctions, level, flag );
}

template < typename ValueType >
void P2Function< ValueType >::setToZero( const uint_t level ) const
{
   vertexDoFFunction_.setToZero( level );
   edgeDoFFunction_.setToZero( level );
}

template < typename ValueType >
void P2Function< ValueType >::swap( const P2Function< ValueType >& other, const uint_t& level, const DoFType& flag ) const
{
   vertexDoFFunction_.swap( other.getVertexDoFFunction(), level, flag );
   edgeDoFFunction_.swap( other.getEdgeDoFFunction(), level, flag );
}

template < typename ValueType >
void P2Function< ValueType >::copyFrom( const P2Function< ValueType >& other, const uint_t& level ) const
{
   vertexDoFFunction_.copyFrom( other.getVertexDoFFunction(), level );
   edgeDoFFunction_.copyFrom( other.getEdgeDoFFunction(), level );
}

template < typename ValueType >
void P2Function< ValueType >::copyFrom( const P2Function< ValueType >&         other,
                                        const uint_t&                          level,
                                        const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                                        const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const
{
   vertexDoFFunction_.copyFrom( other.getVertexDoFFunction(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
   edgeDoFFunction_.copyFrom( other.getEdgeDoFFunction(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
}

template < typename ValueType >
void P2Function< ValueType >::assign( const std::vector< ValueType >&                                               scalars,
                                      const std::vector< std::reference_wrapper< const P2Function< ValueType > > >& functions,
                                      uint_t                                                                        level,
                                      DoFType                                                                       flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              edgeDoFFunctions;

   for ( const P2Function< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
   }

   vertexDoFFunction_.assign( scalars, vertexDoFFunctions, level, flag );
   edgeDoFFunction_.assign( scalars, edgeDoFFunctions, level, flag );
}

template < typename ValueType >
void P2Function< ValueType >::add( ValueType scalar, uint_t level, DoFType flag ) const
{
   vertexDoFFunction_.add( scalar, level, flag );
   edgeDoFFunction_.add( scalar, level, flag );
}

template < typename ValueType >
void P2Function< ValueType >::add( const std::vector< ValueType >&                                               scalars,
                                   const std::vector< std::reference_wrapper< const P2Function< ValueType > > >& functions,
                                   uint_t                                                                        level,
                                   DoFType                                                                       flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              edgeDoFFunctions;

   for ( const P2Function< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
   }

   vertexDoFFunction_.add( scalars, vertexDoFFunctions, level, flag );
   edgeDoFFunction_.add( scalars, edgeDoFFunctions, level, flag );
}

template < typename ValueType >
void P2Function< ValueType >::multElementwise(
    const std::vector< std::reference_wrapper< const P2Function< ValueType > > >& functions,
    uint_t                                                                        level,
    DoFType                                                                       flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              edgeDoFFunctions;

   for ( const P2Function< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
   }

   vertexDoFFunction_.multElementwise( vertexDoFFunctions, level, flag );
   edgeDoFFunction_.multElementwise( edgeDoFFunctions, level, flag );
}

template < typename ValueType >
void P2Function< ValueType >::invertElementwise( uint_t level, DoFType flag, bool workOnHalos ) const
{
   vertexDoFFunction_.invertElementwise( level, flag, workOnHalos );
   edgeDoFFunction_.invertElementwise( level, flag, workOnHalos );
}

template < typename ValueType >
ValueType P2Function< ValueType >::dotGlobal( const P2Function< ValueType >& rhs, const uint_t level, const DoFType& flag ) const
{
   ValueType sum = dotLocal( rhs, level, flag );
   this->startTiming( "Dot (reduce)" );
   walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Dot (reduce)" );
   return sum;
}

template < typename ValueType >
ValueType P2Function< ValueType >::dotLocal( const P2Function< ValueType >& rhs, const uint_t level, const DoFType& flag ) const
{
   auto sum = ValueType( 0 );
   sum += vertexDoFFunction_.dotLocal( rhs.vertexDoFFunction_, level, flag );
   sum += edgeDoFFunction_.dotLocal( rhs.edgeDoFFunction_, level, flag );
   return sum;
}

template < typename ValueType >
ValueType P2Function< ValueType >::sumGlobal( const uint_t level, const DoFType& flag, const bool& absolute ) const
{
   ValueType sum = sumLocal( level, flag, absolute );
   this->startTiming( "Sum (reduce)" );
   walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Sum (reduce)" );
   return sum;
}

template < typename ValueType >
ValueType P2Function< ValueType >::sumLocal( const uint_t level, const DoFType& flag, const bool& absolute ) const
{
   auto sum = ValueType( 0 );
   sum += vertexDoFFunction_.sumLocal( level, flag, absolute );
   sum += edgeDoFFunction_.sumLocal( level, flag, absolute );
   return sum;
}

template < typename ValueType >
void P2Function< ValueType >::prolongateP1ToP2( const hyteg::P1Function< ValueType >& p1Function,
                                                const uint_t&                         level,
                                                const DoFType&                        flag ) const
{
   // Note: 'this' is the dst function - therefore we test this' boundary conditions

   this->startTiming( "Prolongate P1 -> P2" );

   p1Function.template startCommunication< Vertex, Edge >( level );
   p1Function.template startCommunication< Edge, Face >( level );

   for ( const auto& it : this->getStorage()->getVertices() )
   {
      const Vertex& vertex = *it.second;

      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         P2::macrovertex::prolongateP1ToP2< ValueType >( level,
                                                         vertex,
                                                         vertexDoFFunction_.getVertexDataID(),
                                                         edgeDoFFunction_.getVertexDataID(),
                                                         p1Function.getVertexDataID() );
      }
   }

   p1Function.template endCommunication< Vertex, Edge >( level );

   for ( const auto& it : this->getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::prolongateP1ToP2< ValueType >(
             level, edge, vertexDoFFunction_.getEdgeDataID(), edgeDoFFunction_.getEdgeDataID(), p1Function.getEdgeDataID() );
      }
   }

   p1Function.template endCommunication< Edge, Face >( level );

   for ( const auto& it : this->getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::prolongateP1ToP2< ValueType >(
             level, face, vertexDoFFunction_.getFaceDataID(), edgeDoFFunction_.getFaceDataID(), p1Function.getFaceDataID() );
      }
   }

   this->stopTiming( "Prolongate P1 -> P2" );
}

template < typename ValueType >
void P2Function< ValueType >::restrictP2ToP1( const P1Function< ValueType >& p1Function,
                                              const uint_t&                  level,
                                              const DoFType&                 flag ) const
{
   this->startTiming( "Restrict P2 -> P1" );

   vertexDoFFunction_.template startCommunication< Edge, Vertex >( level );
   edgeDoFFunction_.template startCommunication< Edge, Vertex >( level );

   vertexDoFFunction_.template startCommunication< Vertex, Edge >( level );
   edgeDoFFunction_.template startCommunication< Vertex, Edge >( level );

   vertexDoFFunction_.template startCommunication< Face, Edge >( level );
   edgeDoFFunction_.template startCommunication< Face, Edge >( level );

   vertexDoFFunction_.template startCommunication< Edge, Face >( level );
   edgeDoFFunction_.template startCommunication< Edge, Face >( level );

   vertexDoFFunction_.template endCommunication< Edge, Vertex >( level );
   edgeDoFFunction_.template endCommunication< Edge, Vertex >( level );

   for ( const auto& it : this->getStorage()->getVertices() )
   {
      const Vertex& vertex = *it.second;

      const DoFType vertexBC = p1Function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         P2::macrovertex::restrictP2ToP1< ValueType >( level,
                                                       vertex,
                                                       vertexDoFFunction_.getVertexDataID(),
                                                       edgeDoFFunction_.getVertexDataID(),
                                                       p1Function.getVertexDataID() );
      }
   }

   vertexDoFFunction_.template endCommunication< Vertex, Edge >( level );
   edgeDoFFunction_.template endCommunication< Vertex, Edge >( level );

   vertexDoFFunction_.template endCommunication< Face, Edge >( level );
   edgeDoFFunction_.template endCommunication< Face, Edge >( level );

   for ( const auto& it : this->getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      const DoFType edgeBC = p1Function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::restrictP2ToP1< ValueType >(
             level, edge, vertexDoFFunction_.getEdgeDataID(), edgeDoFFunction_.getEdgeDataID(), p1Function.getEdgeDataID() );
      }
   }

   vertexDoFFunction_.template endCommunication< Edge, Face >( level );
   edgeDoFFunction_.template endCommunication< Edge, Face >( level );

   for ( const auto& it : this->getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = p1Function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::restrictP2ToP1< ValueType >(
             level, face, vertexDoFFunction_.getFaceDataID(), edgeDoFFunction_.getFaceDataID(), p1Function.getFaceDataID() );
      }
   }

   this->stopTiming( "Restrict P2 -> P1" );
}

template < typename ValueType >
void P2Function< ValueType >::restrictInjection( uint_t sourceLevel, DoFType flag ) const
{
   for ( const auto& it : this->getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::restrictInjection< ValueType >(
             sourceLevel, face, vertexDoFFunction_.getFaceDataID(), edgeDoFFunction_.getFaceDataID() );
      }
   }

   for ( const auto& it : this->getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::restrictInjection< ValueType >(
             sourceLevel, edge, vertexDoFFunction_.getEdgeDataID(), edgeDoFFunction_.getEdgeDataID() );
      }
   }

   for ( const auto& it : this->getStorage()->getVertices() )
   {
      const Vertex& vertex = *it.second;

      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         P2::macrovertex::restrictInjection< ValueType >(
             sourceLevel, vertex, vertexDoFFunction_.getVertexDataID(), edgeDoFFunction_.getVertexDataID() );
      }
   }
}

template < typename ValueType >
ValueType P2Function< ValueType >::getMaxValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMax = -std::numeric_limits< ValueType >::max();
   localMax      = std::max( localMax, vertexDoFFunction_.getMaxValue( level, flag, false ) );
   localMax      = std::max( localMax, edgeDoFFunction_.getMaxValue( level, flag, false ) );
   walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType P2Function< ValueType >::getMaxMagnitude( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMax = ValueType( 0.0 );
   localMax      = std::max( localMax, vertexDoFFunction_.getMaxMagnitude( level, flag, false ) );
   localMax      = std::max( localMax, edgeDoFFunction_.getMaxMagnitude( level, flag, false ) );

   walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType P2Function< ValueType >::getMinValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMin = std::numeric_limits< ValueType >::max();
   localMin      = std::min( localMin, vertexDoFFunction_.getMinValue( level, flag, false ) );
   localMin      = std::min( localMin, edgeDoFFunction_.getMinValue( level, flag, false ) );
   walberla::mpi::allReduceInplace( localMin, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMin = localMin;
   if ( mpiReduce )
   {
      globalMin = -walberla::mpi::allReduce( -localMin, walberla::mpi::MAX );
   }

   return globalMin;
}

template < typename ValueType >
BoundaryCondition P2Function< ValueType >::getBoundaryCondition() const
{
   WALBERLA_ASSERT_EQUAL( vertexDoFFunction_.getBoundaryCondition(),
                          edgeDoFFunction_.getBoundaryCondition(),
                          "P2Function: boundary conditions of underlying vertex- and edgedof functions differ!" )
   return vertexDoFFunction_.getBoundaryCondition();
}

template < typename ValueType >
void P2Function< ValueType >::setBoundaryCondition( BoundaryCondition bc )
{
   vertexDoFFunction_.setBoundaryCondition( bc );
   edgeDoFFunction_.setBoundaryCondition( bc );
}

template < typename ValueType >
void P2Function< ValueType >::enumerate( uint_t level ) const
{
   this->startTiming( "Enumerate" );

   uint_t counterVertexDoFs = hyteg::numberOfLocalDoFs< VertexDoFFunctionTag >( *( this->getStorage() ), level );
   uint_t counterEdgeDoFs   = hyteg::numberOfLocalDoFs< EdgeDoFFunctionTag >( *( this->getStorage() ), level );

   std::vector< uint_t > vertexDoFsPerRank = walberla::mpi::allGather( counterVertexDoFs );
   std::vector< uint_t > edgeDoFsPerRank   = walberla::mpi::allGather( counterEdgeDoFs );

   ValueType offset = 0;

   for ( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
   {
      offset += static_cast< ValueType >( vertexDoFsPerRank[i] + edgeDoFsPerRank[i] );
   }
   enumerate( level, offset );
   this->stopTiming( "Enumerate" );
}

template < typename ValueType >
void P2Function< ValueType >::enumerate( uint_t level, ValueType& offset ) const
{
   vertexDoFFunction_.enumerate( level, offset );
   edgeDoFFunction_.enumerate( level, offset );
}

template < typename ValueType >
void P2Function< ValueType >::setLocalCommunicationMode(
    const communication::BufferedCommunicator::LocalCommunicationMode& localCommMode )
{
   vertexDoFFunction_.setLocalCommunicationMode( localCommMode );
   edgeDoFFunction_.setLocalCommunicationMode( localCommMode );
}

// ========================
//  explicit instantiation
// ========================
template class P2Function< double >;
// template class P2Function< float >;
template class P2Function< int32_t >;
template class P2Function< int64_t >;

namespace p2function {

void projectMean( const P2Function< real_t >& pressure, const uint_t& level )
{
   if ( pressure.isDummy() )
   {
      return;
   }
   const uint_t numGlobalVertices = numberOfGlobalDoFs< P2FunctionTag >( *pressure.getStorage(), level );
   const real_t sum               = pressure.sumGlobal( level, All );
   pressure.add( -sum / real_c( numGlobalVertices ), level, All );
}

} // namespace p2function

} //namespace hyteg
