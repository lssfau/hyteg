/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "P2PlusBubbleFunction.hpp"

#include <algorithm>

#include "hyteg/geometry/BlendingHelpers.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/p2functionspace/P2MacroCell.hpp"
#include "hyteg/p2functionspace/P2MacroFace.hpp"
#include "hyteg/p2functionspace/P2Multigrid.hpp"
#include "hyteg/p2functionspace/P2TransferOperators.hpp"

namespace hyteg {

using volumedofspace::VolumeDoFFunction;

template < typename ValueType >
P2PlusBubbleFunction< ValueType >::P2PlusBubbleFunction( const std::string&                         name,
                                                         const std::shared_ptr< PrimitiveStorage >& storage )
: Function< P2PlusBubbleFunction< ValueType > >( name, storage )
, vertexDoFFunction_( vertexdof::VertexDoFFunction< ValueType >( name + "_VertexDoF_dummy", storage ) )
, edgeDoFFunction_( EdgeDoFFunction< ValueType >( name + "__EdgeDoF_dummy", storage ) )
, volumeDoFFunction_( VolumeDoFFunction< ValueType >( name + "__VolumeDoF_dummy",
                                                      storage,
                                                      1,
                                                      1,
                                                      1,
                                                      volumedofspace::indexing::VolumeDoFMemoryLayout::AoS ) )
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
}

template < typename ValueType >
P2PlusBubbleFunction< ValueType >::P2PlusBubbleFunction( const std::string&                         name,
                                                         const std::shared_ptr< PrimitiveStorage >& storage,
                                                         uint_t                                     minLevel,
                                                         uint_t                                     maxLevel )
: P2PlusBubbleFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
{}

template < typename ValueType >
P2PlusBubbleFunction< ValueType >::P2PlusBubbleFunction( const std::string&                         name,
                                                         const std::shared_ptr< PrimitiveStorage >& storage,
                                                         uint_t                                     minLevel,
                                                         uint_t                                     maxLevel,
                                                         BoundaryCondition                          boundaryCondition )
: Function< P2PlusBubbleFunction< ValueType > >( name, storage, minLevel, maxLevel )
, vertexDoFFunction_(
      vertexdof::VertexDoFFunction< ValueType >( name + "_VertexDoF", storage, minLevel, maxLevel, boundaryCondition ) )
, edgeDoFFunction_( EdgeDoFFunction< ValueType >( name + "_EdgeDoF", storage, minLevel, maxLevel, boundaryCondition ) )
, volumeDoFFunction_( VolumeDoFFunction< ValueType >( name + "_VolumeDoF",
                                                      storage,
                                                      minLevel,
                                                      maxLevel,
                                                      1,
                                                      volumedofspace::indexing::VolumeDoFMemoryLayout::AoS ) )
{
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      /// one has to use the communicators of the vertexDoF and edgeDoF function to communicate
      /// TODO: find better solution
      communicators_[level] = nullptr;
   }
}

template < typename ValueType >
bool P2PlusBubbleFunction< ValueType >::evaluate( const Point3D& physicalCoords,
                                                  uint_t         level,
                                                  ValueType&     value,
                                                  real_t         searchToleranceRadius ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( physicalCoords );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( value );
      WALBERLA_UNUSED( searchToleranceRadius );
      WALBERLA_ABORT( "P2PlusBubbleFunction< ValueType >::evaluate not implemented for requested template parameter" );
      return false;
   }
   else
   {
      if ( !this->getStorage()->hasGlobalCells() )
      {
         auto [found, faceID, computationalCoords] =
             mapFromPhysicalToComputationalDomain2D( this->getStorage(), physicalCoords, searchToleranceRadius );
         if ( found )
         {
            value = P2::macroface::evaluate( level,
                                             *( this->getStorage()->getFace( faceID ) ),
                                             computationalCoords,
                                             vertexDoFFunction_.getFaceDataID(),
                                             edgeDoFFunction_.getFaceDataID() );
            return true;
         }
      }

      else
      {
         auto [found, cellID, computationalCoords] =
             mapFromPhysicalToComputationalDomain3D( this->getStorage(), physicalCoords, searchToleranceRadius );
         if ( found )
         {
            value = P2::macrocell::evaluate( level,
                                             *( this->getStorage()->getCell( cellID ) ),
                                             computationalCoords,
                                             vertexDoFFunction_.getCellDataID(),
                                             edgeDoFFunction_.getCellDataID() );
            return true;
         }
      }

      // no match found
      return false;
   }

   // will not be reached, but some compilers complain otherwise
   return false;
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::evaluateGradient( const Point3D& physicalCoords, uint_t level, Point3D& gradient ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( physicalCoords );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( gradient );
      WALBERLA_ABORT( "P2PlusBubbleFunction< ValueType >::evaluateGradient not implemented for requested template parameter" );
   }
   else
   {
      // negative value would exclude this alternative feature in finding primitive ID
      real_t searchToleranceRadius = real_c( 1e-12 );

      // Check if 2D or 3D function
      if ( !this->getStorage()->hasGlobalCells() )
      {
         auto [found, faceID, computationalCoords] =
             mapFromPhysicalToComputationalDomain2D( this->getStorage(), physicalCoords, searchToleranceRadius );
         if ( found )
         {
            Face& face = *( this->getStorage()->getFace( faceID ) );

            // evaluate gradient on computational domain
            P2::macroface::evaluateGradient( level,
                                             face,
                                             computationalCoords,
                                             vertexDoFFunction_.getFaceDataID(),
                                             edgeDoFFunction_.getFaceDataID(),
                                             gradient );

            // transform gradient to physical coordinates
            Matrix2r DFinv;
            face.getGeometryMap()->evalDFinv( physicalCoords, DFinv );
            real_t aux0 = gradient[0];
            real_t aux1 = gradient[1];
            gradient[0] = DFinv( 0, 0 ) * aux0 + DFinv( 0, 1 ) * aux1;
            gradient[1] = DFinv( 1, 0 ) * aux0 + DFinv( 1, 1 ) * aux1;

            return;
         }
      }
      else
      {
         WALBERLA_ABORT( "P2PlusBubbleFunction< real_t >::evaluateGradient not implemented for 3D case" );
      }

      WALBERLA_ABORT( "There is no local macro element including a point at the given mapped back coordinates for "
                      << physicalCoords );
   }
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::interpolate( ValueType constant, uint_t level, DoFType flag ) const
{
   vertexDoFFunction_.interpolate( constant, level, flag );
   edgeDoFFunction_.interpolate( constant, level, flag );

   // The meaning of the BubbleDoF is different. We cannot get its value by simply using the supplied
   // constant value. The reason is that the other shape functions will not be zero at the barycenter
   // of the element. Thus, we need to correct for this. In the constant case, the bubble DoF is not
   // needed, as the P2 part already interpolates the constant function exactly.
   volumeDoFFunction_.setToZero( level );
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::interpolate( ValueType constant, uint_t level, BoundaryUID boundaryUID ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   vertexDoFFunction_.interpolate( constant, level, boundaryUID );
   edgeDoFFunction_.interpolate( constant, level, boundaryUID );

   // need to check whether face/cell has the given boundaryUID (which despite the name need not flag a boundary)
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                     uint_t                                              level,
                                                     DoFType                                             flag ) const
{
   if constexpr ( !std::is_floating_point_v< ValueType > )
   {
      WALBERLA_ABORT( "P2PlusBubbleFunction::interpolate() only works with FP-types!" );
   }
   else
   {
      vertexDoFFunction_.interpolate( expr, level, flag );
      edgeDoFFunction_.interpolate( expr, level, flag );

      // For 3D we work on cells and for 2D on faces
      if ( this->storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "Sorry, but 3D interpolation not implemented, yet, for P2PlusBubbleFunction!" );
         // we only perform computations on cell primitives
         // for ( auto& macroIter : this->storage_->getCells() )
         {
            // Cell& cell = *macroIter.second;
            // bubbleData[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayout )] =
         }
      }

      // 2D
      else
      {
         // The meaning of the BubbleDoF is different. We cannot get its value by simply using the supplied
         // constant value. The reason is that the other shape functions will not be zero at the barycenter
         // (x_b,y_b) of the element. We have
         //
         // P_2^+(x_b,y_b) = w_{\text{bubble}} - \frac{1}{9} \Big[ f(x_1,y_1) + f(x_2,y_2) + f(x_3,y_3) \Big]
         //                  + \frac{4}{9} \Big[ f(x_4,y_4) + f(x_5,y_5) + f(x_6,y_6) \Big]
         //                \stackrel{!}{=} f((x_b,y_b)
         //
         // where f() represents the expression, (x_k, y_k) the coordinates of the vertex and edge dof positions.
         // We solve this for the bubble DoF w_{\text{bubble}}
         // we only perform computations on face primitives
         for ( auto& it : this->storage_->getFaces() )
         {
            const auto  faceID = it.first;
            const Face& face   = *it.second;

            ValueType* vertexData = face.getData( vertexDoFFunction_.getFaceDataID() )->getPointer( level );
            ValueType* edgeData   = face.getData( edgeDoFFunction_.getFaceDataID() )->getPointer( level );

            auto       bubbleData = volumeDoFFunction_.dofMemory( faceID, level );
            const auto memLayout  = volumeDoFFunction_.memoryLayout();

            // loop over micro-faces
            for ( const auto& faceType : facedof::allFaceTypes )
            {
               for ( const auto& microFace : facedof::macroface::Iterator( level, faceType, 0 ) )
               {
                  // obtain data indices of dofs associated with vertices and edges of micro-face
                  std::array< uint_t, 3 > vertexDoFIndices;
                  vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, level, vertexDoFIndices );

                  std::array< uint_t, 3 > edgeDoFIndices;
                  edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, faceType, level, edgeDoFIndices );

                  // assemble "correcton term"
                  ValueType vertexDoFsum = ValueType( 0 );
                  vertexDoFsum += vertexData[vertexDoFIndices[0]];
                  vertexDoFsum += vertexData[vertexDoFIndices[1]];
                  vertexDoFsum += vertexData[vertexDoFIndices[2]];

                  ValueType edgeDoFsum = ValueType( 0 );
                  edgeDoFsum += edgeData[edgeDoFIndices[0]];
                  edgeDoFsum += edgeData[edgeDoFIndices[1]];
                  edgeDoFsum += edgeData[edgeDoFIndices[2]];

                  ValueType correction = ( vertexDoFsum - ValueType( 4 ) * edgeDoFsum ) / ValueType( 9 );

                  // obtain location of barycenter mapped to physical domain
                  Point3D barycenter = micromesh::microFaceCenterPosition( this->storage_, faceID, level, microFace, faceType  );

                  // evaluate expression at this center and add the correction
                  bubbleData[volumedofspace::indexing::index( microFace.x(), microFace.y(), faceType, 0, 1, level, memLayout )] =
                      expr( barycenter ) + correction;
               }
            }
         }
      }
   }
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                     uint_t                                              level,
                                                     BoundaryUID                                         boundaryUID ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   vertexDoFFunction_.interpolate( expr, level, boundaryUID );
   edgeDoFFunction_.interpolate( expr, level, boundaryUID );
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::interpolate(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >&    expr,
    const std::vector< std::reference_wrapper< const P2PlusBubbleFunction< ValueType > > >& srcFunctions,
    uint_t                                                                                  level,
    DoFType                                                                                 flag ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              edgeDoFFunctions;

   for ( const P2PlusBubbleFunction< ValueType >& function : srcFunctions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
   }

   vertexDoFFunction_.interpolate( expr, vertexDoFFunctions, level, flag );
   edgeDoFFunction_.interpolate( expr, edgeDoFFunctions, level, flag );
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::setToZero( const uint_t level ) const
{
   vertexDoFFunction_.setToZero( level );
   edgeDoFFunction_.setToZero( level );
   volumeDoFFunction_.setToZero( level );
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::swap( const P2PlusBubbleFunction< ValueType >& other,
                                              const uint_t&                            level,
                                              const DoFType&                           flag ) const
{
   vertexDoFFunction_.swap( other.getVertexDoFFunction(), level, flag );
   edgeDoFFunction_.swap( other.getEdgeDoFFunction(), level, flag );

   // bubble DoFs are always Inner
   if ( ( flag | DoFType::Inner ) == DoFType::Inner )
   {
      volumeDoFFunction_.swap( other.getVolumeDoFFunction(), level );
   }
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::copyFrom( const P2PlusBubbleFunction< ValueType >& other, const uint_t& level ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   vertexDoFFunction_.copyFrom( other.getVertexDoFFunction(), level );
   edgeDoFFunction_.copyFrom( other.getEdgeDoFFunction(), level );
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::copyFrom( const P2PlusBubbleFunction< ValueType >& other,
                                                  const uint_t&                            level,
                                                  const std::map< PrimitiveID, uint_t >&   localPrimitiveIDsToRank,
                                                  const std::map< PrimitiveID, uint_t >&   otherPrimitiveIDsToRank ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   vertexDoFFunction_.copyFrom( other.getVertexDoFFunction(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
   edgeDoFFunction_.copyFrom( other.getEdgeDoFFunction(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::assign(
    const std::vector< ValueType >&                                                         scalars,
    const std::vector< std::reference_wrapper< const P2PlusBubbleFunction< ValueType > > >& functions,
    uint_t                                                                                  level,
    DoFType                                                                                 flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              edgeDoFFunctions;
   std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > >            volumeDoFFunctions;

   for ( const P2PlusBubbleFunction< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
      volumeDoFFunctions.push_back( function.volumeDoFFunction_ );
   }

   vertexDoFFunction_.assign( scalars, vertexDoFFunctions, level, flag );
   edgeDoFFunction_.assign( scalars, edgeDoFFunctions, level, flag );

   // bubble DoFs are always Inner
   if ( ( flag | DoFType::Inner ) == DoFType::Inner )
   {
      volumeDoFFunction_.assign( scalars, volumeDoFFunctions, level );
   }
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::add( ValueType scalar, uint_t level, DoFType flag ) const
{
   vertexDoFFunction_.add( scalar, level, flag );
   edgeDoFFunction_.add( scalar, level, flag );

   // bubble DoFs are always Inner
   if ( ( flag | DoFType::Inner ) == DoFType::Inner )
   {
      volumeDoFFunction_.add( scalar, level );
   }
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::add(
    const std::vector< ValueType >&                                                         scalars,
    const std::vector< std::reference_wrapper< const P2PlusBubbleFunction< ValueType > > >& functions,
    uint_t                                                                                  level,
    DoFType                                                                                 flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              edgeDoFFunctions;
   std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > >            volumeDoFFunctions;

   for ( const P2PlusBubbleFunction< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
      volumeDoFFunctions.push_back( function.volumeDoFFunction_ );
   }

   vertexDoFFunction_.add( scalars, vertexDoFFunctions, level, flag );
   edgeDoFFunction_.add( scalars, edgeDoFFunctions, level, flag );

   // bubble DoFs are always Inner
   if ( ( flag | DoFType::Inner ) == DoFType::Inner )
   {
      volumeDoFFunction_.add( scalars, volumeDoFFunctions, level );
   }
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::multElementwise(
    const std::vector< std::reference_wrapper< const P2PlusBubbleFunction< ValueType > > >& functions,
    uint_t                                                                                  level,
    DoFType                                                                                 flag ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              edgeDoFFunctions;

   for ( const P2PlusBubbleFunction< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
   }

   vertexDoFFunction_.multElementwise( vertexDoFFunctions, level, flag );
   edgeDoFFunction_.multElementwise( edgeDoFFunctions, level, flag );

   // bubble DoFs are always Inner
   if ( ( flag | DoFType::Inner ) == DoFType::Inner )
   {
      // function below now available, yet!
      // volumeDoFFunction_.multElementwise( volumeDoFFunctions, level, flag );
   }
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::invertElementwise( uint_t level, DoFType flag, bool workOnHalos ) const
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   vertexDoFFunction_.invertElementwise( level, flag, workOnHalos );
   edgeDoFFunction_.invertElementwise( level, flag, workOnHalos );
}

template < typename ValueType >
ValueType P2PlusBubbleFunction< ValueType >::dotGlobal( const P2PlusBubbleFunction< ValueType >& rhs,
                                                        const uint_t                             level,
                                                        const DoFType&                           flag ) const
{
   ValueType sum = dotLocal( rhs, level, flag );
   this->startTiming( "Dot (reduce)" );
   walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Dot (reduce)" );
   return sum;
}

template < typename ValueType >
ValueType P2PlusBubbleFunction< ValueType >::dotLocal( const P2PlusBubbleFunction< ValueType >& rhs,
                                                       const uint_t                             level,
                                                       const DoFType&                           flag ) const
{
   auto sum = ValueType( 0 );
   sum += vertexDoFFunction_.dotLocal( rhs.vertexDoFFunction_, level, flag );
   sum += edgeDoFFunction_.dotLocal( rhs.edgeDoFFunction_, level, flag );

   // bubble DoFs are always Inner
   if ( ( flag | DoFType::Inner ) == DoFType::Inner )
   {
      sum += volumeDoFFunction_.dotLocal( rhs.volumeDoFFunction_, level );
   }

   return sum;
}

template < typename ValueType >
ValueType P2PlusBubbleFunction< ValueType >::sumGlobal( const uint_t level, const DoFType& flag, const bool& absolute ) const
{
   ValueType sum = sumLocal( level, flag, absolute );
   this->startTiming( "Sum (reduce)" );
   walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Sum (reduce)" );
   return sum;
}

template < typename ValueType >
ValueType P2PlusBubbleFunction< ValueType >::sumLocal( const uint_t level, const DoFType& flag, const bool& absolute ) const
{
   auto sum = ValueType( 0 );
   sum += vertexDoFFunction_.sumLocal( level, flag, absolute );
   sum += edgeDoFFunction_.sumLocal( level, flag, absolute );

   // bubble DoFs are always Inner
   if ( ( flag | DoFType::Inner ) == DoFType::Inner )
   {
      sum += volumeDoFFunction_.sumLocal( level, absolute );
   }

   return sum;
}

template < typename ValueType >
ValueType P2PlusBubbleFunction< ValueType >::getMaxDoFValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMax = -std::numeric_limits< ValueType >::max();
   localMax      = std::max( localMax, vertexDoFFunction_.getMaxDoFValue( level, flag, false ) );
   localMax      = std::max( localMax, edgeDoFFunction_.getMaxDoFValue( level, flag, false ) );
   localMax      = std::max( localMax, volumeDoFFunction_.getMaxDoFValue( level, false ) );
   walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType P2PlusBubbleFunction< ValueType >::getMaxDoFMagnitude( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMax = ValueType( 0.0 );
   localMax      = std::max( localMax, vertexDoFFunction_.getMaxDoFMagnitude( level, flag, false ) );
   localMax      = std::max( localMax, edgeDoFFunction_.getMaxDoFMagnitude( level, flag, false ) );
   localMax      = std::max( localMax, volumeDoFFunction_.getMaxDoFMagnitude( level, false ) );

   walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType P2PlusBubbleFunction< ValueType >::getMinDoFValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMin = std::numeric_limits< ValueType >::max();
   localMin      = std::min( localMin, vertexDoFFunction_.getMinDoFValue( level, flag, false ) );
   localMin      = std::min( localMin, edgeDoFFunction_.getMinDoFValue( level, flag, false ) );
   localMin      = std::min( localMin, volumeDoFFunction_.getMinDoFValue( level, false ) );

   walberla::mpi::allReduceInplace( localMin, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMin = localMin;
   if ( mpiReduce )
   {
      globalMin = -walberla::mpi::allReduce( -localMin, walberla::mpi::MAX );
   }

   return globalMin;
}

template < typename ValueType >
BoundaryCondition P2PlusBubbleFunction< ValueType >::getBoundaryCondition() const
{
   WALBERLA_ASSERT_EQUAL( vertexDoFFunction_.getBoundaryCondition(),
                          edgeDoFFunction_.getBoundaryCondition(),
                          "P2PlusBubbleFunction: boundary conditions of underlying vertex- and edgedof functions differ!" )
   return vertexDoFFunction_.getBoundaryCondition();
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::setBoundaryCondition( BoundaryCondition bc )
{
   vertexDoFFunction_.setBoundaryCondition( bc );
   edgeDoFFunction_.setBoundaryCondition( bc );
   // BCs not set for the bubble DoFs as not relevant
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::enumerate( uint_t level ) const
{
   this->startTiming( "Enumerate" );

   uint_t counterVertexDoFs = hyteg::numberOfLocalDoFs< VertexDoFFunctionTag >( *( this->getStorage() ), level );
   uint_t counterEdgeDoFs   = hyteg::numberOfLocalDoFs< EdgeDoFFunctionTag >( *( this->getStorage() ), level );
   uint_t counterVolumeDoFs = 0u;

   if ( this->getStorage()->hasGlobalCells() )
   {
      uint_t numCells   = this->getStorage()->getNumberOfLocalCells();
      counterVolumeDoFs = numCells * levelinfo::num_microcells_per_cell( level );
   }
   else
   {
      uint_t numFaces   = this->getStorage()->getNumberOfLocalFaces();
      counterVolumeDoFs = numFaces * levelinfo::num_microfaces_per_face( level );
   }

   std::vector< uint_t > vertexDoFsPerRank = walberla::mpi::allGather( counterVertexDoFs );
   std::vector< uint_t > edgeDoFsPerRank   = walberla::mpi::allGather( counterEdgeDoFs );
   std::vector< uint_t > volumeDoFsPerRank = walberla::mpi::allGather( counterVolumeDoFs );

   ValueType offset = 0;

   for ( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
   {
      offset += static_cast< ValueType >( vertexDoFsPerRank[i] + edgeDoFsPerRank[i] + volumeDoFsPerRank[i] );
   }
   enumerate( level, offset );
   this->stopTiming( "Enumerate" );
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::enumerate( uint_t level, ValueType& offset ) const
{
   vertexDoFFunction_.enumerate( level, offset );
   edgeDoFFunction_.enumerate( level, offset );

   // We need to treat the bubble dofs ourselves, since the VolumeDoFFunction
   // does not offer an enumerate()
   if ( this->storage_->hasGlobalCells() )
   {
      for ( const auto& it : this->storage_->getCells() )
      {
         const auto cellID = it.first;
         const auto cell   = *it.second;

         auto       dofs      = volumeDoFFunction_.dofMemory( cellID, level );
         const auto memLayout = volumeDoFFunction_.memoryLayout();

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
            {
               dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayout )] =
                   offset++;
            }
         }
      }
   }
   else
   {
      for ( const auto& it : this->storage_->getFaces() )
      {
         const auto faceID = it.first;
         const auto face   = *it.second;

         auto       dofs      = volumeDoFFunction_.dofMemory( faceID, level );
         const auto memLayout = volumeDoFFunction_.memoryLayout();

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayout )] = offset++;
            }
         }
      }
   }
}

template < typename ValueType >
void P2PlusBubbleFunction< ValueType >::setLocalCommunicationMode(
    const communication::BufferedCommunicator::LocalCommunicationMode& localCommMode )
{
   WALBERLA_ABORT( "Member function not bubble ready!" );
   vertexDoFFunction_.setLocalCommunicationMode( localCommMode );
   edgeDoFFunction_.setLocalCommunicationMode( localCommMode );
}

// ========================
//  explicit instantiation
// ========================
template class P2PlusBubbleFunction< double >;
template class P2PlusBubbleFunction< float >;
template class P2PlusBubbleFunction< int32_t >;
template class P2PlusBubbleFunction< idx_t >;

} //namespace hyteg
