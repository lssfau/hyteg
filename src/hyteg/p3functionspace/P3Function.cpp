/*
 * Copyright (c) 2026 Marcus Mohr.
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

#include "P3Function.hpp"

#include <algorithm>

#include "hyteg/p3functionspace/P3FunctionPackInfo.hpp"

namespace hyteg {

template < typename ValueType >
P3Function< ValueType >::P3Function( const std::string&                         name,
                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                     uint_t                                     minLevel,
                                     uint_t                                     maxLevel )
: P3Function( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
{}

template < typename ValueType >
P3Function< ValueType >::P3Function( const std::string&                         name,
                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                     uint_t                                     minLevel,
                                     uint_t                                     maxLevel,
                                     BoundaryCondition                          boundaryCondition )
: Function< P3Function< ValueType > >( name, storage, minLevel, maxLevel )
, vertexDoFFunction_(
      vertexdof::VertexDoFFunction< ValueType >( name + "_VertexDoF", storage, minLevel, maxLevel, boundaryCondition ) )
, edgeDoFFunctionBlue_( EdgeDoFFunction< ValueType >( name + "_EdgeDoFBlue", storage, minLevel, maxLevel, boundaryCondition ) )
, edgeDoFFunctionRed_( EdgeDoFFunction< ValueType >( name + "_EdgeDoFRed", storage, minLevel, maxLevel, boundaryCondition ) )
, faceDoFFunction_( volumedofspace::VolumeDoFFunction< ValueType >( name + "_FaceDoF",
                                                                    storage,
                                                                    minLevel,
                                                                    maxLevel,
                                                                    1,
                                                                    volumedofspace::indexing::VolumeDoFMemoryLayout::SoA ) )
{
   if ( storage->hasGlobalCells() )
   {
      WALBERLA_ABORT( "This proof-of-concept implementation only works for 2D! We are missing a FaceDoFFunction for 3D!" );
   }

   std::array< PrimitiveDataID< FunctionMemory< ValueType >, Vertex >, 3 > dataIDsMacroVertex;
   std::array< PrimitiveDataID< FunctionMemory< ValueType >, Edge >, 3 >   dataIDsMacroEdge;
   std::array< PrimitiveDataID< FunctionMemory< ValueType >, Face >, 3 >   dataIDsMacroFace;
   std::array< PrimitiveDataID< FunctionMemory< ValueType >, Cell >, 3 >   dataIDsMacroCell; // currently unused

   dataIDsMacroVertex[0] = vertexDoFFunction_.getVertexDataID();
   dataIDsMacroVertex[1] = edgeDoFFunctionBlue_.getVertexDataID();
   dataIDsMacroVertex[2] = edgeDoFFunctionRed_.getVertexDataID();

   dataIDsMacroEdge[0] = vertexDoFFunction_.getEdgeDataID();
   dataIDsMacroEdge[1] = edgeDoFFunctionBlue_.getEdgeDataID();
   dataIDsMacroEdge[2] = edgeDoFFunctionRed_.getEdgeDataID();

   dataIDsMacroFace[0] = vertexDoFFunction_.getFaceDataID();
   dataIDsMacroFace[1] = edgeDoFFunctionBlue_.getFaceDataID();
   dataIDsMacroFace[2] = edgeDoFFunctionRed_.getFaceDataID();

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      communicators_[level]->addPackInfo( std::make_shared< P3FunctionPackInfo< ValueType > >(
          level, dataIDsMacroVertex, dataIDsMacroEdge, dataIDsMacroFace, dataIDsMacroCell, this->getStorage() ) );
   }
}

template < typename ValueType >
void P3Function< ValueType >::setToZero( const uint_t level ) const
{
   vertexDoFFunction_.setToZero( level );
   edgeDoFFunctionRed_.setToZero( level );
   edgeDoFFunctionBlue_.setToZero( level );
   faceDoFFunction_.setToZero( level );
};

template < typename ValueType >
void P3Function< ValueType >::interpolate( ValueType constant, uint_t level, DoFType flag ) const
{
   vertexDoFFunction_.interpolate( constant, level, flag );
   edgeDoFFunctionBlue_.interpolate( constant, level, flag );
   edgeDoFFunctionRed_.interpolate( constant, level, flag );
   faceDoFHelpers::interpolate( faceDoFFunction_, constant, level, flag );
}

template < typename ValueType >
void P3Function< ValueType >::faceDoFHelpers::interpolate( const volumedofspace::VolumeDoFFunction< ValueType >& function,
                                                           ValueType                                             constant,
                                                           uint_t                                                level,
                                                           DoFType                                               flag )
{
   // NOTE: We pass "flag" here, as this will become important in 3D; for 2D it is not used
   WALBERLA_UNUSED( flag );

   if ( function.getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "No 3D support in P3Function, yet!" );
   }
   else
   {
      for ( auto& it : function.getStorage()->getFaces() )
      {
         const auto  faceID = it.first;
         const auto& face   = *it.second;

         WALBERLA_CHECK_EQUAL( function.getNumScalarsPerPrimitive( faceID ), 1 );

         const auto memLayout = function.memoryLayout();
         auto       dofs      = function.dofMemory( faceID, level );

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

template < typename ValueType >
void P3Function< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                           uint_t                                              level,
                                           DoFType                                             flag ) const
{
   vertexDoFFunction_.interpolate( expr, level, flag );
   edgeDoFFunctionBlue_.interpolate( expr, level, real_c( 1.0 / 3.0 ), flag );
   edgeDoFFunctionRed_.interpolate( expr, level, real_c( 2.0 / 3.0 ), flag );
   faceDoFHelpers::interpolate( faceDoFFunction_, expr, level, flag );
}

template < typename ValueType >
void P3Function< ValueType >::faceDoFHelpers::interpolate( const volumedofspace::VolumeDoFFunction< ValueType >& function,
                                                           const std::function< ValueType( const Point3D& ) >&   expr,
                                                           uint_t                                                level,
                                                           DoFType                                               flag )
{
   // NOTE: We pass "flag" here, as this will become important in 3D; for 2D it is not used
   WALBERLA_UNUSED( flag );

   if ( function.getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "No 3D support in P3Function, yet!" );
   }
   else
   {
      for ( auto& it : function.getStorage()->getFaces() )
      {
         const auto  faceID = it.first;
         const auto& face   = *it.second;

         WALBERLA_CHECK_EQUAL( function.getNumScalarsPerPrimitive( faceID ), 1 );

         const auto memLayout = function.memoryLayout();
         auto       dofs      = function.dofMemory( faceID, level );

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               const Point3D centroid =
                   micromesh::microFaceCenterPosition( function.getStorage(), faceID, level, idxIt, faceType );

               const auto val = expr( Point3D( centroid( 0 ), centroid( 1 ), 0 ) );

               dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, 0, 1, level, memLayout )] = ValueType( val );
            }
         }
      }
   }
}

template < typename ValueType >
void P3Function< ValueType >::multElementwise(
    const std::vector< std::reference_wrapper< const P3Function< ValueType > > >& functions,
    uint_t                                                                        level,
    DoFType                                                                       flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > >      vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >                   edgeDoFFunctionsBlue;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >                   edgeDoFFunctionsRed;
   std::vector< std::reference_wrapper< const volumedofspace::VolumeDoFFunction< ValueType > > > faceDoFFunctions;

   for ( const P3Function< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctionsBlue.push_back( function.edgeDoFFunctionBlue_ );
      edgeDoFFunctionsRed.push_back( function.edgeDoFFunctionRed_ );
      faceDoFFunctions.push_back( function.faceDoFFunction_ );
   }

   vertexDoFFunction_.multElementwise( vertexDoFFunctions, level, flag );
   edgeDoFFunctionBlue_.multElementwise( edgeDoFFunctionsBlue, level, flag );
   edgeDoFFunctionRed_.multElementwise( edgeDoFFunctionsRed, level, flag );

   WALBERLA_ABORT( "Need to implement multelementwise for facedof part!" );
   // faceDoFFunction_.multElementwise( faceDoFFunctions, level, flag );
}

template < typename ValueType >
void P3Function< ValueType >::add( ValueType scalar, uint_t level, DoFType flag ) const
{
   vertexDoFFunction_.add( scalar, level, flag );
   edgeDoFFunctionBlue_.add( scalar, level, flag );
   edgeDoFFunctionRed_.add( scalar, level, flag );
   faceDoFFunction_.add( scalar, level, flag );
}

template < typename ValueType >
void P3Function< ValueType >::add( const std::vector< ValueType >&                                               scalars,
                                   const std::vector< std::reference_wrapper< const P3Function< ValueType > > >& functions,
                                   uint_t                                                                        level,
                                   DoFType                                                                       flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > >      vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >                   edgeDoFFunctionsBlue;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >                   edgeDoFFunctionsRed;
   std::vector< std::reference_wrapper< const volumedofspace::VolumeDoFFunction< ValueType > > > faceDoFFunctions;

   for ( const P3Function< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctionsBlue.push_back( function.edgeDoFFunctionBlue_ );
      edgeDoFFunctionsRed.push_back( function.edgeDoFFunctionRed_ );
      faceDoFFunctions.push_back( function.faceDoFFunction_ );
   }

   vertexDoFFunction_.add( scalars, vertexDoFFunctions, level, flag );
   edgeDoFFunctionBlue_.add( scalars, edgeDoFFunctionsBlue, level, flag );
   edgeDoFFunctionRed_.add( scalars, edgeDoFFunctionsRed, level, flag );

   // Will need a version that supports flag argument for 3D!
   faceDoFFunction_.add( scalars, faceDoFFunctions, level );
}

template < typename ValueType >
void P3Function< ValueType >::swap( const P3Function< ValueType >& other, const uint_t& level, const DoFType& flag ) const
{
   vertexDoFFunction_.swap( other.getVertexDoFFunction(), level, flag );
   edgeDoFFunctionBlue_.swap( other.getEdgeDoFFunctionBlue(), level, flag );
   edgeDoFFunctionRed_.swap( other.getEdgeDoFFunctionRed(), level, flag );
   faceDoFFunction_.swap( other.getFaceDoFFunction(), level );
}

template < typename ValueType >
void P3Function< ValueType >::assign( const std::vector< ValueType >&                                               scalars,
                                      const std::vector< std::reference_wrapper< const P3Function< ValueType > > >& functions,
                                      uint_t                                                                        level,
                                      DoFType                                                                       flag ) const
{
   std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > >      vertexDoFFunctions;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >                   edgeDoFFunctionsBlue;
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >                   edgeDoFFunctionsRed;
   std::vector< std::reference_wrapper< const volumedofspace::VolumeDoFFunction< ValueType > > > faceDoFFunctions;

   for ( const P3Function< ValueType >& function : functions )
   {
      vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
      edgeDoFFunctionsBlue.push_back( function.edgeDoFFunctionBlue_ );
      edgeDoFFunctionsRed.push_back( function.edgeDoFFunctionRed_ );
      faceDoFFunctions.push_back( function.faceDoFFunction_ );
   }

   vertexDoFFunction_.assign( scalars, vertexDoFFunctions, level, flag );
   edgeDoFFunctionBlue_.assign( scalars, edgeDoFFunctionsBlue, level, flag );
   edgeDoFFunctionRed_.assign( scalars, edgeDoFFunctionsRed, level, flag );

   // Will need a version that supports flag argument for 3D!
   faceDoFFunction_.assign( scalars, faceDoFFunctions, level );
}

template < typename ValueType >
ValueType P3Function< ValueType >::dotLocal( const P3Function< ValueType >& rhs, const uint_t level, const DoFType& flag ) const
{
   auto sum = ValueType( 0 );
   sum += vertexDoFFunction_.dotLocal( rhs.vertexDoFFunction_, level, flag );
   sum += edgeDoFFunctionBlue_.dotLocal( rhs.edgeDoFFunctionBlue_, level, flag );
   sum += edgeDoFFunctionRed_.dotLocal( rhs.edgeDoFFunctionRed_, level, flag );

   // Will need a version that supports flag argument for 3D!
   sum += faceDoFFunction_.dotLocal( rhs.faceDoFFunction_, level );

   return sum;
}

template < typename ValueType >
ValueType P3Function< ValueType >::dotGlobal( const P3Function< ValueType >& rhs, const uint_t level, const DoFType& flag ) const
{
   ValueType sum = dotLocal( rhs, level, flag );
   this->startTiming( "Dot (reduce)" );
   walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Dot (reduce)" );
   return sum;
}

template < typename ValueType >
BoundaryCondition P3Function< ValueType >::getBoundaryCondition() const
{
   WALBERLA_ASSERT_EQUAL( vertexDoFFunction_.getBoundaryCondition(),
                          edgeDoFFunctionBlue_.getBoundaryCondition(),
                          "P3Function: boundary conditions of underlying component functions differ!" )
   WALBERLA_ASSERT_EQUAL( vertexDoFFunction_.getBoundaryCondition(),
                          edgeDoFFunctionRed_.getBoundaryCondition(),
                          "P3Function: boundary conditions of underlying component functions differ!" )

   // NOTE: We will need faceDoFFunction_.getBoundaryCondition() for 3D!

   // WALBERLA_ASSERT_EQUAL( vertexDoFFunction_.getBoundaryCondition(),
   //                        faceDoFFunction_.getBoundaryCondition(),
   //                        "P3Function: boundary conditions of underlying component functions differ!" )

   return vertexDoFFunction_.getBoundaryCondition();
}

template < typename ValueType >
void P3Function< ValueType >::setBoundaryCondition( BoundaryCondition bc )
{
   vertexDoFFunction_.setBoundaryCondition( bc );
   edgeDoFFunctionBlue_.setBoundaryCondition( bc );
   edgeDoFFunctionRed_.setBoundaryCondition( bc );

   // NOTE: We will need faceDoFFunction_.setBoundaryCondition() for 3D!
   // faceDoFFunction_.setBoundaryCondition( bc );
}

template < typename ValueType >
void P3Function< ValueType >::copyFrom( const P3Function< ValueType >&         other,
                                        const uint_t&                          level,
                                        const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                                        const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const
{
   vertexDoFFunction_.copyFrom( other.getVertexDoFFunction(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
   edgeDoFFunctionBlue_.copyFrom( other.getEdgeDoFFunctionBlue(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
   edgeDoFFunctionRed_.copyFrom( other.getEdgeDoFFunctionRed(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );

   WALBERLA_ABORT( "Need to implement copyFrom for facedof part!" );
   // faceDoFFunction_.copyFrom( other.getFaceDoFFunction(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
}

template < typename ValueType >
void P3Function< ValueType >::toVector( const P3Function< idx_t >&            numerator,
                                        const std::shared_ptr< VectorProxy >& vec,
                                        uint_t                                level,
                                        DoFType                               flag ) const
{
   vertexDoFFunction_.toVector( numerator.getVertexDoFFunction(), vec, level, flag );
   edgeDoFFunctionBlue_.toVector( numerator.getEdgeDoFFunctionBlue(), vec, level, flag );
   edgeDoFFunctionRed_.toVector( numerator.getEdgeDoFFunctionRed(), vec, level, flag );
   faceDoFFunction_.toVector( numerator.getFaceDoFFunction(), vec, level, flag );
}

template < typename ValueType >
void P3Function< ValueType >::fromVector( const P3Function< idx_t >&            numerator,
                                          const std::shared_ptr< VectorProxy >& vec,
                                          uint_t                                level,
                                          DoFType                               flag ) const
{
   vertexDoFFunction_.fromVector( numerator.getVertexDoFFunction(), vec, level, flag );
   edgeDoFFunctionBlue_.fromVector( numerator.getEdgeDoFFunctionBlue(), vec, level, flag );
   edgeDoFFunctionRed_.fromVector( numerator.getEdgeDoFFunctionRed(), vec, level, flag );
   faceDoFFunction_.fromVector( numerator.getFaceDoFFunction(), vec, level, flag );
};

template < typename ValueType >
ValueType P3Function< ValueType >::getMaxDoFValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMax = -std::numeric_limits< ValueType >::max();
   localMax      = std::max( localMax, vertexDoFFunction_.getMaxDoFValue( level, flag, false ) );
   localMax      = std::max( localMax, edgeDoFFunctionBlue_.getMaxDoFValue( level, flag, false ) );
   localMax      = std::max( localMax, edgeDoFFunctionRed_.getMaxDoFValue( level, flag, false ) );
   // will need to pass a flag below for 3D!
   localMax = std::max( localMax, faceDoFFunction_.getMaxDoFValue( level, false ) );
   walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType P3Function< ValueType >::getMaxDoFMagnitude( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMax = ValueType( 0.0 );
   localMax      = std::max( localMax, vertexDoFFunction_.getMaxDoFMagnitude( level, flag, false ) );
   localMax      = std::max( localMax, edgeDoFFunctionBlue_.getMaxDoFMagnitude( level, flag, false ) );
   localMax      = std::max( localMax, edgeDoFFunctionRed_.getMaxDoFMagnitude( level, flag, false ) );
   // will need to pass a flag below for 3D!
   localMax = std::max( localMax, faceDoFFunction_.getMaxDoFMagnitude( level, false ) );
   walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMax = localMax;
   if ( mpiReduce )
   {
      globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
   }

   return globalMax;
}

template < typename ValueType >
ValueType P3Function< ValueType >::getMinDoFValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   auto localMin = std::numeric_limits< ValueType >::max();
   localMin      = std::min( localMin, vertexDoFFunction_.getMinDoFValue( level, flag, false ) );
   localMin      = std::min( localMin, edgeDoFFunctionBlue_.getMinDoFValue( level, flag, false ) );
   localMin      = std::min( localMin, edgeDoFFunctionRed_.getMinDoFValue( level, flag, false ) );
   // will need to pass a flag below for 3D!
   localMin = std::min( localMin, faceDoFFunction_.getMinDoFValue( level, false ) );
   walberla::mpi::allReduceInplace( localMin, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() );

   ValueType globalMin = localMin;
   if ( mpiReduce )
   {
      globalMin = -walberla::mpi::allReduce( -localMin, walberla::mpi::MAX );
   }

   return globalMin;
}

// ========================
//  explicit instantiation
// ========================
template class P3Function< double >;
template class P3Function< float >;
template class P3Function< int32_t >;
template class P3Function< int64_t >;

} //namespace hyteg
