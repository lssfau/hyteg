/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "N1E1VectorFunction.hpp"

#include "hyteg/geometry/BlendingHelpers.hpp"
#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/n1e1functionspace/N1E1AdditivePackInfo.hpp"
#include "hyteg/n1e1functionspace/N1E1PackInfo.hpp"

#include "N1E1MacroCell.hpp"
#include "N1E1MacroEdge.hpp"
#include "N1E1MacroFace.hpp"

namespace hyteg {
namespace n1e1 {

using walberla::int_c;

template < typename ValueType >
N1E1VectorFunction< ValueType >::N1E1VectorFunction( const std::string&                         name,
                                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                                     const uint_t&                              minLevel,
                                                     const uint_t&                              maxLevel )
: N1E1VectorFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
{}

template < typename ValueType >
N1E1VectorFunction< ValueType >::N1E1VectorFunction( const std::string&                         name,
                                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                                     const uint_t&                              minLevel,
                                                     const uint_t&                              maxLevel,
                                                     const BoundaryCondition&                   boundaryCondition )
: VectorFunction< N1E1VectorFunction< ValueType > >( name, storage, minLevel, maxLevel )
, storage_( storage )
, dofs_( std::make_shared< EdgeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel, boundaryCondition ) )
, boundaryCondition_( boundaryCondition )
{
   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      communicators_[level]->addPackInfo( std::make_shared< N1E1PackInfo< ValueType > >( level,
                                                                                         dofs_->getVertexDataID(),
                                                                                         dofs_->getEdgeDataID(),
                                                                                         dofs_->getFaceDataID(),
                                                                                         dofs_->getCellDataID(),
                                                                                         this->getStorage() ) );
      additiveCommunicators_[level]->addPackInfo( std::make_shared< N1E1AdditivePackInfo< ValueType > >( level,
                                                                                                         dofs_->getVertexDataID(),
                                                                                                         dofs_->getEdgeDataID(),
                                                                                                         dofs_->getFaceDataID(),
                                                                                                         dofs_->getCellDataID(),
                                                                                                         this->getStorage() ) );
   }
}

template < typename ValueType >
bool N1E1VectorFunction< ValueType >::evaluate( const Point3D& physicalCoords,
                                                uint_t         level,
                                                VectorType&    value,
                                                real_t         searchToleranceRadius ) const
{
   WALBERLA_UNUSED( physicalCoords );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( value );
   WALBERLA_UNUSED( searchToleranceRadius );

   WALBERLA_ABORT( "N1E1VectorFunction< ValueType >::evaluate not implemented for requested template parameter" );
}

template <>
bool N1E1VectorFunction< real_t >::evaluate( const Point3D& physicalCoords,
                                             uint_t         level,
                                             VectorType&    value,
                                             real_t         searchToleranceRadius ) const
{
   auto [found, cellID, computationalCoords] =
       mapFromPhysicalToComputationalDomain3D( this->getStorage(), physicalCoords, searchToleranceRadius );

   if ( found )
   {
      value = n1e1::macrocell::evaluate(
          level, *( this->getStorage()->getCell( cellID ) ), computationalCoords, dofs_->getCellDataID() );
   }

   return found;
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::evaluateOnMicroElement( const Point3D&         coordinates,
                                                              uint_t                 level,
                                                              const PrimitiveID&     cellID,
                                                              hyteg::indexing::Index elementIndex,
                                                              celldof::CellType      cellType,
                                                              VectorType&            value ) const
{
   WALBERLA_UNUSED( coordinates );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( cellID );
   WALBERLA_UNUSED( elementIndex );
   WALBERLA_UNUSED( cellType );
   WALBERLA_UNUSED( value );

   WALBERLA_ABORT( "N1E1VectorFunction< ValueType >::evaluateOnMicroElement not implemented for requested template parameter" );
}

template <>
void N1E1VectorFunction< real_t >::evaluateOnMicroElement( const Point3D&               coordinates,
                                                           const uint_t                 level,
                                                           const PrimitiveID&           cellID,
                                                           const hyteg::indexing::Index elementIndex,
                                                           const celldof::CellType      cellType,
                                                           VectorType&                  value ) const
{
   const Cell& cell = *storage_->getCell( cellID );
   value = n1e1::macrocell::evaluateOnMicroElement( level, cell, elementIndex, cellType, coordinates, dofs_->getCellDataID() );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::assign(
    const std::vector< ValueType >&                                                       scalars,
    const std::vector< std::reference_wrapper< const N1E1VectorFunction< ValueType > > >& functions,
    uint_t                                                                                level,
    DoFType                                                                               flag ) const
{
   // the interpolation operator (x ↦ ∫ₑ x·t dΓ) is linear
   // ⇒ we can assign on the DoFs directly

   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > > dofFunctions;

   for ( const N1E1VectorFunction& function : functions )
   {
      dofFunctions.push_back( *function.getDoFs() );
   }

   dofs_->assign( scalars, dofFunctions, level, flag );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::multElementwise(
    const std::vector< std::reference_wrapper< const N1E1VectorFunction< ValueType > > >& functions,
    uint_t                                                                                level,
    DoFType                                                                               flag ) const
{
   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > > dofFunctions;

   for ( const N1E1VectorFunction& function : functions )
   {
      dofFunctions.push_back( *function.getDoFs() );
   }

   dofs_->multElementwise( dofFunctions, level, flag );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::add( VectorType vector, uint_t level, DoFType flag ) const
{
   WALBERLA_UNUSED( vector );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( flag );

   WALBERLA_ABORT( "N1E1VectorFunction< ValueType >::add not implemented for requested template parameter" );
}

template <>
void N1E1VectorFunction< real_t >::add( VectorType vector, uint_t level, DoFType flag ) const
{
   this->startTiming( "Add (vector)" );

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         n1e1::macroedge::add( level, edge, vector, dofs_->getEdgeDataID() );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         n1e1::macroface::add( level, face, vector, dofs_->getFaceDataID() );
      }
   }

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            n1e1::macrocell::add( level, cell, vector, dofs_->getCellDataID() );
         }
      }
   }

   this->stopTiming( "Add (vector)" );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::add(
    const std::vector< ValueType >&                                                       scalars,
    const std::vector< std::reference_wrapper< const N1E1VectorFunction< ValueType > > >& functions,
    size_t                                                                                level,
    DoFType                                                                               flag ) const
{
   // the interpolation operator (x ↦ ∫ₑ x·t dΓ) is linear
   // ⇒ we can add the DoFs directly

   std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > > dofFunctions;

   for ( const N1E1VectorFunction& function : functions )
   {
      dofFunctions.push_back( *function.getDoFs() );
   }

   dofs_->add( scalars, dofFunctions, level, flag );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::interpolate( VectorType constant, uint_t level, DoFType flag ) const
{
   WALBERLA_UNUSED( constant );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( flag );

   WALBERLA_ABORT( "N1E1VectorFunction< ValueType >::interpolate not implemented for requested template parameter" );
}

template <>
void N1E1VectorFunction< real_t >::interpolate( VectorType constant, uint_t level, DoFType flag ) const
{
   this->startTiming( "Interpolate" );

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         n1e1::macroedge::interpolate( level, edge, dofs_->getEdgeDataID(), constant );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         n1e1::macroface::interpolate( level, face, dofs_->getFaceDataID(), constant );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         n1e1::macrocell::interpolate( level, cell, dofs_->getCellDataID(), constant );
      }
   }

   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::interpolate( VectorType constant, uint_t level, BoundaryUID boundaryUID ) const
{
   WALBERLA_UNUSED( constant );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( boundaryUID );

   WALBERLA_ABORT( "N1E1VectorFunction< ValueType >::interpolate not implemented for requested template parameter" );
}

template <>
void N1E1VectorFunction< real_t >::interpolate( VectorType constant, uint_t level, BoundaryUID boundaryUID ) const
{
   this->startTiming( "Interpolate" );

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( edge.getMeshBoundaryFlag() ) == boundaryUID )
      {
         n1e1::macroedge::interpolate( level, edge, dofs_->getEdgeDataID(), constant );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( face.getMeshBoundaryFlag() ) == boundaryUID )
      {
         n1e1::macroface::interpolate( level, face, dofs_->getFaceDataID(), constant );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( cell.getMeshBoundaryFlag() ) == boundaryUID )
      {
         n1e1::macrocell::interpolate( level, cell, dofs_->getCellDataID(), constant );
      }
   }

   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::interpolate( const std::function< VectorType( const Point3D& ) >& expr,
                                                   uint_t                                               level,
                                                   DoFType                                              flag ) const
{
   std::function< VectorType( const Point3D&, const std::vector< VectorType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< VectorType >& ) { return expr( x ); };
   interpolate( exprExtended, {}, level, flag );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::interpolate( const std::function< VectorType( const Point3D& ) >& expr,
                                                   uint_t                                               level,
                                                   BoundaryUID                                          boundaryUID ) const
{
   std::function< VectorType( const Point3D&, const std::vector< VectorType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< VectorType >& ) { return expr( x ); };
   interpolate( exprExtended, {}, level, boundaryUID );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::interpolate(
    const std::function< VectorType( const Point3D&, const std::vector< VectorType >& ) >& expr,
    const std::vector< std::reference_wrapper< const N1E1VectorFunction< ValueType > > >&  srcFunctions,
    uint_t                                                                                 level,
    DoFType                                                                                flag ) const
{
   WALBERLA_UNUSED( expr );
   WALBERLA_UNUSED( srcFunctions );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( flag );

   WALBERLA_ABORT( "N1E1VectorFunction< ValueType >::interpolate not implemented for requested template parameter" );
}

template <>
void N1E1VectorFunction< real_t >::interpolate(
    const std::function< VectorType( const Point3D&, const std::vector< VectorType >& ) >& expr,
    const std::vector< std::reference_wrapper< const N1E1VectorFunction< real_t > > >&     srcFunctions,
    uint_t                                                                                 level,
    DoFType                                                                                flag ) const
{
   this->startTiming( "Interpolate" );

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         n1e1::macroedge::interpolate( level, edge, dofs_->getEdgeDataID(), srcFunctions, expr );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         n1e1::macroface::interpolate( level, face, dofs_->getFaceDataID(), srcFunctions, expr );
      }
   }

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            n1e1::macrocell::interpolate( level, cell, dofs_->getCellDataID(), srcFunctions, expr );
         }
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::interpolate(
    const std::function< VectorType( const Point3D&, const std::vector< VectorType >& ) >& expr,
    const std::vector< std::reference_wrapper< const N1E1VectorFunction< ValueType > > >&  srcFunctions,
    uint_t                                                                                 level,
    BoundaryUID                                                                            boundaryUID ) const
{
   WALBERLA_UNUSED( expr );
   WALBERLA_UNUSED( srcFunctions );
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( boundaryUID );

   WALBERLA_ABORT( "N1E1VectorFunction< ValueType >::interpolate not implemented for requested template parameter" );
}

template <>
void N1E1VectorFunction< real_t >::interpolate(
    const std::function< VectorType( const Point3D&, const std::vector< VectorType >& ) >& expr,
    const std::vector< std::reference_wrapper< const N1E1VectorFunction< real_t > > >&     srcFunctions,
    uint_t                                                                                 level,
    BoundaryUID                                                                            boundaryUID ) const
{
   this->startTiming( "Interpolate" );

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( edge.getMeshBoundaryFlag() ) == boundaryUID )
      {
         n1e1::macroedge::interpolate( level, edge, dofs_->getEdgeDataID(), srcFunctions, expr );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( face.getMeshBoundaryFlag() ) == boundaryUID )
      {
         n1e1::macroface::interpolate( level, face, dofs_->getFaceDataID(), srcFunctions, expr );
      }
   }

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

         if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( cell.getMeshBoundaryFlag() ) == boundaryUID )
         {
            n1e1::macrocell::interpolate( level, cell, dofs_->getCellDataID(), srcFunctions, expr );
         }
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void N1E1VectorFunction< ValueType >::setLocalCommunicationMode(
    const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode )
{
   for ( auto& communicator : communicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
   for ( auto& communicator : additiveCommunicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
}

// ========================
//  explicit instantiation
// ========================
template class N1E1VectorFunction< double >; // TODO real_t?
// template class N1E1VectorFunction< float >;
template class N1E1VectorFunction< int32_t >;
template class N1E1VectorFunction< int64_t >;

} // namespace n1e1
} // namespace hyteg
