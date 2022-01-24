/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/dgfunctionspace/DGFunction.hpp"

#include "hyteg/geometry/Intersection.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {
namespace dg {

template < typename ValueType >
DGFunction< ValueType >::DGFunction( const std::string&                         name,
                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                     uint_t                                     minLevel,
                                     uint_t                                     maxLevel,
                                     const std::shared_ptr< DGBasisInfo >&      basis,
                                     int                                        initialPolyDegree )
: Function< DGFunction< ValueType > >( name, storage, minLevel, maxLevel )
, name_( name )
, storage_( storage )
, minLevel_( minLevel )
, maxLevel_( maxLevel )
, basis_( basis )
{
   if ( storage->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Not implemented for cells yet." )
   }
   else
   {
      for ( auto pid : storage->getFaceIDs() )
      {
         polyDegreesPerPrimitive_[pid] = initialPolyDegree;
      }
   }

   volumeDoFFunction_ =
       std::make_shared< volumedofspace::VolumeDoFFunction< ValueType > >( name,
                                                                           storage,
                                                                           minLevel,
                                                                           maxLevel,
                                                                           basis->numDoFsPerElement( initialPolyDegree ),
                                                                           volumedofspace::VolumeDoFMemoryLayout::SoA );
}

template < typename ValueType >
bool DGFunction< ValueType >::evaluate( const Point3D& coordinates,
                                        uint_t         level,
                                        ValueType&     value,
                                        real_t         searchToleranceRadius ) const
{
   if ( !this->storage_->hasGlobalCells() )
   {
      // 2D

      Point2D coordinates2D( { coordinates[0], coordinates[1] } );

      for ( auto& it : this->getStorage()->getFaces() )
      {
         PrimitiveID faceID = it.first;
         Face&       face   = *it.second;

         const auto polyDegree = polyDegreesPerPrimitive_.at( faceID );
         const auto ndofs      = uint_c( basis_->numDoFsPerElement( polyDegree ) );

         Point2D faceCoodinates0( { face.getCoordinates()[0][0], face.getCoordinates()[0][1] } );
         Point2D faceCoodinates1( { face.getCoordinates()[1][0], face.getCoordinates()[1][1] } );
         Point2D faceCoodinates2( { face.getCoordinates()[2][0], face.getCoordinates()[2][1] } );

         if ( isPointInTriangle( coordinates2D, faceCoodinates0, faceCoodinates1, faceCoodinates2 ) ||
              ( searchToleranceRadius > 0 &&
                circleTriangleIntersection(
                    coordinates2D, searchToleranceRadius, faceCoodinates0, faceCoodinates1, faceCoodinates2 ) ) )
         {
            indexing::Index   elementIndex;
            facedof::FaceType faceType;
            Point2D           localCoordinates;

            volumedofspace::getLocalElementFromCoordinates< ValueType >(
                level, face, coordinates2D, elementIndex, faceType, localCoordinates );

            Eigen::Matrix< real_t, 2, 1 > refPos( localCoordinates[0], localCoordinates[1] );
#if 0
            std::array< Eigen::Matrix< real_t, 2, 1 >, 3 > affineElementVertices;
            auto vertexIndices = facedof::macroface::getMicroVerticesFromMicroFace( elementIndex, faceType );
            for ( uint_t i = 0; i < 3; i++ )
            {
               const auto coord              = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
               affineElementVertices[i]( 0 ) = coord[0];
               affineElementVertices[i]( 1 ) = coord[1];
            }
#endif
            std::vector< real_t > dofs( ndofs );
            for ( uint_t i = 0; i < ndofs; i++ )
            {
               dofs[i] = real_t( volumeDoFFunction_->dof( faceID, elementIndex, i, faceType, level ) );
            }

            real_t value_r;
            basis_->evaluate( polyDegree, refPos, dofs, value_r );

            // Eigen::Matrix< real_t, 2, 1 > coords( coordinates2D[0], coordinates2D[1] );
            // basis_->evaluate( polyDegree, affineElementVertices, coords, dofs, value_r );

            value = ValueType( value_r );
            return true;
         }
      }
   }
   else
   {
      WALBERLA_ABORT( "not implemented" );
#if 0
      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if ( isPointInTetrahedron( coordinates,

                                    cell.getCoordinates()[0],
                                    cell.getCoordinates()[1],
                                    cell.getCoordinates()[2],
                                    cell.getCoordinates()[3],
                                    cell.getFaceInwardNormal( 0 ),
                                    cell.getFaceInwardNormal( 1 ),
                                    cell.getFaceInwardNormal( 2 ),
                                    cell.getFaceInwardNormal( 3 ) ) )
         {
            value = vertexdof::macrocell::evaluate< real_t >( level, cell, coordinates, cellDataID_ );
            return true;
         }
      }

      if ( searchToleranceRadius > 0 )
      {
         for ( auto& it : this->getStorage()->getCells() )
         {
            Cell& cell = *it.second;

            if ( sphereTetrahedronIntersection( coordinates,
                                                searchToleranceRadius,
                                                cell.getCoordinates()[0],
                                                cell.getCoordinates()[1],
                                                cell.getCoordinates()[2],
                                                cell.getCoordinates()[3] ) )
            {
               value = vertexdof::macrocell::evaluate< real_t >( level, cell, coordinates, cellDataID_ );
               return true;
            }
         }
      }
#endif
   }

   return false;
}

template < typename ValueType >
void DGFunction< ValueType >::evaluateLinearFunctional( const std::function< real_t( const Point3D& ) >& f, uint_t level )
{
   if ( storage_->hasGlobalCells() )
   {
      WALBERLA_ABORT( "Linear functional evaluation not implemented." )
   }
   else
   {
      for ( auto& it : this->getStorage()->getFaces() )
      {
         const auto faceID = it.first;
         const auto face   = *it.second;

         const auto degree  = polyDegreesPerPrimitive_.at( faceID );
         const auto numDofs = basis()->numDoFsPerElement( degree );

         auto       dofs      = volumeDoFFunction()->dofMemory( faceID, level );
         const auto memLayout = volumeDoFFunction()->memoryLayout();

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               const std::array< indexing::Index, 3 > vertexIndices =
                   facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
               std::array< Eigen::Matrix< real_t, 2, 1 >, 3 > elementVertices;
               for ( uint_t i = 0; i < 3; i++ )
               {
                  const auto elementVertex = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                  elementVertices[i]( 0 )  = elementVertex[0];
                  elementVertices[i]( 1 )  = elementVertex[1];
               }

               std::vector< real_t > dofValues( numDofs );
               basis()->integrateBasisFunction( degree, elementVertices, f, dofValues );

               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )] =
                      ValueType( dofValues[i] );
               }
            }
         }
      }
   }
}

template < typename ValueType >
void DGFunction< ValueType >::enumerate( uint_t level ) const
{
   ValueType offset = 0;
   this->enumerate( level, offset );
}

template < typename ValueType >
void DGFunction< ValueType >::enumerate( uint_t level, ValueType& offset ) const
{
   if ( storage_->hasGlobalCells() )
   {
      // 3D
      WALBERLA_ABORT( "enumerate() not implemented in 3D." );
   }
   else
   {
      // 2D
      for ( const auto& it : storage_->getFaces() )
      {
         const auto faceID = it.first;
         const auto face   = *it.second;

         const auto degree  = polyDegreesPerPrimitive_.at( faceID );
         const auto numDofs = basis()->numDoFsPerElement( degree );

         auto       dofs      = volumeDoFFunction()->dofMemory( faceID, level );
         const auto memLayout = volumeDoFFunction()->memoryLayout();

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )] =
                      offset++;
               }
            }
         }
      }
   }
}

/// explicit instantiation
template class DGFunction< double >;
template class DGFunction< float >;
template class DGFunction< int32_t >;
template class DGFunction< int64_t >;

} // namespace dg

void createVectorFromFunction( const dg::DGFunction< real_t >&       function,
                               const dg::DGFunction< idx_t >&        numerator,
                               const std::shared_ptr< VectorProxy >& vec,
                               uint_t                                level,
                               DoFType                               flag )
{
   if ( function.getStorage()->hasGlobalCells() )
   {
      // 3D
      WALBERLA_ABORT( "enumerate() not implemented in 3D." );
   }
   else
   {
      // 2D
      for ( const auto& it : function.getStorage()->getFaces() )
      {
         const auto faceID = it.first;
         const auto face   = *it.second;

         const auto degree  = function.polynomialDegree( faceID );
         const auto numDofs = function.basis()->numDoFsPerElement( degree );

         const auto indices   = numerator.volumeDoFFunction()->dofMemory( faceID, level );
         const auto dofs      = function.volumeDoFFunction()->dofMemory( faceID, level );
         const auto memLayout = function.volumeDoFFunction()->memoryLayout();

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  const auto vectorIdx =
                      indices[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )];
                  const auto value =
                      dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )];
                  vec->setValue( vectorIdx, value );
               }
            }
         }
      }
   }
}

void createFunctionFromVector( const dg::DGFunction< real_t >&       function,
                               const dg::DGFunction< idx_t >&        numerator,
                               const std::shared_ptr< VectorProxy >& vec,
                               uint_t                                level,
                               DoFType                               flag )
{
   if ( function.getStorage()->hasGlobalCells() )
   {
      // 3D
      WALBERLA_ABORT( "enumerate() not implemented in 3D." );
   }
   else
   {
      // 2D
      for ( const auto& it : function.getStorage()->getFaces() )
      {
         const auto faceID = it.first;
         const auto face   = *it.second;

         const auto degree  = function.polynomialDegree( faceID );
         const auto numDofs = function.basis()->numDoFsPerElement( degree );

         const auto indices   = numerator.volumeDoFFunction()->dofMemory( faceID, level );
         auto       dofs      = function.volumeDoFFunction()->dofMemory( faceID, level );
         const auto memLayout = function.volumeDoFFunction()->memoryLayout();

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               for ( uint_t i = 0; i < numDofs; i++ )
               {
                  const auto vectorIdx =
                      indices[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )];
                  const auto value = vec->getValue( vectorIdx );
                  dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), faceType, i, numDofs, level, memLayout )] = value;
               }
            }
         }
      }
   }
}

void applyDirichletBC( const dg::DGFunction< idx_t >& numerator, std::vector< idx_t >& mat, uint_t level )
{
   WALBERLA_LOG_WARNING_ON_ROOT( "DGFunction: BCs are not applied to sparse matrix." );
}

} // namespace hyteg