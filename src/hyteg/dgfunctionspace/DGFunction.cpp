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
, volumeDoFFunction_( name,
                      storage,
                      minLevel,
                      maxLevel,
                      basis->numDoFsPerElement( initialPolyDegree ),
                      volumedofspace::VolumeDoFMemoryLayout::SoA )
{}

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

         if ( isPointInTriangle( coordinates2D, faceCoodinates0, faceCoodinates1, faceCoodinates2 ) )
         {
            indexing::Index   elementIndex;
            facedof::FaceType faceType;
            Point2D           localCoordinates;

            volumedofspace::getLocalElementFromCoordinates< ValueType >(
                level, face, coordinates2D, elementIndex, faceType, localCoordinates );

            Eigen::Matrix< real_t, 2, 1 > refPos( localCoordinates[0], localCoordinates[1] );

            std::vector< real_t > dofs( ndofs );
            for ( uint_t i = 0; i < ndofs; i++ )
            {
               dofs[i] = real_t( volumeDoFFunction_.dof( faceID, elementIndex, i, faceType, level ) );
            }

            real_t value_r;
            basis_->evaluate( polyDegree, refPos, dofs, value_r );

            value = ValueType( value_r );
            return true;
         }
      }

      if ( searchToleranceRadius > 0 )
      {
         WALBERLA_ABORT( "not implemented" );
         for ( auto& it : this->getStorage()->getFaces() )
         {
            Face& face = *it.second;

            Point2D faceCoodinates0( { face.getCoordinates()[0][0], face.getCoordinates()[0][1] } );
            Point2D faceCoodinates1( { face.getCoordinates()[1][0], face.getCoordinates()[1][1] } );
            Point2D faceCoodinates2( { face.getCoordinates()[2][0], face.getCoordinates()[2][1] } );

            if ( circleTriangleIntersection(
                     coordinates2D, searchToleranceRadius, faceCoodinates0, faceCoodinates1, faceCoodinates2 ) )
            {
               // value = vertexdof::macroface::evaluate< real_t >( level, face, coordinates, faceDataID_ );
               return true;
            }
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

/// explicit instantiation
template class DGFunction< double >;
template class DGFunction< float >;
template class DGFunction< int32_t >;
template class DGFunction< int64_t >;

} // namespace dg
} // namespace hyteg