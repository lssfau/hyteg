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
   edgeDoFFunctionBlue_.interpolate( expr, level, flag ); // WRONG !!! need to interpolate at real_c( 1.0 / 3.0 ) not 0.5
   edgeDoFFunctionRed_.interpolate( expr, level, flag );  // WRONG !!! need to interpolate at real_c( 2.0 / 3.0 ) not 0.5
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

// ========================
//  explicit instantiation
// ========================
template class P3Function< double >;
template class P3Function< float >;
template class P3Function< int32_t >;
template class P3Function< int64_t >;

} //namespace hyteg
