/*
 * Copyright (c) 2022-2023 Daniel Bauer.
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

#include <initializer_list>
#include <vector>

#include "core/DataTypes.h"
#include "core/debug/Debug.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFOrientation.hpp"
#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/LocalIDMappings.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/n1e1functionspace/N1E1Indexing.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/primitives/Cell.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {
namespace n1e1 {
namespace macrocell {

using indexing::Index;
using walberla::int_c;
using walberla::uint_t;
template < typename ValueType >
using VectorType = typename N1E1VectorFunction< ValueType >::VectorType;

/// Returns the non-unit tangent vector of a micro-edge, where the length of the
/// vector equals the length of the edge.
inline Point3D microEdgeDirection( const uint_t& level, const Cell& cell, const edgedof::EdgeDoFOrientation& orientation )
{
   const real_t  stepFrequency = real_c( 1.0 ) / real_c( levelinfo::num_microedges_per_edge( level ) );
   const Point3D xDir          = ( cell.getCoordinates()[1] - cell.getCoordinates()[0] ) * stepFrequency;
   const Point3D yDir          = ( cell.getCoordinates()[2] - cell.getCoordinates()[0] ) * stepFrequency;
   const Point3D zDir          = ( cell.getCoordinates()[3] - cell.getCoordinates()[0] ) * stepFrequency;

   switch ( orientation )
   {
   case edgedof::EdgeDoFOrientation::X:
      return xDir;
   case edgedof::EdgeDoFOrientation::Y:
      return yDir;
   case edgedof::EdgeDoFOrientation::Z:
      return zDir;
   case edgedof::EdgeDoFOrientation::XY:
      return yDir - xDir;
   case edgedof::EdgeDoFOrientation::XZ:
      return zDir - xDir;
   case edgedof::EdgeDoFOrientation::YZ:
      return zDir - yDir;
   case edgedof::EdgeDoFOrientation::XYZ:
      return zDir - yDir + xDir;
   default:
      WALBERLA_ABORT( "wrong orienation" )
   }
}

/// Evaluates the element local function at the specified coordinates.
/// If `xComp` is outside the specified element, the function is extrapolated.
/// \param xComp Coordinates in computational coordinate system.
/// \returns The function value transformed to the physical coordinate system.
inline VectorType< real_t > evaluateOnMicroElement( const uint_t&                                            level,
                                                    const Cell&                                              cell,
                                                    const indexing::Index                                    elementIndex,
                                                    const celldof::CellType                                  cellType,
                                                    const Point3D&                                           xComp,
                                                    const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dataID )
{
   using ValueType = real_t;

   // get micro vertex and edge dof indices of microcell

   std::array< uint_t, 6 >          edgeDoFIndices;
   std::array< indexing::Index, 4 > microCellIndices = getMicroVerticesFromMicroCell( elementIndex, cellType );
   edgedof::getEdgeDoFDataIndicesFromMicroVerticesFEniCSOrdering( microCellIndices, level, edgeDoFIndices );

   // get local coordinates

   auto microTet0 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[0] );
   auto microTet1 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[1] );
   auto microTet2 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[2] );
   auto microTet3 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[3] );

   auto xLocal = vertexdof::macrocell::detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, xComp );

   auto x = xLocal[0];
   auto y = xLocal[1];
   auto z = xLocal[2];

   // get DoFs

   auto edgedofData = cell.getData( dataID )->getPointer( level );

   auto valueTetE0 = edgedofData[edgeDoFIndices[0]];
   auto valueTetE1 = edgedofData[edgeDoFIndices[1]];
   auto valueTetE2 = edgedofData[edgeDoFIndices[2]];
   auto valueTetE3 = edgedofData[edgeDoFIndices[3]];
   auto valueTetE4 = edgedofData[edgeDoFIndices[4]];
   auto valueTetE5 = edgedofData[edgeDoFIndices[5]];

   // evaluate function in reference tet

   // basis functions N1E1 φ(x, y, z):
   //   at [0.  0.5 0.5]: (     0,   -z  ,    y  )ᵀ
   //   at [0.5 0.  0.5]: (  -z  ,      0,  x    )ᵀ
   //   at [0.5 0.5 0. ]: (-y    ,  x    ,      0)ᵀ
   //   at [0.  0.  0.5]: (   z  ,    z  , -x-y+1)ᵀ
   //   at [0.  0.5 0. ]: ( y    , -x-z+1,    y  )ᵀ
   //   at [0.5 0.  0. ]: (-y-z+1,  x    ,  x    )ᵀ

   VectorType< ValueType > scaleE0 = valueTetE0 * VectorType< ValueType >{ 0, -z, y };
   VectorType< ValueType > scaleE1 = valueTetE1 * VectorType< ValueType >{ -z, 0, x };
   VectorType< ValueType > scaleE2 = valueTetE2 * VectorType< ValueType >{ -y, x, 0 };
   VectorType< ValueType > scaleE3 = valueTetE3 * VectorType< ValueType >{ z, z, -x - y + 1 };
   VectorType< ValueType > scaleE4 = valueTetE4 * VectorType< ValueType >{ y, -x - z + 1, y };
   VectorType< ValueType > scaleE5 = valueTetE5 * VectorType< ValueType >{ -y - z + 1, x, x };

   VectorType< ValueType > localValue = scaleE0 + scaleE1 + scaleE2 + scaleE3 + scaleE4 + scaleE5;

   // transform to affine space (covariant Piola mapping)

   Matrix3r DF;
   cell.getGeometryMap()->evalDF( xComp, DF );

   // TODO precompute and store foctorized A (for each cell type), use also to find xLocal
   Matrix3r A;
   A.row( 0 ) = microTet1 - microTet0;
   A.row( 1 ) = microTet2 - microTet0;
   A.row( 2 ) = microTet3 - microTet0;

   return ( A * DF.transpose() ).fullPivLu().solve( localValue );
}

/// \param xComp Coordinates in computational coordinate system.
/// \returns The function value transformed to the physical coordinate system.
inline VectorType< real_t > evaluate( const uint_t&                                            level,
                                      const Cell&                                              cell,
                                      const Point3D&                                           xComp,
                                      const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dataID )
{
   using ValueType = real_t;

   // find microcell which contains `xComp`
   // note that local coordinates are determined with respect to a different vertex ordering
   // and are therefore wrong!

   indexing::Index   elementIndex;
   celldof::CellType cellType;
   Point3D           wrongLocalCoordinates;

   volumedofspace::getLocalElementFromCoordinates< ValueType >(
       level, cell, xComp, elementIndex, cellType, wrongLocalCoordinates );

   return evaluateOnMicroElement( level, cell, elementIndex, cellType, xComp, dataID );
}

inline void add( const uint_t&                                            level,
                 Cell&                                                    cell,
                 const VectorType< real_t >&                              vector,
                 const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId )
{
   using ValueType = real_t;

   // x ↦ ∫ₑ x·t dΓ, direction = tangent·length
   const ValueType dofScalarX   = vector.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::X ) );
   const ValueType dofScalarY   = vector.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::Y ) );
   const ValueType dofScalarZ   = vector.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::Z ) );
   const ValueType dofScalarXY  = vector.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::XY ) );
   const ValueType dofScalarXZ  = vector.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::XZ ) );
   const ValueType dofScalarYZ  = vector.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::YZ ) );
   const ValueType dofScalarXYZ = vector.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::XYZ ) );

   auto dstData = cell.getData( dstId )->getPointer( level );

   for ( const auto& it : edgedof::macrocell::Iterator( level ) )
   {
      const uint_t idxX  = edgedof::macrocell::xIndex( level, it.x(), it.y(), it.z() );
      const uint_t idxY  = edgedof::macrocell::yIndex( level, it.x(), it.y(), it.z() );
      const uint_t idxZ  = edgedof::macrocell::zIndex( level, it.x(), it.y(), it.z() );
      const uint_t idxXY = edgedof::macrocell::xyIndex( level, it.x(), it.y(), it.z() );
      const uint_t idxXZ = edgedof::macrocell::xzIndex( level, it.x(), it.y(), it.z() );
      const uint_t idxYZ = edgedof::macrocell::yzIndex( level, it.x(), it.y(), it.z() );

      dstData[idxX] += dofScalarX;
      dstData[idxY] += dofScalarY;
      dstData[idxZ] += dofScalarZ;
      dstData[idxXY] += dofScalarXY;
      dstData[idxXZ] += dofScalarXZ;
      dstData[idxYZ] += dofScalarYZ;
   }

   for ( const auto& it : edgedof::macrocell::IteratorXYZ( level ) )
   {
      const uint_t idxXYZ = edgedof::macrocell::xyzIndex( level, it.x(), it.y(), it.z() );
      dstData[idxXYZ] += dofScalarXYZ;
   }
}

inline void
    interpolate( const uint_t&                                                                      level,
                 Cell&                                                                              cell,
                 const PrimitiveDataID< FunctionMemory< real_t >, Cell >&                           cellMemoryId,
                 const std::vector< std::reference_wrapper< const N1E1VectorFunction< real_t > > >& srcFunctions,
                 const std::function< VectorType< real_t >( const Point3D&, const std::vector< VectorType< real_t > >& ) >& expr )
{
   using ValueType = real_t;

   auto                                   cellData = cell.getData( cellMemoryId )->getPointer( level );
   std::vector< VectorType< ValueType > > srcVector( srcFunctions.size() );

   for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) )
   {
      const Point3D                  microVertexPosition = vertexdof::macrocell::coordinateFromIndex( level, cell, it );
      const std::array< Point3D, 6 > dofShiftsFromVertex{
          edgedof::macrocell::xShiftFromVertex( level, cell ),
          edgedof::macrocell::yShiftFromVertex( level, cell ),
          edgedof::macrocell::zShiftFromVertex( level, cell ),
          edgedof::macrocell::xShiftFromVertex( level, cell ) + edgedof::macrocell::yShiftFromVertex( level, cell ),
          edgedof::macrocell::xShiftFromVertex( level, cell ) + edgedof::macrocell::zShiftFromVertex( level, cell ),
          edgedof::macrocell::yShiftFromVertex( level, cell ) + edgedof::macrocell::zShiftFromVertex( level, cell ),
      };

      for ( uint_t dofIdx = 0; dofIdx < dofShiftsFromVertex.size(); ++dofIdx )
      {
         const edgedof::EdgeDoFOrientation edgeOrientation = edgedof::allEdgeDoFOrientationsWithoutXYZ[dofIdx];
         const Point3D                     dofCoordsComp   = microVertexPosition + dofShiftsFromVertex[dofIdx];
         Point3D                           dofCoords;
         cell.getGeometryMap()->evalF( dofCoordsComp, dofCoords );

         for ( uint_t k = 0; k < srcFunctions.size(); ++k )
         {
            srcFunctions[k].get().evaluate( dofCoords, level, srcVector[k] );
         }

         Matrix3r DF;
         cell.getGeometryMap()->evalDF( dofCoordsComp, DF );

         // x ↦ ∫ₑ x·t dΓ, direction = tangent·length
         const ValueType dofValue =
             ( DF.transpose() * expr( dofCoords, srcVector ) ).dot( microEdgeDirection( level, cell, edgeOrientation ) );
         cellData[edgedof::macrocell::index( level, it.x(), it.y(), it.z(), edgeOrientation )] = dofValue;
      }
   }

   for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) )
   {
      const Point3D microVertexPosition  = vertexdof::macrocell::coordinateFromIndex( level, cell, it );
      const Point3D xyzMicroEdgePosition = microVertexPosition + edgedof::macrocell::xShiftFromVertex( level, cell ) +
                                           edgedof::macrocell::yShiftFromVertex( level, cell ) +
                                           edgedof::macrocell::zShiftFromVertex( level, cell );

      Point3D xyzBlend;
      cell.getGeometryMap()->evalF( xyzMicroEdgePosition, xyzBlend );

      for ( uint_t k = 0; k < srcFunctions.size(); ++k )
      {
         srcFunctions[k].get().evaluate( xyzBlend, level, srcVector[k] );
      }

      Matrix3r xyzDF;
      cell.getGeometryMap()->evalDF( xyzMicroEdgePosition, xyzDF );

      const ValueType dofScalarXYZ = ( xyzDF.transpose() * expr( xyzBlend, srcVector ) )
                                         .dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::XYZ ) );

      cellData[edgedof::macrocell::xyzIndex( level, it.x(), it.y(), it.z() )] = dofScalarXYZ;
   }
}

inline void interpolate( const uint_t&                                            level,
                         Cell&                                                    cell,
                         const PrimitiveDataID< FunctionMemory< real_t >, Cell >& cellMemoryId,
                         const VectorType< real_t >&                              constant )
{
   using ValueType = real_t;

   const auto geometryMap = cell.getGeometryMap();
   if ( !geometryMap->isAffine() )
   {
      // If the blending map is not affine, the vector field is not constant in computational space.
      // In that case, we delegate to the non-constant interpolation routine.
      interpolate(
          level, cell, cellMemoryId, {}, [&]( const Point3D&, const std::vector< VectorType< real_t > >& ) { return constant; } );
      return;
   }

   // Geometry map is affine ⇒ DF is constant
   Matrix3r DF;
   geometryMap->evalDF( {}, DF );
   const VectorType< real_t > valComp = DF.transpose() * constant;

   auto cellData = cell.getData( cellMemoryId )->getPointer( level );

   // x ↦ ∫ₑ x·t dΓ, direction = tangent·length
   const ValueType dofScalarX   = valComp.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::X ) );
   const ValueType dofScalarY   = valComp.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::Y ) );
   const ValueType dofScalarZ   = valComp.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::Z ) );
   const ValueType dofScalarXY  = valComp.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::XY ) );
   const ValueType dofScalarXZ  = valComp.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::XZ ) );
   const ValueType dofScalarYZ  = valComp.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::YZ ) );
   const ValueType dofScalarXYZ = valComp.dot( microEdgeDirection( level, cell, edgedof::EdgeDoFOrientation::XYZ ) );

   for ( const auto& it : edgedof::macrocell::Iterator( level, 0 ) )
   {
      cellData[edgedof::macrocell::xIndex( level, it.x(), it.y(), it.z() )]  = dofScalarX;
      cellData[edgedof::macrocell::yIndex( level, it.x(), it.y(), it.z() )]  = dofScalarY;
      cellData[edgedof::macrocell::zIndex( level, it.x(), it.y(), it.z() )]  = dofScalarZ;
      cellData[edgedof::macrocell::xyIndex( level, it.x(), it.y(), it.z() )] = dofScalarXY;
      cellData[edgedof::macrocell::xzIndex( level, it.x(), it.y(), it.z() )] = dofScalarXZ;
      cellData[edgedof::macrocell::yzIndex( level, it.x(), it.y(), it.z() )] = dofScalarYZ;
   }

   for ( const auto& it : edgedof::macrocell::IteratorXYZ( level, 0 ) )
   {
      cellData[edgedof::macrocell::xyzIndex( level, it.x(), it.y(), it.z() )] = dofScalarXYZ;
   }
}

/// Transformation from/to the local edge orientations to/from the edge
/// orientations of the owning primitive.
///
/// The returned matrix takes care of handling mismatching edge orientations. It
/// transforms DoFs from owning primitives into the basis defined by the given
/// cell, or vice versa. It is always diagonal and all entries have magnitude 1.
/// For elements inside the macro-cell it is just the identity.
///
/// In general FEM spaces, neighboring elements must agree on the orientation
/// of shared mesh entities. In n1e1, only the orientations of shared edges
/// are important. They define a local basis for FEM functions. Basis
/// transformations between these spaces allow us to view DoFs from neighboring
/// elements in our local orientation.
/// Recommended reading:
///   Scroggs et al., "Construction of Arbitrary Order Finite Element Degree-
///   of-Freedom Maps on Polygonal and Polyhedral Cell Meshes," 2022, doi:
///   https://doi.org/10.1145/3524456.
///
/// Local operations (e.g. matrix-free operator applications) can work
/// in their local basis while the communication is responsible for the
/// basis transformations. When assembling operators into matrices, these
/// transformations must be "baked into" the matrix, since vectors are assembled
/// locally and our communication routine is not performed during the operator
/// application.
inline Eigen::DiagonalMatrix< real_t, 6 >
    basisTransformation( const uint_t level, const Cell& cell, const Index& microCell, const celldof::CellType cellType )
{
   Eigen::DiagonalMatrix< real_t, 6 > transform;
   transform.setIdentity();

   // Sets of macro primitives owning DoFs shared by the given micro cell.
   std::set< uint_t > owningEdges; // local indices ∈ {0, …, 5}
   std::set< uint_t > owningFaces; // local indices ∈ {0, …, 3}

   // An invalid edge id.
   const Eigen::Index X = 42;

   // The index of the DoF (in the local element matrix/vector) with the same
   // direction as the macro edge with the given local index (HyTeG ordering).
   // Note that every element type misses one of the seven edge types (marked
   // by X).
   //
   // Example:
   //
   //   3                    3                         .3
   //   |\`\.                |\`\.                  ,/' |\
   //   | 5 `\.              | 0 `\.              .1    | 0
   //   |  \   4             |  \   1          ,/'      3  \
   //   3  2 _  `\.     -->  3  2 _  `\.      1------2--|--2
   //   |  /  `-2 `\.        |  /  `-2 `\.      `-5     | /
   //   | 1      `\_`\       | 4      `\_`\        `\_  |4
   //   0------0------1      0------5------1          `-0
   //
   //   Macro-cell           Micro-cell
   //   always WHITE UP      e.g. WHITE UP         BLUE UP
   //   HyTeG edge indices   FEniCS edge indices / DoF indices
   static const std::map< celldof::CellType, std::array< Eigen::Index, 6 > > dofIdx{
       { celldof::CellType::WHITE_UP, { 5, 4, 2, 3, 1, 0 } },
       { celldof::CellType::WHITE_DOWN, { 0, 1, 2, 3, 4, 5 } }, // unused, only for completeness
       { celldof::CellType::BLUE_UP, { 2, 4, 5, 3, X, 0 } },
       { celldof::CellType::BLUE_DOWN, { 2, 1, 0, 3, X, 5 } },
       { celldof::CellType::GREEN_UP, { 0, X, 5, 3, 4, 2 } },
       { celldof::CellType::GREEN_DOWN, { 5, X, 0, 3, 1, 2 } } };

   // For a given micro cell type and local macro face index, gives the local
   // macro edge indices (HyTeG ordering) that are contained in this macro face
   // *and* this micro cell type.
   static const std::map< celldof::CellType, std::map< uint_t, std::set< uint_t > > > edgesOwnedByFace{
       { celldof::CellType::WHITE_UP, { { 0, { 0, 1, 2 } }, { 1, { 0, 3, 4 } }, { 2, { 1, 3, 5 } }, { 3, { 2, 4, 5 } } } },
       { celldof::CellType::WHITE_DOWN, {} }, // unused, only for completeness
       { celldof::CellType::BLUE_UP, { { 0, { 0, 1, 2 } }, { 1, { 3 } }, { 3, { 5 } } } },
       { celldof::CellType::BLUE_DOWN, { { 1, { 0 } }, { 2, { 1, 3, 5 } }, { 3, { 2 } } } },
       { celldof::CellType::GREEN_UP, { { 0, { 2 } }, { 1, { 0, 3, 4 } }, { 2, { 5 } } } },
       { celldof::CellType::GREEN_DOWN, { { 0, { 0 } }, { 2, { 3 } }, { 3, { 2, 4, 5 } } } } };
   // https://xkcd.com/297

   switch ( cellType )
   {
   case celldof::CellType::WHITE_UP:
      owningEdges = indexing::isOnCellEdge( microCell, levelinfo::num_microedges_per_edge( level ) );
      owningFaces = indexing::isOnCellFace( microCell, levelinfo::num_microedges_per_edge( level ) );
      break;
   case celldof::CellType::WHITE_DOWN:
      // All edges are owned by the macro-cell.
      break;
   case celldof::CellType::BLUE_UP:
      owningFaces = indexing::isOnCellFace( microCell + Index{ 1, 0, 0 }, levelinfo::num_microedges_per_edge( level ) );
      break;
   case celldof::CellType::BLUE_DOWN:
      owningFaces = indexing::isOnCellFace( microCell + Index{ 0, 0, 1 }, levelinfo::num_microedges_per_edge( level ) );
      break;
   case celldof::CellType::GREEN_UP:
      owningFaces = indexing::isOnCellFace( microCell, levelinfo::num_microedges_per_edge( level ) );
      break;
   case celldof::CellType::GREEN_DOWN:
      owningFaces = indexing::isOnCellFace( microCell + Index{ 0, 1, 0 }, levelinfo::num_microedges_per_edge( level ) );
      break;
   default:
      WALBERLA_ABORT( "Invalid cell type" );
   }

   for ( const uint_t localEdgeId : owningEdges )
   {
      const uint_t edgeVertex0 = cell.getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeId ).at( 0 );
      const uint_t edgeVertex1 = cell.getEdgeLocalVertexToCellLocalVertexMaps().at( localEdgeId ).at( 1 );
      const real_t sign        = ( edgeVertex1 < edgeVertex0 ) ? real_c( -1 ) : real_c( 1 );
      transform.diagonal()( dofIdx.at( cellType )[localEdgeId] ) = sign;
   }

   for ( const uint_t localFaceId : owningFaces )
   {
      const std::map< uint_t, uint_t > faceLocalVertexToCellLocalVertexMap =
          cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceId );

      for ( const auto& [v0, v1] : std::initializer_list< std::pair< uint_t, uint_t > >{ { 0, 1 }, { 0, 2 }, { 1, 2 } } )
      {
         const uint_t faceVertex0 = faceLocalVertexToCellLocalVertexMap.at( v0 );
         const uint_t faceVertex1 = faceLocalVertexToCellLocalVertexMap.at( v1 );

         const real_t       sign   = ( faceVertex1 < faceVertex0 ) ? real_c( -1 ) : real_c( 1 );
         const uint_t       edgeId = indexing::getCellLocalEdgeIDFromCellLocalVertexIDs( faceVertex0, faceVertex1 );
         const Eigen::Index dof    = dofIdx.at( cellType )[edgeId];

         if ( owningEdges.count( edgeId ) == 0 && edgesOwnedByFace.at( cellType ).at( localFaceId ).count( edgeId ) > 0 )
         {
            WALBERLA_ASSERT_UNEQUAL( dof, X );
            transform.diagonal()( dof ) = sign;
         }
      }
   }

   return transform;
}

} // namespace macrocell
} // namespace n1e1
} // namespace hyteg
