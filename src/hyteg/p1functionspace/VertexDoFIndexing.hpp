/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include <vector>

#include "core/DataTypes.h"

#include "hyteg/StencilDirections.hpp"
#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"
#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"

namespace hyteg {

using walberla::uint_t;
namespace vertexdof {

// ##################
// ### Macro Edge ###
// ##################

namespace macroedge {

uint_t neighborFaceGhostLayerSize( const uint_t& level );

uint_t neighborCellGhostLayerSize( const uint_t& level );

/// Index of a vertex DoF on a macro edge (only access to owned DoFs, no ghost layers).
uint_t index( const uint_t& level, const idx_t& x );

/// 'Enumerates' the inner vertexdofs. x is still the index but must not be located on the boundary.
uint_t innerIndex( const uint_t& level, const idx_t& x );

/// Index of a vertex DoF on a ghost layer of a macro edge.
/// \param level Refinement level
/// \param x index on the macro-edge
/// \param neighbor 0 to access the first neighbor data, 1 to access second neighbor, ...
uint_t indexOnNeighborFace( const uint_t& level, const idx_t& x, const uint_t& neighbor );

uint_t indexOnNeighborCell( const uint_t& level, const idx_t& x, const uint_t& neighbor, const uint_t& numNeighborFaces );

// Stencil access functions

uint_t
    indexFromVertexOnNeighborCell( const uint_t& level, const idx_t& x, const uint_t& cellID, const uint_t& numNeighborFaces );

uint_t indexFromVertexOnNeighborFace( const uint_t& level, const idx_t& x, const uint_t& faceID, const stencilDirection& dir );

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
uint_t indexFromVertex( const uint_t& level, const idx_t& x, const stencilDirection& dir );

uint_t stencilIndexOnEdge( const stencilDirection& dir );

uint_t stencilIndexOnNeighborFace( const stencilDirection dir, const uint_t faceID );

uint_t stencilIndexOnNeighborCell( const uint_t& cellID, const uint_t& numNeighborFaces );

/// Have a look into the documentation to understand the calculations here
/// The west vertices have the same col index as the horizonal edge
uint_t indexFromHorizontalEdge( const uint_t& level, const idx_t& x, const stencilDirection& dir );

/// neighbor arrays from vertex dof
constexpr std::array< stencilDirection, 7 > neighborsWithCenter          = {{hyteg::stencilDirection::VERTEX_C,
                                                                    hyteg::stencilDirection::VERTEX_S,
                                                                    hyteg::stencilDirection::VERTEX_SE,
                                                                    hyteg::stencilDirection::VERTEX_E,
                                                                    hyteg::stencilDirection::VERTEX_N,
                                                                    hyteg::stencilDirection::VERTEX_NW,
                                                                    hyteg::stencilDirection::VERTEX_W}};
constexpr std::array< stencilDirection, 2 > neighborsOnEdgeFromVertexDoF = {
    {hyteg::stencilDirection::VERTEX_E, hyteg::stencilDirection::VERTEX_W}};
constexpr std::array< stencilDirection, 2 > neighborsOnSouthFaceFromVertexDoF = {
    {hyteg::stencilDirection::VERTEX_S, hyteg::stencilDirection::VERTEX_SE}};
constexpr std::array< stencilDirection, 2 > neighborsOnNorthFaceFromVertexDoF = {
    {hyteg::stencilDirection::VERTEX_N, hyteg::stencilDirection::VERTEX_NW}};

/// neighbor arrays need to connect vertex dof and edge dof
constexpr std::array< stencilDirection, 2 > neighborsOnEdgeFromHorizontalEdgeDoF = {
    {hyteg::stencilDirection::VERTEX_E, hyteg::stencilDirection::VERTEX_W}};
constexpr std::array< stencilDirection, 1 > neighborsOnSouthFaceFromHorizontalEdgeDoF = {{hyteg::stencilDirection::VERTEX_SE}};
constexpr std::array< stencilDirection, 1 > neighborsOnNorthFaceFromHorizontalEdgeDoF = {{hyteg::stencilDirection::VERTEX_NW}};

/// Iterator over a vertex DoF macro edge.
/// See \ref EdgeIterator for more information.
class Iterator : public hyteg::indexing::EdgeIterator
{
 public:
   explicit Iterator( const uint_t& level, const uint_t& offsetToCenter = 0, const bool& backwards = false );
};

/// map[neighborCellID][indexOffset] = weight
typedef std::map< uint_t, std::map< indexing::IndexIncrement, real_t > > StencilMap_T;

} // namespace macroedge

// ##################
// ### Macro Face ###
// ##################

namespace macroface {

/// Direct access functions

/// Index of a vertex DoF on a macro face (only access to owned DoFs, no ghost layers).
uint_t index( const uint_t& level, const idx_t& x, const idx_t& y );

/// 'Enumerates' the inner vertexdofs. x and y are still the indices but must not be located on the boundary.
uint_t innerIndex( const uint_t& level, const idx_t& x, const idx_t& y );

/// Index of a vertex DoF on a ghost layer of a macro face.
/// \param level Refinement level
/// \param x x-index on the macro-face
/// \param y y-index on the macro-face
/// \param neighbor 0 or 1 for the respective neighbor
uint_t index( const uint_t& level, const idx_t& x, const idx_t& y, const uint_t& neighbor );

// Stencil access functions

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
uint_t indexFromVertex( const uint_t& level, const idx_t& x, const idx_t& y, const stencilDirection& dir );

/// Have a look into the documentation to understand the calculations here
/// The west vertex has the same col and row index as the horizonal edge
uint_t indexFromHorizontalEdge( const uint_t& level, const idx_t& x, const idx_t& y, const stencilDirection& dir );

constexpr std::array< stencilDirection, 4 > neighborsFromHorizontalEdge = {
    {stencilDirection::VERTEX_SE, stencilDirection::VERTEX_E, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_W}};

/// Have a look into the documentation to understand the calculations here
/// The south west vertex has the same col and row index as the horizonal edge
uint_t indexFromDiagonalEdge( const uint_t& level, const idx_t& x, const idx_t& y, const stencilDirection& dir );

constexpr std::array< stencilDirection, 4 > neighborsFromDiagonalEdge = {
    {stencilDirection::VERTEX_SE, stencilDirection::VERTEX_NE, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_SW}};

/// Have a look into the documentation to understand the calculations here
/// The south vertex has the same col and row index as the horizonal edge
uint_t indexFromVerticalEdge( const uint_t& level, const idx_t& x, const idx_t& y, const stencilDirection& dir );

uint_t indexFromGrayFace( const uint_t& level, const idx_t& x, const idx_t& y, const stencilDirection& dir );

uint_t indexFromBlueFace( const uint_t& level, const idx_t& x, const idx_t& y, const stencilDirection& dir );

constexpr std::array< stencilDirection, 13 > neighborsWithOneNeighborCellWithCenter = {{
    hyteg::stencilDirection::VERTEX_C,
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_TS,
    hyteg::stencilDirection::VERTEX_TSE,
    hyteg::stencilDirection::VERTEX_TSW,
    hyteg::stencilDirection::VERTEX_TNW,
}};

constexpr std::array< stencilDirection, 19 > neighborsWithTwoNeighborCellsWithCenter = {{
    hyteg::stencilDirection::VERTEX_C,   hyteg::stencilDirection::VERTEX_S,   hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_E,   hyteg::stencilDirection::VERTEX_N,   hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,   hyteg::stencilDirection::VERTEX_TC,  hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_TS,  hyteg::stencilDirection::VERTEX_TSE, hyteg::stencilDirection::VERTEX_TSW,
    hyteg::stencilDirection::VERTEX_TNW, hyteg::stencilDirection::VERTEX_BC,  hyteg::stencilDirection::VERTEX_BW,
    hyteg::stencilDirection::VERTEX_BS,  hyteg::stencilDirection::VERTEX_BSE, hyteg::stencilDirection::VERTEX_BSW,
    hyteg::stencilDirection::VERTEX_BNW,
}};

constexpr std::array< stencilDirection, 12 > neighborsWithOneNeighborCellWithoutCenter = {{hyteg::stencilDirection::VERTEX_S,
                                                                                           hyteg::stencilDirection::VERTEX_SE,
                                                                                           hyteg::stencilDirection::VERTEX_E,
                                                                                           hyteg::stencilDirection::VERTEX_N,
                                                                                           hyteg::stencilDirection::VERTEX_NW,
                                                                                           hyteg::stencilDirection::VERTEX_W,
                                                                                           hyteg::stencilDirection::VERTEX_TC,
                                                                                           hyteg::stencilDirection::VERTEX_TW,
                                                                                           hyteg::stencilDirection::VERTEX_TS,
                                                                                           hyteg::stencilDirection::VERTEX_TSE,
                                                                                           hyteg::stencilDirection::VERTEX_TSW,
                                                                                           hyteg::stencilDirection::VERTEX_TNW}};

constexpr std::array< stencilDirection, 18 > neighborsWithTwoNeighborCellsWithoutCenter = {{hyteg::stencilDirection::VERTEX_S,
                                                                                            hyteg::stencilDirection::VERTEX_SE,
                                                                                            hyteg::stencilDirection::VERTEX_E,
                                                                                            hyteg::stencilDirection::VERTEX_N,
                                                                                            hyteg::stencilDirection::VERTEX_NW,
                                                                                            hyteg::stencilDirection::VERTEX_W,
                                                                                            hyteg::stencilDirection::VERTEX_TC,
                                                                                            hyteg::stencilDirection::VERTEX_TW,
                                                                                            hyteg::stencilDirection::VERTEX_TS,
                                                                                            hyteg::stencilDirection::VERTEX_TSE,
                                                                                            hyteg::stencilDirection::VERTEX_TSW,
                                                                                            hyteg::stencilDirection::VERTEX_TNW,
                                                                                            hyteg::stencilDirection::VERTEX_BC,
                                                                                            hyteg::stencilDirection::VERTEX_BW,
                                                                                            hyteg::stencilDirection::VERTEX_BS,
                                                                                            hyteg::stencilDirection::VERTEX_BSE,
                                                                                            hyteg::stencilDirection::VERTEX_BSW,
                                                                                            hyteg::stencilDirection::VERTEX_BNW}};

constexpr std::array< stencilDirection, 7 > neighborsWithCenter    = {{hyteg::stencilDirection::VERTEX_C,
                                                                    hyteg::stencilDirection::VERTEX_S,
                                                                    hyteg::stencilDirection::VERTEX_SE,
                                                                    hyteg::stencilDirection::VERTEX_E,
                                                                    hyteg::stencilDirection::VERTEX_N,
                                                                    hyteg::stencilDirection::VERTEX_NW,
                                                                    hyteg::stencilDirection::VERTEX_W}};
constexpr std::array< stencilDirection, 6 > neighborsWithoutCenter = {{hyteg::stencilDirection::VERTEX_S,
                                                                       hyteg::stencilDirection::VERTEX_SE,
                                                                       hyteg::stencilDirection::VERTEX_E,
                                                                       hyteg::stencilDirection::VERTEX_N,
                                                                       hyteg::stencilDirection::VERTEX_NW,
                                                                       hyteg::stencilDirection::VERTEX_W}};

constexpr std::array< stencilDirection, 4 > neighborsFromVerticalEdge = {
    {stencilDirection::VERTEX_S, stencilDirection::VERTEX_SE, stencilDirection::VERTEX_N, stencilDirection::VERTEX_NW}};

constexpr std::array< stencilDirection, 3 > neighborsFromGrayFace = {
    {stencilDirection::VERTEX_SW, stencilDirection::VERTEX_SE, stencilDirection::VERTEX_NW}};
constexpr std::array< stencilDirection, 3 > neighborsFromBlueFace = {
    {stencilDirection::VERTEX_SE, stencilDirection::VERTEX_NW, stencilDirection::VERTEX_NE}};

// Iterators

/// Iterator over a vertex DoF macro face.
/// See \ref FaceIterator for more information.
class Iterator : public hyteg::indexing::FaceIterator
{
 public:
   explicit Iterator( const uint_t& level, const uint_t& offsetToCenter = 0 );
};

/// Iterator over the border of a vertex DoF macro face.
/// See FaceBorderIterator for more information.
class BoundaryIterator : public hyteg::indexing::FaceBoundaryIterator
{
 public:
   BoundaryIterator( const uint_t&                                 level,
                     const hyteg::indexing::FaceBoundaryDirection& direction,
                     const uint_t&                                 offsetToCenter     = 0,
                     const uint_t&                                 offsetFromVertices = 0 );
};

bool isVertexOnBoundary( const uint_t& level, const hyteg::indexing::Index& idx );

// map[neighborCellID][indexOffset] = weight
typedef std::map< uint_t, std::map< indexing::IndexIncrement, real_t > > StencilMap_T;

} // namespace macroface

// ##################
// ### Macro Cell ###
// ##################

namespace macrocell {

/// Index of a vertex DoF on a macro cell.
uint_t index( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z );

/// Index of neighboring vertices of a vertex DoF specified by the coordinates.
uint_t indexFromVertex( const uint_t& level, const idx_t& x, const idx_t& y, const idx_t& z, const stencilDirection& dir );

/// neighbor arrays from vertex dof
constexpr std::array< stencilDirection, 15 > neighborsWithCenter = {{
    hyteg::stencilDirection::VERTEX_C,
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_TS,
    hyteg::stencilDirection::VERTEX_TSE,
    hyteg::stencilDirection::VERTEX_BC,
    hyteg::stencilDirection::VERTEX_BN,
    hyteg::stencilDirection::VERTEX_BE,
    hyteg::stencilDirection::VERTEX_BNW,
}};

constexpr std::array< stencilDirection, 14 > neighborsWithoutCenter = {{
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_TS,
    hyteg::stencilDirection::VERTEX_TSE,
    hyteg::stencilDirection::VERTEX_BC,
    hyteg::stencilDirection::VERTEX_BN,
    hyteg::stencilDirection::VERTEX_BE,
    hyteg::stencilDirection::VERTEX_BNW,
}};

const std::vector< stencilDirection > neighborsOnFace0WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_TS,
    hyteg::stencilDirection::VERTEX_TSE,
}};

const std::vector< stencilDirection > neighborsOnFace1WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_BC,
    hyteg::stencilDirection::VERTEX_BN,
    hyteg::stencilDirection::VERTEX_BE,
    hyteg::stencilDirection::VERTEX_BNW,
}};

const std::vector< stencilDirection > neighborsOnFace2WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_TS,
    hyteg::stencilDirection::VERTEX_TSE,
    hyteg::stencilDirection::VERTEX_BC,
    hyteg::stencilDirection::VERTEX_BN,
    hyteg::stencilDirection::VERTEX_BE,
}};

const std::vector< stencilDirection > neighborsOnFace3WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_TS,
    hyteg::stencilDirection::VERTEX_BC,
    hyteg::stencilDirection::VERTEX_BN,
    hyteg::stencilDirection::VERTEX_BE,
    hyteg::stencilDirection::VERTEX_BNW,
}};

const std::vector< stencilDirection > neighborsOnEdge0WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_TW,
}};

const std::vector< stencilDirection > neighborsOnEdge1WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_TS,
    hyteg::stencilDirection::VERTEX_TSE,
}};

const std::vector< stencilDirection > neighborsOnEdge2WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_TS,
}};

const std::vector< stencilDirection > neighborsOnEdge3WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_E,
    hyteg::stencilDirection::VERTEX_N,
    hyteg::stencilDirection::VERTEX_TC,
    hyteg::stencilDirection::VERTEX_BC,
    hyteg::stencilDirection::VERTEX_BN,
    hyteg::stencilDirection::VERTEX_BE,
}};

const std::vector< stencilDirection > neighborsOnEdge4WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_NW,
    hyteg::stencilDirection::VERTEX_W,
    hyteg::stencilDirection::VERTEX_TW,
    hyteg::stencilDirection::VERTEX_BC,
    hyteg::stencilDirection::VERTEX_BN,
    hyteg::stencilDirection::VERTEX_BE,
    hyteg::stencilDirection::VERTEX_BNW,
}};

const std::vector< stencilDirection > neighborsOnEdge5WithoutCenter = {{
    hyteg::stencilDirection::VERTEX_S,
    hyteg::stencilDirection::VERTEX_SE,
    hyteg::stencilDirection::VERTEX_TS,
    hyteg::stencilDirection::VERTEX_BC,
    hyteg::stencilDirection::VERTEX_BN,
    hyteg::stencilDirection::VERTEX_BE,
}};

const std::vector< stencilDirection > neighborsOnVertex0WithoutCenter = {
    {hyteg::stencilDirection::VERTEX_E, hyteg::stencilDirection::VERTEX_N, hyteg::stencilDirection::VERTEX_TC}};

const std::vector< stencilDirection > neighborsOnVertex1WithoutCenter = {
    {hyteg::stencilDirection::VERTEX_NW, hyteg::stencilDirection::VERTEX_W, hyteg::stencilDirection::VERTEX_TW}};

const std::vector< stencilDirection > neighborsOnVertex2WithoutCenter = {
    {hyteg::stencilDirection::VERTEX_S, hyteg::stencilDirection::VERTEX_SE, hyteg::stencilDirection::VERTEX_TS}};

const std::vector< stencilDirection > neighborsOnVertex3WithoutCenter = {
    {hyteg::stencilDirection::VERTEX_BC, hyteg::stencilDirection::VERTEX_BN, hyteg::stencilDirection::VERTEX_BE}};

const std::array< std::vector< stencilDirection >, 4 > neighborsOnFaceWithoutCenter = {
    {neighborsOnFace0WithoutCenter, neighborsOnFace1WithoutCenter, neighborsOnFace2WithoutCenter, neighborsOnFace3WithoutCenter}};

const std::array< std::vector< stencilDirection >, 6 > neighborsOnEdgeWithoutCenter = {{neighborsOnEdge0WithoutCenter,
                                                                                        neighborsOnEdge1WithoutCenter,
                                                                                        neighborsOnEdge2WithoutCenter,
                                                                                        neighborsOnEdge3WithoutCenter,
                                                                                        neighborsOnEdge4WithoutCenter,
                                                                                        neighborsOnEdge5WithoutCenter}};

const std::array< std::vector< stencilDirection >, 4 > neighborsOnVertexWithoutCenter = {{neighborsOnVertex0WithoutCenter,
                                                                                          neighborsOnVertex1WithoutCenter,
                                                                                          neighborsOnVertex2WithoutCenter,
                                                                                          neighborsOnVertex3WithoutCenter}};

// Iterators

/// Iterator over a vertex DoF macro cell.
/// See \ref CellIterator for more information.
class Iterator : public hyteg::indexing::CellIterator
{
 public:
   explicit Iterator( const uint_t& level, const uint_t& offsetToCenter = 0 );
};

/// Iterator over the borders (faces) of a macro-cell.
/// See CellBorderIterator for more information.
class BoundaryIterator : public hyteg::indexing::CellBoundaryIterator
{
 public:
   BoundaryIterator( const uint_t& level,
                     const uint_t& vertex0,
                     const uint_t& vertex1,
                     const uint_t& vertex2,
                     const uint_t& offsetToCenter = 0 );
};

/// See \ref indexing::isOnCellFace
std::set< uint_t > isOnCellFace( const indexing::Index& index, const uint_t& level );

/// See \ref indexing::isOnCellEdge
std::set< uint_t > isOnCellEdge( const indexing::Index& index, const uint_t& level );

/// See \ref indexing::isOnCellVertex
std::set< uint_t > isOnCellVertex( const indexing::Index& index, const uint_t& level );

/// map[indexOffset] = weight
typedef std::map< indexing::IndexIncrement, real_t > StencilMap_T;

} // namespace macrocell

// ################
// ### Stencils ###
// ################

/// Returns the logical index offset from a micro-vertex resulting from moving in the passed stencil direction
indexing::IndexIncrement logicalIndexOffsetFromVertex( const stencilDirection& dir );

/// Returns the logical index offset from a micro-vertex resulting from moving in the passed stencil direction.
stencilDirection stencilDirectionFromLogicalOffset( const indexing::IndexIncrement& offset );

uint_t stencilIndexFromVertex( const stencilDirection dir );

uint_t stencilIndexFromHorizontalEdge( const stencilDirection dir );

uint_t stencilIndexFromDiagonalEdge( const stencilDirection dir );

uint_t stencilIndexFromVerticalEdge( const stencilDirection dir );

uint_t stencilIndexFromGrayFace( const stencilDirection& dir );

uint_t stencilIndexFromBlueFace( const stencilDirection& dir );


// ##############
// ### Others ###
// ##############
void getVertexDoFDataIndicesFromMicroFace( const indexing::Index & microFaceIndex,
                                           const facedof::FaceType & faceType,
                                           const uint_t level,
                                           std::array< uint_t, 3 >& vertexDoFIndices );

void getVertexDoFDataIndicesFromMicroFace( const indexing::Index & microFaceIndex,
                                           const facedof::FaceType & faceType,
                                           const uint_t level,
                                           std::vector< uint_t >& vertexDoFIndices );

void getVertexDoFDataIndicesFromMicroCell( const indexing::Index & microCellIndex,
                                           const celldof::CellType & cellType,
                                           const uint_t level,
                                           std::array<uint_t, 4>& vertexDoFIndices );

        void getVertexDoFDataIndicesFromMicroCell( const indexing::Index & microCellIndex,
                                                   const celldof::CellType & cellType,
                                                   const uint_t level,
                                                   std::vector< uint_t >& vertexDoFIndices );

} // namespace vertexdof
} // namespace hyteg
