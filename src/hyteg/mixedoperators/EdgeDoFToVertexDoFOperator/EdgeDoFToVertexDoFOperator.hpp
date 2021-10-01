/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hyteg {

template < class EdgeDoFToVertexDoFForm >
class EdgeDoFToVertexDoFOperator final : public Operator< hyteg::EdgeDoFFunction< real_t >, hyteg::P1Function< real_t > >
{
 public:
   EdgeDoFToVertexDoFOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                               const size_t&                              minLevel,
                               const size_t&                              maxLevel );
   EdgeDoFToVertexDoFOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                               const size_t&                              minLevel,
                               const size_t&                              maxLevel,
                               const EdgeDoFToVertexDoFForm&              form );
   ~EdgeDoFToVertexDoFOperator() = default;

   void apply( const EdgeDoFFunction< real_t >& src,
               const P1Function< double >&      dst,
               uint_t                           level,
               DoFType                          flag,
               UpdateType                       updateType ) const;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const EdgeDoFFunction< idx_t >&             src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const;

   const PrimitiveDataID< StencilMemory< real_t >, Vertex >&                                        getVertexStencilID() const;
   const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex >& getVertexStencil3DID() const;
   const PrimitiveDataID< StencilMemory< real_t >, Edge >&                                          getEdgeStencilID() const;
   const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >&     getEdgeStencil3DID() const;
   const PrimitiveDataID< StencilMemory< real_t >, Face >&                                          getFaceStencilID() const;
   const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face >&     getFaceStencil3DID() const;
   const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >&     getCellStencilID() const;

 private:
   void assembleStencils();

   PrimitiveDataID< StencilMemory< real_t >, Vertex >                                        vertexStencilID_;
   PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex > vertexStencil3DID_;
   PrimitiveDataID< StencilMemory< real_t >, Edge >                                          edgeStencilID_;
   PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >     edgeStencil3DID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >                                          faceStencilID_;
   PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face >     faceStencil3DID_;
   PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >     cellStencilID_;

   EdgeDoFToVertexDoFForm form_;
};

template < typename EdgeDoFToVertexDoFForm >
void assembleEdgeToVertexStencils(
    const std::shared_ptr< PrimitiveStorage >&                                                       storage,
    const uint_t&                                                                                    minLevel,
    const uint_t&                                                                                    maxLevel,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex >& macroVertexStencilID,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >&     macroEdgeStencilID,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face >&     macroFaceStencilID,
    const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >&     macroCellStencilID,
    const EdgeDoFToVertexDoFForm&                                                                    form )
{
   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      ////////////////////
      // Macro-vertices //
      ////////////////////

      for ( const auto& it : storage->getVertices() )
      {
         const auto& vertex = *it.second;
         WALBERLA_ASSERT_GREATER( vertex.getNumNeighborCells(), 0 );

         for ( uint_t neighborCellID = 0; neighborCellID < vertex.getNumNeighborCells(); neighborCellID++ )
         {
            const auto& cell = *( storage->getCell( vertex.neighborCells().at( neighborCellID ) ) );

            const uint_t cellLocalVertexID = cell.getLocalVertexID( vertex.getID() );
            const auto   basisInCell       = algorithms::getMissingIntegersAscending< 1, 4 >( {cellLocalVertexID} );

            const auto vertexAssemblyIndexInCell = indexing::basisConversion(
                indexing::Index( 0, 0, 0 ), basisInCell, {0, 1, 2, 3}, levelinfo::num_microvertices_per_edge( level ) );

            auto& edgeToVertexStencilMemory = vertex.getData( macroVertexStencilID )->getData( level );
            for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
            {
               edgeToVertexStencilMemory[neighborCellID][leafOrientation] =
                   P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                       vertexAssemblyIndexInCell, leafOrientation, cell, level, form );
               WALBERLA_ASSERT_EQUAL( edgeToVertexStencilMemory[neighborCellID][leafOrientation].size(), 1 );
            }
         }
      }

      /////////////////
      // Macro-edges //
      /////////////////

      if ( level >= 1 )
      {
         for ( const auto& it : storage->getEdges() )
         {
            const auto& edge = *it.second;
            WALBERLA_ASSERT_GREATER( edge.getNumNeighborCells(), 0 );

            for ( uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++ )
            {
               const auto& cell = *( storage->getCell( edge.neighborCells().at( neighborCellID ) ) );

               const uint_t cellLocalEdgeID = cell.getLocalEdgeID( edge.getID() );
               const auto   basisInCell     = algorithms::getMissingIntegersAscending< 2, 4 >(
                   {cell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 ),
                    cell.getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 )} );

               const auto vertexAssemblyIndexInCell = indexing::basisConversion(
                   indexing::Index( 1, 0, 0 ), basisInCell, {0, 1, 2, 3}, levelinfo::num_microvertices_per_edge( level ) );

               auto& edgeToVertexStencilMemory = edge.getData( macroEdgeStencilID )->getData( level );
               for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
               {
                  edgeToVertexStencilMemory[neighborCellID][leafOrientation] =
                      P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                          vertexAssemblyIndexInCell, leafOrientation, cell, level, form );
               }
            }
         }
      }

      /////////////////
      // Macro-faces //
      /////////////////

      if ( level >= 1 )
      {
         for ( const auto& it : storage->getFaces() )
         {
            const auto& face = *it.second;

            WALBERLA_ASSERT_GREATER( face.getNumNeighborCells(), 0 );

            for ( uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++ )
            {
               const auto& cell = *( storage->getCell( face.neighborCells().at( neighborCellID ) ) );

               const uint_t                  localFaceID          = cell.getLocalFaceID( face.getID() );
               const std::array< uint_t, 4 > localVertexIDsAtCell = {
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ),
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ),
                   cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ),
                   6 - cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ) -
                       cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ) -
                       cell.getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 )};
               const auto vertexAssemblyIndexInCell = indexing::basisConversion( indexing::Index( 1, 1, 0 ),
                                                                                 localVertexIDsAtCell,
                                                                                 {0, 1, 2, 3},
                                                                                 levelinfo::num_microvertices_per_edge( level ) );

               auto& edgeToVertexStencilMemory = face.getData( macroFaceStencilID )->getData( level );
               for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
               {
                  edgeToVertexStencilMemory[neighborCellID][leafOrientation] =
                      P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                          vertexAssemblyIndexInCell, leafOrientation, cell, level, form );
               }
            }
         }
      }

      /////////////////
      // Macro-cells //
      /////////////////
      if ( level >= 2 )
      {
         for ( const auto& it : storage->getCells() )
         {
            const auto& cell = *it.second;

            auto& edgeToVertexStencilMemory = cell.getData( macroCellStencilID )->getData( level );
            for ( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
            {
               const auto edgeToVertexStencilMap = P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                   indexing::Index( 1, 1, 1 ), leafOrientation, cell, level, form );
               for ( const auto& stencilIt : edgeToVertexStencilMap )
               {
                  edgeToVertexStencilMemory[leafOrientation][stencilIt.first] = stencilIt.second;
               }
            }
         }
      }
   }
}

namespace EdgeDoFToVertexDoF {

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro vertex, be aware that this one entry to large on the boundaries
uint_t macroVertexEdgeDoFToVertexDoFStencilSize( const uint_t& level, const Primitive& primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro edge
uint_t macroEdgeEdgeDoFToVertexDoFStencilSize( const uint_t& level, const Primitive& primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro face
uint_t macroFaceEdgeDoFToVertexDoFStencilSize( const uint_t& level, const Primitive& primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro cell
uint_t macroCellEdgeDoFToVertexDoFStencilSize( const uint_t& level, const Primitive& primitive );

} // namespace EdgeDoFToVertexDoF

typedef EdgeDoFToVertexDoFOperator< P2FenicsForm< hyteg::fenics::NoAssemble, hyteg::fenics::NoAssemble > >
                                                                                       GenericEdgeDoFToVertexDoFOperator;
typedef EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_div_cell_integral_0_otherwise > > EdgeToVertexDivxOperator;
typedef EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_div_cell_integral_1_otherwise > > EdgeToVertexDivyOperator;

typedef EdgeDoFToVertexDoFOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >
    P2ToP1DivxEdgeToVertexOperator;

typedef EdgeDoFToVertexDoFOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >
    P2ToP1DivyEdgeToVertexOperator;

typedef EdgeDoFToVertexDoFOperator< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >
    P2ToP1DivzEdgeToVertexOperator;

} // namespace hyteg
