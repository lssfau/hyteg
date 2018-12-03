#pragma once

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/LevelWiseMemory.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/p2functionspace/generated/p2_div.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_diffusion.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_mass.h"
#include "tinyhhg_core/mixedoperators/generated/p2_to_p1_tet_div_tet.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

namespace hhg {

template< class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class EdgeDoFToVertexDoFOperator : public Operator<hhg::EdgeDoFFunction< real_t >, hhg::P1Function < real_t > >
{
public:

  EdgeDoFToVertexDoFOperator(const std::shared_ptr <PrimitiveStorage> &storage, const size_t & minLevel, const size_t & maxLevel);
  ~EdgeDoFToVertexDoFOperator() final = default;

  void apply_impl(const EdgeDoFFunction <real_t> &src, const P1Function<double> &dst, uint_t level, DoFType flag, UpdateType updateType) const override;

  const PrimitiveDataID< StencilMemory< real_t >, Vertex> &getVertexStencilID() const;
  const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex  > &getVertexStencil3DID() const;
  const PrimitiveDataID< StencilMemory< real_t >, Edge  > &getEdgeStencilID() const;
  const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge  > &getEdgeStencil3DID() const;
  const PrimitiveDataID< StencilMemory< real_t >, Face  > &getFaceStencilID() const;
  const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face  > &getFaceStencil3DID() const;
  const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell  > &getCellStencilID() const;

private:

  void assembleStencils();

  void compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type);

  PrimitiveDataID<StencilMemory< real_t >, Vertex> vertexStencilID_;
  PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex > vertexStencil3DID_;
  PrimitiveDataID<StencilMemory< real_t >, Edge  > edgeStencilID_;
  PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge > edgeStencil3DID_;
  PrimitiveDataID<StencilMemory< real_t >, Face  > faceStencilID_;
  PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face > faceStencil3DID_;
  PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell > cellStencilID_;

};

template< typename UFCOperator3D >
void assembleEdgeToVertexStencils( const std::shared_ptr< PrimitiveStorage > & storage,
                                   const uint_t & minLevel,
                                   const uint_t & maxLevel,
                                   const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex > & macroVertexStencilID,
                                   const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge > & macroEdgeStencilID,
                                   const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face > & macroFaceStencilID,
                                   const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell > & macroCellStencilID  )
{
    UFCOperator3D ufcOperator;

    for( uint_t level = minLevel; level <= maxLevel; ++level )
    {
        ////////////////////
        // Macro-vertices //
        ////////////////////

        for ( const auto &it : storage->getVertices() )
        {
            const auto &vertex = *it.second;
            WALBERLA_ASSERT_GREATER(vertex.getNumNeighborCells(), 0);

            for (uint_t neighborCellID = 0; neighborCellID < vertex.getNumNeighborCells(); neighborCellID++) {
                const auto &cell = *( storage->getCell(vertex.neighborCells().at(neighborCellID)) );

                const uint_t cellLocalVertexID = cell.getLocalVertexID(vertex.getID());
                const auto basisInCell = algorithms::getMissingIntegersAscending<1, 4>({ cellLocalVertexID });

                const auto vertexAssemblyIndexInCell = indexing::basisConversion(
                  indexing::Index(0, 0, 0), basisInCell, {0, 1, 2, 3},
                  levelinfo::num_microvertices_per_edge(level));

                auto& edgeToVertexStencilMemory = vertex.getData( macroVertexStencilID )->getData( level );
                for( const auto& leafOrientation : edgedof::allEdgeDoFOrientationsWithoutXYZ )
                {
                    edgeToVertexStencilMemory[neighborCellID][leafOrientation] =
                    P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                    vertexAssemblyIndexInCell, leafOrientation, cell, level, ufcOperator );
                    WALBERLA_ASSERT_EQUAL( edgeToVertexStencilMemory[neighborCellID][leafOrientation].size(), 1 );
                }
            }
        }

        /////////////////
        // Macro-edges //
        /////////////////

        for ( const auto &it : storage->getEdges() )
        {
            const auto &edge = *it.second;
            WALBERLA_ASSERT_GREATER(edge.getNumNeighborCells(), 0);

            for (uint_t neighborCellID = 0; neighborCellID < edge.getNumNeighborCells(); neighborCellID++) {
                const auto &cell = *( storage->getCell(edge.neighborCells().at(neighborCellID)) );

                const uint_t cellLocalEdgeID = cell.getLocalEdgeID(edge.getID());
                const auto basisInCell = algorithms::getMissingIntegersAscending<2, 4>(
                        {cell.getEdgeLocalVertexToCellLocalVertexMaps().at(cellLocalEdgeID).at(0),
                         cell.getEdgeLocalVertexToCellLocalVertexMaps().at(cellLocalEdgeID).at(1)});

                const auto vertexAssemblyIndexInCell = indexing::basisConversion(
                        indexing::Index(1, 0, 0), basisInCell, {0, 1, 2, 3},
                        levelinfo::num_microvertices_per_edge(level));

                auto& edgeToVertexStencilMemory = edge.getData( macroEdgeStencilID )->getData( level );
                for( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
                {
                    edgeToVertexStencilMemory[neighborCellID][leafOrientation] =
                            P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                                    vertexAssemblyIndexInCell, leafOrientation, cell, level, ufcOperator );
                }
            }
        }

        /////////////////
        // Macro-faces //
        /////////////////

        for (const auto &it : storage->getFaces()) {
            const auto &face = *it.second;

            WALBERLA_ASSERT_GREATER(face.getNumNeighborCells(), 0);

            for (uint_t neighborCellID = 0; neighborCellID < face.getNumNeighborCells(); neighborCellID++) {
                const auto &cell = *(storage->getCell(face.neighborCells().at(neighborCellID)));

                const uint_t localFaceID = cell.getLocalFaceID(face.getID());
                const std::array<uint_t, 4> localVertexIDsAtCell = {
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(0),
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(1),
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(2),
                        6 - cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(0) -
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(1) -
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(2)};
                const auto vertexAssemblyIndexInCell = indexing::basisConversion(
                        indexing::Index(1, 1, 0), localVertexIDsAtCell, {0, 1, 2, 3},
                        levelinfo::num_microvertices_per_edge(level));

                auto& edgeToVertexStencilMemory = face.getData( macroFaceStencilID )->getData( level );
                for( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
                {
                    edgeToVertexStencilMemory[neighborCellID][leafOrientation] =
                            P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                                    vertexAssemblyIndexInCell, leafOrientation, cell, level, ufcOperator );
                }
            }
        }

        /////////////////
        // Macro-cells //
        /////////////////

        for( const auto& it : storage->getCells() )
        {
            const auto &cell = *it.second;

            auto& edgeToVertexStencilMemory = cell.getData( macroCellStencilID )->getData( level );
            for( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
            {
                const auto edgeToVertexStencilMap = P2Elements::P2Elements3D::calculateEdgeToVertexStencilInMacroCell(
                        indexing::Index( 1, 1, 1 ), leafOrientation, cell, level, ufcOperator );
                for( const auto stencilIt : edgeToVertexStencilMap )
                {
                    edgeToVertexStencilMemory[leafOrientation][stencilIt.first] = stencilIt.second;
                }
            }
        }
    }
}


namespace EdgeDoFToVertexDoF {

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro vertex, be aware that this one entry to large on the boundaries
uint_t macroVertexEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro edge
uint_t macroEdgeEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro face
uint_t macroFaceEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level the stencil size is independent of the level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of entries in the stencil on a macro cell
uint_t macroCellEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive );

}/// namespace EdgeDoFToVertexDoF

typedef EdgeDoFToVertexDoFOperator< hhg::fenics::NoAssemble, hhg::fenics::NoAssemble > GenericEdgeDoFToVertexDoFOperator;
typedef EdgeDoFToVertexDoFOperator<p2_div_cell_integral_0_otherwise> EdgeToVertexDivxOperator;
typedef EdgeDoFToVertexDoFOperator<p2_div_cell_integral_1_otherwise> EdgeToVertexDivyOperator;

typedef EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > P2ToP1DivxEdgeToVertexOperator;
typedef EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > P2ToP1DivyEdgeToVertexOperator;
typedef EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > P2ToP1DivzEdgeToVertexOperator;

}/// namespace hhg
