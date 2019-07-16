#pragma once

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/LevelWiseMemory.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFOperatorTypeDefs.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/Algorithms.hpp"
#include "tinyhhg_core/indexing/DistanceCoordinateSystem.hpp"
#include "tinyhhg_core/Levelinfo.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements3D.hpp"
#include "tinyhhg_core/p2functionspace/generated_new/P2FenicsForm.hpp"

namespace hhg{

template< class EdgeDoFForm >
class EdgeDoFOperator : public Operator<hhg::EdgeDoFFunction< real_t >, hhg::EdgeDoFFunction < real_t > >
{
public:

  EdgeDoFOperator(const std::shared_ptr <PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel);
  ~EdgeDoFOperator() final = default;

  void apply(const EdgeDoFFunction< real_t >& src,const  EdgeDoFFunction< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType) const;

  const PrimitiveDataID<StencilMemory< real_t >, Edge  > &getEdgeStencilID() const;
  const PrimitiveDataID<LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge  > &getEdgeStencil3DID() const;
  const PrimitiveDataID<StencilMemory< real_t >, Face  > &getFaceStencilID() const;
  const PrimitiveDataID<LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face  > &getFaceStencil3DID() const;
  const PrimitiveDataID<LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell  > &getCellStencilID() const;

private:

  void assembleStencils();

  PrimitiveDataID<StencilMemory< real_t >, Edge  > edgeStencilID_;
  PrimitiveDataID<LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge  > edgeStencil3DID_;
  PrimitiveDataID<StencilMemory< real_t >, Face  > faceStencilID_;
  PrimitiveDataID<LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face  > faceStencil3DID_;
  PrimitiveDataID<LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell  > cellStencilID_;

  EdgeDoFForm form;
};

/// on edges only one stencil is required since only the horizontal edge DoFs belong to the edge
uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// on face three stencils are needed for horizontal, vertical and diagonal DoFs
uint_t macroFaceEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

uint_t macroCellEdgeDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

template< typename UFCOperator3D >
void assembleEdgeToEdgeStencils( const std::shared_ptr< PrimitiveStorage > & storage,
                                 const uint_t & minLevel,
                                 const uint_t & maxLevel,
                                 const PrimitiveDataID<LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge> & macroEdgeStencilID,
                                 const PrimitiveDataID<LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face> & macroFaceStencilID,
                                 const PrimitiveDataID<LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell> & macroCellStencilID  )
{
    UFCOperator3D ufcOperator;

    for( uint_t level = minLevel; level <= maxLevel; ++level )
    {
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

                const auto edgeAssemblyIndexInCell = indexing::basisConversion(
                        indexing::Index(0, 0, 0), basisInCell, {0, 1, 2, 3}, levelinfo::num_microedges_per_edge(level));

                auto &edgeToEdgeStencilMemory = edge.getData( macroEdgeStencilID )->getData(level);

                for (const auto &leafOrientation : edgedof::allEdgeDoFOrientations) {
                    const auto convertedCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
                            edgedof::EdgeDoFOrientation::X,
                            basisInCell.at(0),
                            basisInCell.at(1),
                            basisInCell.at(2));
                    edgeToEdgeStencilMemory[neighborCellID][convertedCenterOrientation][leafOrientation] =
                            P2Elements::P2Elements3D::calculateEdgeToEdgeStencilInMacroCell(
                                    edgeAssemblyIndexInCell, convertedCenterOrientation, leafOrientation, cell, level,
                                    ufcOperator);
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

                const std::array<edgedof::EdgeDoFOrientation, 3> faceOrientations = {
                        edgedof::EdgeDoFOrientation::X, edgedof::EdgeDoFOrientation::Y,
                        edgedof::EdgeDoFOrientation::XY};

                const uint_t localFaceID = cell.getLocalFaceID(face.getID());
                const std::array<uint_t, 4> localVertexIDsAtCell = {
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(0),
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(1),
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(2),
                        6 - cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(0) -
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(1) -
                        cell.getFaceLocalVertexToCellLocalVertexMaps().at(localFaceID).at(2)};

                const auto edgeAssemblyIndexInCell = indexing::basisConversion(
                        indexing::Index(1, 1, 0), localVertexIDsAtCell, {0, 1, 2, 3},
                        levelinfo::num_microedges_per_edge(level));

                auto &edgeToEdgeStencilMemory = face.getData(macroFaceStencilID)->getData(level);
                for (const auto &centerOrientation : faceOrientations) {
                    for (const auto &leafOrientation : edgedof::allEdgeDoFOrientations) {
                        const auto convertedCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
                                centerOrientation,
                                localVertexIDsAtCell.at(0),
                                localVertexIDsAtCell.at(1),
                                localVertexIDsAtCell.at(2));
                        edgeToEdgeStencilMemory[neighborCellID][convertedCenterOrientation][leafOrientation] =
                                P2Elements::P2Elements3D::calculateEdgeToEdgeStencilInMacroCell(
                                        edgeAssemblyIndexInCell, convertedCenterOrientation, leafOrientation, cell,
                                        level, ufcOperator);
                    }
                }
            }
        }

        /////////////////
        // Macro-cells //
        /////////////////

        for( const auto& it : storage->getCells() )
        {
            const auto &cell = *it.second;

            auto& edgeToEdgeStencilMemory = cell.getData( macroCellStencilID )->getData( level );
            for( const auto& centerOrientation : edgedof::allEdgeDoFOrientations )
            {
                for( const auto& leafOrientation : edgedof::allEdgeDoFOrientations )
                {
                    const auto edgeToEdgeStencilMap = P2Elements::P2Elements3D::calculateEdgeToEdgeStencilInMacroCell(
                            edgedof::macrocell::getInnerIndexByOrientation( centerOrientation ),
                            centerOrientation,
                            leafOrientation,
                            cell,
                            level,
                            ufcOperator );
                    for( const auto stencilIt : edgeToEdgeStencilMap )
                    {
                        edgeToEdgeStencilMemory[centerOrientation][leafOrientation][stencilIt.first] = stencilIt.second;
                    }
                }
            }
        }
    }
}

typedef EdgeDoFOperator< P2FenicsForm< hhg::fenics::NoAssemble, fenics::NoAssemble > > GenericEdgeDoFOperator;

}
