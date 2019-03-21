#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFApply.hpp"
#include "tinyhhg_core/LevelWiseMemory.hpp"
#include "tinyhhg_core/p2functionspace/generated_new/P2FenicsForm.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/p2functionspace/generated/p2_divt.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_diffusion.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_mass.h"
#include "tinyhhg_core/mixedoperators/generated/p1_to_p2_tet_divt_tet.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif


namespace hhg {

template< class VertexDoFToEdgeDoFForm >
class VertexDoFToEdgeDoFOperator : public Operator<P1Function< real_t >, EdgeDoFFunction< real_t > >
{
public:
  VertexDoFToEdgeDoFOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel);
  ~VertexDoFToEdgeDoFOperator() final = default;

  void apply(const P1Function< real_t > & src,const EdgeDoFFunction< real_t > & dst, size_t level, DoFType flag, UpdateType updateType = Replace) const;

  /// since the Vertex does not own any EdgeDoFs only edge, face and cell are needed
  const PrimitiveDataID< StencilMemory< real_t >, Edge> &getEdgeStencilID() const { return edgeStencilID_; }
  const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge> &getEdgeStencil3DID() const { return edgeStencil3DID_; }
  const PrimitiveDataID< StencilMemory< real_t >, Face> &getFaceStencilID() const { return faceStencilID_; }
  const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face> &getFaceStencil3DID() const { return faceStencil3DID_; }
  const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell> &getCellStencilID() const { return cellStencilID_; }

private:
  void assembleStencils();

  PrimitiveDataID< StencilMemory< real_t >, Edge> edgeStencilID_;
  PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge> edgeStencil3DID_;
  PrimitiveDataID< StencilMemory< real_t >, Face> faceStencilID_;
  PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face> faceStencil3DID_;
  PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell> cellStencilID_;

  VertexDoFToEdgeDoFForm form;
};

template< typename UFCOperator3D >
void assembleVertexToEdgeStencils( const std::shared_ptr< PrimitiveStorage > & storage,
                                   const uint_t & minLevel,
                                   const uint_t & maxLevel,
                                   const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge> & macroEdgeStencilID,
                                   const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face> & macroFaceStencilID,
                                   const PrimitiveDataID< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell> & macroCellStencilID  )
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

        auto& vertexToEdgeStencilMemory = edge.getData( macroEdgeStencilID )->getData( level );
        {
          const auto convertedCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
                  edgedof::EdgeDoFOrientation::X, basisInCell.at( 0 ), basisInCell.at( 1 ), basisInCell.at( 2 ) );
          vertexToEdgeStencilMemory[neighborCellID][convertedCenterOrientation] = P2Elements::P2Elements3D::calculateVertexToEdgeStencilInMacroCell(
                  edgeAssemblyIndexInCell, convertedCenterOrientation, cell, level, ufcOperator );
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

        auto& vertexToEdgeStencilMemory = face.getData( macroFaceStencilID )->getData( level );
        for( const auto& centerOrientation : faceOrientations )
        {
          const auto convertedCenterOrientation = edgedof::convertEdgeDoFOrientationFaceToCell(
                  centerOrientation, localVertexIDsAtCell.at(0), localVertexIDsAtCell.at(1), localVertexIDsAtCell.at(2));
          vertexToEdgeStencilMemory[neighborCellID][convertedCenterOrientation] =
                  P2Elements::P2Elements3D::calculateVertexToEdgeStencilInMacroCell(
                          edgeAssemblyIndexInCell, convertedCenterOrientation, cell, level, ufcOperator );
        }
      }
    }

    /////////////////
    // Macro-cells //
    /////////////////

    for( const auto& it : storage->getCells() )
    {
      const auto &cell = *it.second;

      auto& vertexToEdgeStencilMemory = cell.getData( macroCellStencilID )->getData( level );
      for( const auto& centerOrientation : edgedof::allEdgeDoFOrientations )
      {
        const auto vertexToEdgeStencilMap = P2Elements::P2Elements3D::calculateVertexToEdgeStencilInMacroCell(
                edgedof::macrocell::getInnerIndexByOrientation( centerOrientation ),
                centerOrientation,
                cell,
                level,
                ufcOperator );
        for( const auto stencilIt : vertexToEdgeStencilMap )
        {
          vertexToEdgeStencilMemory[centerOrientation][stencilIt.first] = stencilIt.second;
        }
      }
    }
  }
}

namespace VertexDoFToEdgeDoF {

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

/// \param level stencil size is independent of level
/// \param primitive \ref Primitive the memory is allocated on
/// \return number of the stencil entries
uint_t macroCellVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive );

}

typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< hhg::fenics::NoAssemble, fenics::NoAssemble > > GenericVertexDoFToEdgeDoFOperator;
typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise > > VertexToEdgeDivTxOperator;
typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise > > VertexToEdgeDivTyOperator;

typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > > P1ToP2DivTxVertexToEdgeOperator;
typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > > P1ToP2DivTyVertexToEdgeOperator;
typedef VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > > P1ToP2DivTzVertexToEdgeOperator;

}/// namespace hhg
