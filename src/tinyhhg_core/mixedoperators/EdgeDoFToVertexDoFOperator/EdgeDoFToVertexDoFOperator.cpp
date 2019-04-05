#include "EdgeDoFToVertexDoFOperator.hpp"
#include "EdgeDoFToVertexDoFApply.hpp"
#include "generatedKernels/GeneratedKernelsEdgeToVertexMacroFace2D.hpp"
#include "generatedKernels/GeneratedKernelsEdgeToVertexMacroCell3D.hpp"

#include "tinyhhg_core/p2functionspace/P2Elements.hpp"

namespace hhg{

template< class UFCOperator2D, class UFCOperator3D >
EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::EdgeDoFToVertexDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
                                                       const size_t & minLevel,
                                                       const size_t & maxLevel)
  :Operator(storage,minLevel,maxLevel)
{

  using namespace EdgeDoFToVertexDoF;

  auto vertexDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Vertex >>(minLevel_, maxLevel_, macroVertexEdgeDoFToVertexDoFStencilSize);

  auto vertex3DDataHandling   =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex > >(minLevel_, maxLevel_);

  auto edgeDataHandling   =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Edge   >>(minLevel_, maxLevel_, macroEdgeEdgeDoFToVertexDoFStencilSize);

  auto edge3DDataHandling   =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge > >(minLevel_, maxLevel_);
  
  auto faceDataHandling   =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Face   >>(minLevel_, maxLevel_, macroFaceEdgeDoFToVertexDoFStencilSize);

  auto face3DDataHandling   =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face > >(minLevel_, maxLevel_);

  auto cellDataHandling   =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell > >(minLevel_, maxLevel_);

  storage->addVertexData(vertexStencilID_, vertexDataHandling, "VertexDoFToEdgeDoFOperatorVertexStencil");
  storage->addVertexData(vertexStencil3DID_, vertex3DDataHandling, "VertexDoFToEdgeDoFOperatorVertexStencil3D");
  storage->addEdgeData(edgeStencilID_, edgeDataHandling  , "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addEdgeData(edgeStencil3DID_, edge3DDataHandling  , "VertexDoFToEdgeDoFOperatorEdgeStencil3D");
  storage->addFaceData(faceStencilID_, faceDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addFaceData(faceStencil3DID_, face3DDataHandling  , "VertexDoFToEdgeDoFOperatorFaceStencil3D");
  storage->addCellData(cellStencilID_, cellDataHandling  , "VertexDoFToEdgeDoFOperatorCellStencil");

  if ( this->getStorage()->hasGlobalCells() )
  {
    if ( !std::is_same< UFCOperator3D, fenics::NoAssemble >::value )
    {
      assembleEdgeToVertexStencils< UFCOperator3D >( this->getStorage(),
                                                     this->minLevel_,
                                                     this->maxLevel_,
                                                     getVertexStencil3DID(),
                                                     getEdgeStencil3DID(),
                                                     getFaceStencil3DID(),
                                                     getCellStencilID());
    }
  }
  else
  {
    // Only assemble stencils if UFCOperator is specified
    if ( !std::is_same< UFCOperator2D, fenics::NoAssemble >::value )
    {
      assembleStencils();
    }
  }

}

template< class UFCOperator2D, class UFCOperator3D >
void EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::assembleStencils() {
  using namespace P2Elements;

  // Initialize memory for local 6x6 matrices
  Matrix6r local_stiffness_gray;
  Matrix6r local_stiffness_blue;

  // Assemble stencils on all levels
  for (uint_t level = minLevel_; level <= maxLevel_; ++level)
  {

    // Assemble face stencils
    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      // Compute both local stiffness matrices
      compute_local_stiffness(face, level, local_stiffness_gray, fenics::GRAY);
      compute_local_stiffness(face, level, local_stiffness_blue, fenics::BLUE);

      // Assemble edgeToVertex stencil
      real_t* vStencil = storage_->getFace(face.getID())->getData(faceStencilID_)->getPointer(level);
      P2Face::EdgeToVertex::assembleStencil(local_stiffness_gray, local_stiffness_blue, vStencil);
//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Face = {}", PointND<real_t, 12>(&vStencil[0])));
    }

    // Assemble edge stencils
    for (auto& it : storage_->getEdges()) {
      Edge &edge = *it.second;

      // Assemble edgeToVertex
      Face* face = storage_->getFace(edge.neighborFaces()[0]);
      real_t* vStencil = storage_->getEdge(edge.getID())->getData(edgeStencilID_)->getPointer(level);
      compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
      compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
      P2Edge::EdgeToVertex::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true);

      if (edge.getNumNeighborFaces() == 2) {
        face = storage_->getFace(edge.neighborFaces()[1]);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
        P2Edge::EdgeToVertex::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false);
      }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Edge = {}", PointND<real_t, 7>(&vStencil[0])));
    }

    for (auto& it : storage_->getVertices()) {
      Vertex &vertex = *it.second;

      // Assemble EdgeToVertex
      real_t* vStencil = storage_->getVertex(vertex.getID())->getData(vertexStencilID_)->getPointer(level);
      for (auto& faceId : vertex.neighborFaces())
      {
        Face* face = storage_->getFace(faceId);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        P2Vertex::EdgeToVertex::assembleStencil(vertex, *face, local_stiffness_gray, vStencil, storage_);
      }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("edgeToVertex/Vertex = {}", PointND<real_t, 5>(&vStencil[0])));
    }

  }
}

template< class UFCOperator2D, class UFCOperator3D >
void EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type) {
  real_t coords[6];
  fenics::compute_micro_coords(face, level, coords, element_type);
  UFCOperator2D gen;
  gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
}

template< class UFCOperator2D, class UFCOperator3D >
void EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::apply(const EdgeDoFFunction <real_t> &src,
                                                                       const P1Function<double> &dst,
                                                                       uint_t level,
                                                                       DoFType flag,
                                                                       UpdateType updateType) const
{
  using namespace EdgeDoFToVertexDoF;
  this->startTiming( "Apply" );

  ///there might be room for optimization in the communication. i.e. splitting communicate into start and end to overlap comm and calc

  src.communicate< Face, Cell >( level );

  this->timingTree_->start( "Macro-Cell" );

  for ( auto& it : storage_->getCells() )
  {
     Cell& cell = *it.second;

     const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
     if ( testFlag( cellBC, flag ) )
     {
        if ( hhg::globalDefines::useGeneratedKernels && updateType == Add )
        {
           typedef edgedof::EdgeDoFOrientation eo;
           auto                                dstData     = cell.getData( dst.getCellDataID() )->getPointer( level );
           auto                                srcData     = cell.getData( src.getCellDataID() )->getPointer( level );
           auto                                stencilData = cell.getData( cellStencilID_ )->getData( level );
           std::map< eo, uint_t >              firstIdx;
           for ( auto e : edgedof::allEdgeDoFOrientations )
              firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );
           EdgeDoFToVertexDoF::generated::apply_3D_macrocell_edgedof_to_vertexdof_add( &srcData[firstIdx[eo::X]],
                                                                                       &srcData[firstIdx[eo::XY]],
                                                                                       &srcData[firstIdx[eo::XYZ]],
                                                                                       &srcData[firstIdx[eo::XZ]],
                                                                                       &srcData[firstIdx[eo::Y]],
                                                                                       &srcData[firstIdx[eo::YZ]],
                                                                                       &srcData[firstIdx[eo::Z]],
                                                                                       dstData,
                                                                                       stencilData,
                                                                                       static_cast< int64_t >( level ) );
        }
        else
        {
           applyCell( level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
        }
     }
  }

  this->timingTree_->stop( "Macro-Cell" );

  src.communicate<Edge, Face>( level );
  src.communicate<Cell, Face>( level );

  this->timingTree_->start( "Macro-Face" );

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag))
    {
      if ( storage_->hasGlobalCells() )
      {
        applyFace3D(level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
      else if( hhg::globalDefines::useGeneratedKernels )
      {
        real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
        real_t* src_data = face.getData( src.getFaceDataID() )->getPointer( level );
        real_t* dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
        if( updateType == hhg::Replace )
        {
          EdgeDoFToVertexDoF::generated::apply_2D_macroface_edgedof_to_vertexdof_replace( src_data, opr_data, dst_data, static_cast< int64_t >( level ) );
        } else if( updateType == hhg::Add )
        {
          EdgeDoFToVertexDoF::generated::apply_2D_macroface_edgedof_to_vertexdof_add( src_data, opr_data, dst_data, static_cast< int64_t >( level ) );
        }
      }
      else
      {
        applyFace(level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }
  }

  this->timingTree_->stop( "Macro-Face" );

  src.communicate<Face, Edge>( level );

  this->timingTree_->start( "Macro-Edge" );

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag))
    {
      if ( storage_->hasGlobalCells() )
      {
        applyEdge3D( level, edge, *getStorage(), edgeStencil3DID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
      }
      else
      {
        applyEdge( level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
      }
    }
  }

  this->timingTree_->stop( "Macro-Edge" );

  src.communicate<Edge, Vertex>( level );

  this->timingTree_->start( "Macro-Vertex" );

  for (auto& it : storage_->getVertices()) {
    Vertex& vertex = *it.second;

    const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
    if (testFlag(vertexBC, flag))
    {
      if ( storage_->hasGlobalCells() )
      {
        applyVertex3D( level, vertex, *getStorage(), vertexStencil3DID_, src.getVertexDataID(), dst.getVertexDataID(), updateType );
      }
      else
      {
        applyVertex( level, vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), updateType );
      }
    }
  }

  this->timingTree_->stop( "Macro-Vertex" );

  this->stopTiming( "Apply" );
}

template< class UFCOperator2D, class UFCOperator3D >
const PrimitiveDataID<StencilMemory< real_t >, Vertex > &EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::getVertexStencilID() const {
  return vertexStencilID_;
}

template< class UFCOperator2D, class UFCOperator3D >
const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex > &EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::getVertexStencil3DID() const {
  return vertexStencil3DID_;
}

template< class UFCOperator2D, class UFCOperator3D >
const PrimitiveDataID<StencilMemory< real_t >, Edge > &EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::getEdgeStencilID() const {
  return edgeStencilID_;
}

template< class UFCOperator2D, class UFCOperator3D >
const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge > &EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::getEdgeStencil3DID() const {
  return edgeStencil3DID_;
}

template< class UFCOperator2D, class UFCOperator3D >
const PrimitiveDataID<StencilMemory< real_t >, Face > &EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::getFaceStencilID() const {
  return faceStencilID_;
}

template< class UFCOperator2D, class UFCOperator3D >
const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face > &EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::getFaceStencil3DID() const {
  return faceStencil3DID_;
}

template< class UFCOperator2D, class UFCOperator3D >
const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell > &EdgeDoFToVertexDoFOperator< UFCOperator2D, UFCOperator3D >::getCellStencilID() const {
  return cellStencilID_;
}

namespace EdgeDoFToVertexDoF {
////////// Stencil sizes //////////
uint_t macroVertexEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive ) {
  WALBERLA_UNUSED(level);
  return primitive.getNumNeighborEdges() + primitive.getNumNeighborFaces();
}

uint_t macroEdgeEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive ) {
  WALBERLA_UNUSED(level);
  return 2 + 5 * primitive.getNumNeighborFaces();
}

uint_t macroFaceEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive ) {
  WALBERLA_UNUSED(level);
  WALBERLA_UNUSED(primitive);
  return 12;
}

uint_t macroCellEdgeDoFToVertexDoFStencilSize(const uint_t &level, const Primitive & primitive ) {
  WALBERLA_UNUSED(level);
  WALBERLA_UNUSED(primitive);
  return 7 * 27;
}

}/// namespace EdgeDoFToVertexDoF

template class EdgeDoFToVertexDoFOperator< hhg::fenics::NoAssemble, hhg::fenics::NoAssemble >;
template class EdgeDoFToVertexDoFOperator< hhg::fenics::NoAssemble, hhg::fenics::UndefinedAssembly >;
template class EdgeDoFToVertexDoFOperator<p2_div_cell_integral_0_otherwise>;
template class EdgeDoFToVertexDoFOperator<p2_div_cell_integral_1_otherwise>;

template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise >;
template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_mass_cell_integral_0_otherwise >;
template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_pspg_tet_cell_integral_0_otherwise >;

template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_div_tet_cell_integral_0_otherwise >;
template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_div_tet_cell_integral_1_otherwise >;
template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise >;

template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_0_otherwise >;
template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_1_otherwise >;
template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise >;

template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise >;
template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise >;
template class EdgeDoFToVertexDoFOperator< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise >;

}/// namespace hhg
