#include "EdgeDoFToVertexDoFOperator.hpp"
#include "EdgeDoFToVertexDoFApply.hpp"
#include "generatedKernels/GeneratedKernelsEdgeToVertexMacroFace2D.hpp"

#include "tinyhhg_core/p2functionspace/variablestencil/P2VariableStencilCommon.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements.hpp"

#include "tinyhhg_core/mixedoperators/P2ToP1FenicsForm.hpp"

namespace hhg{

template< class EdgeDoFToVertexDoFForm >
EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::EdgeDoFToVertexDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage,
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
    if ( form.assemble3D() )
    {
      WALBERLA_ABORT("Not implemented.");
//      assembleEdgeToVertexStencils< UFCOperator3D >( this->getStorage(),
//                                                     this->minLevel_,
//                                                     this->maxLevel_,
//                                                     getVertexStencil3DID(),
//                                                     getEdgeStencil3DID(),
//                                                     getFaceStencil3DID(),
//                                                     getCellStencilID());
    }
  }
  else
  {
    // Only assemble stencils if UFCOperator is specified
    if ( form.assemble2D() )
    {
      assembleStencils();
    }
  }

}

template< class EdgeDoFToVertexDoFForm >
void EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::assembleStencils() {
  using namespace P2Elements;

  // Assemble stencils on all levels
  for (uint_t level = minLevel_; level <= maxLevel_; ++level)
  {

    // Assemble face stencils
    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

       uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );

       Point3D x( face.coords[0] );
       real_t  h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

       Point3D d0 = h * ( face.coords[1] - face.coords[0] );
       Point3D d2 = h * ( face.coords[2] - face.coords[0] );

       form.geometryMap = face.getGeometryMap();

       Point3D dirS  = -1.0 * d2;
       Point3D dirSE = d0 - 1.0 * d2;
       Point3D dirE  = d0;
       Point3D dirW  = -1.0 * d0;
       Point3D dirNW = -1.0 * d0 + d2;
       Point3D dirN  = d2;

       real_t* vStencil = storage_->getFace(face.getID())->getData(faceStencilID_)->getPointer(level);

       P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
           form, {x, x + dirW, x + dirS}, P2Elements::P2Face::elementSW_reord, vStencil );
       P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
           form, {x, x + dirS, x + dirSE}, P2Elements::P2Face::elementS_reord, vStencil );
       P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
           form, {x, x + dirSE, x + dirE}, P2Elements::P2Face::elementSE_reord, vStencil );
       P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
           form, {x, x + dirE, x + dirN}, P2Elements::P2Face::elementNE_reord, vStencil );
       P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
           form, {x, x + dirN, x + dirNW}, P2Elements::P2Face::elementN_reord, vStencil );
       P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
           form, {x, x + dirNW, x + dirW}, P2Elements::P2Face::elementNW_reord, vStencil );
    }

    // Assemble edge stencils
    for (auto& it : storage_->getEdges()) {
      Edge &edge = *it.second;
      real_t *vStencil = storage_->getEdge(edge.getID())->getData(edgeStencilID_)->getPointer(level);

      size_t rowsize = levelinfo::num_microvertices_per_edge(level);

      Face *faceS = storage_->getFace(edge.neighborFaces()[0]);
      Face *faceN = nullptr;
      uint_t s_south = faceS->vertex_index(edge.neighborVertices()[0]);
      uint_t e_south = faceS->vertex_index(edge.neighborVertices()[1]);
      uint_t o_south = faceS->vertex_index(faceS->get_vertex_opposite_to_edge(edge.getID()));

      real_t h = 1.0 / (walberla::real_c(rowsize - 1));

      Point3D dS_se = h * (faceS->coords[e_south] - faceS->coords[s_south]);
      Point3D dS_so = h * (faceS->coords[o_south] - faceS->coords[s_south]);
      Point3D dS_oe = h * (faceS->coords[e_south] - faceS->coords[o_south]);

      Point3D dir_S = -1.0 * dS_oe;
      Point3D dir_E = dS_se;
      Point3D dir_SE = dS_so;
      Point3D dir_W = -1.0 * dS_se;

      Point3D x = edge.getCoordinates()[0];
      Point3D dx = h * edge.getDirection();
      x += dx;

      uint_t s_north, e_north, o_north;
      Point3D dir_N;
      Point3D dir_NW;

      if (edge.getNumNeighborFaces() == 2) {
        faceN = storage_->getFace(edge.neighborFaces()[1]);
        s_north = faceN->vertex_index(edge.neighborVertices()[0]);
        e_north = faceN->vertex_index(edge.neighborVertices()[1]);
        o_north = faceN->vertex_index(faceN->get_vertex_opposite_to_edge(edge.getID()));

        Point3D dN_so = h * (faceN->coords[o_north] - faceN->coords[s_north]);
        Point3D dN_oe = h * (faceN->coords[e_north] - faceN->coords[o_north]);

        dir_N = dN_so;
        dir_NW = -1.0 * dN_oe;
      }

     // assemble south
     form.geometryMap = faceS->getGeometryMap();
     P2::variablestencil::assembleEdgeToVertexStencil<EdgeDoFToVertexDoFForm>(form, {x, x + dir_W, x + dir_S},
                                                         P2Elements::P2Face::elementSW_reord, vStencil);
     P2::variablestencil::assembleEdgeToVertexStencil<EdgeDoFToVertexDoFForm>(form, {x, x + dir_S, x + dir_SE},
                                                         P2Elements::P2Face::elementS_reord, vStencil);
     P2::variablestencil::assembleEdgeToVertexStencil<EdgeDoFToVertexDoFForm>(
         form, {x, x + dir_SE, x + dir_E}, P2Elements::P2Face::elementSE_reord, vStencil);

     if (edge.getNumNeighborFaces() == 2) {
       form.geometryMap = faceN->getGeometryMap();
       P2::variablestencil::assembleEdgeToVertexStencil<EdgeDoFToVertexDoFForm>(
           form, {x, x + dir_E, x + dir_N}, P2Elements::P2Face::elementNE_reord, vStencil);
       P2::variablestencil::assembleEdgeToVertexStencil<EdgeDoFToVertexDoFForm>(
           form, {x, x + dir_N, x + dir_NW}, P2Elements::P2Face::elementN_reord, vStencil);
       P2::variablestencil::assembleEdgeToVertexStencil<EdgeDoFToVertexDoFForm>(
           form, {x, x + dir_NW, x + dir_W}, P2Elements::P2Face::elementNW_reord, vStencil);
        }
    }

    for (auto& it : storage_->getVertices()) {
      Vertex &vertex = *it.second;

      // Assemble EdgeToVertex
      real_t* vStencil = storage_->getVertex(vertex.getID())->getData(vertexStencilID_)->getPointer(level);

       uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

       Point3D x;
       Point3D d0;
       Point3D d2;

       real_t h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

       uint_t neighborId = 0;
       for( auto& faceId : vertex.neighborFaces() )
       {
          Face* face       = storage_->getFace( faceId );
          form.geometryMap = face->getGeometryMap();

          uint_t                     v_i       = face->vertex_index( vertex.getID() );
          std::vector< PrimitiveID > adj_edges = face->adjacent_edges( vertex.getID() );

          x = face->coords[v_i];
          d0 =
              ( face->coords[face->vertex_index( storage_->getEdge( adj_edges[0] )->get_opposite_vertex( vertex.getID() ) )] - x ) * h;
          d2 =
              ( face->coords[face->vertex_index( storage_->getEdge( adj_edges[1] )->get_opposite_vertex( vertex.getID() ) )] - x ) * h;

          Point3D matrixRow;
          form.integrateEdgeToVertex( {{x, x + d0, x + d2}}, matrixRow );

          uint_t i = 1;
          // iterate over adjacent edges
          for( auto& edgeId : adj_edges )
          {
             uint_t edge_idx = vertex.edge_index( edgeId );
             vStencil[edge_idx] += matrixRow[3 - i];
             i += 1;
          }

          walberla::uint_t face_idx = vertex.getNumNeighborEdges() + vertex.face_index( face->getID() );
          vStencil[face_idx] += matrixRow[0];

          ++neighborId;
       }
    }
  }
}

template< class EdgeDoFToVertexDoFForm >
void EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::apply(const EdgeDoFFunction <real_t> &src,
                                                                       const P1Function<double> &dst,
                                                                       uint_t level,
                                                                       DoFType flag,
                                                                       UpdateType updateType) const
{
  using namespace EdgeDoFToVertexDoF;
  this->startTiming( "Apply" );

  ///there might be room for optimization in the communication. i.e. splitting communicate into start and end to overlap comm and calc

  src.communicate<Face, Cell>( level );

  for (auto& it : storage_->getCells()) {
    Cell& cell = *it.second;

    const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
    if (testFlag(cellBC, flag))
    {
      applyCell(level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType);
    }
  }

  src.communicate<Edge, Face>( level );
  src.communicate<Cell, Face>( level );

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


  src.communicate<Face, Edge>( level );

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

  src.communicate<Edge, Vertex>( level );

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
  this->stopTiming( "Apply" );
}

template< class EdgeDoFToVertexDoFForm >
const PrimitiveDataID<StencilMemory< real_t >, Vertex > &EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getVertexStencilID() const {
  return vertexStencilID_;
}

template< class EdgeDoFToVertexDoFForm >
const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex > &EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getVertexStencil3DID() const {
  return vertexStencil3DID_;
}

template< class EdgeDoFToVertexDoFForm >
const PrimitiveDataID<StencilMemory< real_t >, Edge > &EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getEdgeStencilID() const {
  return edgeStencilID_;
}

template< class EdgeDoFToVertexDoFForm >
const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge > &EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getEdgeStencil3DID() const {
  return edgeStencil3DID_;
}

template< class EdgeDoFToVertexDoFForm >
const PrimitiveDataID<StencilMemory< real_t >, Face > &EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getFaceStencilID() const {
  return faceStencilID_;
}

template< class EdgeDoFToVertexDoFForm >
const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face > &EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getFaceStencil3DID() const {
  return faceStencil3DID_;
}

template< class EdgeDoFToVertexDoFForm >
const PrimitiveDataID<LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell > &EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getCellStencilID() const {
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

template class EdgeDoFToVertexDoFOperator< P2FenicsForm< hhg::fenics::NoAssemble, hhg::fenics::NoAssemble > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_div_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_div_cell_integral_1_otherwise > >;

template class EdgeDoFToVertexDoFOperator< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2ToP1FenicsForm< fenics::NoAssemble,                     p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;

}/// namespace hhg
