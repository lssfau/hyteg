#include "VertexDoFToEdgeDoFOperator.hpp"

#include "tinyhhg_core/p2functionspace/P2Elements.hpp"
#include "generatedKernels/GeneratedKernelsVertexToEdgeMacroFace2D.hpp"
#include "generatedKernels/GeneratedKernelsVertexToEdgeMacroCell3D.hpp"

namespace hhg {

template< class UFCOperator2D, class UFCOperator3D >
VertexDoFToEdgeDoFOperator< UFCOperator2D, UFCOperator3D >::VertexDoFToEdgeDoFOperator(const std::shared_ptr<PrimitiveStorage> &storage, size_t minLevel, size_t maxLevel)
  : Operator(storage, minLevel, maxLevel) {
  /// since the Vertex does not own any EdgeDoFs only edge and face are needed

  auto edgeDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Edge >>(minLevel_,
                                                                        maxLevel_,
                                                                        VertexDoFToEdgeDoF::macroEdgeVertexDoFToEdgeDoFStencilSize);

  auto edge3DDataHandling =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge > >( minLevel_, maxLevel_ );

  auto faceDataHandling =
    std::make_shared< MemoryDataHandling<StencilMemory<real_t>, Face >>(minLevel_,
                                                                        maxLevel_,
                                                                        VertexDoFToEdgeDoF::macroFaceVertexDoFToEdgeDoFStencilSize);

  auto face3DDataHandling =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face > >( minLevel_, maxLevel_ );

  auto cellDataHandling =
    std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell > >( minLevel_, maxLevel_ );

  storage->addEdgeData(edgeStencilID_, edgeDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil");
  storage->addEdgeData(edgeStencil3DID_, edge3DDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil3D");
  storage->addFaceData(faceStencilID_, faceDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil");
  storage->addFaceData(faceStencil3DID_, face3DDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil3D");
  storage->addCellData(cellStencilID_, cellDataHandling, "VertexDoFToEdgeDoFOperatorCellStencil");

  if ( this->getStorage()->hasGlobalCells() )
  {
    if ( !std::is_same< UFCOperator3D, fenics::NoAssemble >::value )
    {
      assembleVertexToEdgeStencils< UFCOperator3D >( this->getStorage(),
                                                     this->minLevel_,
                                                     this->maxLevel_,
                                                     getEdgeStencil3DID(),
                                                     getFaceStencil3DID(),
                                                     getCellStencilID() );
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
void VertexDoFToEdgeDoFOperator< UFCOperator2D, UFCOperator3D >::assembleStencils() {
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

      // Assemble vertexToEdge stencil
      real_t * vStencil = storage_->getFace(face.getID())->getData(faceStencilID_)->getPointer(level);
      P2Face::VertexToEdge::assembleStencil(local_stiffness_gray, local_stiffness_blue, vStencil);
//      WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Face = {}", PointND<real_t, 12>(&vStencil[0])));
    }

    // Assemble edge stencils
    for (auto& it : storage_->getEdges()) {
      Edge &edge = *it.second;

      // Assemble vertexToEdge stencil
      Face* face = storage_->getFace(edge.neighborFaces()[0]);
      real_t* vStencil = storage_->getEdge(edge.getID())->getData(edgeStencilID_)->getPointer(level);
      compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
      compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
      P2Edge::VertexToEdge::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true);

      if (edge.getNumNeighborFaces() == 2) {
        face = storage_->getFace(edge.neighborFaces()[1]);
        compute_local_stiffness(*face, level, local_stiffness_gray, fenics::GRAY);
        compute_local_stiffness(*face, level, local_stiffness_blue, fenics::BLUE);
        P2Edge::VertexToEdge::assembleStencil(edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false);
      }

//      WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Edge = {}", PointND<real_t, 4>(&vStencil[0])));
    }

  }
}

template< class UFCOperator2D, class UFCOperator3D >
void VertexDoFToEdgeDoFOperator< UFCOperator2D, UFCOperator3D >::compute_local_stiffness(const Face &face, size_t level, Matrix6r& local_stiffness, fenics::ElementType element_type) {
  real_t coords[6];
  fenics::compute_micro_coords(face, level, coords, element_type);
  UFCOperator2D gen;
  gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
}

template < class UFCOperator2D, class UFCOperator3D >
void VertexDoFToEdgeDoFOperator< UFCOperator2D, UFCOperator3D >::apply( const P1Function< real_t >&      src,
                                                                        const EdgeDoFFunction< real_t >& dst,
                                                                        size_t                           level,
                                                                        DoFType                          flag,
                                                                        UpdateType                       updateType ) const
{
  this->startTiming( "Apply" );
  ///the order of communication is crucial here.
  ///first the vertex dofs on the macro vertex need to be communicated to the edge since they are needed on the edge and the face
  src.communicate<Vertex, Edge>( level );
  ///secondly the vertex dofs on the macro edge are communicated to the face passing on the vertex dof from the macro vertex
  src.communicate<Edge, Face>( level );
  src.communicate< Face, Cell >( level );
  src.communicate< Cell, Face >( level );
  ///lastly the vertex dofs on the macro face are communicated to the edge which also contain vertex dofs which are located on neighboring edges
  src.startCommunication<Face, Edge>( level );

  this->timingTree_->start( "Macro-Cell" );

  for (auto& it : storage_->getCells()) {
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
          VertexDoFToEdgeDoF::generated::apply_3D_macrocell_vertexdof_to_edgedof_add( &dstData[firstIdx[eo::X]],
                                                                                      &dstData[firstIdx[eo::XY]],
                                                                                      &dstData[firstIdx[eo::XYZ]],
                                                                                      &dstData[firstIdx[eo::XZ]],
                                                                                      &dstData[firstIdx[eo::Y]],
                                                                                      &dstData[firstIdx[eo::YZ]],
                                                                                      &dstData[firstIdx[eo::Z]],
                                                                                      srcData,
                                                                                      static_cast< int64_t >( level ),
                                                                                      stencilData );
       }
       else
       {
          VertexDoFToEdgeDoF::applyCell( level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
       }
    }
  }

  this->timingTree_->stop( "Macro-Cell" );

  this->timingTree_->start( "Macro-Face" );

  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;

    const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
    if (testFlag(faceBC, flag))
    {
      if ( storage_->hasGlobalCells() )
      {
        VertexDoFToEdgeDoF::applyFace3D( level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
      }
      else if( hhg::globalDefines::useGeneratedKernels )
      {
        real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
        real_t* vertexToDiagonalEdgeStencil   = &opr_data[4];
        real_t* vertexToHorizontalEdgeStencil = &opr_data[0];
        real_t* vertexToVerticalEdgeStencil   = &opr_data[8];
        real_t* src_data = face.getData( src.getFaceDataID() )->getPointer( level );
        real_t* dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
        if( updateType == hhg::Replace )
        {
          VertexDoFToEdgeDoF::generated::apply_2D_macroface_vertexdof_to_edgedof_replace( dst_data, src_data, vertexToDiagonalEdgeStencil, vertexToHorizontalEdgeStencil, vertexToVerticalEdgeStencil, static_cast< int64_t  >( level ) );
        } else if( updateType == hhg::Add )
        {
          VertexDoFToEdgeDoF::generated::apply_2D_macroface_vertexdof_to_edgedof_add( dst_data, src_data, vertexToDiagonalEdgeStencil, vertexToHorizontalEdgeStencil, vertexToVerticalEdgeStencil, static_cast< int64_t  >( level ) );
        }
      }
      else
      {
        VertexDoFToEdgeDoF::applyFace( level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
      }
    }
  }

  this->timingTree_->stop( "Macro-Face" );

  src.endCommunication<Face, Edge>( level );

  this->timingTree_->start( "Macro-Edge" );

  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;

    const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
    if (testFlag(edgeBC, flag))
    {
      if ( storage_->hasGlobalCells() )
      {
        VertexDoFToEdgeDoF::applyEdge3D( level, edge, *getStorage(), edgeStencil3DID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
      }
      else
      {
        VertexDoFToEdgeDoF::applyEdge( level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
      }
    }
  }

  this->timingTree_->stop( "Macro-Edge" );

  this->stopTiming( "Apply" );
}

namespace VertexDoFToEdgeDoF {


uint_t macroEdgeVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  return 2 + primitive.getNumNeighborFaces();
}

uint_t macroFaceVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 4 + 4 + 4;
}

uint_t macroCellVertexDoFToEdgeDoFStencilSize(const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  WALBERLA_UNUSED( primitive );
  return 7 * 27;
}
}

template class VertexDoFToEdgeDoFOperator< hhg::fenics::NoAssemble, hhg::fenics::NoAssemble >;
template class VertexDoFToEdgeDoFOperator< hhg::fenics::NoAssemble, hhg::fenics::UndefinedAssembly >;
template class VertexDoFToEdgeDoFOperator<p2_divt_cell_integral_0_otherwise>;
template class VertexDoFToEdgeDoFOperator<p2_divt_cell_integral_1_otherwise>;

template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise >;
template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_mass_cell_integral_0_otherwise >;
template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_pspg_tet_cell_integral_0_otherwise >;

template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_div_tet_cell_integral_0_otherwise >;
template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_div_tet_cell_integral_1_otherwise >;
template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise >;

template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_0_otherwise >;
template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_1_otherwise >;
template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise >;

template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise >;
template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise >;
template class VertexDoFToEdgeDoFOperator< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise >;

}/// namespace hhg
