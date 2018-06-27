
#pragma once

#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "P1DataHandling.hpp"

#include "tinyhhg_core/p1functionspace/generated_new/P1FormLaplace.hpp"
#include "tinyhhg_core/p1functionspace/generated_new/P1FormMass.hpp"
#include "tinyhhg_core/p1functionspace/generated_new/P1FormEpsilon.hpp"
#include "tinyhhg_core/p1functionspace/generated_new/P1FormDivT.hpp"
#include "tinyhhg_core/p1functionspace/generated_new/P1FormDiv.hpp"
#include "tinyhhg_core/p1functionspace/generated_new/P1FormPSPG.hpp"

#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"

#include <tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp>
#include <tinyhhg_core/p1functionspace/blending/VertexDoFBlending.hpp>

namespace hhg
{

template<class P1Form>
class P1BlendingOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  P1BlendingOperator(const std::shared_ptr< PrimitiveStorage > & storage,
                        size_t minLevel,
                        size_t maxLevel)
    : Operator(storage, minLevel, maxLevel)
  {
  }

  ~P1BlendingOperator()
  {
  }

private:
  void apply_impl(P1Function< real_t >& src, P1Function< real_t >& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    src.communicate< Vertex, Edge>( level );
    src.communicate< Edge, Face>( level );
    src.communicate< Face, Edge>( level );
    src.communicate< Edge, Vertex>( level );

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if (testFlag(vertexBC, flag))
      {
        vertexdof::blending::macrovertex::applyBlending< real_t, P1Form >(level, vertex, form, storage_, src.getVertexDataID(), dst.getVertexDataID(), updateType);
      }
    }

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if (testFlag(edgeBC, flag))
      {
        vertexdof::blending::macroedge::applyBlending< real_t, P1Form >(level, edge, form, storage_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }
    }

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if (testFlag(faceBC, flag))
      {
        vertexdof::blending::macroface::applyBlending< real_t, P1Form >(level, face, form, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }

  }

  void smooth_gs_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag)
  {
    // start pulling vertex halos
    dst.startCommunication<Edge, Vertex>( level );

    // start pulling edge halos
    dst.startCommunication<Face, Edge>( level );

    // end pulling vertex halos
    dst.endCommunication<Edge, Vertex>( level );

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if (testFlag(vertexBC, flag))
      {
        vertexdof::blending::macrovertex::smoothGSBlending(level, vertex, form, storage_, dst.getVertexDataID(), rhs.getVertexDataID());
      }
    }

    dst.startCommunication<Vertex, Edge>( level );

    // end pulling edge halos
    dst.endCommunication<Face, Edge>( level );

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if (testFlag(edgeBC, flag))
      {
        vertexdof::blending::macroedge::smoothGSBlending<real_t>(level, edge, form, storage_, dst.getEdgeDataID(), rhs.getEdgeDataID());
      }
    }

    dst.endCommunication<Vertex, Edge>( level );

    dst.startCommunication<Edge, Face>( level );

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if (testFlag(faceBC, flag))
      {
        vertexdof::blending::macroface::smoothGSBlending<real_t>(level, face, form, dst.getFaceDataID(), rhs.getFaceDataID());
      }
    }

    dst.endCommunication<Edge, Face>( level );
  }

  void smooth_jac_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, P1Function< real_t >& tmp, size_t level, DoFType flag)
  {
    // start pulling vertex halos
    tmp.startCommunication<Edge, Vertex>( level );

    // start pulling edge halos
    tmp.startCommunication<Face, Edge>( level );

    // end pulling vertex halos
    tmp.endCommunication<Edge, Vertex>( level );

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if (testFlag(vertexBC, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Vertex::smooth_jac(vertex, vertexLocalMatrixID_, dst.getVertexDataID(), rhs.getVertexDataID(), tmp.getVertexDataID(), level);
      }
    }

    dst.startCommunication<Vertex, Edge>( level );

    // end pulling edge halos
    tmp.endCommunication<Face, Edge>( level );

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if (testFlag(edgeBC, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Edge::smooth_jac(level, edge, edgeLocalMatrixID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), tmp.getEdgeDataID());
      }
    }

    dst.endCommunication<Vertex, Edge>( level );

    dst.startCommunication<Edge, Face>( level );

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if (testFlag(faceBC, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Face::smooth_jac(level, face, faceLocalMatrixID_, dst.getFaceDataID(), rhs.getFaceDataID(), tmp.getFaceDataID());
      }
    }

    dst.endCommunication<Edge, Face>( level );
  }

#ifdef HHG_BUILD_WITH_PETSC
  void createMatrix_impl(P1Function< real_t >& src, P1Function< real_t >& dst, Mat& mat, size_t level, DoFType flag)
  {
    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if (testFlag(vertexBC, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Vertex::saveOperator(vertex, vertexLocalMatrixID_, src.getVertexDataID(), dst.getVertexDataID(), mat, level);
      }
    }

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if (testFlag(edgeBC, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Edge::saveOperator(level, edge, edgeLocalMatrixID_, src.getEdgeDataID(), dst.getEdgeDataID(), mat);
      }
    }

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if (testFlag(faceBC, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Face::saveOperator(level, face, faceLocalMatrixID_, src.getFaceDataID(), dst.getFaceDataID(), mat);
      }
    }
  }
#endif

  P1Form form;
};

typedef P1BlendingOperator<P1Form_laplace> P1BlendingLaplaceOperator;
typedef P1BlendingOperator<P1Form_mass> P1BlendingMassOperator;

typedef P1BlendingOperator<P1Form_epsilon_11> P1BlendingEpsilonOperator_11;
typedef P1BlendingOperator<P1Form_epsilon_12> P1BlendingEpsilonOperator_12;
typedef P1BlendingOperator<P1Form_epsilon_21> P1BlendingEpsilonOperator_21;
typedef P1BlendingOperator<P1Form_epsilon_22> P1BlendingEpsilonOperator_22;

typedef P1BlendingOperator<P1Form_divT_1> P1BlendingDivTOperator_1;
typedef P1BlendingOperator<P1Form_divT_2> P1BlendingDivTOperator_2;

typedef P1BlendingOperator<P1Form_div_1> P1BlendingDivOperator_1;
typedef P1BlendingOperator<P1Form_div_2> P1BlendingDivOperator_2;

typedef P1BlendingOperator<P1Form_pspg> P1BlendingPSPGOperator;

}


