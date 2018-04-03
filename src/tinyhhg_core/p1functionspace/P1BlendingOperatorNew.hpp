
#pragma once

#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"
#include "P1DataHandling.hpp"

#include "tinyhhg_core/fenics/fenics.hpp"

#include "tinyhhg_core/p1functionspace/generated_new/P1FormLaplace.hpp"

#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"

#include <tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp>
#include <tinyhhg_core/p1functionspace/VertexDoFBlendingNew.hpp>

namespace hhg
{

template<class P1Form>
class P1BlendingOperatorNew : public Operator< P1Function< real_t >, P1Function< real_t > >
{
public:
  P1BlendingOperatorNew(const std::shared_ptr< PrimitiveStorage > & storage,
                        size_t minLevel,
                        size_t maxLevel)
    : Operator(storage, minLevel, maxLevel)
  {
  }

  ~P1BlendingOperatorNew()
  {
  }

private:
  void apply_impl(P1Function< real_t >& src, P1Function< real_t >& dst, size_t level, DoFType flag, UpdateType updateType = Replace)
  {
    // start pulling vertex halos
    src.getCommunicator(level)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    src.getCommunicator(level)->startCommunication<Face, Edge>();

    src.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
//        vertexdof::blendingnew::macrovertex::applyBlending< real_t >(level, vertex, storage_, src.getVertexDataID(), dst.getVertexDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    src.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
//        vertexdof::blending::macroedge::applyBlending< real_t >(level, edge, storage_, edgeLocalMatrixIDs_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        vertexdof::blendingnew::macroface::applyBlending< real_t, P1Form >(level, face, form, src.getFaceDataID(), dst.getFaceDataID(), updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_gs_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag)
  {
    // start pulling vertex halos
    dst.getCommunicator(level)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    dst.getCommunicator(level)->startCommunication<Face, Edge>();

    // end pulling vertex halos
    dst.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
//        vertexdof::blending::macrovertex::smooth_gs_blending(level, vertex, storage_, vertexLocalMatrixIDs_, dst.getVertexDataID(), rhs.getVertexDataID());
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    dst.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
//        vertexdof::blending::macroedge::smoothGSBlending<real_t>(level, edge, storage_, edgeLocalMatrixIDs_, dst.getEdgeDataID(), rhs.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
//        vertexdof::blending::macroface::smoothGSBlending<real_t>(level, face, faceLocalMatrixIDs_, dst.getFaceDataID(), rhs.getFaceDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

  void smooth_jac_impl(P1Function< real_t >& dst, P1Function< real_t >& rhs, P1Function< real_t >& tmp, size_t level, DoFType flag)
  {
    // start pulling vertex halos
    tmp.getCommunicator(level)->startCommunication<Edge, Vertex>();

    // start pulling edge halos
    tmp.getCommunicator(level)->startCommunication<Face, Edge>();

    // end pulling vertex halos
    tmp.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Vertex::smooth_jac(vertex, vertexLocalMatrixID_, dst.getVertexDataID(), rhs.getVertexDataID(), tmp.getVertexDataID(), level);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    // end pulling edge halos
    tmp.getCommunicator(level)->endCommunication<Face, Edge>();

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Edge::smooth_jac(level, edge, edgeLocalMatrixID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), tmp.getEdgeDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Face::smooth_jac(level, face, faceLocalMatrixID_, dst.getFaceDataID(), rhs.getFaceDataID(), tmp.getFaceDataID());
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

#ifdef HHG_BUILD_WITH_PETSC
  void createMatrix_impl(P1Function< real_t >& src, P1Function< real_t >& dst, Mat& mat, size_t level, DoFType flag)
  {
    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Vertex::saveOperator(vertex, vertexLocalMatrixID_, src.getVertexDataID(), dst.getVertexDataID(), mat, level);
      }
    }

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Edge::saveOperator(level, edge, edgeLocalMatrixID_, src.getEdgeDataID(), dst.getEdgeDataID(), mat);
      }
    }

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        WALBERLA_ABORT("To be implemented")
//        P1Face::saveOperator(level, face, faceLocalMatrixID_, src.getFaceDataID(), dst.getFaceDataID(), mat);
      }
    }
  }
#endif

  P1Form form;
};

typedef P1BlendingOperatorNew<P1FormLaplace> P1BlendingLaplaceOperatorNew;

}


