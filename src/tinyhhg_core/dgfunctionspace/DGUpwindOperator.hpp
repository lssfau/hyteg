#pragma once

#include <tinyhhg_core/Operator.hpp>

#include <array>
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg
{

template<class VelocityBaseType>
class DGUpwindOperator : public Operator< DGFunction< real_t >, DGFunction< real_t > >
{
  typedef std::array<std::shared_ptr<VelocityBaseType>, 2> VelocityType;

public:
  DGUpwindOperator(const std::shared_ptr< PrimitiveStorage > & storage, const VelocityType& velocity, size_t minLevel, size_t maxLevel)
    : Operator(storage, minLevel, maxLevel), velocity_(velocity)
  {
  }

  ~DGUpwindOperator()
  {
  }

private:

  void apply_impl(DGFunction< real_t >& src, DGFunction< real_t >& dst, uint_t level, DoFType flag, UpdateType updateType = Replace)
  {
    // start pulling edge halos
    src.getCommunicator(level)->startCommunication<Face, Edge>();

    // end pulling edge halos
    src.getCommunicator(level)->endCommunication<Face, Edge>();

    // start pulling vertex halos
    src.getCommunicator(level)->startCommunication<Edge, Vertex>();

    // end pulling vertex halos
    src.getCommunicator(level)->endCommunication<Edge, Vertex>();

    for(auto velocityComponent : velocity_) {
      velocityComponent->getCommunicator(level)->template startCommunication<Edge, Vertex>();
    }


    for(auto velocityComponent : velocity_) {
      velocityComponent->getCommunicator(level)->template startCommunication<Face, Edge>();
    }


    for(auto velocityComponent : velocity_) {
      velocityComponent->getCommunicator(level)->template endCommunication<Edge, Vertex>();
    }

    for (auto& it : storage_->getVertices()) {
      Vertex& vertex = *it.second;

      if (testFlag(vertex.getDoFType(), flag))
      {
        DGVertex::upwind< real_t >(level, vertex, storage_, src.getVertexDataID(), dst.getVertexDataID(), std::array<PrimitiveDataID< FunctionMemory< real_t >, Vertex>, 2>{{velocity_[0]->getVertexDataID(), velocity_[1]->getVertexDataID()}}, updateType);
      }
    }

    dst.getCommunicator(level)->startCommunication<Vertex, Edge>();

    for(auto velocityComponent : velocity_) {
      velocityComponent->getCommunicator(level)->template endCommunication<Face, Edge>();
    }

    for (auto& it : storage_->getEdges()) {
      Edge& edge = *it.second;

      if (testFlag(edge.getDoFType(), flag))
      {
        DGEdge::upwind< real_t >(level, edge, storage_, src.getEdgeDataID(), dst.getEdgeDataID(), std::array<PrimitiveDataID< FunctionMemory< real_t >, Edge>, 2>{{velocity_[0]->getEdgeDataID(), velocity_[1]->getEdgeDataID()}}, updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Vertex, Edge>();

    dst.getCommunicator(level)->startCommunication<Edge, Face>();

    for (auto& it : storage_->getFaces()) {
      Face& face = *it.second;

      if (testFlag(face.type, flag))
      {
        DGFace::upwind<real_t>(level, face, storage_, src.getFaceDataID(), dst.getFaceDataID(), std::array<PrimitiveDataID< FunctionMemory< real_t >, Face>, 2>{{velocity_[0]->getFaceDataID(), velocity_[1]->getFaceDataID()}}, updateType);
      }
    }

    dst.getCommunicator(level)->endCommunication<Edge, Face>();
  }

private:
  VelocityType velocity_;
};

}


