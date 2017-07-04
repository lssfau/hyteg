#pragma once

#include "tinyhhg_core/communication/PackInfo.hpp"

namespace hhg {
namespace communication {

class P1BubblePackInfo : public PackInfo {

public:


    P1BubblePackInfo(uint_t level,
                     PrimitiveDataID<VertexP1BubbleFunctionMemory, Vertex> dataIDVertex,
                     PrimitiveDataID<EdgeP1BubbleFunctionMemory  , Edge>   dataIDEdge,
                     PrimitiveDataID<FaceP1BubbleFunctionMemory  , Face>   dataIDFace)
        : level_(level),
          dataIDVertex_(dataIDVertex),
          dataIDEdge_(dataIDEdge),
          dataIDFace_(dataIDFace)
    {

    }

    virtual void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer);
    virtual void unpackVertexFromEdge(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer);
    virtual void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver);

    virtual void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer);
    virtual void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer);
    virtual void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver);

    virtual void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer);
    virtual void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer);
    virtual void communicateLocalEdgeToFace(const Edge *sender, Face *receiver);

    virtual void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer);
    virtual void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer);
    virtual void communicateLocalFaceToEdge(const Face *sender, Edge *receiver);


private:
    uint_t level_;
    PrimitiveDataID<VertexP1BubbleFunctionMemory, Vertex> dataIDVertex_;
    PrimitiveDataID<EdgeP1BubbleFunctionMemory  , Edge>   dataIDEdge_;
    PrimitiveDataID<FaceP1BubbleFunctionMemory  , Face>   dataIDFace_;
};

} //namespace communication
} //namespace hhg