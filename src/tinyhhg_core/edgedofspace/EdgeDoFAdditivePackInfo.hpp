#pragma once

#include "core/DataTypes.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include "tinyhhg_core/communication/DoFSpacePackInfo.hpp"
#include "tinyhhg_core/boundary/BoundaryConditions.hpp"

namespace hhg {

using walberla::uint_t;

template< typename ValueType >
class FunctionMemory;
template < typename DataType, typename PrimitiveType >
class PrimitiveDataID;

class Vertex;
class Edge;
class Face;
class Cell;
class PrimitiveStorage;
class PrimitiveID;

template < typename ValueType >
class EdgeDoFAdditivePackInfo : public communication::DoFSpacePackInfo< ValueType >
{
 public:
   EdgeDoFAdditivePackInfo( uint_t                                         level,
                    PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dataIDVertex,
                    PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dataIDEdge,
                    PrimitiveDataID< FunctionMemory< ValueType >, Face >   dataIDFace,
                    PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dataIDCell,
                    std::weak_ptr< PrimitiveStorage >                      storage ,
                    BoundaryCondition                                      boundaryCondition,
                    DoFType                                                boundaryTypeToSkip ) :
     communication::DoFSpacePackInfo< ValueType >( level, dataIDVertex, dataIDEdge, dataIDFace, dataIDCell, storage ),
     boundaryCondition_( boundaryCondition ), boundaryTypeToSkip_( boundaryTypeToSkip )
   {}

  void packVertexForEdge(const Vertex *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackEdgeFromVertex(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalVertexToEdge(const Vertex *sender, Edge *receiver) const override;

  void packEdgeForVertex(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackVertexFromEdge(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalEdgeToVertex(const Edge *sender, Vertex *receiver) const override;

  void packEdgeForFace(const Edge *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackFaceFromEdge(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalEdgeToFace(const Edge *sender, Face *receiver) const override;

  void packFaceForEdge(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackEdgeFromFace(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalFaceToEdge(const Face *sender, Edge *receiver) const override;

  void packFaceForVertex(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackVertexFromFace(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalFaceToVertex(const Face *sender, Vertex *receiver) const override;

  void packFaceForCell(const Face *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackCellFromFace(Cell *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalFaceToCell(const Face *sender, Cell *receiver) const override;

  void packCellForFace(const Cell *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackFaceFromCell(Face *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalCellToFace(const Cell *sender, Face *receiver) const override;

  void packCellForEdge(const Cell *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackEdgeFromCell(Edge *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalCellToEdge(const Cell *sender, Edge *receiver) const override;

  void packCellForVertex(const Cell *sender, const PrimitiveID &receiver, walberla::mpi::SendBuffer &buffer) const override;

  void unpackVertexFromCell(Vertex *receiver, const PrimitiveID &sender, walberla::mpi::RecvBuffer &buffer) const override;

  void communicateLocalCellToVertex(const Cell *sender, Vertex *receiver) const override;

 private:
   using communication::DoFSpacePackInfo< ValueType >::level_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDVertex_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDEdge_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDFace_;
   using communication::DoFSpacePackInfo< ValueType >::dataIDCell_;
   using communication::DoFSpacePackInfo< ValueType >::storage_;

   BoundaryCondition boundaryCondition_;
   DoFType           boundaryTypeToSkip_;
};

} //namespace hhg
