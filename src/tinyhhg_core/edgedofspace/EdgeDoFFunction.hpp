#pragma once

#include "core/mpi/all.h"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include "tinyhhg_core/edgedofspace/EdgeDoFMacroCell.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroFace.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFPackInfo.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/boundary/BoundaryConditions.hpp"

namespace hhg {

using walberla::real_c;

namespace edgedof {
///@name Size Functions
///@{

inline uint_t edgeDoFMacroVertexFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  return primitive.getNumNeighborEdges() + primitive.getNumNeighborFaces();
}

inline uint_t edgeDoFMacroEdgeFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  return levelinfo::num_microedges_per_edge( level ) + primitive.getNumNeighborFaces() * ( 3 * ( levelinfo::num_microedges_per_edge( level ) ) - 1 )
    + primitive.getNumNeighborCells() * levelinfo::num_microedges_per_edge( level );
}

inline uint_t edgeDoFMacroFaceFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( primitive );
  ///"inner/own" points on the face
  uint_t innerDofs = 3 * ( ( ( levelinfo::num_microedges_per_edge( level ) + 1 ) * levelinfo::num_microedges_per_edge( level ) ) / 2 );

  ///ghost points on one adjacent tet
  uint_t GhostDoFsOneSide = 0;
  if(primitive.getNumNeighborCells() != 0){
    /// points in the "white up" tets
    GhostDoFsOneSide += 3 * levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microedges_per_edge( level ) );
    /// points from the xyz edge
    GhostDoFsOneSide += levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microedges_per_edge( level ) - 1 );
    /// points on the parallel face inside the tet
    GhostDoFsOneSide += 3 * levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microedges_per_edge( level ) - 1 );
  }

  return innerDofs + primitive.getNumNeighborCells() * GhostDoFsOneSide;
}

inline uint_t edgeDoFMacroCellFunctionMemorySize( const uint_t & level, const Primitive & primitive )
{
  WALBERLA_UNUSED( primitive );
  return   6 * ( levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) )
         + ( levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) - 1 ) );
}

///@}

}// namespace edgedof

template< typename ValueType >
class EdgeDoFFunction : public Function< EdgeDoFFunction< ValueType > >
{
public:

  EdgeDoFFunction( const std::string & name, const std::shared_ptr< PrimitiveStorage > & storage ) :
    Function< EdgeDoFFunction< ValueType > >( name, storage ),
    vertexDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Vertex >, Vertex >() ),
    edgeDataID_(   storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Edge >,   Edge >() ),
    faceDataID_(   storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Face >,   Face >() )
  {}

  EdgeDoFFunction( const std::string & name, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & minLevel, const uint_t & maxLevel ) :
    EdgeDoFFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create012BC() )
  {}

  EdgeDoFFunction( const std::string & name, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & minLevel, const uint_t & maxLevel, const BoundaryCondition & boundaryCondition ) :
      Function< EdgeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel ), boundaryCondition_( boundaryCondition )
  {
    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Vertex >> vertexDataHandling =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Vertex >>(minLevel, maxLevel, edgedof::edgeDoFMacroVertexFunctionMemorySize);

    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Edge >> edgeDataHandling   =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Edge   >>(minLevel, maxLevel, edgedof::edgeDoFMacroEdgeFunctionMemorySize);

    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Face >> faceDataHandling   =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Face   >>(minLevel, maxLevel, edgedof::edgeDoFMacroFaceFunctionMemorySize);

    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Cell >> cellDataHandling   =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Cell   >>(minLevel, maxLevel, edgedof::edgeDoFMacroCellFunctionMemorySize);


    storage->addVertexData( vertexDataID_, vertexDataHandling, name );
    storage->addEdgeData(   edgeDataID_,   edgeDataHandling,   name );
    storage->addFaceData(   faceDataID_,   faceDataHandling,   name );
    storage->addCellData(   cellDataID_,   cellDataHandling,   name );

    for (uint_t level = minLevel; level <= maxLevel; ++level) {
      //communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
      communicators_[level]->addPackInfo(
        std::make_shared<EdgeDoFPackInfo<ValueType> >(level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage()));
    }
  }

  inline void
  assign( const std::vector< ValueType > scalars, const std::vector< EdgeDoFFunction< ValueType >* > functions,
          uint_t level, DoFType flag = All );

  inline void
  add( const std::vector< ValueType > scalars, const std::vector< EdgeDoFFunction< ValueType >* > functions,
       uint_t level, DoFType flag = All );

  /// Interpolates a given expression to a EdgeDoFFunction

  inline void
  interpolate( const ValueType& constant, uint_t level, DoFType flag = All );

  inline void
  interpolate( const std::function< ValueType( const Point3D & ) >& expr,
                          uint_t level, DoFType flag = All);

  inline void
  interpolateExtended( const std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                       const std::vector<EdgeDoFFunction<ValueType>*> srcFunctions,
                      uint_t level,
                      DoFType flag = All);

  inline real_t
  dotLocal( EdgeDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All );

  inline void enumerate( uint_t level );

  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex>   & getVertexDataID() const { return vertexDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Edge>   & getEdgeDataID()   const { return edgeDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Face>   & getFaceDataID()   const { return faceDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Cell>   & getCellDataID()   const { return cellDataID_; }


  inline ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true );

  inline BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }

  template< typename SenderType, typename ReceiverType >
  inline void startCommunication( const uint_t & level ) const
  {
    if ( isDummy() ) { return; }
    communicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
  }

  template< typename SenderType, typename ReceiverType >
  inline void endCommunication( const uint_t & level ) const
  {
    if ( isDummy() ) { return; }
    communicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
  }

  template< typename SenderType, typename ReceiverType >
  inline void communicate( const uint_t & level ) const
  {
    if ( isDummy() ) { return; }
    communicators_.at( level )->template communicate< SenderType, ReceiverType >();
  }

  inline void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode )
  {
    if ( isDummy() ) { return; }
    for ( auto & communicator : communicators_ )
    {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
    }
  }

   using Function< EdgeDoFFunction< ValueType > >::isDummy;

private:

   inline void enumerate( uint_t level, ValueType& offset );

   using Function< EdgeDoFFunction< ValueType > >::communicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge > edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face > faceDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell > cellDataID_;

   BoundaryCondition boundaryCondition_;

   /// friend P2Function for usage of enumerate
   friend class P2Function< ValueType >;
};


}// namespace hhg
