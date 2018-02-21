#pragma once

#include "core/mpi/all.h"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include "EdgeDoFMacroFace.hpp"
#include "EdgeDoFMacroEdge.hpp"
#include "EdgeDoFPackInfo.hpp"


namespace hhg {

namespace edgedof {
///@name Size Functions
///@{

inline uint_t edgeDoFMacroVertexFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( level );
  return 2 * primitive.getNumNeighborEdges();
}

inline uint_t edgeDoFMacroEdgeFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  return levelinfo::num_microedges_per_edge( level ) + primitive.getNumNeighborFaces() * ( 3 * ( levelinfo::num_microedges_per_edge( level ) ) - 1 );
}

inline uint_t edgeDoFMacroFaceFunctionMemorySize( const uint_t &level, const Primitive & primitive )
{
  WALBERLA_UNUSED( primitive );
  return 3 * ( ( ( levelinfo::num_microedges_per_edge( level ) + 1 ) * levelinfo::num_microedges_per_edge( level ) ) / 2 );
}

///@}

}// namespace edgedof

template< typename ValueType >
class EdgeDoFFunction : public Function< EdgeDoFFunction< ValueType > >
{
public:

  EdgeDoFFunction( const std::string & name, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & minLevel, const uint_t & maxLevel ) :
      Function< EdgeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
  {
    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Vertex >> vertexDataHandling =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Vertex >>(minLevel, maxLevel, edgedof::edgeDoFMacroVertexFunctionMemorySize);

    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Edge >> edgeDataHandling   =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Edge   >>(minLevel, maxLevel, edgedof::edgeDoFMacroEdgeFunctionMemorySize);

    std::shared_ptr<MemoryDataHandling<FunctionMemory< ValueType >, Face >> faceDataHandling   =
        std::make_shared< MemoryDataHandling<FunctionMemory< ValueType >, Face   >>(minLevel, maxLevel, edgedof::edgeDoFMacroFaceFunctionMemorySize);


    storage->addVertexData( vertexDataID_, vertexDataHandling, name );
    storage->addEdgeData(   edgeDataID_,   edgeDataHandling,   name );
    storage->addFaceData(   faceDataID_,   faceDataHandling,   name );

    for (uint_t level = minLevel; level <= maxLevel; ++level) {
      //communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
      communicators_[level]->addPackInfo(
        std::make_shared<EdgeDoFPackInfo<ValueType> >(level, vertexDataID_, edgeDataID_, faceDataID_, this->getStorage()));
    }
  }

  inline void
  assign( const std::vector< ValueType > scalars, const std::vector< EdgeDoFFunction< ValueType >* > functions,
          uint_t level, DoFType flag = All );

  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex>   & getVertexDataID() const { return vertexDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Edge>   & getEdgeDataID()   const { return edgeDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Face>   & getFaceDataID()   const { return faceDataID_; }

private:

    using Function< EdgeDoFFunction< ValueType > >::communicators_;

    /// Interpolates a given expression to a EdgeDoFFunction
    inline void
    interpolate_impl(std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                                             const std::vector<EdgeDoFFunction<ValueType>*> srcFunctions,
                                             uint_t level,
                                             DoFType flag = All);

    inline void
    add_impl( const std::vector< ValueType > scalars, const std::vector< EdgeDoFFunction< ValueType >* > functions,
              uint_t level, DoFType flag = All );

    inline real_t
    dot_impl( EdgeDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All );

    inline void
    prolongate_impl( uint_t sourceLevel, DoFType flag = All );

    inline void
    prolongateQuadratic_impl( uint_t sourceLevel, DoFType flag = All );

    inline void
    restrict_impl( uint_t sourceLevel, DoFType flag = All );

    inline void
    enumerate_impl( uint_t level, uint_t& num );

    PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
    PrimitiveDataID< FunctionMemory< ValueType >, Edge > edgeDataID_;
    PrimitiveDataID< FunctionMemory< ValueType >, Face > faceDataID_;
};

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::interpolate_impl(std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                                                           const std::vector<EdgeDoFFunction<ValueType>*> srcFunctions,
                                                           uint_t level,
                                                           DoFType flag)
{
  // Collect all source IDs in a vector
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>>   srcEdgeIDs;
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>>   srcFaceIDs;

  for (auto& function : srcFunctions)
  {
    srcEdgeIDs.push_back(function->edgeDataID_);
    srcFaceIDs.push_back(function->faceDataID_);
  }

  for ( auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( edge.getDoFType(), flag ) )
    {
      edgedof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
    }
  }

  communicators_[ level ]->template startCommunication< Edge, Face >();

  for ( auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( face.type, flag ) )
    {
      edgedof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
    }
  }

  communicators_[ level ]->template endCommunication< Edge, Face >();
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::assign(const std::vector<ValueType> scalars, const std::vector<EdgeDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge >>     srcEdgeIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face >>     srcFaceIDs;

  for ( auto& function : functions )
  {
      srcEdgeIDs.push_back( function->edgeDataID_);
      srcFaceIDs.push_back( function->faceDataID_ );
  }

  for ( auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( edge.getDoFType(), flag ) )
    {
      edgedof::macroedge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
    }
  }

  communicators_[ level ]->template startCommunication< Edge, Face >();

  for ( auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( face.type, flag ) )
    {
      edgedof::macroface::assign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
    }
  }

  communicators_[ level ]->template endCommunication< Edge, Face >();

}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::add_impl(const std::vector<ValueType> scalars, const std::vector<EdgeDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge >>     srcEdgeIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face >>     srcFaceIDs;

  for ( auto& function : functions )
  {
      srcEdgeIDs.push_back( function->edgeDataID_);
      srcFaceIDs.push_back( function->faceDataID_ );
  }

  for ( auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( edge.getDoFType(), flag ) )
    {
      edgedof::macroedge::add< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
    }
  }

  communicators_[ level ]->template startCommunication< Edge, Face >();

  for ( auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( face.type, flag ) )
    {
      edgedof::macroface::add< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
    }
  }

  communicators_[ level ]->template endCommunication< Edge, Face >();
}

template< typename ValueType >
inline real_t EdgeDoFFunction< ValueType >::dot_impl(EdgeDoFFunction< ValueType >& rhs, size_t level, DoFType flag)
{
  real_t scalarProduct =  0.0 ;

  for ( auto & it : this->getStorage()->getEdges() )
  {
    Edge & edge = *it.second;

    if ( testFlag( edge.getDoFType(), flag ) )
    {
      scalarProduct += edgedof::macroedge::dot< ValueType >( level, edge, edgeDataID_, rhs.edgeDataID_ );
    }
  }

  for ( auto & it : this->getStorage()->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( face.type, flag ) )
    {
      scalarProduct += edgedof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
    }
  }

  walberla::mpi::allReduceInplace( scalarProduct, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
  return scalarProduct;
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::prolongate_impl(size_t sourceLevel, DoFType flag)
{
  WALBERLA_ASSERT( false, "To be implemented..." );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::prolongateQuadratic_impl(size_t sourceLevel, DoFType flag)
{
  WALBERLA_ASSERT( false, "To be implemented..." );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::restrict_impl(size_t sourceLevel, DoFType flag)
{
  WALBERLA_ASSERT( false, "To be implemented..." );
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::enumerate_impl(uint_t level, uint_t& num)
{
  for (auto& it : this->getStorage()->getEdges()) {
    Edge& edge = *it.second;
    edgedof::macroedge::enumerate< ValueType >(level, edge, edgeDataID_, num);
  }

  communicators_[level]->template startCommunication<Edge, Face>();


  for (auto& it : this->getStorage()->getFaces()) {
    Face& face = *it.second;
    edgedof::macroface::enumerate< ValueType >(level, face, faceDataID_, num);
  }

  communicators_[level]->template endCommunication<Edge, Face>();

  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();

}


}// namespace hhg
