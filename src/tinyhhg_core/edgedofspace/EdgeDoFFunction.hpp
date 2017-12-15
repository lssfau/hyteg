#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"

#include "tinyhhg_core/communication/BufferedCommunication.hpp"

#include "EdgeDoFMemory.hpp"
#include "EdgeDoFDataHandling.hpp"
#include "EdgeDoFMacroFace.hpp"
#include "EdgeDoFMacroEdge.hpp"
#include "EdgeDoFPackInfo.hpp"


namespace hhg {

template< typename ValueType >
class EdgeDoFFunction : public Function< EdgeDoFFunction< ValueType > >
{
public:

  EdgeDoFFunction( const std::string & name, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & minLevel, const uint_t & maxLevel ) :
      Function< EdgeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
  {
    auto edgeDoFMacroVertexFunctionMemoryDataHandling = std::make_shared< EdgeDoFMacroVertexFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
    auto edgeDoFMacroEdgeFunctionMemoryDataHandling   = std::make_shared< EdgeDoFMacroEdgeFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
    auto edgeDoFMacroFaceFunctionMemoryDataHandling   = std::make_shared< EdgeDoFMacroFaceFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );

    storage->addVertexData( vertexDataID_, edgeDoFMacroVertexFunctionMemoryDataHandling, name );
    storage->addEdgeData(   edgeDataID_,   edgeDoFMacroEdgeFunctionMemoryDataHandling,   name );
    storage->addFaceData(   faceDataID_,   edgeDoFMacroFaceFunctionMemoryDataHandling,   name );

    for (uint_t level = minLevel; level <= maxLevel; ++level) {
      //communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
      communicators_[level]->addPackInfo(
        std::make_shared<EdgeDoFPackInfo<ValueType> >(level, vertexDataID_, edgeDataID_, faceDataID_, storage_));
    }
  }

  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex>   & getVertexDataID() const { return vertexDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Edge>   & getEdgeDataID()   const { return edgeDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Face>   & getFaceDataID()   const { return faceDataID_; }

private:

    using Function< EdgeDoFFunction< ValueType > >::storage_;
    using Function< EdgeDoFFunction< ValueType > >::communicators_;

    /// Interpolates a given expression to a EdgeDoFFunction
    inline void
    interpolate_impl(std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                                             const std::vector<EdgeDoFFunction<ValueType>*> srcFunctions,
                                             uint_t level,
                                             DoFType flag = All);

    inline void
    assign_impl( const std::vector< ValueType > scalars, const std::vector< EdgeDoFFunction< ValueType >* > functions,
                 uint_t level, DoFType flag = All );

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
  WALBERLA_LOG_WARNING_ON_ROOT( "Interpolate not fully implemented!" );

  // Collect all source IDs in a vector
//  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> srcVertexIDs;
//  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>>   srcEdgeIDs;
  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Face>>   srcFaceIDs;

  for (auto& function : srcFunctions)
  {
//    srcVertexIDs.push_back(function->vertexDataID_);
//    srcEdgeIDs.push_back(function->edgeDataID_);
    srcFaceIDs.push_back(function->faceDataID_);
  }

  for ( auto & it : storage_->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( face.type, flag ) )
    {
      edgedof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
    }
  }
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::assign_impl(const std::vector<ValueType> scalars, const std::vector<EdgeDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
  WALBERLA_LOG_WARNING_ON_ROOT( "Assign not fully implemented!" );

  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Vertex >> srcVertexIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge >>     srcEdgeIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face >>     srcFaceIDs;

  for ( auto& function : functions )
  {
      srcVertexIDs.push_back( function->vertexDataID_ );
      srcEdgeIDs.push_back( function->edgeDataID_);
      srcFaceIDs.push_back( function->faceDataID_ );
  }

  for ( auto & it : storage_->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( face.type, flag ) )
    {
      edgedof::macroface::assign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
    }
  }
}

template< typename ValueType >
inline void EdgeDoFFunction< ValueType >::add_impl(const std::vector<ValueType> scalars, const std::vector<EdgeDoFFunction< ValueType >*> functions, size_t level, DoFType flag)
{
  WALBERLA_LOG_WARNING_ON_ROOT( "Add not fully implemented!" );

  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Vertex >> srcVertexIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Edge >>     srcEdgeIDs;
  std::vector<PrimitiveDataID< FunctionMemory< ValueType >, Face >>     srcFaceIDs;

  for ( auto& function : functions )
  {
      srcVertexIDs.push_back( function->vertexDataID_ );
      srcEdgeIDs.push_back( function->edgeDataID_);
      srcFaceIDs.push_back( function->faceDataID_ );
  }

  for ( auto & it : storage_->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( face.type, flag ) )
    {
      edgedof::macroface::add< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
    }
  }
}

template< typename ValueType >
inline real_t EdgeDoFFunction< ValueType >::dot_impl(EdgeDoFFunction< ValueType >& rhs, size_t level, DoFType flag)
{
  WALBERLA_LOG_WARNING_ON_ROOT( "Add not fully implemented!" );

  real_t scalarProduct =  0.0 ;

  for ( auto & it : storage_->getFaces() )
  {
    Face & face = *it.second;

    if ( testFlag( face.type, flag ) )
    {
      scalarProduct += edgedof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
    }
  }

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
  for (auto& it : storage_->getEdges()) {
    Edge& edge = *it.second;
    edgedof::macroedge::enumerate< ValueType >(level, edge, edgeDataID_, num);
  }

  communicators_[level]->template startCommunication<Edge, Face>();


  for (auto& it : storage_->getFaces()) {
    Face& face = *it.second;
    edgedof::macroface::enumerate< ValueType >(level, face, faceDataID_, num);
  }

  communicators_[level]->template endCommunication<Edge, Face>();

  communicators_[level]->template startCommunication<Face, Edge>();
  communicators_[level]->template endCommunication<Face, Edge>();

  communicators_[level]->template startCommunication<Edge, Vertex>();
  communicators_[level]->template endCommunication<Edge, Vertex>();

}

}
