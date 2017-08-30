
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

#include <algorithm>
#include <iomanip>
#include <set>

namespace hhg {

SetupPrimitiveStorage::SetupPrimitiveStorage( const MeshInfo & meshInfo, const uint_t & numberOfProcesses ) :
    numberOfProcesses_( numberOfProcesses )
{
  WALBERLA_ASSERT_GREATER( numberOfProcesses_, 0, "Number of processes must be positive" );

  // since the MeshInfo IDs of the vertices do not necessarily
  // match the primitive IDs of the vertices in the SetupStorage, we need an assignment
  std::map< uint_t, PrimitiveID > meshVertexIDToPrimitiveID;

  // Adding vertices to storage
  const MeshInfo::VertexContainer vertices = meshInfo.getVertices();
  for ( const auto & it : vertices )
  {
    const MeshInfo::Vertex meshInfoVertex = it.second;

    PrimitiveID vertexID = generatePrimitiveID();

    meshVertexIDToPrimitiveID[ meshInfoVertex.getID() ] = vertexID;

    Point3D coordinates( meshInfoVertex.getCoordinates() );
    vertices_[ vertexID.getID() ] = std::make_shared< Vertex >( vertexID, coordinates );

    // All to root by default
    primitiveIDToTargetRankMap_[ vertexID.getID() ] = 0;
  }

  // Adding edges to storage
  const MeshInfo::EdgeContainer edges = meshInfo.getEdges();
  for ( const auto & it : edges )
  {
    const MeshInfo::Edge meshInfoEdge = it.second;

    PrimitiveID edgeID = generatePrimitiveID();

    WALBERLA_ASSERT_EQUAL( meshInfoEdge.getVertices().size(), 2, "Edges are expected to have two vertices." );
    PrimitiveID vertexID0 = meshVertexIDToPrimitiveID[ meshInfoEdge.getVertices().at( 0 )  ];
    PrimitiveID vertexID1 = meshVertexIDToPrimitiveID[ meshInfoEdge.getVertices().at( 1 ) ];

    DoFType dofType = meshInfoEdge.getDoFType();

    std::array<Point3D, 2> coords;

    coords[0] = vertices_[ vertexID0.getID() ]->getCoordinates();
    coords[1] = vertices_[ vertexID1.getID() ]->getCoordinates();

    WALBERLA_ASSERT_EQUAL( edges_.count( edgeID.getID() ), 0 );
    WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID0.getID() ), 1 );
    WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID1.getID() ), 1 );
    edges_[ edgeID.getID() ] = std::make_shared< Edge >( edgeID, vertexID0, vertexID1, dofType, coords);

    // All to root by default
    primitiveIDToTargetRankMap_[ edgeID.getID() ] = 0;

    // Adding edge ID as neighbor to SetupVertices
    vertices_[ vertexID0.getID() ]->addEdge( edgeID );
    vertices_[ vertexID1.getID() ]->addEdge( edgeID );
  }

  // Adding faces to storage
  const MeshInfo::FaceContainer faces = meshInfo.getFaces();
  for ( const auto & it : faces )
  {
    const MeshInfo::Face meshInfoFace = it.second;

    PrimitiveID faceID = generatePrimitiveID();

    WALBERLA_ASSERT_EQUAL( meshInfoFace.getVertices().size(), 3, "Only supporting triangle faces." );
    PrimitiveID vertexID0 = meshVertexIDToPrimitiveID[ meshInfoFace.getVertices().at( 0 ) ];
    PrimitiveID vertexID1 = meshVertexIDToPrimitiveID[ meshInfoFace.getVertices().at( 1 ) ];
    PrimitiveID vertexID2 = meshVertexIDToPrimitiveID[ meshInfoFace.getVertices().at( 2 ) ];

    WALBERLA_ASSERT_EQUAL( faces_.count( faceID.getID() ), 0 );
    WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID0.getID() ), 1 );
    WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID1.getID() ), 1 );
    WALBERLA_ASSERT_EQUAL( vertices_.count( vertexID2.getID() ), 1 );

    PrimitiveID edgeID0;
    PrimitiveID edgeID1;
    PrimitiveID edgeID2;

    const bool foundEdge0 = findEdgeByVertexIDs( vertexID0, vertexID1, edgeID0 );
    const bool foundEdge1 = findEdgeByVertexIDs( vertexID1, vertexID2, edgeID1 );
    const bool foundEdge2 = findEdgeByVertexIDs( vertexID2, vertexID0, edgeID2 );

    WALBERLA_CHECK( foundEdge0 && foundEdge1 && foundEdge2, "Could not successfully construct faces from MeshInfo" );

    WALBERLA_ASSERT_EQUAL( edges_.count( edgeID0.getID() ), 1 );
    WALBERLA_ASSERT_EQUAL( edges_.count( edgeID1.getID() ), 1 );
    WALBERLA_ASSERT_EQUAL( edges_.count( edgeID2.getID() ), 1 );

    // Edge Orientation
    std::array< int, 3 > edgeOrientation;

    PrimitiveID edge0Vertex0 = edges_[ edgeID0.getID() ]->getVertexID0();
    PrimitiveID edge0Vertex1 = edges_[ edgeID0.getID() ]->getVertexID1();
    PrimitiveID edge1Vertex0 = edges_[ edgeID1.getID() ]->getVertexID0();
    PrimitiveID edge1Vertex1 = edges_[ edgeID1.getID() ]->getVertexID1();
    PrimitiveID edge2Vertex0 = edges_[ edgeID2.getID() ]->getVertexID0();
    PrimitiveID edge2Vertex1 = edges_[ edgeID2.getID() ]->getVertexID1();

    if (edge0Vertex1 == edge1Vertex0 && edge1Vertex1 == edge2Vertex0 && edge2Vertex1 == edge0Vertex0)
    {
      edgeOrientation = {{1, 1, 1}};
    }
    else if (edge0Vertex1 == edge1Vertex0 && edge1Vertex1 == edge2Vertex1 && edge2Vertex0 == edge0Vertex0)
    {
      edgeOrientation = {{1, 1, -1}};
    }
    else if (edge0Vertex1 == edge1Vertex1 && edge1Vertex0 == edge2Vertex0 && edge2Vertex1 == edge0Vertex0)
    {
      edgeOrientation = {{1, -1, 1}};
    }
    else if (edge0Vertex1 == edge1Vertex1 && edge1Vertex0 == edge2Vertex1 && edge2Vertex0 == edge0Vertex0)
    {
      edgeOrientation = {{1, -1, -1}};
    }
    else if (edge0Vertex0 == edge1Vertex0 && edge1Vertex1 == edge2Vertex0 && edge2Vertex1 == edge0Vertex1)
    {
      edgeOrientation = {{-1, 1, 1}};
    }
    else if (edge0Vertex0 == edge1Vertex0 && edge1Vertex1 == edge2Vertex1 && edge2Vertex0 == edge0Vertex1)
    {
      edgeOrientation = {{-1, 1, -1}};
    }
    else if (edge0Vertex0 == edge1Vertex1 && edge1Vertex0 == edge2Vertex0 && edge2Vertex1 == edge0Vertex1)
    {
      edgeOrientation = {{-1, -1, 1}};
    }
    else if (edge0Vertex0 == edge1Vertex1 && edge1Vertex0 == edge2Vertex1 && edge2Vertex0 == edge0Vertex1)
    {
      edgeOrientation = {{-1, -1, -1}};
    }

    // Corner coordinates
    std::array< Point3D, 3 > coordinates;
    std::array< PrimitiveID, 3 > vertexIDs;

    if (edgeOrientation[0] == 1)
    {
      coordinates[0] = vertices_[ edge0Vertex0.getID() ]->getCoordinates();
      coordinates[1] = vertices_[ edge0Vertex1.getID() ]->getCoordinates();

      vertexIDs[0] = edge0Vertex0.getID();
      vertexIDs[1] = edge0Vertex1.getID();
    }
    else
    {
      coordinates[0] = vertices_[ edge0Vertex1.getID() ]->getCoordinates();
      coordinates[1] = vertices_[ edge0Vertex0.getID() ]->getCoordinates();

      vertexIDs[0] = edge0Vertex1.getID();
      vertexIDs[1] = edge0Vertex0.getID();
    }

    if (edgeOrientation[1] == 1)
    {
      coordinates[2] = vertices_[ edge1Vertex1.getID() ]->getCoordinates();

      vertexIDs[2] = edge1Vertex1.getID();
    }
    else
    {
      coordinates[2] = vertices_[ edge1Vertex0.getID() ]->getCoordinates();

      vertexIDs[2] = edge1Vertex0.getID();
    }

    faces_[ faceID.getID() ] = std::shared_ptr< Face >( new Face( faceID, vertexIDs, {{edgeID0, edgeID1, edgeID2}}, edgeOrientation, coordinates ) );

    // All to root by default
    primitiveIDToTargetRankMap_[ faceID.getID() ] = 0;

    // Adding face ID to vertices as neighbors
    vertices_[vertexIDs[0].getID()]->addFace(faceID);
    vertices_[vertexIDs[1].getID()]->addFace(faceID);
    vertices_[vertexIDs[2].getID()]->addFace(faceID);

    // Adding face ID to edges as neighbors
    edges_[ edgeID0.getID() ]->addFace( faceID );
    edges_[ edgeID1.getID() ]->addFace( faceID );
    edges_[ edgeID2.getID() ]->addFace( faceID );
  }

  const MeshInfo::CellContainer cells = meshInfo.getCells();
  for ( const auto & it : cells )
  {
    const MeshInfo::Cell meshInfoCell = it.second;

    PrimitiveID cellID = generatePrimitiveID();

    WALBERLA_ASSERT_EQUAL( meshInfoCell.getVertices().size(), 4, "Only supporting tetrahedron cells." );

    PrimitiveID vertexID0 = meshVertexIDToPrimitiveID[ meshInfoCell.getVertices().at( 0 ) ];
    PrimitiveID vertexID1 = meshVertexIDToPrimitiveID[ meshInfoCell.getVertices().at( 1 ) ];
    PrimitiveID vertexID2 = meshVertexIDToPrimitiveID[ meshInfoCell.getVertices().at( 2 ) ];
    PrimitiveID vertexID3 = meshVertexIDToPrimitiveID[ meshInfoCell.getVertices().at( 3 ) ];

    PrimitiveID edgeID0;
    PrimitiveID edgeID1;
    PrimitiveID edgeID2;
    PrimitiveID edgeID3;
    PrimitiveID edgeID4;
    PrimitiveID edgeID5;

    PrimitiveID faceID0;
    PrimitiveID faceID1;
    PrimitiveID faceID2;
    PrimitiveID faceID3;

    const bool foundEdge0 = findEdgeByVertexIDs( vertexID0, vertexID1, edgeID0 );
    const bool foundEdge1 = findEdgeByVertexIDs( vertexID0, vertexID2, edgeID1 );
    const bool foundEdge2 = findEdgeByVertexIDs( vertexID0, vertexID3, edgeID2 );
    const bool foundEdge3 = findEdgeByVertexIDs( vertexID1, vertexID2, edgeID3 );
    const bool foundEdge4 = findEdgeByVertexIDs( vertexID1, vertexID3, edgeID4 );
    const bool foundEdge5 = findEdgeByVertexIDs( vertexID2, vertexID3, edgeID5 );

    const bool foundFace0 = findFaceByVertexIDs( vertexID0, vertexID1, vertexID2, faceID0 );
    const bool foundFace1 = findFaceByVertexIDs( vertexID0, vertexID1, vertexID3, faceID1 );
    const bool foundFace2 = findFaceByVertexIDs( vertexID0, vertexID2, vertexID3, faceID2 );
    const bool foundFace3 = findFaceByVertexIDs( vertexID1, vertexID2, vertexID3, faceID3 );

    WALBERLA_CHECK( foundEdge0 && foundEdge1 && foundEdge2 && foundEdge3 && foundEdge4 && foundEdge5, "Could not successfully construct cell from MeshInfo." );
    WALBERLA_CHECK( foundFace0 && foundFace1 && foundFace2 && foundFace3,                             "Could not successfully construct cell from MeshInfo." );

    std::vector< PrimitiveID > cellVertices = {{ vertexID0, vertexID1, vertexID2, vertexID3 }};
    std::vector< PrimitiveID > cellEdges    = {{ edgeID0, edgeID1, edgeID2, edgeID3, edgeID4, edgeID5 }};
    std::vector< PrimitiveID > cellFaces    = {{ faceID0, faceID1, faceID2, faceID3 }};

    for ( const auto & id : cellVertices ) { WALBERLA_ASSERT( vertexExists( id ) ); vertices_[ id.getID() ]->addCell( cellID ); }
    for ( const auto & id : cellEdges    ) { WALBERLA_ASSERT(   edgeExists( id ) );    edges_[ id.getID() ]->addCell( cellID ); }
    for ( const auto & id : cellFaces    ) { WALBERLA_ASSERT(   faceExists( id ) );    faces_[ id.getID() ]->addCell( cellID ); }

    cells_[ cellID.getID() ] = std::make_shared< Cell >( cellID, cellVertices, cellEdges, cellFaces );
  }

  for (auto& it : edges_) {
    Edge& edge = *it.second;

    if (testFlag(edge.dofType_, hhg::NeumannBoundary)) {

      for (auto& itv : edge.neighborVertices()) {
        vertices_[itv.getID()]->dofType_ = hhg::NeumannBoundary;
      }

    }
  }

  for (auto& it : edges_) {
    Edge& edge = *it.second;

    if (testFlag(edge.dofType_, hhg::DirichletBoundary)) {

      for (auto& itv : edge.neighborVertices()) {
        vertices_[itv.getID()]->dofType_ = hhg::DirichletBoundary;
      }

    }
  }
}

const Primitive * SetupPrimitiveStorage::getPrimitive( const PrimitiveID & id ) const
{
  if ( vertexExists( id ) ) { return getVertex( id ); }
  if ( edgeExists( id ) )   { return getEdge( id ); }
  if ( faceExists( id ) )   { return getFace( id ); }
  if ( cellExists( id ) )   { return getCell( id ); }
  return nullptr;
}

bool SetupPrimitiveStorage::findEdgeByVertexIDs( const PrimitiveID & vertexID0, const PrimitiveID & vertexID1, PrimitiveID & edge ) const
{
  if ( vertices_.count( vertexID0.getID() ) == 0 || vertices_.count( vertexID1.getID() ) == 0 )
  {
    return false;
  }

  for ( auto it = edges_.begin(); it != edges_.end(); it++ )
  {
    if (   ( it->second->getVertexID0() == vertexID0 && it->second->getVertexID1() == vertexID1 )
        || ( it->second->getVertexID0() == vertexID1 && it->second->getVertexID1() == vertexID0 ) )
    {
      edge = PrimitiveID( it->first );
      return true;
    }
  }

  return false;
}

bool SetupPrimitiveStorage::findFaceByVertexIDs( const PrimitiveID & vertexID0, const PrimitiveID & vertexID1, const PrimitiveID & vertexID2, PrimitiveID & faceID ) const
{
  if ( !vertexExists( vertexID0 ) || !vertexExists( vertexID1 ) || !vertexExists( vertexID2 ) )
  {
    return false;
  }

  for ( const auto & it : faces_ )
  {
    auto face = it.second;
    if (   ( face->getVertexID0() == vertexID0 && face->getVertexID1() == vertexID1 && face->getVertexID2() == vertexID2 )
        || ( face->getVertexID0() == vertexID0 && face->getVertexID1() == vertexID2 && face->getVertexID2() == vertexID1 )
        || ( face->getVertexID0() == vertexID1 && face->getVertexID1() == vertexID0 && face->getVertexID2() == vertexID2 )
        || ( face->getVertexID0() == vertexID1 && face->getVertexID1() == vertexID2 && face->getVertexID2() == vertexID0 )
        || ( face->getVertexID0() == vertexID2 && face->getVertexID1() == vertexID0 && face->getVertexID2() == vertexID1 )
        || ( face->getVertexID0() == vertexID2 && face->getVertexID1() == vertexID1 && face->getVertexID2() == vertexID0 )
       )
    {
      faceID = PrimitiveID( it.first );
      return true;
    }
  }

  return false;
}


void SetupPrimitiveStorage::assembleRankToSetupPrimitivesMap( RankToSetupPrimitivesMap & rankToSetupPrimitivesMap ) const
{
  rankToSetupPrimitivesMap.clear();

  PrimitiveMap setupPrimitives;
  getSetupPrimitives( setupPrimitives );
  for ( uint_t rank = 0; rank < numberOfProcesses_; rank++ )
  {
    rankToSetupPrimitivesMap[ rank ] = std::vector< PrimitiveID::IDType >();
    for ( auto setupPrimitive : setupPrimitives )
    {
      if ( rank == getTargetRank( setupPrimitive.first ) )
      {
	      rankToSetupPrimitivesMap[ rank ].push_back( setupPrimitive.first );
      }
    }
  }

  WALBERLA_ASSERT_LESS_EQUAL( rankToSetupPrimitivesMap.size(), numberOfProcesses_ );
}

uint_t SetupPrimitiveStorage::getNumberOfEmptyProcesses() const
{
  uint_t numberOfEmptyProcesses = 0;
  RankToSetupPrimitivesMap rankToSetupPrimitivesMap;
  assembleRankToSetupPrimitivesMap( rankToSetupPrimitivesMap );
  for ( auto const & rankToSetupPrimitives : rankToSetupPrimitivesMap )
  {
    if ( rankToSetupPrimitives.second.size() == 0 )
    {
      numberOfEmptyProcesses++;
    }
  }
  return numberOfEmptyProcesses;
}

uint_t SetupPrimitiveStorage::getNumberOfPrimitives() const
{
  return vertices_.size() + edges_.size() + faces_.size() + cells_.size();
}

uint_t SetupPrimitiveStorage::getMinPrimitivesPerRank() const
{
  uint_t minNumberOfPrimitives = std::numeric_limits< uint_t >::max();
  RankToSetupPrimitivesMap rankToSetupPrimitivesMap;
  assembleRankToSetupPrimitivesMap( rankToSetupPrimitivesMap );
  for ( auto const & rankToSetupPrimitives : rankToSetupPrimitivesMap )
  {
    minNumberOfPrimitives = std::min( rankToSetupPrimitives.second.size(), minNumberOfPrimitives );
  }
  return minNumberOfPrimitives;
}

uint_t SetupPrimitiveStorage::getMaxPrimitivesPerRank() const
{
  uint_t maxNumberOfPrimitives = 0;
  RankToSetupPrimitivesMap rankToSetupPrimitivesMap;
  assembleRankToSetupPrimitivesMap( rankToSetupPrimitivesMap );
  for ( auto const & rankToSetupPrimitives : rankToSetupPrimitivesMap )
  {
    maxNumberOfPrimitives = std::max( rankToSetupPrimitives.second.size(), maxNumberOfPrimitives );
  }
  return maxNumberOfPrimitives;
}

real_t SetupPrimitiveStorage::getAvgPrimitivesPerRank() const
{
  return real_t( getNumberOfPrimitives() ) / real_t( numberOfProcesses_ );
}


void SetupPrimitiveStorage::toStream( std::ostream & os ) const
{
  os << "SetupPrimitiveStorage:\n";

  os << " - Processes (overall): " << std::setw(10) << numberOfProcesses_ << "\n";
  os << " - Processes (empty)  : " << std::setw(10) << getNumberOfEmptyProcesses() << "\n";

  os << " - Number of...\n"
     << "   +  Vertices: " << std::setw(10) << vertices_.size() << "\n"
     << "   +     Edges: " << std::setw(10) << edges_.size() << "\n"
     << "   +     Faces: " << std::setw(10) << faces_.size() << "\n";

  os << " - Primitives per process...\n"
     << "   +      min: " << std::setw(10) << getMinPrimitivesPerRank() << "\n"
     << "   +      max: " << std::setw(10) << getMaxPrimitivesPerRank() << "\n"
     << "   +      avg: " << std::setw(10) << getAvgPrimitivesPerRank() << "\n";

#ifndef NDEBUG
  os << "\n";
  os << "Vertices:   ID | Target Rank | Position  | Neighbor Edges \n"
     << "---------------------------------------------------------\n";
  for ( auto it = vertices_.begin(); it != vertices_.end(); it++ )
  {
    Point3D coordinates = it->second->getCoordinates();
    os << "          " << std::setw(4) << it->first << " | "
       << std::setw(11) << getTargetRank( it->first ) << " | "
       << coordinates << " | ";
    for ( const auto & neighborEdgeID : it->second->getHigherDimNeighbors() )
    {
      os << neighborEdgeID.getID() << " ";
    }
    os << "\n";

  }
  os << "\n";

  os << "Edges:      ID | Target Rank | VertexID_0 | VertexID_1 | DoF Type             | Neighbor Faces \n"
     << "----------------------------------------------------------------------------------------------\n";
  for ( auto it = edges_.begin(); it != edges_.end(); it++ )
  {
    os << "          " << std::setw(4) << it->first << " | "
       << std::setw(11) << getTargetRank( it->first ) << " | "
       << std::setw(10) << it->second->getVertexID0().getID() << " | "
       << std::setw(10) << it->second->getVertexID1().getID() << " | "
       << std::setw(20) << it->second->getDoFType() << " | ";
    for ( const auto & neighborFaceID : it->second->getHigherDimNeighbors() )
    {
      os << neighborFaceID.getID() << " ";
    }
        os << "\n";
  }
  os << "\n";

  os << "Faces:      ID | Target Rank | EdgeID_0 | EdgeID_1 | EdgeID_2\n"
     << "-------------------------------------------------------------\n";
  for ( auto it = faces_.begin(); it != faces_.end(); it++ )
  {
    os << "          " << std::setw(4) << it->first << " | "
       << std::setw(11) << getTargetRank( it->first ) << " | "
       << std::setw(8) << it->second->getEdgeID0().getID() << " | "
       << std::setw(8) << it->second->getEdgeID1().getID() << " | "
       << std::setw(8) << it->second->getEdgeID2().getID() << "\n";
  }
#endif
}


void SetupPrimitiveStorage::getSetupPrimitives( PrimitiveMap & setupPrimitiveMap ) const
{
  setupPrimitiveMap.clear();

  setupPrimitiveMap.insert( beginVertices(), endVertices() );
  setupPrimitiveMap.insert( beginEdges(), endEdges() );
  setupPrimitiveMap.insert( beginFaces(), endFaces() );
  setupPrimitiveMap.insert( beginCells(), endCells() );

  WALBERLA_ASSERT_EQUAL( setupPrimitiveMap.size(), vertices_.size() + edges_.size() + faces_.size() + cells_.size() );
}


PrimitiveID SetupPrimitiveStorage::generatePrimitiveID() const
{
  PrimitiveID newID( getNumberOfPrimitives() );
  WALBERLA_ASSERT( !primitiveExists( newID ) );
  return newID;
}

}
