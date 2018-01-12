#include "vtkwriter.hpp"
#include "levelinfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"

namespace hhg
{

using walberla::real_c;

static void writeXMLHeader( std::ostream & output )
{
  WALBERLA_ROOT_SECTION()
  {
    output << "<?xml version=\"1.0\"?>\n";
    output << "<VTKFile type=\"UnstructuredGrid\">\n";
    output << " <UnstructuredGrid>\n";
  }
}

static void writeXMLFooter( std::ostream & output )
{
  WALBERLA_ROOT_SECTION()
  {
    output << " </UnstructuredGrid>\n";
    output << "</VTKFile>\n";
  }
}

static void writePieceHeader( std::ostream & output, const uint_t & numberOfPoints, const uint_t & numberOfCells )
{
  output << "<Piece "
         << "NumberOfPoints=\"" << numberOfPoints << "\" "
         << "NumberOfCells=\"" << numberOfCells << "\""
         << ">\n";
}

static void writePieceFooter( std::ostream & output )
{
  output << "</Piece>\n";
}


const std::map< VTKOutput::DoFType, std::string > VTKOutput::DoFTypeToString_ =
{
  { DoFType::VERTEX,          "VertexDoF" },
  { DoFType::EDGE_HORIZONTAL, "HorizontalEdgeDoF" },
  { DoFType::EDGE_VERTICAL,   "VerticalEdgeDoF" },
  { DoFType::EDGE_DIAGONAL,   "DiagonalEdgeDoF" },
  { DoFType::DG,              "DGDoF" },
};


std::string VTKOutput::fileNameExtension( const VTKOutput::DoFType & dofType, const uint_t & level, const uint_t & timestep ) const
{
  return fmt::format( "_{}_level{}_ts{}", DoFTypeToString_.at( dofType ), std::to_string( level ), std::to_string( timestep ) );
}


void VTKOutput::writePointsForMicroVertices( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level ) const
{
  output << "<Points>\n";
  output << "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";

  for (auto& it : storage->getFaces()) {
    Face &face = *it.second;

    size_t rowsize = levelinfo::num_microvertices_per_edge( level );
    Point3D x, x0;

    x0 = face.coords[0];

    Point3D d0 = (face.coords[1] - face.coords[0]) / (real_c(rowsize)-1);
    Point3D d2 = (face.coords[2] - face.coords[0]) / (real_c(rowsize)-1);

    size_t inner_rowsize = rowsize;

    for (size_t i = 0; i < rowsize; ++i)
    {
      x = x0;
      x += real_c(i) * d2;

      for (size_t j = 0; j < inner_rowsize; ++j)
      {
        output << std::scientific << x[0] << " " << x[1] << " " << x[2] << " ";
        x += d0;
      }

      --inner_rowsize;
    }
  }

  output << "\n</DataArray>\n";
  output << "</Points>\n";
}

void VTKOutput::writePointsForMicroEdges( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage,
                                          const uint_t & level, const VTKOutput::DoFType & dofType ) const
{
  WALBERLA_ASSERT(    dofType == VTKOutput::DoFType::EDGE_HORIZONTAL
                   || dofType == VTKOutput::DoFType::EDGE_VERTICAL
                   || dofType == VTKOutput::DoFType::EDGE_DIAGONAL );

  output << "<Points>\n";
  output << "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";

  for ( const auto & it : storage->getFaces() )
  {
    Face &face = *it.second;

    const Point3D faceBottomLeftCoords  = face.coords[0];
    const Point3D faceBottomRightCoords = face.coords[1];
    const Point3D faceTopLeftCoords     = face.coords[2];

    const Point3D horizontalMicroEdgeOffset = ( ( faceBottomRightCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( level ) ) ) * 0.5;
    const Point3D verticalMicroEdgeOffset   = ( ( faceTopLeftCoords     - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( level ) ) ) * 0.5;

    switch ( dofType )
    {
    case DoFType::EDGE_HORIZONTAL:
    {
      for ( const auto & itIdx : indexing::edgedof::macroface::Iterator( level, 0 ) )
      {
        const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( real_c( itIdx.col() * 2 + 1 ) * horizontalMicroEdgeOffset + real_c( itIdx.row() * 2     ) * verticalMicroEdgeOffset );
        output << horizontalMicroEdgePosition[0] << " " << horizontalMicroEdgePosition[1] << " " << horizontalMicroEdgePosition[2] << "\n";
      }
      break;
    }
    case DoFType::EDGE_VERTICAL:
    {
      for ( const auto & itIdx : indexing::edgedof::macroface::Iterator( level, 0 ) )
      {
        const Point3D verticalMicroEdgePosition   = faceBottomLeftCoords + ( real_c( itIdx.col() * 2     ) * horizontalMicroEdgeOffset + real_c( itIdx.row() * 2 + 1 ) * verticalMicroEdgeOffset );
        output << verticalMicroEdgePosition[0]   << " " << verticalMicroEdgePosition[1]   << " " << verticalMicroEdgePosition[2]   << "\n";
      }
      break;
    }
    case DoFType::EDGE_DIAGONAL:
    {
      for ( const auto & itIdx : indexing::edgedof::macroface::Iterator( level, 0 ) )
      {
        const Point3D horizontalMicroEdgePosition = faceBottomLeftCoords + ( real_c( itIdx.col() * 2 + 1 ) * horizontalMicroEdgeOffset + real_c( itIdx.row() * 2     ) * verticalMicroEdgeOffset );
        const Point3D diagonalMicroEdgePosition   = horizontalMicroEdgePosition + verticalMicroEdgeOffset;
        output << diagonalMicroEdgePosition[0]   << " " << diagonalMicroEdgePosition[1]   << " " << diagonalMicroEdgePosition[2]   << "\n";
      }
      break;
    }
    default:
      WALBERLA_ABORT( "Bad DoF type in VTK output for edge DoFs" );
      break;
    }
  }

  output << "\n</DataArray>\n";
  output << "</Points>\n";
}

void VTKOutput::writeCells( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & faceWidth ) const
{
  const uint_t numberOfCells = (((faceWidth - 1) * faceWidth) / 2) + (((faceWidth - 2) * (faceWidth - 1)) / 2);

  output << "<Cells>\n";
  output << "<DataArray type=\"Int32\" Name=\"connectivity\">\n";

  // connectivity
  size_t offset = 0;

  for (auto & it : storage->getFaces()) {
    //TODO is it really unused?
    WALBERLA_UNUSED(it);
    size_t rowsize = faceWidth - 1;
    size_t inner_rowsize = rowsize;

    for (size_t i = 0; i < rowsize; ++i)
    {
      for (size_t j = 0; j < inner_rowsize-1; ++j)
      {
        output << offset << " " << offset + 1 << " " << offset + inner_rowsize + 1 << " ";
        output << offset + 1 << " " << offset + inner_rowsize + 2 << " " << offset + inner_rowsize + 1 << " ";
        ++offset;
      }

      output << offset << " " << offset + 1 << " " << offset + inner_rowsize + 1 << " ";

      offset += 2;
      --inner_rowsize;
    }

    ++offset;
  }

  output << "\n</DataArray>\n";
  output << "<DataArray type=\"Int32\" Name=\"offsets\">\n";

  // offsets
  offset = 3;
  for (auto& it : storage->getFaces()) {
    WALBERLA_UNUSED(it);

    for (size_t i = 0; i < numberOfCells; ++i)
    {
      output << offset << " ";
      offset += 3;
    }
  }

  output << "\n</DataArray>\n";
  output << "<DataArray type=\"UInt8\" Name=\"types\">\n";

  // cell types
  for (auto& it : storage->getFaces()) {
    WALBERLA_UNUSED(it);
    for (size_t i = 0; i < numberOfCells; ++i)
    {
      output << "5 ";
    }
  }

  output << "\n</DataArray>\n";
  output << "</Cells>\n";
}

void VTKOutput::writeP1( std::ostream & output, const uint_t & level ) const
{
  if ( p1Functions_.size() == 0 )
  {
    return;
  }

  auto & storage = p1Functions_[0]->getStorage();

  const uint_t numberOfPoints = storage->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level );
  const uint_t numberOfCells  = storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level );

  writePieceHeader( output, numberOfPoints, numberOfCells );

  writePointsForMicroVertices( output, storage, level );

  writeCells( output, storage, levelinfo::num_microvertices_per_edge( level ) );

  output << "<PointData>\n";

  for ( const auto & function : p1Functions_ )
  {
    output << "<DataArray type=\"Float64\" Name=\"" << function->getFunctionName() <<  "\" NumberOfComponents=\"1\">\n";

    for ( const auto & it : storage->getFaces() )
    {
      const Face &face = *it.second;

      size_t len = levelinfo::num_microvertices_per_face( level );
      output << std::scientific;

      for ( size_t i = 0; i < len; ++i )
      {
        output << face.getData( function->getFaceDataID() )->getPointer( level )[i] << " ";
      }
    }
    output << "\n</DataArray>\n";
  }

  output << "</PointData>\n";

  writePieceFooter( output );
}


void VTKOutput::writeEdgeDoFs( std::ostream & output, const uint_t & level, const VTKOutput::DoFType & dofType ) const
{
  WALBERLA_ASSERT(    dofType == VTKOutput::DoFType::EDGE_HORIZONTAL
                   || dofType == VTKOutput::DoFType::EDGE_VERTICAL
                   || dofType == VTKOutput::DoFType::EDGE_DIAGONAL );

  if ( edgeDoFFunctions_.size() == 0 )
  {
    return;
  }

  auto & storage = edgeDoFFunctions_[0]->getStorage();

  const uint_t numberOfPoints = storage->getNumberOfLocalFaces() * levelinfo::num_microedges_per_face( level ) / 3;
  const uint_t faceWidth = levelinfo::num_microedges_per_edge( level );
  const uint_t numberOfCells = storage->getNumberOfLocalFaces() * ((((faceWidth - 1) * faceWidth) / 2) + (((faceWidth - 2) * (faceWidth - 1)) / 2));

  writePieceHeader( output, numberOfPoints, numberOfCells );

  writePointsForMicroEdges( output, storage, level, dofType );

  output << "<PointData>\n";

  for ( const auto & function : edgeDoFFunctions_ )
  {
    output << "<DataArray type=\"Float64\" Name=\"" << function->getFunctionName() <<  "\" NumberOfComponents=\"1\">\n";

    for ( const auto & it : storage->getFaces() )
    {
      const Face & face = *it.second;

      output << std::scientific;

      switch ( dofType )
      {
      case VTKOutput::DoFType::EDGE_HORIZONTAL:
      {
        for ( const auto & itIdx : indexing::edgedof::macroface::Iterator( level ) )
        {
          output << face.getData( function->getFaceDataID() )->getPointer( level )[ vtkDetail::horizontalEdgeOnMacroFaceIndex( level, itIdx.col(), itIdx.row() ) ] << "\n";
        }
        break;
      }
      case VTKOutput::DoFType::EDGE_VERTICAL:
      {
        for ( const auto & itIdx : indexing::edgedof::macroface::Iterator( level ) )
        {
          output << face.getData( function->getFaceDataID() )->getPointer( level )[ vtkDetail::verticalEdgeOnMacroFaceIndex( level, itIdx.col(), itIdx.row() ) ] << "\n";
        }
        break;
      }
      case VTKOutput::DoFType::EDGE_DIAGONAL:
      {
        for ( const auto & itIdx : indexing::edgedof::macroface::Iterator( level ) )
        {
          output << face.getData( function->getFaceDataID() )->getPointer( level )[ vtkDetail::diagonalEdgeOnMacroFaceIndex( level, itIdx.col(), itIdx.row() ) ] << "\n";
        }
        break;
      }
      default:
        WALBERLA_ABORT( "Bad DoF type in VTK output for edge DoFs" );
        break;
      }

    }
    output << "\n</DataArray>\n";
  }

  output << "</PointData>\n";

  writeCells( output, storage, faceWidth );

  writePieceFooter( output );

}

void VTKOutput::writeDGDoFs( std::ostream & output, const uint_t & level ) const
{
  if ( dgFunctions_.size() == 0 )
  {
    return;
  }

  auto & storage = dgFunctions_[0]->getStorage();

  const uint_t numberOfPoints = storage->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level );
  const uint_t numberOfCells  = storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level );

  writePieceHeader( output, numberOfPoints, numberOfCells );

  writePointsForMicroVertices( output, storage, level );

  writeCells( output, storage, levelinfo::num_microvertices_per_edge( level ) );

  output << "<CellData>";

  for ( const auto & function : dgFunctions_ )
  {
    output << "<DataArray type=\"Float64\" Name=\"" << function->getFunctionName() << "\" NumberOfComponents=\"1\">\n";
    for ( const auto & it : storage->getFaces() )
    {
      const Face & face = *it.second;

      uint_t rowsize = levelinfo::num_microvertices_per_edge( level );
      uint_t inner_rowsize = rowsize;
      output << std::scientific;

      uint_t idx;

      for ( size_t j = 0; j < rowsize - 1; ++j )
      {
        for ( size_t i = 0; i < inner_rowsize - 2; ++i )
        {
          idx = vtkDetail::bubbleGrayFaceIndex( level, i, j, stencilDirection::CELL_GRAY_C );
          output << face.getData( function->getFaceDataID() )->getPointer( level )[idx] << " ";
          idx = vtkDetail::bubbleBlueFaceIndex( level, i, j, stencilDirection::CELL_BLUE_C );
          output << face.getData( function->getFaceDataID() )->getPointer( level )[idx] << " ";
        }
        idx = vtkDetail::bubbleGrayFaceIndex( level, inner_rowsize - 2, j, stencilDirection::CELL_GRAY_C );
        output << face.getData( function->getFaceDataID() )->getPointer( level )[idx] << " ";
        --inner_rowsize;
      }
    }
    output << "\n</DataArray>\n";
  }

  output << "\n</CellData>\n";

  writePieceFooter( output );
}


void VTKOutput::writeDoFByType( std::ostream & output, const uint_t & level, const VTKOutput::DoFType & dofType ) const
{
  switch ( dofType )
  {
  case DoFType::VERTEX:
    writeP1( output, level );
    break;
  case DoFType::EDGE_HORIZONTAL:
  case DoFType::EDGE_VERTICAL:
  case DoFType::EDGE_DIAGONAL:
    writeEdgeDoFs( output, level, dofType );
    break;
  case DoFType::DG:
    writeDGDoFs( output, level );
    break;
  default:
    WALBERLA_ABORT( "[VTK] DoFType not supported!" );
    break;
  }
}

uint_t VTKOutput::getNumRegisteredFunctions( const VTKOutput::DoFType & dofType ) const
{
  switch ( dofType )
  {
  case DoFType::VERTEX:
    return p1Functions_.size();
  case DoFType::EDGE_HORIZONTAL:
  case DoFType::EDGE_VERTICAL:
  case DoFType::EDGE_DIAGONAL:
    return edgeDoFFunctions_.size();
  case DoFType::DG:
    return dgFunctions_.size();
    break;
  default:
    WALBERLA_ABORT( "[VTK] DoFType not supported!" );
    return 0;
  }
}


void VTKOutput::write( const uint_t & level, const uint_t & timestep ) const
{
  if ( writeFrequency_ > 0 && timestep % writeFrequency_ == 0 )
  {
    const std::vector< VTKOutput::DoFType > dofTypes = { DoFType::VERTEX, DoFType::EDGE_HORIZONTAL, DoFType::EDGE_VERTICAL, DoFType::EDGE_DIAGONAL, DoFType::DG };

    for ( const auto & dofType : dofTypes )
    {
      if ( getNumRegisteredFunctions( dofType ) > 0 )
      {
        const std::string completeFilePath( fmt::format( "{}/{}{}.vtu", dir_, filename_, fileNameExtension( dofType, level, timestep ) ) );

        WALBERLA_LOG_PROGRESS_ON_ROOT( "[VTK] Writing output to " << completeFilePath );

        std::ostringstream output;

        writeXMLHeader( output );

        writeDoFByType( output, level, dofType );

        walberla::mpi::writeMPITextFile( completeFilePath, output.str() );

        WALBERLA_ROOT_SECTION()
        {
          std::ofstream pvtu_file;
          pvtu_file.open( completeFilePath.c_str(), std::ofstream::out | std::ofstream::app );
          WALBERLA_CHECK( !!pvtu_file, "[VTKWriter] Error opening file: " << completeFilePath );
          writeXMLFooter( pvtu_file );
          pvtu_file.close();
        }
      }
    }
  }

}

}
