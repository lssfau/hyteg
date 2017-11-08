
#pragma once

#include "core/mpi/MPITextFile.h"

#include "levelinfo.hpp"
#include "Function.hpp"
#include "p1functionspace/P1Function.hpp"
#include "p1functionspace/P1Memory.hpp"

#include "dgfunctionspace/DGFunction.hpp"

#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/indexing/EdgeDoFIndexing.hpp"

#include <string>

namespace hhg
{

namespace vtkDetail{
SPECIALIZE(uint_t,BubbleFace::indexFaceFromGrayFace,bubbleGrayFaceIndex)
SPECIALIZE(uint_t,BubbleFace::indexFaceFromBlueFace,bubbleBlueFaceIndex)

SPECIALIZE(uint_t, indexing::edgedof::macroface::horizontalIndex, horizontalEdgeOnMacroFaceIndex)
SPECIALIZE(uint_t, indexing::edgedof::macroface::verticalIndex,   verticalEdgeOnMacroFaceIndex)
SPECIALIZE(uint_t, indexing::edgedof::macroface::diagonalIndex,   diagonalEdgeOnMacroFaceIndex)
}

using walberla::uint_t;
using walberla::uint_c;
using walberla::real_t;
using walberla::real_c;
////FIXME this typedef can be remove when we move into walberla namespace
typedef walberla::uint64_t uint64_t;

class VTKOutput
{
public:

  /// \param writeFrequency outputs VTK in the specified frequency
  VTKOutput( const std::string & dir, const std::string & filename, const uint_t & writeFrequency = 1 ) :
     dir_( dir ), filename_( filename ), writeFrequency_( writeFrequency )
  {}

  void add( const P1Function     < real_t > * function ) { p1Functions_.push_back( function ); } ;
  void add( const EdgeDoFFunction< real_t > * function ) { edgeDoFFunctions_.push_back( function ); };
  void add( const BubbleFunction < real_t > * function ) { bubbleFunctions_.push_back( function ); };
  void add( const DGFunction     < real_t > * function ) { dgFunctions_.push_back( function ); };

  void add( const std::shared_ptr< P1Function     < real_t > > & function ) { p1Functions_.push_back( function.get() ); } ;
  void add( const std::shared_ptr< EdgeDoFFunction< real_t > > & function ) { edgeDoFFunctions_.push_back( function.get() ); };
  void add( const std::shared_ptr< BubbleFunction < real_t > > & function ) { bubbleFunctions_.push_back( function.get() ); };
  void add( const std::shared_ptr< DGFunction     < real_t > > & function ) { dgFunctions_.push_back( function.get() ); };

  /// Writes the VTK output only if writeFrequency > 0 and timestep % writeFrequency == 0
  /// Appends the time step to the filename.
  void write( const uint_t & level, const uint_t & timestep = 0 ) const;

private:

  enum class DoFType
  {
    VERTEX,
    EDGE_HORIZONTAL,
    EDGE_VERTICAL,
    EDGE_DIAGONAL,
    DG
  };

  static const std::map< VTKOutput::DoFType, std::string > DoFTypeToString_;

  void writeDoFByType( std::ostream & output, const uint_t & level, const VTKOutput::DoFType & dofType ) const;
  uint_t getNumRegisteredFunctions( const VTKOutput::DoFType & dofType ) const;

  void writeP1      ( std::ostream & output, const uint_t & level                                     ) const;
  void writeEdgeDoFs( std::ostream & output, const uint_t & level, const VTKOutput::DoFType & dofType ) const;
  void writeDGDoFs  ( std::ostream & output, const uint_t & level                                     ) const;

  std::string fileNameExtension( const VTKOutput::DoFType & dofType, const uint_t & level, const uint_t & timestep ) const;

  void writeHeader       ( std::ostringstream & output, const uint_t & numberOfPoints, const uint_t & numberOfCells ) const;
  void writeFooterAndFile( std::ostringstream & output, const std::string & completeFilePath ) const;

  void writePointsForMicroVertices( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level ) const;
  void writePointsForMicroEdges   ( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level, const VTKOutput::DoFType & dofType ) const;

  void writeCells( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & faceWidth ) const;

  std::string dir_;
  std::string filename_;

  uint_t writeFrequency_;

  std::vector< const P1Function     < real_t > * > p1Functions_;
  std::vector< const EdgeDoFFunction< real_t > * > edgeDoFFunctions_;
  std::vector< const BubbleFunction < real_t > * > bubbleFunctions_;
  std::vector< const DGFunction     < real_t > * > dgFunctions_;

};




template< typename ContinuousFunctionType, typename DiscontinuousFunctionType>
void VTKWriter(std::vector<const Function<ContinuousFunctionType> *> functionsC,
               std::vector<const Function<DiscontinuousFunctionType> *> functionsD, const uint_t level, const std::string &dir,
               const std::string &filename)
{
  uint_t rk = uint_c(walberla::mpi::MPIManager::instance()->rank());

  std::string pvtu_filename(fmt::format("{}/{}.vtu", dir, filename));
  WALBERLA_LOG_INFO_ON_ROOT("[VTKWriter] Writing functions on level " << level << " to '" << pvtu_filename << "'");

  auto& storage = functionsC[0]->getStorage();

  // prepare local piece
  std::ostringstream output;

  size_t num_faces = storage->getNumberOfLocalFaces();

  if (rk == 0) {
    output << "<?xml version=\"1.0\"?>\n";
    output << "<VTKFile type=\"UnstructuredGrid\">\n";
    output << "<UnstructuredGrid>\n";
  }
  output << "<Piece NumberOfPoints=\"" << num_faces * levelinfo::num_microvertices_per_face(level) << "\" NumberOfCells=\"" << num_faces * levelinfo::num_microfaces_per_face(level) << "\">\n";
  output << "<Points>\n";
  output << "<DataArray type=\"Float64\" NumberOfComponents=\"3\">\n";

  // write out coords
  for (auto& it : storage->getFaces()) {
    Face &face = *it.second;

    size_t rowsize = levelinfo::num_microvertices_per_edge(level);
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
  output << "<Cells>\n";
  output << "<DataArray type=\"Int32\" Name=\"connectivity\">\n";

  // connectivity
  size_t offset = 0;

  for (auto& it : storage->getFaces()) {
    //TODO is it really unused?
    WALBERLA_UNUSED(it);
    size_t rowsize = levelinfo::num_microvertices_per_edge(level) - 1;
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

    for (size_t i = 0; i < levelinfo::num_microfaces_per_face(level); ++i)
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
    for (size_t i = 0; i < levelinfo::num_microfaces_per_face(level); ++i)
    {
      output << "5 ";
    }
  }

  output << "\n</DataArray>\n";
  output << "</Cells>\n";
  output << "<PointData>\n";

  // point data
  for (const auto function : functionsC)
  {
    output << "<DataArray type=\"Float64\" Name=\"" << function->getFunctionName() <<  "\" NumberOfComponents=\"1\">\n";
    for (auto& it : storage->getFaces()) {
      Face &face = *it.second;

      size_t len = levelinfo::num_microvertices_per_face(level);
      output << std::scientific;

      const P1Function< real_t >* p1Function = dynamic_cast<const P1Function< real_t >*>(function);

      for (size_t i = 0; i < len; ++i)
      {
        output << face.getData(p1Function->getFaceDataID())->getPointer(level)[i] << " ";
      }
    }
    output << "\n</DataArray>\n";
  }

  output << "</PointData>\n";
  output << "<CellData>";

  // cell data
  for (const auto function : functionsD)
  {
    output << "<DataArray type=\"Float64\" Name=\"" << function->getFunctionName() <<  "\" NumberOfComponents=\"1\">\n";
    for (auto& it : storage->getFaces()) {
      Face &face = *it.second;

      uint_t rowsize = levelinfo::num_microvertices_per_edge(level);
      uint_t inner_rowsize = rowsize;
      output << std::scientific;

      auto dgFunction = dynamic_cast<const DGFunction< real_t >*>(function);
      if (dgFunction == nullptr) {
        WALBERLA_ABORT("Function is of wrong type");
      }

      uint_t idx;

      for (size_t j = 0; j < rowsize - 1; ++j) {
        for (size_t i = 0; i < inner_rowsize - 2; ++i) {
          idx = vtkDetail::bubbleGrayFaceIndex(level,i, j, stencilDirection::CELL_GRAY_C);
          output << face.getData(dgFunction->getFaceDataID())->getPointer(level)[idx] << " ";
          idx = vtkDetail::bubbleBlueFaceIndex(level,i, j, stencilDirection::CELL_BLUE_C);
          output << face.getData(dgFunction->getFaceDataID())->getPointer(level)[idx] << " ";
        }
        idx = vtkDetail::bubbleGrayFaceIndex(level,inner_rowsize-2, j, stencilDirection::CELL_GRAY_C);
        output << face.getData(dgFunction->getFaceDataID())->getPointer(level)[idx] << " ";
        --inner_rowsize;
      }
    }
    output << "\n</DataArray>\n";
  }

  output << "\n</CellData>\n";
  output << "</Piece>\n";

  walberla::mpi::writeMPITextFile(pvtu_filename, output.str());

  if (rk == 0)
  {
    std::ofstream pvtu_file;
    pvtu_file.open(pvtu_filename.c_str(), std::ofstream::out | std::ofstream::app);

    if(!pvtu_file)
    {
      WALBERLA_ABORT("[VTKWriter] Error opening file: " << pvtu_filename);
    }

    pvtu_file << " </UnstructuredGrid>\n";
    pvtu_file << "</VTKFile>\n";
    pvtu_file.close();
  }

}

}

