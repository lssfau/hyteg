
#pragma once

#include <tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp>
#include "core/mpi/MPITextFile.h"

#include "levelinfo.hpp"
#include "Function.hpp"
#include "p1functionspace/P1Function.hpp"

#include "dgfunctionspace/DGFunction.hpp"

#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFIndexing.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"

#include <string>

namespace hhg
{

namespace vtkDetail
{
SPECIALIZE( uint_t, vertexdof::macrocell::index, vertexDoFOnMacroCellIndex )
SPECIALIZE(uint_t, vertexdof::macroface::indexFromVertex, vertexDoFOnMacroFaceIndex)

SPECIALIZE(uint_t,BubbleFace::indexFaceFromGrayFace,bubbleGrayFaceIndex)
SPECIALIZE(uint_t,BubbleFace::indexFaceFromBlueFace,bubbleBlueFaceIndex)

SPECIALIZE(uint_t, edgedof::macroface::horizontalIndex, horizontalEdgeOnMacroFaceIndex)
SPECIALIZE(uint_t, edgedof::macroface::verticalIndex,   verticalEdgeOnMacroFaceIndex)
SPECIALIZE(uint_t, edgedof::macroface::diagonalIndex,   diagonalEdgeOnMacroFaceIndex)
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
     dir_( dir ), filename_( filename ), writeFrequency_( writeFrequency ), write2D_( true )
  {}

  /// Writes only macro-faces.
  void set2D() { write2D_ = true; }
  /// Writes only macro-cells.
  void set3D() { write2D_ = false; }

  void add( const P1Function     < real_t > * function ) { p1Functions_.push_back( function ); } ;
  void add( const EdgeDoFFunction< real_t > * function ) { edgeDoFFunctions_.push_back( function ); };
  void add( const BubbleFunction < real_t > * function ) { bubbleFunctions_.push_back( function ); };
  void add( const DGFunction     < real_t > * function ) { dgFunctions_.push_back( function ); };

  void add( const P2Function     < real_t > * function ) { p2Functions_.push_back( function );
                                                           p1Functions_.push_back( function->getVertexDoFFunction().get() );
                                                           edgeDoFFunctions_.push_back( function->getEdgeDoFFunction().get() ); }

  void add( const std::shared_ptr< P1Function     < real_t > > & function ) { p1Functions_.push_back( function.get() ); } ;
  void add( const std::shared_ptr< EdgeDoFFunction< real_t > > & function ) { edgeDoFFunctions_.push_back( function.get() ); };
  void add( const std::shared_ptr< BubbleFunction < real_t > > & function ) { bubbleFunctions_.push_back( function.get() ); };
  void add( const std::shared_ptr< DGFunction     < real_t > > & function ) { dgFunctions_.push_back( function.get() ); };

  void add( const std::shared_ptr< P2Function     < real_t > > & function ) { p2Functions_.push_back( function.get() );
                                                                              p1Functions_.push_back( function->getVertexDoFFunction().get() );
                                                                              edgeDoFFunctions_.push_back( function->getEdgeDoFFunction().get() ); }

  /// Writes the VTK output only if writeFrequency > 0 and timestep % writeFrequency == 0.
  /// Therefore always writes output if timestep is 0.
  /// Appends the time step to the filename.
  /// Note: files will be overwritten if called twice with the same time step!
  void write( const uint_t & level, const uint_t & timestep = 0 ) const;

private:

  enum class DoFType
  {
    VERTEX,
    EDGE_HORIZONTAL,
    EDGE_VERTICAL,
    EDGE_DIAGONAL,
    DG,
    P2
  };

  static const std::map< VTKOutput::DoFType, std::string > DoFTypeToString_;

  void writeDoFByType( std::ostream & output, const uint_t & level, const VTKOutput::DoFType & dofType ) const;
  uint_t getNumRegisteredFunctions( const VTKOutput::DoFType & dofType ) const;

  void writeP1      ( std::ostream & output, const uint_t & level                                     ) const;
  void writeEdgeDoFs( std::ostream & output, const uint_t & level, const VTKOutput::DoFType & dofType ) const;
  void writeDGDoFs  ( std::ostream & output, const uint_t & level                                     ) const;
  void writeP2      ( std::ostream & output, const uint_t & level                                     ) const;

  std::string fileNameExtension( const VTKOutput::DoFType & dofType, const uint_t & level, const uint_t & timestep ) const;

  void writeHeader       ( std::ostringstream & output, const uint_t & numberOfPoints, const uint_t & numberOfCells ) const;
  void writeFooterAndFile( std::ostringstream & output, const std::string & completeFilePath ) const;

  void writePointsForMicroVertices( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level ) const;
  void writePointsForMicroEdges   ( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level, const VTKOutput::DoFType & dofType ) const;

  void writeVertexDoFData( std::ostream & output, const vertexdof::VertexDoFFunction< real_t > * function,
                           const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level ) const;
  void writeEdgeDoFData  ( std::ostream & output, const EdgeDoFFunction< real_t > * function,
                           const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & level, const DoFType & dofType ) const;

  void writeCells( std::ostream & output, const std::shared_ptr< PrimitiveStorage > & storage, const uint_t & faceWidth ) const;

  void syncAllFunctions( const uint_t & level ) const;

  std::string dir_;
  std::string filename_;

  uint_t writeFrequency_;

  bool write2D_;

  std::vector< const P1Function     < real_t > * > p1Functions_;
  std::vector< const EdgeDoFFunction< real_t > * > edgeDoFFunctions_;
  std::vector< const BubbleFunction < real_t > * > bubbleFunctions_;
  std::vector< const DGFunction     < real_t > * > dgFunctions_;

  std::vector< const P2Function     < real_t > * > p2Functions_;

};


}

