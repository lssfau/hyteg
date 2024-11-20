/*
 * Copyright (c) 2017-2024 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"

#include "core/Format.hpp"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKDGWriter.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::int32_c;
using walberla::real_c;
using walberla::uint32_c;
using walberla::vtk::typeToString;

VTKOutput::VTKOutput( std::string                                dir,
                      std::string                                filename,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      const uint_t&                              writeFrequency )
: dir_( std::move( dir ) )
, filename_( std::move( filename ) )
, writeFrequency_( writeFrequency )
, write2D_( true )
, storage_( storage )
, vtkDataFormat_( vtk::DataFormat::ASCII )
{
   // set output to 3D if storage contains cells
   if ( storage->hasGlobalCells() )
   {
      set3D();
   }
}

const std::map< vtk::DoFType, std::string > VTKOutput::DoFTypeToString_ = {
    { vtk::DoFType::VERTEX, "VertexDoF" },
    { vtk::DoFType::EDGE_X, "XEdgeDoF" },
    { vtk::DoFType::EDGE_Y, "YEdgeDoF" },
    { vtk::DoFType::EDGE_Z, "ZEdgeDoF" },
    { vtk::DoFType::EDGE_XY, "XYEdgeDoF" },
    { vtk::DoFType::EDGE_XZ, "XZEdgeDoF" },
    { vtk::DoFType::EDGE_YZ, "YZEdgeDoF" },
    { vtk::DoFType::EDGE_XYZ, "XYZEdgeDoF" },
    { vtk::DoFType::DG, "DGDoF" },
    { vtk::DoFType::P2, "P2" },
    { vtk::DoFType::P2_PLUS_BUBBLE, "P2+Bubble" },
    { vtk::DoFType::N1E1, "N1E1" },
    { vtk::DoFType::P1DGE, "P1DGE" },
};

std::string VTKOutput::fileNameExtension( const vtk::DoFType& dofType, const uint_t& level, const uint_t& timestep ) const
{
   return walberla::format( "_%s_level%u_ts%u", VTKOutput::DoFTypeToString_.at( dofType ).c_str(), level, timestep );
}

void VTKOutput::writeDoFByType( std::ostream& output, const uint_t& level, const vtk::DoFType& dofType ) const
{
   switch ( dofType )
   {
   case vtk::DoFType::VERTEX:
      VTKP1Writer::write( *this, output, level );
      break;
   case vtk::DoFType::EDGE_X:
   case vtk::DoFType::EDGE_Y:
   case vtk::DoFType::EDGE_Z:
   case vtk::DoFType::EDGE_XY:
   case vtk::DoFType::EDGE_XZ:
   case vtk::DoFType::EDGE_YZ:
   case vtk::DoFType::EDGE_XYZ:
      VTKEdgeDoFWriter::write( *this, output, level, dofType );
      break;
   case vtk::DoFType::DG:
      VTKDGWriter::write( *this, output, level );
      break;
   case vtk::DoFType::P2:
      VTKP2Writer::write( *this, output, level );
      break;
   case vtk::DoFType::P2_PLUS_BUBBLE:
      VTKP2PlusBubbleWriter::write( *this, output, level );
      break;
   case vtk::DoFType::N1E1:
      VTKN1E1Writer::write( *this, output, level );
      break;
   case vtk::DoFType::P1DGE:
      VTKP1DGEWriter::write( *this, output, level );
      break;
   default:
      WALBERLA_ABORT( "[VTK] DoFType not supported!" );
      break;
   }
}

uint_t VTKOutput::getNumRegisteredFunctions( const vtk::DoFType& dofType ) const
{
   switch ( dofType )
   {
   case vtk::DoFType::VERTEX:
      return feFunctionRegistry_.getP1Functions().size() + feFunctionRegistry_.getP1VectorFunctions().size();
   case vtk::DoFType::EDGE_X:
   case vtk::DoFType::EDGE_Y:
   case vtk::DoFType::EDGE_Z:
   case vtk::DoFType::EDGE_XY:
   case vtk::DoFType::EDGE_XZ:
   case vtk::DoFType::EDGE_YZ:
   case vtk::DoFType::EDGE_XYZ:
      return feFunctionRegistry_.getEdgeDoFFunctions().size();
   case vtk::DoFType::DG:
      return feFunctionRegistry_.getDGFunctions().size() + feFunctionRegistry_.getDGVectorFunctions().size();
   case vtk::DoFType::P2:
      return feFunctionRegistry_.getP2Functions().size() + feFunctionRegistry_.getP2VectorFunctions().size();
   case vtk::DoFType::P2_PLUS_BUBBLE:
      return feFunctionRegistry_.getP2PlusBubbleFunctions().size(); // + feFunctionRegistry_.getP2PlusBubbleVectorFunctions().size();
   case vtk::DoFType::P1DGE:
      return feFunctionRegistry_.getEGFunctions().size();
      break;
   case vtk::DoFType::N1E1:
      return feFunctionRegistry_.getN1E1VectorFunctions().size();
      break;
   default:
      WALBERLA_ABORT( "[VTK] DoFType not supported!" );
      return 0;
   }
}

void VTKOutput::write( const uint_t level, const uint_t timestep )
{
   storage_->getTimingTree()->start( "VTK write" );

   if ( writeFrequency_ > 0 && timestep % writeFrequency_ == 0 )
   {
      micromesh::communicate( storage_, level );
      bool excludeDG = true;
      communication::syncRegisteredFunctions( feFunctionRegistry_, level, excludeDG, communication::syncDirection_t::LOW2HIGH );

      const std::vector< vtk::DoFType > dofTypes2D = { vtk::DoFType::VERTEX,
                                                       vtk::DoFType::EDGE_X,
                                                       vtk::DoFType::EDGE_Y,
                                                       vtk::DoFType::EDGE_XY,
                                                       vtk::DoFType::DG,
                                                       vtk::DoFType::P2,
                                                       vtk::DoFType::P2_PLUS_BUBBLE,
                                                       vtk::DoFType::P1DGE };

      const std::vector< vtk::DoFType > dofTypes3D = { vtk::DoFType::VERTEX,
                                                       vtk::DoFType::EDGE_X,
                                                       vtk::DoFType::EDGE_Y,
                                                       vtk::DoFType::EDGE_Z,
                                                       vtk::DoFType::EDGE_XY,
                                                       vtk::DoFType::EDGE_XZ,
                                                       vtk::DoFType::EDGE_YZ,
                                                       vtk::DoFType::EDGE_XYZ,
                                                       vtk::DoFType::DG,
                                                       vtk::DoFType::P2,
                                                       vtk::DoFType::P1DGE,
                                                       vtk::DoFType::N1E1 };

      auto dofTypes = write2D_ ? dofTypes2D : dofTypes3D;

      for ( const auto& dofType : dofTypes )
      {
         if ( getNumRegisteredFunctions( dofType ) > 0 )
         {
            const std::string completeFilePath = walberla::format(
                "%s/%s%s.vtu", dir_.c_str(), filename_.c_str(), fileNameExtension( dofType, level, timestep ).c_str() );
            //( fmt::format( "{}/{}{}.vtu", dir_, filename_, fileNameExtension( dofType, level, timestep ) ) );

            std::ostringstream output;

            vtk::writeXMLHeader( output );

            writeDoFByType( output, level, dofType );

            walberla::mpi::writeMPITextFile( completeFilePath, output.str() );

            WALBERLA_ROOT_SECTION()
            {
               std::ofstream pvtu_file;
               pvtu_file.open( completeFilePath.c_str(), std::ofstream::out | std::ofstream::app );
               WALBERLA_CHECK( !!pvtu_file, "[VTKWriter] Error opening file: " << completeFilePath );
               vtk::writeXMLFooter( pvtu_file );
               pvtu_file.close();
            }
         }
      }
   }

   storage_->getTimingTree()->stop( "VTK write" );
}

// -------------------------
//  Explicit Instantiations
// -------------------------
template void VTKOutput::add( const GenericFunction< double >& function );
template void VTKOutput::add( const GenericFunction< float >& function );
template void VTKOutput::add( const GenericFunction< int32_t >& function );
template void VTKOutput::add( const GenericFunction< int64_t >& function );

} // namespace hyteg
