/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "hyteg/dataexport/VTKOutput.hpp"

#include "core/Format.hpp"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/facedofspace/FaceDoFIndexing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

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
   /// set output to 3D is storage contains cells
   if ( storage->hasGlobalCells() )
   {
      set3D();
   }
}

void VTKOutput::add( const P2Function< real_t >& function )
{
   p2Functions_.push_back( function );
}

void VTKOutput::add( const P1VectorFunction< real_t >& function )
{
   p1VecFunctions_.push_back( function );
}

void VTKOutput::add( const P2VectorFunction< real_t >& function )
{
   p2VecFunctions_.push_back( function );
}

void VTKOutput::add( const P1Function< real_t >& function )
{
   p1Functions_.push_back( function );
}

void VTKOutput::add( const EdgeDoFFunction< real_t >& function )
{
   edgeDoFFunctions_.push_back( function );
}

void VTKOutput::add( const DGFunction< real_t >& function )
{
   dgFunctions_.push_back( function );
}

void VTKOutput::add( const P1StokesFunction< real_t >& function )
{
   add( function.uvw );
   add( function.p );
}

void VTKOutput::add( const P2P1TaylorHoodFunction< real_t >& function )
{
   add( function.uvw );
   add( function.p );
}

void VTKOutput::add( const BlockFunction< real_t >& function )
{
   for ( uint_t k = 0; k < function.getNumberOfBlocks(); k++ )
   {
      add( function[k] );
   }
}

void VTKOutput::add( const GenericFunction< real_t >& function )
{
   bool matchFound = false;
   switch ( function.getFunctionKind() )
   {
   case functionTraits::P1_FUNCTION:
      matchFound = tryUnwrapAndAdd< FunctionWrapper< P1Function< real_t > > >( function );
      break;

   case functionTraits::P2_FUNCTION:
      matchFound = tryUnwrapAndAdd< FunctionWrapper< P2Function< real_t > > >( function );
      break;

   case functionTraits::P1_VECTOR_FUNCTION:
      matchFound = tryUnwrapAndAdd< FunctionWrapper< P1VectorFunction< real_t > > >( function );
      break;

   case functionTraits::P2_VECTOR_FUNCTION:
      matchFound = tryUnwrapAndAdd< FunctionWrapper< P2VectorFunction< real_t > > >( function );
      break;

   case functionTraits::EDGE_DOF_FUNCTION:
      matchFound = tryUnwrapAndAdd< FunctionWrapper< EdgeDoFFunction< real_t > > >( function );
      break;

   case functionTraits::DG_FUNCTION:
      matchFound = tryUnwrapAndAdd< FunctionWrapper< DGFunction< real_t > > >( function );
      break;

   default:
      matchFound = false;
   }

   if ( !matchFound )
   {
      WALBERLA_ABORT( "VTKOutput: Failed to add GenericFunction object!" );
   }
}

const std::map< vtk::DoFType, std::string > VTKOutput::DoFTypeToString_ = {
    {vtk::DoFType::VERTEX, "VertexDoF"},
    {vtk::DoFType::EDGE_X, "XEdgeDoF"},
    {vtk::DoFType::EDGE_Y, "YEdgeDoF"},
    {vtk::DoFType::EDGE_Z, "ZEdgeDoF"},
    {vtk::DoFType::EDGE_XY, "XYEdgeDoF"},
    {vtk::DoFType::EDGE_XZ, "XZEdgeDoF"},
    {vtk::DoFType::EDGE_YZ, "YZEdgeDoF"},
    {vtk::DoFType::EDGE_XYZ, "XYZEdgeDoF"},
    {vtk::DoFType::DG, "DGDoF"},
    {vtk::DoFType::P2, "P2"},
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
      VTKDGDoFWriter::write( *this, output, level );
      break;
   case vtk::DoFType::P2:
      VTKP2Writer::write( *this, output, level );
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
      return p1Functions_.size() + p1VecFunctions_.size();
   case vtk::DoFType::EDGE_X:
   case vtk::DoFType::EDGE_Y:
   case vtk::DoFType::EDGE_Z:
   case vtk::DoFType::EDGE_XY:
   case vtk::DoFType::EDGE_XZ:
   case vtk::DoFType::EDGE_YZ:
   case vtk::DoFType::EDGE_XYZ:
      return edgeDoFFunctions_.size();
   case vtk::DoFType::DG:
      return dgFunctions_.size();
      break;
   case vtk::DoFType::P2:
      return p2Functions_.size() + p2VecFunctions_.size();
      break;
   default:
      WALBERLA_ABORT( "[VTK] DoFType not supported!" );
      return 0;
   }
}

void VTKOutput::write( const uint_t& level, const uint_t& timestep ) const
{
   if ( level <= 1 )
   {
      return;
   }

   storage_->getTimingTree()->start( "VTK write" );

   if ( writeFrequency_ > 0 && timestep % writeFrequency_ == 0 )
   {
      syncAllFunctions( level );

      const std::vector< vtk::DoFType > dofTypes2D = {vtk::DoFType::VERTEX,
                                                      vtk::DoFType::EDGE_X,
                                                      vtk::DoFType::EDGE_Y,
                                                      vtk::DoFType::EDGE_XY,
                                                      vtk::DoFType::DG,
                                                      vtk::DoFType::P2};

      const std::vector< vtk::DoFType > dofTypes3D = {vtk::DoFType::VERTEX,
                                                      vtk::DoFType::EDGE_X,
                                                      vtk::DoFType::EDGE_Y,
                                                      vtk::DoFType::EDGE_Z,
                                                      vtk::DoFType::EDGE_XY,
                                                      vtk::DoFType::EDGE_XZ,
                                                      vtk::DoFType::EDGE_YZ,
                                                      vtk::DoFType::EDGE_XYZ,
                                                      vtk::DoFType::DG,
                                                      vtk::DoFType::P2};

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

void VTKOutput::syncAllFunctions( const uint_t& level ) const
{
   for ( const auto& function : p1Functions_ )
   {
      hyteg::communication::syncFunctionBetweenPrimitives< hyteg::P1Function< real_t > >( function, level );
   }

   for ( const auto& function : p1VecFunctions_ )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }

   for ( const auto& function : p2Functions_ )
   {
      hyteg::communication::syncP2FunctionBetweenPrimitives( function, level );
   }

   for ( const auto& function : p2VecFunctions_ )
   {
      hyteg::communication::syncVectorFunctionBetweenPrimitives( function, level );
   }

   for ( const auto& function : edgeDoFFunctions_ )
   {
      hyteg::communication::syncFunctionBetweenPrimitives< hyteg::EdgeDoFFunction< real_t > >( function, level );
   }

   for ( const auto& function : dgFunctions_ )
   {
      function.communicate< Vertex, Edge >( level );
      function.communicate< Edge, Face >( level );
   }
}

} // namespace hyteg
