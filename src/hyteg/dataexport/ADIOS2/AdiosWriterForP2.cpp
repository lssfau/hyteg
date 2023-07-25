/*
 * Copyright (c) 2023 Marcus Mohr.
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
#include <adios2.h>

#include "hyteg/dataexport/ADIOS2/AdiosHelperFunctions.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriterForP1.hpp"
#include "hyteg/dataexport/FEFunctionRegistry.hpp"
#include "hyteg/dataexport/FEFunctionWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKMeshWriter.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

AdiosWriterForP2::AdiosWriterForP2( adios2::ADIOS&                             adios,
                                    std::string&                               filePath,
                                    std::string&                               fileBaseName,
                                    const std::string&                         engineType,
                                    uint_t                                     level,
                                    const std::shared_ptr< PrimitiveStorage >& storage )
: storage_( storage )
, level_( level )
{
   // create the name of the output file
   std::stringstream tag;
   tag << "-P2_level" << level_ << ".bp";
   fileName_ = filePath + "/" + fileBaseName + tag.str();

   // create our own writer
   io_ = adios.DeclareIO( fileName_ );

   // set the type of engine
   io_.SetEngine( engineType );
}

void AdiosWriterForP2::write( const FEFunctionRegistry& registry, uint_t timestep )
{
   // working on faces/cells some processes might not have data to export
   if ( !adiosHelpers::mpiProcessHasMacrosOfHighestDimension( storage_ ) )
   {
      // however, open() and close() are collective operations
      if ( !firstWriteCompleted_ )
      {
         // open output file
         engine_              = io_.Open( fileName_, adios2::Mode::Write );
         firstWriteCompleted_ = true;
      }

      // also appears to be collective for BP
      engine_.BeginStep();
      engine_.EndStep();

      return;
   }

   // Assemble list of names of node-centred and cell-centred funcions (latter is empty)
   std::vector< std::string > emptyList;
   std::vector< std::string > p2FunctionList;
   registry.extractFunctionNames( p2FunctionList, functionTraits::P2_FUNCTION );
   registry.extractFunctionNames( p2FunctionList, functionTraits::P2_VECTOR_FUNCTION );

   // ------------------
   //  First write only
   // ------------------
   if ( !firstWriteCompleted_ )
   {
      // open output file
      engine_ = io_.Open( fileName_, adios2::Mode::Write );

      // export the mesh information
      engine_.BeginStep();
      writeMesh( p2FunctionList );

      // add meta data on simulation software
      adiosHelpers::generateSoftwareMetaData( io_ );

      // define ADIOS variables to be associated with our functions
      defineVariables< real_t >( registry );

      // done
      firstWriteCompleted_ = true;
   }

   // ----------------------------
   //  First and Follow-up Writes
   // ----------------------------
   else
   {
      engine_.BeginStep();
   }

   // set time-step info
   adiosHelpers::putTimeStepInfo( io_, engine_, timestep );
   for ( const auto& func : registry.getP2Functions().getFunctions< real_t >() )
   {
      scheduleScalarFunctionForExport( func );
   }

   for ( const auto& func : registry.getP2VectorFunctions().getFunctions< real_t >() )
   {
      scheduleVectorFunctionForExport( func );
   }

   engine_.EndStep();
}

void AdiosWriterForP2::writeMesh( const std::vector< std::string >& p2FunctionList )
{
   // store the generic XML file as attribute as part of the binary data
   std::vector< std::string > emptyList;
   std::string                vtkMetaInfo = adiosHelpers::generateVTKMetaInfo( p2FunctionList, emptyList );
   io_.DefineAttribute( "vtk.xml", vtkMetaInfo, "", "" );

   // determine basic entity counts
   const uint_t dim         = storage_->hasGlobalCells() ? 3u : 2u;
   uint_t       numElements = 0u;
   uint_t       numVertices = 0u;
   if ( dim == 3 )
   {
      numVertices = storage_->getNumberOfLocalCells() * levelinfo::num_microvertices_per_cell( level_ + 1 );
      numElements = storage_->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( level_ );
   }
   else
   {
      numVertices = storage_->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level_ + 1 );
      numElements = storage_->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level_ );
   }

   // store entity counts as scalar variables
   adios2::Variable< uint32_t > varNumberOfNodes =
       io_.DefineVariable< uint32_t >( "NumberOfVertices", { adios2::LocalValueDim } );
   adios2::Variable< uint32_t > varNumberOfFaces =
       io_.DefineVariable< uint32_t >( "NumberOfElements", { adios2::LocalValueDim } );

   engine_.Put( varNumberOfNodes, static_cast< uint32_t >( numVertices ) );
   engine_.Put( varNumberOfFaces, static_cast< uint32_t >( numElements ) );

   // store the type of elements in the mesh as scalar variable
   // 22u = "QuadraticTriangle", 24u = "QuadraticTetrahedron"
   using elementMarker_t                               = uint32_t;
   elementMarker_t                     elementMarker   = dim == 2 ? 22u : 24u;
   adios2::Variable< elementMarker_t > varElementTypes = io_.DefineVariable< uint32_t >( "types" );
   engine_.Put( varElementTypes, elementMarker );

   // store node coordinates (always x,y,z even for 2D!)
   adios2::Variable< real_t >                varVertices = io_.DefineVariable< real_t >( "vertices", {}, {}, { numVertices, 3 } );
   adios2::Variable< real_t >::Span          vertices    = engine_.Put< real_t >( varVertices );
   AdiosWriter::StreamAccessBuffer< real_t > vertexStream( vertices, varVertices.Count() );
   VTKMeshWriter::writePointsForMicroVertices( dim == 2, vertexStream, storage_, level_ + 1, false );

   // store element connectivity
   adios2::Variable< uint64_t > varConnectivity =
       io_.DefineVariable< uint64_t >( "connectivity", {}, {}, { numElements, dim == 2 ? 7u : 11u } );
   adios2::Variable< uint64_t >::Span connectivity = engine_.Put< uint64_t >( varConnectivity );

   if ( dim == 3 )
   {
      AdiosWriter::StreamAccessBuffer< uint64_t, 11 > connectivityStream( connectivity, varConnectivity.Count() );
      VTKMeshWriter::writeElementNodeAssociationP2Tetrahedrons( connectivityStream, storage_, level_ );
   }
   else
   {
      AdiosWriter::StreamAccessBuffer< uint64_t, 7 > connectivityStream( connectivity, varConnectivity.Count() );
      VTKMeshWriter::writeElementNodeAssociationP2Triangles( connectivityStream, storage_, level_ );
   }
}

template < typename value_t >
void AdiosWriterForP2::defineVariables( const FEFunctionRegistry& registry )
{
   auto checkAndDefine = [this]( std::string funcName, uint_t funcDim ) {
      adios2::Variable< value_t > varDoFData = io_.InquireVariable< value_t >( funcName );
      if ( !varDoFData )
      {
         uint_t extent{ 0 };
         if ( !storage_->hasGlobalCells() )
         {
            extent = levelinfo::num_microvertices_per_face( level_ + 1 ) * storage_->getNumberOfLocalFaces();
         }
         else
         {
            extent = levelinfo::num_microvertices_per_cell( level_ + 1 ) * storage_->getNumberOfLocalCells();
         }
         if ( funcDim == 1 )
         {
            io_.DefineVariable< value_t >( funcName, {}, {}, { extent } );
         }
         else
         {
            // dofs for vector-valued functions need to be collected together in node-based fashion
            io_.DefineVariable< value_t >( funcName, {}, {}, { extent, funcDim } );
         }
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( "ADIOS2 variable '" << funcName << "' already exists!" );
         WALBERLA_ABORT( "So geht das nicht!" );
      }
   };

   std::vector< P2Function< value_t > > p2Funcs = registry.getP2Functions().getFunctions< value_t >();
   for ( const auto& func : p2Funcs )
   {
      checkAndDefine( func.getFunctionName(), func.getDimension() );
   }

   std::vector< P2VectorFunction< value_t > > p2VecFuncs = registry.getP2VectorFunctions().getFunctions< value_t >();
   for ( const auto& func : p2VecFuncs )
   {
      checkAndDefine( func.getFunctionName(), func.getDimension() );
   }
}

template < typename value_t >
void AdiosWriterForP2::scheduleScalarFunctionForExport( const P2Function< value_t >& func )
{
   adios2::Variable< value_t > varDoFData = io_.InquireVariable< value_t >( func.getFunctionName() );
   if ( !varDoFData )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "ADIOS2 variable '" << func.getFunctionName() << "' does not exists!" );
      WALBERLA_ABORT( "So geht das nicht!" );
   }
   else
   {
      typename adios2::Variable< value_t >::Span dofData = engine_.Put< value_t >( varDoFData );
      AdiosWriter::StreamAccessBuffer< value_t > dataStream( dofData, varDoFData.Count() );
      VTKP2Writer::writeP2FunctionData( !storage_->hasGlobalCells(), dataStream, func, storage_, level_ );
   }
}

template < typename value_t >
void AdiosWriterForP2::scheduleVectorFunctionForExport( const P2VectorFunction< value_t >& func )
{
   adios2::Variable< value_t > varDoFData = io_.InquireVariable< value_t >( func.getFunctionName() );
   if ( !varDoFData )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "ADIOS2 variable '" << func.getFunctionName() << "' does not exists!" );
      WALBERLA_ABORT( "So geht das nicht!" );
   }
   else
   {
      typename adios2::Variable< value_t >::Span dofData = engine_.Put< value_t >( varDoFData );
      AdiosWriter::StreamAccessBuffer< value_t > dataStream( dofData, varDoFData.Count() );
      VTKP2Writer::writeP2VectorFunctionData( !storage_->hasGlobalCells(), dataStream, func, storage_, level_ );
   }
}

} // namespace hyteg
