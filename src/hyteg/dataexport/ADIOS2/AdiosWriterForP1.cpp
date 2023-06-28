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
#include "hyteg/dataexport/ADIOS2/AdiosWriterForP1.hpp"

#include <adios2.h>

#include "hyteg/dataexport/ADIOS2/AdiosHelperFunctions.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/FEFunctionRegistry.hpp"
#include "hyteg/dataexport/FEFunctionWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKMeshWriter.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

AdiosWriterForP1::AdiosWriterForP1( std::string&                               filePath,
                                    std::string&                               fileBaseName,
                                    std::string&                               configFile,
                                    uint_t                                     level,
                                    const std::shared_ptr< PrimitiveStorage >& storage )
: storage_( storage )
, level_( level )
{
   // setup central ADIOS2 interface object
#ifdef WALBERLA_BUILD_WITH_MPI
   if ( configFile.empty() )
   {
      adios_ = adios2::ADIOS( walberla::MPIManager::instance()->comm() );
   }
   else
   {
      adios_ = adios2::ADIOS( configFile, walberla::MPIManager::instance()->comm() );
   }
#else
   if ( configFile.empty() )
   {
      adios_ = adios2::ADIOS;
   }
   else
   {
      adios_ = adios2::ADIOS( configFile );
   }
#endif

   // create our writer
   io_ = adios_.DeclareIO( "AdiosWriter" );

   // we will use the BP format, but the most recent "BP5" format
   // is unsupported by ParaView 5.11.1
   io_.SetEngine( "BP4" );

   // create the name of the output file
   std::stringstream tag;
   tag << "-P1_level" << level_ << ".bp";
   fileName_ = filePath + "/" + fileBaseName + tag.str();
}

void AdiosWriterForP1::write( const FEFunctionRegistry& registry, uint_t timestep )
{
   // Assemble list of names of node-centred and cell-centred funcions (latter is empty)
   std::vector< std::string > emptyList;
   std::vector< std::string > p1FunctionList;
   registry.extractFunctionNames( p1FunctionList, functionTraits::P1_FUNCTION );
   registry.extractFunctionNames( p1FunctionList, functionTraits::P1_VECTOR_FUNCTION );

   // ------------------
   //  First write only
   // ------------------
   if ( !firstWriteCompleted_ )
   {
      // open output file
      engine_ = io_.Open( fileName_, adios2::Mode::Write );

      // export the mesh information
      engine_.BeginStep();
      writeMesh( p1FunctionList );

      // for testing
      adios2::Variable< real_t > varTime = io_.DefineVariable< real_t >( "TIME" );
      engine_.Put( varTime, real_c( 0 ) );

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

   for ( const auto& func : registry.getP1Functions().getFunctions< real_t >() )
   {
      scheduleScalarFunctionForExport( func );
   }

   for ( const auto& func : registry.getP1VectorFunctions().getFunctions< real_t >() )
   {
      scheduleVectorFunctionForExport( func );
   }

   engine_.EndStep();

   WALBERLA_LOG_WARNING_ON_ROOT( "AdiosWriterForP1::write() not fully functional, yet!" );
}

void AdiosWriterForP1::writeMesh( const std::vector< std::string >& p1FunctionList ) const
{
   // store the generic XML file as attribute as part of the binary data
   std::vector< std::string > emptyList;
   std::string                vtkMetaInfo = adiosHelpers::generateVTKMetaInfo( p1FunctionList, emptyList );
   io_.DefineAttribute( "vtk.xml", vtkMetaInfo, "", "" );

   // determine basic entity counts
   const uint_t dim         = storage_->hasGlobalCells() ? 3u : 2u;
   uint_t       numElements = 0u;
   uint_t       numVertices = 0u;
   if ( dim == 3 )
   {
      numVertices = storage_->getNumberOfLocalCells() * levelinfo::num_microvertices_per_cell( level_ );
      numElements = storage_->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( level_ );
   }
   else
   {
      numVertices = storage_->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level_ );
      numElements = storage_->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level_ );
   }
   WALBERLA_LOG_INFO_ON_ROOT( "writeMesh: numVertices = " << numVertices );
   WALBERLA_LOG_INFO_ON_ROOT( "writeMesh: numElements = " << numElements );

   // store entity counts as scalar variables
   adios2::Variable< uint32_t > varNumberOfNodes =
       io_.DefineVariable< uint32_t >( "NumberOfVertices", { adios2::LocalValueDim } );
   adios2::Variable< uint32_t > varNumberOfFaces =
       io_.DefineVariable< uint32_t >( "NumberOfElements", { adios2::LocalValueDim } );

   engine_.Put( varNumberOfNodes, static_cast< uint32_t >( numVertices ) );
   engine_.Put( varNumberOfFaces, static_cast< uint32_t >( numElements ) );

   // store the type of elements in the mesh as scalar variable
   // 5u = "Triangle", 10u = "Tetrahedron"
   using elementMarker_t                               = uint32_t;
   elementMarker_t                     elementMarker   = dim == 2 ? 5u : 10u;
   adios2::Variable< elementMarker_t > varElementTypes = io_.DefineVariable< uint32_t >( "types" );
   engine_.Put( varElementTypes, elementMarker );

   // store node coordinates (always x,y,z even for 2D!)
   adios2::Variable< real_t >                varVertices = io_.DefineVariable< real_t >( "vertices", {}, {}, { numVertices, 3 } );
   adios2::Variable< real_t >::Span          vertices    = engine_.Put< real_t >( varVertices );
   AdiosWriter::StreamAccessBuffer< real_t > vertexStream( vertices, varVertices.Count() );
   VTKMeshWriter::writePointsForMicroVertices( dim == 2, vertexStream, storage_, level_, false );

   // store element connectivity
   adios2::Variable< uint64_t > varConnectivity =
       io_.DefineVariable< uint64_t >( "connectivity", {}, {}, { numElements, dim + 2 } );
   adios2::Variable< uint64_t >::Span connectivity = engine_.Put< uint64_t >( varConnectivity );

   if ( dim == 3 )
   {
      AdiosWriter::StreamAccessBuffer< uint64_t, 5 > connectivityStream( connectivity, varConnectivity.Count() );
      VTKMeshWriter::writeElementNodeAssociationP1Tetrahedrons(
          connectivityStream, storage_, levelinfo::num_microvertices_per_edge( level_ ), false );
   }
   else
   {
      AdiosWriter::StreamAccessBuffer< uint64_t, 4 > connectivityStream( connectivity, varConnectivity.Count() );
      VTKMeshWriter::writeElementNodeAssociationP1Triangles(
          connectivityStream, storage_, levelinfo::num_microvertices_per_edge( level_ ), false );
   }
}

template < typename value_t >
void AdiosWriterForP1::defineVariables( const FEFunctionRegistry& registry ) const
{
   auto checkAndDefine = [this]( std::string funcName, uint_t funcDim ) {
      adios2::Variable< real_t > varDoFData = io_.InquireVariable< value_t >( funcName );
      if ( !varDoFData )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Defining ADIOS2 variable '" << funcName << "'" );
         uint_t extent{ 0 };
         if ( !storage_->hasGlobalCells() )
         {
            extent = levelinfo::num_microvertices_per_face( level_ ) * storage_->getNumberOfLocalFaces();
         }
         else
         {
            extent = levelinfo::num_microvertices_per_cell( level_ ) * storage_->getNumberOfLocalCells();
         }
         WALBERLA_LOG_INFO_ON_ROOT( " -> extent = " << extent );
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

   std::vector< P1Function< value_t > > p1Funcs = registry.getP1Functions().getFunctions< value_t >();
   for ( const auto& func : p1Funcs )
   {
      checkAndDefine( func.getFunctionName(), func.getDimension() );
   }

   std::vector< P1VectorFunction< value_t > > p1VecFuncs = registry.getP1VectorFunctions().getFunctions< value_t >();
   for ( const auto& func : p1VecFuncs )
   {
      checkAndDefine( func.getFunctionName(), func.getDimension() );
   }
}

template < typename value_t >
void AdiosWriterForP1::scheduleScalarFunctionForExport( const P1Function< value_t >& func ) const
{
   adios2::Variable< real_t > varDoFData = io_.InquireVariable< value_t >( func.getFunctionName() );
   if ( !varDoFData )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "ADIOS2 variable '" << func.getFunctionName() << "' does not exists!" );
      WALBERLA_ABORT( "So geht das nicht!" );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Putting ADIOS2 variable '" << func.getFunctionName() << "'" );
      adios2::Variable< real_t >::Span dofData = engine_.Put< value_t >( varDoFData );
      if ( !storage_->hasGlobalCells() )
      {
         uint_t blockSize = levelinfo::num_microvertices_per_face( level_ );
         uint_t currentPos{ 0 };
         for ( const auto& it : func.getStorage()->getFaces() )
         {
            const Face& face = *it.second;
            // assuming here that the Span is a contiguous memory buffer
            memcpy(
                &dofData[currentPos], face.getData( func.getFaceDataID() )->getPointer( level_ ), blockSize * sizeof( value_t ) );
            currentPos += blockSize;
         }
         WALBERLA_LOG_INFO_ON_ROOT( "Put " << currentPos << " data items!" );
      }
      else
      {
         uint_t blockSize = levelinfo::num_microvertices_per_cell( level_ );
         uint_t currentPos{ 0 };
         for ( const auto& it : func.getStorage()->getCells() )
         {
            const Cell& cell = *it.second;
            // assuming here that the Span is a contiguous memory buffer
            memcpy(
                &dofData[currentPos], cell.getData( func.getCellDataID() )->getPointer( level_ ), blockSize * sizeof( value_t ) );
            currentPos += blockSize;
         }
      }
   }
}

template < typename value_t >
void AdiosWriterForP1::scheduleVectorFunctionForExport( const P1VectorFunction< value_t >& func ) const
{
   adios2::Variable< real_t > varDoFData = io_.InquireVariable< value_t >( func.getFunctionName() );
   if ( !varDoFData )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "ADIOS2 variable '" << func.getFunctionName() << "' does not exists!" );
      WALBERLA_ABORT( "So geht das nicht!" );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Putting ADIOS2 variable '" << func.getFunctionName() << "'" );
      adios2::Variable< real_t >::Span dofData = engine_.Put< value_t >( varDoFData );
      if ( !storage_->hasGlobalCells() )
      {
         uint_t currentPos{ 0 };
         for ( const auto& it : storage_->getFaces() )
         {
            const Face& face = *it.second;
            for ( uint_t i = 0; i < levelinfo::num_microvertices_per_face( level_ ); ++i )
            {
               for ( uint_t idx = 0; idx < func.getDimension(); ++idx )
               {
                  dofData[currentPos++] = face.getData( func[idx].getFaceDataID() )->getPointer( level_ )[i];
               }
            }
         }
      }
      else
      {
        uint_t      currentPos{ 0 };
         for ( const auto& it : storage_->getCells() )
         {
            const Cell& cell      = *it.second;
            const auto  cellData0 = cell.getData( func[0].getCellDataID() )->getPointer( level_ );
            const auto  cellData1 = cell.getData( func[1].getCellDataID() )->getPointer( level_ );
            const auto  cellData2 = cell.getData( func[2].getCellDataID() )->getPointer( level_ );

            for ( const auto& idxIt : vertexdof::macrocell::Iterator( level_ ) )
            {
               dofData[currentPos++] = cellData0[vertexdof::macrocell::index( level_, idxIt.x(), idxIt.y(), idxIt.z() )];
               dofData[currentPos++] = cellData1[vertexdof::macrocell::index( level_, idxIt.x(), idxIt.y(), idxIt.z() )];
               dofData[currentPos++] = cellData2[vertexdof::macrocell::index( level_, idxIt.x(), idxIt.y(), idxIt.z() )];
            }
         }
      }
   }
}

} // namespace hyteg
