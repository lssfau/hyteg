/*
* Copyright (c) 2017-2024 Nils Kohl.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/math/Constants.h"
#include "core/timing/all.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dg1functionspace/DG1Function.hpp"
#include "hyteg/egfunctionspace/EGFunction.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg {

using walberla::math::pi;

const std::string P0_SCALAR   = "p0-scalar";
const std::string P1_SCALAR   = "p1-scalar";
const std::string P2_SCALAR   = "p2-scalar";
const std::string DG1_SCALAR  = "dg1-scalar";
const std::string P1_VECTOR   = "p1-vector";
const std::string P2_VECTOR   = "p2-vector";
const std::string EG_VECTOR   = "eg-vector";
const std::string N1E1_VECTOR = "n1e1-vector";

const std::vector< std::string > functionStrings =
    { P0_SCALAR, P1_SCALAR, P2_SCALAR, DG1_SCALAR, P1_VECTOR, P2_VECTOR, EG_VECTOR, N1E1_VECTOR };

const auto gradient = []( const Point3D& x ) { return x[0] + x[1] + x[2]; };

template < typename FunctionType >
static void interpolateAndWriteVTK( const std::string&                         directory,
                                    uint_t                                     level,
                                    const std::shared_ptr< PrimitiveStorage >& storage,
                                    const FunctionType&                        function )
{
   function.interpolate( gradient, level );
   VTKOutput vtkOutput( directory, function.getFunctionName(), storage );
   vtkOutput.add( function );
   vtkOutput.write( level );
}

static void exportFunctions( const std::string& functionType, uint_t level, const std::shared_ptr< PrimitiveStorage >& storage )
{
   uint_t minLevel = level;
   uint_t maxLevel = level;

   std::string directory = "VTKOutputPictureNormTest-Fresh";

   std::string functionName = "u";

   if ( functionType == P0_SCALAR )
   {
      interpolateAndWriteVTK( directory, level, storage, P0Function< real_t >( functionName, storage, minLevel, maxLevel ) );
   }
   else if ( functionType == P1_SCALAR )
   {
      interpolateAndWriteVTK( directory, level, storage, P1Function< real_t >( functionName, storage, minLevel, maxLevel ) );
   }
   else if ( functionType == P2_SCALAR )
   {
      interpolateAndWriteVTK( directory, level, storage, P2Function< real_t >( functionName, storage, minLevel, maxLevel ) );
   }
   else if ( functionType == DG1_SCALAR )
   {
      interpolateAndWriteVTK( directory, level, storage, DG1Function< real_t >( functionName, storage, minLevel, maxLevel ) );
   }
   else if ( functionType == P1_VECTOR )
   {
      interpolateAndWriteVTK(
          directory, level, storage, P1VectorFunction< real_t >( functionName, storage, minLevel, maxLevel ) );
   }
   else if ( functionType == P2_VECTOR )
   {
      interpolateAndWriteVTK(
          directory, level, storage, P2VectorFunction< real_t >( functionName, storage, minLevel, maxLevel ) );
   }
   else if ( functionType == EG_VECTOR )
   {
      interpolateAndWriteVTK( directory, level, storage, EGFunction< real_t >( functionName, storage, minLevel, maxLevel ) );
   }
   else if ( functionType == N1E1_VECTOR )
   {
      WALBERLA_CHECK( storage->hasGlobalCells(), "Cannot create N1E1 function if the domain dimension is not 3." )
      interpolateAndWriteVTK(
          directory, level, storage, n1e1::N1E1VectorFunction< real_t >( functionName, storage, minLevel, maxLevel ) );
   }
}

static std::shared_ptr< PrimitiveStorage > storageFromMesh( const std::string& meshFile )
{
   WALBERLA_LOG_INFO_ON_ROOT( "Opening mesh file " << prependHyTeGMeshDir( meshFile ) );
   MeshInfo                            mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( meshFile ) );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );
   return storage;
}

static std::shared_ptr< PrimitiveStorage > storageCurvedShell( int dim )
{
   MeshInfo mesh = MeshInfo::emptyMeshInfo();

   if ( dim == 2 )
   {
      mesh = MeshInfo::meshAnnulus( 0.5, 1.0, MeshInfo::CRISS, 12, 2 );
   }
   else if ( dim == 3 )
   {
      mesh = MeshInfo::meshSphericalShell( 2, 2, 0.5, 1.0 );
   }
   else
   {
      WALBERLA_ABORT( "Invalid dimension." )
   }

   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   if ( dim == 2 )
   {
      AnnulusMap::setMap( setupStorage );
   }
   else
   {
      IcosahedralShellMap::setMap( setupStorage );
   }

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   return storage;
}

static std::shared_ptr< PrimitiveStorage > storageParametric( int dim, uint_t level, uint_t mappingDegree )
{
   MeshInfo mesh = MeshInfo::emptyMeshInfo();

   if ( dim == 2 )
   {
      mesh = MeshInfo::meshUnitSquare( 1 );
   }
   else if ( dim == 3 )
   {
      mesh = MeshInfo::meshCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 2 ), 1, 1, 1 );
   }
   else
   {
      WALBERLA_ABORT( "Invalid dimension." )
   }

   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   std::vector< std::function< real_t( const Point3D& ) > > curveFunc = {
       [&]( const Point3D& x ) { return x[0] + 0.1 * sin( 2 * pi * 1.5 * x[1] ); },
       [&]( const Point3D& x ) { return x[1] + 0.1 * sin( 2 * pi * x[0] ); },
       [&]( const Point3D& x ) { return x[2] + 0.1 * sin( 2 * pi * 1.1 * x[0] ); },
   };

   const auto microMesh = std::make_shared< micromesh::MicroMesh >( storage, level, level, mappingDegree, dim );

   micromesh::interpolateAndCommunicate( *microMesh, curveFunc, level );
   storage->setMicroMesh( microMesh );

   return storage;
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   const uint_t level = 2;

   const auto functionFlag      = "--function";
   const auto meshFlag          = "--mesh";
   const auto blendingShellFlag = "--curved-shell";
   const auto parametricFlag    = "--parametric";

   std::string                                functionParam;
   std::shared_ptr< hyteg::PrimitiveStorage > storage;

   if ( std::find_if( argv, argv + argc, []( const std::string& s ) { return s == "--help" || s == "-h"; } ) != argv + argc )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Generated VTK output of various functions to be used for VTK output regression testing in the "
                                 "'picture-norm'." )
      WALBERLA_LOG_INFO_ON_ROOT( "Usage:" )
      WALBERLA_LOG_INFO_ON_ROOT( "  -h/--help:             show this message and exit" )
      WALBERLA_LOG_INFO_ON_ROOT( "  " << functionFlag << " <function>: function type to output, any of" )
      for ( const auto& s : hyteg::functionStrings )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "      " << s )
      }
      WALBERLA_LOG_INFO_ON_ROOT( "  " << meshFlag << " <mesh-file>:    file path to a mesh file, e.g., \"2D/penta_5el.msh\"" )
      WALBERLA_LOG_INFO_ON_ROOT( "  " << blendingShellFlag
                                      << " <dim>:  generates a blended annulus (dim == 2) or thick spherical shell (dim == 3)" )
      WALBERLA_LOG_INFO_ON_ROOT( "  " << parametricFlag << " <dim>:  generates a parametrically (P1) curved domain" )
      return EXIT_SUCCESS;
   }

   // KISS parameter search
   for ( int i = 0; i < argc; i++ )
   {
      const std::string prm( argv[i] );
      if ( prm == functionFlag )
      {
         WALBERLA_CHECK_GREATER( argc, i + 1, "Found " << functionFlag << "but no argument." )
         functionParam = std::string( argv[i + 1] );

         WALBERLA_CHECK( std::find( hyteg::functionStrings.begin(), hyteg::functionStrings.end(), functionParam ) !=
                             hyteg::functionStrings.end(),
                         "Invalid function. Run --help for list." )

         WALBERLA_LOG_INFO_ON_ROOT( "Function: " << functionParam )
         break;
      }
   }

   for ( int i = 0; i < argc; i++ )
   {
      const std::string prm( argv[i] );
      if ( prm == meshFlag )
      {
         if ( storage )
         {
            WALBERLA_ABORT( "Decide: either mesh file, shell, or parametric!" )
         }
         WALBERLA_CHECK_GREATER( argc, i + 1, "Found " << meshFlag << "but no argument." )
         const auto meshFile = std::string( argv[i + 1] );
         WALBERLA_LOG_INFO_ON_ROOT( "Mesh file: " << meshFile )
         storage = hyteg::storageFromMesh( meshFile );
      }

      if ( prm == blendingShellFlag )
      {
         if ( storage )
         {
            WALBERLA_ABORT( "Decide: either mesh file, shell, or parametric!" )
         }
         WALBERLA_CHECK_GREATER( argc, i + 1, "Found " << blendingShellFlag << "but no argument." )
         const auto dim = std::stoi( argv[i + 1] );
         WALBERLA_LOG_INFO_ON_ROOT( "Curved shell, dim: " << dim )
         storage = hyteg::storageCurvedShell( dim );
      }

      if ( prm == parametricFlag )
      {
         if ( storage )
         {
            WALBERLA_ABORT( "Decide: either mesh file, shell, or parametric!" )
         }

         WALBERLA_LOG_INFO_ON_ROOT( functionParam );

         if ( functionParam != hyteg::P1_SCALAR && functionParam != hyteg::P1_VECTOR && functionParam != hyteg::P2_SCALAR &&
              functionParam != hyteg::P2_VECTOR )
         {
            WALBERLA_ABORT(
                "Parametric mappings only supported here for P1 and P2 scalar + vector functions (but could be extended probably)." )
         }

         WALBERLA_CHECK_GREATER( argc, i + 1, "Found " << parametricFlag << "but no argument." )
         const auto dim = std::stoi( argv[i + 1] );
         WALBERLA_LOG_INFO_ON_ROOT( "Parametric map, dim: " << dim )

         uint_t mappingDegree = 0;
         if ( functionParam == hyteg::P1_SCALAR || functionParam == hyteg::P1_VECTOR )
         {
            mappingDegree = 1;
         }
         else if ( functionParam == hyteg::P2_SCALAR || functionParam == hyteg::P2_VECTOR )
         {
            mappingDegree = 2;
         }

         storage = hyteg::storageParametric( dim, level, mappingDegree );
      }
   }

   hyteg::exportFunctions( functionParam, level, storage );

   return EXIT_SUCCESS;
}
