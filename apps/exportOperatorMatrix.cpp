/*
 * Copyright (c) 2017-2025 Christoph Schwarzmeier, Dominik Thoennes, Marcus Mohr.
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
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>
#include <sstream>
#include <string>

// Primitive management
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

// Function spaces
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

// P1 Operators
#include "constant_stencil_operator/P1ConstantOperator.hpp"

// P2 Operators
#include "constant_stencil_operator/P2ConstantOperator.hpp"

// P2PlusBubble Operators
#include "hyteg_operators/operators/diffusion/P2PlusBubbleElementwiseDiffusion.hpp"

// Mixed Operators
#include "hyteg/composites/CCRStokesOperator.hpp"

#include "mixed_operator/P2P1TaylorHoodStokesOperator.hpp"

// PETSc interface
#include "hyteg/petsc/PETScExportOperatorMatrix.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

typedef enum
{
   CCRSTOKES,
   P1MASS,
   P1LAPLACE,
   P2MASS,
   P2LAPLACE,
   P2PLUSBUBBLELAPLACE,
   P2P1STOKES,
   P2EDGEMASS
} operatorTag;

typedef struct
{
   operatorTag oprEnum;
   std::string matName;
   bool        elimDirichletBC;
} oprInfo;

std::map< std::string, oprInfo > oprMap = { { "CCRStokes", { CCRSTOKES, "StokesOpCCR", true } },
                                            { "P1Mass", { P1MASS, "MassOpP1", false } },
                                            { "P1Diff", { P1LAPLACE, "DiffOpP1", true } },
                                            { "P2Mass", { P2MASS, "MassOpP2", false } },
                                            { "P2Diff", { P2LAPLACE, "DiffOpP2", true } },
                                            { "P2PlusBubbleDiff", { P2PLUSBUBBLELAPLACE, "DiffOpP2PlusBubble", true } },
                                            { "P2P1Stokes", { P2P1STOKES, "StokesOpP2P1", true } },
                                            { "P2EdgeMass", { P2EDGEMASS, "MassOpP2_EdgeDoFs", false } } };

void showUsage()
{
   std::stringstream mesg;

   mesg << "Please specify the following two parameters in the given order:\n\n"
        << "  <level>               on which level do you want the operator to be set up?\n"
        << "  <operator>            which operator do you want to export?\n"
        << "  <format>              either Matlab or MatrixMarket\n\n"
        << "Optionally also give\n\n"
        << "  <name of Gmsh file>   if none is given unit square will be meshed with two triangles\n\n"
        << " Choices available for <operator> are\n";

   for ( auto it = oprMap.begin(); it != oprMap.end(); ++it )
   {
      mesg << " - " << it->first << "\n";
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" << mesg.str() );
}

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );

   WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "--------- Matrix Export of Operators ---------" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "----------------------------------------------" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "" );

   // Process command-line
   if ( argc < 4 || argc > 5 )
   {
      showUsage();
      WALBERLA_ABORT( "\n" );
   }

   bool useMeshFile = ( argc == 5 ) ? true : false;

   uint_t      level   = static_cast< uint_t >( std::stoul( argv[1] ) );
   std::string oprName = ( argv[2] );
   if ( oprMap.find( oprName ) == oprMap.end() )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Sorry, but '" << oprName << "' does not seem to be a valid choice." );
      WALBERLA_ABORT( "\n" );
   }
   operatorTag oprTag  = oprMap.at( oprName ).oprEnum;
   std::string matName = oprMap.at( oprName ).matName;
   bool        elim    = oprMap[oprName].elimDirichletBC;

   // determine output format and filename postfix
   PetscViewerFormat format   = PETSC_VIEWER_ASCII_MATLAB;
   std::string       fileName = oprName;
   std::string       fmtOpt   = argv[3];
   if ( fmtOpt == "Matlab" )
   {
      format = PETSC_VIEWER_ASCII_MATLAB;
      fileName.append( ".m" );
   }
   else if ( fmtOpt == "MatrixMarket" )
   {
      format = PETSC_VIEWER_ASCII_MATRIXMARKET;
      fileName.append( ".mtx" );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Format option '" << fmtOpt << "' not valid!\n" );
      showUsage();
      WALBERLA_ABORT( "\n" );
   }

   // Mesh generation
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();
   if ( useMeshFile )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Generating mesh from file " << argv[4] );
      meshInfo = MeshInfo::fromGmshFile( argv[4] );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Generating criss mesh on unit square" );
      Point2D cornerLL( 0.0, 0.0 );
      Point2D cornerUR( 1.0, 1.0 );
      meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CRISS, 1, 1 );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Should elimination of Dirichlet DoFs be symmetric?
   bool symm = true;

   // Should petsc::exportOperator() be verbose?
   bool verb = true;

   // Operator creation and export
   switch ( oprTag )
   {
      // --------------
      //  P1 operators
      // --------------
   case P1LAPLACE: {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Laplace operator for P1 elements" );
      hyteg::P1ConstantLaplaceOperator opr( storage, level, level );
      petsc::exportOperator< P1ConstantLaplaceOperator >( opr, fileName, matName, format, storage, level, elim, symm, verb );
   }
   break;

   case P1MASS: {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Mass operator for P1 elements" );
      hyteg::P1ConstantMassOperator opr( storage, level, level );
      petsc::exportOperator< P1ConstantMassOperator >( opr, fileName, matName, format, storage, level, elim, symm, verb );
   }
   break;

      // --------------
      //  P2 operators
      // --------------
   case P2LAPLACE: {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Laplace operator for P2 elements" );
      hyteg::P2ConstantLaplaceOperator opr( storage, level, level );
      petsc::exportOperator< P2ConstantLaplaceOperator >( opr, fileName, matName, format, storage, level, elim, symm, verb );
   }
   break;

   case P2MASS: {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Mass operator for P2 elements" );
      hyteg::P2ConstantMassOperator opr( storage, level, level );
      petsc::exportOperator< P2ConstantMassOperator >( opr, fileName, matName, format, storage, level, elim, symm, verb );
   }
   break;

   case P2EDGEMASS: {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Mass operator for P2 elements (EdgeDoFs only)" );
      typedef EdgeDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
                            P2EdgeDoFMassOperator;
      P2EdgeDoFMassOperator opr( storage, level, level );
      petsc::exportOperator< P2EdgeDoFMassOperator >( opr, fileName, matName, format, storage, level, elim, symm, verb );
   }
   break;

      // ------------------------
      //  P2PlusBubble operators
      // ------------------------
   case P2PLUSBUBBLELAPLACE: {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Laplace operator for P2PlusBubble elements" );
      using P2PlusBubbleElementwiseDiffusion = hyteg::operatorgeneration::P2PlusBubbleElementwiseDiffusion;
      P2PlusBubbleElementwiseDiffusion opr( storage, level, level );
      petsc::exportOperator< P2PlusBubbleElementwiseDiffusion >( opr, fileName, matName, format, storage, level, elim, symm, verb );
   }
   break;

      // -----------------
      //  Mixed operators
      // -----------------
   case P2P1STOKES: {
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting Stokes Operator for P2-P1 element" );
      hyteg::P2P1TaylorHoodStokesOperator opr( storage, level, level );
      petsc::exportOperator< P2P1TaylorHoodStokesOperator >( opr, fileName, matName, format, storage, level, elim, symm, verb );
   }
   break;

   case CCRSTOKES: {
     WALBERLA_LOG_INFO_ON_ROOT( "Exporting Stokes Operator for Conforming Crouzeix-Raviart element" );
     hyteg::CCRStokesOperator opr( storage, level, level );
     petsc::exportOperator< CCRStokesOperator >( opr, fileName, matName, format, storage, level, elim, symm, verb );
   }
   break;
   }

   WALBERLA_LOG_DEVEL_ON_ROOT( "----------------------------------------------" );

   return 0;
}
