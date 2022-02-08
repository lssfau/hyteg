/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "core/Environment.h"
#include "core/math/Constants.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScExportOperatorMatrix.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

void testMass()
{
   using namespace dg;

   // MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 1, 1 );
   MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );
   // MeshInfo                            meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t minLevel = 2;
   const uint_t maxLevel = 5;

   // std::function< real_t( const Point3D& ) > f = []( const Point3D& x ) { return sqrt( x[0] * x[0] + x[1] * x[1] ) < 0.3 ? 1 : 0; };
   // std::function< real_t( const Point3D& ) > f = []( const Point3D& x ) { return tanh( 8 * x[0] ); };
   std::function< real_t( const Point3D& ) > f = []( const Point3D& x ) { return sin( 2 * pi * x[0] ); };

   auto basis = std::make_shared< DGBasisLinearLagrange_Example >();

   DGFunction< real_t > u( "u", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel, basis, 1 );

   tmp.evaluateLinearFunctional( f, maxLevel );

   auto       mass = std::make_shared< DGMassForm_Example >();
   DGOperator M( storage, minLevel, maxLevel, mass );

   std::string                     fileName = "/tmp/mass.m";
   PETScSparseMatrix< DGOperator > mat;
   DGFunction< idx_t >             numeratorSrc( "numeratorSrc", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< idx_t >             numeratorDst( "numeratorDst", storage, minLevel, maxLevel, basis, 1 );
   numeratorSrc.enumerate( maxLevel );
   numeratorDst.enumerate( maxLevel );
   mat.createMatrixFromOperator( M, maxLevel, numeratorSrc, numeratorDst );
   mat.print( fileName, false, PETSC_VIEWER_ASCII_MATLAB );

   PETScCGSolver< DGOperator > solver( storage, maxLevel, numeratorSrc );
   solver.solve( M, u, tmp, maxLevel );

   VTKOutput vtkOutput( "../../output", "DGSmokeTestMass", storage );

   vtkOutput.add( u );
   vtkOutput.add( tmp );

   vtkOutput.write( maxLevel );
}

void testDiffusion( uint_t maxLevel )
{
   using namespace dg;

   MeshInfo meshInfo = MeshInfo::meshFaceChain( 1 );

   auto ff = []( const Point3D& p ) {
      Matrix3r mat;
      mat( 0, 0 ) = .34;
      mat( 0, 1 ) = .678;
      mat( 0, 2 ) = 0;

      mat( 1, 0 ) = .282;
      mat( 1, 1 ) = .35353;
      mat( 1, 2 ) = 0;

      mat( 2, 0 ) = 0;
      mat( 2, 1 ) = 0;
      mat( 2, 2 ) = 0;

      return mat.mul( p );
   };

   // meshInfo.applyCoordinateMap( ff );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   const uint_t minLevel = 2;

   std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
      return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * ( x[0] + x[1] - 1 ) );
   };

   std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
      return 4 * pi * pi * ( -2 * sin( 4 * pi * ( x[0] + x[1] ) ) + sin( 4 * pi * x[0] ) + sin( 4 * pi * x[1] ) );
   };

   auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();
   auto laplaceForm = std::make_shared< DGDiffusionForm_Example >();
   auto massForm    = std::make_shared< DGMassForm_Example >();

   DGFunction< real_t > u( "u", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > f( "f", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > sol( "sol", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > tmp( "tmp", storage, minLevel, maxLevel, basis, 1 );
   DGFunction< real_t > err( "err", storage, minLevel, maxLevel, basis, 1 );

   DGFunction< idx_t > numerator( "numerator", storage, minLevel, maxLevel, basis, 1 );
   numerator.enumerate( maxLevel );

   DGOperator A( storage, minLevel, maxLevel, laplaceForm );
   DGOperator M( storage, minLevel, maxLevel, massForm );

   std::string                     fileName = "/tmp/diffusion.m";
   PETScSparseMatrix< DGOperator > mat;
   mat.createMatrixFromOperator( A, maxLevel, numerator, numerator );
   mat.print( fileName, false, PETSC_VIEWER_ASCII_MATLAB );

   f.evaluateLinearFunctional( rhsFunc, maxLevel );
   tmp.evaluateLinearFunctional( solFunc, maxLevel );

   PETScCGSolver< DGOperator > solverM( storage, maxLevel, numerator );
   solverM.solve( M, sol, tmp, maxLevel );
   PETScCGSolver< DGOperator > solverA( storage, maxLevel, numerator );
   solverA.solve( A, u, f, maxLevel );

   err.assign( { 1.0, -1.0 }, { u, sol }, maxLevel );
   auto discrL2 = sqrt( err.dotGlobal( err, maxLevel ) / real_c( numberOfGlobalDoFs( u, maxLevel ) ) );

   WALBERLA_LOG_INFO_ON_ROOT( "level " << maxLevel << ": L2 error " << discrL2 );

   VTKOutput vtkOutput( "../../output", "DGSmokeTestLaplace", storage );

   vtkOutput.add( u );
   vtkOutput.add( err );
   vtkOutput.add( sol );

   vtkOutput.write( maxLevel );
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   // hyteg::testMass();

   hyteg::testDiffusion( 3 );
   // hyteg::testDiffusion( 4 );
   // hyteg::testDiffusion( 5 );
   hyteg::testDiffusion( 6 );

   return EXIT_SUCCESS;
}
