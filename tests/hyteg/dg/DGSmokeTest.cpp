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

void test()
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

   VTKOutput vtkOutput( "../../output", "DGSmokeTest", storage );

   vtkOutput.add( u );
   vtkOutput.add( tmp );

   vtkOutput.write( maxLevel );
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   hyteg::test();

   return EXIT_SUCCESS;
}
