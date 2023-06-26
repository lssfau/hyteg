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

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGVectorFunction.hpp"
#include "hyteg/dgfunctionspace/DGVectorOperators.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

real_t testHomDirichlet( uint_t level, uint_t degree )
{
   using namespace dg;

   MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_16el.msh" );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   std::function< real_t( const Point3D& ) > solFunc0 = []( const Point3D& x ) {
     return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] );
   };

   std::function< real_t( const Point3D& ) > solFunc1 = []( const Point3D& x ) {
     return sin( 3 * pi * x[0] ) * sin( 3 * pi * x[1] );
   };

   std::function< real_t( const Point3D& ) > rhsFunc0 = []( const Point3D& x ) {
     return 8 * pi * pi * ( sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) );
   };

   std::function< real_t( const Point3D& ) > rhsFunc1 = []( const Point3D& x ) {
     return 2 * 9 * pi * pi * ( sin( 3 * pi * x[0] ) * sin( 3 * pi * x[1] ) );
   };

   auto basis       = std::make_shared< DGBasisLinearLagrange_Example >();

   DGVectorFunction< real_t > u( "u", storage, level, level, basis, degree );
   DGVectorFunction< real_t > f( "f", storage, level, level, basis, degree );
   DGVectorFunction< real_t > sol( "sol", storage, level, level, basis, degree );
   DGVectorFunction< real_t > tmp( "tmp", storage, level, level, basis, degree );
   DGVectorFunction< real_t > err( "err", storage, level, level, basis, degree );

   DGVectorFunction< idx_t > numerator( "numerator", storage, level, level, basis, degree );
   numerator.enumerate( level );

   DGVectorMassOperator< real_t >    M(storage, level, level);
   DGVectorLaplaceOperator< real_t > A(storage, level, level);

   f[0].evaluateLinearFunctional( rhsFunc0, level );
   f[1].evaluateLinearFunctional( rhsFunc1, level );

   tmp[0].evaluateLinearFunctional( solFunc0, level );
   tmp[1].evaluateLinearFunctional( solFunc1, level );

   PETScCGSolver< DGVectorMassOperator< real_t > > solverM( storage, level, numerator );
   solverM.solve( M, sol, tmp, level );

   PETScCGSolver< DGVectorLaplaceOperator< real_t > > solverA( storage, level, numerator, 1e-12, 1e-12, 10000 );
   solverA.solve( A, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   auto discrL2 = sqrt( err.dotGlobal( err, level ) / real_c( numberOfGlobalDoFs( u, level ) ) );

   WALBERLA_LOG_INFO(discrL2);

   VTKOutput vtk( "../../output/", "DGVectorPoisson2DConvergenceTest_testDirichlet", storage );
   vtk.add( u );
   vtk.add( sol );
   vtk.add( err );
   vtk.add( f );
   vtk.write( level );

   return discrL2;
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   {
      const uint_t degree = 1;
      const uint_t level = 4;
      auto err_level= hyteg::testHomDirichlet( level, degree );
      auto err_level_next = hyteg::testHomDirichlet( level+1, degree );
      auto rate = err_level / err_level_next;

      WALBERLA_CHECK_LESS(rate, 4.1);
      WALBERLA_CHECK_GREATER(rate, 3.9);
   }

   return EXIT_SUCCESS;
}
