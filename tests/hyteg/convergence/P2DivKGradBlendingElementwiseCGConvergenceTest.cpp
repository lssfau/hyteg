/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
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
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg//geometry/IcosahedralShellMap.hpp"
#include "hyteg//geometry/AnnulusMap.hpp"
//#include "hyteg/petsc/PETScMinResSolver.hpp"
//#include "hyteg/petsc/PETScLUSolver.hpp"
//#include "hyteg/petsc/PETScManager.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void P2DivKGradBlendingElementwiseCGTest()
{
   const uint_t level = 2;
   MeshInfo meshInfo = MeshInfo::meshSphericalShell( 3, 3, 0.5, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   IcosahedralShellMap::setMap( setupStorage );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "P2DivKGradBlendingElementwiseCGTest_domain" );

   hyteg::P2ElementwiseDivKGradBlendingOperator L( storage, level, level );

   hyteg::P2Function< real_t > r( "r", storage, level, level );
   hyteg::P2Function< real_t > f( "f", storage, level, level );
   hyteg::P2Function< real_t > u( "u", storage, level, level );
   hyteg::P2Function< real_t > u_exact( "u_exact", storage, level, level );
   hyteg::P2Function< real_t > err( "err", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const hyteg::Point3D& ) > rhs   = []( const hyteg::Point3D& ) { return 0; };
   std::function< real_t( const hyteg::Point3D& ) > ones  = []( const hyteg::Point3D& ) { return 1.0; };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   communication::syncP2FunctionBetweenPrimitives( u, level ); 


   auto solver = hyteg::CGSolver< hyteg::P2ElementwiseDivKGradBlendingOperator >( storage, level, level );
   solver.setPrintInfo( true );
   // auto solver = std::make_shared< PETScLUSolver< hyteg::P2ElementwiseBlendingLaplaceOperator > >( storage, level );

   err.assign( {1.0, -1.0}, {u, u_exact}, level );

   hyteg::VTKOutput vtkOutput( "../../output", "P2DivKGradBlendingElementwiseCGTest", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( f );
   vtkOutput.add( r );
   vtkOutput.add( err );
   vtkOutput.write( level, 0 );

   solver.solve( L, u, f, level );

   err.assign( {1.0, -1.0}, {u, u_exact}, level );

   vtkOutput.write( level, 1 );

   const auto npoints      = real_c( numberOfGlobalDoFs<P2FunctionTag>( *storage, level ) );
   const real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / real_c( npoints ) );


   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err << " (level " << level << ", " << "unknowns incl boundary: " << npoints << ")" );
   WALBERLA_CHECK_LESS( discr_l2_err, 3.5e-05 );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

//   hyteg::PETScManager manager( &argc, &argv );
   hyteg::P2DivKGradBlendingElementwiseCGTest();
}