/*
* Copyright (c) 2017-2024 Andreas Burkhart, Nils Kohl.
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
#include <core/config/Config.h>
#include <core/mpi/MPIManager.h>
#include <core/math/Constants.h>

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/solvertemplates/StokesFSGMGSolverTemplate.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "hyteg/numerictools/L2Space.hpp"

#include "terraneo/operators/P2P1StokesOperatorWithProjection.hpp"
#include "terraneo/operators/P2StokesABlockWithProjection.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

using namespace hyteg;

using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   real_t   rMin = 0.5, rMax = 1.0;
   MeshInfo meshInfo = MeshInfo::meshSphericalShell( 3, 2, 0.5, 1.0 );

   SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   IcosahedralShellMap::setMap( setupStorage );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   BoundaryCondition bcVelocity;
   bcVelocity.createDirichletBC( "DirichletOuter", {MeshInfo::hollowFlag::flagOuterBoundary} );
   // bcVelocity.createDirichletBC( "DirichletInner", {MeshInfo::hollowFlag::flagInnerBoundary} );
   bcVelocity.createFreeslipBC( "FreeslipInner", {MeshInfo::hollowFlag::flagInnerBoundary} );

   uint_t minLevel = 2U, maxLevel = 3U;

   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > fStrong( "fStrong", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > res( "res", storage, minLevel, maxLevel, bcVelocity );

   P2Function< real_t > mu( "mu", storage, minLevel, maxLevel );
   P1Function< real_t > muInv( "muInv", storage, minLevel, maxLevel );

   uint_t lMax = 5U;

   terraneo::SphericalHarmonicsTool sphTool( lMax );

   // Some coefficient and its inverse.
   auto visc    = [&]( const Point3D& x ) { return 2 + std::sin( x[0] ) * std::cos( x[1] ); };
   // auto visc    = [&]( const Point3D& x ) { return 1.0; };
   auto viscInv = [&]( const Point3D& x ) { return 1.0 / visc( x ); };
   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      mu.interpolate( visc, level, All );
      muInv.interpolate( viscInv, level, All );
   }

   std::function< void( const Point3D&, Point3D& ) > normalsFS = []( const Point3D& x, Point3D& n ) { n = -x / x.norm(); };

   real_t                                    rhsScale = 1.0;
   std::function< real_t( const Point3D& ) > fRHS     = [&]( const Point3D& x ) {
      real_t r = x.norm();
      return rhsScale * std::sin( walberla::math::pi * ( r - rMin ) / ( rMax - rMin ) ) *
             sphTool.shconvert_eval( 4, 2, x[0], x[1], x[2] );
   };

   fStrong.uvw().interpolate( fRHS, maxLevel, All );

   auto projectionOperator = std::make_shared< P2ProjectNormalOperator >( storage, minLevel, maxLevel, normalsFS );

   auto stokesOperatorFS = std::make_shared< P2P1StokesFullIcosahedralShellMapOperatorFS >(
       storage, minLevel, maxLevel, mu, muInv, *projectionOperator, bcVelocity );


   P2ElementwiseBlendingVectorMassOperator vecMassOperator(storage, minLevel, maxLevel);

   vecMassOperator.apply(fStrong.uvw(), f.uvw(), maxLevel, All);

   auto stokesSolverTest =
       solvertemplates::temporary::stokesGMGFSSolver< P2P1StokesFullIcosahedralShellMapOperatorFS, P2ProjectNormalOperator >(
           storage, minLevel, maxLevel, stokesOperatorFS, projectionOperator, bcVelocity );

   projectionOperator->project(f, maxLevel, FreeslipBoundary);
   stokesSolverTest->solve( *stokesOperatorFS, u, f, maxLevel );

   stokesOperatorFS->apply( u, res, maxLevel, Inner | NeumannBoundary | FreeslipBoundary );
   res.assign( { 1.0, -1.0 }, { res, f }, maxLevel );

   real_t unscaledFinalResiduum = std::sqrt( res.dotGlobal( res, maxLevel, Inner | FreeslipBoundary ) );

   real_t unscaledResidualEpsilon = 1e-3;

   WALBERLA_LOG_INFO_ON_ROOT( "Final residual: " << unscaledFinalResiduum );
   WALBERLA_CHECK_LESS( unscaledFinalResiduum, unscaledResidualEpsilon );

   // AdiosWriter adiosWriter( "../../output", "solverTest", storage );
   // adiosWriter.add( u );
   // adiosWriter.add( f );
   // adiosWriter.write( maxLevel, 0U );

   return EXIT_SUCCESS;
}
