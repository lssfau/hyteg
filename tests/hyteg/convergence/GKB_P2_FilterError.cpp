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
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GKBSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

using walberla::math::pi;

void setRightBFSBoundaryNeumannPoiseuille( SetupPrimitiveStorage& setupStorage, const real_t & channelLength )
{
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   const real_t eps = 0.001;

   for ( const auto& it : setupStorage.getVertices() )
   {
      if ( std::fabs( it.second->getCoordinates()[0] - 1.0 ) < eps && it.second->getCoordinates()[1] > -1.0 + eps &&
           it.second->getCoordinates()[1] < 1.0 - eps )
      {
         setupStorage.setMeshBoundaryFlag( it.first, 2 );
      }
   }

   for ( const auto& it : setupStorage.getEdges() )
   {
      const auto edgeCoordinates = it.second->getCoordinates();
      if ( std::fabs( edgeCoordinates[0][0] - channelLength/2 ) < eps && std::fabs( edgeCoordinates[1][0] - channelLength/2 ) < eps )
      {
         setupStorage.setMeshBoundaryFlag( it.first, 2 );
      }
   }
}



void runBenchmark(const uint_t & minlevel, const uint_t & maxlevel,  const uint_t & channelLength)
{
  
   WALBERLA_LOG_INFO_ON_ROOT( "Poiseuille flow benchmark with channel length " << channelLength);


   /////////////////////////////////////////////////////////////////////////// Domain setup /////////////////////////////////////////////////////////
   //create a Rectangle as mesh with 4 triangles
   real_t halfLength = static_cast<double>(channelLength)/2;
   auto meshInfo = MeshInfo::meshRectangle( Point2D( {-halfLength, -1} ), Point2D( {halfLength, 1} ), MeshInfo::CRISSCROSS, channelLength, 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setRightBFSBoundaryNeumannPoiseuille( setupStorage, channelLength );
   
   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   auto                                      storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   writeDomainPartitioningVTK( storage, "../../output", "P2P1ChannelTest" );

   ////////////////////////////////////////////////////////////////////////// Problem and Functions //////////////////////////////////////////////////
   P2P1TaylorHoodFunction< real_t > r( "r", storage, minlevel,maxlevel );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, minlevel, maxlevel );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, minlevel, maxlevel );
   P2P1TaylorHoodFunction< real_t > Au( "Au", storage, minlevel, maxlevel );
   P2P1TaylorHoodFunction< real_t > u_exact( "u_exact", storage, minlevel, maxlevel );
   P2P1TaylorHoodFunction< real_t > err( "err", storage, minlevel, maxlevel );
   hyteg::P2P1TaylorHoodStokesOperator A( storage, minlevel, maxlevel );


      const auto setUVelocityBC = [halfLength]( const Point3D& x ) -> real_t {
         if ( x[0] < -halfLength + 1e-8 )
         {
             return real_c( (1 - x[1] * x[1])); //*real_c(sin(((x[1]+1)/2)*pi)) );
         }
         else
         {
            return real_c( 0 );
         }
      };

     
      const auto solutionU = []( const Point3D& x ) -> real_t { return real_c( 1 - x[1] * x[1] );};//*real_c(sin((fabs(x[1])/2 + 1/2)*pi)); };

      const auto solutionP = [channelLength]( const Point3D& x ) -> real_t {  return real_c( -2 * x[0] + channelLength); }; // normalize x

      u_exact.uvw[0].interpolate( solutionU, maxlevel );
      u_exact.p.interpolate( solutionP, maxlevel );

   u.uvw[0].interpolate( setUVelocityBC, maxlevel, DirichletBoundary );
     

    real_t localDoFs1 =static_cast<real_t>(hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxlevel ));
   


   /////////////////////////////////////////// SOLVER SETUP GKB FAILS ///////////////////////////////////////////////////////////
   GKBSolver_P2P1THOP_NO_AL GKB_HOUSE_solver( 
      storage, 
      maxlevel,
      CGSolver<P2ConstantVectorLaplaceOperator>(storage, maxlevel, maxlevel), 
      std::numeric_limits< PetscInt >::max(), 
      1e-6
   );



   GKB_HOUSE_solver.solve( A, u, f, maxlevel );
   

 
}


int main( int argc, char* argv[] ) {
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

  PETScManager petscManager( &argc, &argv );
  //const uint_t & minlevel, const uint_t & maxlevel, const uint_t & solverType, const uint_t & channelLength
   runBenchmark(2, 2, 2);

   return 0;
}
