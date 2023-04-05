#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/preconditioners/JacobiPreconditioner.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   uint_t numProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "./GMRESparam.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );
   parameters.listParameters();

   const uint_t minLevel       = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel       = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t maxKrylowDim   = parameters.getParameter< uint_t >( "maxKrylowDim" );
   const uint_t restartLength  = parameters.getParameter< uint_t >( "restartLength" );
   const real_t arnoldiTOL     = parameters.getParameter< real_t >( "arnoldiTOL" );
   const real_t approxTOL      = parameters.getParameter< real_t >( "approxTOL" );
   const real_t doubleOrthoTOL = parameters.getParameter< real_t >( "doubleOrthoTOL" );

   hyteg::MeshInfo meshInfo = hyteg::MeshInfo::meshRectangle(
       hyteg::Point2D(  0.0, 0.0  ), hyteg::Point2D(  1.0, 1.0  ), hyteg::MeshInfo::CRISS, 20, 20 );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, numProcesses );
   hyteg::loadbalancing::roundRobin( setupStorage, numProcesses );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P1Function< real_t > residual( "residual", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > f( "f", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > laplaceTimesFunction( "laplaceTimesFunction", storage, minLevel, maxLevel );

   std::function< real_t( const hyteg::Point3D& ) > boundaryConditions = []( const hyteg::Point3D& x ) {
      if ( x[0] <= 0.1 || x[1] <= 0.1 || x[0] >= 0.9 || x[1] >= 0.9 )
      {
         return ( 0.25 - std::pow( x[0] - 0.5, 2 ) );
      }
      WALBERLA_ABORT( "point is not on the boundary" );
   };

   u.interpolate( boundaryConditions, maxLevel, hyteg::DirichletBoundary );

   hyteg::P1ConstantLaplaceOperator L( storage, minLevel, maxLevel );
   L.computeInverseDiagonalOperatorValues();

   auto preCondi =
       std::make_shared< hyteg::JacobiPreconditioner< hyteg::P1ConstantLaplaceOperator > >( storage, minLevel, maxLevel, 2 );
   hyteg::GMRESSolver gmresSolver = hyteg::GMRESSolver< hyteg::P1ConstantLaplaceOperator >(
       storage, minLevel, maxLevel, maxKrylowDim, restartLength, arnoldiTOL, approxTOL, doubleOrthoTOL, preCondi );
   gmresSolver.solve( L, u, f, maxLevel );

   L.apply( u, laplaceTimesFunction, maxLevel, hyteg::Inner );
   residual.assign( { 1.0, -1.0 }, { f, laplaceTimesFunction }, maxLevel, hyteg::Inner );
   real_t residualEuclideanNorm = std::sqrt( residual.dotGlobal( residual, maxLevel, hyteg::Inner ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Euclidean norm of residual: " << residualEuclideanNorm );

   if ( parameters.getParameter< bool >( "vtkOutput" ) )
   {
      hyteg::VTKOutput vtkOutput( ".", "GMRESApp", storage );
      vtkOutput.add( u );
      vtkOutput.add( residual );
      vtkOutput.add( f );
      vtkOutput.write( maxLevel );
   }
}
