#include "core/DataTypes.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/composites/P1BlendingStokesOperator.hpp"
#include "tinyhhg_core/composites/P1PolynomialBlendingStokesOperator.hpp"
#include "tinyhhg_core/composites/P1CoefficientStokesOperator.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/geometry/CircularMap.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/UzawaSolver.hpp"
#include "tinyhhg_core/format.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName = "../data/meshes/unitsquare_with_circular_hole_neumann.msh";

   hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{{0.5, 0.5, 0}};
   real_t  circleRadius = 0.25;

   for( auto it = setupStorage.beginFaces(); it != setupStorage.endFaces(); ++it )
   {
      Face& face = *it->second;

      if( face.hasBoundaryEdge() )
      {
         Edge& edge = *setupStorage.getEdge( face.edgesOnBoundary[0] );

         if( ( edge.getCoordinates()[0] - circleCenter ).norm() < 0.4 )
         {
            setupStorage.setGeometryMap( edge.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
            setupStorage.setGeometryMap( face.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
         }
      }
   }

   hhg::loadbalancing::roundRobin( setupStorage );

   size_t minLevel = 2;
   size_t maxLevel = 4;
   size_t maxiter  = 5000;
   const uint_t coarseMaxiter = 200;
   const real_t mg_tolerance = 1e-9;
   const uint_t maxOuterIter = 100;
   const uint_t interpolationLevel = 3;
   const uint_t polyDegree = 4;

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::P1Function< real_t >                    tmp( "tmp", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              r( "r", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              f( "f", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              u( "u", storage, minLevel, maxLevel );
   std::shared_ptr< hhg::P1Function< real_t > > coefficient =
       std::make_shared< hhg::P1Function< real_t > >( "coeff", storage, minLevel, maxLevel );

   typedef hhg::P1BlendingStokesOperator SolveOperator;
   SolveOperator L( storage, minLevel, maxLevel );

//   typedef hhg::P1PolynomialBlendingStokesOperator SolveOperator;
//   SolveOperator L( storage, minLevel, maxLevel, interpolationLevel );
//   L.interpolateStencils(polyDegree);
//   L.useDegree(polyDegree);

   std::function< real_t( const hhg::Point3D& ) > coeff  = []( const hhg::Point3D& x ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > zero   = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones   = []( const hhg::Point3D& ) { return 1.0; };
   std::function< real_t( const hhg::Point3D& ) > inflow = []( const hhg::Point3D& x ) {
      if( x[0] < 1e-4 )
      {
         return 4.0 * x[1] * ( 1.0 - x[1] );
      }

      return 0.0;
   };

   u.u.interpolate( inflow, maxLevel, hhg::DirichletBoundary );
   u.v.interpolate( zero, maxLevel, hhg::DirichletBoundary );
   coefficient->interpolate( coeff, maxLevel );

//   auto solver = hhg::MinResSolver< hhg::P1StokesFunction< real_t >, SolveOperator >( storage, minLevel, maxLevel );
//   solver.solve( L, u, f, r, maxLevel, 1e-8, maxiter, hhg::Inner | hhg::NeumannBoundary, true );

   typedef hhg::UzawaSolver<hhg::P1StokesFunction<real_t>, SolveOperator> Solver;
   auto solver = Solver(storage, minLevel, maxLevel);

   WALBERLA_LOG_INFO_ON_ROOT("Starting Uzawa cycles");
   WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6s|%10s|%10s|%10s|%10s","iter","abs_res","rel_res","conv","Time"));

   L.apply(u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary);
   r.assign({1.0, -1.0}, { &f, &r }, maxLevel, hhg::Inner | hhg::NeumannBoundary);
   real_t begin_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner | hhg::NeumannBoundary));
   real_t abs_res_old = begin_res;
   real_t rel_res = 1.0;

   WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e",0,begin_res, rel_res, begin_res/abs_res_old, 0));

   real_t totalTime = real_c(0.0);
   real_t averageConvergenceRate = real_c(0.0);
   const uint_t convergenceStartIter = 3;

   uint_t outer;
   for (outer = 0; outer < maxOuterIter; ++outer) {
      auto start = walberla::timing::getWcTime();
      solver.solve(L, u, f, r, maxLevel, 1e-6, coarseMaxiter, hhg::Inner | hhg::NeumannBoundary, Solver::CycleType::VCYCLE, true);
      auto end = walberla::timing::getWcTime();
//      hhg::vertexdof::projectMean(u->p, *tmp, maxLevel);


      L.apply(u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary);

      r.assign({1.0, -1.0}, { &f, &r }, maxLevel, hhg::Inner | hhg::NeumannBoundary);
      real_t abs_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner | hhg::NeumannBoundary));
      rel_res = abs_res / begin_res;
      WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e",outer+1,abs_res, rel_res, abs_res/abs_res_old, end-start));
      totalTime += end-start;

      if (outer >= convergenceStartIter) {
         averageConvergenceRate += abs_res/abs_res_old;
      }

      abs_res_old = abs_res;

      if (rel_res < mg_tolerance)
      {
         break;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << std::scientific << totalTime);
   WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(outer+1-convergenceStartIter));

   // u_u*iHat + u_v*jHat
   hhg::VTKOutput vtkOutput( "../output", "stokes_stab_varcoeff" );
   vtkOutput.add( &u.u );
   vtkOutput.add( &u.v );
   vtkOutput.add( &u.p );
   vtkOutput.write( maxLevel, 0 );
   return EXIT_SUCCESS;
}
