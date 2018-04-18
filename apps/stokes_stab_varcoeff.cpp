#include "core/DataTypes.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/p1functionspace/P1BlendingOperatorNew.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/composites/P1BlendingStokesOperator.hpp"
#include "tinyhhg_core/composites/P1PolynomialBlendingStokesOperator.hpp"
#include "tinyhhg_core/composites/P1CoefficientStokesOperator.hpp"
#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/geometry/CircularMap.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
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

   std::string meshFileName = "../data/meshes/unitsquare_with_circular_hole.msh";

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
   size_t maxLevel = 5;
   size_t maxiter  = 5000;
   const uint_t coarseMaxiter = 200;
   const real_t mg_tolerance = 1e-9;
   const uint_t maxOuterIter = 30;
   const uint_t interpolationLevel = 4;
   const uint_t polyDegree = 12;

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   hhg::P1Function< real_t >                    one( "one", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              tmp( "tmp", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              tmp2( "tmp", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              r( "r", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              f( "f", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              u( "u", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              u_exact( "u_exact", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              err( "err", storage, minLevel, maxLevel );
   hhg::P1StokesFunction< real_t >              Lu( "Lu", storage, minLevel, maxLevel );

   typedef hhg::P1BlendingStokesOperator SolveOperator;
   SolveOperator L( storage, minLevel, maxLevel );

//   typedef hhg::P1PolynomialBlendingStokesOperator SolveOperator;
//   SolveOperator L( storage, minLevel, maxLevel, interpolationLevel );
//   L.interpolateStencils(polyDegree);
//   L.useDegree(polyDegree);

   P1BlendingMassOperatorNew M(storage, minLevel, maxLevel);

   std::function< real_t( const hhg::Point3D& ) > zeros   = []( const hhg::Point3D& ) { return 0.0; };
   std::function< real_t( const hhg::Point3D& ) > ones   = []( const hhg::Point3D& ) { return 1.0; };
   std::function<real_t(const hhg::Point3D&)> exact_u = [](const hhg::Point3D& x) { return sin(x[0])*cos(x[1])/(x[0] + 1); };
   std::function<real_t(const hhg::Point3D&)> exact_v = [](const hhg::Point3D& x) { return -((x[0] + 1)*cos(x[0]) - sin(x[0]))*sin(x[1])/pow(x[0] + 1, 2); };
   std::function<real_t(const hhg::Point3D&)> exact_p = [](const hhg::Point3D& x) { return pow(x[0], 2)*pow(x[1], 3); };
   std::function<real_t(const hhg::Point3D&)> rhs_u = [](const hhg::Point3D& x) { return (2*x[0]*pow(x[1], 3)*pow(x[0] + 1, 3) + 2.0*pow(x[0] + 1, 2)*sin(x[0])*cos(x[1]) + (4.0*x[0] + 4.0)*cos(x[0])*cos(x[1]) - 2.0*((x[0] + 1)*cos(x[0]) + sin(x[0]))*cos(x[1]))/pow(x[0] + 1, 3); };
   std::function<real_t(const hhg::Point3D&)> rhs_v = [](const hhg::Point3D& x) { return (3*pow(x[0], 2)*pow(x[1], 2)*pow(x[0] + 1, 4) + 2.0*pow(x[0] + 1, 2)*(-(x[0] + 1)*cos(x[0]) + 2*sin(x[0]))*sin(x[1]) + 6.0*((x[0] + 1)*cos(x[0]) - sin(x[0]))*sin(x[1]))/pow(x[0] + 1, 4); };

   u.u.interpolate( exact_u, maxLevel, hhg::DirichletBoundary );
   u.v.interpolate( exact_v, maxLevel, hhg::DirichletBoundary );

   u_exact.u.interpolate( exact_u, maxLevel, hhg::All );
   u_exact.v.interpolate( exact_v, maxLevel, hhg::All );
   u_exact.p.interpolate( exact_p, maxLevel, hhg::All );

   hhg::vertexdof::projectMean(u_exact.p, tmp.p, maxLevel);

   // Integrate RHS for u
   tmp.u.interpolate(rhs_u, maxLevel, hhg::All);
   M.apply(tmp.u, f.u, maxLevel, hhg::All);

   // Integrate RHS for v
   tmp.v.interpolate(rhs_v, maxLevel, hhg::All);
   M.apply(tmp.v, f.v, maxLevel, hhg::All);

   // Apply compatibility projection
   one.interpolate(zeros, maxLevel, hhg::All);
   one.interpolate(ones, maxLevel, hhg::DirichletBoundary);
   L.div_x.apply(u_exact.u, tmp2.p, maxLevel, hhg::DirichletBoundary, Replace);
   L.div_y.apply(u_exact.v, tmp2.p, maxLevel, hhg::DirichletBoundary, Add);
   real_t corr = one.dot(tmp2.p, maxLevel, hhg::DirichletBoundary);
   M.apply(one, tmp2.p, maxLevel, hhg::DirichletBoundary, Replace);
   real_t volume = one.dot(tmp2.p, maxLevel, hhg::DirichletBoundary);
   tmp2.p.assign({corr/volume}, {&one}, maxLevel, hhg::DirichletBoundary);
   M.apply(tmp2.p, f.p, maxLevel, hhg::DirichletBoundary);

   one.interpolate(ones, maxLevel, hhg::All);
   real_t npoints = one.dot( one, maxLevel );

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
      hhg::vertexdof::projectMean(u.p, tmp.p, maxLevel);


      L.apply(u, r, maxLevel, hhg::Inner | hhg::NeumannBoundary);

      r.assign({1.0, -1.0}, { &f, &r }, maxLevel, hhg::Inner | hhg::NeumannBoundary);
      real_t abs_res = std::sqrt(r.dot(r, maxLevel, hhg::Inner | hhg::NeumannBoundary));
      rel_res = abs_res / begin_res;
      WALBERLA_LOG_INFO_ON_ROOT(hhg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e",outer+1,abs_res, rel_res, abs_res/abs_res_old, end-start));
      totalTime += end-start;

      if (abs_res/abs_res_old > 0.95) {
         break;
      }

      if (outer >= convergenceStartIter) {
         averageConvergenceRate += abs_res/abs_res_old;
      }

      abs_res_old = abs_res;

      if (rel_res < mg_tolerance)
      {
         break;
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << totalTime);
   WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(outer+1-convergenceStartIter));
   WALBERLA_LOG_INFO_ON_ROOT("Dofs: " << 3 * npoints);

   err.assign( {1.0, -1.0}, {&u, &u_exact}, maxLevel );
   real_t discr_u_l2_err = std::sqrt( (err.u.dot( err.u, maxLevel ) + err.v.dot( err.v, maxLevel )) / (2*npoints) );
   real_t discr_p_l2_err = std::sqrt( (err.p.dot( err.p, maxLevel )) / (npoints) );

   WALBERLA_LOG_INFO_ON_ROOT("velocity_err = " << std::scientific << discr_u_l2_err);
   WALBERLA_LOG_INFO_ON_ROOT("pressure_err = " << std::scientific << discr_p_l2_err);

   // u_u*iHat + u_v*jHat
   hhg::VTKOutput vtkOutput( "../output", "stokes_stab_varcoeff" );
   vtkOutput.add( &u.u );
   vtkOutput.add( &u.v );
   vtkOutput.add( &u.p );
   vtkOutput.add( &u_exact.u );
   vtkOutput.add( &u_exact.v );
   vtkOutput.add( &u_exact.p );
   vtkOutput.add( &err.u );
   vtkOutput.add( &err.v );
   vtkOutput.add( &err.p );
   vtkOutput.add( &r.u );
   vtkOutput.add( &r.v );
   vtkOutput.add( &r.p );
   vtkOutput.add( &f.u );
   vtkOutput.add( &f.v );
   vtkOutput.add( &f.p );
   vtkOutput.write( maxLevel, 0 );
   return EXIT_SUCCESS;
}
