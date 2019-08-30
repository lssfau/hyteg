#include "core/DataTypes.h"
#include <core/config/Create.h>
#include <hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp>

#include "hyteg/VTKWriter.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/composites/P1StokesOperator.hpp"
#include "hyteg/composites/P1PolynomialBlendingStokesOperator.hpp"
#include "hyteg/composites/P1CoefficientStokesOperator.hpp"
#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/Format.hpp"

using walberla::real_t;
using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr<walberla::config::Config> cfg;

   if (argc == 1) {
      walberla::shared_ptr<walberla::config::Config> cfg_(new walberla::config::Config);
      cfg_->readParameterFile("../../data/param/stokes_blending.prm");
      cfg = cfg_;
   } else {
      cfg = walberla::config::create(argc, argv);
   }
   WALBERLA_LOG_INFO_ON_ROOT("config = " << *cfg);
   walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

   const std::string meshFileName = parameters.getParameter<std::string>("meshFileName");
   const uint_t minLevel = parameters.getParameter<uint_t>("minLevel");
   const uint_t maxLevel = parameters.getParameter<uint_t>("maxLevel");
   const uint_t polyDegree = parameters.getParameter<uint_t>("polyDegree");
   const uint_t interpolationLevel = parameters.getParameter<uint_t>("interpolationLevel");
   const uint_t coarseMaxiter =  parameters.getParameter<uint_t>("coarseMaxiter");
   const uint_t maxOuterIter = parameters.getParameter<uint_t>("maxOuterIter");

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{{parameters.getParameter<real_t>("circleCenterX"), parameters.getParameter<real_t>("circleCenterY"), 0}};
   real_t  circleRadius = parameters.getParameter<real_t>("circleRadius");

   for( const auto & it : setupStorage.getFaces() )
   {
      Face& face = *(it.second);

      std::vector< PrimitiveID > neighborEdgesOnBoundary = face.neighborEdges();
      std::remove_if( neighborEdgesOnBoundary.begin(), neighborEdgesOnBoundary.end(),
                      [ &setupStorage ]( const PrimitiveID & id ){ return !setupStorage.onBoundary( id ); } );

      if( neighborEdgesOnBoundary.size() > 0 )
      {
         Edge& edge = *setupStorage.getEdge( neighborEdgesOnBoundary[0] );

         if( ( edge.getCoordinates()[0] - circleCenter ).norm() < 0.4 )
         {
            setupStorage.setGeometryMap( edge.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
            setupStorage.setGeometryMap( face.getID(),
                                         std::make_shared< CircularMap >( face, setupStorage, circleCenter, circleRadius ) );
         }
      }
   }

   hyteg::loadbalancing::roundRobin( setupStorage );


   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   hyteg::P1Function< real_t >                    one( "one", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t >              tmp( "tmp", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t >              tmp2( "tmp", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t >              r( "r", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t >              f( "f", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t >              u( "u", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t >              u_exact( "u_exact", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t >              err( "err", storage, minLevel, maxLevel );
   hyteg::P1StokesFunction< real_t >              Lu( "Lu", storage, minLevel, maxLevel );

   //typedef hyteg::P1PolynomialBlendingStokesOperator SolveOperator;
   auto start = walberla::timing::getWcTime();
   hyteg::P1PolynomialBlendingStokesOperator L( storage, minLevel, maxLevel, interpolationLevel );
   L.interpolateStencils(polyDegree);
   L.useDegree(polyDegree);
   auto end = walberla::timing::getWcTime();
   real_t setupTime = end-start;

   P1BlendingMassOperator M(storage, minLevel, maxLevel);

   std::function< real_t( const hyteg::Point3D& ) > zeros   = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > ones   = []( const hyteg::Point3D& ) { return 1.0; };
   std::function<real_t(const hyteg::Point3D&)> exact_u = [](const hyteg::Point3D& x) { return sin(x[0])*cos(x[1])/(x[0] + 1); };
   std::function<real_t(const hyteg::Point3D&)> exact_v = [](const hyteg::Point3D& x) { return -((x[0] + 1)*cos(x[0]) - sin(x[0]))*sin(x[1])/pow(x[0] + 1, 2); };
   std::function<real_t(const hyteg::Point3D&)> exact_p = [](const hyteg::Point3D& x) { return pow(x[0], 2)*pow(x[1], 3); };
   std::function<real_t(const hyteg::Point3D&)> rhs_u = [](const hyteg::Point3D& x) { return (2*x[0]*pow(x[1], 3)*pow(x[0] + 1, 3) + 2.0*pow(x[0] + 1, 2)*sin(x[0])*cos(x[1]) + (4.0*x[0] + 4.0)*cos(x[0])*cos(x[1]) - 2.0*((x[0] + 1)*cos(x[0]) + sin(x[0]))*cos(x[1]))/pow(x[0] + 1, 3); };
   std::function<real_t(const hyteg::Point3D&)> rhs_v = [](const hyteg::Point3D& x) { return (3*pow(x[0], 2)*pow(x[1], 2)*pow(x[0] + 1, 4) + 2.0*pow(x[0] + 1, 2)*(-(x[0] + 1)*cos(x[0]) + 2*sin(x[0]))*sin(x[1]) + 6.0*((x[0] + 1)*cos(x[0]) - sin(x[0]))*sin(x[1]))/pow(x[0] + 1, 4); };

   u.u.interpolate( exact_u, maxLevel, hyteg::DirichletBoundary );
   u.v.interpolate( exact_v, maxLevel, hyteg::DirichletBoundary );

   u_exact.u.interpolate( exact_u, maxLevel, hyteg::All );
   u_exact.v.interpolate( exact_v, maxLevel, hyteg::All );
   u_exact.p.interpolate( exact_p, maxLevel, hyteg::All );

   hyteg::vertexdof::projectMean(u_exact.p, maxLevel);

   // Integrate RHS for u
   tmp.u.interpolate(rhs_u, maxLevel, hyteg::All);
   M.apply(tmp.u, f.u, maxLevel, hyteg::All);

   // Integrate RHS for v
   tmp.v.interpolate(rhs_v, maxLevel, hyteg::All);
   M.apply(tmp.v, f.v, maxLevel, hyteg::All);

   // Apply compatibility projection
   one.interpolate(zeros, maxLevel, hyteg::All);
   one.interpolate(ones, maxLevel, hyteg::DirichletBoundary);
   L.div_x.apply(u_exact.u, tmp2.p, maxLevel, hyteg::DirichletBoundary, Replace);
   L.div_y.apply(u_exact.v, tmp2.p, maxLevel, hyteg::DirichletBoundary, Add);
   real_t corr = one.dotGlobal(tmp2.p, maxLevel, hyteg::DirichletBoundary);
   M.apply(one, tmp2.p, maxLevel, hyteg::DirichletBoundary, Replace);
   real_t volume = one.dotGlobal(tmp2.p, maxLevel, hyteg::DirichletBoundary);
   tmp2.p.assign({corr/volume}, {one}, maxLevel, hyteg::DirichletBoundary);
   M.apply(tmp2.p, f.p, maxLevel, hyteg::DirichletBoundary);

   one.interpolate(ones, maxLevel, hyteg::All);
   real_t npoints = one.dotGlobal( one, maxLevel );

   auto smoother = std::make_shared< hyteg::UzawaSmoother< hyteg::P1PolynomialBlendingStokesOperator > >(
      storage, minLevel, maxLevel, 0.3 );
   auto coarseGridSolver = std::make_shared< hyteg::MinResSolver< hyteg::P1PolynomialBlendingStokesOperator > >( storage, minLevel, minLevel, coarseMaxiter );
   auto restrictionOperator = std::make_shared< hyteg::P1P1StokesToP1P1StokesRestriction>();
   auto prolongationOperator = std::make_shared< hyteg::P1P1StokesToP1P1StokesProlongation >();

   auto solver = hyteg::GeometricMultigridSolver< hyteg::P1PolynomialBlendingStokesOperator >(
      storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2, 2 );

//   typedef hyteg::P1P1StokesToP1P1StokesRestriction  RestrictionOperator;
//   typedef hyteg::P1P1StokesToP1P1StokesProlongation ProlongationOperator;
//   typedef hyteg::MinResSolver< hyteg::P1StokesFunction<real_t>, SolveOperator > MinResSolver_T;
//
//   MinResSolver_T minResSolver( storage, minLevel, maxLevel );
//   RestrictionOperator restrictionOperator;
//   ProlongationOperator prolongationOperator;
//
//   typedef hyteg::UzawaSolver<hyteg::P1StokesFunction<real_t>, SolveOperator, MinResSolver_T, RestrictionOperator, ProlongationOperator, true> Solver;
//   auto solver = Solver(storage, minResSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2, 2);

   WALBERLA_LOG_INFO_ON_ROOT("Starting Uzawa cycles");
   WALBERLA_LOG_INFO_ON_ROOT(hyteg::format("%6s|%10s|%10s|%10s|%10s","iter","abs_res","rel_res","conv","Time"));

   L.apply(u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary);
   r.assign({1.0, -1.0}, { f, r }, maxLevel, hyteg::Inner | hyteg::NeumannBoundary);
   real_t begin_res = std::sqrt(r.dotGlobal(r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary));
   real_t abs_res_old = begin_res;
   real_t rel_res = 1.0;

   WALBERLA_LOG_INFO_ON_ROOT(hyteg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e",0,begin_res, rel_res, begin_res/abs_res_old, 0));

   real_t solveTime = real_c(0.0);
   real_t averageConvergenceRate = real_c(0.0);
   const uint_t convergenceStartIter = 3;

   uint_t outer;
   for (outer = 0; outer < maxOuterIter; ++outer) {
      start = walberla::timing::getWcTime();
      solver.solve( L, u, f, maxLevel );
      end = walberla::timing::getWcTime();
      hyteg::vertexdof::projectMean(u.p, maxLevel);


      L.apply(u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary);

      r.assign({1.0, -1.0}, { f, r }, maxLevel, hyteg::Inner | hyteg::NeumannBoundary);
      real_t abs_res = std::sqrt(r.dotGlobal(r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary));
      rel_res = abs_res / begin_res;
      WALBERLA_LOG_INFO_ON_ROOT(hyteg::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e",outer+1,abs_res, rel_res, abs_res/abs_res_old, end-start));
      solveTime += end-start;

      if (abs_res/abs_res_old > 0.95) {
         break;
      }

      if (outer >= convergenceStartIter) {
         averageConvergenceRate += abs_res/abs_res_old;
      }

      abs_res_old = abs_res;

//      if (rel_res < mg_tolerance)
//      {
//         break;
//      }
   }

//   WALBERLA_CHECK_LESS( outer, maxOuterIter );

   WALBERLA_LOG_INFO_ON_ROOT("Setup time: " << setupTime);
   WALBERLA_LOG_INFO_ON_ROOT("Solve time: " << solveTime);
   WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << setupTime + solveTime);
   WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(outer+1-convergenceStartIter));
   WALBERLA_LOG_INFO_ON_ROOT("Dofs: " << 3 * npoints);

   err.assign( {1.0, -1.0}, {u, u_exact}, maxLevel );
   real_t discr_u_l2_err = std::sqrt( (err.u.dotGlobal( err.u, maxLevel ) + err.v.dotGlobal( err.v, maxLevel )) / (2*npoints) );
   real_t discr_p_l2_err = std::sqrt( (err.p.dotGlobal( err.p, maxLevel )) / (npoints) );

   WALBERLA_LOG_INFO_ON_ROOT("velocity_err = " << std::scientific << discr_u_l2_err);
   WALBERLA_LOG_INFO_ON_ROOT("pressure_err = " << std::scientific << discr_p_l2_err);

   WALBERLA_CHECK_LESS( discr_u_l2_err, 2e-4 );
   WALBERLA_CHECK_LESS( discr_p_l2_err, 2.1e-2 );

   // u_u*iHat + u_v*jHat
//   hyteg::VTKOutput vtkOutput( "../output", "stokes_stab_varcoeff" );
//   vtkOutput.add( u.u );
//   vtkOutput.add( u.v );
//   vtkOutput.add( u.p );
//   vtkOutput.add( u_exact.u );
//   vtkOutput.add( u_exact.v );
//   vtkOutput.add( u_exact.p );
//   vtkOutput.add( err.u );
//   vtkOutput.add( err.v );
//   vtkOutput.add( err.p );
//   vtkOutput.add( r.u );
//   vtkOutput.add( r.v );
//   vtkOutput.add( r.p );
//   vtkOutput.add( f.u );
//   vtkOutput.add( f.v );
//   vtkOutput.add( f.p );
//   vtkOutput.write( maxLevel, 0 );
   return EXIT_SUCCESS;
}
