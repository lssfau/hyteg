#include <core/Environment.h>
#include <core/config/Config.h>

#include "core/timing/Timer.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/geometry/AffineMap.hpp"
#include "tinyhhg_core/geometry/CircularMap.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/p1functionspace/P1BlendingOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1BlendingOperatorNew.hpp"
#include "tinyhhg_core/p1functionspace/P1PolynomialBlendingOperatorNew.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hhg;

int main(int argc, char* argv[])
{

   /// create enviroment
   walberla::Environment walberlaEnv(argc, argv);
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   const size_t level = 4;
   const uint_t maxiter = 10000;

   /// read mesh file and create storage
   MeshInfo meshInfo = MeshInfo::fromGmshFile( "../data/meshes/unitsquare_with_circular_hole.msh" );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   Point3D circleCenter{{0.5, 0.5, 0}};
   real_t circleRadius = 0.25;

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
   std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

//   typedef hhg::P1BlendingLaplaceOperator SolveOperator;

   typedef hhg::P1BlendingLaplaceOperatorNew SolveOperator;
   SolveOperator L(storage, level, level);

//   typedef hhg::P1PolynomialBlendingLaplaceOperatorNew SolveOperator;
//   const uint_t interpolationLevel = level;
//   const uint_t polyDegree = 12;
//   SolveOperator L(storage, level, level, interpolationLevel);
//   L.interpolateStencils(polyDegree);
//   L.useDegree(polyDegree);

   P1BlendingMassOperatorNew M(storage, level, level);

   auto x = std::make_shared<hhg::P1Function< real_t > >("x", storage, level, level);
   auto y = std::make_shared<hhg::P1Function< real_t > >("y", storage, level, level);
   auto u = std::make_shared<hhg::P1Function< real_t > >("u", storage, level, level);
   auto u_exact = std::make_shared<hhg::P1Function< real_t > >("u_exact", storage, level, level);
   auto f = std::make_shared<hhg::P1Function< real_t > >("f", storage, level, level);
   auto r = std::make_shared<hhg::P1Function< real_t > >("r", storage, level, level);
   auto err = std::make_shared<hhg::P1Function< real_t > >("err", storage, level, level);
   auto helper = std::make_shared<hhg::P1Function< real_t > >("helper", storage, level, level);

   std::function<real_t(const hhg::Point3D&)> ones  = [](const hhg::Point3D&) { return 1.0; };
   std::function<real_t(const hhg::Point3D&)> zeros  = [](const hhg::Point3D&) { return 0.0; };

   helper->interpolate(ones, level);
   real_t npoints = helper->dot(*helper, level);

   std::function<real_t(const hhg::Point3D&)> tmp_x = [&](const hhg::Point3D& x_) {
     return x_[0];
   };

   std::function<real_t(const hhg::Point3D&)> tmp_y = [&](const hhg::Point3D& x_) {
     return x_[1];
   };

   std::function<real_t(const hhg::Point3D&)> exact = [](const hhg::Point3D& x) { return sin(x[0])*cos(x[1])/(x[0]*x[1] + 1); };
   std::function<real_t(const hhg::Point3D&)> rhs = [](const hhg::Point3D& x) { return 2*(-(pow(x[0], 2) + pow(x[1], 2))*sin(x[0])*cos(x[1]) + pow(x[0]*x[1] + 1, 2)*sin(x[0])*cos(x[1]) + (x[0]*x[1] + 1)*(-x[0]*sin(x[0])*sin(x[1]) + x[1]*cos(x[0])*cos(x[1])))/pow(x[0]*x[1] + 1, 3); };

   x->interpolate(tmp_x, level, hhg::All);
   y->interpolate(tmp_y, level, hhg::All);

   u->interpolate(exact, level, hhg::DirichletBoundary);
   u_exact->interpolate(exact, level);
   helper->interpolate(rhs, level, hhg::All);
   M.apply(*helper, *f, level, hhg::All);

   auto solver = hhg::CGSolver<hhg::P1Function< real_t >, SolveOperator>(storage, level, level);
   solver.solve(L, *u, *f, *r, level, 1e-10, maxiter, hhg::Inner, true);

   err->assign({1.0, -1.0}, {u.get(), u_exact.get()}, level, hhg::All);


   real_t discr_l2_err = std::sqrt(err->dot(*err, level) / npoints);
   WALBERLA_LOG_INFO_ON_ROOT("discrete L2 error = " << discr_l2_err);

   VTKOutput vtkOutput("../output", "cg_P1_blending");
   vtkOutput.add(x.get());
   vtkOutput.add(y.get());
   vtkOutput.add(u.get());
   vtkOutput.add(u_exact.get());
   vtkOutput.add(err.get());
   vtkOutput.write(level);

   return 0;
}
