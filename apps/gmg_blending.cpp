/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#define FP_FAST_FMA
#define FP_FAST_FMAF
#define FP_FAST_FMAL

#include <core/timing/Timer.h>
#include <core/Environment.h>
#include <core/config/Create.h>

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p1functionspace/P1SurrogateOperator.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1QuadraticProlongation.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "core/Format.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

using namespace hyteg;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  walberla::shared_ptr<walberla::config::Config> cfg;

  if (argc == 1) {
    walberla::shared_ptr<walberla::config::Config> cfg_(new walberla::config::Config);
    cfg_->readParameterFile("../data/param/gmg_blending.prm");
    cfg = cfg_;
  } else {
    cfg = walberla::config::create(argc, argv);
  }
  WALBERLA_LOG_INFO_ON_ROOT("config = " << *cfg);
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  uint_t level_H = parameters.getParameter<uint_t>("level_h_coarse");
  uint_t level_h = parameters.getParameter<uint_t>("level_h_fine");
  const uint_t minLevel = level_H;
  const uint_t maxLevel = level_h;
  const uint_t maxMemoryLevel = maxLevel + 1;
  const uint_t maxPolyDegree = parameters.getParameter<uint_t>("maxPolyDegree");
  const uint_t interpolationLevel = parameters.getParameter<uint_t>("interpolationLevel");
  const uint_t max_outer_iter =  parameters.getParameter<uint_t>("max_outer_iter");
  const uint_t max_cg_iter =  parameters.getParameter<uint_t>("max_cg_iter");
  const real_t mg_tolerance = parameters.getParameter<real_t>("mg_tolerance");
  const real_t coarse_tolerance = parameters.getParameter<real_t>("coarse_tolerance");
  const bool polynomialOperator = parameters.getParameter<bool>("polynomialOperator");

  MeshInfo meshInfo = MeshInfo::fromGmshFile(parameters.getParameter<std::string>("meshFilename"));
  // MeshInfo meshInfo = MeshInfo::meshUnitSquare(level_H);
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  Point3D circleCenter{{0.5, 0.5, 0}};
  real_t circleRadius = 0.25;

  for( const auto & it : setupStorage.getFaces() )
  {
     Face& face = *(it.second);

     std::vector< PrimitiveID > neighborEdgesOnBoundary = face.neighborEdges();
     neighborEdgesOnBoundary.erase(
      std::remove_if( neighborEdgesOnBoundary.begin(), neighborEdgesOnBoundary.end(),
                     [ &setupStorage ]( const PrimitiveID & id ){ return !setupStorage.onBoundary( id ); } )
      , neighborEdgesOnBoundary.end());
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

  if (polynomialOperator) {
    WALBERLA_LOG_INFO_ON_ROOT("Polynomial Operator enabled");
  } else {
    WALBERLA_LOG_INFO_ON_ROOT("Polynomial Operator disabled");
  }

  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

#ifdef WALBERLA_BUILD_WITH_PARMETIS
//  loadbalancing::distributed::parmetis( *storage );
#endif

  hyteg::P1Function< real_t > r("r", storage, minLevel, maxMemoryLevel);
//  hyteg::P1Function< real_t > r_fe("r_fe", storage, minLevel, maxMemoryLevel);
  hyteg::P1Function< real_t > f("f", storage, minLevel, maxMemoryLevel);
  hyteg::P1Function< real_t > u("u", storage, minLevel, maxMemoryLevel);
//  hyteg::P1Function< real_t > u_fe("u_fe", storage, minLevel, maxMemoryLevel);
  hyteg::P1Function< real_t > Lu("Lu", storage, minLevel, maxMemoryLevel);
  hyteg::P1Function< real_t > u_exact("u_exact", storage, minLevel, maxMemoryLevel);
  hyteg::P1Function< real_t > err("err", storage, minLevel, maxMemoryLevel);
//  hyteg::P1Function< real_t > err_est("err_est", storage, minLevel, maxMemoryLevel);
  hyteg::P1Function< real_t > npoints_helper("npoints_helper", storage, minLevel, maxMemoryLevel);
  hyteg::P1Function< real_t > tmp("tmp", storage, minLevel, maxMemoryLevel);
//  hyteg::P1Function< real_t > tmp_fe("tmp_fe", storage, minLevel, maxMemoryLevel);
  auto coordX = std::make_shared<hyteg::P1Function<real_t>>("x", storage, minLevel, maxMemoryLevel);
  auto coordY = std::make_shared<hyteg::P1Function<real_t>>("y", storage, minLevel, maxMemoryLevel);

  typedef hyteg::P1BlendingMassOperator MassOperator;
  typedef hyteg::P1BlendingLaplaceOperator SolveOperatorNodal;
  typedef hyteg::P1SurrogateLaplaceOperator SolveOperatorPoly;

  std::function<real_t(const hyteg::Point3D&)> exact = [](const hyteg::Point3D& x) { return sin(x[0])*sinh(x[1]); };
  // std::function<real_t(const hyteg::Point3D&)> rhs = [](const hyteg::Point3D& x) { return -2*(x[0] + 1)*cos(x[0])*sinh(x[1]) - 3*sin(x[0])*cosh(x[1]); };
  std::function<real_t(const hyteg::Point3D&)> rhs = [](const hyteg::Point3D&) { return 0.0; };

  std::function<real_t(const hyteg::Point3D&)> zeros = [](const hyteg::Point3D&) { return 0.0; };
  std::function<real_t(const hyteg::Point3D&)> ones  = [](const hyteg::Point3D&) { return 1.0; };

  std::function<real_t(const hyteg::Point3D&)> xExpr = [](const hyteg::Point3D& x) { return x[0]; };
  std::function<real_t(const hyteg::Point3D&)> yExpr = [](const hyteg::Point3D& x) { return x[1]; };

  coordX->interpolate(xExpr, maxLevel, hyteg::All);
  coordY->interpolate(yExpr, maxLevel, hyteg::All);

  WALBERLA_LOG_INFO_ON_ROOT("Interpolating u");
  u.interpolate(exact, maxLevel, hyteg::DirichletBoundary);
  u.interpolate(exact, maxMemoryLevel, hyteg::DirichletBoundary);

//  u_fe.interpolate(exact, maxLevel, hyteg::DirichletBoundary);
//  u_fe.interpolate(exact, maxMemoryLevel, hyteg::DirichletBoundary);

  WALBERLA_LOG_INFO_ON_ROOT("Setting up operators");
  MassOperator M(storage, minLevel, maxMemoryLevel);
  std::shared_ptr<SolveOperatorPoly> Lpoly;
  std::shared_ptr<SolveOperatorNodal> L;
  uint_t useDegree = maxPolyDegree;

  auto start = walberla::timing::getWcTime();
  L = std::make_shared<SolveOperatorNodal>(storage, minLevel, maxMemoryLevel);

  if (polynomialOperator) {
    Lpoly = std::make_shared<SolveOperatorPoly>(storage, minLevel, maxLevel);

    Lpoly->interpolateStencils(useDegree, interpolationLevel);

//    real_t polyError34 = Lpoly->lInfinityError(3, 4, 5);
//    real_t polyError45 = Lpoly->lInfinityError(4, 5, 5);
//    real_t polyError56 = Lpoly->lInfinityError(5, 6, 5);
//    real_t polyError67 = Lpoly->lInfinityError(6, 7, 5);
//    real_t polyError78 = Lpoly->lInfinityError(7, 8, 5);
//    real_t polyError89 = Lpoly->lInfinityError(8, 9, 5);
//    real_t polyError910 = Lpoly->lInfinityError(9, 10, 5);
//    real_t polyError1011 = Lpoly->lInfinityError(10, 11, 5);
//    real_t polyError1112 = Lpoly->lInfinityError(11, 12, 5);

//    WALBERLA_LOG_INFO_ON_ROOT("polyError34 = " << polyError34);
//    WALBERLA_LOG_INFO_ON_ROOT("polyError45 = " << polyError45);
//    WALBERLA_LOG_INFO_ON_ROOT("polyError56 = " << polyError56);
    //    WALBERLA_LOG_INFO_ON_ROOT("polyError67 = " << polyError67);
    //    WALBERLA_LOG_INFO_ON_ROOT("polyError78 = " << polyError78);
    //    WALBERLA_LOG_INFO_ON_ROOT("polyError89 = " << polyError89);
    //    WALBERLA_LOG_INFO_ON_ROOT("polyError910 = " << polyError910);
    //    WALBERLA_LOG_INFO_ON_ROOT("polyError1011 = " << polyError1011);
    //    WALBERLA_LOG_INFO_ON_ROOT("polyError1112 = " << polyError1112);
  }
  auto   end       = walberla::timing::getWcTime();
  double setupTime = end - start;

  WALBERLA_LOG_INFO_ON_ROOT( "Interpolating exact function" );
  u_exact.interpolate( exact, maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( "Integrating rhs" );
  tmp.interpolate( rhs, maxLevel );
  M.apply( tmp, f, maxLevel, hyteg::All );

  npoints_helper.interpolate( ones, maxLevel );
  real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel );

  //  npoints_helper.interpolate(ones, interpolationLevel);
  //  real_t npointsCoarse = npoints_helper.dotGlobal(npoints_helper, interpolationLevel);

  auto coarseLaplaceSolverPoly =
      std::make_shared< hyteg::CGSolver< SolveOperatorPoly > >( storage, minLevel, minLevel, max_cg_iter, coarse_tolerance );
  auto coarseLaplaceSolverNodal =
      std::make_shared< hyteg::CGSolver< SolveOperatorNodal > >( storage, minLevel, minLevel, max_cg_iter, coarse_tolerance );

  auto restrictionOperator           = std::make_shared< hyteg::P1toP1LinearRestriction<> >();
  auto prolongationOperator          = std::make_shared< hyteg::P1toP1LinearProlongation<> >();
  auto quadraticProlongationOperator = std::make_shared< hyteg::P1toP1QuadraticProlongation >();

  auto smootherPoly  = std::make_shared< hyteg::GaussSeidelSmoother< SolveOperatorPoly > >();
  auto smootherNodal = std::make_shared< hyteg::GaussSeidelSmoother< SolveOperatorNodal > >();

  GeometricMultigridSolver< SolveOperatorPoly > laplaceSolverPoly(
      storage, smootherPoly, coarseLaplaceSolverPoly, restrictionOperator, prolongationOperator, minLevel, maxMemoryLevel, 2, 2 );
  GeometricMultigridSolver< SolveOperatorNodal > laplaceSolverNodal( storage,
                                                                     smootherNodal,
                                                                     coarseLaplaceSolverNodal,
                                                                     restrictionOperator,
                                                                     prolongationOperator,
                                                                     minLevel,
                                                                     maxMemoryLevel,
                                                                     2,
                                                                     2 );

  WALBERLA_LOG_INFO_ON_ROOT( "Starting V cycles" );
  WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%10s|%10s|%10s|%10s|%10s|%10s|%10s",
                                               "iter",
                                               "abs_res",
                                               "rel_res",
                                               "conv",
                                               "L2-error",
                                               "est. L2",
                                               "Cycle-Time",
                                               "Est-Time" ) );

  real_t rel_res = 1.0;

  if ( polynomialOperator )
  {
    Lpoly->apply( u, Lu, maxLevel, hyteg::Inner );
  }
  else
  {
    L->apply( u, Lu, maxLevel, hyteg::Inner );
  }
  r.assign( { 1.0, -1.0 }, { f, Lu }, maxLevel, hyteg::Inner );

  real_t begin_res   = std::sqrt( r.dotGlobal( r, maxLevel, hyteg::Inner ) );
  real_t abs_res_old = begin_res;

  err.assign( { 1.0, -1.0 }, { u, u_exact }, maxLevel );
  real_t discr_l2_err = std::sqrt( err.dotGlobal( err, maxLevel ) / npoints );

  // Estimating discretization error
  quadraticProlongationOperator->prolongate( u, maxLevel, hyteg::Inner);
  r.interpolate(zeros, maxMemoryLevel, hyteg::All);
//  L->applyPartial(u, r, maxMemoryLevel, interpolationLevel, hyteg::Inner);
//  tmp.interpolate(zeros, maxMemoryLevel, hyteg::All);
//  L->smooth_gs(tmp, r, maxMemoryLevel, hyteg::Inner);
//  real_t estL2Error = std::sqrt(r.dotGlobal(r, maxMemoryLevel) / npointsCoarse);
//  real_t estL2ErrorOld = estL2Error;
  real_t estL2Error = 0;

  WALBERLA_LOG_INFO_ON_ROOT(walberla::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res/abs_res_old, discr_l2_err,estL2Error,0.0));

  double solveTime = 0.0;
  real_t averageConvergenceRate = real_c(0.0);
  const uint_t convergenceStartIter = 3;

  uint_t i = 0;
  for (; i < max_outer_iter; ++i)
  {
    start = walberla::timing::getWcTime();
    if (polynomialOperator) {
      laplaceSolverPoly.solve(*Lpoly, u, f, maxLevel);
    } else {
      laplaceSolverNodal.solve(*L, u, f, maxLevel);
    }
    end = walberla::timing::getWcTime();
    double vCycleTime = end - start;

    start = walberla::timing::getWcTime();
    // Estimating discretization error
    quadraticProlongationOperator->prolongate( u, maxLevel, hyteg::Inner);
    r.interpolate(zeros, maxMemoryLevel, hyteg::All);
//    L->applyPartial(u, r, maxMemoryLevel, interpolationLevel, hyteg::Inner);
//    tmp.interpolate(zeros, maxMemoryLevel, hyteg::All);
//    L->smooth_gs(tmp, r, maxMemoryLevel, hyteg::Inner);
//    estL2Error = std::sqrt(r.dotGlobal(r, maxMemoryLevel) / npointsCoarse);
    end = walberla::timing::getWcTime();
    double estimatorTime = end - start;
    if (polynomialOperator) {
      Lpoly->apply(u, Lu, maxLevel, hyteg::Inner);
    } else {
      L->apply(u, Lu, maxLevel, hyteg::Inner);
    }
    r.assign({1.0, -1.0}, { f, Lu }, maxLevel, hyteg::Inner);
    real_t abs_res = std::sqrt(r.dotGlobal(r, maxLevel, hyteg::Inner));
    rel_res = abs_res / begin_res;
    err.assign({1.0, -1.0}, { u, u_exact }, maxLevel);
    discr_l2_err = std::sqrt(err.dotGlobal(err, maxLevel) / npoints);

    WALBERLA_LOG_INFO_ON_ROOT(walberla::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", i+1, abs_res, rel_res, abs_res/abs_res_old, discr_l2_err, estL2Error, vCycleTime, estimatorTime));

#if 0
    if (polynomialOperator) {
      if (estL2Error / estL2ErrorOld > 0.99 && useDegree < 12) {

        if (updatedDegree && abs_res/abs_res_old <= 1.0) {
          WALBERLA_LOG_INFO_ON_ROOT("Increasing polynomial had no effect, finishing");
          break;
        }

        WALBERLA_LOG_INFO_ON_ROOT("Increasing polynomial degree to " << useDegree + 1);
        ++useDegree;
        Lpoly->interpolateStencils(useDegree, interpolationLevel);
        updatedDegree = true;
      } else {
        updatedDegree = false;
      }
    } else {
      if (estL2Error / estL2ErrorOld > 0.99) {
        WALBERLA_LOG_INFO_ON_ROOT("Error estimator converged");
        break;
      }
    }
    solveTime += vCycleTime + estimatorTime;
#endif
    solveTime += vCycleTime;

    if (i >= convergenceStartIter) {
      averageConvergenceRate += abs_res/abs_res_old;
    }

//    estL2ErrorOld = estL2Error;
    abs_res_old = abs_res;

    if (rel_res < mg_tolerance)
    {
      break;
    }
  }

  WALBERLA_LOG_INFO_ON_ROOT("L(H): " << level_H);
  WALBERLA_LOG_INFO_ON_ROOT("L(h): " << level_h);
  if (polynomialOperator) {
    WALBERLA_LOG_INFO_ON_ROOT("Polynomial degree: " << maxPolyDegree);
    WALBERLA_LOG_INFO_ON_ROOT("LSQP level: " << interpolationLevel);
  }
  WALBERLA_LOG_INFO_ON_ROOT("Setup time: " << std::defaultfloat << setupTime);
  WALBERLA_LOG_INFO_ON_ROOT("Solve time " << std::defaultfloat << solveTime);
  WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << std::defaultfloat << setupTime + solveTime);
  WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(i-convergenceStartIter));
  WALBERLA_LOG_INFO_ON_ROOT("L^2 error: " << std::scientific << discr_l2_err);
  WALBERLA_LOG_INFO_ON_ROOT("DoFs: " << (uint_t) npoints);

  if (parameters.getParameter<bool>("vtkOutput")) {
    hyteg::VTKOutput vtkOutput("../output", "gmg_blending", storage);
    vtkOutput.add(u);
//    vtkOutput.add(&u_fe);
    vtkOutput.add(err);
//    vtkOutput.add(&err_est);
    vtkOutput.add(r);
//    vtkOutput.add(&r_fe);
//    vtkOutput.add(&tmp);
//    vtkOutput.add(&tmp_fe);
    vtkOutput.write(maxLevel, 0);
//    vtkOutput.write(maxMemoryLevel, 0);
  }

  return 0;
}
