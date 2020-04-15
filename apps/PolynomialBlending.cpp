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
#include "hyteg/p1functionspace/P1PolynomialBlendingOperator.hpp"

#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VariableOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2PolynomialBlendingOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"

#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"

#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "core/Format.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;

#define PI 3.141592653589793238

using namespace hyteg;


typedef std::function<real_t(const hyteg::Point3D&)> c_function;

// FE-spaces

struct P1Space
{
  typedef hyteg::P1Function<real_t> Function;
  typedef hyteg::P1BlendingMassOperator Mass;
  typedef hyteg::P1BlendingLaplaceOperator Laplace;
  typedef hyteg::P1toP1LinearRestriction Restriction;
  typedef hyteg::P1toP1LinearProlongation Prolongation;

  static std::shared_ptr<Laplace> laplaceOperator(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel)
  {
    return std::make_shared<Laplace>(storage, minLevel, maxLevel);
  }
};

struct P1Space_LSQP
{
  static uint_t polyDegree;
  static uint_t interpolationLevel;
  typedef hyteg::P1Function<real_t> Function;
  typedef hyteg::P1BlendingMassOperator Mass;
  typedef hyteg::P1PolynomialBlendingLaplaceOperator Laplace;
  typedef hyteg::P1toP1LinearRestriction Restriction;
  typedef hyteg::P1toP1LinearProlongation Prolongation;

  static std::shared_ptr<Laplace> laplaceOperator(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel)
  {
    return std::make_shared<Laplace>(storage, minLevel, maxLevel, interpolationLevel, polyDegree);
  }
};
uint_t P1Space_LSQP::polyDegree = 0;
uint_t P1Space_LSQP::interpolationLevel = 0;

struct P2Space
{
  typedef hyteg::P2Function<real_t> Function;
  typedef hyteg::P2BlendingMassOperator Mass;
  typedef hyteg::P2BlendingLaplaceOperator Laplace;
  typedef hyteg::P2toP2QuadraticRestriction Restriction;
  typedef hyteg::P2toP2QuadraticProlongation Prolongation;

  static std::shared_ptr<Laplace> laplaceOperator(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel)
  {
    return std::make_shared<Laplace>(storage, minLevel, maxLevel);
  }
};

struct P2Space_LSQP
{
  static uint_t polyDegree;
  static uint_t interpolationLevel;
  typedef hyteg::P2Function<real_t> Function;
  typedef hyteg::P2BlendingMassOperator Mass;
  typedef hyteg::P2PolynomialBlendingLaplaceOperator Laplace;
  typedef hyteg::P2toP2QuadraticRestriction Restriction;
  typedef hyteg::P2toP2QuadraticProlongation Prolongation;

  static std::shared_ptr<Laplace> laplaceOperator(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel)
  {
    return std::make_shared<Laplace>(storage, minLevel, maxLevel, interpolationLevel, polyDegree);
  }
};
uint_t P2Space_LSQP::polyDegree = 0;
uint_t P2Space_LSQP::interpolationLevel = 0;

struct P2Space_const
{
  typedef hyteg::P2Function<real_t> Function;
  typedef hyteg::P2ConstantMassOperator Mass;
  typedef hyteg::P2ConstantLaplaceOperator Laplace;
  typedef hyteg::P2toP2QuadraticRestriction Restriction;
  typedef hyteg::P2toP2QuadraticProlongation Prolongation;

  static std::shared_ptr<Laplace> laplaceOperator(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel)
  {
    return std::make_shared<Laplace>(storage, minLevel, maxLevel);
  }
};

struct P2Space_elementwise
{
  typedef hyteg::P2Function<real_t> Function;
  typedef hyteg::P2ElementwiseBlendingMassOperator Mass;
  typedef hyteg::P2ElementwiseBlendingLaplaceOperator Laplace;
  typedef hyteg::P2toP2QuadraticRestriction Restriction;
  typedef hyteg::P2toP2QuadraticProlongation Prolongation;

  static std::shared_ptr<Laplace> laplaceOperator(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel)
  {
    return std::make_shared<Laplace>(storage, minLevel, maxLevel);
  }
};

// compare variable and constant P2 operators (debug purpose)
void compareP2(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel
           , const uint_t max_outer_iter, const uint_t max_cg_iter, const real_t coarse_tolerance, const bool vtk
           , c_function& exact, c_function& boundary, c_function& rhs)
{
  // define functions and operators

  c_function zeros = [](const hyteg::Point3D&) {return 0.0;};
  c_function ones  = [](const hyteg::Point3D&) {return 1.0;};
  c_function xExpr = [](const hyteg::Point3D & x) {return x[0];};
  c_function yExpr = [](const hyteg::Point3D & x) {return x[1];};

  typename P2Space::Function r("r", storage, minLevel, maxLevel);
  typename P2Space::Function f_var("f_var", storage, minLevel, maxLevel);
  typename P2Space::Function f_const("f_const", storage, minLevel, maxLevel);
  typename P2Space::Function f_ew("f_ew", storage, minLevel, maxLevel);
  typename P2Space::Function u_var("u_var", storage, minLevel, maxLevel);
  typename P2Space::Function u_const("u_const", storage, minLevel, maxLevel);
  typename P2Space::Function u_ew("u_ew", storage, minLevel, maxLevel);
  typename P2Space::Function L_var_u("L_var_u", storage, minLevel, maxLevel);
  typename P2Space::Function L_const_u("L_const_u", storage, minLevel, maxLevel);
  typename P2Space::Function L_ew_u("L_ew_u", storage, minLevel, maxLevel);
  typename P2Space::Function u_exact("u_exact", storage, minLevel, maxLevel);
  typename P2Space::Function err_c_M("err_c_M", storage, minLevel, maxLevel);
  typename P2Space::Function err_c_L("err_c_L", storage, minLevel, maxLevel);
  typename P2Space::Function err_c_solver("err_c_solver", storage, minLevel, maxLevel);
  typename P2Space::Function err_ew_M("err_ew_M", storage, minLevel, maxLevel);
  typename P2Space::Function err_ew_L("err_ew_L", storage, minLevel, maxLevel);
  typename P2Space::Function err_ew_solver("err_ew_solver", storage, minLevel, maxLevel);
  typename P2Space::Function err_var("err_var", storage, minLevel, maxLevel);
  typename P2Space::Function err_const("err_const", storage, minLevel, maxLevel);
  typename P2Space::Function err_ew("err_ew", storage, minLevel, maxLevel);
  typename P2Space::Function one("1", storage, minLevel, maxLevel);
  typename P2Space::Function tmp("tmp", storage, minLevel, maxLevel);
  auto coordX = std::make_shared<typename P2Space::Function>("x", storage, minLevel, maxLevel);
  auto coordY = std::make_shared<typename P2Space::Function>("y", storage, minLevel, maxLevel);

  coordX->interpolate(xExpr, maxLevel, hyteg::All);
  coordY->interpolate(yExpr, maxLevel, hyteg::All);

  one.interpolate(ones, maxLevel);

  u_const.interpolate(boundary, maxLevel, hyteg::DirichletBoundary);
  u_var.interpolate(boundary, maxLevel, hyteg::DirichletBoundary);
  u_ew.interpolate(boundary, maxLevel, hyteg::DirichletBoundary);

  typename P2Space::Mass M_var(storage, minLevel, maxLevel);
  typename P2Space_const::Mass M_const(storage, minLevel, maxLevel);
  typename P2Space_elementwise::Mass M_ew(storage, minLevel, maxLevel);

  auto L_var = P2Space::laplaceOperator(storage, minLevel, maxLevel);
  auto L_const = P2Space_const::laplaceOperator(storage, minLevel, maxLevel);
  auto L_ew = P2Space_elementwise::laplaceOperator(storage, minLevel, maxLevel);

  u_exact.interpolate(exact, maxLevel);
  tmp.interpolate(rhs, maxLevel);

  M_const.apply(tmp, f_const, maxLevel, hyteg::All);
  M_var.apply(tmp, f_var, maxLevel, hyteg::All);
  M_ew.apply(tmp, f_ew, maxLevel, hyteg::All);
  err_c_M.assign({1.0, -1.0}, {f_const, f_var}, maxLevel);
  err_ew_M.assign({1.0, -1.0}, {f_ew, f_var}, maxLevel);

  L_var->apply(u_var, L_var_u, maxLevel, hyteg::Inner);
  L_const->apply(u_const, L_const_u, maxLevel, hyteg::Inner);
  L_ew->apply(u_ew, L_ew_u, maxLevel, hyteg::Inner);
  err_c_L.assign({1.0, -1.0}, {L_const_u, L_var_u}, maxLevel);
  err_ew_L.assign({1.0, -1.0}, {L_ew_u, L_var_u}, maxLevel);

  real_t npoints = one.dotGlobal(one, maxLevel);

  // define solver

  auto cg_var = std::make_shared<hyteg::CGSolver<typename P2Space::Laplace>>(storage, minLevel, minLevel, max_cg_iter, coarse_tolerance);
  auto cg_const = std::make_shared<hyteg::CGSolver<typename P2Space_const::Laplace>>(storage, minLevel, minLevel, max_cg_iter, coarse_tolerance);
  auto restrictionOperator = std::make_shared<typename P2Space::Restriction>();
  auto prolongationOperator = std::make_shared<typename P2Space::Prolongation>();
  auto smoother_var = std::make_shared<hyteg::GaussSeidelSmoother<typename P2Space::Laplace>>();
  auto smoother_const = std::make_shared<hyteg::GaussSeidelSmoother<typename P2Space_const::Laplace>>();
  GeometricMultigridSolver<typename P2Space::Laplace> laplaceSolver_var(storage, smoother_var, cg_var, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2);
  GeometricMultigridSolver<typename P2Space_const::Laplace> laplaceSolver_const(storage, smoother_const, cg_const, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2);

  auto CGsolver_var = std::make_shared<hyteg::CGSolver<typename P2Space::Laplace>>(storage, maxLevel, maxLevel, max_cg_iter, coarse_tolerance);
  auto CGsolver_const = std::make_shared<hyteg::CGSolver<typename P2Space_const::Laplace>>(storage, maxLevel, maxLevel, max_cg_iter, coarse_tolerance);
  auto CGsolver_ew = std::make_shared<hyteg::CGSolver<typename P2Space_elementwise::Laplace>>(storage, maxLevel, maxLevel, max_cg_iter, coarse_tolerance);

  CGsolver_var->solve(*L_var,u_var,f_var,maxLevel);
  CGsolver_const->solve(*L_const,u_const,f_const,maxLevel);
  CGsolver_ew->solve(*L_ew,u_ew,f_ew,maxLevel);

  for (uint_t k = 0; k < max_outer_iter; ++k)
  {
  //   // smoother_const->solve(*L_const, u_const, f_const, maxLevel);
  //   // smoother_var->solve(*L_var, u_var, f_var, maxLevel);
  //   laplaceSolver_const.solve(*L_const, u_const, f_const, maxLevel);
  //   laplaceSolver_var.solve(*L_var, u_var, f_var, maxLevel);
  }
  err_c_solver.assign({1.0, -1.0}, {u_const, u_var}, maxLevel);
  err_ew_solver.assign({1.0, -1.0}, {u_ew, u_var}, maxLevel);

  // error
  WALBERLA_LOG_INFO_ON_ROOT(walberla::format("error Mass: const: %10.3e,  elwise: %10.3e", std::sqrt(err_c_M.dotGlobal(err_c_M, maxLevel)),std::sqrt(err_ew_M.dotGlobal(err_ew_M, maxLevel))));
  WALBERLA_LOG_INFO_ON_ROOT(walberla::format("error Lapl: const: %10.3e,  elwise: %10.3e", std::sqrt(err_c_L.dotGlobal(err_c_L, maxLevel)),std::sqrt(err_ew_L.dotGlobal(err_ew_L, maxLevel))));
  WALBERLA_LOG_INFO_ON_ROOT(walberla::format("error Solver: const: %10.3e,  elwise: %10.3e", std::sqrt(err_c_solver.dotGlobal(err_c_solver, maxLevel)),std::sqrt(err_ew_solver.dotGlobal(err_ew_solver, maxLevel))));

  err_var.assign({1.0, -1.0}, {u_var, u_exact}, maxLevel);
  real_t L2_err_var = std::sqrt(err_var.dotGlobal(err_var, maxLevel) / npoints);
  err_const.assign({1.0, -1.0}, {u_const, u_exact}, maxLevel);
  real_t L2_err_const = std::sqrt(err_const.dotGlobal(err_const, maxLevel) / npoints);
  err_ew.assign({1.0, -1.0}, {u_ew, u_exact}, maxLevel);
  real_t L2_err_ew = std::sqrt(err_ew.dotGlobal(err_ew, maxLevel) / npoints);

  WALBERLA_LOG_INFO_ON_ROOT(walberla::format("discrete L2 error: var: %10.3e,  const: %10.3e,  elwise: %10.3e",L2_err_var,L2_err_const,L2_err_ew ));


  if (vtk)
  {
    hyteg::VTKOutput vtkOutput("../output", "PolynomialBlending", storage);
    vtkOutput.add(err_c_M);
    vtkOutput.add(err_c_L);
    vtkOutput.add(err_ew_M);
    vtkOutput.add(err_ew_L);
    vtkOutput.add(err_var);
    vtkOutput.add(err_const);
    vtkOutput.add(err_ew);
    vtkOutput.write(maxLevel, 0);
  }
}

template<typename FE>
void solve(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel
           , const uint_t max_outer_iter, const uint_t max_cg_iter, const real_t mg_tolerance, const real_t coarse_tolerance, const bool vtk
           , c_function& exact, c_function& boundary, c_function& rhs)
{
  // define functions and operators

  c_function zeros = [](const hyteg::Point3D&) {return 0.0;};
  c_function ones  = [](const hyteg::Point3D&) {return 1.0;};
  c_function xExpr = [](const hyteg::Point3D & x) {return x[0];};
  c_function yExpr = [](const hyteg::Point3D & x) {return x[1];};

  typename FE::Function r("r", storage, minLevel, maxLevel);
  typename FE::Function f("f", storage, minLevel, maxLevel);
  typename FE::Function u("u", storage, minLevel, maxLevel);
  typename FE::Function Lu("Lu", storage, minLevel, maxLevel);
  typename FE::Function u_exact("u_exact", storage, minLevel, maxLevel);
  typename FE::Function err("err", storage, minLevel, maxLevel);
  typename FE::Function npoints_helper("npoints_helper", storage, minLevel, maxLevel);
  typename FE::Function tmp("tmp", storage, minLevel, maxLevel);
  auto coordX = std::make_shared<typename FE::Function>("x", storage, minLevel, maxLevel);
  auto coordY = std::make_shared<typename FE::Function>("y", storage, minLevel, maxLevel);

  coordX->interpolate(xExpr, maxLevel, hyteg::All);
  coordY->interpolate(yExpr, maxLevel, hyteg::All);

  WALBERLA_LOG_INFO_ON_ROOT("Interpolating boundary");
  u.interpolate(boundary, maxLevel, hyteg::DirichletBoundary);

  WALBERLA_LOG_INFO_ON_ROOT("Setting up operators");
  typename FE::Mass M(storage, minLevel, maxLevel);

  auto start = walberla::timing::getWcTime();

  auto L = FE::laplaceOperator(storage, minLevel, maxLevel);

  auto end = walberla::timing::getWcTime();
  real_t setupTime = end - start;

  WALBERLA_LOG_INFO_ON_ROOT("Interpolating exact function");
  u_exact.interpolate(exact, maxLevel);
  WALBERLA_LOG_INFO_ON_ROOT("Integrating rhs");
  tmp.interpolate(rhs, maxLevel);
  M.apply(tmp, f, maxLevel, hyteg::All);

  npoints_helper.interpolate(ones, maxLevel);
  real_t npoints = npoints_helper.dotGlobal(npoints_helper, maxLevel);


  // define solver

  auto coarseLaplaceSolver = std::make_shared<hyteg::CGSolver<typename FE::Laplace>>(storage, minLevel, minLevel, max_cg_iter, coarse_tolerance);
  auto restrictionOperator = std::make_shared<typename FE::Restriction>();
  auto prolongationOperator = std::make_shared<typename FE::Prolongation>();
  auto smoother = std::make_shared<hyteg::GaussSeidelSmoother<typename FE::Laplace>>();

  GeometricMultigridSolver<typename FE::Laplace> laplaceSolver(storage, smoother, coarseLaplaceSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 2, 2);


  // initial state

  WALBERLA_LOG_INFO_ON_ROOT("Starting V cycles");
  WALBERLA_LOG_INFO_ON_ROOT(walberla::format("%6s|%10s|%10s|%10s|%10s|%10s", "iter", "abs_res", "rel_res", "conv", "L2-error", "Cycle-Time"));

  err.assign({1.0, -1.0}, {u, u_exact}, maxLevel);
  real_t discr_l2_err = std::sqrt(err.dotGlobal(err, maxLevel) / npoints);

  L->apply(u, Lu, maxLevel, hyteg::Inner);
  r.assign({1.0, -1.0}, {f, Lu}, maxLevel, hyteg::Inner);
  real_t begin_res = std::sqrt(r.dotGlobal(r, maxLevel, hyteg::Inner));
  real_t abs_res_old = begin_res;
  real_t rel_res = 1.0;

  WALBERLA_LOG_INFO_ON_ROOT(walberla::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", 0, begin_res, rel_res, begin_res / abs_res_old, discr_l2_err,  0.0));

  real_t solveTime = real_c(0.0);
  real_t averageConvergenceRate = real_c(0.0);
  const uint_t convergenceStartIter = 3;

  uint_t i = 0;

  for (; i < max_outer_iter; ++i)
  {
    // solve

    start = walberla::timing::getWcTime();

    laplaceSolver.solve(*L, u, f, maxLevel);
    end = walberla::timing::getWcTime();
    real_t vCycleTime = end - start;

    // compute residual

    L->apply(u, Lu, maxLevel, hyteg::Inner);
    r.assign({1.0, -1.0}, { f, Lu }, maxLevel, hyteg::Inner);
    real_t abs_res = std::sqrt(r.dotGlobal(r, maxLevel, hyteg::Inner));
    rel_res = abs_res / begin_res;

    // compute error

    err.assign({1.0, -1.0}, { u, u_exact }, maxLevel);
    discr_l2_err = std::sqrt(err.dotGlobal(err, maxLevel) / npoints);

    WALBERLA_LOG_INFO_ON_ROOT(walberla::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e|%10.3e", i + 1, abs_res, rel_res, abs_res / abs_res_old, discr_l2_err, vCycleTime));

    solveTime += vCycleTime;

    // compute convergence rate

    if (i >= convergenceStartIter)
    {
      averageConvergenceRate += abs_res / abs_res_old;
    }

    abs_res_old = abs_res;

    // stopping criterion

    if (rel_res < mg_tolerance)
    {
      break;
    }
  }

  WALBERLA_LOG_INFO_ON_ROOT("Setup time: " << std::defaultfloat << setupTime);
  WALBERLA_LOG_INFO_ON_ROOT("Solve time " << std::defaultfloat << solveTime);
  WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << std::defaultfloat << setupTime + solveTime);
  WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(i - convergenceStartIter));
  WALBERLA_LOG_INFO_ON_ROOT("L^2 error: " << std::scientific << discr_l2_err);
  WALBERLA_LOG_INFO_ON_ROOT("DoFs: " << (uint_t) npoints);

  if (vtk)
  {
    hyteg::VTKOutput vtkOutput("../output", "PolynomialBlending", storage);
    vtkOutput.add(u);
    vtkOutput.add(err);
    vtkOutput.add(r);
    vtkOutput.write(maxLevel, 0);
  }
}


int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel(walberla::logging::Logging::PROGRESS);
  walberla::MPIManager::instance()->useWorldComm();

  walberla::shared_ptr<walberla::config::Config> cfg;

  if (argc == 1)
  {
    walberla::shared_ptr<walberla::config::Config> cfg_(new walberla::config::Config);
    cfg_->readParameterFile("../data/param/PolynomialBlending.prm");
    cfg = cfg_;
  }
  else
  {
    cfg = walberla::config::create(argc, argv);
  }

  WALBERLA_LOG_INFO_ON_ROOT("config = " << *cfg);
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  const uint_t minLevel = parameters.getParameter<uint_t>("level_h_coarse");
  const uint_t maxLevel = parameters.getParameter<uint_t>("level_h_fine");
  const uint_t FE_space = parameters.getParameter<uint_t>("FE_space");

  const bool polynomialOperator = parameters.getParameter<bool>("polynomialOperator");
  const bool blending = parameters.getParameter<bool>("blending");
  const bool annulus = parameters.getParameter<bool>("annulus");
  const uint_t nX = parameters.getParameter<uint_t>("nX");
  const uint_t nY = parameters.getParameter<uint_t>("nY");
  const uint_t maxPolyDegree = parameters.getParameter<uint_t>("maxPolyDegree");
  // const uint_t maxInterpolationLevel = parameters.getParameter<uint_t>("interpolationLevel");
  // const uint_t interpolationLevel = std::min(maxLevel, maxInterpolationLevel);
  const uint_t interpolationLevel = parameters.getParameter<uint_t>("interpolationLevel");

  const uint_t max_outer_iter =  parameters.getParameter<uint_t>("max_outer_iter");
  const uint_t max_cg_iter =  parameters.getParameter<uint_t>("max_cg_iter");
  const real_t mg_tolerance = parameters.getParameter<real_t>("mg_tolerance");
  const real_t coarse_tolerance = parameters.getParameter<real_t>("coarse_tolerance");

  const bool vtk = parameters.getParameter<bool>("vtkOutput");

 if (polynomialOperator)
  {
    WALBERLA_LOG_INFO_ON_ROOT("Polynomial Operator enabled");
  }
  else
  {
    WALBERLA_LOG_INFO_ON_ROOT("Polynomial Operator disabled");
  }

  if (annulus)
  {
    WALBERLA_LOG_INFO_ON_ROOT("Geometry: Annulus");
    if (blending)
    {
      WALBERLA_LOG_INFO_ON_ROOT("Geometry blending enabled");
    }
    else
    {
      WALBERLA_LOG_INFO_ON_ROOT("Geometry blending disabled");
    }
  }
  else
  {
    WALBERLA_LOG_INFO_ON_ROOT("Geometry: Rectangle");
  }

  // define functions and domain

  /// case rectangle
  c_function exact = [](const hyteg::Point3D & x) {return sin(PI*x[0])*sin(PI* x[1]);};
  c_function boundary = [](const hyteg::Point3D &) {return 0;};
  c_function rhs = [](const hyteg::Point3D & x) {return 2*PI*PI*sin(PI*x[0])*sin(PI* x[1]);};

  MeshInfo meshInfo = MeshInfo::meshRectangle(Point2D({0.0, 0.0}), Point2D({1.0, 1.0}), MeshInfo::CRISS, nX, nY);

  /// case annulus
  if (annulus)
  {
    Point3D circleCenter{{0, 0, 0}};
    real_t innerRadius = 1.0;
    real_t outerRadius = 2.0;
    real_t middle = (outerRadius + innerRadius) / 2.0;

    exact = [](const hyteg::Point3D & x) {return sin(x[0] * x[0] + x[1] * x[1]);};
    boundary = [circleCenter, innerRadius, outerRadius, middle](const hyteg::Point3D & x) {return ((x - circleCenter).norm() < middle) ? sin(innerRadius * innerRadius) : sin(outerRadius * outerRadius);};
    rhs = [](const hyteg::Point3D & x) {return - 4 * cos(x[0] * x[0] + x[1] * x[1]) + 4 * (x[0] * x[0] + x[1] * x[1]) * sin(x[0] * x[0] + x[1] * x[1]);};

    meshInfo = MeshInfo::meshAnnulus(innerRadius, outerRadius, MeshInfo::CRISS, nX, nY);
  }

  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);

  if (annulus && blending)
    AnnulusMap::setMap(setupStorage);

  hyteg::loadbalancing::roundRobin(setupStorage);
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);


  if (parameters.getParameter<bool>("compare_P2"))
  {
    compareP2(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, coarse_tolerance, vtk, exact, boundary, rhs);
  }
  else
  {
    if (polynomialOperator)
    {
      if (FE_space == 1)
      {
        P1Space_LSQP::interpolationLevel = interpolationLevel;
        P1Space_LSQP::polyDegree = maxPolyDegree;
        solve<P1Space_LSQP>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
      }
      else if (FE_space == 2)
      {
        P2Space_LSQP::interpolationLevel = interpolationLevel;
        P2Space_LSQP::polyDegree = maxPolyDegree;
        solve<P2Space_LSQP>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
      }
      else
      {
        WALBERLA_ABORT("The desired FE space is not implemented");
      }
    }
    else
    {
      if (FE_space == 1)
      {
        solve<P1Space>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
      }
      else if (FE_space == 2)
      {
        solve<P2Space>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
      }
      else
      {
        WALBERLA_ABORT("The desired FE space is not implemented");
      }
    }
  }
  return 0;
}
