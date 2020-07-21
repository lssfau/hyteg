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

// ==================== FE-spaces ====================

enum ElementType
{
  P1    = 1,
  P2    = 2
};

enum StencilType
{
  NONE      = -1,
  // CONST     = 0,
  VARIABLE  = 1,
  LSQP      = 2
};

template <StencilType T>
struct OperatorHandler
{
  static uint_t interpolationLevel;

  template <class OP>
  static std::shared_ptr<OP> make_shared(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel)
  {
    return std::make_shared<OP>(storage, minLevel, maxLevel);
  }

  template <class OP>
  static real_t setup(std::shared_ptr<OP>, const uint_t) {return 0;}
};

template<>
uint_t OperatorHandler<StencilType::LSQP>::interpolationLevel = 0;

template<>
template<class OP>
std::shared_ptr<OP> OperatorHandler<StencilType::LSQP>::make_shared(
                            std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel)
{
  return std::make_shared<OP>(storage, minLevel, maxLevel, interpolationLevel);
}

template<>
template <class OP>
real_t OperatorHandler<StencilType::LSQP>::setup(std::shared_ptr<OP> op, const uint_t polyDegree)
{
  auto start = walberla::timing::getWcTime();
  op->interpolateStencils(polyDegree);
  auto end = walberla::timing::getWcTime();
  op->useDegree(polyDegree);
  return (end - start);
}

template <ElementType P, StencilType T>
struct FE_Space{};

// template <>
// struct FE_Space<ElementType::P1, StencilType::NONE>
// {
//   using Function = hyteg::P1Function<real_t>;
//   using Restriction = hyteg::P1toP1LinearRestriction;
//   using Prolongation = hyteg::P1toP1LinearProlongation;

//   using Mass = hyteg::P1BlendingMassOperator;

//   using LaplaceCONST = hyteg::P1ConstantLaplaceOperator;
//   using LaplaceVAR = hyteg::P1BlendingLaplaceOperator;
//   using LaplaceLSQP = hyteg::P1PolynomialBlendingLaplaceOperator;
// };

template <>
struct FE_Space<ElementType::P2, StencilType::NONE>
{
  using Function = hyteg::P2Function<real_t>;
  using Restriction = hyteg::P2toP2QuadraticRestriction;
  using Prolongation = hyteg::P2toP2QuadraticProlongation;

  using Mass = hyteg::P2ConstantMassOperator;

  using LaplaceCONST = hyteg::P2ConstantLaplaceOperator;
  using LaplaceVAR = hyteg::P2BlendingLaplaceOperator;
  using LaplaceLSQP = hyteg::P2PolynomialBlendingLaplaceOperator;

  using divKgradVAR = hyteg::P2divKgradOperator;
  using divKgradLSQP = hyteg::P2PolynomialDivKgradOperator;
};

template <ElementType P>
using P_Space = FE_Space<P, StencilType::NONE>;

template<ElementType P>
struct FE_Space<P, StencilType::VARIABLE> : public P_Space<P>, public OperatorHandler<StencilType::VARIABLE>
{
  using typename P_Space<P>::Function;
  using typename P_Space<P>::Restriction;
  using typename P_Space<P>::Prolongation;

  using typename P_Space<P>::Mass;
  using Laplace = typename P_Space<P>::LaplaceVAR;
  using divKgrad = typename P_Space<P>::divKgradVAR;
};

template<ElementType P>
struct FE_Space<P,StencilType::LSQP> : public P_Space<P>, public OperatorHandler<StencilType::LSQP>
{
  using typename P_Space<P>::Function;
  using typename P_Space<P>::Restriction;
  using typename P_Space<P>::Prolongation;

  using typename P_Space<P>::Mass;
  using Laplace = typename P_Space<P>::LaplaceLSQP;
  using divKgrad = typename P_Space<P>::divKgradLSQP;

  static void setInterpolationLevel(uint_t level)
  {
    interpolationLevel = level;
  };
};

// template<ElementType P>
// struct FE_Space<P,StencilType::CONST> : public P_Space<P>, public OperatorHandler<StencilType::CONST>
// {
//   using typename P_Space<P>::Function;
//   using typename P_Space<P>::Restriction;
//   using typename P_Space<P>::Prolongation;

//   using typename P_Space<P>::Mass;
//   using Laplace = typename P_Space<P>::LaplaceCONST;
// };

// ================================================

template <ElementType P, StencilType T>
void solveTmpl(std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel
           , const uint_t max_outer_iter, const uint_t max_cg_iter, const real_t mg_tolerance, const real_t coarse_tolerance, const bool vtk
           , c_function& exact, c_function& boundary, c_function& rhs, const uint_t minPolyDegree = 0, const uint_t maxPolyDegree = 0)
{
  using FE = FE_Space<P,T>;

  // define functions and operators
  typename FE::Function r("r", storage, minLevel, maxLevel);
  typename FE::Function f("f", storage, minLevel, maxLevel);
  typename FE::Function u("u", storage, minLevel, maxLevel);
  typename FE::Function Lu("Lu", storage, minLevel, maxLevel);
  typename FE::Function u_exact("u_exact", storage, minLevel, maxLevel);
  typename FE::Function err("err", storage, minLevel, maxLevel);
  typename FE::Function tmp("tmp", storage, minLevel, maxLevel);

  // rhs
  typename FE::Mass M(storage, minLevel, maxLevel);
  tmp.interpolate(rhs, maxLevel);
  M.apply(tmp, f, maxLevel, hyteg::All);

  // operator
  auto L = FE::template make_shared<typename FE::divKgrad>(storage, minLevel, maxLevel);

  // exact solution
  u_exact.interpolate(exact, maxLevel);

  // boundary
  u.interpolate(boundary, maxLevel, hyteg::DirichletBoundary);

    // define solver
  auto coarseGridSolver = std::make_shared<hyteg::CGSolver<typename FE::divKgrad>>(storage, minLevel, minLevel, max_cg_iter, coarse_tolerance);
  auto restrictionOperator = std::make_shared<typename FE::Restriction>();
  auto prolongationOperator = std::make_shared<typename FE::Prolongation>();
  auto smoother = std::make_shared<hyteg::GaussSeidelSmoother<typename FE::divKgrad>>();
  GeometricMultigridSolver<typename FE::divKgrad> GMGSolver(storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3);

  // solve for each polynomial degree
  for (uint_t polyDegree = minPolyDegree; polyDegree <= maxPolyDegree; ++polyDegree)
  {
    if (T == StencilType::LSQP)
    {
      WALBERLA_LOG_INFO_ON_ROOT(walberla::format("Apply LSQ-fit with q = %d", polyDegree));
    }
    real_t setupTime = FE::setup(L, polyDegree);

    // reset u
    u.interpolate([](const hyteg::Point3D&) {return 0.0;}, maxLevel, hyteg::Inner);

    // solve iteratively
    uint_t iter = 0;
    real_t res = 0, res_old, discr_l2_err;
    real_t vCycleTime, solveTime = 0;
    real_t averageConvergenceRate = 0;
    const uint_t convergenceStartIter = 3;

    WALBERLA_LOG_INFO_ON_ROOT("Starting V cycles");
    WALBERLA_LOG_INFO_ON_ROOT(walberla::format("%6s|%10s|%10s|%10s|%10s", "iter", "res", "conv", "L2-error", "Cycle-Time"));

    while(1)
    {
       // compute error
      err.assign({1.0, -1.0}, {u, u_exact}, maxLevel);
      M.apply(err, tmp, maxLevel, hyteg::All, Replace);
      discr_l2_err = std::sqrt(err.dotGlobal(tmp, maxLevel));

      // compute residual
      res_old = res;
      L->apply(u, Lu, maxLevel, hyteg::Inner, Replace);
      r.assign({1.0, -1.0}, {f, Lu}, maxLevel, hyteg::Inner);
      M.apply(r, tmp, maxLevel, hyteg::All, Replace);
      res = std::sqrt(r.dotGlobal(tmp, maxLevel, hyteg::Inner));

      // compute convergence rate
      real_t convRate = res / res_old;
      if (iter >= convergenceStartIter)
      {
        averageConvergenceRate += convRate;
      }

      WALBERLA_LOG_INFO_ON_ROOT(walberla::format("%6d|%10.3e|%10.3e|%10.3e|%10.3e", iter, res, convRate, discr_l2_err, vCycleTime));

      // stopping criterion
      if (++iter > max_outer_iter || (iter > convergenceStartIter && convRate > 0.999))
      {
        WALBERLA_LOG_INFO_ON_ROOT("Ending multigrid without reaching desired tolerance!");
        break;
      }
      if (res < mg_tolerance)
      {
        WALBERLA_LOG_INFO_ON_ROOT("Multigrid converged!");
        break;
      }


      // solve
      auto start = walberla::timing::getWcTime();
      GMGSolver.solve(*L, u, f, maxLevel);
      auto end = walberla::timing::getWcTime();
      vCycleTime = end - start;
      solveTime += vCycleTime;
    }

    WALBERLA_LOG_INFO_ON_ROOT("Setup time: " << std::defaultfloat << setupTime);
    WALBERLA_LOG_INFO_ON_ROOT("Solve time " << std::defaultfloat << solveTime);
    WALBERLA_LOG_INFO_ON_ROOT("Time to solution: " << std::defaultfloat << setupTime + solveTime);
    WALBERLA_LOG_INFO_ON_ROOT("Avg. convergence rate: " << std::scientific << averageConvergenceRate / real_c(iter - convergenceStartIter));
    WALBERLA_LOG_INFO_ON_ROOT("L^2 error: " << std::scientific << discr_l2_err);
    // WALBERLA_LOG_INFO_ON_ROOT("DoFs: " << (uint_t) npoints);

    if (vtk)
    {
      std::string name = "Polynomial_divKgrad_P" + std::to_string(P);
      switch (T)
      {
        // case StencilType::CONST:
        //   name += "const";
        //   break;

        case StencilType::VARIABLE:
          name += "var";
          break;

        case StencilType::LSQP:
          name += "poly" + std::to_string(polyDegree);
          break;
      }

      hyteg::VTKOutput vtkOutput("../output", name, storage);
      vtkOutput.add(u_exact);
      vtkOutput.add(u);
      vtkOutput.add(err);
      vtkOutput.add(r);
      vtkOutput.write(maxLevel, 0);
    }

    if (T != StencilType::LSQP) break;
    WALBERLA_LOG_INFO_ON_ROOT("=======================================================");
  }
}

template <ElementType P>
void solve(const StencilType T, const uint_t interpolationLevel, std::shared_ptr<PrimitiveStorage> storage, const uint_t minLevel, const uint_t maxLevel
           , const uint_t max_outer_iter, const uint_t max_cg_iter, const real_t mg_tolerance, const real_t coarse_tolerance, const bool vtk
           , c_function& exact, c_function& boundary, c_function& rhs, const uint_t minPolyDegree, const uint_t maxPolyDegree)
{
  switch (T)
  {
    case VARIABLE:
      WALBERLA_LOG_INFO_ON_ROOT("Operatortype: Variable Stencil");
      solveTmpl<P, VARIABLE>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
      break;

    case LSQP:
      WALBERLA_LOG_INFO_ON_ROOT("Operatortype: Surrogate Polynomial Stencil");
      WALBERLA_LOG_INFO_ON_ROOT("Interpolation level: " << interpolationLevel);

      FE_Space<P,LSQP>::setInterpolationLevel(interpolationLevel);
      solveTmpl<P, LSQP>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs, minPolyDegree, maxPolyDegree);
      break;

    default:
      WALBERLA_ABORT("The desired Operator Type is not supported!");
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
    cfg_->readParameterFile("../data/param/VariableCoeffLaplace_Surrogates.prm");
    cfg = cfg_;
  }
  else
  {
    cfg = walberla::config::create(argc, argv);
  }

  // read parameter file
  // WALBERLA_LOG_INFO_ON_ROOT("config = " << *cfg);
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  const uint_t minLevel = parameters.getParameter<uint_t>("level_h_coarse");
  const uint_t maxLevel = parameters.getParameter<uint_t>("level_h_fine");

  const int elTypeInt = parameters.getParameter<int>("elementType");
  const ElementType elType = ElementType(elTypeInt);

  const int opTypeInt = parameters.getParameter<int>("operatorType");
  const StencilType opType = StencilType(opTypeInt);

  const uint_t minPolyDegree = parameters.getParameter<uint_t>("minPolyDegree");
  const uint_t maxPolyDegree = parameters.getParameter<uint_t>("maxPolyDegree");
  const uint_t maxInterpolationLevel = parameters.getParameter<uint_t>("interpolationLevel");
  const uint_t interpolationLevel = std::min(maxLevel, maxInterpolationLevel);

  const uint_t max_outer_iter =  parameters.getParameter<uint_t>("max_outer_iter");
  const uint_t max_cg_iter =  parameters.getParameter<uint_t>("max_cg_iter");
  const real_t mg_tolerance = parameters.getParameter<real_t>("mg_tolerance");
  const real_t coarse_tolerance = parameters.getParameter<real_t>("coarse_tolerance");

  const bool vtk = parameters.getParameter<bool>("vtkOutput");

  const bool discontinuousK = parameters.getParameter<bool>("discontinuous_k");

  // define functions and domain
  if (discontinuousK)
  {
    WALBERLA_LOG_INFO_ON_ROOT("k discontinuous (main-singularity)");
  }
  else
  {
    WALBERLA_LOG_INFO_ON_ROOT("smooth k, single macro-element");
  }
  // case: smooth k, domain = triangle
  c_function exact = [](const hyteg::Point3D & x) {return sin(x[0])*sinh(x[1]);};
  c_function boundary = exact;
  c_function rhs = [](const hyteg::Point3D & x) {return - cos(x[0])*cos(x[0]) * sinh(x[1]);};
  c_function k = [](const hyteg::Point3D & x) {return 2.0 +  sin(x[0]);};
  MeshInfo meshInfo = MeshInfo::fromGmshFile("../data/meshes/tri_1el.msh");

  // case discontinuous k, domain = square (main-singularity)
  if (discontinuousK)
  {
    k = [](const hyteg::Point3D & x) {return (x[0]*x[1] < 0)? 1 : 5;};
    rhs = [](const hyteg::Point3D &) {return 0;};
    real_t a[4] = {0.4472135955,-0.7453559925,-0.9441175905,-2.401702643};
    real_t b[4] = {1.0,2.333333333,0.55555555555,-0.4814814814};
    real_t delta = 0.5354409456;
    exact = [a,b,delta](const hyteg::Point3D & x){
      // domain partition
      int i;
      if (x[1] >= 0)
      {
        i = (x[0] > 0)? 0 : 1;
      }
      else
      {
        i = (x[0] < 0)? 2 : 3;
      }

      // radius
      real_t r = std::sqrt(x[0]*x[0]+x[1]*x[1]);

      // angle
      real_t theta;
      if (std::abs(x[0]) < 1E-16)
      {
        theta = (x[1] > 0)? PI/2. : -PI/2.;
      }
      else
      {
        theta = std::atan(x[1]/x[0]);
        if (x[0] < 0)
          theta += PI;
      }
      if (theta < 0) theta += 2*PI;

      return std::pow(r,delta) * (a[i]*sin(delta*theta) + b[i]*cos(delta*theta));
    };
    boundary = exact;

    meshInfo = MeshInfo::meshRectangle(Point2D({-1.0, -1.0}), Point2D({1.0, 1.0}), MeshInfo::CRISS, 1, 1);
  }

  P2Form_divKgrad::callback = k;

  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);

  hyteg::loadbalancing::roundRobin(setupStorage);
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  // size_t N = storage->getFaces().size();
  WALBERLA_LOG_INFO_ON_ROOT("Refinement levels: " << minLevel << "->" << maxLevel);
  // WALBERLA_LOG_INFO_ON_ROOT("Number of Microfaces:" << N * (1 << (2*maxLevel)));
  WALBERLA_LOG_INFO_ON_ROOT("h_L = 1/" << (1 << maxLevel));

  // choose FE space
  switch (elType)
  {
    case ElementType::P1:
      WALBERLA_ABORT("Operator not implemented for P1 elements!")
      // WALBERLA_LOG_INFO_ON_ROOT("Element Type: P1");
      // solve<ElementType::P1>(opType, interpolationLevel, storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs, minPolyDegree, maxPolyDegree);
      break;

    case ElementType::P2:
      WALBERLA_LOG_INFO_ON_ROOT("Element Type: P2");
      solve<ElementType::P2>(opType, interpolationLevel, storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs, minPolyDegree, maxPolyDegree);
      break;

    default:
      WALBERLA_ABORT("The desired FE space is not supported!");
  }

  return 0;
}
