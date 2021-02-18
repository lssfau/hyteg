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
#include "hyteg/p1functionspace/P1ConstantOperator_new.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator_new.hpp"
#include "hyteg/p1functionspace/P1PolynomialBlendingOperator.hpp"

#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VariableOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2SurrogateOperator.hpp"
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

using walberla::math::pi;

using namespace hyteg;


typedef std::function<real_t(const hyteg::Point3D&)> c_function;

// cartesian to polar coordinates
c_function radius = [](const hyteg::Point3D& x) {return std::hypot(x[0], x[1]);};
c_function angle = [](const hyteg::Point3D& x) {return std::atan2(x[1], x[0]);};

// ==================== FE-spaces ====================

enum ElementType
{
  P1    = 1,
  P2    = 2
};

enum StencilType
{
  NONE      = -1,
  CONST     = 0,
  VARIABLE  = 1,
  LSQP      = 2,
  CONST_NEW = 3,
  VARIABLE_NEW   = 4
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

template <>
struct FE_Space<ElementType::P1, StencilType::NONE>
{
  using Function = hyteg::P1Function<real_t>;
  using Restriction = hyteg::P1toP1LinearRestriction;
  using Prolongation = hyteg::P1toP1LinearProlongation;

  using Mass = hyteg::P1ConstantMassOperator; //! todo use new form for blending mass

  using LaplaceCONST = hyteg::P1ConstantLaplaceOperator;
  using LaplaceCONST_NEW = hyteg::P1ConstantLaplaceOperator_new;
  using LaplaceVAR = hyteg::P1BlendingLaplaceOperator;
  using LaplaceVAR_NEW = hyteg::P1BlendingLaplaceOperator_new;
  using LaplaceLSQP = hyteg::P1PolynomialBlendingLaplaceOperator;
};

template <>
struct FE_Space<ElementType::P2, StencilType::NONE>
{
  using Function = hyteg::P2Function<real_t>;
  using Restriction = hyteg::P2toP2QuadraticRestriction;
  using Prolongation = hyteg::P2toP2QuadraticProlongation;

  using Mass = hyteg::P2ConstantMassOperator; //! todo use new form for blending mass

  using LaplaceCONST = hyteg::P2ConstantLaplaceOperator;
  using LaplaceCONST_NEW = hyteg::P2ConstantLaplaceOperator; //todo new p2 operator not implemented yet
  using LaplaceVAR = hyteg::P2BlendingLaplaceOperator;
  using LaplaceVAR_NEW = hyteg::P2BlendingLaplaceOperator;//todo new p2 operator not implemented yet
  using LaplaceLSQP = hyteg::P2SurrogateLaplaceOperator;
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
};

template<ElementType P>
struct FE_Space<P,StencilType::VARIABLE_NEW> : public P_Space<P>, public OperatorHandler<StencilType::VARIABLE_NEW>
{
  using typename P_Space<P>::Function;
  using typename P_Space<P>::Restriction;
  using typename P_Space<P>::Prolongation;

  using typename P_Space<P>::Mass;
  using Laplace = typename P_Space<P>::LaplaceVAR_NEW;
};

template<ElementType P>
struct FE_Space<P,StencilType::LSQP> : public P_Space<P>, public OperatorHandler<StencilType::LSQP>
{
  using typename P_Space<P>::Function;
  using typename P_Space<P>::Restriction;
  using typename P_Space<P>::Prolongation;

  using typename P_Space<P>::Mass;
  using Laplace = typename P_Space<P>::LaplaceLSQP;

  static void setInterpolationLevel(uint_t level)
  {
    interpolationLevel = level;
  };
};

template<ElementType P>
struct FE_Space<P,StencilType::CONST> : public P_Space<P>, public OperatorHandler<StencilType::CONST>
{
  using typename P_Space<P>::Function;
  using typename P_Space<P>::Restriction;
  using typename P_Space<P>::Prolongation;

  using typename P_Space<P>::Mass;
  using Laplace = typename P_Space<P>::LaplaceCONST;
};

template<ElementType P>
struct FE_Space<P,StencilType::CONST_NEW> : public P_Space<P>, public OperatorHandler<StencilType::CONST_NEW>
{
  using typename P_Space<P>::Function;
  using typename P_Space<P>::Restriction;
  using typename P_Space<P>::Prolongation;

  using typename P_Space<P>::Mass;
  using Laplace = typename P_Space<P>::LaplaceCONST_NEW;
};



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
  auto L = FE::template make_shared<typename FE::Laplace>(storage, minLevel, maxLevel);
  // auto L_compare = FE_Space<P, VARIABLE>::template make_shared<typename FE_Space<P, VARIABLE>::Laplace>(storage, minLevel, maxLevel);

  // exact solution
  u_exact.interpolate(exact, maxLevel);

  // boundary
  u.interpolate(boundary, maxLevel, hyteg::DirichletBoundary);

    // define solver
  auto coarseGridSolver = std::make_shared<hyteg::CGSolver<typename FE::Laplace>>(storage, minLevel, minLevel, max_cg_iter, coarse_tolerance);
  auto restrictionOperator = std::make_shared<typename FE::Restriction>();
  auto prolongationOperator = std::make_shared<typename FE::Prolongation>();
  auto smoother = std::make_shared<hyteg::GaussSeidelSmoother<typename FE::Laplace>>();
  GeometricMultigridSolver<typename FE::Laplace> GMGSolver(storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3);

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
    real_t vCycleTime = 0, solveTime = 0;
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
      // L_compare->apply(u, Lu, maxLevel, hyteg::Inner, Replace);
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
      // for (uint_t i = 0; i < max_cg_iter; ++i)
        // smoother->solve(*L, u, f, maxLevel);
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
      std::string name = "PolynomialBlending_P" + std::to_string(P);
      switch (T)
      {
        case StencilType::CONST:
          name += "const";
          break;

        case StencilType::CONST_NEW:
          name += "const_NEW";
          break;

        case StencilType::VARIABLE:
          name += "var";
          break;

        case StencilType::VARIABLE_NEW:
          name += "var_NEW";
          break;

        case StencilType::LSQP:
          name += "poly" + std::to_string(polyDegree);
          break;
      }

      hyteg::VTKOutput vtkOutput("output", name, storage);
      vtkOutput.add(u);
      vtkOutput.add(err);
      vtkOutput.add(r);
      vtkOutput.add(u_exact);
      vtkOutput.add(f);
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
    case CONST:
      WALBERLA_LOG_INFO_ON_ROOT("Operatortype: Constant Stencil");
      solveTmpl<P, CONST>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
      break;

    case CONST_NEW:
      WALBERLA_LOG_INFO_ON_ROOT("Operatortype: NEW Constant Stencil");
      solveTmpl<P, CONST_NEW>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
      break;

    case VARIABLE:
      WALBERLA_LOG_INFO_ON_ROOT("Operatortype: Variable Stencil");
      solveTmpl<P, VARIABLE>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
      break;

    case VARIABLE_NEW:
      WALBERLA_LOG_INFO_ON_ROOT("Operatortype: NEW Variable Stencil");
      solveTmpl<P, VARIABLE_NEW>(storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs);
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


void showGrid(std::shared_ptr<PrimitiveStorage> storage, const uint_t interpolationlevel)
{
  P1Function<real_t> grid("grid", storage, interpolationlevel, interpolationlevel);

  size_t rowsize = levelinfo::num_microvertices_per_edge( interpolationlevel );
  auto edgeId = grid.getEdgeDataID();
  auto vtxId = grid.getVertexDataID();

  for ( auto& it : storage->getEdges() )
  {
      Edge& edge = *it.second;
      auto data = edge.getData( edgeId )->getPointer( interpolationlevel );

      for( size_t i = 1; i < rowsize - 1; ++i )
      {
        data[vertexdof::macroedge::indexFromVertex( interpolationlevel, i, stencilDirection::VERTEX_C )] = 1;
      }
  }

  for ( auto& it : storage->getVertices() )
  {
      Vertex& vertex = *it.second;
      auto data = vertex.getData( vtxId )->getPointer( interpolationlevel );
      data[0] = 2;
  }

  std::string name = "Annulus_macro";
  hyteg::VTKOutput vtkOutput("../output", name, storage);
  vtkOutput.add(grid);
  vtkOutput.write(interpolationlevel, 0);
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
    // cfg_->readParameterFile("../data/param/PolynomialBlending.prm");
    cfg_->readParameterFile("../hyteg/data/param/PolynomialBlending.prm");
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

  const bool blending = parameters.getParameter<bool>("blending");
  const bool annulus = parameters.getParameter<bool>("annulus");
  uint_t nX = parameters.getParameter<uint_t>("nX");
  uint_t nY = parameters.getParameter<uint_t>("nY");
  uint_t nZ = parameters.getParameter<uint_t>("nZ");
  const uint_t macroLevel = parameters.getParameter<uint_t>("macro_refinement");
  nX = nX << macroLevel;
  nY = nY << macroLevel;
  nZ = nZ << macroLevel;

  const uint_t minPolyDegree = parameters.getParameter<uint_t>("minPolyDegree");
  const uint_t maxPolyDegree = parameters.getParameter<uint_t>("maxPolyDegree");
  const uint_t maxInterpolationLevel = parameters.getParameter<uint_t>("interpolationLevel");
  const uint_t interpolationLevel = std::min(maxLevel, maxInterpolationLevel);

  const uint_t max_outer_iter =  parameters.getParameter<uint_t>("max_outer_iter");
  const uint_t max_cg_iter =  parameters.getParameter<uint_t>("max_cg_iter");
  const real_t mg_tolerance = parameters.getParameter<real_t>("mg_tolerance");
  const real_t coarse_tolerance = parameters.getParameter<real_t>("coarse_tolerance");

  const bool vtk = parameters.getParameter<bool>("vtkOutput");
  const bool show_macrogrid = parameters.getParameter<bool>("show_macrogrid");


  if (annulus)
  {
    if (nZ > 0)
    {
      // todo
      WALBERLA_ABORT("spherical shell not supported yet!");
    }

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
    if (nZ == 0)
    {
      WALBERLA_LOG_INFO_ON_ROOT("Geometry: Rectangle");
    }
    else
    {
      WALBERLA_LOG_INFO_ON_ROOT("Geometry: Cube");
    }
  }


  // define functions and domain

  /// case rectangle
  c_function exact = [](const hyteg::Point3D & x) {return sin(pi*x[0])*sin(pi* x[1]);};
  c_function boundary = [](const hyteg::Point3D &) {return 0;};
  c_function rhs = [](const hyteg::Point3D & x) {return 2*pi*pi*sin(pi*x[0])*sin(pi* x[1]);};

  MeshInfo meshInfo = MeshInfo::meshRectangle(Point2D({0.0, 0.0}), Point2D({1.0, 1.0}), MeshInfo::CRISS, nX, nY);

  if (nZ > 0)
  {
    exact = [](const hyteg::Point3D& x) { return sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]); };
    rhs = [](const hyteg::Point3D& x) { return 3*pi*pi*sin(pi*x[0])*sin(pi*x[1])*sin(pi*x[2]); };

    meshInfo = MeshInfo::meshCuboid(Point3D({0.0,0.0,0.0}), Point3D({1.0,1.0,1.0}), nX, nY, nZ);
  }

  /// case annulus
  Point3D circleCenter{{0, 0, 0}};
  real_t rMin = pi;
  real_t rMax = 2.0*pi;
  real_t middle = (rMax + rMin) / 2.0;
  std::function<real_t(const real_t, const real_t)> exact_polar = [](const real_t r, const real_t phi) { return sin(2*r) * sin(4*phi); };
  if (annulus)
  {

    exact = [&](const hyteg::Point3D& x) { return exact_polar(radius(x), angle(x)); };
    boundary = [&](const hyteg::Point3D & x) {
        real_t r = (radius(x) < middle) ? rMin : rMax;
        return exact_polar(r, angle(x));
    };
    rhs = [&](const hyteg::Point3D& x) {
      real_t r = radius(x);
      real_t T1 = 8 * x[0]*x[1] * (x[0]*x[0] - x[1]*x[1]);
      real_t T2 = pow(r,8)*(-r*cos(2*r) + (2*r*r + 8)*sin(2*r));
      real_t T3 = pow(r,14);
      return T1*T2/T3;
    };

    meshInfo = MeshInfo::meshAnnulus(rMin, rMax, MeshInfo::CRISS, nX, nY);
  }

  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);

  if (annulus && blending)
    AnnulusMap::setMap(setupStorage);

  hyteg::loadbalancing::roundRobin(setupStorage);
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  // WALBERLA_LOG_INFO_ON_ROOT(storage->getGlobalInfo());
  WALBERLA_LOG_INFO_ON_ROOT("Refinement levels: " << minLevel << "->" << maxLevel);

  // choose FE space
  switch (elType)
  {
    case ElementType::P1:
      WALBERLA_LOG_INFO_ON_ROOT("Element Type: P1");
      solve<ElementType::P1>(opType, interpolationLevel, storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs, minPolyDegree, maxPolyDegree);
      break;

    case ElementType::P2:
      WALBERLA_LOG_INFO_ON_ROOT("Element Type: P2");
      solve<ElementType::P2>(opType, interpolationLevel, storage, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk, exact, boundary, rhs, minPolyDegree, maxPolyDegree);
      break;

    default:
      WALBERLA_ABORT("The desired FE space is not supported!");
  }

  if (show_macrogrid)
    showGrid(storage, maxInterpolationLevel);

  return 0;
}
