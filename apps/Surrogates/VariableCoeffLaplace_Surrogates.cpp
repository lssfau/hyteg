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

#include <core/Environment.h>
#include <core/config/Create.h>
#include <core/timing/Timer.h>

#include "core/Format.hpp"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/CircularMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2SurrogateOperator.hpp"
#include "hyteg/p2functionspace/P2VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"
#include "constantStencilOperator/P1ConstantOperator.hpp"

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
  static double setup(std::shared_ptr<OP>, const uint_t) {return 0;}
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

template <>
template < class OP >
double OperatorHandler< StencilType::LSQP >::setup( std::shared_ptr< OP > op, const uint_t polyDegree )
{
  auto start = walberla::timing::getWcTime();
  op->interpolateStencils( polyDegree );
  auto end = walberla::timing::getWcTime();
  op->useDegree( polyDegree );
  return ( end - start );
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
  using LaplaceLSQP = hyteg::P2SurrogateLaplaceOperator;

  using divKgradVAR = hyteg::P2divKgradOperator;
  using divKgradLSQP = hyteg::P2SurrogateDivKgradOperator;
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
    real_t       res        = 0, res_old, discr_l2_err;
    double       vCycleTime = 0, solveTime = 0;
    real_t       averageConvergenceRate = 0;
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

      hyteg::VTKOutput vtkOutput("../../output", name, storage);
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

// write diagonal entries of stiffness matrix to vtk-file
void showStencilFunction(std::shared_ptr<PrimitiveStorage> storage, const uint_t level, const real_t alpha)
{
  typedef stencilDirection SD;

  P2Function<real_t> stencil("stencil", storage, level, level);

  stencil.add( real_c( nan( "" ) ), level, All );

  auto vtxID  = stencil.getVertexDoFFunction().getFaceDataID();
  auto edgeID = stencil.getEdgeDoFFunction().getFaceDataID();

  for (auto& faceit : storage->getFaces())
  {
    hyteg::Face& face = *faceit.second;

    real_t* vtxData = face.getData(vtxID)->getPointer(level);
    real_t* edgeData = face.getData(edgeID)->getPointer(level);

    P2Form_divKgrad form;
    form.setGeometryMap(face.getGeometryMap());

    Point3D x0( face.getCoordinates()[0] ), x;
    real_t  h = real_c( 1.0 ) / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );

    Point3D d0 = h * (face.getCoordinates()[1] - face.getCoordinates()[0]);
    Point3D d2 = h * (face.getCoordinates()[2] - face.getCoordinates()[0]);

    // directions
    const Point3D dirS  = -d2;
    const Point3D dirSE = d0 - d2;
    const Point3D dirE  = d0;
    const Point3D dirW  = -dirE;
    const Point3D dirNW = -dirSE;
    const Point3D dirN  = -dirS;
    const Point3D dirNE = dirN + dirE;

    // stencil entries
    std::array<real_t, P2::NumStencilentries2D::VtV> VtVStencil;
    std::array<real_t, P2::NumStencilentries2D::EtV> EtVStencil;
    std::array<real_t, P2::NumStencilentries2D::VtE> VtEStencil;
    std::array<real_t, P2::NumStencilentries2D::EtE> EtEStencil;

    // loop over all DOFs
    for (const auto& it : hyteg::edgedof::macroface::Iterator(level, 0))
    {
      x = x0 + walberla::real_c(it.y()) * d2 + walberla::real_c(it.x()) * d0;

      P2::variablestencil::macroface::assembleStencil(form, x, dirS, dirSE, dirE, dirN, dirNW, dirW, dirNE,
                                                      VtVStencil, EtVStencil, VtEStencil, EtEStencil);
      idx_t i = it.x();
      idx_t j = it.y();

      // VERTEX DoF
      if (!vertexdof::macroface::isVertexOnBoundary(level, it))
      {
        vtxData[vertexdof::macroface::indexFromVertex(level, i, j, SD::VERTEX_C)] =
          VtVStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)];
      }

      // HORIZONTAL EDGE DoF
      if (!edgedof::macroface::isHorizontalEdgeOnBoundary(level, it))
      {
        edgeData[edgedof::macroface::indexFromHorizontalEdge(level, i, j, SD::EDGE_HO_C)] =
          EtEStencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_HO_C)];
      }

      // VERTICAL EDGE DoF
      if (!edgedof::macroface::isVerticalEdgeOnBoundary(level, it))
      {
        edgeData[edgedof::macroface::indexFromVerticalEdge(level, i, j, SD::EDGE_VE_C)] =
          EtEStencil[edgedof::stencilIndexFromVerticalEdge(SD::EDGE_VE_C)];
      }

      // DIAGONAL EDGE DoF
      if (!edgedof::macroface::isDiagonalEdgeOnBoundary(level, it))
      {
        edgeData[edgedof::macroface::indexFromDiagonalEdge(level, i, j, SD::EDGE_DI_C)] =
          EtEStencil[edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_DI_C)];
      }
    }
  }

  std::string name = "divKgrad_stencil_alpha=" + std::to_string(int(alpha));
  hyteg::VTKOutput vtkOutput("../../output", name, storage);
  vtkOutput.add(stencil);
  vtkOutput.write(level, 0);
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
    cfg_->readParameterFile("../../data/param/VariableCoeffLaplace_Surrogates.prm");
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
  const bool show_stencil = parameters.getParameter<bool>("show_stencil");

  const real_t alpha = parameters.getParameter<real_t>("alpha");
  const real_t phi = parameters.getParameter<real_t>("phi");

  const uint_t n_el = parameters.getParameter<uint_t>("n_el");

  // define functions and domain
  // case: smooth k, domain = triangle
  // c_function exact = [](const hyteg::Point3D & x) {return sin(x[0])*sinh(x[1]);};
  // c_function boundary = exact;
  // c_function rhs = [](const hyteg::Point3D & x) {return - cos(x[0])*cos(x[0]) * sinh(x[1]);};
  // c_function k = [](const hyteg::Point3D & x) {return 2.0 +  sin(x[0]);};
  c_function exact = [phi](const hyteg::Point3D& x) { return sin(phi*PI*x[0])*sinh(PI*x[1]); };
  c_function boundary = exact;
  c_function k = [alpha](const hyteg::Point3D& x) { return tanh(alpha*(x[0] - 0.5)) + 2; };
  c_function rhs = [phi,alpha](const hyteg::Point3D& x) {
     real_t t0 = tanh( alpha * ( x[0] - real_c( 0.5 ) ) );
     real_t t1 = phi * alpha * ( t0 * t0 - real_c( 1 ) ) * cos( phi * PI * x[0] );
     real_t t2 = ( phi * phi - real_c( 1 ) ) * PI * ( t0 + real_c( 2 ) ) * sin( phi * PI * x[0] );
     return PI * ( t1 + t2 ) * sinh( PI * x[1] );
  };

  std::string msh = "../../data/meshes/tri_" + std::to_string(n_el) + "el.msh";
  MeshInfo meshInfo = MeshInfo::fromGmshFile(msh);
  P2Form_divKgrad::callback = k;

  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  setupStorage.setMeshBoundaryFlagsOnBoundary(1, 0, true);

  hyteg::loadbalancing::roundRobin(setupStorage);
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);

  WALBERLA_LOG_INFO_ON_ROOT("Refinement levels: " << minLevel << "->" << maxLevel);
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

  if (show_stencil)
    showStencilFunction(storage, maxInterpolationLevel, alpha);

  return 0;
}
