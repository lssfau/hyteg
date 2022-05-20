/*
 * Copyright (c) 2017-2021 Nils Kohl.
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

/**
 * \page 07_IsoviscousConvectionAnnulus A full app simulating isoviscous convection on an annulus.
 *
 * \dontinclude tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp
 *
 * \brief In this tutorial we will set up a complete app that solves a coupled system
 * of the Stokes equations and the advection-diffusion equation.
 *
 * \section IsoviscousConvectionAnnulus-equations Governing equations
 *
 * We are in this tutorial considering a 2D simulation of isoviscous convection on an annular domain.
 * The convection is modeled using the Boussinesq-approximation, by coupling the Stokes equations to an advection-diffusion
 * equation. The system is given as:
 *
 * \f{align*}{
 *  -\Delta u + \nabla p &= \text{Ra}\, c\, g \\
 *  \nabla u &= 0 \\
 *  \frac{\partial}{\partial t}c + u \cdot \nabla c - \nabla \cdot \kappa \nabla c &= q
 * \f}
 *
 * where \f$u\f$ is the velocity of the fluid, \f$p\f$ the pressure, \f$c\f$ the temperature,
 * \f$q\f$ internal heating, \f$\text{Ra}\f$ the dimensionless Rayleigh number, \f$\kappa\f$ a diffusivity constant
 * and \f$g\f$ the gravitational force.
 *
 * For simplicity, we set \f$\kappa = 1\f$  and steer the vigor of convection through the Rayleigh number.
 * Also we assume no internal heating.
 *
 *
 * \section IsoviscousConvectionAnnulus-domain Domain
 *
 * The annulus domain is created in four steps.
 * 1. We define the unstructured mesh through a MeshInfo object. There are some pre-defined mesh generators, for example for the
 *    annulus domain. For external unstructured tetrahedral or triangular meshes, a gmsh reader is implemented.
 * 2. The HyTeG primitive data structures are created from the mesh via a SetupPrimitiveStorage object.
 * 3. A geometry mapping is applied to map the primitives to the physical domain.
 * 4. The primitives are distributed among the parallel processes via construction of the final PrimitiveStorage object.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Domain setup
 *
 * \section IsoviscousConvectionAnnulus-functionsandoperators Function spaces and operators
 *
 * We now define the function spaces and operators that we need in our app.
 *
 * In HyTeG, there are multiple implementations of the same linear operator.
 * This way we exploit properties of the operator to achieve maximum performance of the underlying compute kernels.
 * For example, operators that arise from bilinear forms with constant coefficients are implemented by patchwise constant
 * stencils.
 *
 * A certain set of operators is already implemented and the user does not have to care about assembly of corresponding
 * stencils and/or (local) stiffness matrices.
 *
 * Functions (== elements of a finite element space) are represented by corresponding types in HyTeG.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Function and operator typedefs
 *
 * Since we focus multigrid solvers, functions and operators are constructed over a hierarchy of levels.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Function setup
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Operator setup
 *
 * \section IsoviscousConvectionAnnulus-initialcondition Initial and boundary conditions
 *
 * Analytical functions can be interpolated (think 'sampled') into the finite element functions by a call to interpolate.
 * Here we define an initial (and boundary) temperature function for the mantle.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Initial temperature
 *
 * \section IsoviscousConvectionAnnulus-solvers Solvers
 *
 * We solve the non-linear system by splitting it into three components that are solved in an alternating fashion.
 *
 * First we construct the iterative solver for the Stokes equation.
 *
 * We employ a monolithic multigrid solver with inexact-Uzawa relaxation.
 * The solver is composed by combining different components such as grid transfer operators, smoothers, and coarse
 * grid solvers. To allow for simple composition in C++ we make use std::shared_ptr<>.
 *
 * Most components are derived from a base class Solver, and templated with the linear operator
 * so that e.g. the same GeometricMultigridSolver implementation is reused also for other discretizations or equations.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Stokes solver setup
 *
 * For the unsteady diffusion, we require a linear solver for the operator and use a wrapper that simplifies time stepping.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Diffusion solver setup
 *
 * The advection is treated by a Lagrangian solver that employs tracer particles.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Advection solver setup
 *
 * To visualize the results, we employ VTK. To output functions in VTK format, we prepare a VTKOutput object and add the
 * functions of interest.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp VTK
 *
 * Before solving the Stokes equation, we need to assemble the right-hand side according to the Boussinesq-approximation.
 * Note that the temperature enters the RHS in strong form if we simply set \f$ f := Ra * c * n \f$.
 * So we need to premultiply with the finite element mass matrix first.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp RHS
 *
 * For an initial velocity field != 0 we solve the Stokes equation before starting the actual time stepping.
 * For simplicity, we apply a fixed number of v-cycles. In practical applications we could also steer this
 * via a residual threshold.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Stokes solve initial
 *
 * We compute the variable time-step size via a CFL condition.
 * Therefore we compute the maximum velocity.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Max velocity
 *
 * \section IsoviscousConvectionAnnulus-output Output
 *
 * VTK output is pretty straightforward ...
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp VTK write
 *
 * Simulations are worthless without storing relevant results in some way.
 * This can be done e.g. via a SQLite interface. In parallel runs, it is advisable to only write into the database from
 * one process (e.g. root). If necessary, relevant, distributed data should be reduced first.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp DB
 *
 * \section IsoviscousConvectionAnnulus-timeloop Time loop
 *
 * As discussed, we compute the time-step size via a CFL condition.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp CFL
 *
 * In each time step, we now first call the advection solver.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Advection
 *
 * Then the diffusion solver.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Diffusion
 *
 * Update the right-hand side.
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp RHS update
 *
 * And then solve the Stokes equation (again with a fixed number of iterations for simplicity).
 *
 * \snippet tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp Stokes
 *
 * \section IsoviscousConvectionAnnulus-results Results
 *
 * The resulting temperature field after a simulation over 1000 time-steps is shown below.
 * Plumes are rising from the heated inner boundary outwards to the outer boundary.
 * Once a temperature gradient has developed at the outer boundary plumes of cool material are
 * sinking in the opposite direction. After about 800 time-steps a chaotic behavior sets in.
 *
 * \htmlonly
  <center>
  <table>
  <tr>
  <td><img src="convection_annulus.0000.png" width="100%"/><center>initial condition</center></td>
  <td><img src="convection_annulus.0025.png" width="100%"/><center>after 25 time-steps</center></td>
  <td><img src="convection_annulus.0050.png" width="100%"/><center>after 50 time-steps</center></td>
  </tr>
  <tr>
  <td><img src="convection_annulus.0075.png" width="100%"/><center>after 75 time-steps</center></td>
  <td><img src="convection_annulus.0100.png" width="100%"/><center>after 100 time-steps</center></td>
  <td><img src="convection_annulus.0200.png" width="100%"/><center>after 200 time-steps</center></td>
  </tr>
  <tr>
  <td><img src="convection_annulus.0500.png" width="100%"/><center>after 500 time-steps</center></td>
  <td><img src="convection_annulus.0800.png" width="100%"/><center>after 800 time-steps</center></td>
  <td><img src="convection_annulus.1000.png" width="100%"/><center>after 1000 time-steps</center></td>
  </tr>
  </table>
  </center>
  \endhtmlonly
 *
 *
 * \section IsoviscousConvectionAnnulus-fullApp Full Application
 * \include tutorials/07_IsoviscousConvectionAnnulus/IsoviscousConvectionAnnulus.cpp
*/

#include <cmath>
#include <core/Environment.h>

#include "core/DataTypes.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/SQL.hpp"
#include "hyteg/dataexport/TimingOutput.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/UzawaSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesVelocityBlockBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"

namespace hyteg {

using walberla::int_c;
using walberla::real_t;
using walberla::math::pi;

struct DomainInfo
{
   real_t rMin = 0;
   real_t rMax = 0;
   uint_t nTan = 0;
   uint_t nRad = 0;
};

struct SolverInfo
{
   uint_t stokesMaxNumIterations = 10;
   uint_t uzawaInnerIterations   = 10;
   uint_t uzawaPreSmooth         = 6;
   uint_t uzawaPostSmooth        = 6;
   real_t uzawaOmega             = 0.3;

   uint_t diffusionMaxNumIterations           = 10000;
   real_t diffusionAbsoluteResidualUTolerance = 10000;
};

/// Calculates and returns
///
///     ||u||_L2 = sqrt( u^T M u )
///
template < typename FunctionType, typename MassOperator >
real_t normL2( const FunctionType& u, const FunctionType& tmp, const MassOperator& M, const uint_t& level, const DoFType& flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   return std::sqrt( u.dotGlobal( tmp, level, flag ) );
}

template < typename FunctionType, typename MassOperator >
real_t normL2Squared( const FunctionType& u,
                      const FunctionType& tmp,
                      const MassOperator& M,
                      const uint_t&       level,
                      const DoFType&      flag )
{
   tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   return u.dotGlobal( tmp, level, flag );
}

template < typename StokesFunction, typename StokesOperator, typename MassOperatorVelocity, typename MassOperatorPressure >
void calculateStokesResiduals( const StokesOperator&       A,
                               const MassOperatorVelocity& Mu,
                               const MassOperatorPressure& Mp,
                               const StokesFunction&       x,
                               const StokesFunction&       b,
                               uint_t                      level,
                               StokesFunction&             r,
                               StokesFunction&             tmp,
                               real_t&                     residual )
{
   tmp.interpolate( 0, level, All );
   r.interpolate( 0, level, All );

   A.apply( x, tmp, level, Inner );
   r.assign( { 1.0, -1.0 }, { b, tmp }, level, Inner );
   residual = normL2Squared( r.uvw()[0], tmp.uvw()[0], Mu, level, Inner );
   residual += normL2Squared( r.uvw()[1], tmp.uvw()[1], Mu, level, Inner );
   residual += normL2Squared( r.p(), tmp.p(), Mp, level, Inner );
   residual = std::sqrt( residual );
}

/// [Domain setup]
std::shared_ptr< PrimitiveStorage > createPrimitiveStorage( DomainInfo domainInfo )
{
   // A MeshInfo object represents an unstructured mesh.
   // In HyTeG this mesh usually acts as a coarse approximation to the domain.
   // The elements of this coarse mesh are then refined multiple times in a later step.
   MeshInfo meshInfo =
       MeshInfo::meshAnnulus( domainInfo.rMin, domainInfo.rMax, MeshInfo::CRISS, domainInfo.nTan, domainInfo.nRad );

   // A SetupPrimitiveStorage is constructed to convert the unstructured mesh to a graph
   // of macro-primitives that is used to efficiently parallelize the mesh data structure.
   // Interface primitives are inserted to allow for communication between macro-elements.
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // To ensure that the macro-primitives are evenly distributed, a load balancing algorithm
   // assigns process ranks to the individual macro-primitives.
   loadbalancing::roundRobinVolume( *setupStorage );

   // Boundary conditions (e.g. Dirichlet, Neumann) are assigned in two steps.
   // First a boundary flag (integer) is assigned to each primitive.
   // This helper function assigns 1 to all macro-primitives on the boundary and 0 to all remaining primitives.
   // Then a boundary condition is assigned to each integer. The default assignment is
   //  0 -> inner
   //  1 -> Dirichlet
   //  2 -> Neumann
   // Since we employ Dirichlet conditions in this tutorial, this is all we need to do.
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   // Since we project the computational domain to the annulus we need to apply a mapping.
   // There are some predefined mappings such as the annulus map and a map for the spherical shell.
   AnnulusMap::setMap( *setupStorage );

   // Finally the distributed PrimitiveStorage object is constructed.
   // From now on this object is used to handle the domain information.
   // It only contains information of process-local primitives and their direct neighborhood.
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

   // Note that the MeshInfo and SetupPrimitiveStorage objects are destroyed when we leave the scope of this function.

   return storage;
}
/// [Domain setup]

void runBenchmark( real_t      cflMax,
                   real_t      rayleighNumber,
                   uint_t      maxNumTimeSteps,
                   uint_t      minLevel,
                   uint_t      maxLevel,
                   SolverInfo  solverInfo,
                   DomainInfo  domainInfo,
                   real_t      simulationTime,
                   bool        vtk,
                   std::string outputDirectory,
                   std::string outputBaseName )
{
   auto storage = createPrimitiveStorage( domainInfo );

   if ( vtk )
   {
      writeDomainPartitioningVTK( storage, outputDirectory, outputBaseName + "_domain" );
   }

   const uint_t unknownsTemperature = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );
   const uint_t unknownsStokes      = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );

   const real_t hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );
   const real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   const real_t diffusivity                 = 1.0;
   const real_t internalHeating             = 0.0;
   const real_t initialTemperatureSteepness = 10.0;

   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark name: " << outputBaseName )
   WALBERLA_LOG_INFO_ON_ROOT( " - time discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + max CFL:                                      " << cflMax )
   WALBERLA_LOG_INFO_ON_ROOT( "   + simulation time:                              " << simulationTime )
   WALBERLA_LOG_INFO_ON_ROOT( " - space discretization: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + rMin:                                         " << domainInfo.rMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + rMax:                                         " << domainInfo.rMax )
   WALBERLA_LOG_INFO_ON_ROOT( "   + nTan:                                         " << domainInfo.nTan )
   WALBERLA_LOG_INFO_ON_ROOT( "   + nRad:                                         " << domainInfo.nRad )
   WALBERLA_LOG_INFO_ON_ROOT( "   + dimensions:                                   " << ( storage->hasGlobalCells() ? "3" : "2" ) )
   WALBERLA_LOG_INFO_ON_ROOT( "   + maxLevel:                                        " << maxLevel )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns temperature, including boundary:     " << unknownsTemperature )
   WALBERLA_LOG_INFO_ON_ROOT( "   + unknowns Stokes, including boundary:          " << unknownsStokes )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_min:                                        " << hMin )
   WALBERLA_LOG_INFO_ON_ROOT( "   + h_max:                                        " << hMax )
   WALBERLA_LOG_INFO_ON_ROOT( " - benchmark settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Rayleigh number:                              " << rayleighNumber )
   WALBERLA_LOG_INFO_ON_ROOT( "   + internal heating:                             " << internalHeating )
   WALBERLA_LOG_INFO_ON_ROOT( " - app settings: " )
   WALBERLA_LOG_INFO_ON_ROOT( "   + Uzawa omega:                                  " << solverInfo.uzawaOmega )
   WALBERLA_LOG_INFO_ON_ROOT( "   + output directory:                             " << outputDirectory )
   WALBERLA_LOG_INFO_ON_ROOT( "   + output base name:                             " << outputBaseName )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );
   WALBERLA_LOG_INFO_ON_ROOT( "" );



   /// [Function and operator typedefs]
   typedef P2P1TaylorHoodFunction< real_t > StokesFunction;
   typedef P2Function< real_t >             ScalarFunction;

   typedef P2P1ElementwiseBlendingStokesOperator  StokesOperator;
   typedef P2ElementwiseBlendingLaplaceOperator   LaplaceOperator;
   typedef P2ElementwiseBlendingMassOperator      MassOperatorVelocity;
   typedef P1ConstantMassOperator                 MassOperatorPressure;
   typedef P2ElementwiseUnsteadyDiffusionOperator UnsteadyDiffusionOperator;
   /// [Function and operator typedefs]

   /// [Function setup]
   ScalarFunction c( "c", storage, minLevel, maxLevel );
   /// [Function setup]

   ScalarFunction cOld( "cOld", storage, minLevel, maxLevel );
   ScalarFunction cTmp( "cTmp", storage, minLevel, maxLevel );
   ScalarFunction cTmp2( "cTmp2", storage, minLevel, maxLevel );
   ScalarFunction q( "q", storage, minLevel, maxLevel );

   StokesFunction u( "u", storage, minLevel, maxLevel );
   StokesFunction uLast( "uLast", storage, minLevel, maxLevel );
   StokesFunction f( "f", storage, minLevel, maxLevel );
   StokesFunction outwardNormal( "outwardNormal", storage, minLevel, maxLevel );
   StokesFunction stokesTmp( "stokesTmp", storage, minLevel, maxLevel );
   StokesFunction stokesResidual( "stokesResidual", storage, minLevel, maxLevel );

   ScalarFunction uMagnitudeSquared( "uMagnitudeSquared", storage, minLevel, maxLevel );

   ScalarFunction uTmp( "uTmp", storage, minLevel, maxLevel );
   ScalarFunction uTmp2( "uTmp2", storage, minLevel, maxLevel );

   /// [Operator setup]
   auto                 A = std::make_shared< StokesOperator >( storage, minLevel, maxLevel );
   LaplaceOperator      L( storage, minLevel, maxLevel );

   MassOperatorVelocity MVelocity( storage, minLevel, maxLevel );
   MassOperatorPressure MPressure( storage, minLevel, maxLevel );

   UnsteadyDiffusionOperator diffusionOperator(
       storage, minLevel, maxLevel, 1.0, diffusivity, DiffusionTimeIntegrator::ImplicitEuler );
   /// [Operator setup]

   /// [Initial temperature]
   std::function< real_t( const Point3D& ) > initialTemperature = [&]( const Point3D& x ) {
      auto radius = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      return std::exp( -initialTemperatureSteepness * ( ( radius - domainInfo.rMin ) / ( domainInfo.rMax - domainInfo.rMin ) ) );
   };

   c.interpolate( initialTemperature, maxLevel, All );
   /// [Initial temperature]

   std::function< real_t( const Point3D& ) > normalX = []( const Point3D& x ) { return x[0] / x.norm(); };
   std::function< real_t( const Point3D& ) > normalY = []( const Point3D& x ) { return x[1] / x.norm(); };

   outwardNormal.uvw().interpolate( {normalX, normalY}, maxLevel );

   A->computeAndStoreLocalElementMatrices();
   L.computeAndStoreLocalElementMatrices();
   MVelocity.computeAndStoreLocalElementMatrices();

   /// [Stokes solver setup]

   // Solvers in HyTeG are constructed by plugging together individual solver components.
   // Most components are derived from the Solver class.
   std::shared_ptr< Solver< StokesOperator > > stokesSolver;

   // For the inexact Uzawa relaxation, we need a smoother for the A11-block of the discretized Stokes system.
   // We construct a scalar, weighted Jacobi smoother. Note that this component is derived from the Solver class.
   std::shared_ptr< Solver< typename StokesOperator::VelocityOperator_T > > smoother =
       std::make_shared< WeightedJacobiSmoother< typename StokesOperator::VelocityOperator_T > >(
           storage, minLevel, maxLevel, 0.66 );

   // The relaxation on the viscous block is now built using the scalar velocity smoother.
   // It is applied to both components (x, and y).
   auto uzawaVelocityPreconditioner =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, smoother );

   // We plug the relaxation into the full Uzawa relaxation instance.
   auto uzawaSmoother = std::make_shared< UzawaSmoother< StokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, solverInfo.uzawaOmega, Inner, solverInfo.uzawaInnerIterations );

   // To solve on the coarse grid, we use a standard MinRes solver.
   // Some predefined solvers are available under the solvertemplates namespace.
   auto coarseGridSolver = solvertemplates::stokesMinResSolver< StokesOperator >( storage, minLevel, 1e-8, 1000 );

   // The monolithic MG solver requires grid transfer of both velocity and pressure.
   auto prolongationOperator = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto restrictionOperator  = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );

   // Finally, all components are plugged together to build our solver for the Stokes equation.
   auto multigridSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >( storage,
                                                                                          uzawaSmoother,
                                                                                          coarseGridSolver,
                                                                                          restrictionOperator,
                                                                                          prolongationOperator,
                                                                                          minLevel,
                                                                                          maxLevel,
                                                                                          solverInfo.uzawaPreSmooth,
                                                                                          solverInfo.uzawaPostSmooth,
                                                                                          2,
                                                                                          CycleType::VCYCLE );

   stokesSolver = multigridSolver;
   /// [Stokes solver setup]

   /// [Diffusion solver setup]

   // The diffusive term is solved using a standard CG solver.
   auto diffusionLinearSolver = std::make_shared< CGSolver< UnsteadyDiffusionOperator > >(
       storage, minLevel, maxLevel, solverInfo.diffusionMaxNumIterations, solverInfo.diffusionAbsoluteResidualUTolerance );

   // To simplify the time-stepping procedure, a wrapper is available.
   // We split the construction into two parts to be able easily to exchange the linear solver.
   UnsteadyDiffusion< ScalarFunction, UnsteadyDiffusionOperator, LaplaceOperator, MassOperatorVelocity > diffusionSolver(
       storage, minLevel, maxLevel, diffusionLinearSolver );
   /// [Diffusion solver setup]

   /// [Advection solver setup]
   // For the advection term we employ a characteristic Galerkin solver that employs tracer particles to transport
   // a scalar quantity (temperature in this case). MMOC stands for modified method of characteristics.
   MMOCTransport< ScalarFunction > transport( storage, minLevel, maxLevel, TimeSteppingScheme::RK4 );
   /// [Advection solver setup]

   real_t timeTotal = 0;
   real_t vMax      = 0;
   real_t vRms      = 0;
   real_t residual  = 0;

   real_t timeStepTotal = 0;
   real_t timeStokes    = 0;
   real_t timeMMOC      = 0;
   real_t timeDiffusion = 0;
   real_t timeVTK       = 0;

   /// [VTK]
   VTKOutput vtkOutput( outputDirectory, outputBaseName, storage );
   vtkOutput.add( c );
   vtkOutput.add( u );
   vtkOutput.add( f );
   /// [VTK]

   uint_t            timeStep = 0;
   walberla::WcTimer timeStepTimer;
   walberla::WcTimer localTimer;

   timeStepTimer.start();

   /// [RHS]
   MVelocity.apply( c, f.uvw()[0], maxLevel, All );
   MVelocity.apply( c, f.uvw()[1], maxLevel, All );

   f.uvw().multElementwise( { f.uvw(), outwardNormal.uvw() }, maxLevel );
   f.uvw().assign( { rayleighNumber }, { f.uvw() }, maxLevel, All );
   /// [RHS]

   calculateStokesResiduals( *A, MVelocity, MPressure, u, f, maxLevel, stokesResidual, stokesTmp, residual );

   localTimer.start();

   /// [Stokes solve initial]
   for ( uint_t i = 0; i < solverInfo.stokesMaxNumIterations; i++ )
   {
      stokesSolver->solve( *A, u, f, maxLevel );
   }
   /// [Stokes solve initial]

   localTimer.end();
   timeStokes = localTimer.last();

   uLast.assign( { 1.0 }, { u }, maxLevel, All );

   calculateStokesResiduals( *A, MVelocity, MPressure, u, f, maxLevel, stokesResidual, stokesTmp, residual );

   /// [Max velocity]
   vMax = velocityMaxMagnitude( u.uvw(), uTmp, uMagnitudeSquared, maxLevel, All );
   /// [Max velocity]

   localTimer.start();
   if ( vtk )
   {
      /// [VTK write]
      vtkOutput.write( maxLevel );
      /// [VTK write]
   }
   localTimer.end();
   timeVTK = localTimer.last();

   timeStepTimer.end();
   timeStepTotal = timeStepTimer.last();

   WALBERLA_LOG_INFO_ON_ROOT(
       " timestep |           dt |   time total | velocity RMS | velocity max magnitude |  res. Stokes |  total | Stokes |   MMOC |   diff |    VTK |" )
   WALBERLA_LOG_INFO_ON_ROOT(
       "----------+--------------+--------------+--------------+------------------------+--------------+--------+--------+--------+--------+--------+" )
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( " %8s | %12s | %12.8f | %12.4f | %22.4f | %12.5e | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f |",
                         "initial",
                         "-",
                         timeTotal,
                         vRms,
                         vMax,
                         residual,
                         timeStepTotal,
                         timeStokes,
                         timeMMOC,
                         timeDiffusion,
                         timeVTK ) )

   /// [DB]
   FixedSizeSQLDB db( outputDirectory + "/" + outputBaseName + ".db" );

   db.setVariableEntry( "ts", uint_c( timeStep ) );
   db.setVariableEntry( "sim_time", timeTotal );
   db.setVariableEntry( "v_max", vMax );
   db.setVariableEntry( "dt", real_c( 0 ) );
   db.writeRowOnRoot();
   /// [DB]

   while ( timeTotal < simulationTime && timeStep < maxNumTimeSteps )
   {
      timeStepTimer.start();

      timeStep++;

      // new time step size

      /// [CFL]
      vMax      = velocityMaxMagnitude( u.uvw(), uTmp, uTmp2, maxLevel, All );
      real_t dt = ( cflMax / vMax ) * hMin;
      /// [CFL]

      // energy

      // advection

      localTimer.start();

      /// [Advection]
      transport.step( c, u.uvw(), u.uvw(), maxLevel, All, dt, 1, true );
      /// [Advection]
      localTimer.end();
      timeMMOC = localTimer.last();

      // diffusion

      diffusionOperator.setDt( dt );
      diffusionOperator.getOperator().computeAndStoreLocalElementMatrices();

      cOld.assign( { 1.0 }, { c }, maxLevel, All );

      localTimer.start();

      /// [Diffusion]
      diffusionSolver.step( diffusionOperator, L, MVelocity, c, cOld, q, q, maxLevel, Inner );
      /// [Diffusion]
      localTimer.end();
      timeDiffusion += localTimer.last();

      // Stokes

      /// [RHS update]
      MVelocity.apply( c, f.uvw()[0], maxLevel, All );
      MVelocity.apply( c, f.uvw()[1], maxLevel, All );

      f.uvw().multElementwise( { f.uvw(), outwardNormal.uvw() }, maxLevel );
      f.uvw().assign( { rayleighNumber }, { f.uvw() }, maxLevel, All );
      /// [RHS update]

      calculateStokesResiduals( *A, MVelocity, MPressure, u, f, maxLevel, stokesResidual, stokesTmp, residual );

      localTimer.start();
      /// [Stokes]
      for ( uint_t i = 0; i < solverInfo.stokesMaxNumIterations; i++ )
      {
         stokesSolver->solve( *A, u, f, maxLevel );
      }
      /// [Stokes]
      localTimer.end();
      timeStokes = localTimer.last();

      calculateStokesResiduals( *A, MVelocity, MPressure, u, f, maxLevel, stokesResidual, stokesTmp, residual );

      timeTotal += dt;

      vMax = velocityMaxMagnitude( u.uvw(), uTmp, uMagnitudeSquared, maxLevel, All );

      localTimer.start();
      if ( vtk )
      {
         vtkOutput.write( maxLevel, timeStep );
      }
      localTimer.end();
      timeVTK = localTimer.last();

      timeStepTimer.end();
      timeStepTotal = timeStepTimer.last();

      db.setVariableEntry( "ts", uint_c( timeStep ) );
      db.setVariableEntry( "sim_time", timeTotal );
      db.setVariableEntry( "v_max", vMax );
      db.setVariableEntry( "dt", dt );
      db.writeRowOnRoot();

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( " %8d | %12.5e | %12.8f | %12.4f | %22.4f | %12.5e | %6.2f | %6.2f | %6.2f | %6.2f | %6.2f |",
                            timeStep,
                            dt,
                            timeTotal,
                            vRms,
                            vMax,
                            residual,
                            timeStepTotal,
                            timeStokes,
                            timeMMOC,
                            timeDiffusion,
                            timeVTK ) )
   }
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./parameters.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   hyteg::DomainInfo domainInfo;
   hyteg::SolverInfo solverInfo;

   domainInfo.rMin = mainConf.getParameter< hyteg::real_t >( "rMin" );
   domainInfo.rMax = mainConf.getParameter< hyteg::real_t >( "rMax" );
   domainInfo.nTan = mainConf.getParameter< hyteg::uint_t >( "nTan" );
   domainInfo.nRad = mainConf.getParameter< hyteg::uint_t >( "nRad" );

   const hyteg::real_t cflMax          = mainConf.getParameter< hyteg::real_t >( "cflMax" );
   const hyteg::real_t rayleighNumber  = mainConf.getParameter< hyteg::real_t >( "rayleighNumber" );
   const uint_t        minLevel        = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t        level           = mainConf.getParameter< uint_t >( "level" );
   const hyteg::real_t simulationTime  = mainConf.getParameter< hyteg::real_t >( "simulationTime" );
   const hyteg::uint_t maxNumTimeSteps = mainConf.getParameter< hyteg::uint_t >( "maxNumTimeSteps" );

   solverInfo.stokesMaxNumIterations = mainConf.getParameter< uint_t >( "stokesMaxNumIterations" );
   solverInfo.uzawaOmega             = mainConf.getParameter< real_t >( "uzawaOmega" );
   solverInfo.uzawaInnerIterations   = mainConf.getParameter< uint_t >( "uzawaInnerIterations" );
   solverInfo.uzawaPreSmooth         = mainConf.getParameter< uint_t >( "uzawaPreSmooth" );
   solverInfo.uzawaPostSmooth        = mainConf.getParameter< uint_t >( "uzawaPostSmooth" );

   solverInfo.diffusionMaxNumIterations           = mainConf.getParameter< uint_t >( "diffusionMaxNumIterations" );
   solverInfo.diffusionAbsoluteResidualUTolerance = mainConf.getParameter< real_t >( "diffusionAbsoluteResidualUTolerance" );

   const std::string outputDirectory = mainConf.getParameter< std::string >( "outputDirectory" );
   const std::string outputBaseName  = mainConf.getParameter< std::string >( "outputBaseName" );
   const bool        vtk             = mainConf.getParameter< bool >( "vtk" );

   hyteg::runBenchmark( cflMax,
                        rayleighNumber,
                        maxNumTimeSteps,
                        minLevel,
                        level,
                        solverInfo,
                        domainInfo,
                        simulationTime,
                        vtk,
                        outputDirectory,
                        outputBaseName );
}
