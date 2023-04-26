/*
 * Copyright (c) 2023 Ponsuganth Ilangovan P,.
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
 * \page 11_CouetteFlow Tutorial to implement simple Couette flow in HyTeG.
 *
 * \dontinclude tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp
 *
 * \brief This tutorial demonstrates the implementation of a Couette flow problem on a circular annulus with different tangential velocities on the inner and outer boundaries of the annulus. On polar coordinates the \f$ u_r = 0\f$ everywhere. After implementation we also perform a convergence analysis to show that FEM solution from HyTeG converges to the analytical solution of the problem with expected order of convergence.
 *
 * \section task Task
 * 
 * <img src="11_CouetteFlow_Outline.png" width="30%" />
 * 
 * The problem that needs to be solved is the Stokes equation on the annulus which is given as below
 * \f[
   \begin{align*}
      -\Delta\mathbf{u}_c + \nabla p &= 0\\[2ex]
      \nabla\cdot\mathbf{u}_c &= 0
   \end{align*}
 * \f]
 * Here \f$ \mathbf{u}_c \f$ denotes the equations are written in Cartesian coordinates which we will denote henceforth with just \f$ \mathbf{u} \f$, however for this problem, it is easy to specify the boundary conditions in polar coordinates as \f$\mathbf{u}_\phi(r = r_{min}) = u_{inner},\quad \mathbf{u}_\phi(r = r_{max}) = u_{outer}\f$
 * 
 * The finite element weak form of the boundary value problem is to find \f$ \mathbf{u}_h \in \mathbf{V}_h \f$ and \f$ p_h \in Q_h \f$ such that,
 * \f[
 * \begin{align*}
 *    a(\mathbf{u}_h, \mathbf{v}_h) + b(p_h, \mathbf{v}_h) &= \mathbf{0}\quad\quad\forall\mathbf{v}_h \in \mathbf{V}_h\\[2ex]
 *    c(\mathbf{u}_h, q_h) &= 0\quad\quad\forall q_h \in Q_h
 * \end{align*}
 * \f]
 * where \f$ \mathbf{V}_h \f$ and \f$ Q_h \f$ are appropriate function spaces for our P2P1 Taylor Hood elements which is an inf-sup stable pair for velocity and pressure. We don't have a Neumann integral because the problem we have to solve, specifies Dirichlet boundary condition on all of the boundaries.
 * 
 * The terms that we have to implement in HyTeG are as follows,
 * \f[
 * \begin{align*}
 *    a(\mathbf{u}_h, \mathbf{v}_h) &= \int_{\Omega} \nabla\mathbf{u}_h : \nabla\mathbf{v}_h d\Omega\\[2ex]
 *    b(p_h, \mathbf{v}_h) &= \int_{\Omega} \nabla p \cdot \mathbf{v}_h d\Omega\\[2ex]
 *    c(\mathbf{u}_h, q_h) &= \int_{\Omega} (\nabla\cdot\mathbf{u}_h)q_h d\Omega
 * \end{align*}
 * \f]
 * 
 * \section CouetteImplementation Implementation
 * 
 * Let us start with our implementation in HyTeG,
 * 
 * \subsection env Environment
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp Create Environment
 * 
 * Here we set up the MPI environment with walberla functions and also the PETSc manager if needed.
 * 
 * \subsection readPrm Read parameter file
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp Read prm file
 * 
 * \subsection setupMesh Setup the annulus mesh
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp annulus mesh
 * HyTeG provides a function `MeshInfo::meshAnnulus` to create an Annulus mesh with an inner and outer radius. Here we also specify the type of mesh CRISS or CROSS which decides the orientation of triangles and also the number of initial refinements in the radial and tangential direction. 
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp setup setupstorage
 * Then the `SetupPrimitiveStorage` object is set with the annulus mesh and the MPI manager which sets up the distributed data structure in the MPI processes to store the mesh
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp blending map
 * When the mesh is refined upto the maximum level from an initial coarse mesh, the curvature of the annulus will not be captured by the fine mesh. The appropriate blending map (`AnnulusMap` here) is set and the mesh is mapped accordingly. This map is also essential for blending compatible operators which takes this into account in the finite element integration.
 * 
 * \subsection BCs Boundary conditions
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp boundary cond
 * 
 * \subsection FEOps FE Operators
 * 
 * \subsection gmg Multigrid solver
 * 
 * \subsection err Error calculation
 *
 * \subsection res Results
 */

#include <iostream>

#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"

#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"

#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"

#include "hyteg/geometry/AnnulusMap.hpp"

#include "hyteg/types/PointND.hpp"
#include "hyteg/MeshQuality.hpp"

#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"

#include "core/math/Constants.h"
#include "core/logging/Logging.h"

#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"

#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"


#include <iomanip>

using walberla::real_t;

std::function< real_t( const hyteg::Point3D& ) > radius = []( const hyteg::Point3D& x ) {
   return std::sqrt( x[0] * x[0] + x[1] * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > angle = []( const hyteg::Point3D& x ) { return std::atan2( x[1], x[0] ); };

template < typename FunctionType, typename MassOperator >
real_t normL2( const FunctionType& u, const FunctionType& tmp, const MassOperator& M, const uint_t& level, const hyteg::DoFType& flag )
{
   // tmp.interpolate( 0, level );
   M.apply( u, tmp, level, flag );
   return std::sqrt( u.dotGlobal( tmp, level, flag ) );
}

real_t convAnalysis(const walberla::Config::BlockHandle& mainConf, const uint_t minLevel, const uint_t maxLevel, real_t& hMin)
{
   /* Setup Mesh */
/// [annulus mesh]
   const real_t rMin = real_c(0.5), rMax = real_c(1.5);
   hyteg::MeshInfo meshInfo = hyteg::MeshInfo::meshAnnulus(rMin, rMax, hyteg::MeshInfo::CRISS, mainConf.getParameter<uint_t>("annulusNTan"), mainConf.getParameter<uint_t>("annulusNRad"));
/// [annulus mesh]

/// [setup setupstorage]
   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
/// [setup setupstorage]

/// [blending map]
   hyteg::AnnulusMap::setMap(*setupStorage);
/// [blending map]

/// [setup storage]
   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 3 );
/// [setup storage]

/// [boundary cond]
   hyteg::BoundaryCondition bcVelocity;
   bcVelocity.createDirichletBC("DirichletInnerAndOuter", {hyteg::MeshInfo::hollowFlag::flagInnerBoundary, hyteg::MeshInfo::hollowFlag::flagOuterBoundary});
/// [boundary cond]

   /* Value of minimum cell length */

   hMin = hyteg::MeshQuality::getMinimalEdgeLength( storage, maxLevel );

   /* VTK Output */

   hyteg::VTKOutput vtkOutput("./output", "SimpleAnnulus", storage);

/// [declare p2p1funcs]
   hyteg::P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel + 1, bcVelocity );
   hyteg::P2P1TaylorHoodFunction< real_t > rhs( "rhs", storage, minLevel, maxLevel, bcVelocity );
   hyteg::P2P1TaylorHoodFunction< real_t > residual( "res", storage, minLevel, maxLevel, bcVelocity );
/// [declare p2p1funcs]

   vtkOutput.add( u );
   vtkOutput.add( rhs );

   real_t uInner = real_c(4.0);
   real_t uOuter = real_c(-8.0);

   real_t rMean = (rMin + rMax) / real_c(2.0);

   std::function< real_t( const hyteg::Point3D& ) > boundaryConditionsX = [uInner, uOuter, rMean]( const hyteg::Point3D& x ) 
   {
      real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
      
      real_t sin_theta = x[1]/r;

      if(r < rMean)
      {
         return uInner * sin_theta;
      }
      else if(r > rMean)
      {
         return uOuter * sin_theta;
      }
      else return real_c(0.0);
   };

   std::function< real_t( const hyteg::Point3D& ) > boundaryConditionsY = [uInner, uOuter, rMean]( const hyteg::Point3D& x ) 
   {
      real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
      
      real_t cos_theta = -x[0]/r;

      if(r < rMean)
      {
         return uInner * cos_theta;
      }
      else if(r > rMean)
      {
         return uOuter * cos_theta;
      }
      else return real_c(0.0);
   };

   u.uvw().interpolate({boundaryConditionsX, boundaryConditionsY}, maxLevel, hyteg::DirichletBoundary);

   typedef hyteg::P2P1ElementwiseBlendingStokesOperator       StokesOperator;
   typedef hyteg::P2ElementwiseBlendingVectorMassOperator     MassOperator;

   StokesOperator stokesOperator( storage, minLevel, maxLevel );

   const uint_t coarseIter     = mainConf.getParameter<uint_t>("coarseIter");
   const uint_t fineIter       = mainConf.getParameter<uint_t>("fineIter");

   /* Multigrid solver */

   const uint_t uzawaPre       = 10;
   const uint_t uzawaPost      = 10;
   const uint_t innerJacSmooth = 4;
   const real_t uzawaOmega     = real_c( 0.37 );
   const real_t jacobiOmega    = real_c( 0.66 );
   const uint_t numIterations  = 10;

   auto coarseGridSolver = hyteg::solvertemplates::stokesMinResSolver< StokesOperator >( storage, minLevel, real_c( 1e-14 ), coarseIter );

   auto fineGridSolver = hyteg::solvertemplates::stokesMinResSolver< StokesOperator >( storage, maxLevel, real_c( 1e-14 ), fineIter );

   auto restriction    = std::make_shared< hyteg::P2P1StokesToP2P1StokesRestriction >( true );
   auto prolongation   = std::make_shared< hyteg::P2P1StokesToP2P1StokesProlongation >();
   auto jacobiSmoother = std::make_shared< hyteg::WeightedJacobiSmoother< StokesOperator::VelocityOperator_T > >(
       storage, minLevel, maxLevel, jacobiOmega );
   auto uzawaVelocityPreconditioner =
       std::make_shared< hyteg::StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, jacobiSmoother );
   auto uzawaSmoother = std::make_shared< hyteg::UzawaSmoother< StokesOperator > >(
       storage, uzawaVelocityPreconditioner, minLevel, maxLevel, uzawaOmega, hyteg::Inner | hyteg::NeumannBoundary, innerJacSmooth );
   auto gmgSolver = std::make_shared< hyteg::GeometricMultigridSolver< StokesOperator > >(
       storage, uzawaSmoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, uzawaPre, uzawaPost, 2 );

   /* Solve the system */

   for( uint_t i = 0; i < numIterations; ++i )
   {
      gmgSolver->solve( stokesOperator, u, rhs, maxLevel );
   }
   
   /* Calculate residual */

   hyteg::P2P1TaylorHoodFunction<real_t> Au("Au", storage, minLevel, maxLevel);

   stokesOperator.apply( u, Au, maxLevel, hyteg::Inner );

   residual.assign( {real_c(1.0), real_c(-1.0)}, {rhs, Au}, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

   vtkOutput.add( residual );

   /* Calculate error with analytical solution */

   hyteg::P2P1TaylorHoodFunction< real_t > u_analytical( "u_analytical", storage, minLevel, maxLevel + 1, bcVelocity );
   hyteg::P2P1TaylorHoodFunction< real_t > error( "error", storage, minLevel, maxLevel + 1, bcVelocity );
   hyteg::P2VectorFunction< real_t > tmp("tmp", storage, minLevel, maxLevel + 1, bcVelocity);
   MassOperator M(storage, minLevel, maxLevel + 1);

   std::cout << std::setprecision(14);

   std::function< real_t( const hyteg::Point3D& ) > analyticalU = [uInner, uOuter, rMin, rMax]( const hyteg::Point3D& x )
   {
      real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
      real_t theta = std::atan2(x[1], x[0]);

      real_t c2 = (uOuter - (uInner * rMax / rMin)) * rMax * rMin * rMin / ((rMin * rMin) - (rMax * rMax));
      real_t c1 = (uInner / rMin) - (c2 / (rMin * rMin));
      real_t u_theta = c1 * r + c2 * (real_c(1.0)/r);

      return u_theta * std::sin(theta);
   };

   std::function< real_t( const hyteg::Point3D& ) > analyticalV = [uInner, uOuter, rMin, rMax]( const hyteg::Point3D& x )
   {
      real_t r = std::sqrt(x[0] * x[0] + x[1] * x[1]);
      real_t theta = std::atan2(x[1], x[0]);

      real_t c2 = (uOuter - (uInner * rMax / rMin)) * rMax * rMin * rMin / ((rMin * rMin) - (rMax * rMax));
      real_t c1 = (uInner / rMin) - (c2 / (rMin * rMin));
      real_t u_theta = c1 * r + c2 * (real_c(1.0)/r);

      return real_c(-1.0) * u_theta * std::cos(theta);
   };

   u_analytical.uvw().interpolate({analyticalU, analyticalV}, maxLevel, hyteg::All);

   hyteg::P2P1StokesToP2P1StokesProlongation StokesProlongation;

   StokesProlongation.prolongate(u, maxLevel, hyteg::All);
   StokesProlongation.prolongate(u_analytical, maxLevel, hyteg::All);
   

   error.assign({real_c(1.0), real_c(-1.0)}, {u, u_analytical}, maxLevel);
   error.assign({real_c(1.0), real_c(-1.0)}, {u, u_analytical}, maxLevel + 1);

   vtkOutput.add(u_analytical);
   vtkOutput.add(error);

   /* Write all out */

   vtkOutput.write(maxLevel);

   hyteg::VTKOutput vtkOutputHigher("./output", "SimpleAnnulus", storage);

   vtkOutputHigher.add(u_analytical);
   vtkOutputHigher.add(error);
   vtkOutputHigher.write(maxLevel + 1);

   real_t errorUV = normL2(error.uvw(), tmp, M, maxLevel + 1, hyteg::Inner);

   real_t residualNorm = std::sqrt(residual.uvw().dotGlobal(residual.uvw(), maxLevel, hyteg::All));

   uint_t nDOF = error.getNumberOfGlobalDoFs(maxLevel);

   real_t errorEnergyNorm = std::sqrt(error.uvw().dotGlobal(error.uvw(), maxLevel, hyteg::Inner) / real_c(nDOF));

   WALBERLA_ROOT_SECTION()
   {
      std::cout << "Residual = " << residualNorm << std::endl << "errorUV = " << errorUV << std::endl << "U energy norm = " << errorEnergyNorm << std::endl;
   }

   return errorUV;
}

/// [Create Environment]
int main(int argc, char* argv[])
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif
/// [Create Environment]

   /* check if a config was given on command line or load default file otherwise */

/// [Read prm file]
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./CouetteFlow.prm" );
   }
   else
   {
      cfg = env.config();
   }
/// [Read prm file]

   /* Parameters */
   const walberla::Config::BlockHandle mainConf    = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   real_t hMin = real_c(0.0);

   std::ofstream file("error_analysis.txt", std::ofstream::app);

   const uint minLevel = mainConf.getParameter<uint>("minLevel");
   const uint maxLevel = mainConf.getParameter<uint>("maxLevel");

   for(uint_t level = 1U; level < maxLevel; level++)
   {
      real_t l2Error = convAnalysis(mainConf, minLevel, level, hMin);
      WALBERLA_MPI_BARRIER()
      WALBERLA_ROOT_SECTION()
      {
         std::cout << walberla::format( "Value of minimum edge length and L2 Error for maxLevel %d = %10.10f, %10.10f", level, hMin, l2Error) << std::endl;
         file << walberla::format( "%10.10f, %10.10f", hMin, l2Error) << std::endl;
      }
   }

   file.close();

   return 0;
}