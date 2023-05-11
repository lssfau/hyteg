/*
 * Copyright (c) 2023 Ponsuganth Ilangovan P, Marcus Mohr.
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
 * \page 11_CouetteFlow Tutorial to implement simple Couette flow
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
 * The weak finite element formulation of the boundary value problem is to find \f$ \mathbf{u}_h \in \mathbf{V}_h \f$ and \f$ p_h \in Q_h \f$ such that,
 * \f[
 * \begin{align*}
 *    a(\mathbf{u}_h, \mathbf{v}_h) + b(p_h, \mathbf{v}_h) &= \mathbf{0}\quad\quad\forall\mathbf{v}_h \in \mathbf{V}_h\\[2ex]
 *    b(q_h, \mathbf{u}_h) &= 0\quad\quad\forall q_h \in Q_h
 * \end{align*}
 * \f]
 * where \f$ \mathbf{V}_h \f$ and \f$ Q_h \f$ are appropriate function spaces for our P2P1 Taylor Hood elements which is an inf-sup stable pair for velocity and pressure.
 * 
 * The terms that we have to implement in HyTeG are as follows,
 * \f[
 * \begin{align*}
 *    a(\mathbf{u}_h, \mathbf{v}_h) &= \int_{\Omega} \nabla\mathbf{u}_h : \nabla\mathbf{v}_h d\Omega\\[2ex]
 *    b(p_h, \mathbf{v}_h) &= \int_{\Omega} (\nabla\cdot\mathbf{v}_h)p_h d\Omega
 * \end{align*}
 * \f]
 * 
 * \section CouetteImplementation Implementation
 * 
 * Let us start with our implementation in HyTeG,
 * 
 * Here we set up the MPI environment with walberla functions and also the PETSc manager if needed.
 * 
 * \subsection env Environment
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp Create Environment
 * 
 * \subsection readPrm Read parameter file
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp Read prm file
 * 
 * HyTeG provides a function `MeshInfo::meshAnnulus` to create an Annulus mesh with an inner and outer radius. Here we also specify the type of mesh CRISS or CROSS which decides the orientation of triangles and also the number of initial refinements in the radial and tangential direction. 
 * 
 * \subsection setupMesh Setup the annulus mesh
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp annulus mesh
 * 
 * Then the `SetupPrimitiveStorage` object is set with the annulus mesh and the MPI manager which sets up the distributed data structure in the MPI processes to store the mesh
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp setup setupstorage
 * 
 * When the mesh is refined upto the maximum level from an initial coarse mesh, the curvature of the annulus will not be captured by the fine mesh. The appropriate blending map (`AnnulusMap` here) is set and the mesh is mapped accordingly. This map is also essential for blending compatible operators which takes this into account in the finite element integration.
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp blending map
 * 
 * \subsection BCs Boundary conditions
 * 
 * The `meshAnnulus` function automatically tags the inner and outer boundary of the Annulus, and hence those flags are used to set the required boundary condition, which in our case we set the Dirichlet boundary condition on the inner and outer boundaries of the annulus. Now this boundary condition object `bcVelocity` must be passed to the FE functions which needs to use this boundary condition.
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp boundary cond
 * 
 * Here `u` is the P2P1 pair FE function which HyTeG provides to hold the velocity and pressure values in a same class. The value of the velocity `u.uvw()` boundary condition (BC) is set using the vector interpolate function by passing the appropriate lambda functions which sets the BC values on the Dirichlet boundary and on the level specified.
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp bcValues
 * 
 * The lambda functions `boundaryConditionsX` and `boundaryConditionsY` are used to set the boundary values for velocity. In our example, the velocities are radial in the boundaries and hence we calculate the appropriate cartesian components and impose them on the boundary. The lambda functions that are passed will receive the coordinates of the boundary point and the user has to check it's position and impose the correct values by returning the same.
 * 
 * \subsection FEfns FE Functions
 * 
 * Here we define the finite element functions needed for our FE computations and error norm calculations. These functions store the DOF values in a distributed fashion on the MPI processes that are running the program.
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp declare p2p1funcs
 * 
 * \subsection FEOps FE Operators
 * 
 * To implement \f$ a(\mathbf{u}_h, \mathbf{v}_h) , b(p_h, \mathbf{v}_h), c(\mathbf{u}_h, q_h)\f$, HyTeG provides two main operators. One is the `P2P1TaylorHoodStokesOperator` but it does not use the blending map in the finite element computations. 
 * 
 * The other operator `P2P1ElementwiseBlendingStokesOperator` is compatible with blending maps and hence is a more accurate way to implement the Stokes problem for the annulus.
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp StokesOperator
 * 
 * \subsection gmg Multigrid solver
 * 
 * Here we set up the Geometric multigrid solver which requires the following,
 * - smoother \f$ \Rightarrow\f$ we use the Uzawa smoother with Jacobi smoothing for velocity.
 * - coarse grid solver \f$ \Rightarrow\f$ for this we use the pressure preconditioned MINRES solver which is wrapped as a solver template under `solvertemplates::stokesMinResSolver`.
 * - restriction and prolongation operators \f$ \Rightarrow\f$ the `P2P1StokesToP2P1StokesRestriction(Prolongation)` uses quadratic interpolation for prolongation of velocity, linear one for prolongation of pressure and the transpose of prolongation for restriction
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp gmg
 * 
 * We then call the GMG solver required number of times to solve the system on the level from which we want the V cycle to start.
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp gmgSolve
 * 
 * \subsection err Error calculation
 * 
 * The computed finite element solution is projected to a higher level and then the error is computed with the analytical solution on that level. To calculate the L2 norm of the error we would have to use the mass operator which performs the integration of the L2 norm in the FE space.
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp calcErr
 * 
 * This function basically computes the L2 norm on the finite element space,
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp normL2
 * 
 * \f[
 * \begin{equation*}
 *    ||u||_2 = \left(\int_{\Omega} u^2 d\Omega\right)^{\frac{1}{2}}
 * \end{equation*}
 * \f]
 * In Einstein summation notation,
 * \f[
 * \begin{align*}
 *    ||u_h||_2 &= \left(\sum_{e}^{}\int_{\Omega^e} {}_eu_i^h\ \phi^e_i\ \phi^e_j\ {}_eu_j^h\ d\Omega^e\right)^{\frac{1}{2}}\\
 *    ||u_h||_2 &= \left(\mathbf{U}^\top \mathbf{M} \mathbf{U}\right)^{\frac{1}{2}}
 * \end{align*}
 * \f]
 * 
 * Hence we define the appropriate mass operator for our FE space which also supports blending for the annulus to compute the norm of the error.
 * 
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp massOpErr
 * \snippet tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp calcErrNorm
 *
 * \subsection res Results
 * 
 * ## Convergence plot for P2P1 Taylor Hood elements
 * <img src="11_CouetteFlow_Convergence.png" width="75%"/>
 * 
 * ## Velocity magnitude contour plot on the annulus
 * <img src="11_CouetteFlow_Velocity.png" width="50%"/> 
 * 
 * ## Velocity magnitude contour plot on the annulus
 * <img src="11_CouetteFlow_Pressure.png" width="50%"/> 
 * 
 * \section fullCode Code without comments
 * \include tutorials/11_CouetteFlow/11_CouetteFlow2D.cpp
 */

#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/operators/VectorMassOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

using walberla::real_t;

using namespace hyteg;

/// [normL2]
template < typename FunctionType, typename MassOperator >
real_t normL2( const FunctionType& u, const FunctionType& tmp, const MassOperator& M, uint_t level, DoFType flag )
{
   M.apply( u, tmp, level, flag );
   return std::sqrt( u.dotGlobal( tmp, level, flag ) );
}
/// [normL2]

real_t convAnalysis( const walberla::Config::BlockHandle& mainConf, uint_t minLevel, uint_t maxLevel, real_t& hMin )
{
   /// [annulus mesh]
   const real_t rMin = real_c( 0.5 ), rMax = real_c( 1.5 );
   MeshInfo     meshInfo = MeshInfo::meshAnnulus( rMin,
                                              rMax,
                                              MeshInfo::CRISS,
                                              mainConf.getParameter< uint_t >( "annulusNTan" ),
                                              mainConf.getParameter< uint_t >( "annulusNRad" ) );
   /// [annulus mesh]

   /// [setup setupstorage]
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   /// [setup setupstorage]

   /// [blending map]
   AnnulusMap::setMap( *setupStorage );
   /// [blending map]

   /// [setup storage]
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 3 );
   /// [setup storage]

   /// [boundary cond]
   BoundaryCondition bcVelocity;
   bcVelocity.createDirichletBC( "DirichletInnerAndOuter",
                                 { MeshInfo::hollowFlag::flagInnerBoundary, MeshInfo::hollowFlag::flagOuterBoundary } );
   /// [boundary cond]

   hMin = MeshQuality::getMinimalEdgeLength( storage, maxLevel );

   VTKOutput vtkOutput( "./output", "SimpleAnnulus", storage );

   /// [declare p2p1funcs]
   P2P1TaylorHoodFunction< real_t > u( "u", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > rhs( "rhs", storage, minLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > residual( "res", storage, maxLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > Au( "Au", storage, maxLevel, maxLevel );
   P2P1TaylorHoodFunction< real_t > uAnalytical( "uAnalytical", storage, maxLevel, maxLevel, bcVelocity );
   P2P1TaylorHoodFunction< real_t > error( "error", storage, minLevel, maxLevel + 1, bcVelocity );
   P2VectorFunction< real_t >       tmp( "tmp", storage, minLevel, maxLevel + 1, bcVelocity );
   /// [declare p2p1funcs]

   vtkOutput.add( u );
   vtkOutput.add( rhs );

   real_t uInner = real_c( 4.0 );
   real_t uOuter = real_c( -8.0 );

   real_t rMean = ( rMin + rMax ) / real_c( 2.0 );

   std::function< real_t( const Point3D& ) > boundaryConditionsX = [uInner, uOuter, rMean]( const Point3D& x ) {
      real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] );

      real_t sin_theta = x[1] / r;

      if ( r < rMean )
      {
         return uInner * sin_theta;
      }
      else if ( r > rMean )
      {
         return uOuter * sin_theta;
      }
      else
         return real_c( 0.0 );
   };

   std::function< real_t( const Point3D& ) > boundaryConditionsY = [uInner, uOuter, rMean]( const Point3D& x ) {
      real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] );

      real_t cos_theta = -x[0] / r;

      if ( r < rMean )
      {
         return uInner * cos_theta;
      }
      else if ( r > rMean )
      {
         return uOuter * cos_theta;
      }
      else
         return real_c( 0.0 );
   };

   /// [bcValues]
   u.uvw().interpolate( { boundaryConditionsX, boundaryConditionsY }, maxLevel, DirichletBoundary );
   /// [bcValues]

   /// [StokesOperator]
   typedef P2P1ElementwiseBlendingStokesOperator StokesOperator;
   StokesOperator                                stokesOperator( storage, minLevel, maxLevel );
   /// [StokesOperator]

   const uint_t coarseIter = mainConf.getParameter< uint_t >( "coarseIter" );

   const uint_t uzawaPre       = 10;
   const uint_t uzawaPost      = 10;
   const uint_t innerJacSmooth = 4;
   const real_t uzawaOmega     = real_c( 0.37 );
   const real_t jacobiOmega    = real_c( 0.66 );
   const uint_t numIterations  = 10;

   /// [gmg]
   auto coarseGridSolver =
       solvertemplates::stokesMinResSolver< StokesOperator >( storage, minLevel, real_c( 1e-14 ), coarseIter );

   auto restriction    = std::make_shared< P2P1StokesToP2P1StokesRestriction >( true );
   auto prolongation   = std::make_shared< P2P1StokesToP2P1StokesProlongation >();
   auto jacobiSmoother = std::make_shared< WeightedJacobiSmoother< StokesOperator::VelocityOperator_T > >(
       storage, minLevel, maxLevel, jacobiOmega );
   auto uzawaVelocitySmoother =
       std::make_shared< StokesVelocityBlockBlockDiagonalPreconditioner< StokesOperator > >( storage, jacobiSmoother );
   auto uzawaSmoother = std::make_shared< UzawaSmoother< StokesOperator > >(
       storage, uzawaVelocitySmoother, minLevel, maxLevel, uzawaOmega, Inner | NeumannBoundary, innerJacSmooth );
   auto gmgSolver = std::make_shared< GeometricMultigridSolver< StokesOperator > >(
       storage, uzawaSmoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, uzawaPre, uzawaPost, 2 );

   /// [gmg]

   /// [gmgSolve]
   for ( uint_t i = 0; i < numIterations; ++i )
   {
      gmgSolver->solve( stokesOperator, u, rhs, maxLevel );
   }
   /// [gmgSolve]

   stokesOperator.apply( u, Au, maxLevel, Inner );

   residual.assign( { real_c( 1.0 ), real_c( -1.0 ) }, { rhs, Au }, maxLevel, Inner | NeumannBoundary );

   vtkOutput.add( residual );

   /// [massOpErr]
   typedef P2ElementwiseBlendingVectorMassOperator MassOperator;

   MassOperator M( storage, minLevel, maxLevel + 1 );
   /// [massOpErr]

   std::function< real_t( const Point3D& ) > analyticalU = [uInner, uOuter, rMin, rMax]( const Point3D& x ) {
      real_t r     = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t theta = std::atan2( x[1], x[0] );

      real_t c2      = ( uOuter - ( uInner * rMax / rMin ) ) * rMax * rMin * rMin / ( ( rMin * rMin ) - ( rMax * rMax ) );
      real_t c1      = ( uInner / rMin ) - ( c2 / ( rMin * rMin ) );
      real_t u_theta = c1 * r + c2 * ( real_c( 1.0 ) / r );

      return u_theta * std::sin( theta );
   };

   std::function< real_t( const Point3D& ) > analyticalV = [uInner, uOuter, rMin, rMax]( const Point3D& x ) {
      real_t r     = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t theta = std::atan2( x[1], x[0] );

      real_t c2      = ( uOuter - ( uInner * rMax / rMin ) ) * rMax * rMin * rMin / ( ( rMin * rMin ) - ( rMax * rMax ) );
      real_t c1      = ( uInner / rMin ) - ( c2 / ( rMin * rMin ) );
      real_t u_theta = c1 * r + c2 * ( real_c( 1.0 ) / r );

      return real_c( -1.0 ) * u_theta * std::cos( theta );
   };

   uAnalytical.uvw().interpolate( { analyticalU, analyticalV }, maxLevel, All );

   /// [calcErr]
   P2P1StokesToP2P1StokesProlongation StokesProlongation;

   error.assign( { real_c( 1.0 ), real_c( -1.0 ) }, { u, uAnalytical }, maxLevel );
   StokesProlongation.prolongate( error, maxLevel, All );
   /// [calcErr]

   vtkOutput.add( uAnalytical );
   vtkOutput.add( error );

   vtkOutput.write( maxLevel );

   VTKOutput vtkOutputHigher( "./output", "SimpleAnnulus", storage );

   vtkOutputHigher.add( error );
   vtkOutputHigher.write( maxLevel + 1 );

   /// [calcErrNorm]
   real_t errorUV = normL2( error.uvw(), tmp, M, maxLevel + 1, All );
   /// [calcErrNorm]

   real_t residualNorm = std::sqrt( residual.uvw().dotGlobal( residual.uvw(), maxLevel, All ) );

   WALBERLA_ROOT_SECTION()
   {
      std::cout << "Residual = " << residualNorm << std::endl << "errorUV = " << errorUV << std::endl;
   }

   return errorUV;
}

/// [Create Environment]
int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   /// [Create Environment]

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

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   real_t hMin = real_c( 0.0 );

   std::ofstream file( "error_analysis.txt", std::ofstream::app );

   const uint minLevel = mainConf.getParameter< uint >( "minLevel" );
   const uint maxLevel = mainConf.getParameter< uint >( "maxLevel" );

   for ( uint_t level = 1U; level < maxLevel; level++ )
   {
      real_t l2Error = convAnalysis( mainConf, minLevel, level, hMin );
      WALBERLA_MPI_BARRIER()
      WALBERLA_ROOT_SECTION()
      {
         std::cout << walberla::format(
                          "Value of minimum edge length and L2 Error for maxLevel %d = %10.10f, %10.10f", level, hMin, l2Error )
                   << std::endl;
         file << walberla::format( "%10.10f, %10.10f", hMin, l2Error ) << std::endl;
      }
   }

   file.close();

   return 0;
}
