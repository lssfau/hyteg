/*
 * Copyright (c) 2024 Ponsuganth Ilangovan P
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg_operators/operators/k_mass/P2ToP1ElementwiseKMass.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "terraneo/operators/TransportOperatorStd.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"

/**
 * \page GB.02_KingCompressible Tutorial GB.02 - Compressible Benchmark
 *
 * \dontinclude tutorials/geo-benchmarks/GB.02_KingCompressible/GB.02_KingCompressible.cpp
 * 
 * \brief We have demonstrated the incompressible Stokes-Transport benchmark (Blankenbach) in the previous tutorial. 
 * In this one we will implement a compressible benchmark from Geophysical literature, specifically the one from
 * <a href="https://doi.org/10.1111/j.1365-246X.2009.04413.x">A community benchmark for 2-D Cartesian compressible convection in the Earth's mantle,
 * King et al., 2010, GJI</a>.
 * It involves solving the respective compressible equations on unit square and calculation of Nusselt numbers and 
 * velocity RMS values to verify the code.
 * 
 * \section GB02-KingCompressible-GoverningEquations Model and Equations
 * 
 * The governing equations are basically the Stokes and energy equations,
 * 
 * \f{align*}{
 *  -\nabla\cdot\tau + \nabla p &= \text{Ra}T \\
 *  \nabla\cdot(\rho u) &= 0 \\
 *  \frac{\partial T}{\partial t} + u \cdot \nabla T - \nabla \cdot \left(\frac{\kappa}{\bar{\rho}C_p} 
 *  \nabla T\right) &= \frac{\alpha \text{Di}}{C_p} (u\cdot g) T + \frac{\text{Di}}{\text{Ra}C_p}\tau:\epsilon
 * \f}
 * 
 * The initial conditions for temperature is prescribed as,
 * 
 * \f{equation*}{
 *  T(x, y, t = 0) = (1-y) + A\cos{\pi x}\sin{\pi y}
 * \f}
 * 
 * where $x$ and $y$ are the coordinates of the unit square with origin at bottom-left. 
 * Freeslip boundary conditions is imposed on all four walls of the square.
 * For the temperature field, a Dirichlet boundary condition of T = 0 and T = 1 is imposed 
 * at the top and bottom and zero flux on the sides. This induces a single convection cell 
 * in the square and eventually reaches a steady state. Once this state is reached, 
 * the Nusselt number values are calculated and verified. 
 * As mentioned in the article mentioned above, if we start the model at a higher Rayleigh number 
 * we may not converge to a single cell solution. Hence we start at lower Rayleigh numbers 
 * and store the checkpoints, which we then use to initialise the model at a higher Rayleigh number.
 * 
 * \section GB02-KingCompressible-Domain Domain
 * 
 * Through the geometry module in HyTeG, a rectangle mesh can be created with the `meshRectangle`, 
 * which is used to create the unit square mesh. In addition to creating the macro mesh with required subdivisions, 
 * the boundary nodes must be marked properly which will be used downstream when defining the boundary conditions
 * for the finite element functions, which is what is done in the following snippet.
 * 
 * \snippet{trimleft} this SetupStorageAndMarkings
 * 
 * \section GB02-KingCompressible-BCs Boundary Conditions
 * 
 * The boundary conditions for the Temperature and Velocity are defined as follows. For temperature, a Dirichlet is defined
 * on the top and bottom and zero flux on the side walls.
 * 
 * \snippet{trimleft} this BoundaryConditionsTemperature
 * 
 * For the velocity freeslip boundary conditions must be imposed on all four walls. For this, 
 * as we have straight boundaries here, the axis aligned coordinates can be seperately defined a Dirichlet 
 * and Neumann condition which inturn will impose the freeslip condition (no-normal flow).
 * 
 * \snippet{trimleft} this BoundaryConditionsVelocity
 * 
 * \section GB02-KingCompressible-OpStokes Operators -- Stokes
 * 
 * The Stokes operator can be set up similar to previous tutorials, but here, to treat the compressibility constraint, 
 * we expand the mass conservation with chain rule and treat some terms explicitly to 
 * make sure, we have the same symmetric system matrix.
 * 
 * \snippet{trimleft} this StokesFrozenVelocity
 * 
 * \section GB02-KingCompressible-OpEnergy Operators -- Energy
 * 
 * For treating the advection part of the energy equation, we use the method of modified characteristics (MMOC).
 * Hence this advection operator can be defined with
 * 
 * \snippet{trimleft} this MMOCForTransport
 * 
 * As the operator splitting approach is used, this must be taken care in the time stepping algorithm. We use an implicit Euler
 * scheme for timestepping the diffusion and other terms of the energy equations, while the MMOC is used to step for advection. 
 * The `apply` function of the `P2TransportTALAOperator` applies the time discretized form of the weak form considering only 
 * the terms other than advection.
 * 
 * \section GB02-KingCompressible-SolverStokes Solver -- Stokes
 *  
 * The boundary conditions are set with the correct Dirichlet/Neumann flag to ensure that no-normal flow is applied on the walls. 
 * This then helps in creating a single convection cell in the unit square. We start the model at a Rayleigh number of
 * Ra \f$ = 10^4 \f$, then we store a checkpoint with this solution. This is then used to initialise the model at 
 * Ra \f$ = 10^5 \f$, and then eventually to even higher values for testing. But here we use the values of Ra \f$ = 10^5 \f$
 * and Di \f$ = 0.25 \f$ to calculate the Nusselt numbers for comparison.
 * 
 * Minres or PETSc direct solver can be used to compute the solution to the Stokes system.
 * 
 * \snippet{trimleft} this StokesSolverLambdaFunction
 * 
 * \section GB02-KingCompressible-SolverEnergy Solver -- Energy
 * 
 * The timestep for the energy equation solver is calculated according to the CFL condition using the velocity from the 
 * Stokes solution. This is then used for the MMOC solver to step the advection, then also used for the energy equation.
 * Here we use GMRES for the energy solver as the system is nicely symmetric as advection is absent in the linear system.
 * 
 * \snippet{trimleft} this TransportSolverLambdaFunction
 * 
 * \section GB02-KingCompressible-Timestepping Timestepping
 * 
 * First the temperature \f$ T(x, t_0) \f$ is initialized to the prescribed field according to the above equation.
 * 
 * \snippet{trimleft} this TemperatureInitialization
 * 
 * Then the Stokes system is solved to get the initial velocity field, else the program can also be started from a checkpoint.
 * 
 * \snippet{trimleft} this InitialSolveOrCheckpoint
 * 
 * This is then used to step the transport solver, both advection and the remaining energy parts, with 
 * MMOC and implicit Euler respectively. With this we compute the temperature field at time \f$ t_1 \f$. Then the Stokes system is
 * solved again to obtain the corresponding velocity field at time \f$ t_1 \f$ and so on.
 * 
 * 
 * \section GB02-KingCompressible-Results Results
 * 
 * <img src="GB.02_KingCompressibleSquare.png" width="50%" />
 * 
 */

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;

using namespace terraneo;

using P2ToP1ElementwiseKMass = operatorgeneration::P2ToP1ElementwiseKMass;
// using P1ToP2ElementwiseKMass = operatorgeneration::P1ToP2ElementwiseKMass;

using P2TransportTALAOperator = terraneo::P2TransportOperator;

namespace hyteg {

const real_t boundaryMarkerThreshold = 1e-6;

std::function< bool( const Point3D& ) > bottomMarker = []( const Point3D& x ) {
   if ( std::abs( x[1] ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > rightMarker = []( const Point3D& x ) {
   if ( std::abs( x[0] - 1.0 ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > leftMarker = []( const Point3D& x ) {
   if ( std::abs( x[0] ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > topMarker = []( const Point3D& x ) {
   if ( std::abs( x[1] - 1.0 ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > cornersMarker = []( const Point3D& x ) {
   if ( ( topMarker( x ) && rightMarker( x ) ) || ( bottomMarker( x ) && rightMarker( x ) ) ||
        ( topMarker( x ) && leftMarker( x ) ) || ( bottomMarker( x ) && leftMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

using P2P1StokesOperator = operatorgeneration::P2P1StokesFullOperator;

enum BoundaryMarkers
{
   Bottom = 23,
   Right,
   Left,
   Top,
   Corners
};

struct ParameterContainer
{
   bool verbose = true;

   real_t rMin = 1.22;
   real_t rMax = 2.22;

   uint_t maxTimeSteps      = 1000;
   uint_t vtkWriteFrequency = 1U;

   bool MMOC             = true;
   bool SUPG             = false;
   bool compressible     = true;
   bool adiabaticHeating = true;
   bool shearHeating     = true;

   real_t Ra          = 1e5;
   real_t Di          = 0.5;
   real_t T0          = 0.091;
   real_t diffusivity = 1.0;
   real_t cflMax      = 0.75;
   real_t AiniPerturb = 0.1;

   real_t rho0       = 1.0;
   real_t alpha      = 1.0;
   real_t cpr        = 1.0;
   real_t cvr        = 1.0;
   real_t grueneisen = 1.0;
   real_t alphabar   = 1.0;
   real_t cpbar      = 1.0;
   real_t chibar     = 1.0;
   real_t k_         = 1.0;

   real_t minresRelTol = 1e-4;
   real_t minresAbsTol = 1e-8;
   real_t gmresTol     = 1e-5;
   uint_t minresIter   = 1000U;
   uint_t gmresIter    = 1000U;

   uint_t nsCalcFreq = 10U;
};

class P2TransportTimesteppingOperator : public Operator< P2Function< real_t >, P2Function< real_t > >
{
 public:
   P2TransportTimesteppingOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                    const uint_t                               minLevel,
                                    const uint_t                               maxLevel,
                                    real_t                                     k )
   : Operator( storage, minLevel, maxLevel )
   , diffusionOperator( storage, minLevel, maxLevel )
   , massOperator( storage, minLevel, maxLevel )
   , k_( k )
   {}

   ///[TransportOperatorApply]
   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      diffusionOperator.apply( src, dst, level, flag, updateType );
      dst.assign( { k_ * dt }, { dst }, level, flag );
      massOperator.apply( src, dst, level, flag, Add );
   }
   ///[TransportOperatorApply]

   void setDt( real_t dt_ ) { dt = dt_; }

 private:
   P2ElementwiseBlendingLaplaceOperator diffusionOperator;
   P2ElementwiseBlendingMassOperator    massOperator;

   real_t dt = 0.01;
   real_t k_ = 1.0;
};

class TALASimulation
{
 public:
   TALASimulation( const walberla::Config::BlockHandle& mainConf_,
                   std::shared_ptr< PrimitiveStorage >  storage_,
                   uint_t                               minLevel_,
                   uint_t                               maxLevel_ )
   : mainConf( mainConf_ )
   , storage( storage_ )
   , minLevel( minLevel_ )
   , maxLevel( maxLevel_ )
   , vecMassOperator( storage, minLevel_, maxLevel_ )
   , massOperator( storage, minLevel_, maxLevel_ )
   , massOperatorP1( storage, minLevel_, maxLevel_ )
   , transport( storage, minLevel_, maxLevel_, TimeSteppingScheme::RK4 )
   {
      params.rMin = mainConf.getParameter< real_t >( "rMin" );
      params.rMax = mainConf.getParameter< real_t >( "rMax" );

      endTime             = mainConf.getParameter< real_t >( "simulationTime" );
      params.maxTimeSteps = mainConf.getParameter< uint_t >( "maxTimeSteps" );

      params.vtkWriteFrequency = mainConf.getParameter< uint_t >( "vtkWriteFrequency" );

      params.AiniPerturb = mainConf.getParameter< real_t >( "AiniPerturb" );

      params.Ra = mainConf.getParameter< real_t >( "RayleighNumber" );
      params.Di = mainConf.getParameter< real_t >( "DissipationNumber" );

      params.MMOC             = mainConf.getParameter< bool >( "MMOC" );
      params.SUPG             = mainConf.getParameter< bool >( "SUPG" );
      params.compressible     = mainConf.getParameter< bool >( "compressible" );
      params.adiabaticHeating = mainConf.getParameter< bool >( "adiabaticHeating" );
      params.shearHeating     = mainConf.getParameter< bool >( "shearHeating" );

      params.minresIter   = mainConf.getParameter< uint_t >( "stokesMinresIter" );
      params.minresRelTol = mainConf.getParameter< real_t >( "stokesMinresTol" );

      params.nsCalcFreq = mainConf.getParameter< uint_t >( "nsCalcFreq" );

      params.gmresIter = mainConf.getParameter< uint_t >( "transportGmresIter" );
      params.gmresTol  = mainConf.getParameter< real_t >( "transportGmresTol" );

      normalsFS = [this]( const Point3D& x, Point3D& nx ) {
         if ( rightMarker( x ) )
         {
            nx[0] = 1.0;
            nx[1] = 0.0;
         }
         else if ( leftMarker( x ) )
         {
            nx[0] = -1.0;
            nx[1] = 0.0;
         }
         else if ( topMarker( x ) )
         {
            nx[0] = 0.0;
            nx[1] = 1.0;
         }
         else if ( bottomMarker( x ) )
         {
            nx[0] = 0.0;
            nx[1] = -1.0;
         }
         else
         {
            WALBERLA_LOG_INFO_ON_ROOT( "Probably shouldn't be here!" );
         }
      };

      tempDevBC = [this]( const Point3D& x ) {
         if ( topMarker( x ) )
         {
            return 0.0;
         }
         else if ( bottomMarker( x ) )
         {
            return 1.0;
         }
         return 0.0;
      };

      /// [TemperatureInitialization]
      tempIni = [this]( const Point3D& x ) {
         return ( 1 - x[1] ) + params.AiniPerturb * std::cos( walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[1] );
      };
      /// [TemperatureInitialization]

      TRefFunc = [this]( const Point3D& x ) { return params.T0 * std::exp( ( 1 - x[1] ) * params.Di ) - params.T0; };

      BoundaryCondition bcTemp, bcVelocity, bcVelocityX, bcVelocityY;

      /// [BoundaryConditionsTemperature]
      bcTemp.createDirichletBC( "DirichletBottomAndTop",
                                { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Corners } );
      bcTemp.createNeumannBC( "NeumannLeftAndRight", { BoundaryMarkers::Left, BoundaryMarkers::Right } );
      /// [BoundaryConditionsTemperature]

      bcVelocity.createAllInnerBC();
      // bcVelocity.createDirichletBC(
      //     "AllDirichlet", { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Left, BoundaryMarkers::Right } );
      bcVelocity.createFreeslipBC(
          "AllFreeslip", { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Left, BoundaryMarkers::Right } );
      // bcVelocity.createFreeslipBC( "DirichletCorners", { BoundaryMarkers::Corners } );
      bcVelocity.createDirichletBC( "DirichletCorners", { BoundaryMarkers::Corners } );

      /// [BoundaryConditionsVelocity]
      bcVelocityX.createAllInnerBC();
      bcVelocityX.createDirichletBC( "LRDirichlet", { BoundaryMarkers::Left, BoundaryMarkers::Right, BoundaryMarkers::Corners } );
      bcVelocityX.createNeumannBC( "TBNeumann", { BoundaryMarkers::Top, BoundaryMarkers::Bottom } );

      bcVelocityY.createAllInnerBC();
      bcVelocityY.createNeumannBC( "LRNeumann", { BoundaryMarkers::Left, BoundaryMarkers::Right } );
      bcVelocityY.createDirichletBC( "TBDirichlet", { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Corners } );
      /// [BoundaryConditionsVelocity]

      TDev     = std::make_shared< P2Function< real_t > >( "TDev", storage_, minLevel_, maxLevel_, bcTemp );
      TDevPrev = std::make_shared< P2Function< real_t > >( "TDevPrev", storage_, minLevel_, maxLevel_, bcTemp );
      TDevInt  = std::make_shared< P2Function< real_t > >( "TDevInt", storage_, minLevel_, maxLevel_, bcTemp );
      TRef     = std::make_shared< P2Function< real_t > >( "TRef", storage_, minLevel_, maxLevel_, bcTemp );
      TRefDev  = std::make_shared< P2Function< real_t > >( "TRefDev", storage_, minLevel_, maxLevel_, bcTemp );
      TRes     = std::make_shared< P2Function< real_t > >( "TRes", storage_, minLevel_, maxLevel_, bcTemp );
      TRhs     = std::make_shared< P2Function< real_t > >( "TRhs", storage_, minLevel_, maxLevel_, bcTemp );
      rhoP2    = std::make_shared< P2Function< real_t > >( "rhoP2", storage_, minLevel_, maxLevel_, bcTemp );
      rhoInvP2 = std::make_shared< P2Function< real_t > >( "rhoInvP2", storage_, minLevel_, maxLevel_, bcTemp );
      viscP2   = std::make_shared< P2Function< real_t > >( "viscP2", storage_, minLevel_, maxLevel_ );

      Ttemp = std::make_shared< P2Function< real_t > >( "Ttemp", storage_, minLevel_, maxLevel_ );

      zero = std::make_shared< P2Function< real_t > >( "zero", storage_, minLevel_, maxLevel_, bcTemp );

      gradRhoByRhoP2 = std::make_shared< P2VectorFunction< real_t > >( "gradRhoByRhoP2", storage_, minLevel_, maxLevel_ );

      u     = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "u", storage_, minLevel_, maxLevel_, bcVelocity );
      uRes  = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRes", storage_, minLevel_, maxLevel_, bcVelocity );
      uPrev = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPrev", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhs  = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhs", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhsStrong =
          std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhsStrong", storage_, minLevel_, maxLevel_, bcVelocity );
      uTemp = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uTemp", storage_, minLevel_, maxLevel_, bcVelocity );

      u->uvw().setBoundaryCondition( bcVelocityX, 0U );
      u->uvw().setBoundaryCondition( bcVelocityY, 1U );

      uRhs->uvw().setBoundaryCondition( bcVelocityX, 0U );
      uRhs->uvw().setBoundaryCondition( bcVelocityY, 1U );

      transportOp = std::make_shared< P2TransportTimesteppingOperator >( storage_, minLevel_, maxLevel_, params.diffusivity );

      transportTALAOp = std::make_shared< P2TransportTALAOperator >( storage_, minLevel_, maxLevel_ );

      gradRhoByRhoP2->component( 0U ).interpolate( 0.0, maxLevel_, All );

      if ( params.compressible )
      {
         rhoFunc = [this]( const Point3D& x ) { return params.rho0 * std::exp( ( 1 - x[1] ) * params.alpha / params.Di ); };
         gradRhoByRhoP2->component( 1U ).interpolate( params.alpha / params.Di, maxLevel_, All );
      }
      else
      {
         rhoFunc = [this]( const Point3D& ) { return 1.0; };
         gradRhoByRhoP2->component( 1U ).interpolate( 0.0, maxLevel_, All );
      }

      rhoInvFunc = [this]( const Point3D& x ) { return 1.0 / rhoFunc( x ); };

      gradRhoByRhoX =
          std::make_shared< P2ToP1ElementwiseKMass >( storage_, minLevel_, maxLevel_, gradRhoByRhoP2->component( 0U ) );
      gradRhoByRhoY =
          std::make_shared< P2ToP1ElementwiseKMass >( storage_, minLevel_, maxLevel_, gradRhoByRhoP2->component( 1U ) );

      transportTALAOp->setVelocity( u );
      transportTALAOp->setViscosity( viscP2 );
      transportTALAOp->setTemperature( TDev );

      invGravityField = std::make_shared< P2VectorFunction< real_t > >( "invGravityField", storage_, minLevel_, maxLevel_ );

      diffusionTermCoeff    = std::make_shared< P2Function< real_t > >( "diffusionTermCoeff", storage_, minLevel_, maxLevel_ );
      adiabaticTermCoeff    = std::make_shared< P2Function< real_t > >( "adiabaticTermCoeff", storage_, minLevel_, maxLevel_ );
      shearHeatingTermCoeff = std::make_shared< P2Function< real_t > >( "shearHeatingTermCoeff", storage_, minLevel_, maxLevel_ );
      constEnergyCoeff      = std::make_shared< P2Function< real_t > >( "constEnergyCoeff", storage_, minLevel_, maxLevel_ );
      surfTempCoeff         = std::make_shared< P2Function< real_t > >( "surfTempCoeff", storage_, minLevel_, maxLevel_ );

      rhoP2->interpolate( rhoFunc, maxLevel_, All );
      rhoInvP2->interpolate( rhoInvFunc, maxLevel_, All );

      diffusivityCoeffFunc = [this]( const Point3D& x ) { return params.k_ / ( params.cpbar * rhoFunc( x ) ); };

      adiabaticCoeffFunc = [this]( const Point3D& ) { return params.alphabar * params.Di / params.cpbar; };

      constEnergyCoeffFunc = [this]( const Point3D& ) { return 0.0; };

      surfTempCoeffFunc = [this]( const Point3D& ) { return 0.0; };

      invGravityX = [this]( const Point3D& ) { return 0.0; };
      invGravityY = [this]( const Point3D& ) { return 1.0; };

      diffusionTermCoeff->interpolate( params.k_ / params.cpbar, maxLevel_, All );
      diffusionTermCoeff->multElementwise( { *diffusionTermCoeff, *rhoInvP2 }, maxLevel_, All );

      adiabaticTermCoeff->interpolate( params.alphabar * params.Di / params.cpbar, maxLevel_, All );

      shearHeatingTermCoeff->interpolate( params.Di / ( params.Ra * params.cpbar ), maxLevel_, All );
      shearHeatingTermCoeff->multElementwise( { *shearHeatingTermCoeff, *rhoInvP2 }, maxLevel_, All );

      constEnergyCoeff->interpolate( 0.0, maxLevel_, All );
      surfTempCoeff->interpolate( 0.0, maxLevel_, All );

      TRef->interpolate( TRefFunc, maxLevel_, All );

      invGravityField->component( 0u ).interpolate( 0.0, maxLevel_, All );
      invGravityField->component( 1u ).interpolate( 1.0, maxLevel_, All );

      transportTALAOp->setInvGravity( { std::make_shared< std::function< real_t( const Point3D& ) > >( invGravityX ),
                                        std::make_shared< std::function< real_t( const Point3D& ) > >( invGravityY ) } );
      transportTALAOp->setShearHeatingCoeff( shearHeatingTermCoeff );

      transportTALAOp->setDiffusivityCoeff(
          std::make_shared< std::function< real_t( const Point3D& ) > >( diffusivityCoeffFunc ) );
      transportTALAOp->setAdiabaticCoeff( std::make_shared< std::function< real_t( const Point3D& ) > >( adiabaticCoeffFunc ) );
      transportTALAOp->setConstEnergyCoeff(
          std::make_shared< std::function< real_t( const Point3D& ) > >( constEnergyCoeffFunc ) );
      transportTALAOp->setSurfTempCoeff( std::make_shared< std::function< real_t( const Point3D& ) > >( surfTempCoeffFunc ) );
      transportTALAOp->setReferenceTemperature( std::make_shared< std::function< real_t( const Point3D& ) > >( TRefFunc ) );

      transportTALAOp->setTALADict( {
          { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, !params.MMOC },
          { terraneo::TransportOperatorTermKey::DIFFUSION_TERM, true },
          { terraneo::TransportOperatorTermKey::SHEAR_HEATING_TERM, params.shearHeating },
          { terraneo::TransportOperatorTermKey::ADIABATIC_HEATING_TERM, params.adiabaticHeating },
          { terraneo::TransportOperatorTermKey::INTERNAL_HEATING_TERM, false },
          { terraneo::TransportOperatorTermKey::SUPG_STABILISATION, params.SUPG },
      } );

      transportTALAOp->initializeOperators();

      projectionOperator = std::make_shared< P2ProjectNormalOperator >( storage_, minLevel_, maxLevel_, normalsFS );

      viscP2->interpolate( 1.0, maxLevel_, All );
      /// [StokesOperator]
      stokesOperator = std::make_shared< P2P1StokesOperator >( storage_, minLevel_, maxLevel_, *viscP2 );
      /// [StokesOperator]

      params.cflMax = mainConf.getParameter< real_t >( "cflMax" );

      params.verbose     = mainConf.getParameter< bool >( "verbose" );
      stokesMinresSolver = std::make_shared< MinResSolver< P2P1StokesOperator > >(
          storage_, minLevel_, maxLevel_, params.minresIter, params.minresRelTol );
      stokesMinresSolver->setPrintInfo( params.verbose );

#if defined( HYTEG_BUILD_WITH_PETSC )
      stokesDirectSolver = std::make_shared< PETScLUSolver< P2P1StokesOperator > >( storage_, maxLevel_ );
#endif

      transportGmresSolver = std::make_shared< GMRESSolver< P2TransportTimesteppingOperator > >(
          storage_, minLevel_, maxLevel_, params.gmresIter, params.gmresTol, params.gmresTol );
      transportGmresSolver->setPrintInfo( params.verbose );

      transportTALAGmresSolver = std::make_shared< GMRESSolver< P2TransportTALAOperator > >(
          storage_, minLevel_, maxLevel_, params.gmresIter, params.gmresTol, params.gmresTol );
      transportTALAGmresSolver->setPrintInfo( true );

      std::string outputFilename = mainConf.getParameter< std::string >( "outputFilename" );
      std::string outputPath     = mainConf.getParameter< std::string >( "outputPath" );

      if ( params.SUPG )
      {
         outputFilename.append( "_SUPG" );
      }
      else if ( params.MMOC )
      {
         outputFilename.append( "_MMOC" );
      }
      else
      {
         outputFilename.append( "_noMMOC" );
      }

      vtkOutput = std::make_shared< VTKOutput >( outputPath, outputFilename, storage );

      adiosXmlConfig = mainConf.getParameter< std::string >( "adiosXmlConfig" );
      adios2Output   = std::make_shared< AdiosWriter >( outputPath, outputFilename, adiosXmlConfig, storage );

      adios2Output->add( *u );
      adios2Output->add( *TDev );
      adios2Output->add( *TRef );
      adios2Output->add( *TRefDev );
      adios2Output->add( *diffusionTermCoeff );
      adios2Output->add( *rhoP2 );

      storeCheckpoint     = mainConf.getParameter< bool >( "storeCheckpoint" );
      startFromCheckpoint = mainConf.getParameter< bool >( "startFromCheckpoint" );

      storeCheckpointFreq = mainConf.getParameter< uint_t >( "storeCheckpointFreq" );

      cpPath     = mainConf.getParameter< std::string >( "cpPath" );
      cpFilename = mainConf.getParameter< std::string >( "cpFilename" );

      cpStartFilename = mainConf.getParameter< std::string >( "cpStartFilename" );
   }

   void solveU();
   void solveT();
   void step();
   void solve();
   void writeVTK( uint_t timestep = 0 ) { adios2Output->write( maxLevel, timestep ); }

 private:
   const walberla::Config::BlockHandle& mainConf;

   std::shared_ptr< PrimitiveStorage > storage;
   uint_t                              minLevel, maxLevel;

   std::shared_ptr< P2Function< real_t > > TDev, TDevPrev, TDevInt, TRef, TRefDev, TRhs, TRes, rhoP2, rhoInvP2, zero, viscP2,
       Ttemp;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > u, uRes, uPrev, uPrevIter, uRhs, uRhsStrong, uTemp;
   std::shared_ptr< P2VectorFunction< real_t > >       gradRhoByRhoP2;
   std::shared_ptr< P2VectorFunction< real_t > >       invGravityField;

   std::shared_ptr< P2Function< real_t > > diffusionTermCoeff, adiabaticTermCoeff, shearHeatingTermCoeff, constEnergyCoeff,
       surfTempCoeff;

   std::shared_ptr< P2P1StokesOperator > stokesOperator;

   std::shared_ptr< P2ProjectNormalOperator > projectionOperator;

   P2ElementwiseBlendingVectorMassOperator vecMassOperator;
   P2ElementwiseBlendingMassOperator       massOperator;
   P1ElementwiseBlendingMassOperator       massOperatorP1;

   std::shared_ptr< P2ToP1ElementwiseKMass > gradRhoByRhoX;
   std::shared_ptr< P2ToP1ElementwiseKMass > gradRhoByRhoY;

   std::shared_ptr< P2TransportTimesteppingOperator > transportOp;

   std::shared_ptr< P2TransportTALAOperator > transportTALAOp;

   MMOCTransport< P2Function< real_t > > transport;

   // Solvers
#if defined( HYTEG_BUILD_WITH_PETSC )
   std::shared_ptr< PETScLUSolver< P2P1StokesOperator > > stokesDirectSolver;
#endif
   std::shared_ptr< MinResSolver< P2P1StokesOperator > > stokesMinresSolver;

   std::shared_ptr< GMRESSolver< P2TransportTimesteppingOperator > > transportGmresSolver;
   std::shared_ptr< GMRESSolver< P2TransportTALAOperator > >         transportTALAGmresSolver;

   ParameterContainer params;

   uint_t iTimeStep = 0U;

   real_t simulationTime = 0.0, endTime = 1.0;

   std::function< real_t( const Point3D& ) > tempIni, tempDevBC, rhoFunc, rhoInvFunc, TRefFunc;

   std::function< real_t( const Point3D& ) > diffusivityCoeffFunc;
   std::function< real_t( const Point3D& ) > adiabaticCoeffFunc;
   std::function< real_t( const Point3D& ) > shearHeatingCoeffFunc;
   std::function< real_t( const Point3D& ) > constEnergyCoeffFunc;
   std::function< real_t( const Point3D& ) > surfTempCoeffFunc;

   std::function< real_t( const Point3D& ) > invGravityX;
   std::function< real_t( const Point3D& ) > invGravityY;

   std::function< void( const Point3D&, Point3D& ) > normalsFS;

   // Output

   std::shared_ptr< VTKOutput > vtkOutput;

   uint_t      storeCheckpointFreq = 1000U;
   bool        storeCheckpoint = false, startFromCheckpoint = false;
   std::string cpFilename, cpPath, cpStartFilename, adiosXmlConfig;

   std::shared_ptr< AdiosWriter > adios2Output;

   std::shared_ptr< AdiosCheckpointExporter > adios2CheckpointExporter;
   std::shared_ptr< AdiosCheckpointImporter > adios2CheckpointImporter;
};

void TALASimulation::solveU()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING STOKES SOLVER" ) );

   uRhsStrong->uvw().component( 0U ).interpolate( 0.0, maxLevel, All );
   uRhsStrong->uvw().component( 1U ).interpolate( params.Ra * params.alphabar, maxLevel, All );

   TRefDev->assign( { 1.0, -1.0 }, { *TDev, *TRef }, maxLevel, All );
   uRhsStrong->uvw().component( 0 ).multElementwise( { uRhsStrong->uvw().component( 0 ), *TRefDev }, maxLevel, All );
   uRhsStrong->uvw().component( 1 ).multElementwise( { uRhsStrong->uvw().component( 1 ), *TRefDev }, maxLevel, All );

   vecMassOperator.apply( uRhsStrong->uvw(), uRhs->uvw(), maxLevel, All );

   /// [StokesFrozenVelocity]
   gradRhoByRhoY->apply( u->uvw().component( 1U ), uRhs->p(), maxLevel, All );
   uRhs->p().assign( { -1.0 }, { uRhs->p() }, maxLevel, All );
   /// [StokesFrozenVelocity]

   projectionOperator->project( *uRhs, maxLevel, FreeslipBoundary );

   u->uvw().interpolate( 0.0, maxLevel, DirichletBoundary );

   bool directSolverFlag = mainConf.getParameter< bool >( "directSolver" );

   /// [StokesSolverLambdaFunction]
   if ( directSolverFlag )
   {
#if defined( HYTEG_BUILD_WITH_PETSC )
      stokesDirectSolver->solve( *stokesOperator, *u, *uRhs, maxLevel );
#else
      WALBERLA_ABORT( "Direct solver requested but PETSc not compiled!" );
#endif
   }
   else
   {
      stokesMinresSolver->solve( *stokesOperator, *u, *uRhs, maxLevel );
   }
   /// [StokesSolverLambdaFunction]

   vertexdof::projectMean( u->p(), maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STOKES SOLVER DONE!" ) );
}

/// [TransportSolverLambdaFunction]
void TALASimulation::solveT()
{
   transportTALAOp->calculateTimestep( params.cflMax );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING TRANSPORT SOLVER with dt = %2.6e", transportTALAOp->timestep ) );

   /// [MMOCForTransport]
   if ( params.MMOC )
   {
      transportTALAOp->stepMMOC( maxLevel );
   }
   /// [MMOCForTransport]

   TDev->interpolate( tempDevBC, maxLevel, DirichletBoundary );
   transportTALAOp->applyRHS( *TRhs, maxLevel, All );

   TDev->interpolate( tempDevBC, maxLevel, DirichletBoundary );

   transportTALAGmresSolver->solve( *transportTALAOp, *TDev, *TRhs, maxLevel );

   transportTALAOp->incrementTimestep();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER DONE!" ) );
}
/// [TransportSolverLambdaFunction]

void TALASimulation::step()
{
   real_t vMax = u->uvw().getMaxComponentMagnitude( maxLevel, All );
   real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   real_t Pe = hMax * vMax / ( 4 * params.k_ );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Peclet number = %f", Pe ) );

   solveT();

   solveU();

   TDevPrev->assign( { 1.0 }, { *TDev }, maxLevel, All );

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );
}

void TALASimulation::solve()
{
   TDev->interpolate( tempIni, maxLevel, Inner | NeumannBoundary );
   TDev->interpolate( tempDevBC, maxLevel, DirichletBoundary );

   TDevPrev->assign( { 1.0 }, { *TDev }, maxLevel, All );

   /// [InitialSolveOrCheckpoint]
   if ( startFromCheckpoint )
   {
      adios2CheckpointImporter = std::make_shared< AdiosCheckpointImporter >( cpPath, cpStartFilename, adiosXmlConfig );

      adios2CheckpointImporter->restoreFunction( u->uvw() );
      adios2CheckpointImporter->restoreFunction( *TDev );
   }
   else
   {
      solveU();
   }
   /// [InitialSolveOrCheckpoint]

   TRefDev->assign( { 1.0, 1.0 }, { *TRef, *TDev }, maxLevel, All );

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );

   writeVTK( iTimeStep );

   iTimeStep++;

   while ( simulationTime < endTime && iTimeStep < params.maxTimeSteps )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nStarting step %d at time = %f!\n", iTimeStep, simulationTime ) );

      step();

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Step done!" ) );

      iTimeStep++;

      simulationTime += transportTALAOp->timestep;

      TRefDev->assign( { -1.0, 1.0 }, { *TRef, *TDev }, maxLevel, All );

      if ( iTimeStep % params.vtkWriteFrequency == 0 )
      {
         writeVTK( iTimeStep );
      }

      if ( iTimeStep % params.nsCalcFreq == 0 )
      {
         real_t deltaT           = TRefDev->getMaxDoFValue( maxLevel ) - TRefDev->getMinDoFValue( maxLevel );
         real_t nusseltNumberTop = nusseltcalc::calculateNusseltNumber2D( *TDev, maxLevel, 0.01, 1e-6, 101 );
         real_t velocityRMSValue = nusseltcalc::velocityRMS( *u, *uTemp, massOperator, 1, 1, maxLevel );

         uint_t nVelocityDoFs    = numberOfGlobalDoFs( *u, maxLevel );
         uint_t nTemperatureDoFs = numberOfGlobalDoFs( *TDev, maxLevel );

         WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
             "nusseltNumberTop = %4.7e\nnusseltNumberBot = not_calculated\ndeltaT = %4.7e\nvelocityRMS = %4.7e\nnVelocityDoFs = %u\nnTemperatureDoFs = %u\nnDoFs = %u",
             nusseltNumberTop,
             deltaT,
             velocityRMSValue,
             nVelocityDoFs,
             nTemperatureDoFs,
             nVelocityDoFs + nTemperatureDoFs ) );
      }

      if ( iTimeStep % storeCheckpointFreq == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Storing Checkpoint!" );

         adios2CheckpointExporter = std::make_shared< AdiosCheckpointExporter >( adiosXmlConfig );

         adios2CheckpointExporter->registerFunction( u->uvw(), minLevel, maxLevel );
         adios2CheckpointExporter->registerFunction( *TDev, minLevel, maxLevel );

         adios2CheckpointExporter->storeCheckpoint( cpPath, cpFilename );
      }
   }
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#if defined( HYTEG_BUILD_WITH_PETSC )
   PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./GB.02_KingCompressible.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   const uint_t nx = mainConf.getParameter< uint_t >( "nx" );
   const uint_t ny = mainConf.getParameter< uint_t >( "ny" );

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   auto meshInfo = hyteg::MeshInfo::meshRectangle( Point2D( 0.0, 0.0 ), Point2D( 1.0, 1.0 ), hyteg::MeshInfo::CRISS, nx, ny );

   /// [SetupStorageAndMarkings]
   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Bottom, bottomMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Left, leftMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Right, rightMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Top, topMarker );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( BoundaryMarkers::Corners, cornersMarker );

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1 );
   /// [SetupStorageAndMarkings]

   uint_t nMacroFaces      = storage->getNumberOfGlobalFaces();
   uint_t nMacroPrimitives = storage->getNumberOfGlobalPrimitives();

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "\n\nMacroFaces = %d, nMacroPrimitives = %d\n", nMacroFaces, nMacroPrimitives ) );

   TALASimulation simulation( mainConf, storage, minLevel, maxLevel );

   simulation.solve();

   return 0;
}