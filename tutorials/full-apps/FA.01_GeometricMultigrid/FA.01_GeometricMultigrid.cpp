/*
 * Copyright (c) 2017-2024 Dominik Bartuschat, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

/**
 * \page FA.01_GeometricMultigrid Tutorial FA.01 - Geometric multigrid
 *
 * \dontinclude tutorials/full-apps/FA.01_GeometricMultigrid/FA.01_GeometricMultigrid.cpp
 *
 * \brief In this tutorial we will set up a complete app using P1 elements which will perform geometric
 * multigrid to solve the Laplace equation.
 *
 * \section FA01-FullAppP1GMG-setup Setup MPI
 *
 * At we first we create a walberla Environment which handles mpi init and set the communicator for all
 * prozesses to WorldComm
 *
 * \snippet{trimleft} this Create Environment
 *
 * \section FA01-FullAppP1GMG-parameters Set Parameters
 *
 * One way to set the parameters for the simulation is to use a parameter file.
 *
 * We need to create a walberla Config and read a parameter file. In this case the file is in the same
 * directiory as the executable.
 *
 * Now we can extract blocks from the config and access the desired parameters in the block.
 *
 * See walberla documentation for more details:
 * http://www.walberla.net/doxygen/classwalberla_1_1config_1_1Config.html
 *
 * \snippet{trimleft} this Get Parameters
 *
 *
 * \section FA01-FullAppP1GMG-storage Primitive Storage
 *
 * In this step we create a fully distributed PrimitiveStorage from a mesh file
 *
 * This creation can be split into more steps for more flexibility like using other load
 * balancing techniques
 *
 * \snippet{trimleft} this Primitive Storage
 *
 * \section FA01-FullAppP1GMG-functionSpaces Create P1 Function Spaces
 *
 * Here we allocate the necessary P1 Functions.
 *
 * In the constructor we need to specify:
 *
 * - a name which is needed for the vtk output
 * - the storage on which the data should be allocated
 * - the smallest level and
 * - the highest level which should be allocated
 *
 * Be aware that the highest level is included
 *
 * \snippet{trimleft} this Function Spaces
 *
 * \section FA01-FullAppP1GMG-boundaries Create functions for boundary conditions
 *
 * To set the boundary conditions we create a function using a lamba function. The functions needs to return
 * a real_t and take and hyteg::Point3D as an argument.
 *
 * In this case we transfer the Cartesian to polar coordinates and use the angle to apply some sine pattern.
 *
 * The interpolate function is used to evaluate the function at all desired points.
 * We need to call the function with:
 *
 * - the function which should be evaluated
 * - the desired level
 * - on which points the function should be evaluated
 *
 * In our case we only evaluate the function on Dirichlet boundary conditions. Other options would be
 * for example hyteg::Inner for all points no on a boundary or hyteg::NeumannBondary for points which are
 * specified as Neumann boundaries
 *
 * \snippet{trimleft} this Boundary Conditions
 *
 * \section FA01-FullAppP1GMG-solver Solver
 *
 * Now we will setup the solvers.
 *
 * As a coarse grid solver we use a CG solver which works on a P1Function and uses a P1LaplaceOperator.
 * They are specified as template parameters. In the constructor we need to specify the storage, the min
 * and the max level.
 *
 * For the geometric multigrid solver we need to specify the function (P1Function), the operator
 * (P1LaplaceOperator), the restriction operator (P1toP1LinearRestriction), the prolongation operator
 * (P1toP1LinearProlongation) and the coarse grid solver (CoarseSolver)
 * as template arguments.
 * The storage, the coarse grid solver and the levels are needed as parameters.
 *
 * Furthermore we need to create an instance of the P1LaplaceOperator
 *
 * \snippet{trimleft} this Solvers
 *
 * \section FA01-FullAppP1GMG-multigrid Perform Multigrid
 *
 * In order to perform the multigrid iterations we call the function solve.
 *
 * The parameters for the function call are:
 * - operator
 * - function to work on
 * - the right-hand side of the equation
 * - an additional function which will be used as a residual
 * - highest level in the multigrid hierarchy
 * - tolerance for the coarse-grid solver
 * - maximal iterations for the coarse-grid solver
 * - on which points to work
 * - cycle type (V or W)
 *
 * Be aware that it is not assured that the residual functin will contain the up to date residual after
 * the function call. Also the multigrid will always solve the problem down to the smallest level of
 * the function
 *
 * When calling the solve function with CycleType::VCycle one cycle is performed per call.
 * We use a for loop to perform the desired number of V cycles
 *
 * \snippet{trimleft} this Multigrid
 *
 * \section FA01-FullAppP1GMG-Calculate Residual
 *
 * To validate the result we calculate the residual by first writing the result of multiplying the
 * Laplace operator with the function into the laplaceTimesFunction function.
 *
 * Next we subtract this from the right-hand side and calculate the Euclidean norm.
 *
 * We can use the walberla logging functionality to print out the norm of the residual
 *
 * See: http://www.walberla.net/doxygen/classwalberla_1_1logging_1_1Logging.html
 *
 * \snippet{trimleft} this Residual
 *
 * \section FA01-FullAppP1GMG-output Write VTK Output
 *
 * If set in the parameter file the results are written onto disc in the vtk file format.
 *
 * First we create an VTKOutput object which needs to path to write the files to (".") and a name for the
 * output ("FullAppP1GMG").
 *
 * Next we register all functions that should be displayed by using the add function.
 *
 * By specifiying a certain level as paramter in the write function, the desired level will be written to
 * disk
 *
 * \snippet{trimleft} this Output
 *
 * \section FA01-FullAppP1GMG-fullApp Full Application
 * \include tutorials/full-apps/FA.01_GeometricMultigrid/FA.01_GeometricMultigrid.cpp
 *
 */

int main( int argc, char** argv )
{
   /// [Create Environment]
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   /// [Create Environment]

   /// [Get Parameters]
   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "./FA.01_GeometricMultigrid.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );
   parameters.listParameters();

   const uint_t      minLevel         = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel         = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t      max_outer_iter   = parameters.getParameter< uint_t >( "max_outer_iter" );
   const uint_t      max_coarse_iter  = parameters.getParameter< uint_t >( "max_coarse_iter" );
   const real_t      coarse_tolerance = parameters.getParameter< real_t >( "coarse_tolerance" );
   const std::string meshFile         = parameters.getParameter< std::string >( "mesh" );
   /// [Get Parameters]

   /// [Primitive Storage]
   std::shared_ptr< hyteg::PrimitiveStorage > storage = hyteg::PrimitiveStorage::createFromGmshFile( hyteg::prependHyTeGMeshDir( meshFile ) );
   /// [Primitive Storage]

   /// [Function Spaces]
   hyteg::P1Function< real_t > residual( "residual", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > rightHandSide( "rightHandSide", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > function( "function", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > laplaceTimesFunction( "laplaceTimesFunction", storage, minLevel, maxLevel );
   /// [Function Spaces]

   /// [Boundary Conditions]
   std::function< real_t( const hyteg::Point3D& ) > boundaryConditions = []( const hyteg::Point3D& x ) {
      real_t radius = sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t angle  = std::atan2( x[1], x[0] );
      if ( radius < 1.1 )
      {
         return std::sin( 10 * angle );
      }
      else if ( radius > 1.9 )
      {
         return std::sin( 5 * angle );
      }
      WALBERLA_ABORT( "point is not on the boundary" );
   };

   function.interpolate( boundaryConditions, maxLevel, hyteg::DirichletBoundary );
   /// [Boundary Conditions]

   /// [Solvers]
   auto smoother         = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1ConstantLaplaceOperator > >();
   auto coarseGridSolver = std::make_shared< hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator > >(
       storage, minLevel, minLevel, max_coarse_iter, real_c(0), coarse_tolerance );
   auto restrictionOperator  = std::make_shared< hyteg::P1toP1LinearRestriction<> >();
   auto prolongationOperator = std::make_shared< hyteg::P1toP1LinearProlongation<> >();

   auto multiGridSolver = hyteg::GeometricMultigridSolver< hyteg::P1ConstantLaplaceOperator >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

   hyteg::P1ConstantLaplaceOperator laplaceOperator( storage, minLevel, maxLevel );
   /// [Solvers]

   /// [Multigrid]
   for ( uint_t i = 0; i < max_outer_iter; ++i )
   {
      multiGridSolver.solve( laplaceOperator, function, rightHandSide, maxLevel );
   }
   /// [Multigrid]

   /// [Residual]
   laplaceOperator.apply( function, laplaceTimesFunction, maxLevel, hyteg::Inner );
   residual.assign( { 1.0, -1.0 }, { rightHandSide, laplaceTimesFunction }, maxLevel, hyteg::Inner );
   real_t residualEuclideanNorm = std::sqrt( residual.dotGlobal( residual, maxLevel, hyteg::Inner ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Euclidean norm of residual: " << residualEuclideanNorm )
   /// [Residual]

   /// [Output]
   if ( parameters.getParameter< bool >( "vtkOutput" ) )
   {
      hyteg::VTKOutput vtkOutput( ".", "FullAppP1GMG", storage );
      vtkOutput.add( function );
      vtkOutput.add( residual );
      vtkOutput.add( rightHandSide );
      vtkOutput.write( maxLevel );
   }
   /// [Output]
}
