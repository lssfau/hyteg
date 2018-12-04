#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

/**
 * \page 05_FullAppP1GMG Creating a full app using P1 elements
 *
 * \dontinclude tutorials/05_FullAppP1GMG/P1GMG.cpp
 *
 * \brief In this tutorial we will set up a complete app using P1 elements which will perform geometric
 * multigrid to solve the Laplace equation.
 *
 * \section setup Setup MPI
 *
 * At we first we create a walberla Environment which handles mpi init and set the communicator for all
 * prozesses to WorldComm
 *
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Create Environment
 *
 * \section parameters Set Parameters
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
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Get Parameters
 *
 *
 * \section storage Primitive Storage
 *
 * In this step we create a fully distributed PrimitiveStorage from a mesh file
 *
 * This creation can be split into more steps for more flexibility like using other load
 * balancing techniques
 *
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Primitive Storage
 *
 * \section functionSpaces Create P1 Function Spaces
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
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Function Spaces
 *
 * \section boundaries Create functions for boundary conditions
 *
 * To set the boundary conditions we create a function using a lamba function. The functions needs to return
 * a real_t and take and hhg::Point3D as an argument.
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
 * for example hhg::Inner for all points no on a boundary or hhg::NeumannBondary for points which are
 * specified as Neumann boundaries
 *
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Boundary Conditions
 *
 * \section solver Solver
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
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Solvers
 *
 * \section multigrid Perform Multigrid
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
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Multigrid
 *
 * \section Calculate Residual
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
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Residual
 *
 * \section output Write VTK Output
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
 * \snippet tutorials/05_FullAppP1GMG/P1GMG.cpp Output
 *
 * \section fullApp Full Application
 * \include tutorials/05_FullAppP1GMG/P1GMG.cpp
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
   cfg->readParameterFile( "./P1GMG.prm" );
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
   std::shared_ptr< hhg::PrimitiveStorage > storage = hhg::PrimitiveStorage::createFromGmshFile( meshFile );
   /// [Primitive Storage]

   /// [Function Spaces]
   hhg::P1Function< real_t > residual( "residual", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > rightHandSide( "rightHandSide", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > function( "function", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > laplaceTimesFunction( "laplaceTimesFunction", storage, minLevel, maxLevel );
   /// [Function Spaces]

   /// [Boundary Conditions]
   std::function< real_t( const hhg::Point3D& ) > boundaryConditions = []( const hhg::Point3D& x ) {
      real_t radius = sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t angle  = std::atan2( x[1], x[0] );
      if( radius < 1.1 )
      {
         return std::sin( 10 * angle );
      } else if( radius > 1.9 )
      {
         return std::sin( 5 * angle );
      }
      WALBERLA_ABORT("point is not on the boundary");
   };

   function.interpolate( boundaryConditions, maxLevel, hhg::DirichletBoundary );
   /// [Boundary Conditions]

   /// [Solvers]
   typedef hhg::P1toP1LinearRestriction RestrictionOperator;
   typedef hhg::P1toP1LinearProlongation ProlongationOperator;
   typedef hhg::CGSolver< hhg::P1Function< real_t >, hhg::P1ConstantLaplaceOperator > CoarseSolver;

   RestrictionOperator restrictionOperator;
   ProlongationOperator prolongationOperator;
   auto coarseSolver = std::make_shared< CoarseSolver >( storage, minLevel, maxLevel );

   typedef hhg::GMultigridSolver< hhg::P1Function< real_t >, hhg::P1ConstantLaplaceOperator, CoarseSolver, RestrictionOperator, ProlongationOperator > GMG;
   GMG multiGridSolver( storage, coarseSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

   hhg::P1ConstantLaplaceOperator laplaceOperator( storage, minLevel, maxLevel );
   /// [Solvers]

   /// [Multigrid]
   for( uint_t i = 0; i < max_outer_iter; ++i )
   {
      multiGridSolver.solve( laplaceOperator,
                             function,
                             rightHandSide,
                             residual,
                             maxLevel,
                             coarse_tolerance,
                             max_coarse_iter,
                             hhg::Inner,
                             GMG::CycleType::VCYCLE );
   }
   /// [Multigrid]

   /// [Residual]
   laplaceOperator.apply( function, laplaceTimesFunction, maxLevel, hhg::Inner );
   residual.assign( {1.0, -1.0}, {&rightHandSide, &laplaceTimesFunction}, maxLevel, hhg::Inner );
   real_t residualEuclideanNorm = std::sqrt( residual.dotGlobal( residual, maxLevel, hhg::Inner ) );
   WALBERLA_LOG_INFO_ON_ROOT( "Euclidean norm of residual: " << residualEuclideanNorm )
   /// [Residual]

   /// [Output]
   if( parameters.getParameter< bool >( "vtkOutput" ) )
   {
      hhg::VTKOutput vtkOutput(".", "FullAppP1GMG", storage);
      vtkOutput.add( &function );
      vtkOutput.add( &residual );
      vtkOutput.add( &rightHandSide );
      vtkOutput.write( maxLevel );
   }
   /// [Output]
}
