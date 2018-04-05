#include <cmath>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultiGrid.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

//using namespace hhg;

int main( int argc, char** argv )
{
   /// [Create Environment]
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   /// [Create Environment]

   /// [Get Parameters]
   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "../data/param/GMG_P1.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t      minLevel         = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t      maxLevel         = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t      max_outer_iter   = parameters.getParameter< uint_t >( "max_outer_iter" );
   const uint_t      max_coarse_iter  = parameters.getParameter< uint_t >( "max_coarse_iter" );
   const real_t      mg_tolerance     = parameters.getParameter< real_t >( "mg_tolerance" );
   const real_t      coarse_tolerance = parameters.getParameter< real_t >( "coarse_tolerance" );
   const std::string meshFile         = parameters.getParameter< std::string >( "mesh" );
   /// [Get Parameters]

   /// [Create Primitive Storage]
   std::shared_ptr< hhg::PrimitiveStorage > storage = hhg::PrimitiveStorage::createFromGmshFile( meshFile );
   /// [Create Primitive Storage]

   /// [Create Function Spaces]
   hhg::P1Function< real_t > residuum( "residuum", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > rightHandSide( "rightHandSide", storage, minLevel, maxLevel );
   hhg::P1Function< real_t > function( "function", storage, minLevel, maxLevel );
   /// [Create Function Spaces]

   /// [Create Functions]
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
   };
   /// [Create Functions]

   /// [Set Boundary Conditions]
   function.interpolate( boundaryConditions, maxLevel, hhg::DirichletBoundary );
   /// [Set Boundary Conditions]

   /// [Setup Solvers]
   typedef hhg::CGSolver< hhg::P1Function< real_t >, hhg::P1LaplaceOperator >                       CoarseSolver;
   typedef hhg::GMultigridSolver< hhg::P1Function< real_t >, hhg::P1LaplaceOperator, CoarseSolver > GMG;

   hhg::P1LaplaceOperator laplaceOperator( storage, minLevel, maxLevel );
   auto                   coarseSolver = std::make_shared< CoarseSolver >( storage, minLevel, maxLevel );
   GMG                    multiGridSolver( storage, coarseSolver, minLevel, maxLevel );
   /// [Setup Solvers]

   for( uint_t i = 0; i < max_outer_iter; ++i )
   {
      multiGridSolver.solve( laplaceOperator,
                             function,
                             rightHandSide,
                             residuum,
                             maxLevel,
                             coarse_tolerance,
                             max_coarse_iter,
                             hhg::Inner,
                             GMG::CycleType::VCYCLE );
   }

   if (parameters.getParameter<bool>("vtkOutput")) {
      VTKOutput vtkOutput("../output", "gmg_P2_h_refinement");
      vtkOutput.add(&function);
      vtkOutput.add(&residuum);
      vtkOutput.add(&rightHandSide);
      vtkOutput.write(maxLevel);
   }
}