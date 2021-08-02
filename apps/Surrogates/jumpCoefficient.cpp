/*
 * Copyright (c) 2021 Benjamin Mann.
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

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator_new.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1ScaledSurrogateOperator.hpp"
#include "hyteg/p1functionspace/P1SurrogateOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator_new.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2SurrogateOperator.hpp"
#include "hyteg/p2functionspace/P2VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

using function = std::function< real_t( const hyteg::Point3D& ) >;

enum OpType
{
   VARIABLE         = 0,
   SURROGATE        = 1,
   SCALED_SURROGATE = 2
};

// print stencils for all points
template < class A_t >
void printStencils( std::shared_ptr< PrimitiveStorage > storage, const A_t& A )
{
   if ( storage->hasGlobalCells() )
   { // 3D version
      for ( auto& cellID : storage->getCellIDs() )
      {
         auto cellptr  = storage->getCell( cellID );
         auto stencils = A.computeStencilsForCell( *cellptr );

         for ( auto& [node, stencil] : stencils )
         {
            std::cout << node[0] << " " << node[1] << " " << node[2];

            for ( auto& [idx, weight] : stencil )
            {
               std::cout << " " << weight;
            }

            std::cout << "\n";
         }

         std::cout << "\n\n";
      }
   }
   else
   { // 2D version
      for ( auto& faceID : storage->getFaceIDs() )
      {
         auto faceptr  = storage->getFace( faceID );
         auto stencils = A.computeStencilsForFace( *faceptr );

         for ( auto& [node, stencil] : stencils )
         {
            std::cout << node[0] << " " << node[1];

            for ( auto& weight : stencil )
            {
               std::cout << " " << weight;
            }

            std::cout << "\n";
         }

         std::cout << "\n\n";
      }
   }
}

// print ||A-Aq||_max
template < class A_t >
void printError( const A_t& A, uint_t l_min, uint_t l_max )
{
   double maxerr = 0;
   double minerr = 1e100;
   std::cout << "lvl | ||A-Aq||_cell for all cells\n";

   for ( uint_t lvl = l_min; lvl <= l_max; ++lvl )
   {
      auto err = A.computeSurrogateError( lvl );
      std::cout << lvl << " | [";
      for ( auto& el : err )
      {
         std::cout << el << " ";
         if ( lvl == l_max )
         {
            if ( el > maxerr )
               maxerr = el;
            if ( el < minerr )
               minerr = el;
         }
      }
      std::cout << "]\n";
   }

   auto err = A.computeSurrogateError( l_max );
   std::cout << "refine:\n\'";
   for ( uint_t i = 0; i < err.size(); ++i )
   {
      if ( err[i] > ( minerr + maxerr ) / 2 )
      {
         std::cout << i;
         if ( i < err.size() - 1 )
            std::cout << " ";
      }
   }
   std::cout << "\'\n";
}

// solve Au=b(f)
template < class A_t, class M_t, class FE_t, class R_t, class P_t >
void solve( std::shared_ptr< PrimitiveStorage > storage,
            const A_t&                          A,
            const function&                     u_anal,
            const function&                     f,
            uint_t                              l_min,
            uint_t                              l_max,
            uint_t                              mg_iter,
            uint_t                              cg_iter,
            real_t                              mg_tol,
            real_t                              cg_tol,
            bool                                vtk )
{
   // define functions and operators
   FE_t r( "r", storage, l_min, l_max );
   FE_t b( "b", storage, l_min, l_max );
   FE_t u( "u", storage, l_min, l_max );
   FE_t Au( "Au", storage, l_min, l_max );
   FE_t u_exact( "u_exact", storage, l_min, l_max );
   FE_t err( "err", storage, l_min, l_max );
   FE_t tmp( "tmp", storage, l_min, l_max );

   // rhs
   M_t M( storage, l_min, l_max );
   tmp.interpolate( f, l_max );
   M.apply( tmp, b, l_max, hyteg::All );

   // exact solution
   u_exact.interpolate( u_anal, l_max );

   // boundary
   u.interpolate( u_anal, l_max, hyteg::DirichletBoundary );
   u.interpolate( []( const hyteg::Point3D& ) { return 0.0; }, l_max, hyteg::Inner );

   // define solver
   auto cg           = std::make_shared< hyteg::CGSolver< A_t > >( storage, l_min, l_min, cg_iter, cg_tol );
   auto restriction  = std::make_shared< R_t >();
   auto prolongation = std::make_shared< P_t >();
   auto smoother     = std::make_shared< hyteg::GaussSeidelSmoother< A_t > >();

   GeometricMultigridSolver< A_t > gmg( storage, smoother, cg, restriction, prolongation, l_min, l_max, 3, 3 );

   // solve iteratively
   uint_t       iter       = 0;
   real_t       res        = 0, res_old, discr_l2_err;
   real_t       vCycleTime = 0, solveTime = 0;
   real_t       averageConvergenceRate = 0;
   const uint_t convergenceStartIter   = 3;

   WALBERLA_LOG_INFO_ON_ROOT( "Starting V cycles" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%10s|%10s|%10s|%10s", "iter", "res", "conv", "L2-error", "Cycle-Time" ) );

   while ( 1 )
   {
      // compute error
      err.assign( { 1.0, -1.0 }, { u, u_exact }, l_max );
      M.apply( err, tmp, l_max, hyteg::All, Replace );
      discr_l2_err = std::sqrt( err.dotGlobal( tmp, l_max ) );

      // compute residual
      res_old = res;
      A.apply( u, Au, l_max, hyteg::Inner, Replace );
      r.assign( { 1.0, -1.0 }, { b, Au }, l_max, hyteg::Inner );
      M.apply( r, tmp, l_max, hyteg::All, Replace );
      res = std::sqrt( r.dotGlobal( tmp, l_max, hyteg::Inner ) );

      // compute convergence factor
      real_t convRate = res / res_old;
      if ( iter >= convergenceStartIter )
      {
         averageConvergenceRate += convRate;
      }

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e", iter, res, convRate, discr_l2_err, vCycleTime ) );

      // stopping criterion
      if ( ++iter > mg_iter || ( iter > convergenceStartIter && convRate > 0.999 ) )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Ending multigrid without reaching desired tolerance!" );
         break;
      }
      if ( res < mg_tol )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Multigrid converged!" );
         break;
      }

      // solve
      auto start = walberla::timing::getWcTime();
      gmg.solve( A, u, b, l_max );
      auto end   = walberla::timing::getWcTime();
      vCycleTime = end - start;
      solveTime += vCycleTime;
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Solve time " << std::defaultfloat << solveTime );
   WALBERLA_LOG_INFO_ON_ROOT( "Avg. convergence rate: " << std::scientific
                                                        << averageConvergenceRate / real_c( iter - convergenceStartIter ) );
   WALBERLA_LOG_INFO_ON_ROOT( "L^2 error: " << std::scientific << discr_l2_err );
   // WALBERLA_LOG_INFO_ON_ROOT("DoFs: " << (uint_t) npoints);

   if ( vtk )
   {
      auto        Nel  = ( storage->hasGlobalCells() ) ? storage->getNumberOfGlobalCells() : storage->getNumberOfGlobalFaces();
      std::string name = "JumpCoefficient_Nel" + std::to_string( Nel );

      hyteg::VTKOutput vtkOutput( "output", name, storage );
      vtkOutput.setVTKDataFormat( VTKOutput::VTK_DATA_FORMAT::BINARY );
      vtkOutput.add( u );
      vtkOutput.add( err );
      vtkOutput.add( r );
      vtkOutput.add( u_exact );
      vtkOutput.add( b );
      vtkOutput.write( l_max, 0 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "=======================================================" );
}

int main( int argc, char** argv )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg;

   if ( argc == 1 )
   {
      walberla::shared_ptr< walberla::config::Config > cfg_( new walberla::config::Config );
      cfg_->readParameterFile( "../../hyteg/data/param/jumpCoefficient.prm" );
      cfg = cfg_;
   }
   else
   {
      cfg = walberla::config::create( argc, argv );
   }

   // read parameter file
   WALBERLA_LOG_INFO_ON_ROOT( "config = " << *cfg );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t minLevel = parameters.getParameter< uint_t >( "level_h_coarse" );
   const uint_t maxLevel = parameters.getParameter< uint_t >( "level_h_fine" );

   const int    discretization = parameters.getParameter< int >( "discretization" );
   const int    tmp            = parameters.getParameter< int >( "operatorType" );
   const OpType opType         = OpType( tmp );
   const uint_t polyDegree     = parameters.getParameter< uint_t >( "polyDegree" );

   const real_t alpha    = parameters.getParameter< real_t >( "alpha" );
   const real_t x_jump_0 = parameters.getParameter< real_t >( "x_jump_0" );
   const real_t x_jump_1 = parameters.getParameter< real_t >( "x_jump_1" );
   const real_t k_min    = parameters.getParameter< real_t >( "k_min" );
   const real_t k_max    = parameters.getParameter< real_t >( "k_max" );

   const uint_t max_outer_iter   = parameters.getParameter< uint_t >( "max_outer_iter" );
   const uint_t max_cg_iter      = parameters.getParameter< uint_t >( "max_cg_iter" );
   const real_t mg_tolerance     = parameters.getParameter< real_t >( "mg_tolerance" );
   const real_t coarse_tolerance = parameters.getParameter< real_t >( "coarse_tolerance" );

   const bool vtk            = parameters.getParameter< bool >( "vtkOutput" );
   const bool print_stencils = parameters.getParameter< bool >( "print_stencils" );

   const std::string msh = parameters.getParameter< std::string >( "mesh" );

   // domain
   MeshInfo meshInfo = MeshInfo::fromGmshFile( msh );
   // MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { 0.0, 0.0 } ), Point2D( { 1.0, 1.0 } ), MeshInfo::CRISS, 1, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::stringstream ss;
   setupStorage.toStream( ss, true );
   WALBERLA_LOG_INFO_ON_ROOT( ss.str() );

   // correct boundary for mesh with hanging nodes
   // constexpr real_t minH = 1e-10;
   // std::function< bool( const hyteg::Point3D& ) > leftBnd = [=]( const hyteg::Point3D& x ) {
   //    return x[0] < minH;
   // };
   // std::function< bool( const hyteg::Point3D& ) > rightBnd = [=]( const hyteg::Point3D& x ) {
   //    return x[0] > (1-minH);
   // };
   // std::function< bool( const hyteg::Point3D& ) > downBnd = [=]( const hyteg::Point3D& x ) {
   //    return x[1] < minH;
   // };
   // std::function< bool( const hyteg::Point3D& ) > upBnd = [=]( const hyteg::Point3D& x ) {
   //    return x[1] > (1-minH);
   // };
   // setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, leftBnd, true );
   // setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, rightBnd, true );
   // setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, downBnd, true );
   // setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, upBnd, true );

   hyteg::loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // std::cerr << "\nfaces: ";
   // for (auto& idx : storage->getFaces())
   // {
   //    std::cerr << " " << idx.second->getMeshBoundaryFlag();
   // }
   // std::cerr << "\nedges: ";
   // for (auto& idx : storage->getEdges())
   // {
   //    std::cerr << " " << idx.second->getMeshBoundaryFlag();
   // }
   // std::cerr << "\nvertices: ";
   // for (auto& idx : storage->getVertices())
   // {
   //    std::cerr << " " << idx.second->getMeshBoundaryFlag();
   // }
   // std::cerr << "\n";

   // coefficient function
   function k = [=]( const hyteg::Point3D& x ) {
      // jump at { (x,y) | x + (x_jump_0 - x_jump_1)y = x_jump_0 }
      real_t x_jump_y = x_jump_0 * ( 1.0 - x[1] ) + x_jump_1 * x[1];
      return ( x[0] < x_jump_y ) ? k_min : k_max;
   };

   // analytic solution
   real_t du_dx_min = 1.0 / k_min;
   real_t du_dx_max = 1.0 / k_max;
   real_t du_dy_min = ( x_jump_0 - x_jump_1 ) / k_min;
   real_t du_dy_max = ( x_jump_0 - x_jump_1 ) / k_max;
   real_t u_jump    = du_dx_min * x_jump_0;

   function u = [=]( const hyteg::Point3D& x ) {
      // kink at { (x,y) | x + (x_jump_0 - x_jump_1)y = x_jump_0 }
      real_t x_jump_y = x_jump_0 * ( 1.0 - x[1] ) + x_jump_1 * x[1];
      if ( x[0] < x_jump_y )
      {
         return du_dx_min * x[0] + du_dy_min * x[1];
      }
      else
      {
         return u_jump + du_dx_max * ( x[0] - x_jump_0 ) + du_dy_max * x[1];
      }
   };

   // rhs
   function f = []( const hyteg::Point3D& ) { return 0.0; };

   // test with smooth solution and coeff
   if ( alpha > 0.0 )
   {
      // real_t phi              = 6;
      real_t scaling_factor_u = 1;

      // u = [=]( const hyteg::Point3D& x ) { return scaling_factor_u * sin( phi * M_PI * x[0] ) * sinh( M_PI * x[1] ); };
      // k = [=]( const hyteg::Point3D& x ) { return tanh( alpha * ( x[0] - 0.5 ) ) + 2; };
      // f = [=]( const hyteg::Point3D& x ) {
      //    real_t t0 = tanh( alpha * ( x[0] - 0.5 ) );
      //    real_t t1 = phi * alpha * ( t0 * t0 - 1 ) * cos( phi * M_PI * x[0] );
      //    real_t t2 = ( phi * phi - 1 ) * M_PI * ( t0 + 2 ) * sin( phi * M_PI * x[0] );
      //    return scaling_factor_u * M_PI * ( t1 + t2 ) * sinh( M_PI * x[1] );
      // };

      u = [=]( const hyteg::Point3D& x ) {
         return scaling_factor_u * ( -2 * alpha * x[0] + log( 2 * exp( alpha * ( x[0] - 0.5 ) ) + 1 ) ) / alpha;
      };
      k = [=]( const hyteg::Point3D& x ) {
         return ( 2 * exp( alpha * ( x[0] - 0.5 ) ) + 1 ) / ( exp( alpha * ( x[0] - 0.5 ) ) + 1 );
      };
   }

   if ( discretization == 1 )
   {
      using M_t = P1ConstantMassOperator_new;
      using FE  = hyteg::P1Function< real_t >;
      using R_t = hyteg::P1toP1LinearRestriction;
      using P_t = hyteg::P1toP1LinearProlongation;

      using A_form = forms::p1_div_k_grad_affine_q1;
      A_form form( k, k );

      P1VariableOperator_new< A_form > A1( storage, minLevel, maxLevel, form );
      P1SurrogateOperator< A_form >    A1q( storage, minLevel, maxLevel, form );
      A1q.interpolateStencils( polyDegree, maxLevel );
      P1ScaledSurrogateOperator< A_form > A1qs( storage, minLevel, maxLevel, form );
      A1qs.interpolateStencils( polyDegree, maxLevel );

      switch ( opType )
      {
      case VARIABLE:
         solve< P1VariableOperator_new< A_form >, M_t, FE, R_t, P_t >(
             storage, A1, u, f, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk );
         if ( print_stencils )
            printStencils< P1VariableOperator_new< A_form > >( storage, A1 );
         break;
      case SURROGATE:
         solve< P1SurrogateOperator< A_form >, M_t, FE, R_t, P_t >(
             storage, A1q, u, f, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk );
         printError< P1SurrogateOperator< A_form > >( A1q, minLevel, maxLevel );
         if ( print_stencils )
            printStencils< P1SurrogateOperator< A_form > >( storage, A1q );
         break;
      case SCALED_SURROGATE:
         solve< P1ScaledSurrogateOperator< A_form >, M_t, FE, R_t, P_t >(
             storage, A1qs, u, f, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk );
         if ( print_stencils )
            printStencils< P1ScaledSurrogateOperator< A_form > >( storage, A1qs );
         break;
      default:
         WALBERLA_ABORT( "The desired Operator Type is not supported!" );
      }
   }
   else if ( discretization == 2 )
   {
      using M_t = P2ConstantMassOperator;
      using FE  = hyteg::P2Function< real_t >;
      using R_t = hyteg::P2toP2QuadraticRestriction;
      using P_t = hyteg::P2toP2QuadraticProlongation;

      P2Form_divKgrad::callback = k;

      P2divKgradOperator          A2( storage, minLevel, maxLevel );
      P2SurrogateDivKgradOperator A2q( storage, minLevel, maxLevel, maxLevel, polyDegree );

      switch ( opType )
      {
      case VARIABLE:
         solve< P2divKgradOperator, M_t, FE, R_t, P_t >(
             storage, A2, u, f, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk );
         break;
      case SURROGATE:
         solve< P2SurrogateDivKgradOperator, M_t, FE, R_t, P_t >(
             storage, A2q, u, f, minLevel, maxLevel, max_outer_iter, max_cg_iter, mg_tolerance, coarse_tolerance, vtk );
         break;
         // case SCALED_SURROGATE:
         //    break;

      default:
         WALBERLA_ABORT( "The desired Operator Type is not supported!" );
      }
   }
   else
   {
      WALBERLA_ABORT( "The desired discretization order is not supported!" );
   }
}
