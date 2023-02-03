/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Dominik Thoennes, Nils Kohl, Marcus Mohr.
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
#include <core/math/Constants.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dgfunctionspace_old/DG0P1UpwindOperator.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1IntegrateDG.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   timingTree->start( "Global" );
   std::string meshFileName = "../data/meshes/flow_around_cylinder.msh";

   real_t viscosity = 1e-4;

   bool   neumann  = true;
   uint_t minLevel = 2;
   uint_t maxLevel = 3;

   real_t time              = 0.0;
   real_t inflowBuildupTime = 0.0;
   // real_t endTime           = 5.0;
   real_t endTime         = 0.1;
   uint_t iter            = 0;
   uint_t max_cg_iter     = 50;
   uint_t outerIterations = 2;

   hyteg::MeshInfo              meshInfo = hyteg::MeshInfo::fromGmshFile( meshFileName );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

#ifdef WALBERLA_BUILD_WITH_PARMETIS
   loadbalancing::distributed::parmetis( *storage );
#endif

   //  const real_t minimalEdgeLength = hyteg::MeshQuality::getMinimalEdgeLength(storage, maxLevel);
   //  real_t dt = 0.025 * minimalEdgeLength;
   real_t dt      = 5e-6;
   real_t dt_plot = 0.005;

   uint_t plotModulo = uint_c( std::ceil( dt_plot / dt ) );

   WALBERLA_LOG_INFO_ON_ROOT( "dt = " << dt );

   std::function< real_t( const hyteg::Point3D& ) > bc_x = [&time, &inflowBuildupTime]( const hyteg::Point3D& x ) {
      const real_t U_m = 5.0;

      if ( x[0] < 1e-8 )
      {
         real_t velocity = 4.0 * U_m * x[1] * ( 0.41 - x[1] ) / ( 0.41 * 0.41 );
         real_t damping;

         if ( time < inflowBuildupTime )
         {
            damping = 0.5 * ( 1.0 + std::cos( walberla::math::pi * ( time / inflowBuildupTime - 1.0 ) ) );
         }
         else
         {
            damping = 1.0;
         }

         return damping * velocity;
      }
      else
      {
         return real_c( 0.0 );
      }
   };

   std::function< real_t( const hyteg::Point3D& ) > bc_y = []( const hyteg::Point3D& ) { return 0.0; };

   std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return 0.0; };
   std::function< real_t( const hyteg::Point3D& ) > one  = []( const hyteg::Point3D& ) { return 1.0; };

   // hyteg::P1Function< real_t > u( "u", storage, minLevel, maxLevel );
   // hyteg::P1Function< real_t > v( "v", storage, minLevel, maxLevel );
   hyteg::P1VectorFunction< real_t > velocity( "velocity", storage, minLevel, maxLevel );

   hyteg::P1Function< real_t > p( "p", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > p_rhs( "p_rhs", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > p_res( "p_res", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > tmp( "tmp", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > tmp2( "tmp2", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > res( "res", storage, minLevel, maxLevel );
   hyteg::P1Function< real_t > ones( "ones", storage, minLevel, maxLevel );

   auto u_dg = std::make_shared< hyteg::DGFunction_old< real_t > >( "u_dg", storage, minLevel, maxLevel );
   auto v_dg = std::make_shared< hyteg::DGFunction_old< real_t > >( "v_dg", storage, minLevel, maxLevel );

   auto u_dg_old = std::make_shared< hyteg::DGFunction_old< real_t > >( "u_dg", storage, minLevel, maxLevel );
   auto v_dg_old = std::make_shared< hyteg::DGFunction_old< real_t > >( "v_dg", storage, minLevel, maxLevel );

   hyteg::P1ConstantLaplaceOperator A( storage, minLevel, maxLevel );
   hyteg::P1ConstantLaplaceOperator Ascaled( storage, minLevel, maxLevel );

   // Scale Laplace operator with viscosity
   Ascaled.scale( viscosity );

   hyteg::P1DivxOperator          div_x( storage, minLevel, maxLevel );
   hyteg::P1DivyOperator          div_y( storage, minLevel, maxLevel );
   hyteg::P1DivTxOperator         divT_x( storage, minLevel, maxLevel );
   hyteg::P1DivTyOperator         divT_y( storage, minLevel, maxLevel );
   hyteg::P1LumpedInvMassOperator invDiagMass( storage, minLevel, maxLevel );

   // std::array< hyteg::P1Function< real_t >, 2 > velocity{ u,v };
   // hyteg::DGUpwindOperator< hyteg::P1Function< real_t > > N( storage, velocity, minLevel, maxLevel );
   hyteg::DG0P1UpwindOperator N( storage, velocity, minLevel, maxLevel );

   typedef hyteg::CGSolver< hyteg::P1ConstantLaplaceOperator > CoarseSolver;
   auto coarseLaplaceSolver  = std::make_shared< CoarseSolver >( storage, minLevel, minLevel, max_cg_iter );
   auto smoother             = std::make_shared< hyteg::GaussSeidelSmoother< P1ConstantLaplaceOperator > >();
   auto restrictionOperator  = std::make_shared< P1toP1LinearRestriction >();
   auto prolongationOperator = std::make_shared< P1toP1LinearProlongation >();

   typedef GeometricMultigridSolver< hyteg::P1ConstantLaplaceOperator > LaplaceSover;
   LaplaceSover                                                         laplaceSolver(
       storage, smoother, coarseLaplaceSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

   // u.interpolate( bc_x, maxLevel, hyteg::DirichletBoundary );
   // v.interpolate( bc_y, maxLevel, hyteg::DirichletBoundary );
   velocity.interpolate( { bc_x, bc_y }, hyteg::DirichletBoundary );
   p.interpolate( zero, maxLevel - 1, hyteg::NeumannBoundary );
   ones.interpolate( real_t( 1 ), maxLevel, hyteg::All );

   u_dg->projectP1( velocity[0], maxLevel, hyteg::All );
   v_dg->projectP1( velocity[1], maxLevel, hyteg::All );

   hyteg::VTKOutput vtkOutput( "../output", "test", storage, plotModulo );
   vtkOutput.add( velocity );
   vtkOutput.add( p );
   vtkOutput.write( maxLevel, iter );
   ++iter;

   while ( time < endTime )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "time = " << time );
      time += dt;
      velocity[0].interpolate( bc_x, maxLevel, hyteg::DirichletBoundary );

      u_dg_old->projectP1( velocity[0], maxLevel, hyteg::All );
      v_dg_old->projectP1( velocity[1], maxLevel, hyteg::All );

      N.apply( *u_dg_old, *u_dg, maxLevel, hyteg::All, Replace );
      N.apply( *v_dg_old, *v_dg, maxLevel, hyteg::All, Replace );

      // Predict u
      // tmp.integrateDG( *u_dg, ones, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      P1IntegrateDG( *u_dg, ones, tmp, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

      Ascaled.apply( velocity[0], tmp, maxLevel, hyteg::Inner | hyteg::NeumannBoundary, Add );
      invDiagMass.apply( tmp, tmp2, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      velocity[0].add( { -dt }, { tmp2 }, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

      // Predict v
      // tmp.integrateDG( *v_dg, ones, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      P1IntegrateDG( *v_dg, ones, tmp, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      Ascaled.apply( velocity[1], tmp, maxLevel, hyteg::Inner | hyteg::NeumannBoundary, Add );
      invDiagMass.apply( tmp, tmp2, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      velocity[1].add( { -dt }, { tmp2 }, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

      // Solve p
      p.interpolate( zero, maxLevel - 1, hyteg::NeumannBoundary );
      div_x.apply( velocity[0], p_rhs, maxLevel, hyteg::Inner | hyteg::DirichletBoundary, Replace );
      div_y.apply( velocity[1], p_rhs, maxLevel, hyteg::Inner | hyteg::DirichletBoundary, Add );

      restrictionOperator->restrict( p_rhs, maxLevel, hyteg::Inner | hyteg::DirichletBoundary );

      if ( !neumann )
      {
         hyteg::vertexdof::projectMean( p_rhs, maxLevel - 1 );
      }

      for ( uint_t outer = 0; outer < outerIterations; ++outer )
      {
         laplaceSolver.solve( A, p, p_rhs, maxLevel - 1 );
      }

      if ( !neumann )
      {
         hyteg::vertexdof::projectMean( p, maxLevel - 1 );
      }

      prolongationOperator->prolongate( p, maxLevel - 1, hyteg::Inner | hyteg::DirichletBoundary );

      // Correct velocities
      divT_x.apply( p, tmp, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      invDiagMass.apply( tmp, tmp2, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      velocity[0].add( { -1.0 }, { tmp2 }, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

      divT_y.apply( p, tmp, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      invDiagMass.apply( tmp, tmp2, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );
      velocity[1].add( { -1.0 }, { tmp2 }, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

      vtkOutput.write( maxLevel, iter );
      ++iter;

      u_dg_old.swap( u_dg );
      v_dg_old.swap( v_dg );
   }

   timingTree->stop( "Global" );
   auto reduced_tt = timingTree->getReduced();
   WALBERLA_LOG_INFO_ON_ROOT( reduced_tt );

   return 0;
}
