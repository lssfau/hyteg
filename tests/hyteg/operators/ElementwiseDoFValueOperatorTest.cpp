/*
* Copyright (c) 2025 Andreas Burkhart
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
#include <iostream>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/geometry/SphericalCoordsMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2FullViscousTDependentOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "elementwise_dof_value_operator/generated/p1_k_mass_visc_centroid_blending_q4_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_adiabatic_heat_u_p2_supg_blending_q5_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_advection_u_p2_blending_q5_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_advection_u_p2_supg_blending_q4_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_diffusion_u_p2_supg_blending_q3_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_mass_u_p2_supg_blending_q5_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_shear_heat_T_p2_dep_eta_blending_q6_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_shear_heat_T_p2_dep_eta_supg_blending_q5_ElementwiseOperator.hpp"
#include "elementwise_dof_value_operator/generated/p2_to_p1_grad_rho_rho_p2_blending_q5_ElementwiseOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;
using namespace hyteg;

int main( int argc, char** argv )
{
   // Create walberla & MPI environment
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   if ( sizeof( real_t ) < 8 )
   {
      // The tolerances in this test are designed for at least double precision.
      WALBERLA_LOG_INFO_ON_ROOT( "Single precision or lower detected. Aborting test." );

      return EXIT_SUCCESS;
   }

   // This test requires a large amount of DoFs / high number of refinements, hence in debug mode
   // we reduce the refinement level. This is not optimal but some compromise has to be made here.

   {
      // Adiabatic Heating SUPG Operator

      // Extract the required parameters
      uint_t minLevel  = 4;
      uint_t maxLevel  = 4;
      real_t tolerance = 2e-4;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 2;
         maxLevel  = 2;
         tolerance = 3e-2;
      }

      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/unitsquare_with_circular_hole.msh",
                                              "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el_offcenter.msh",
                                              "../../../data/meshes/3D/pyramid_4el_offcenter.msh" };

      for ( size_t i = 0; i < 2; i++ )
      {
         if ( i )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "--- With blending ---" );
         }

         for ( auto s : Meshes2D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               PolarCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, minLevel, maxLevel );
            P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
            P2Function< real_t >       res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta_Di_rho_alpha = []( const hyteg::Point3D& x,
                                                                                              real_t                u_abs ) {
               WALBERLA_UNUSED( u_abs );
               WALBERLA_UNUSED( x );
               return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) );
            };

            std::function< real_t( const hyteg::Point3D& ) > gX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0];
            };

            std::function< real_t( const hyteg::Point3D& ) > gY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[1];
            };

            std::function< real_t( const hyteg::Point3D& ) > gZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 0;
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 ) - 5;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[0], 2 ) + x[0] * x[1] + 2;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 );
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 0;
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * std::pow( x[1], 3 ) * std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) ) *
                      ( std::pow( x[1], 2 ) - 5 ) * ( x[0] * ( std::pow( x[0], 2 ) + x[0] * x[1] + 2 ) + std::pow( x[1], 3 ) );
            };

            u.interpolate( { initX, initY }, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_adiabatic_heat_u_p2_supg_blending_q5_ElementwiseOperator HeatOp( storage,
                                                                                minLevel,
                                                                                maxLevel,
                                                                                u.component( 0 ),
                                                                                u.component( 1 ),
                                                                                u.component( 1 ),
                                                                                gX,
                                                                                gY,
                                                                                gZ,
                                                                                delta_Di_rho_alpha );
            P2ElementwiseBlendingMassOperator                           P2Mass( storage, minLevel, maxLevel );

            HeatOp.apply( T, res, maxLevel, All, Replace );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }

         for ( auto s : Meshes3D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               SphericalCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, minLevel, maxLevel );
            P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
            P2Function< real_t >       res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta_Di_rho_alpha = []( const hyteg::Point3D& x,
                                                                                              real_t                u_abs ) {
               WALBERLA_UNUSED( u_abs );
               WALBERLA_UNUSED( x );
               return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) * std::exp( -x[2] );
            };

            std::function< real_t( const hyteg::Point3D& ) > gX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0];
            };

            std::function< real_t( const hyteg::Point3D& ) > gY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[1];
            };

            std::function< real_t( const hyteg::Point3D& ) > gZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[2];
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * std::pow( x[2], 2 ) + 93;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] + 2 * x[1] - x[2];
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0] * x[1] + std::pow( x[2], 2 );
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -4 * x[2] * ( 2 * std::pow( x[2], 2 ) + 93 ) * ( x[0] * x[1] - std::pow( x[2], 2 ) ) *
                      std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                      ( x[0] * ( x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12 ) + x[1] * ( x[0] * x[1] + 2 * x[1] - x[2] ) -
                        x[2] * ( x[0] * x[1] - std::pow( x[2], 2 ) ) ) *
                      std::exp( -x[2] );
            };

            u.interpolate( { initX, initY, initZ }, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_adiabatic_heat_u_p2_supg_blending_q5_ElementwiseOperator HeatOp( storage,
                                                                                minLevel,
                                                                                maxLevel,
                                                                                u.component( 0 ),
                                                                                u.component( 1 ),
                                                                                u.component( 2 ),
                                                                                gX,
                                                                                gY,
                                                                                gZ,
                                                                                delta_Di_rho_alpha );
            P2ElementwiseBlendingMassOperator                           P2Mass( storage, minLevel, maxLevel );

            HeatOp.apply( T, res, maxLevel, All, Replace );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }
      }
   }

   {
      // Adiabatic Heating Operator

      // Extract the required parameters
      uint_t minLevel  = 4;
      uint_t maxLevel  = 4;
      real_t tolerance = 1e-4;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 2;
         maxLevel  = 2;
         tolerance = 3e-3;
      }

      // No Scaling
      std::function< real_t( const Point3D& ) > adiabaticHeatingScaling = [=]( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return real_c( 1.0 );
      };

      // Gravity
      std::function< real_t( const Point3D& ) > gX = [&]( const Point3D& x ) { return -x[0] / x.norm(); };
      std::function< real_t( const Point3D& ) > gY = [&]( const Point3D& x ) { return -x[1] / x.norm(); };
      std::function< real_t( const Point3D& ) > gZ = [&]( const Point3D& x ) { return -x[2] / x.norm(); };

      //Meshes
      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el.msh" };

      for ( auto s : Meshes2D )
      {
         // Init setup storage
         MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > T_( "T", storage, minLevel, maxLevel );
         P2Function< real_t > O_( "O", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > T = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] + std::pow( x[1], 2 ) - real_c( 5 );
         };

         std::function< real_t( const hyteg::Point3D& ) > O = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::pow( x[0], 2 ) + x[0] * x[1] - 3 * x[1] + 3;
         };

         std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::pow( x[0], 2 ) + x[0] * x[1] + real_c( 2 );
         };

         std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] + std::pow( x[1], 2 );
         };

         u.interpolate( { uX, uY }, maxLevel, hyteg::All );
         T_.interpolate( T, maxLevel, hyteg::All );
         O_.interpolate( O, maxLevel, hyteg::All );

         // operators
         auto AdiaOp = std::make_shared< p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator >( storage,
                                                                                                   minLevel,
                                                                                                   maxLevel,
                                                                                                   u.component( 0 ),
                                                                                                   u.component( 1 ),
                                                                                                   u.component( 0 ),
                                                                                                   gX,
                                                                                                   gY,
                                                                                                   gZ,
                                                                                                   adiabaticHeatingScaling );

         AdiaOp->apply( T_, res, maxLevel, All, Replace );

         real_t integVal0    = O_.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( -20.7952031306044 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      for ( auto s : Meshes3D )
      {
         // Init setup storage
         MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > T_( "T", storage, minLevel, maxLevel );
         P2Function< real_t > O_( "O", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > T = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] - x[1] * x[2] + real_c( 9 );
         };

         std::function< real_t( const hyteg::Point3D& ) > O = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[2] + x[1] * x[2] + x[1] - 3 * x[2] + 4;
         };

         std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                   ( x[0] * x[1] + x[2] - real_c( 12 ) );
         };

         std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                   ( x[0] * x[2] + real_c( 2 ) * x[1] - x[2] );
         };

         std::function< real_t( const hyteg::Point3D& ) > uZ = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return ( -x[0] * x[1] + std::pow( x[2], 2 ) ) *
                   std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) );
         };

         u.interpolate( { uX, uY, uZ }, maxLevel, hyteg::All );
         T_.interpolate( T, maxLevel, hyteg::All );
         O_.interpolate( O, maxLevel, hyteg::All );

         // operators
         auto AdiaOp = std::make_shared< p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator >( storage,
                                                                                                   minLevel,
                                                                                                   maxLevel,
                                                                                                   u.component( 0 ),
                                                                                                   u.component( 1 ),
                                                                                                   u.component( 2 ),
                                                                                                   gX,
                                                                                                   gY,
                                                                                                   gZ,
                                                                                                   adiabaticHeatingScaling );

         AdiaOp->apply( T_, res, maxLevel, All, Replace );

         real_t integVal0    = O_.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( -160.500694444444 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "Blending Test: " );
      {
         // Init setup storage
         real_t rCMB     = real_c( 1.0 );
         real_t rSurface = real_c( 2.0 );

         MeshInfo              meshInfo = MeshInfo::meshAnnulus( rCMB, rSurface, MeshInfo::meshFlavour::CRISS, 25, 4 );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         AnnulusMap::setMap( setupStorage );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > T_( "T", storage, minLevel, maxLevel );
         P2Function< real_t > O_( "O", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > T = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] + std::pow( x[1], 2 ) - real_c( 5 );
         };

         std::function< real_t( const hyteg::Point3D& ) > O = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::pow( x[0], 2 ) + x[0] * x[1] - 3 * x[1] + 3;
         };

         std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::pow( x[0], 2 ) + x[0] * x[1] + real_c( 2 );
         };

         std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] + std::pow( x[1], 2 );
         };

         u.interpolate( { uX, uY }, maxLevel, hyteg::All );
         T_.interpolate( T, maxLevel, hyteg::All );
         O_.interpolate( O, maxLevel, hyteg::All );

         // operators
         auto AdiaOp = std::make_shared< p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator >( storage,
                                                                                                   minLevel,
                                                                                                   maxLevel,
                                                                                                   u.component( 0 ),
                                                                                                   u.component( 1 ),
                                                                                                   u.component( 0 ),
                                                                                                   gX,
                                                                                                   gY,
                                                                                                   gZ,
                                                                                                   adiabaticHeatingScaling );

         AdiaOp->apply( T_, res, maxLevel, All, Replace );

         real_t integVal0    = O_.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( 298.989860831646 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      {
         // Init setup storage
         real_t rCMB     = real_c( 1.0 );
         real_t rSurface = real_c( 2.0 );

         // create the spherical shell mesh
         MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 3, 2, rCMB, rSurface );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         IcosahedralShellMap::setMap( setupStorage );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > T_( "T", storage, minLevel, maxLevel );
         P2Function< real_t > O_( "O", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > T = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] - x[1] * x[2] + real_c( 9 );
         };

         std::function< real_t( const hyteg::Point3D& ) > O = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[2] + x[1] * x[2] + x[1] - 3 * x[2] + 4;
         };

         std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                   ( x[0] * x[1] + x[2] - real_c( 12 ) );
         };

         std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                   ( x[0] * x[2] + real_c( 2 ) * x[1] - x[2] );
         };

         std::function< real_t( const hyteg::Point3D& ) > uZ = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return ( -x[0] * x[1] + std::pow( x[2], 2 ) ) *
                   std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) );
         };

         u.interpolate( { uX, uY, uZ }, maxLevel, hyteg::All );
         T_.interpolate( T, maxLevel, hyteg::All );
         O_.interpolate( O, maxLevel, hyteg::All );

         // operators
         auto AdiaOp = std::make_shared< p2_adiabatic_heat_u_p2_blending_q6_ElementwiseOperator >( storage,
                                                                                                   minLevel,
                                                                                                   maxLevel,
                                                                                                   u.component( 0 ),
                                                                                                   u.component( 1 ),
                                                                                                   u.component( 2 ),
                                                                                                   gX,
                                                                                                   gY,
                                                                                                   gZ,
                                                                                                   adiabaticHeatingScaling );

         AdiaOp->apply( T_, res, maxLevel, All, Replace );

         real_t integVal0    = O_.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( 613.159099500636 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }
   }

   {
      // Advection SUPG Operator
      // Extract the required parameters
      uint_t minLevel  = 5;
      uint_t maxLevel  = 5;
      real_t tolerance = 1e-4;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 3;
         maxLevel  = 3;
         tolerance = 1e-2;
      }

      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/unitsquare_with_circular_hole.msh",
                                              "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el_offcenter.msh",
                                              "../../../data/meshes/3D/pyramid_4el_offcenter.msh" };

      for ( size_t i = 0; i < 2; i++ )
      {
         if ( i )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "--- With blending ---" );
         }

         for ( auto s : Meshes2D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               PolarCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, minLevel, maxLevel );
            P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
            P2Function< real_t >       res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta_rho_C_p = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) ) * std::exp( x[0] * x[1] );
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] + std::pow( x[1], 2 ) - 5;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[0], 2 ) + x[0] * x[1] + 2;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 );
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 0;
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) ) *
                      std::pow( std::pow( x[0], 2 ) + x[0] * x[1] + 2 * std::pow( x[1], 3 ) + 2, 2 ) * std::exp( x[0] * x[1] );
            };

            u.interpolate( { initX, initY }, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_advection_u_p2_supg_blending_q4_ElementwiseOperator AdOp(
                storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 1 ), delta_rho_C_p );
            P2ElementwiseBlendingMassOperator P2Mass( storage, minLevel, maxLevel );

            AdOp.apply( T, res, maxLevel, All );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }

         for ( auto s : Meshes3D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               SphericalCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, minLevel, maxLevel );
            P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
            P2Function< real_t >       res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta_rho_C_p = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                      std::exp( x[0] * x[1] + x[2] );
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] - x[1] * x[2] + 2 * std::pow( x[2], 2 ) + 93;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] + 2 * x[1] - x[2];
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0] * x[1] + std::pow( x[2], 2 );
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                      std::pow( x[1] * ( x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12 ) +
                                    ( x[0] - x[2] ) * ( x[0] * x[1] + 2 * x[1] - x[2] ) +
                                    ( x[1] - 4 * x[2] ) * ( x[0] * x[1] - std::pow( x[2], 2 ) ),
                                2 ) *
                      std::exp( x[0] * x[1] + x[2] );
            };

            u.interpolate( { initX, initY, initZ }, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_advection_u_p2_supg_blending_q4_ElementwiseOperator AdOp(
                storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 2 ), delta_rho_C_p );
            P2ElementwiseBlendingMassOperator P2Mass( storage, minLevel, maxLevel );

            AdOp.apply( T, res, maxLevel, All );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }
      }
   }

   {
      // Advection Operator

      // Extract the required parameters
      uint_t minLevel  = 4;
      uint_t maxLevel  = 4;
      real_t tolerance = 1e-7;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 2;
         maxLevel  = 2;
         tolerance = 1e-5;
      }

      // No Scaling
      std::function< real_t( const Point3D& ) > advectionScaling = [=]( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return real_c( 1.0 );
      };

      // Meshes
      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el.msh" };

      for ( auto s : Meshes2D )
      {
         // Init setup storage
         MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > T_( "T", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > T = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[0] + x[1] * x[1] + x[0] + real_c( 2.0 );
         };

         std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] + real_c( 2.0 );
         };
         std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[1];
         };

         u.interpolate( { initX, initY }, maxLevel, hyteg::All );
         T_.interpolate( T, maxLevel, hyteg::All );

         O.interpolate( 1, maxLevel, hyteg::All );

         // operators

         auto AdvecOp = std::make_shared< p2_advection_u_p2_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 0 ), advectionScaling );

         AdvecOp->apply( T_, res, maxLevel, All, Replace );

         real_t integVal0    = O.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( 61 ) / real_c( 12 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      for ( auto s : Meshes3D )
      {
         // Init setup storage
         MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > T_( "T", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > T = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[0] + real_c( 2.0 );
         };

         std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] + real_c( 2.0 );
         };
         std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[1];
         };
         std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[2];
         };

         u.interpolate( { initX, initY, initZ }, maxLevel, hyteg::All );
         T_.interpolate( T, maxLevel, hyteg::All );

         O.interpolate( 1, maxLevel, hyteg::All );
         // operators
         auto AdvecOp = std::make_shared< p2_advection_u_p2_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 2 ), advectionScaling );

         AdvecOp->apply( T_, res, maxLevel, All, Replace );

         real_t integVal0    = O.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( 65 ) / real_c( 12 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "Blending Test: " );
      {
         // Init setup storage
         real_t rCMB     = 1.0;
         real_t rSurface = 2.0;

         MeshInfo              meshInfo = MeshInfo::meshAnnulus( rCMB, rSurface, MeshInfo::meshFlavour::CRISS, 25, 4 );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         AnnulusMap::setMap( setupStorage );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > T_( "T", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > T = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[0] + x[1] * x[1] + x[0] + real_c( 2.0 );
         };

         std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] + real_c( 2.0 );
         };
         std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[1];
         };

         u.interpolate( { initX, initY }, maxLevel, hyteg::All );
         T_.interpolate( T, maxLevel, hyteg::All );

         O.interpolate( 1, maxLevel, hyteg::All );

         // operators
         auto AdvecOp = std::make_shared< p2_advection_u_p2_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 0 ), advectionScaling );

         AdvecOp->apply( T_, res, maxLevel, All, Replace );

         real_t integVal0    = O.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( 6 ) * pi;

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      {
         // Init setup storage
         real_t rCMB     = 1.0;
         real_t rSurface = 2.0;

         // create the spherical shell mesh
         MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 3, 2, rCMB, rSurface );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         IcosahedralShellMap::setMap( setupStorage );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > T_( "T", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > T = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[0] + real_c( 2.0 );
         };

         std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] + real_c( 2.0 );
         };
         std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[1];
         };
         std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[2];
         };

         u.interpolate( { initX, initY, initZ }, maxLevel, hyteg::All );
         T_.interpolate( T, maxLevel, hyteg::All );

         O.interpolate( 1, maxLevel, hyteg::All );

         // operators
         auto AdvecOp = std::make_shared< p2_advection_u_p2_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 2 ), advectionScaling );

         AdvecOp->apply( T_, res, maxLevel, All, Replace );

         real_t integVal0    = O.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = ( real_c( 56 ) * pi ) / real_c( 3 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }
   }

   {
      // Diffusion SUPG Test

      // Extract the required parameters
      uint_t minLevel       = 5;
      uint_t maxLevel       = 5;
      real_t tolerance      = 5e-4;
      uint_t maxLevelShell  = maxLevel - 2;
      real_t toleranceShell = 8.5e-3;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel       = 3;
         maxLevel       = 3;
         maxLevelShell  = maxLevel - 1;
         tolerance      = 7e-3;
         toleranceShell = 6e-2;
      }

      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/unitsquare_with_circular_hole.msh",
                                              "../../../data/meshes/2D/quad_4el_offcenter.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el_offcenter.msh",
                                              "../../../data/meshes/3D/pyramid_4el_offcenter.msh" };

      for ( size_t i = 0; i < 2; i++ )
      {
         if ( i )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "--- With blending ---" );
         }

         for ( auto s : Meshes2D )
         {
            WALBERLA_LOG_INFO_ON_ROOT( s );
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               PolarCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, minLevel, maxLevel );
            P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
            P2Function< real_t >       res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return 1.0;
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 ) - 5;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[0], 2 ) + x[0] * x[1] + 2;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 );
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 0;
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -4 * std::pow( x[1], 3 );
            };

            u.interpolate( { initX, initY }, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_diffusion_u_p2_supg_blending_q3_ElementwiseOperator DiffOp(
                storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 1 ), delta );
            P2ElementwiseBlendingMassOperator P2Mass( storage, minLevel, maxLevel );

            DiffOp.apply( T, res, maxLevel, All );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }

         for ( auto s : Meshes3D )
         {
            WALBERLA_LOG_INFO_ON_ROOT( s );
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               SphericalCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, minLevel, maxLevel );
            P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
            P2Function< real_t >       res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return 1.0;
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * std::pow( x[2], 2 ) + 93;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] + 2 * x[1] - x[2];
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0] * x[1] + std::pow( x[2], 2 );
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 16 * x[2] * ( x[0] * x[1] - std::pow( x[2], 2 ) );
            };

            u.interpolate( { initX, initY, initZ }, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_diffusion_u_p2_supg_blending_q3_ElementwiseOperator DiffOp(
                storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 2 ), delta );
            P2ElementwiseBlendingMassOperator P2Mass( storage, minLevel, maxLevel );

            DiffOp.apply( T, res, maxLevel, All );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }

         {
            WALBERLA_LOG_INFO_ON_ROOT( "SphericalShell" );
            std::vector< real_t > layers = { 0.25, 0.625, 1.0 };

            // Init setup storage
            MeshInfo              meshInfo = hyteg::MeshInfo::meshSphericalShell( 5, layers );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               IcosahedralShellMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, maxLevelShell, maxLevelShell );
            P2Function< real_t >       O( "O", storage, maxLevelShell, maxLevelShell );
            P2VectorFunction< real_t > u( "u", storage, maxLevelShell, maxLevelShell );
            P2Function< real_t >       res( "res", storage, maxLevelShell, maxLevelShell );

            P2Function< real_t > Ctrl( "Ctrl", storage, maxLevelShell, maxLevelShell );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, maxLevelShell, maxLevelShell );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return 1.0;
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::exp( x[0] + x[2] );
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] + 2 * x[1] - x[2];
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0] * x[1] + std::pow( x[2], 2 );
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return ( 6 * x[0] - 4 * std::pow( x[2], 2 ) + 24 ) * std::exp( 2 * x[0] + 2 * x[2] );
            };

            u.interpolate( { initX, initY, initZ }, maxLevelShell, hyteg::All );

            T.interpolate( initT, maxLevelShell, hyteg::All );
            O.interpolate( 1, maxLevelShell, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevelShell, hyteg::All );

            // operators
            p2_diffusion_u_p2_supg_blending_q3_ElementwiseOperator DiffOp(
                storage, maxLevelShell, maxLevelShell, u.component( 0 ), u.component( 1 ), u.component( 2 ), delta );
            P2ElementwiseBlendingMassOperator P2Mass( storage, maxLevelShell, maxLevelShell );

            DiffOp.apply( T, res, maxLevelShell, All );
            P2Mass.apply( Ctrl, CtrlRes, maxLevelShell, All );

            real_t integVal0    = T.dotGlobal( res, maxLevelShell, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevelShell, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), toleranceShell );
         }
      }
   }

   {
      // Full Stokes with temperature dependent viscosity

      // Extract the required parameters
      uint_t minLevel  = 4;
      uint_t maxLevel  = 4;
      real_t tolerance = 1e-6;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 2;
         maxLevel  = 2;
         tolerance = 2e-4;
      }

      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/unitsquare_with_circular_hole.msh",
                                              "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el.msh", "../../../data/meshes/3D/pyramid_4el.msh" };

      for ( size_t i = 0; i < 2; i++ )
      {
         if ( i )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "--- With blending ---" );
         }

         for ( auto s : Meshes2D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               PolarCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, minLevel, maxLevel );
            P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 ) + x[2] - 5;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[0], 2 ) + x[0] * x[1] + 2;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 );
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) ) *
                      ( 5 * std::pow( x[0], 2 ) * std::pow( x[1], 2 ) - 25 * std::pow( x[0], 2 ) -
                        4 * x[0] * std::pow( x[1], 3 ) + 20 * x[0] * x[1] + std::pow( x[1], 4 ) - 5 * std::pow( x[1], 2 ) );
            };

            std::function< real_t( const hyteg::Point3D&, real_t ) > etaT = []( const hyteg::Point3D& x, real_t temp ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( temp );
               real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
               return radius * temp;
            };

            u.interpolate( { initX, initY }, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            P2ElementwiseBlendingFullViscousTDependentOperator FVOp( storage, minLevel, maxLevel, T, etaT );
            P2ElementwiseBlendingMassOperator                  P2Mass( storage, minLevel, maxLevel );

            FVOp.apply( u, res, maxLevel, All, Replace );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = u.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }

         for ( auto s : Meshes3D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               SphericalCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t >       T( "T", storage, minLevel, maxLevel );
            P2Function< real_t >       O( "O", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );
            P2VectorFunction< real_t > res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * std::pow( x[2], 2 ) + 93;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] + 2 * x[1] - x[2];
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0] * x[1] + std::pow( x[2], 2 );
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return ( 1.0 / 3.0 ) * std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                      ( 10 * std::pow( x[0], 2 ) + 2 * x[0] * x[1] - 8 * x[0] * x[2] + 34 * x[0] + 10 * std::pow( x[1], 2 ) -
                        20 * x[1] * x[2] - 32 * x[1] + 28 * std::pow( x[2], 2 ) + 8 * x[2] + 79 ) /
                      ( 2 * std::pow( x[2], 2 ) + 93 );
            };

            std::function< real_t( const hyteg::Point3D&, real_t ) > etaT = []( const hyteg::Point3D& x, real_t temp ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( temp );
               real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
               return radius / temp;
            };

            u.interpolate( { initX, initY, initZ }, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            P2ElementwiseBlendingFullViscousTDependentOperator FVOp( storage, minLevel, maxLevel, T, etaT );
            P2ElementwiseBlendingMassOperator                  P2Mass( storage, minLevel, maxLevel );

            FVOp.apply( u, res, maxLevel, All, Replace );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = u.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }
      }
   }

   {
      // GradRhoRho Test

      // Extract the required parameters
      uint_t minLevel  = 4;
      uint_t maxLevel  = 4;
      real_t tolerance = 1e-5;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 2;
         maxLevel  = 2;
         tolerance = 3e-4;
      }

      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el.msh" };

      for ( auto s : Meshes2D )
      {
         // Init setup storage
         MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P1Function< real_t >       O( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P1Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > rho_fem( "rho", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > rho = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::exp( x[0] + x[1] );
         };

         std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] + real_c( 2.0 );
         };
         std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[1];
         };

         u.interpolate( { initX, initY }, maxLevel, hyteg::All );
         rho_fem.interpolate( rho, maxLevel, hyteg::All );

         O.interpolate( 1, maxLevel, hyteg::All );

         // operators
         auto GradRhoRhoOp0_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_0_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );
         auto GradRhoRhoOp1_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_1_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );

         GradRhoRhoOp0_->apply( u.component( 0 ), res, maxLevel, All, Replace );
         GradRhoRhoOp1_->apply( u.component( 1 ), res, maxLevel, All, Add );

         real_t integVal0    = O.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( real_c( -31 ) / real_c( 12 ) );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      for ( auto s : Meshes3D )
      {
         // Init setup storage
         MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P1Function< real_t >       O( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P1Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > rho_fem( "rho", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > rho = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::exp( x[0] + x[1] + x[2] );
         };

         std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] + real_c( 2.0 );
         };
         std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[1];
         };
         std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[2];
         };

         u.interpolate( { initX, initY, initZ }, maxLevel, hyteg::All );
         rho_fem.interpolate( rho, maxLevel, hyteg::All );

         O.interpolate( 1, maxLevel, hyteg::All );
         // operators
         auto GradRhoRhoOp0_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_0_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );
         auto GradRhoRhoOp1_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_1_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );
         auto GradRhoRhoOp2_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_2_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );

         GradRhoRhoOp0_->apply( u.component( 0 ), res, maxLevel, All, Replace );
         GradRhoRhoOp1_->apply( u.component( 1 ), res, maxLevel, All, Add );
         GradRhoRhoOp2_->apply( u.component( 2 ), res, maxLevel, All, Add );

         real_t integVal0    = O.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( -17 ) / real_c( 6 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "Blending Test: " );
      {
         // Init setup storage
         real_t rCMB     = 1.0;
         real_t rSurface = 2.0;

         MeshInfo              meshInfo = MeshInfo::meshAnnulus( rCMB, rSurface, MeshInfo::meshFlavour::CRISS, 25, 4 );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         AnnulusMap::setMap( setupStorage );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P1Function< real_t >       O( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P1Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > rho_fem( "rho", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > rho = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::exp( x[0] + x[1] );
         };

         std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] + real_c( 2.0 );
         };
         std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[1];
         };

         u.interpolate( { initX, initY }, maxLevel, hyteg::All );
         rho_fem.interpolate( rho, maxLevel, hyteg::All );

         O.interpolate( 1, maxLevel, hyteg::All );

         // operators
         auto GradRhoRhoOp0_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_0_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );
         auto GradRhoRhoOp1_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_1_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );

         GradRhoRhoOp0_->apply( u.component( 0 ), res, maxLevel, All, Replace );
         GradRhoRhoOp1_->apply( u.component( 1 ), res, maxLevel, All, Add );

         real_t integVal0    = O.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( -39 ) / real_c( 4 ) * pi;

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      {
         // Init setup storage
         real_t rCMB     = 1.0;
         real_t rSurface = 2.0;

         // create the spherical shell mesh
         MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 3, 2, rCMB, rSurface );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         IcosahedralShellMap::setMap( setupStorage );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P1Function< real_t >       O( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P1Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > rho_fem( "rho", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D& ) > rho = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::exp( x[0] + x[1] + x[2] );
         };

         std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] + real_c( 2.0 );
         };
         std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[1];
         };
         std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[2];
         };

         u.interpolate( { initX, initY, initZ }, maxLevel, hyteg::All );
         rho_fem.interpolate( rho, maxLevel, hyteg::All );

         O.interpolate( 1, maxLevel, hyteg::All );

         // operators
         auto GradRhoRhoOp0_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_0_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );
         auto GradRhoRhoOp1_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_1_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );
         auto GradRhoRhoOp2_ = std::make_shared< p2_to_p1_grad_rho_rho_p2_2_blending_q5_ElementwiseOperator >(
             storage, minLevel, maxLevel, rho_fem );

         GradRhoRhoOp0_->apply( u.component( 0 ), res, maxLevel, All, Replace );
         GradRhoRhoOp1_->apply( u.component( 1 ), res, maxLevel, All, Add );
         GradRhoRhoOp2_->apply( u.component( 2 ), res, maxLevel, All, Add );

         real_t integVal0    = O.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = ( real_c( -404 ) * pi ) / real_c( 15 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }
   }

   {
      // Mass SUPG Operator

      // Extract the required parameters
      uint_t minLevel  = 4;
      uint_t maxLevel  = 4;
      real_t tolerance = 2e-5;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 2;
         maxLevel  = 2;
         tolerance = 5e-3;
      }

      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/unitsquare_with_circular_hole.msh",
                                              "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el.msh", "../../../data/meshes/3D/pyramid_4el.msh" };

      for ( size_t i = 0; i < 2; i++ )
      {
         if ( i )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "--- With blending ---" );
         }

         for ( auto s : Meshes2D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               PolarCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t > T( "T", storage, minLevel, maxLevel );
            P2Function< real_t > M( "M", storage, minLevel, maxLevel );
            P2Function< real_t > O( "O", storage, minLevel, maxLevel );
            P2Function< real_t > uX( "uX", storage, minLevel, maxLevel );
            P2Function< real_t > uY( "uY", storage, minLevel, maxLevel );
            P2Function< real_t > res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return 1.0;
            };

            std::function< real_t( const hyteg::Point3D& ) > initM = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] + x[1] + 4;
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 ) - 5;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[0], 2 ) + x[0] * x[1] + 2;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 );
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 0;
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * std::pow( x[1], 3 ) * ( x[0] * x[1] + x[1] + 4 );
            };

            uX.interpolate( initX, maxLevel, hyteg::All );
            uY.interpolate( initY, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            M.interpolate( initM, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_mass_u_p2_supg_blending_q5_ElementwiseOperator Op( storage, minLevel, maxLevel, uX, uY, uX, delta );
            P2ElementwiseBlendingMassOperator                 P2Mass( storage, minLevel, maxLevel );

            Op.apply( M, res, maxLevel, All, Replace );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }

         for ( auto s : Meshes3D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               SphericalCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t > T( "T", storage, minLevel, maxLevel );
            P2Function< real_t > M( "M", storage, minLevel, maxLevel );
            P2Function< real_t > O( "O", storage, minLevel, maxLevel );
            P2Function< real_t > uX( "uX", storage, minLevel, maxLevel );
            P2Function< real_t > uY( "uY", storage, minLevel, maxLevel );
            P2Function< real_t > uZ( "uZ", storage, minLevel, maxLevel );
            P2Function< real_t > res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return 1.0;
            };

            std::function< real_t( const hyteg::Point3D& ) > initM = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * x[0] * x[2] + 43 * x[0] + x[1] * x[2];
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * std::pow( x[2], 2 ) + 93;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] + 2 * x[1] - x[2];
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0] * x[1] + std::pow( x[2], 2 );
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -4 * x[2] * ( x[0] * x[1] - std::pow( x[2], 2 ) ) * ( 2 * x[0] * x[2] + 43 * x[0] + x[1] * x[2] );
            };

            uX.interpolate( initX, maxLevel, hyteg::All );
            uY.interpolate( initY, maxLevel, hyteg::All );
            uZ.interpolate( initZ, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            M.interpolate( initM, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_mass_u_p2_supg_blending_q5_ElementwiseOperator Op( storage, minLevel, maxLevel, uX, uY, uZ, delta );
            P2ElementwiseBlendingMassOperator                 P2Mass( storage, minLevel, maxLevel );

            Op.apply( M, res, maxLevel, All, Replace );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }
      }
   }

   {
      // Shear Heating SUPG Operator

      // Extract the required parameters
      uint_t minLevel  = 4;
      uint_t maxLevel  = 4;
      real_t tolerance = 1e-5;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 2;
         maxLevel  = 2;
         tolerance = 2e-3;
      }

      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/unitsquare_with_circular_hole.msh",
                                              "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el.msh", "../../../data/meshes/3D/pyramid_4el.msh" };

      for ( size_t i = 0; i < 2; i++ )
      {
         if ( i )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "--- With blending ---" );
         }

         for ( auto s : Meshes2D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               PolarCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t > T( "T", storage, minLevel, maxLevel );
            P2Function< real_t > O( "O", storage, minLevel, maxLevel );
            P2Function< real_t > uX( "uX", storage, minLevel, maxLevel );
            P2Function< real_t > uY( "uY", storage, minLevel, maxLevel );
            P2Function< real_t > res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta_Di_Ra = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return 1.0;
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 ) - 5;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[0], 2 ) + x[0] * x[1] + 2;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return std::pow( x[1], 2 );
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 0;
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * std::pow( x[1], 3 ) * std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) ) *
                      ( std::pow( x[0], 2 ) + std::pow( 2 * x[0] - x[1], 2 ) ) * ( std::pow( x[1], 2 ) - 5 );
            };

            std::function< real_t( const hyteg::Point3D&, real_t ) > f_fem = []( const hyteg::Point3D& x, real_t temp ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( temp );
               real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
               return radius * temp;
            };

            uX.interpolate( initX, maxLevel, hyteg::All );
            uY.interpolate( initY, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_shear_heat_T_p2_dep_eta_supg_blending_q5_ElementwiseOperator Op(
                storage, minLevel, maxLevel, uX, uY, uX, T, f_fem, delta_Di_Ra );
            P2ElementwiseBlendingMassOperator P2Mass( storage, minLevel, maxLevel );

            Op.apply( O, res, maxLevel, All, Replace );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }

         for ( auto s : Meshes3D )
         {
            // Init setup storage
            MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
            SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

            // Geometry map
            if ( i )
            {
               SphericalCoordsMap::setMap( setupStorage );
            }

            // Create storage
            std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

            // Create functions
            P2Function< real_t > T( "T", storage, minLevel, maxLevel );
            P2Function< real_t > O( "O", storage, minLevel, maxLevel );
            P2Function< real_t > uX( "uX", storage, minLevel, maxLevel );
            P2Function< real_t > uY( "uY", storage, minLevel, maxLevel );
            P2Function< real_t > uZ( "uZ", storage, minLevel, maxLevel );
            P2Function< real_t > res( "res", storage, minLevel, maxLevel );

            P2Function< real_t > Ctrl( "Ctrl", storage, minLevel, maxLevel );
            P2Function< real_t > CtrlRes( "CtrlRes", storage, minLevel, maxLevel );

            std::function< real_t( const hyteg::Point3D&, real_t ) > delta_Di_Ra = []( const hyteg::Point3D& x, real_t u_abs ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( u_abs );
               return 1.0;
            };

            std::function< real_t( const hyteg::Point3D& ) > initT = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return 2 * std::pow( x[2], 2 ) + 93;
            };

            std::function< real_t( const hyteg::Point3D& ) > initX = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] - 3 * x[0] + std::pow( x[2], 2 ) - 12;
            };
            std::function< real_t( const hyteg::Point3D& ) > initY = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return x[0] * x[1] + 2 * x[1] - x[2];
            };
            std::function< real_t( const hyteg::Point3D& ) > initZ = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -x[0] * x[1] + std::pow( x[2], 2 );
            };

            std::function< real_t( const hyteg::Point3D& ) > ctrl = []( const hyteg::Point3D& x ) {
               WALBERLA_UNUSED( x );
               return -4 * x[2] * ( x[0] * x[1] - std::pow( x[2], 2 ) ) *
                      std::sqrt( std::pow( x[0], 2 ) + std::pow( x[1], 2 ) + std::pow( x[2], 2 ) ) *
                      ( 9 * std::pow( x[0] + 1, 2 ) + 9 * std::pow( x[0] + x[1], 2 ) + 9 * std::pow( x[1] - 2 * x[2], 2 ) +
                        2 * std::pow( x[0] - 2 * x[1] + 2 * x[2] + 8, 2 ) + 2 * std::pow( x[0] + x[1] - 4 * x[2] - 1, 2 ) +
                        2 * std::pow( 2 * x[0] - x[1] - 2 * x[2] + 7, 2 ) ) /
                      ( 18 * std::pow( x[2], 2 ) + 837 );
            };

            std::function< real_t( const hyteg::Point3D&, real_t ) > f_fem = []( const hyteg::Point3D& x, real_t temp ) {
               WALBERLA_UNUSED( x );
               WALBERLA_UNUSED( temp );
               real_t radius = std::sqrt( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] );
               return radius / temp;
            };

            uX.interpolate( initX, maxLevel, hyteg::All );
            uY.interpolate( initY, maxLevel, hyteg::All );
            uZ.interpolate( initZ, maxLevel, hyteg::All );

            T.interpolate( initT, maxLevel, hyteg::All );
            O.interpolate( 1, maxLevel, hyteg::All );
            Ctrl.interpolate( ctrl, maxLevel, hyteg::All );

            // operators
            p2_shear_heat_T_p2_dep_eta_supg_blending_q5_ElementwiseOperator Op(
                storage, minLevel, maxLevel, uX, uY, uZ, T, f_fem, delta_Di_Ra );
            P2ElementwiseBlendingMassOperator P2Mass( storage, minLevel, maxLevel );

            Op.apply( O, res, maxLevel, All, Replace );
            P2Mass.apply( Ctrl, CtrlRes, maxLevel, All );

            real_t integVal0    = T.dotGlobal( res, maxLevel, All );
            real_t integValCtrl = O.dotGlobal( CtrlRes, maxLevel, All );

            WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
            WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
            WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
            WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
         }
      }
   }

   {
      // Shear Heating Operator

      // Extract the required parameters
      uint_t minLevel  = 5;
      uint_t maxLevel  = 5;
      real_t tolerance = 1e-4;

      WALBERLA_DEBUG_SECTION()
      {
         minLevel  = 3;
         maxLevel  = 3;
         tolerance = 3e-3;
      }

      // Meshes
      std::vector< std::string > Meshes2D = { "../../../data/meshes/2D/quad_4el.msh" };
      std::vector< std::string > Meshes3D = { "../../../data/meshes/3D/cube_6el.msh" };

      for ( auto s : Meshes2D )
      {
         // Init setup storage
         MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2Function< real_t >       T_( "T", storage, minLevel, maxLevel );
         P2Function< real_t >       O_( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D&, real_t ) > eta = []( const hyteg::Point3D& x, real_t temp ) {
            WALBERLA_UNUSED( x );
            WALBERLA_UNUSED( temp );
            return real_c( 5 ) * std::pow( x[0], 2 ) + real_c( 5 ) * std::pow( x[1], 2 );
         };

         std::function< real_t( const hyteg::Point3D& ) > rho = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return ( x[0] + x[1] + real_c( 5 ) );
         };

         std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return -std::pow( x[1], 2 ) - real_c( 5 ) * x[1];
         };

         std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::pow( x[0], 2 ) + real_c( 5 ) * x[0];
         };

         std::function< real_t( const hyteg::Point3D& ) > O = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] + std::pow( x[1], 2 ) - real_c( 5 );
         };

         u.interpolate( { uX, uY }, maxLevel, hyteg::All );

         O_.interpolate( O, maxLevel, hyteg::All );

         // operators
         auto ShearOp = std::make_shared< p2_shear_heat_T_p2_dep_eta_blending_q6_ElementwiseOperator >(
             storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 0 ), T_, eta, rho );

         ShearOp->apply( O_, res, maxLevel, All, Replace );

         real_t integVal0    = O_.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( 15265 ) / real_c( 63 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      for ( auto s : Meshes3D )
      {
         // Init setup storage
         MeshInfo              meshInfo = MeshInfo::fromGmshFile( s );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2Function< real_t >       T_( "T", storage, minLevel, maxLevel );
         P2Function< real_t >       O_( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D&, real_t ) > eta = []( const hyteg::Point3D& x, real_t temp ) {
            WALBERLA_UNUSED( x );
            WALBERLA_UNUSED( temp );
            return real_c( 5 ) * std::pow( x[0], 2 ) + real_c( 5 ) * std::pow( x[1], 2 ) + real_c( 5 ) * std::pow( x[2], 2 ) +
                   real_c( 1 );
         };

         std::function< real_t( const hyteg::Point3D& ) > rho = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return ( x[0] + x[1] + x[2] + real_c( 7 ) );
         };

         std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[2] + std::pow( x[2], 2 ) + real_c( 2 ) * x[2];
         };

         std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::pow( x[0], 2 ) + real_c( 2 ) * x[0] * x[2];
         };
         std::function< real_t( const hyteg::Point3D& ) > uZ = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return real_c( 2 ) * x[0] - std::pow( x[1], 2 ) - x[1];
         };
         std::function< real_t( const hyteg::Point3D& ) > O = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] - x[0] + std::pow( x[2], 2 ) - real_c( 12 );
         };

         u.interpolate( { uX, uY, uZ }, maxLevel, hyteg::All );

         O_.interpolate( O, maxLevel, hyteg::All );

         // operators
         auto ShearOp = std::make_shared< p2_shear_heat_T_p2_dep_eta_blending_q6_ElementwiseOperator >(
             storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 2 ), T_, eta, rho );

         ShearOp->apply( O_, res, maxLevel, All, Replace );

         real_t integVal0    = O_.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( 319923.293551587 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "Blending Test: " );
      {
         // Init setup storage
         real_t rCMB     = real_c( 1.0 );
         real_t rSurface = real_c( 2.0 );

         MeshInfo              meshInfo = MeshInfo::meshAnnulus( rCMB, rSurface, MeshInfo::meshFlavour::CRISS, 25, 4 );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         AnnulusMap::setMap( setupStorage );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2Function< real_t >       T_( "T", storage, minLevel, maxLevel );
         P2Function< real_t >       O_( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D, real_t ) > eta = []( const hyteg::Point3D& x, real_t temp ) {
            WALBERLA_UNUSED( x );
            WALBERLA_UNUSED( temp );
            return real_c( 5 ) * std::pow( x[0], 2 ) + real_c( 5 ) * std::pow( x[1], 2 );
         };

         std::function< real_t( const hyteg::Point3D& ) > rho = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return ( x[0] + x[1] + real_c( 5 ) );
         };

         std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return -std::pow( x[1], 2 ) - real_c( 5 ) * x[1];
         };

         std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::pow( x[0], 2 ) + real_c( 5 ) * x[0];
         };

         std::function< real_t( const hyteg::Point3D& ) > O = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] + std::pow( x[1], 2 ) - real_c( 5 );
         };

         u.interpolate( { uX, uY }, maxLevel, hyteg::All );

         O_.interpolate( O, maxLevel, hyteg::All );

         // operators
         auto ShearOp = std::make_shared< p2_shear_heat_T_p2_dep_eta_blending_q6_ElementwiseOperator >(
             storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 0 ), T_, eta, rho );

         ShearOp->apply( O_, res, maxLevel, All, Replace );

         real_t integVal0    = O_.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = -real_c( 12750 ) * real_c( M_SQRT2 ) *
                                   ( -real_c( 3 ) / real_c( 32 ) * real_c( M_SQRT2 ) +
                                     ( real_c( 1 ) / real_c( 8 ) ) * real_c( M_SQRT2 ) * real_c( M_PI ) ) -
                               real_c( 19125 ) / real_c( 8 ) + real_c( 31485 ) * real_c( M_PI );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }

      {
         // Init setup storage
         real_t rCMB     = real_c( 0.2 );
         real_t rSurface = real_c( 0.4 );

         // create the spherical shell mesh
         MeshInfo              meshInfo = MeshInfo::meshSphericalShell( 3, 2, rCMB, rSurface );
         SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

         IcosahedralShellMap::setMap( setupStorage );

         // Create storage
         std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

         // Create functions
         P2Function< real_t >       T_( "T", storage, minLevel, maxLevel );
         P2Function< real_t >       O_( "O", storage, minLevel, maxLevel );
         P2VectorFunction< real_t > u( "u", storage, minLevel, maxLevel );

         P2Function< real_t > res( "res", storage, minLevel, maxLevel );
         P2Function< real_t > rho_( "rho", storage, minLevel, maxLevel );
         P2Function< real_t > eta_( "eta", storage, minLevel, maxLevel );

         std::function< real_t( const hyteg::Point3D&, real_t ) > eta = []( const hyteg::Point3D& x, real_t temp ) {
            WALBERLA_UNUSED( x );
            return real_c( 5 ) * std::pow( x[0], 2 ) + real_c( 5 ) * std::pow( x[1], 2 ) + real_c( 5 ) * std::pow( x[2], 2 ) +
                   real_c( 1 );
         };

         std::function< real_t( const hyteg::Point3D& ) > rho = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return ( x[0] + x[1] + x[2] + real_c( 7 ) );
         };

         std::function< real_t( const hyteg::Point3D& ) > uX = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[1] * x[2] + std::pow( x[2], 2 ) + real_c( 2 ) * x[2];
         };

         std::function< real_t( const hyteg::Point3D& ) > uY = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return std::pow( x[0], 2 ) + real_c( 2 ) * x[0] * x[2];
         };
         std::function< real_t( const hyteg::Point3D& ) > uZ = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return real_c( 2 ) * x[0] - std::pow( x[1], 2 ) - x[1];
         };
         std::function< real_t( const hyteg::Point3D& ) > O = []( const hyteg::Point3D& x ) {
            WALBERLA_UNUSED( x );
            return x[0] * x[1] - x[0] + std::pow( x[2], 2 ) - real_c( 12 );
         };

         u.interpolate( { uX, uY, uZ }, maxLevel, hyteg::All );

         O_.interpolate( O, maxLevel, hyteg::All );

         // operators
         auto ShearOp = std::make_shared< p2_shear_heat_T_p2_dep_eta_blending_q6_ElementwiseOperator >(
             storage, minLevel, maxLevel, u.component( 0 ), u.component( 1 ), u.component( 2 ), T_, eta, rho );

         ShearOp->apply( O_, res, maxLevel, All, Replace );

         real_t integVal0    = O_.dotGlobal( res, maxLevel, All );
         real_t integValCtrl = real_c( 6504.77196038933 );

         WALBERLA_LOG_INFO_ON_ROOT( integVal0 );
         WALBERLA_LOG_INFO_ON_ROOT( integValCtrl );
         WALBERLA_LOG_INFO_ON_ROOT( abs( integVal0 - integValCtrl ) );
         WALBERLA_CHECK_LESS( abs( integVal0 - integValCtrl ), tolerance );
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Finish" );

   return EXIT_SUCCESS;
}