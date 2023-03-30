/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include <core/Environment.h>
#include <core/Format.hpp>
#include <core/config/Create.h>
#include <core/math/Constants.h>
#include <core/mpi/Broadcast.h>
#include <core/timing/Timer.h>

#include "hyteg/adaptiverefinement/mesh.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/petsc/PETScVector.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/controlflow/AgglomerationWrapper.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

#define R_min 1.0
#define R_max 2.0

using Mass    = P1ConstantMassOperator;
using Laplace = P1ConstantLaplaceOperator;
// using Laplace      = P1ElementwiseLaplaceOperator;
using DivkGradForm = forms::p1_div_k_grad_affine_q3;
using DivkGrad     = P1VariableOperator< DivkGradForm >;

// representation of a model problem
struct ModelProblem
{
   enum Type
   {
      // constant coefficient
      DIRAC,
      DIRAC_REGULARIZED,
      SLIT,
      REENTRANT_CORNER,
      WAVES,
      __VARIABLE_COEFFICIENT,
      // variable coefficient
      JUMP_COEFF,
      JUMP_COEFF_REGULARIZED,
      SPIKE,
      __NOT_AVAILABLE
   };

   ModelProblem( uint_t mp_type, uint_t spacial_dim )
   : type( Type( mp_type ) )
   , dim( spacial_dim )
   , sigma( 0.0 )
   , offset_err( 0.0 )
   , ignore_offset_for_refinement( false )
   , n_w( 1.0 )
   , p_w( 1.0 )
   , k_lo( 1.0 )
   , k_hi( 1.0 )
   , k_spike( 1.0 )
   {
      if ( mp_type >= Type::__NOT_AVAILABLE )
      {
         WALBERLA_ABORT( "Invalid argument for model problem: Must be less than " << Type::__NOT_AVAILABLE );
      }
      if ( mp_type == Type::__VARIABLE_COEFFICIENT )
      {
         WALBERLA_ABORT( "Invalid argument for model problem: " << Type::__VARIABLE_COEFFICIENT << " undefined!" );
      }
      if ( spacial_dim < 2 || 3 < spacial_dim )
      {
         WALBERLA_ABORT( "Invalid argument for spacial dimension: Must be either 2 or 3" );
      }
   }

   std::string name() const
   {
      std::stringstream ss;

      switch ( type )
      {
      case DIRAC:
         ss << "dirac"; // << "_" << offset_err;
         if ( ignore_offset_for_refinement )
            ss << "_fullyRefined";
         break;

      case DIRAC_REGULARIZED:
         ss << "dirac_regularized"; // << "_" << sigma;
         break;

      case SLIT:
         ss << "slit";
         break;

      case REENTRANT_CORNER:
         ss << "reentrant_corner";
         break;

      case WAVES:
         ss << "waves";
         break;

      case JUMP_COEFF:
         ss << "jump_coeff";
         break;

      case JUMP_COEFF_REGULARIZED:
         ss << "jump_coeff_regularized";
         break;

      case SPIKE:
         ss << "spike";
         break;

      default:
         break;
      }

      ss << "_" << dim << "D";

      return ss.str();
   }

   void set_params( real_t sig,
                    real_t offset,
                    bool   refine_all,
                    uint_t n_waves,
                    uint_t squeeze_waves,
                    real_t kLo,
                    real_t kHi,
                    real_t kSpike )
   {
      if ( type == DIRAC )
      {
         offset_err                   = offset;
         ignore_offset_for_refinement = refine_all;
      }
      if ( type == DIRAC_REGULARIZED )
      {
         sigma = sig;
      }
      if ( type == WAVES )
      {
         n_w = real_c( n_waves );
         p_w = real_c( squeeze_waves );
      }
      if ( type == JUMP_COEFF )
      {
         k_lo = kLo;
         k_hi = kHi;
      }
      if ( type == JUMP_COEFF_REGULARIZED )
      {
         sigma = sig;
      }
      if ( type == SPIKE )
      {
         k_spike = kSpike;
      }
      // todo: add parameters when adding new problems
   }

   bool constant_coefficient() const { return type < __VARIABLE_COEFFICIENT; }

   auto coeff() const
   {
      std::function< real_t( const Point3D& ) > k;

      if ( constant_coefficient() )
      {
         k = [=]( const Point3D& ) { return 1.0; };
      }
      if ( type == JUMP_COEFF )
      {
         k = [=]( const hyteg::Point3D& x ) {
            // jump at { (x,y) | x + (x_jump_0 - x_jump_1)y = x_jump_0 }
            real_t x_jump_y = x_jump_0 * ( real_c( 1.0 ) - x[1] ) + x_jump_1 * x[1];
            return ( x[0] < x_jump_y ) ? k_lo : k_hi;
         };
      }
      if ( type == JUMP_COEFF_REGULARIZED )
      {
         k = [=]( const hyteg::Point3D& x ) {
            return ( 2 * exp( sigma * ( x[0] - 0.5 ) ) + 1 ) / ( exp( sigma * ( x[0] - 0.5 ) ) + 1 );
         };
      }
      if ( type == SPIKE )
      {
         k = [=]( const hyteg::Point3D& x ) { return exp( k_spike * ( x[0] - x[0] * x[0] ) * ( x[1] - x[1] * x[1] ) ); };
      }
      // todo: add required coefficients

      return k;
   }

   void init( std::shared_ptr< PrimitiveStorage > storage,
              P1Function< real_t >&               u,
              P1Function< real_t >&               f,
              P1Function< real_t >&               u_h,
              P1Function< real_t >&               b,
              uint_t                              lvl,
              uint_t                              additional_lvl ) const
   {
      using walberla::math::one_div_root_two;
      using walberla::math::pi;

      // setup rhs function f and analytic solution u of pde
      std::function< real_t( const Point3D& ) > _f;
      std::function< real_t( const Point3D& ) > _u;

      if ( type == DIRAC_REGULARIZED )
      {
         const real_t sqrt_2pi_pow_3 = std::pow( 2.0 * pi, 3.0 / 2.0 ); //sqrt(2π)^3
         const real_t sig2           = sigma * sigma;
         const real_t sig3           = sig2 * sigma;

         if ( dim == 2 )
         {
            WALBERLA_ABORT( "regularized dirac not implemented for 2d" );
         }
         else
         {
            _f = [=]( const Point3D& x ) -> real_t {
               return std::exp( -0.5 * x.squaredNorm() / sig2 ) / ( sqrt_2pi_pow_3 * sig3 );
            };

            _u = [this, sqrt_2pi_pow_3]( const Point3D& x ) -> real_t {
               real_t r = x.norm();
               if ( r <= 0.0 )
                  return 1.0 / ( sqrt_2pi_pow_3 * sigma );
               else
                  return std::erf( one_div_root_two * r / sigma ) / ( 4.0 * pi * r );
            };
         }
      }

      if ( type == DIRAC )
      {
         _f = [=]( const Point3D& x ) -> real_t {
            return ( x.norm() < 1e-100 ) ? std::numeric_limits< real_t >::infinity() : 0.0;
         };

         if ( dim == 2 )
         {
            _u = [=]( const Point3D& x ) -> real_t { return -std::log( x.norm() ) / ( 2.0 * pi ); };
         }
         else
         {
            _u = [=]( const Point3D& x ) -> real_t { return 1.0 / ( 4.0 * pi * x.norm() ); };
         }
      }

      if ( type == SLIT || type == REENTRANT_CORNER )
      {
         const real_t alpha = ( type == SLIT ) ? 1.0 / 2.0 : 2.0 / 3.0;

         if ( type == REENTRANT_CORNER && dim == 3 )
         {
            WALBERLA_ABORT( "" << name() << " not implemented for 3d" );
         }

         _f = []( const Point3D& ) -> real_t { return 0.0; };
         _u = [=]( const Point3D& x ) -> real_t {
            auto r   = x.norm();
            auto phi = std::atan2( x[1], x[0] ); // only valid in 2D!
            if ( phi < 0 )
               phi += 2 * pi;

            return std::pow( r, alpha ) * std::sin( alpha * phi );
         };
      }

      if ( type == WAVES )
      {
         if ( dim == 3 )
         {
            WALBERLA_ABORT( "" << name() << " not implemented for 3d" );
         }

         const real_t rot1 = real_c( 0.1 );
         const real_t rot2 = real_c( 0.7 );

         _u = [=]( const hyteg::Point3D& x ) {
            return ( exp( p_w * ( -pow( x[0], 2 ) + x[0] ) * ( -pow( x[1], 2 ) + x[1] ) ) - 1 ) *
                   sin( 2 * n_w * pi * ( x[0] + x[1] ) );
         };
         _f = [=]( const hyteg::Point3D& x ) {
            auto x0  = n_w * pi;
            auto x1  = 2 * x0 * ( x[0] + x[1] );
            auto x2  = sin( x1 );
            auto x3  = -pow( x[1], 2 ) + x[1];
            auto x4  = -pow( x[0], 2 ) + x[0];
            auto x5  = p_w * x4;
            auto x6  = exp( x3 * x5 );
            auto x7  = x2 * x6;
            auto x8  = 2 * x7;
            auto x9  = p_w * x3;
            auto x10 = 1 - 2 * x[0];
            auto x11 = 4 * x0 * x6 * cos( x1 );
            auto x12 = 1 - 2 * x[1];
            auto x13 = pow( p_w, 2 ) * x7;
            return 8 * pow( n_w, 2 ) * pow( pi, 2 ) * x2 * ( x6 - 1 ) - pow( x10, 2 ) * x13 * pow( x3, 2 ) - x10 * x11 * x9 -
                   x11 * x12 * x5 - pow( x12, 2 ) * x13 * pow( x4, 2 ) + x5 * x8 + x8 * x9;
         };
      }

      if ( type == JUMP_COEFF )
      {
         if ( dim == 3 )
         {
            WALBERLA_ABORT( "" << name() << " not implemented for 3d" );
         }

         real_t du_dx_min = real_c( 1.0 ) / k_lo;
         real_t du_dx_max = real_c( 1.0 ) / k_hi;
         real_t du_dy_min = ( x_jump_0 - x_jump_1 ) / k_lo;
         real_t du_dy_max = ( x_jump_0 - x_jump_1 ) / k_hi;
         real_t u_jump    = du_dx_min * x_jump_0;

         _u = [=]( const hyteg::Point3D& x ) {
            // kink at { (x,y) | x + (x_jump_0 - x_jump_1)y = x_jump_0 }
            real_t x_jump_y = x_jump_0 * ( real_c( 1.0 ) - x[1] ) + x_jump_1 * x[1];
            if ( x[0] < x_jump_y )
            {
               return du_dx_min * x[0] + du_dy_min * x[1];
            }
            else
            {
               return u_jump + du_dx_max * ( x[0] - x_jump_0 ) + du_dy_max * x[1];
            }
         };

         _f = []( const Point3D& ) -> real_t { return 0.0; };
      }
      if ( type == JUMP_COEFF_REGULARIZED )
      {
         if ( dim == 3 )
         {
            WALBERLA_ABORT( "" << name() << " not implemented for 3d" );
         }
         const auto scaling_factor_u = real_c( 1.0 );
         // real_t du_dx_min = real_c( 1.0 ) / k_lo;
         // real_t du_dx_max = real_c( 1.0 ) / k_hi;
         // real_t du_dy_min = ( x_jump_0 - x_jump_1 ) / k_lo;
         // real_t du_dy_max = ( x_jump_0 - x_jump_1 ) / k_hi;
         // real_t u_jump    = du_dx_min * x_jump_0;

         _u = [=]( const hyteg::Point3D& x ) {
            return scaling_factor_u * ( -2 * sigma * x[0] + log( 2 * exp( sigma * ( x[0] - 0.5 ) ) + 1 ) ) / sigma;
         };

         _f = []( const Point3D& ) -> real_t { return 0.0; };
      }
      if ( type == SPIKE )
      {
         if ( dim == 3 )
         {
            WALBERLA_ABORT( "" << name() << " not implemented for 3d" );
         }

         _u = [=]( const hyteg::Point3D& x ) {
            return real_c( 1.0 ) - exp( -k_spike * ( x[0] - x[0] * x[0] ) * ( x[1] - x[1] * x[1] ) );
         };

         _f = [=]( const hyteg::Point3D& x ) {
            return real_c( 2.0 ) * k_spike * ( ( x[0] - x[0] * x[0] ) + ( x[1] - x[1] * x[1] ) );
         };
      }

      // interpolate rhs of pde
      f.interpolate( _f, lvl );
      // construct rhs of linear system s.th. b_i = ∫fφ_i
      if ( type == DIRAC )
      {
         // ∫δ_0 φ_i = φ_i(0)
         auto _b = [=]( const Point3D& x ) -> real_t { return ( _f( x ) > 0.0 ) ? 1.0 : 0.0; };
         b.interpolate( _b, lvl );
      }
      else
      {
         // ∫f φ_i = [Mf]_i
         Mass M( storage, lvl, lvl );
         M.apply( f, b, lvl, All );
      }

      // interpolate analytic solution
      u.interpolate( _u, lvl );
      if ( additional_lvl )
      {
         u.interpolate( _u, lvl + additional_lvl );
      }

      // initialize u_h
      u_h.setToZero( lvl );
      u_h.interpolate( _u, lvl, DirichletBoundary );
      if ( additional_lvl )
      {
         u_h.setToZero( lvl + additional_lvl );
         u_h.interpolate( _u, lvl + additional_lvl, DirichletBoundary );
      }
   }

   void apply_error_filter( const P1Function< real_t >& e, P1Function< real_t >& e_f, uint_t lvl ) const
   {
      e.communicate< Vertex, Edge >( lvl );
      e.communicate< Edge, Face >( lvl );
      e.communicate< Face, Cell >( lvl );

      auto filter = [this]( const Point3D& x, const std::vector< real_t >& val ) -> real_t {
         return ( x.norm() >= offset_err ) ? val[0] : 0.0;
      };
      e_f.interpolate( filter, { e }, lvl );

      e_f.communicate< Vertex, Edge >( lvl );
      e_f.communicate< Edge, Face >( lvl );
      e_f.communicate< Face, Cell >( lvl );
   }

   // general setting
   const Type   type;
   const uint_t dim;
   // parameters
   real_t sigma;      // standard deviation
   real_t offset_err; // offset from singularity for error computation
   bool   ignore_offset_for_refinement;
   real_t n_w;            // number of waves
   real_t p_w;            // 'squeeze' parameter of waves
   real_t k_lo, k_hi;     // values of jump coefficient
   real_t x_jump_0 = 0.4; // position of jump
   real_t x_jump_1 = 0.6;
   real_t k_spike; // severity of spike
};

SetupPrimitiveStorage domain( const ModelProblem& problem, uint_t N, const std::string& inputMsh )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( inputMsh != "" )
   {
      meshInfo = MeshInfo::fromGmshFile( inputMsh );
   }
   else if ( problem.dim == 2 )
   {
      if ( problem.type == ModelProblem::REENTRANT_CORNER )
      {
         if ( N % 2 != 0 || N < 2 || N > 10 )
         {
            WALBERLA_ABORT( "Initial resolution for reentrant corner must be an even number with 2 <= N <= 10!" );
         }
         uint_t n_el     = N * N * 3 / 2;
         auto   filename = "data/meshes/LShape_" + std::to_string( n_el ) + "el.msh";
         meshInfo        = MeshInfo::fromGmshFile( filename );
      }
      else if ( problem.type == ModelProblem::WAVES || problem.type == ModelProblem::JUMP_COEFF ||
                problem.type == ModelProblem::SPIKE )
      {
         // (0,1)^2
         Point2D n( { 1, 1 } );
         meshInfo = MeshInfo::meshRectangle( { 0, 0 }, n, MeshInfo::CRISS, N, N );
      }
      else
      {
         // (-1,1)^2
         Point2D n( { 1, 1 } );
         meshInfo = MeshInfo::meshRectangle( -n, n, MeshInfo::CRISS, N, N );
      }
   }
   else
   {
      if ( problem.type == ModelProblem::SLIT )
      {
         if ( N % 2 != 0 )
         {
            WALBERLA_ABORT( "Initial resolution for slit domain must be an even number!" );
         }
         Point3D n( { 1, 1, 0.1 } );
         meshInfo = MeshInfo::meshCuboid( -n, n, N, N, 1 );
      }
      else
      {
         Point3D n( { 1, 1, 1 } );
         meshInfo = MeshInfo::meshCuboid( -n, n, N, N, N );
      }
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // setup boundary conditions

   // dirichlet boundary conditions on entire boundary
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   if ( problem.type == ModelProblem::SLIT )
   {
      if ( problem.dim == 3 )
      {
         // neumann boundary on top and bottom of plate
         setupStorage.setMeshBoundaryFlagsOnBoundary( 2, 0, true );
         auto onBoundary = []( const Point3D& x ) -> bool {
            return std::abs( x[0] ) >= 1 - 1e-15 || //
                   std::abs( x[1] ) >= 1 - 1e-15;
         };
         setupStorage.setMeshBoundaryFlagsByCentroidLocation( 1, onBoundary, true );
      }
      // dirichlet boundary on slit
      auto onBoundary = []( const Point3D& x ) -> bool { return std::abs( x[1] ) < 1e-50 && x[0] >= 0.0; };
      setupStorage.setMeshBoundaryFlagsByVertexLocation( 1, onBoundary, true );
   }

   return setupStorage;
}

// solve problem with current refinement and return list of elementwise squared errors of local elements
template < class A_t >
adaptiveRefinement::ErrorVector solve( adaptiveRefinement::Mesh&                mesh,
                                       const ModelProblem&                      problem,
                                       uint_t                                   u0,
                                       std::shared_ptr< P1Function< real_t > >& u_old,
                                       uint_t                                   l_interpolate,
                                       uint_t                                   l_min,
                                       uint_t                                   l_max,
                                       uint_t                                   max_iter,
                                       real_t                                   tol,
                                       real_t                                   cg_tol,
                                       std::string                              vtkname,
                                       bool                                     writePartitioning,
                                       bool                                     writeMeshfile,
                                       uint_t                                   refinement_step,
                                       uint_t                                   accurate_error,
                                       bool                                     l2_error_each_iteration = true )
{
   auto l_err = l_max + accurate_error;

   double t0, t1;
   double t_loadbalancing = 0;
   if ( u0 == 0 || u_old == nullptr )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "* apply Loadbalancing" );
      // if u0 is initialized with zero, we apply loadbalancing before creating the storage.
      // else, we first interpolate u before applying loadbalancing
      t0 = walberla::timing::getWcTime();
      mesh.loadbalancing( adaptiveRefinement::Loadbalancing::ROUND_ROBIN );
      t1              = walberla::timing::getWcTime();
      t_loadbalancing = t1 - t0;
   }
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "* create PrimitiveStorage ..." );
   t0           = walberla::timing::getWcTime();
   auto storage = mesh.make_storage();
   t1           = walberla::timing::getWcTime();
   printCurrentMemoryUsage();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> Time spent to create storage: %12.3e", t1 - t0 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( "* setup system ..." );

   // operators
   t0 = walberla::timing::getWcTime();
   std::shared_ptr< A_t > A;
   auto                   M = std::make_shared< Mass >( storage, l_min, l_err );
   auto                   P = std::make_shared< P1toP1LinearProlongation<> >();
   auto                   R = std::make_shared< P1toP1LinearRestriction<> >();

   if constexpr ( std::is_same_v< A_t, DivkGrad > )
   {
      auto k = problem.coeff();

      A = std::make_shared< A_t >( storage, l_min, l_max + 1, DivkGradForm( k, k ) );
   }
   else
   {
      A = std::make_shared< A_t >( storage, l_min, l_err );
   }

   // FE functions
   auto                 u = std::make_shared< P1Function< real_t > >( "u_h", storage, l_min, l_err );
   P1Function< real_t > f( "f", storage, l_max, l_max );
   P1Function< real_t > b( "b(f)", storage, l_min, l_max );
   P1Function< real_t > r( "r", storage, l_max, l_max );
   P1Function< real_t > u_anal( "u_anal", storage, l_max, l_err );
   P1Function< real_t > err( "e_h", storage, l_max, l_err );
   P1Function< real_t > tmp( "tmp", storage, l_min, l_err );
   P1Function< real_t > err_f( "e_h_filtered", storage, l_max, l_err );

   // global DoF
   tmp.interpolate( []( const Point3D& ) { return 1.0; }, l_max, Inner | NeumannBoundary );
   auto n_dof = uint_t( tmp.dotGlobal( tmp, l_max ) );
   // cg-dof
   tmp.interpolate( []( const Point3D& ) { return 1.0; }, l_min, Inner | NeumannBoundary );
   auto n_dof_coarse = uint_t( tmp.dotGlobal( tmp, l_min ) );

   t0          = walberla::timing::getWcTime();
   auto t_init = t1 - t0;
   WALBERLA_LOG_INFO_ON_ROOT( " -> number of global DoF: " << n_dof );

   // computation of residual and L2 error
   t0 = walberla::timing::getWcTime();
   err.setToZero( l_max );
   err.setToZero( l_err );
   tmp.setToZero( l_max );
   tmp.setToZero( l_err );
   double t_residual = 0, t_error = 0;
   auto   compute_residual = [&]() -> real_t {
      auto my_t0 = walberla::timing::getWcTime();
      A->apply( *u, tmp, l_max, Inner | NeumannBoundary, Replace );
      r.assign( { 1.0, -1.0 }, { b, tmp }, l_max, Inner | NeumannBoundary );
      M->apply( r, tmp, l_max, Inner | NeumannBoundary, Replace );
      auto norm_r = std::sqrt( r.dotGlobal( tmp, l_max, Inner | NeumannBoundary ) );
      auto my_t1  = walberla::timing::getWcTime();
      t_residual += my_t1 - my_t0;
      return norm_r;
   };
   auto compute_L2error = [&]() -> real_t {
      auto my_t0 = walberla::timing::getWcTime();
      for ( uint_t lvl = l_max; lvl < l_err; ++lvl )
      {
         P->prolongate( *u, lvl, Inner | NeumannBoundary );
      }
      err.assign( { 1.0, -1.0 }, { *u, u_anal }, l_err, Inner | NeumannBoundary );
      problem.apply_error_filter( err, err_f, l_err );
      M->apply( err_f, tmp, l_err, Inner | NeumannBoundary, Replace );
      auto norm_e = std::sqrt( err_f.dotGlobal( tmp, l_err, Inner | NeumannBoundary ) );
      auto my_t1  = walberla::timing::getWcTime();
      t_error += my_t1 - my_t0;
      return norm_e;
   };

   // initialize analytic solution, rhs and boundary values
   problem.init( storage, u_anal, f, *u, b, l_max, accurate_error );

   t1 = walberla::timing::getWcTime();
   t_init += t1 - t0;

   double t_interpolate = 0;

   // initialize u_h
   if ( u0 == 1 && u_old != nullptr )
   {
      // interpolate old solution to new mesh
      t0 = walberla::timing::getWcTime();
      WALBERLA_LOG_INFO_ON_ROOT( " -> initialize u=u_old" );
      u->interpolate( *u_old, l_max, l_interpolate );
      t1            = walberla::timing::getWcTime();
      t_interpolate = t1 - t0;
      // apply loadbalancing
      WALBERLA_LOG_INFO_ON_ROOT( " -> apply Loadbalancing" );
      t0                 = walberla::timing::getWcTime();
      auto migrationInfo = mesh.loadbalancing( adaptiveRefinement::Loadbalancing::ROUND_ROBIN );
      storage->migratePrimitives( migrationInfo );
      t1              = walberla::timing::getWcTime();
      t_loadbalancing = t1 - t0;
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( " -> initialize u=0" );
   }

   // solver
   t0 = walberla::timing::getWcTime();
   // smoother
   // auto smoother         = std::make_shared< WeightedJacobiSmoother< A_t > >(storage, l_min, l_max, 0.66);
   auto smoother = std::make_shared< GaussSeidelSmoother< A_t > >();
   // coarse grid solver
   auto cgIter = std::max( uint_t( 5 ), n_dof_coarse );
#ifdef HYTEG_BUILD_WITH_PETSC
   auto cgs = std::make_shared< PETScCGSolver< A_t > >( storage, l_min, 1e-30, cg_tol, cgIter );
#else
   auto cgs = std::make_shared< CGSolver< A_t > >( storage, l_min, l_min, cgIter, cg_tol );
#endif
   // multigrid
   GeometricMultigridSolver< A_t > gmg( storage, smoother, cgs, R, P, l_min, l_max, 3, 3 );
   t1 = walberla::timing::getWcTime();
   t_init += t1 - t0;

   WALBERLA_LOG_INFO_ON_ROOT( " -> Time spent to ...  " );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> %20s: %12.3e", "initialize problem", t_init ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> %20s: %12.3e", "interpolate u0", t_interpolate ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> %20s: %12.3e", "loadbalancing", t_loadbalancing ) );

   // solve
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "* solve system ..." );
   WALBERLA_LOG_INFO_ON_ROOT( " -> run multigrid solver" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6s |%18s |%13s", "k", "||r_k||_L2", "||e_k||_L2" ) );

   // initial residual
   real_t norm_r = compute_residual();

   uint_t iter    = 0;
   double t_solve = 0;
   while ( norm_r > tol && iter < max_iter )
   {
      if ( l2_error_each_iteration )
      {
         auto eL2 = compute_L2error();
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6d |%18.3e |%13.3e", iter, norm_r, eL2 ) );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6d |%18.3e |", iter, norm_r ) );
      }

      ++iter;
      t0 = walberla::timing::getWcTime();
      gmg.solve( *A, *u, b, l_max );
      t1 = walberla::timing::getWcTime();
      t_solve += t1 - t0;
      norm_r = compute_residual();
   }

   auto eL2 = compute_L2error();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6d |%18.3e |%13.3e", iter, norm_r, eL2 ) );

   // compute elementwise error
   t0 = walberla::timing::getWcTime();
   adaptiveRefinement::ErrorVector err_2_elwise_loc;

   P1Function< real_t >& err_for_refinement = ( problem.ignore_offset_for_refinement ) ? err : err_f;

   if ( problem.dim == 3 )
   {
      for ( auto& [id, cell] : storage->getCells() )
      {
         real_t err_2_cell = vertexdof::macrocell::dot< real_t >(
             l_err, *cell, err_for_refinement.getCellDataID(), err_for_refinement.getCellDataID(), 0 );

         // scale squared error by cell-volume
         std::array< Point3D, 3 + 1 > vertices;
         uint_t                       i = 0;
         for ( auto& vid : cell->neighborVertices() )
         {
            vertices[i] = storage->getVertex( vid )->getCoordinates();
            ++i;
         }
         err_2_cell *= adaptiveRefinement::Simplex3::volume( vertices );

         err_2_elwise_loc.push_back( { err_2_cell, id } );
      }
   }
   else // dim == 2
   {
      for ( auto& [id, face] : storage->getFaces() )
      {
         real_t err_2_face = vertexdof::macroface::dot< real_t >(
             l_err, *face, err_for_refinement.getFaceDataID(), err_for_refinement.getFaceDataID(), 0 );

         // scale squared error by face-volume
         std::array< Point3D, 2 + 1 > vertices;
         uint_t                       i = 0;
         for ( auto& vid : face->neighborVertices() )
         {
            vertices[i] = storage->getVertex( vid )->getCoordinates();
            ++i;
         }
         err_2_face *= adaptiveRefinement::Simplex2::volume( vertices );

         err_2_elwise_loc.push_back( { err_2_face, id } );
      }
   }

   t1                  = walberla::timing::getWcTime();
   auto t_elwise_error = t1 - t0;

   WALBERLA_LOG_INFO_ON_ROOT( " -> Time spent to ...  " );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> %20s: %12.3e", "system solve", t_solve ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> %20s: %12.3e", "compute residual", t_residual ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> %20s: %12.3e", "compute error", t_error ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> %20s: %12.3e", "compute elwise err", t_elwise_error ) );

   // export to vtk
   if ( vtkname != "" )
   {
      // coefficient
      P1Function< real_t > k( "k", storage, l_max, l_max );
      k.interpolate( problem.coeff(), l_max );

      // error
      for ( uint_t lvl = l_err; lvl > l_max; --lvl )
      {
         R->restrict( err, lvl, Inner | NeumannBoundary );
      }

      // squared error
      P1Function< real_t > err_2( "err^2", storage, l_min, l_max );
      err_2.multElementwise( { err, err }, l_max, All );

      // write vtkfile
      VTKOutput vtkOutput( "output", vtkname, storage );
      vtkOutput.setVTKDataFormat( vtk::DataFormat::BINARY );
      vtkOutput.add( u_anal );
      vtkOutput.add( k );
      vtkOutput.add( f );
      vtkOutput.add( *u );
      vtkOutput.add( b );
      vtkOutput.add( err );
      vtkOutput.add( err_2 );
      vtkOutput.write( l_max, refinement_step );
   }

   if ( writePartitioning )
   {
      std::map< std::string, std::map< PrimitiveID, real_t > > realData;

      writeDomainPartitioningVTK(
          *storage, "output", vtkname + "_partitioning_vertices_ts" + std::to_string( refinement_step ), VTK_VERTEX, realData );
      writeDomainPartitioningVTK(
          *storage, "output", vtkname + "_partitioning_edges_ts" + std::to_string( refinement_step ), VTK_LINE, realData );
      writeDomainPartitioningVTK(
          *storage, "output", vtkname + "_partitioning_faces_ts" + std::to_string( refinement_step ), VTK_TRIANGLE, realData );
      if ( storage->hasGlobalCells() )
      {
         writeDomainPartitioningVTK(
             *storage, "output", vtkname + "_partitioning_cells_ts" + std::to_string( refinement_step ), VTK_TETRA, realData );
      }
   }

   if ( writeMeshfile )
   {
      mesh.exportMesh( "output/" + vtkname + "_mesh_" + std::to_string( refinement_step ) + ".msh" );
   }

   // update u_old
   u_old.swap( u );

   return err_2_elwise_loc;
}

void solve_for_each_refinement( const SetupPrimitiveStorage& setupStorage,
                                const ModelProblem&          problem,
                                uint_t                       n_ref,
                                uint_t                       n_el_max,
                                real_t                       p_ref,
                                uint_t                       l_min,
                                uint_t                       l_max,
                                uint_t                       l_final,
                                uint_t                       u0,
                                uint_t                       max_iter,
                                real_t                       tol,
                                uint_t                       max_iter_final,
                                real_t                       tol_final,
                                real_t                       cg_tol,
                                uint_t                       accurate_error,
                                std::string                  vtkname,
                                bool                         writePartitioning,
                                bool                         writeMeshfile )
{
   // construct adaptive mesh
   adaptiveRefinement::Mesh mesh( setupStorage );
   // store solution from previous refinement as initial guess
   std::shared_ptr< P1Function< real_t > > u_old = nullptr;

   uint_t refinement = 0;
   while ( 1 )
   {
      adaptiveRefinement::ErrorVector local_errors;
      if ( problem.constant_coefficient() )
         local_errors = solve< Laplace >( mesh,
                                          problem,
                                          u0,
                                          u_old,
                                          l_max,
                                          l_min,
                                          l_max,
                                          max_iter,
                                          tol,
                                          cg_tol,
                                          vtkname,
                                          writePartitioning,
                                          writeMeshfile,
                                          refinement,
                                          accurate_error );
      else
         local_errors = solve< DivkGrad >( mesh,
                                           problem,
                                           u0,
                                           u_old,
                                           l_max,
                                           l_min,
                                           l_max,
                                           max_iter,
                                           tol,
                                           cg_tol,
                                           vtkname,
                                           writePartitioning,
                                           writeMeshfile,
                                           refinement,
                                           accurate_error );

      if ( refinement >= n_ref )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "" );
         WALBERLA_LOG_INFO_ON_ROOT( "* maximum number of refinements!" );
         break;
      }

      auto n_el_old = mesh.n_elements();
      ++refinement;
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "* apply refinement " << refinement );
      WALBERLA_LOG_INFO_ON_ROOT( " -> n_el_old = " << n_el_old );

      // refinement strategy
      const auto p_idx     = std::min( uint_t( std::ceil( real_t( n_el_old ) * p_ref ) ), n_el_old - 1 );
      auto       criterion = [&]( const adaptiveRefinement::ErrorVector& err_global, uint_t i ) -> bool {
         if ( i == 0 )
         {
            WALBERLA_LOG_INFO_ON_ROOT( " -> min_i err_i = " << err_global.back().first );
            WALBERLA_LOG_INFO_ON_ROOT( " -> max_i err_i = " << err_global.front().first );
            WALBERLA_LOG_INFO_ON_ROOT( " -> refining all elements i where err_i >= " << 0.9 * err_global[p_idx].first );
         }
         return err_global[i].first >= 0.9 * err_global[p_idx].first;
      };

      // apply refinement
      auto t0    = walberla::timing::getWcTime();
      auto ratio = mesh.refineRG( local_errors, criterion, n_el_max );
      auto t1    = walberla::timing::getWcTime();
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> Time spent for refinement: %12.3e", t1 - t0 ) );
      WALBERLA_LOG_INFO_ON_ROOT( " -> n_el_new = " << mesh.n_elements() );
      // compute mesh quality
      auto v_mean       = mesh.volume() / real_t( mesh.n_elements() );
      auto v_minmax     = mesh.min_max_volume();
      auto a_minmax     = mesh.min_max_angle();
      auto a_meanminmax = mesh.mean_min_max_angle();
      WALBERLA_LOG_INFO_ON_ROOT( " -> el volume (min, max, mean): " << v_minmax.first << ", " << v_minmax.second << ", "
                                                                    << v_mean );
      WALBERLA_LOG_INFO_ON_ROOT( " -> angle (min, max) over all elements: " << a_minmax.first << ", " << a_minmax.second );
      WALBERLA_LOG_INFO_ON_ROOT( " -> angle (min, max) mean over all elements: " << a_meanminmax.first << ", "
                                                                                 << a_meanminmax.second );

      if ( ratio < 0.95 )
      {
         WALBERLA_LOG_INFO_ON_ROOT(
             "* refinement could not be applied to all required elements\n  without exceeding the max number of elements!" );
         break;
      }
   }

   if ( l_final > l_max || max_iter_final > max_iter || tol_final < tol )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "* final solve with " );
      if ( l_final > l_max )
         WALBERLA_LOG_INFO_ON_ROOT( "    higher micro-level" );
      if ( max_iter_final > max_iter )
         WALBERLA_LOG_INFO_ON_ROOT( "    greater number of iterations" );
      if ( tol_final < tol )
         WALBERLA_LOG_INFO_ON_ROOT( "    smaller tolerance for residual" );

      if ( problem.constant_coefficient() )
      {
         solve< Laplace >( mesh,
                           problem,
                           u0,
                           u_old,
                           l_max,
                           l_min,
                           l_final,
                           max_iter_final,
                           tol_final,
                           cg_tol,
                           vtkname,
                           writePartitioning,
                           writeMeshfile,
                           refinement,
                           accurate_error );
      }
      else
      {
         solve< DivkGrad >( mesh,
                            problem,
                            u0,
                            u_old,
                            l_max,
                            l_min,
                            l_final,
                            max_iter_final,
                            tol_final,
                            cg_tol,
                            vtkname,
                            writePartitioning,
                            writeMeshfile,
                            refinement,
                            accurate_error );
      }

      if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
      {
         std::cout << "\nTimingtree final iteration:\n" << *( u_old->getStorage()->getTimingTree() ) << "\n";
      }
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg;

   if ( argc == 1 )
   {
      walberla::shared_ptr< walberla::config::Config > cfg_( new walberla::config::Config );
      cfg_->readParameterFile( "../../hyteg/data/param/adaptiveRefinement.prm" );
      cfg = cfg_;
   }
   else
   {
      cfg = walberla::config::create( argc, argv );
   }

   // read parameters
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t mp            = parameters.getParameter< uint_t >( "modelProblem" );
   const uint_t dim           = parameters.getParameter< uint_t >( "dim" );
   const real_t sigma         = parameters.getParameter< real_t >( "sigma", 0.0 );
   const real_t offset        = parameters.getParameter< real_t >( "offset", 0.0 );
   const bool   refine_all    = parameters.getParameter< bool >( "refine_all", false );
   const uint_t n_waves       = parameters.getParameter< uint_t >( "n_waves", 1u );
   const uint_t squeeze_waves = parameters.getParameter< uint_t >( "squeeze_waves", 1u );
   const real_t kLo           = parameters.getParameter< real_t >( "k_low", 1.0 );
   const real_t kHi           = parameters.getParameter< real_t >( "k_high", 1.0 );
   const real_t kSpike        = parameters.getParameter< real_t >( "k_spike", 1.0 );

   // todo: add parameters when adding new problems

   const uint_t N_default =
       ( ModelProblem::Type( mp ) == ModelProblem::SLIT || ModelProblem::Type( mp ) == ModelProblem::REENTRANT_CORNER ) ? 2 : 1;
   const uint_t N             = parameters.getParameter< uint_t >( "initial_resolution", N_default );
   const uint_t n_refinements = parameters.getParameter< uint_t >( "n_refinements" );
   const uint_t n_el_max      = parameters.getParameter< uint_t >( "n_el_max", std::numeric_limits< uint_t >::max() );
   const real_t p_refinement  = parameters.getParameter< real_t >( "percentile" );

   const uint_t l_min   = parameters.getParameter< uint_t >( "cg_level", 0 );
   const uint_t l_max   = parameters.getParameter< uint_t >( "microlevel" );
   const uint_t l_final = parameters.getParameter< uint_t >( "microlevel_final", l_max );

   const uint_t u0             = parameters.getParameter< uint_t >( "initial_guess", 0 );
   const uint_t max_iter       = parameters.getParameter< uint_t >( "n_iterations" );
   const uint_t max_iter_final = parameters.getParameter< uint_t >( "n_iterations_final", max_iter );
   const real_t tol            = parameters.getParameter< real_t >( "tolerance" );
   const real_t tol_final      = parameters.getParameter< real_t >( "tolerance_final", tol );
   const real_t cg_tol         = parameters.getParameter< real_t >( "cg_tolerance", tol );

   std::string vtkname           = parameters.getParameter< std::string >( "vtkName", "" );
   std::string inputmesh         = parameters.getParameter< std::string >( "initialMesh", "" );
   const bool  writePartitioning = parameters.getParameter< bool >( "writeDomainPartitioning", false );
   const bool  writeMeshfile     = parameters.getParameter< bool >( "writeMeshfile", false );

   const uint_t accurate_error = parameters.getParameter< uint_t >( "compute_error_on_higher_lvl", 0 );

#ifdef HYTEG_BUILD_WITH_PETSC
   PETScManager petscManager( &argc, &argv );
#endif

   // setup model problem
   ModelProblem problem( mp, dim );
   if ( vtkname == "auto" )
   {
      vtkname = problem.name();
   }
   problem.set_params( sigma, offset, refine_all, n_waves, squeeze_waves, kLo, kHi, kSpike );

   // setup domain
   auto setupStorage = domain( problem, N, inputmesh );

   // print parameters
   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d (%s)", "model problem", mp, problem.name().c_str() ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "spacial dimensions", dim ) );
   if ( ModelProblem::Type( mp ) == ModelProblem::DIRAC_REGULARIZED )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e", "sigma", sigma ) );
      if ( sigma <= real_t( 0.0 ) )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Sigma should be greater than zero to prevent singularity in the solution at (0,0,0)!" )
      }
   }
   if ( ModelProblem::Type( mp ) == ModelProblem::DIRAC )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e", "offset", offset ) );
      if ( refine_all )
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %s", "", "(but refine elements within offset)" ) );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %s", "", "(also no refinement within offset)" ) );
      }
      if ( offset <= real_t( 0.0 ) )
      {
         WALBERLA_LOG_WARNING_ON_ROOT(
             "Choosing an offset of zero will produce NaNs in the L2-error due to the singularity at (0,0,0)!" )
      }
   }
   if ( ModelProblem::Type( mp ) == ModelProblem::WAVES )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2d", "n_waves", n_waves ) );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2d", "squeeze", squeeze_waves ) );
   }
   if ( ModelProblem::Type( mp ) == ModelProblem::JUMP_COEFF )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e", "k_low", kLo ) );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e", "k_high", kHi ) );
   }
   if ( ModelProblem::Type( mp ) == ModelProblem::SPIKE )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e", "k_spike", kSpike ) );
   }
   if ( ModelProblem::Type( mp ) == ModelProblem::JUMP_COEFF_REGULARIZED )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e", "sigma", sigma ) );
   }
   // todo: add output for new parameters when adding problems
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "initial resolution", N ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "number of refinements", n_refinements ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "max. number of coarse elements", n_el_max ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e", "proportion to refine per step", p_refinement ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %s", "initial guess", ( u0 ) ? "interpolated" : "zero" ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d / %d", "level (min/max)", l_min, l_max ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "max iterations", max_iter ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e / %2.1e", "tolerance (CG/MG)", cg_tol, tol ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       " %30s: %d", "additional final step", ( l_final > l_max || max_iter_final > max_iter || tol_final < tol ) ) );
   if ( l_final > l_max || max_iter_final > max_iter || tol_final < tol )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "level (final solve)", l_final ) );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "max iterations (final solve)", max_iter_final ) );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e", "tolerance (final solve)", tol_final ) );
   }
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "write vtk output", ( vtkname != "" ) ) );
   if ( vtkname != "" )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %s", "vtk name", vtkname.c_str() ) );
   }
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "write domain partitioning", writePartitioning ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "export final mesh as file", writeMeshfile ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "use higher level for L2 error", accurate_error ) );

   // solve/refine iteratively
   solve_for_each_refinement( setupStorage,
                              problem,
                              n_refinements,
                              n_el_max,
                              p_refinement,
                              l_min,
                              l_max,
                              l_final,
                              u0,
                              max_iter,
                              tol,
                              max_iter_final,
                              tol_final,
                              cg_tol,
                              accurate_error,
                              vtkname,
                              writePartitioning,
                              writeMeshfile );

   return 0;
}
