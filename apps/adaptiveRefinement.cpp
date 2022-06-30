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
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

#define R_min 1.0
#define R_max 2.0

using Mass         = P1ConstantMassOperator;
using Laplace      = P1ConstantLaplaceOperator;
using DivkGradForm = forms::p1_div_k_grad_affine_q3;
using DivkGrad     = P1VariableOperator< DivkGradForm >;

// representation of a model problem
struct ModelProblem
{
   enum Type
   {
      DIRAC,
      DIRAC_REGULARIZED,
      SLIT,
      REENTRANT_CORNER,
      NOT_AVAILABLE
   };

   ModelProblem( uint_t mp_type, uint_t spacial_dim )
   : type( Type( mp_type ) )
   , dim( spacial_dim )
   , sigma( 0.0 )
   , offset_err( 0.0 )
   , ignore_offset_for_refinement( false )
   {
      if ( mp_type >= Type::NOT_AVAILABLE )
      {
         WALBERLA_ABORT( "Invalid argument for model problem: Must be less than " << Type::NOT_AVAILABLE );
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

      default:
         break;
      }

      ss << "_" << dim << "D";

      return ss.str();
   }

   void set_params( real_t sig, real_t offset, bool refine_all )
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
      // todo: add parameters when adding new problems
   }

   bool constant_coefficient() const
   {
      // todo: adapt when adding setting with non-constant coefficient
      return type < NOT_AVAILABLE;
   }

   auto coeff() const
   {
      std::function< real_t( const Point3D& ) > k;

      if ( constant_coefficient() )
      {
         k = [=]( const Point3D& ) { return 1.0; };
      }
      // todo: add required coefficients

      return k;
   }

   void init( std::shared_ptr< PrimitiveStorage > storage,
              P1Function< real_t >&               u,
              P1Function< real_t >&               f,
              P1Function< real_t >&               u_h,
              P1Function< real_t >&               b,
              uint_t                              lvl ) const
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
            _f = [=]( const Point3D& x ) -> real_t { return std::exp( -0.5 * x.normSq() / sig2 ) / ( sqrt_2pi_pow_3 * sig3 ); };
            _u = [=]( const Point3D& x ) -> real_t {
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
      u.interpolate( _u, lvl + 1 );

      // initialize u_h
      u_h.setToZero( lvl );
      u_h.setToZero( lvl + 1 );
      u_h.interpolate( _u, lvl, DirichletBoundary );
      u_h.interpolate( _u, lvl + 1, DirichletBoundary );
   }

   void apply_error_filter( const P1Function< real_t >& e, P1Function< real_t >& e_f, uint_t lvl ) const
   {
      e.communicate< Vertex, Edge >( lvl );
      e.communicate< Edge, Face >( lvl );
      e.communicate< Face, Cell >( lvl );

      auto filter = [=]( const Point3D& x, const std::vector< real_t >& val ) -> real_t {
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
};

SetupPrimitiveStorage domain( const ModelProblem& problem, uint_t N )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( problem.dim == 2 )
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
      else
      {
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
                                       std::string                              vtkname,
                                       bool                                     writePartitioning,
                                       uint_t                                   refinement_step,
                                       bool                                     l2_error_each_iteration = true )
{
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "* create PrimitiveStorage ..." );
   auto storage = mesh.make_storage();
   printCurrentMemoryUsage();

   WALBERLA_LOG_INFO_ON_ROOT( "* solve system ..." );

   // operators
   std::shared_ptr< A_t > A;
   auto                   M = std::make_shared< Mass >( storage, l_min, l_max + 1 );
   auto                   P = std::make_shared< P1toP1LinearProlongation >();
   auto                   R = std::make_shared< P1toP1LinearRestriction >();

   if constexpr ( std::is_same_v< A_t, DivkGrad > )
   {
      auto k = problem.coeff();

      A = std::make_shared< A_t >( storage, l_min, l_max + 1, DivkGradForm( k, k ) );
   }
   else
   {
      A = std::make_shared< A_t >( storage, l_min, l_max + 1 );
   }

   // FE functions
   auto                 u = std::make_shared< P1Function< real_t > >( "u_h", storage, l_min, l_max + 1 );
   P1Function< real_t > f( "f", storage, l_max, l_max );
   P1Function< real_t > b( "b(f)", storage, l_min, l_max );
   P1Function< real_t > r( "r", storage, l_max, l_max );
   P1Function< real_t > u_anal( "u_anal", storage, l_max, l_max + 1 );
   P1Function< real_t > err( "e_h", storage, l_max, l_max + 1 );
   P1Function< real_t > tmp( "tmp", storage, l_min, l_max + 1 );
   P1Function< real_t > err_f( "e_h_filtered", storage, l_max + 1, l_max + 1 );

   // global DoF
   tmp.interpolate( []( const Point3D& ) { return 1.0; }, l_max, Inner | NeumannBoundary );
   auto n_dof = uint_t( tmp.dotGlobal( tmp, l_max ) );
   WALBERLA_LOG_INFO_ON_ROOT( " -> number of global DoF: " << n_dof );

   // computation of residual and L2 error
   err.setToZero( l_max );
   err.setToZero( l_max + 1 );
   tmp.setToZero( l_max );
   tmp.setToZero( l_max + 1 );
   auto compute_residual = [&]() -> real_t {
      A->apply( *u, tmp, l_max, Inner | NeumannBoundary, Replace );
      r.assign( { 1.0, -1.0 }, { b, tmp }, l_max, Inner | NeumannBoundary );
      M->apply( r, tmp, l_max, Inner | NeumannBoundary, Replace );
      return std::sqrt( r.dotGlobal( tmp, l_max, Inner | NeumannBoundary ) );
   };
   auto compute_L2error = [&]() -> real_t {
      P->prolongate( *u, l_max, Inner | NeumannBoundary );
      err.assign( { 1.0, -1.0 }, { *u, u_anal }, l_max + 1, Inner | NeumannBoundary );
      problem.apply_error_filter( err, err_f, l_max + 1 );
      M->apply( err_f, tmp, l_max + 1, Inner | NeumannBoundary, Replace );
      return std::sqrt( err_f.dotGlobal( tmp, l_max + 1, Inner | NeumannBoundary ) );
   };

   // initialize analytic solution, rhs and boundary values
   problem.init( storage, u_anal, f, *u, b, l_max );

   // initialize u_h
   if ( u0 == 1 && u_old != nullptr )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " -> initialize u=u_old" );
      auto u_init = [&]( const Point3D& x ) -> real_t {
         real_t ux = 0.0;
         if ( !u_old->evaluate( x, l_interpolate, ux, 1e-12 ) )
         {
            WALBERLA_LOG_INFO( "     interpolation failed at x = " << x );
         }
         return ux;
      };
      u->interpolate( u_init, l_max, Inner );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( " -> initialize u=0" );
   }

   // apply loadbalancing
   auto migrationInfo = mesh.loadbalancing();
   storage->migratePrimitives( migrationInfo );

   // initial residual
   real_t norm_r = compute_residual();

   // solver
   auto coarseGridSolver = std::make_shared< CGSolver< A_t > >( storage, l_min, l_min, std::max( max_iter, n_dof ), tol * 1e-1 );
   auto smoother         = std::make_shared< GaussSeidelSmoother< A_t > >();

   GeometricMultigridSolver< A_t > gmg( storage, smoother, coarseGridSolver, R, P, l_min, l_max, 3, 3 );

   // solve
   WALBERLA_LOG_INFO_ON_ROOT( " -> run multigrid solver" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6s |%18s |%13s", "k", "||r_k||_L2", "||e_k||_L2" ) );

   uint_t iter = 0;
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
      gmg.solve( *A, *u, b, l_max );
      norm_r = compute_residual();
   }

   auto eL2 = compute_L2error();
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6d |%18.3e |%13.3e", iter, norm_r, eL2 ) );

   // compute elementwise error
   adaptiveRefinement::ErrorVector err_2_elwise_loc;

   P1Function< real_t >& err_for_refinement = ( problem.ignore_offset_for_refinement ) ? err : err_f;

   if ( problem.dim == 3 )
   {
      for ( auto& [id, cell] : storage->getCells() )
      {
         real_t err_2_cell = vertexdof::macrocell::dot< real_t >(
             l_max + 1, *cell, err_for_refinement.getCellDataID(), err_for_refinement.getCellDataID(), 0 );

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
             l_max + 1, *face, err_for_refinement.getFaceDataID(), err_for_refinement.getFaceDataID(), 0 );

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

   // export to vtk
   if ( vtkname != "" )
   {
      // coefficient
      P1Function< real_t > k( "k", storage, l_max, l_max );
      k.interpolate( problem.coeff(), l_max );

      // error
      R->restrict( err, l_max + 1, Inner | NeumannBoundary );

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

   // update u_old
   if ( u0 == 1 )
   {
      u_old.swap( u );
   }

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
                                std::string                  vtkname,
                                bool                         writePartitioning )
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
         local_errors = solve< Laplace >(
             mesh, problem, u0, u_old, l_max, l_min, l_max, max_iter, tol, vtkname, writePartitioning, refinement );
      else
         local_errors = solve< DivkGrad >(
             mesh, problem, u0, u_old, l_max, l_min, l_max, max_iter, tol, vtkname, writePartitioning, refinement );

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
      auto ratio = mesh.refineRG( local_errors, criterion, n_el_max );
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
         solve< Laplace >(
             mesh, problem, u0, u_old, l_max, l_min, l_final, max_iter_final, tol_final, vtkname, writePartitioning, refinement );
      }
      else
      {
         solve< DivkGrad >(
             mesh, problem, u0, u_old, l_max, l_min, l_final, max_iter_final, tol_final, vtkname, writePartitioning, refinement );
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
   WALBERLA_LOG_INFO_ON_ROOT( "config = " << *cfg );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t mp         = parameters.getParameter< uint_t >( "modelProblem" );
   const uint_t dim        = parameters.getParameter< uint_t >( "dim" );
   const real_t sigma      = parameters.getParameter< real_t >( "sigma", 0.0 );
   const real_t offset     = parameters.getParameter< real_t >( "offset", 0.0 );
   const bool   refine_all = parameters.getParameter< bool >( "refine_all", false );
   // todo: add parameters when adding new problems

   const uint_t N             = parameters.getParameter< uint_t >( "initial_resolution", 1 );
   const uint_t n_refinements = parameters.getParameter< uint_t >( "n_refinements" );
   const uint_t n_el_max      = parameters.getParameter< uint_t >( "n_el_max", std::numeric_limits< uint_t >::max() );
   const real_t p_refinement  = parameters.getParameter< real_t >( "percentile" );

   const uint_t l_min   = 0;
   const uint_t l_max   = parameters.getParameter< uint_t >( "microlevel" );
   const uint_t l_final = parameters.getParameter< uint_t >( "microlevel_final", l_max );

   const uint_t u0             = parameters.getParameter< uint_t >( "initial_guess", 0 );
   const uint_t max_iter       = parameters.getParameter< uint_t >( "n_iterations" );
   const uint_t max_iter_final = parameters.getParameter< uint_t >( "n_iterations_final", max_iter );
   const real_t tol            = parameters.getParameter< real_t >( "tolerance" );
   const real_t tol_final      = parameters.getParameter< real_t >( "tolerance_final", tol );

   std::string vtkname           = parameters.getParameter< std::string >( "vtkName", "" );
   const bool  writePartitioning = parameters.getParameter< bool >( "writeDomainPartitioning", false );

   // setup model problem
   ModelProblem problem( mp, dim );
   WALBERLA_LOG_INFO_ON_ROOT( "Model Problem: " << problem.name() );
   problem.set_params( sigma, offset, refine_all );
   // setup domain
   auto setupStorage = domain( problem, N );
   // if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
   //    setupStorage.toStream( std::cout, true );
   // solve/refine iteratively
   if ( vtkname == "auto" )
   {
      vtkname = problem.name();
   }
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
                              vtkname,
                              writePartitioning );

   return 0;
}
