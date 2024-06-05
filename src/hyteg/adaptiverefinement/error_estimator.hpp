/*
 * Copyright (c) 2024 Benjamin Mann
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

#pragma once

#include <limits>

#include "hyteg_operators/operators/mass/P1ElementwiseMass.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

#include "mesh.hpp"

namespace hyteg {
namespace adaptiveRefinement {

/// @brief ErrorEstimator based on CITATION.
///         ! Requires using FMG as solver !
template < class FE_function >
class ErrorEstimator
{
   using valueType = typename FE_function::valueType;

 public:
   enum Norm
   {
      L2
      // todo add other norms
   };

   /// @brief ErrorEstimator based on CITATION
   /// @param u solution vector the estimator shall correspond to
   /// @param j_max eta_j will be computed for all j=1,...,j_max.
   ///              If j_max=0, only local indicator eta_T will be computed
   /// @param norm norm the estimate is based on
   ErrorEstimator( const FE_function& u, uint_t j_max, Norm norm = L2 )
   : _storage( u.getStorage() )
   , _l_min( u.getMinLevel() )
   , _l_max( u.getMaxLevel() )
   , _j_max( j_max )
   , _l_min_fmg( _l_min )
   , _l_max_fmg( _l_max )
   , _u( u )
   , _eta( j_max + 1 )
   , _theta( j_max )
   , _C1( j_max )
   , _C2( j_max )
   , _fmg_called( false )
   , _estimate_called( false )
   {
      // we require one more level than L-j to compute η_j
      if ( _l_max - _l_min < j_max + 1 )
      {
         WALBERLA_ABORT( "ErrorEstimator requires (j_max + 1) levels of MG hierarchy!" )
      }
      else
      {
         _l_min = _l_max - _j_max - 1;
      }

      _w  = std::make_unique< FE_function >( "w", _storage, _l_min, _l_max );
      _Mw = std::make_unique< FE_function >( "Mw", _storage, _l_max, _l_max );

      if constexpr ( std::is_same_v< FE_function, P1Function< valueType > > )
      {
         _P = std::make_unique< P1toP1LinearProlongation<> >();
         if ( norm == L2 )
         {
            _M = std::make_unique< operatorgeneration::P1ElementwiseMass >( _storage, _l_max, _l_max );
         }
      }
      else
      {
         // todo: implement P2, etc.
         WALBERLA_ABORT( "ErrorEstimator not implemented for desired FE space" )
      }
   }

   /// @brief postProlongateCallback required to be used in FullMultigridSolver
   /// @param lvl
   void fmg_callback( uint_t lvl )
   {
      // reset flag
      _estimate_called = false;
      _fmg_called      = true;

      // copy P * u[lvl] to w[lvl+1]
      if ( _l_min <= lvl && lvl < _l_max )
      {
         _w->copyFrom( _u, lvl + 1 );
      }

      // store fmg levels to check later
      if ( !_fmg_called )
      {
         _l_min_fmg = lvl;
      }
      _l_max_fmg = lvl;
   }

   /// @brief trigger computation of error estimates. Call AFTER fmg solver.
   void estimate()
   {
      if ( !_fmg_called )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "No data from FMG available. Maybe ErrorEstimator::fmg_callback() hasn't been used." )
         return;
      }
      if ( _l_min_fmg > _l_min )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Not enough coarse solutions from FMG available." )
         return;
      }
      if ( _l_max_fmg < _l_max )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "No fine solution at u.getMaxLevel() available from FMG." )
         return;
      }
      for ( uint_t j = 1; j <= _j_max + 1; ++j )
      {
         auto lvl = _l_max - j; // we require all levels L-j_max-1, ..., L-1

         // prolongate u_lvl to finest level (the first prolongations happens within FMG)
         for ( uint_t k = lvl + 1; k < _l_max; ++k )
            _P->prolongate( *_w, k, All );
         // substract fine grid solution to obtain estimate for e_lvl
         _w->add( { -1 }, { _u }, _l_max, All );
         _w->interpolate( 0, _l_max, DirichletBoundary );
         // (after M.apply() (ELEMENTWISE operator), volume-primitives only contain local contributions. Additive interface data only stored on interfaces)
         _M->apply( *_w, *_Mw, _l_max, Inner | NeumannBoundary, Replace );

         // compute squared global and local norms w*M*w
         if ( _storage->hasGlobalCells() ) // 3d
         {
            auto& e  = _w->getCellDataID();
            auto& Me = _Mw->getCellDataID();
            for ( auto& [id, cell] : _storage->getCells() )
            {
               real_t eMe = 0.0;
               if constexpr ( std::is_same_v< FE_function, P1Function< valueType > > )
               {
                  eMe = vertexdof::macrocell::dot< valueType >( _l_max, *cell, e, Me, 0 );
               }
               // todo implement for P2, etc.
               _eta[j - 1] += eMe;
               if ( j == 1 )
                  _eta_T_sq.push_back( { eMe, id } );
            }
         }
         else // 2d
         {
            auto& e  = _w->getFaceDataID();
            auto& Me = _Mw->getFaceDataID();
            for ( auto& [id, face] : _storage->getFaces() )
            {
               real_t eMe = 0.0;
               if constexpr ( std::is_same_v< FE_function, P1Function< valueType > > )
               {
                  eMe = vertexdof::macroface::dot< valueType >( _l_max, *face, e, Me, 0 );
               }
               // todo implement for P2, etc.
               _eta[j - 1] += eMe;
               if ( j == 1 )
                  _eta_T_sq.push_back( { eMe, id } );
            }
         }
      }

      _fmg_called      = false;
      _estimate_called = true;

      if ( _j_max > 0 ) // if only local estimate needed, we omit global communication
      {
         compute_eta();
      }
   }

   /// @brief get global estimate eta_j to ||u - u_h||
   /// @param j
   /// @return eta_j = (theta_j)^j * ||e_L - e_(L-j)||
   valueType eta( uint_t j ) const
   {
      if ( check_j_exists( j ) )
      {
         return _eta[j - 1];
      }
      else
      {
         return std::numeric_limits< valueType >::quiet_NaN();
      }
   }

   /// @brief get convergence estimate theta_j
   /// @param j
   /// @return theta_j = ||e_L - e_(L-j)|| / ||e_L - e_(L-j-1)||
   valueType theta( uint_t j ) const
   {
      if ( check_j_exists( j ) )
      {
         return _theta[j - 1];
      }
      else
      {
         return std::numeric_limits< valueType >::quiet_NaN();
      }
   }

   /// @brief get list of squared error estimates for each macro cell
   /// @return [(eta_T)^2, id_T] for all T owned by this process
   const ErrorVector& eta_T_sq() const
   {
      if ( !_estimate_called )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "No estimates have been computed!" );
      }
      return _eta_T_sq;
   }

 private:
   /// @brief trigger computation of eta, theta, C1, C2
   void compute_eta()
   {
      for ( auto& e : _eta )
      {
         walberla::mpi::allReduceInplace( e, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
         e = std::sqrt( e );
      }

      for ( uint_t i = 0; i < _j_max; ++i )
      {
         // j=i+1 (zero indexed arrays)
         _theta[i] = _eta[i] / _eta[i + 1];
      }
      for ( uint_t i = 0; i < _j_max; ++i )
      {
         // j=i+1 (zero indexed arrays)

         auto theta_p_j  = std::pow( _theta[i], i + 1 ); // θ^(j)
         auto theta_p_j1 = theta_p_j * _theta[i];        // θ^(j+1)

         // error estimate
         _eta[i] = _eta[i] * theta_p_j;

         // error constants
         _C1[i] = std::pow( 1.0 - theta_p_j, i + 2 ) / std::pow( 1.0 + theta_p_j1, i + 1 );
         _C2[i] = std::pow( 1.0 + theta_p_j, i + 2 ) / std::pow( 1.0 - theta_p_j1, i + 1 );
      }
   }

   /// @brief Check whether values for level j have been computed
   /// @param j
   /// @return true if values exist, false otherwise
   bool check_j_exists( uint_t j ) const
   {
      if ( j < 1 || j > _j_max || !_estimate_called )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( walberla::format( "Estimate not computed for level j=%d!", j ) );
         return false;
      }
      return true;
   }

   std::shared_ptr< PrimitiveStorage >                     _storage;
   uint_t                                                  _l_min, _l_max, _j_max, _l_min_fmg, _l_max_fmg;
   const FE_function&                                      _u;      // ref to solution vector
   std::unique_ptr< FE_function >                          _w, _Mw; // vectors for temporary values
   std::unique_ptr< Operator< FE_function, FE_function > > _M;      // spd operator for norm
   std::unique_ptr< ProlongationOperator< FE_function > >  _P;
   std::vector< valueType >                                _eta, _theta, _C1, _C2;
   adaptiveRefinement::ErrorVector                         _eta_T_sq;
   bool                                                    _fmg_called, _estimate_called;
};

} // namespace adaptiveRefinement
} // namespace hyteg
