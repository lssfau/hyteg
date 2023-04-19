/*
 * Copyright (c) 2023 Benjamin Mann.
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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

class Undefined
{};

/// @brief Class representing an L2 space
/// @tparam DiscretizationType   Functiontype defining a discrete subspace
/// @tparam CodomainType         Type defining the codomain of the elements of this space (usually R or R^3)
template < class DiscretizationType = Undefined, typename CodomainType = real_t >
class L2Space
{
 public:
   static constexpr uint_t DEFAULT = 0xdef;

   /// @brief Construct a L2 space corresponding to the domain discretized by storage at a certain grid level
   /// @param storage   PrimitiveStorage object
   /// @param lvl       Default grid level to work on
   /// @param q         Default order of quadrature rule to use for computing integrals.
   L2Space( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t& lvl, const uint_t& q )
   : _storage( storage )
   , _lvl( lvl )
   , _q( q )
   {}

   /// @brief Set default order of quadrature rule
   /// @param q
   void setQuad( const uint_t q ) { _q = q; }

   /// @brief Set default working grid level
   /// @param lvl
   void setLvl( const uint_t lvl ) { _lvl = lvl; }

   /// @brief Compute the L2 norm ||u||_0
   /// @param u function in L2
   /// @param q use quadrature of order q instead of the chosen default value
   /// @param lvl operate on level lvl instead of chosen default value
   /// @return ||u||_0
   real_t norm( const std::function< CodomainType( const Point3D& ) >& u, uint_t q = DEFAULT, uint_t lvl = DEFAULT ) const
   {
      return std::sqrt( this->dot( u, u, q, lvl ) );
   }

   /// @brief Compute the L2 norm ||u||_0
   /// @param u      function in L2
   /// @param u2_T   map to store partial result (u,u)_L2(T) on each macro element T
   /// @param q      use quadrature of order q instead of the chosen default value
   /// @param lvl    operate on level lvl instead of chosen default value
   /// @return ||u||_0
   real_t norm( const std::function< CodomainType( const Point3D& ) >& u,
                std::map< PrimitiveID, real_t >&                       u2_T,
                uint_t                                                 q   = DEFAULT,
                uint_t                                                 lvl = DEFAULT ) const
   {
      return std::sqrt( this->dot( u, u, u2_T, q, lvl ) );
   }

   /// @brief Compute the L2 norm ||u||_0
   /// @param u function in L2
   /// @param q use quadrature of order q instead of the chosen default value
   /// @param lvl operate on level lvl instead of chosen default value
   /// @return ||u||_0
   real_t norm( const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& u,
                uint_t                                                                     q   = DEFAULT,
                uint_t                                                                     lvl = DEFAULT ) const
   {
      return std::sqrt( this->dot( u, u, q, lvl ) );
   }

   /// @brief Compute the L2 norm ||u||_0
   /// @param u function in L2
   /// @param u2_T   map to store partial result (u,u)_L2(T) on each macro element T
   /// @param q use quadrature of order q instead of the chosen default value
   /// @param lvl operate on level lvl instead of chosen default value
   /// @return ||u||_0
   real_t norm( const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& u,
                std::map< PrimitiveID, real_t >&                                           u2_T,
                uint_t                                                                     q   = DEFAULT,
                uint_t                                                                     lvl = DEFAULT ) const
   {
      return std::sqrt( this->dot( u, u, u2_T, q, lvl ) );
   }

   /// @brief Compute an approximation to the L2 dot product (u,v)_0
   /// @param u      function in L2
   /// @param v      function in L2
   /// @param q      use quadrature of order q instead of the chosen default value
   /// @param lvl    operate on level lvl instead of chosen default value
   /// @return (u,v)_0
   real_t dot( const std::function< CodomainType( const Point3D& ) >& u,
               const std::function< CodomainType( const Point3D& ) >& v,
               uint_t                                                 q   = DEFAULT,
               uint_t                                                 lvl = DEFAULT ) const
   {
      std::map< PrimitiveID, real_t > _;
      return this->dot( u, v, _, q, lvl );
   }

   /// @brief Compute an approximation to the L2 dot product (u,v)_0
   /// @param u      function in L2
   /// @param v      function in L2
   /// @param uv_T   map to store partial result of dot product on each macro element T
   /// @param q      use quadrature of order q instead of the chosen default value
   /// @param lvl    operate on level lvl instead of chosen default value
   /// @return (u,v)_0
   real_t dot( const std::function< CodomainType( const Point3D& ) >& u,
               const std::function< CodomainType( const Point3D& ) >& v,
               std::map< PrimitiveID, real_t >&                       uv_T,
               uint_t                                                 q   = DEFAULT,
               uint_t                                                 lvl = DEFAULT ) const
   {
      auto uid = [&]( const Point3D& x, const PrimitiveID& ) { return u( x ); };
      auto vid = [&]( const Point3D& x, const PrimitiveID& ) { return v( x ); };
      return this->dot( uid, vid, uv_T, q, lvl );
   }

   /// @brief Compute an approximation to the L2 dot product (u,v)_0
   /// @param u      function in L2
   /// @param v      function in L2
   /// @param q      use quadrature of order q instead of the chosen default value
   /// @param lvl    operate on level lvl instead of chosen default value
   /// @return (u,v)_0
   real_t dot( const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& u,
               const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& v,
               uint_t                                                                     q   = DEFAULT,
               uint_t                                                                     lvl = DEFAULT ) const
   {
      std::map< PrimitiveID, real_t > _;
      return this->dot( u, v, _, q, lvl );
   }

   /// @brief Compute an approximation to the L2 dot product (u,v)_0
   /// @param u      function in L2
   /// @param v      function in L2
   /// @param uv_T   map to store partial result of dot product on each macro element T
   /// @param q      use quadrature of order q instead of the chosen default value
   /// @param lvl    operate on level lvl instead of chosen default value
   /// @return (u,v)_0
   real_t dot( const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& u,
               const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& v,
               std::map< PrimitiveID, real_t >&                                           uv_T,
               uint_t                                                                     q   = DEFAULT,
               uint_t                                                                     lvl = DEFAULT ) const;

   /// @brief Compute b_i = ∫ φ_i f for all basis functions φ_i of the discrete subspace
   /// @param f      L2 function
   /// @param b      output vector to store the values of the integral
   /// @param q      use quadrature of order q instead of the chosen default value
   /// @param lvl    operate on level lvl instead of chosen default value
   void dot( const std::function< CodomainType( const Point3D& ) >& f,
             DiscretizationType&                                    b,
             uint_t                                                 q   = DEFAULT,
             uint_t                                                 lvl = DEFAULT ) const;

 private:
   /// @brief Integrate a function over a given macro element
   /// @param T            primitive over which f shall be integrated
   /// @param f            integrand
   /// @param q            use quadrature of order q instead of the chosen default value
   /// @param lvl          operate on level lvl instead of chosen default value
   /// @return ∫_T f(x) dx
   template < class PrimitiveType >
   real_t integrate( const PrimitiveType& T, const std::function< real_t( const Point3D& ) >& f, uint_t q, uint_t lvl ) const;

   /// @brief Compute b_i = ∫ φ_i f for all basis functions φ_i of the discrete subspace
   /// @tparam Op          Type of operator fitting for the FE space, e.g P1VariableOperator
   /// @tparam LinearForm  Form defining quadrature rule to evaluate the integrals
   /// @param f      L2 function
   /// @param b      output vector to store the values of the integral
   /// @param lvl    operate on level lvl instead of chosen default value
   template < template < class > class Op, class LinearForm >
   void dot( const std::function< real_t( const Point3D& ) >& f, const DiscretizationType& b, uint_t lvl ) const;

   std::shared_ptr< PrimitiveStorage > _storage; // storage corresponding to domain discretization
   uint_t                              _lvl;     // grid level to work on
   uint_t                              _q;       // order of quadrature rule used for computing integrals
};

} // namespace hyteg
