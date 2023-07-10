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
/// @tparam Quad                 Order of quadrature rule used to compute integrals on each micro element
/// @tparam DiscretizationType   Functiontype defining a discrete subspace
/// @tparam CodomainType         Type defining the codomain of the elements of this space (usually R or R^3)
template < uint_t Quad, class DiscretizationType = Undefined, typename CodomainType = real_t >
class L2Space
{
 public:
   /// @brief Construct a L2 space corresponding to the domain discretized by storage at a certain grid level
   /// @param storage   PrimitiveStorage object
   /// @param lvl       grid level to work on
   L2Space( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t& lvl )
   : _storage( storage )
   , _lvl( lvl )
   {}

   /// @brief Set working grid level
   /// @param lvl
   void setLvl( const uint_t lvl ) { _lvl = lvl; }

   /// @brief Compute the L2 norm ||u||_0
   /// @param u function in L2
   /// @return ||u||_0
   real_t norm( const std::function< CodomainType( const Point3D& ) >& u ) const { return std::sqrt( this->dot( u, u ) ); }

   /// @brief Compute the L2 norm ||u||_0
   /// @param u      function in L2
   /// @param u2_T   map to store partial result (u,u)_L2(T) on each macro element T
   /// @return ||u||_0
   real_t norm( const std::function< CodomainType( const Point3D& ) >& u, std::map< PrimitiveID, real_t >& u2_T ) const
   {
      return std::sqrt( this->dot( u, u, u2_T ) );
   }

   /// @brief Compute the L2 norm ||u||_0
   /// @param u function in L2
   /// @return ||u||_0
   real_t norm( const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& u ) const
   {
      return std::sqrt( this->dot( u, u ) );
   }

   /// @brief Compute the L2 norm ||u||_0
   /// @param u function in L2
   /// @param u2_T   map to store partial result (u,u)_L2(T) on each macro element T
   /// @return ||u||_0
   real_t norm( const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& u,
                std::map< PrimitiveID, real_t >&                                           u2_T ) const
   {
      return std::sqrt( this->dot( u, u, u2_T ) );
   }

   /// @brief Compute an approximation to the L2 dot product (u,v)_0
   /// @param u      function in L2
   /// @param v      function in L2
   /// @return (u,v)_0
   real_t dot( const std::function< CodomainType( const Point3D& ) >& u,
               const std::function< CodomainType( const Point3D& ) >& v ) const
   {
      std::map< PrimitiveID, real_t > _;
      return this->dot( u, v, _ );
   }

   /// @brief Compute an approximation to the L2 dot product (u,v)_0
   /// @param u      function in L2
   /// @param v      function in L2
   /// @param uv_T   map to store partial result of dot product on each macro element T
   /// @return (u,v)_0
   real_t dot( const std::function< CodomainType( const Point3D& ) >& u,
               const std::function< CodomainType( const Point3D& ) >& v,
               std::map< PrimitiveID, real_t >&                       uv_T ) const
   {
      auto uid = [&]( const Point3D& x, const PrimitiveID& ) { return u( x ); };
      auto vid = [&]( const Point3D& x, const PrimitiveID& ) { return v( x ); };
      return this->dot( uid, vid, uv_T );
   }

   /// @brief Compute an approximation to the L2 dot product (u,v)_0
   /// @param u      function in L2
   /// @param v      function in L2
   /// @return (u,v)_0
   real_t dot( const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& u,
               const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& v ) const
   {
      std::map< PrimitiveID, real_t > _;
      return this->dot( u, v, _ );
   }

   /// @brief Compute an approximation to the L2 dot product (u,v)_0
   /// @param u      function in L2
   /// @param v      function in L2
   /// @param uv_T   map to store partial result of dot product on each macro element T
   /// @return (u,v)_0
   real_t dot( const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& u,
               const std::function< CodomainType( const Point3D&, const PrimitiveID& ) >& v,
               std::map< PrimitiveID, real_t >&                                           uv_T ) const;

   /// @brief Compute b_i = ∫ φ_i f for all basis functions φ_i of the discrete subspace
   /// @param f      L2 function
   /// @param b      output vector to store the values of the integral
   void dot( const std::function< CodomainType( const Point3D& ) >& f, DiscretizationType& b ) const;

 private:
   /// @brief Integrate a function over a given macro element
   /// @param primitive over which f shall be integrated
   /// @param f  integrand
   /// @return ∫_T f(x) dx
   template < class PrimitiveType >
   real_t integrate( const PrimitiveType& primitive, const std::function< real_t( const Point3D& ) >& f ) const;

   /// @brief Integrate a function over a given macro element
   /// @tparam QuadratureRule
   /// @tparam PrimitiveType
   /// @param primitive  primitive over which f shall be integrated
   /// @param f  integrand
   /// @return ∫_T f(x) dx
   template < class QuadratureRule, class PrimitiveType >
   real_t integrate( const PrimitiveType& primitive, const std::function< real_t( const Point3D& ) >& f ) const;

   /// @brief Compute b_i = ∫ φ_i f for all basis functions φ_i of the discrete subspace
   /// @tparam Op          Type of operator fitting for the FE space, e.g P1VariableOperator
   /// @tparam LinearForm  Form defining quadrature rule to evaluate the integrals
   /// @param f      L2 function
   /// @param b      output vector to store the values of the integral
   template < template < class > class Op, class LinearForm >
   void dot( const std::function< real_t( const Point3D& ) >& f, const DiscretizationType& b ) const;

   std::shared_ptr< PrimitiveStorage > _storage; // storage corresponding to domain discretization
   uint_t                              _lvl;     // grid level to work on
};

} // namespace hyteg
