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

#include "core/Abort.h"
#include "core/math/KahanSummation.h"

#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p0/p0_linear_form_blending_q5.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_linear_form_blending_q5.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

class Undefined
{};

/// @brief Class representing an L2 space
/// @tparam DiscretizationType   Functiontype defining a discrete subspace
/// @tparam CodomainType         Type defining the codomain of the elements of this space (usually R or R^3)
template < typename DiscretizationType = Undefined, typename CodomainType = real_t >
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
               uint_t                                                                     lvl = DEFAULT ) const
   {
      if ( q == DEFAULT )
         q = _q;
      if ( lvl == DEFAULT )
         lvl = _lvl;

      // inner product of the CodomainType space
      auto innerProduct = [&]( const CodomainType& ux, const CodomainType& vx ) -> real_t {
         if constexpr ( std::is_same_v< CodomainType, real_t > )
         {
            return ux * vx;
         }
         if constexpr ( std::is_same_v< CodomainType, Point3D > )
         {
            return ux.dot( vx );
         }

         WALBERLA_ABORT( "L2 dot product not implemented for CodomainType" )
      };

      real_t localsum = 0.0;

      if ( _storage->hasGlobalCells() )
      {
         // 3D: integrate over all cells
         std::vector< PrimitiveID > cellIDs = _storage->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : localsum )
#endif
         for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
         {
            auto id   = cellIDs[uint_c( i )];
            auto uv   = [&]( const Point3D& x ) { return innerProduct( u( x, id ), v( x, id ) ); };
            auto part = integrate<>( *_storage->getCell( id ), uv, q, lvl );
            uv_T[id]  = part;
            localsum += part;
         }
      }
      else
      {
         // 2D: integrate over all faces
         std::vector< PrimitiveID > faceIDs = _storage->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : localsum )
#endif
         for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
         {
            auto id   = faceIDs[uint_c( i )];
            auto uv   = [&]( const Point3D& x ) { return innerProduct( u( x, id ), v( x, id ) ); };
            auto part = integrate<>( *_storage->getFace( id ), uv, q, lvl );
            uv_T[id]  = part;
            localsum += part;
         }
      }

      // compute global sum
      walberla::mpi::allReduceInplace( localsum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
      return localsum;
   }

   /// @brief Compute b_i = ∫ φ_i f for all basis functions φ_i of the discrete subspace
   /// @param f      L2 function
   /// @param b      output vector to store the values of the integral
   /// @param q      use quadrature of order q instead of the chosen default value
   /// @param lvl    operate on level lvl instead of chosen default value
   void dot( const std::function< CodomainType( const Point3D& ) >& f,
             DiscretizationType&                                    b,
             uint_t                                                 q   = DEFAULT,
             uint_t                                                 lvl = DEFAULT ) const;
   // {
   //    if ( q == DEFAULT )
   //       q = _q;
   //    if ( lvl == DEFAULT )
   //       lvl = _lvl;

   //    if constexpr ( std::is_same_v< CodomainType, real_t > )
   //    {
   //       if constexpr ( std::is_same_v< FE, P1Function< real_t > > )
   //       {
   //          switch ( q )
   //          {
   //          case 5:
   //             return this->applyLinearForm< OpP1 >( P1Q5( f, f ), b, lvl );
   //          }
   //       }
   //       if constexpr ( std::is_same_v< FE, P2Function< real_t > > )
   //       {
   //          // todo: implement for P2 functions, ...
   //       }
   //    }
   //    if constexpr ( std::is_same_v< CodomainType, Point3D > )
   //    {
   //       // todo: implement for vector functions
   //    }
   //    WALBERLA_ABORT( "(u,f)_L2 not implemented for selected combination of CodomainType, FE discretization and quadrature rule" );
   // }

 private:
   /// @brief Integrate a function over a given macro element
   /// @param T            primitive over which f shall be integrated
   /// @param f            integrand
   /// @param q            use quadrature of order q instead of the chosen default value
   /// @param lvl          operate on level lvl instead of chosen default value
   /// @return ∫_T f(x) dx
   template < class PrimitiveType >
   real_t integrate( const PrimitiveType& T, const std::function< real_t( const Point3D& ) >& f, uint_t q, uint_t lvl ) const
   {
      // setup quadrature rule
      std::shared_ptr< P0FormHyTeG > quad;
      switch ( q )
      {
      case 5:
         quad = std::make_shared< forms::p0_linear_form_blending_q5 >( f, f );
         break;

      default:
         WALBERLA_ABORT( "Quadrature rule with q = " << q << " for L2 dot product not implemented!" )
         break;
      }
      quad->setGeometryMap( T.getGeometryMap() );

      // compute the integral
      walberla::math::KahanAccumulator< real_t > integral;       // value of integral over T;
      Matrixr< 1, 1 >                            integral_micro; // value of integral over micro element

      if constexpr ( std::is_same_v< PrimitiveType, Cell > )
      {
         // loop over micro-cells
         for ( const auto& cType : celldof::allCellTypes )
         {
            for ( const auto& micro : celldof::macrocell::Iterator( lvl, cType, 0 ) )
            {
               // computational coordinates of micro cell
               auto                     verts = celldof::macrocell::getMicroVerticesFromMicroCell( micro, cType );
               std::array< Point3D, 4 > coords;
               for ( uint_t k = 0; k < 4; ++k )
               {
                  coords[k] = vertexdof::macrocell::coordinateFromIndex( lvl, T, verts[k] );
               }

               // integral over micro cell
               integral_micro.setZero();
               quad->integrateAll( coords, integral_micro );

               integral += integral_micro( 0, 0 );
            }
         }
      }
      else
      {
         // loop over micro-faces
         for ( const auto& fType : facedof::allFaceTypes )
         {
            for ( const auto& micro : facedof::macroface::Iterator( lvl, fType, 0 ) )
            {
               // computational coordinates of micro face
               auto                     verts = facedof::macroface::getMicroVerticesFromMicroFace( micro, fType );
               std::array< Point3D, 3 > coords;
               for ( uint_t k = 0; k < 3; ++k )
               {
                  coords[k] = vertexdof::macroface::coordinateFromIndex( lvl, T, verts[k] );
               }

               // integral over micro face
               integral_micro.setZero();
               quad->integrateAll( coords, integral_micro );

               integral += integral_micro( 0, 0 );
            }
         }
      }

      return integral.get();
   }

   /// @brief Compute b_i = ∫ φ_i f for all basis functions φ_i of the discrete subspace
   /// @tparam Op          Type of operator fitting for the FE space, e.g P1VariableOperator
   /// @tparam LinearForm  Form defining quadrature rule to evaluate the integrals
   /// @param f      L2 function
   /// @param b      output vector to store the values of the integral
   /// @param lvl    operate on level lvl instead of chosen default value
   template < template < class > class Op, class LinearForm >
   void dot( const std::function< real_t( const Point3D& ) >& f, const DiscretizationType& b, uint_t lvl ) const
   {
      // apply linear form
      LinearForm       form( f, f );
      Op< LinearForm > _b( _storage, lvl, lvl, form );
      _b.computeDiagonalOperatorValues();
      b.copyFrom( *_b.getDiagonalValues(), lvl );

      // free memory of diagonal
      // !this should not be necessary -> maybe fix destructor of FE Function!
      for ( const auto& it : _storage->getVertices() )
      {
         _b.getDiagonalValues()->deleteMemory( lvl, *( it.second ) );
      }
      for ( const auto& it : _storage->getEdges() )
      {
         _b.getDiagonalValues()->deleteMemory( lvl, *( it.second ) );
      }
      for ( const auto& it : _storage->getFaces() )
      {
         _b.getDiagonalValues()->deleteMemory( lvl, *( it.second ) );
      }
      for ( const auto& it : _storage->getCells() )
      {
         _b.getDiagonalValues()->deleteMemory( lvl, *( it.second ) );
      }
   }

   std::shared_ptr< PrimitiveStorage > _storage; // storage corresponding to domain discretization
   uint_t                              _lvl;     // grid level to work on
   uint_t                              _q;       // order of quadrature rule used for computing integrals
};

template <>
void L2Space< P1Function< real_t > >::dot( const std::function< real_t( const Point3D& ) >& f,
                                           P1Function< real_t >&                            b,
                                           uint_t                                           q,
                                           uint_t                                           lvl ) const
{
   if ( q == DEFAULT )
      q = _q;
   if ( lvl == DEFAULT )
      lvl = _lvl;

   using namespace forms;

   switch ( q )
   {
   case 5:
      return dot< P1VariableOperator, p1_linear_form_blending_q5 >( f, b, lvl );
   default:
      WALBERLA_ABORT( "(v_i,f)_L2 not implemented for P1 with selected quadrature rule" );
   }
}

} // namespace hyteg
