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

#include "hyteg/numerictools/L2Space.hpp"

#include "core/Abort.h"
#include "core/math/KahanSummation.h"

#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p0/p0_linear_form_blending_q5.hpp"
#include "hyteg/forms/form_hyteg_generated/p0/p0_linear_form_blending_q7.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_linear_form_blending_q5.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_linear_form_blending_q7.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

namespace hyteg {

template < uint_t Quad, class Discr, typename Codomain >
real_t L2Space< Quad, Discr, Codomain >::dot( const std::function< Codomain( const Point3D&, const PrimitiveID& ) >& u,
                                              const std::function< Codomain( const Point3D&, const PrimitiveID& ) >& v,
                                              std::map< PrimitiveID, real_t >&                                       uv_T ) const
{
   // inner product of the Codomain space, u(x)⋅v(x)
   auto innerProduct = [&]( const Codomain& ux, const Codomain& vx ) -> real_t {
      if constexpr ( std::is_same_v< Codomain, real_t > )
      {
         return ux * vx;
      }
      if constexpr ( std::is_same_v< Codomain, Point3D > )
      {
         return ux.dot( vx );
      }

      WALBERLA_ABORT( "L2 dot product not implemented for Codomain" )
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
         auto part = integrate( *_storage->getCell( id ), uv );
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
         auto part = integrate( *_storage->getFace( id ), uv );
         uv_T[id]  = part;
         localsum += part;
      }
   }

   // compute global sum
   walberla::mpi::allReduceInplace( localsum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   return localsum;
}

template < uint_t Quad, class Discr, typename Codomain >
template < class QuadratureRule, class PrimitiveType >
real_t L2Space< Quad, Discr, Codomain >::integrate( const PrimitiveType&                             T,
                                                    const std::function< real_t( const Point3D& ) >& f ) const
{
   QuadratureRule q( f, f );
   q.setGeometryMap( T.getGeometryMap() );

   // compute the integral
   walberla::math::KahanAccumulator< real_t > integral;       // value of integral over T;
   Matrixr< 1, 1 >                            integral_micro; // value of integral over micro element

   if constexpr ( std::is_same_v< PrimitiveType, Cell > )
   {
      // loop over micro-cells
      for ( const auto& cType : celldof::allCellTypes )
      {
         for ( const auto& micro : celldof::macrocell::Iterator( _lvl, cType, 0 ) )
         {
            // computational coordinates of micro cell
            auto                     verts = celldof::macrocell::getMicroVerticesFromMicroCell( micro, cType );
            std::array< Point3D, 4 > coords;
            for ( uint_t k = 0; k < 4; ++k )
            {
               coords[k] = vertexdof::macrocell::coordinateFromIndex( _lvl, T, verts[k] );
            }

            // integral over micro cell
            integral_micro.setZero();
            q.integrateAll( coords, integral_micro );

            integral += integral_micro( 0, 0 );
         }
      }
   }
   else
   {
      // loop over micro-faces
      for ( const auto& fType : facedof::allFaceTypes )
      {
         for ( const auto& micro : facedof::macroface::Iterator( _lvl, fType, 0 ) )
         {
            // computational coordinates of micro face
            auto                     verts = facedof::macroface::getMicroVerticesFromMicroFace( micro, fType );
            std::array< Point3D, 3 > coords;
            for ( uint_t k = 0; k < 3; ++k )
            {
               coords[k] = vertexdof::macroface::coordinateFromIndex( _lvl, T, verts[k] );
            }

            // integral over micro face
            integral_micro.setZero();
            q.integrateAll( coords, integral_micro );

            integral += integral_micro( 0, 0 );
         }
      }
   }

   return integral.get();
}

template < uint_t Quad, class Discr, typename Codomain >
template < template < class > class Op, class LinearForm >
void L2Space< Quad, Discr, Codomain >::dot( const std::function< real_t( const Point3D& ) >& f, const Discr& b ) const
{
   // apply linear form
   LinearForm       form( f, f );
   Op< LinearForm > _b( _storage, _lvl, _lvl, form );
   _b.computeDiagonalOperatorValues();
   b.copyFrom( *_b.getDiagonalValues(), _lvl );

   // free memory of diagonal
   // !this should not be necessary -> maybe fix destructor of FE Function!
   for ( const auto& it : _storage->getVertices() )
   {
      _b.getDiagonalValues()->deleteMemory( _lvl, *( it.second ) );
   }
   for ( const auto& it : _storage->getEdges() )
   {
      _b.getDiagonalValues()->deleteMemory( _lvl, *( it.second ) );
   }
   for ( const auto& it : _storage->getFaces() )
   {
      _b.getDiagonalValues()->deleteMemory( _lvl, *( it.second ) );
   }
   for ( const auto& it : _storage->getCells() )
   {
      _b.getDiagonalValues()->deleteMemory( _lvl, *( it.second ) );
   }
}

// === template specializations ===

//!!! WHEN ADDING NEW QUADRATURE RULES OR DISCRETIZATIONS ALWAYS
// * add an appropriate P0 quadrature rule to integrate()
// * implement the corresponding specialization of L2Space<?,?,?>::dot()
//!!!

// P0 quadratre rules for integrating arbitrary functions over the domain
template < uint_t Quad, class Discr, typename Codomain >
template < class PrimitiveType >
real_t L2Space< Quad, Discr, Codomain >::integrate( const PrimitiveType&                             T,
                                                    const std::function< real_t( const Point3D& ) >& f ) const
{
   switch ( Quad )
   {
   case 5:
      return integrate< forms::p0_linear_form_blending_q5 >( T, f );
   case 7:
      return integrate< forms::p0_linear_form_blending_q7 >( T, f );
   default:
      WALBERLA_ABORT( "Quadrature rule not implemented" )
   }
}

// no discretization specified; Implementaion of dot(std::function, Discr) not required!
template class L2Space< 5, Undefined, real_t >;
template class L2Space< 7, Undefined, real_t >;
template class L2Space< 5, Undefined, Point3D >;
template class L2Space< 7, Undefined, Point3D >;

// P1
template <>
void L2Space< 5, P1Function< real_t > >::dot( const std::function< real_t( const Point3D& ) >& f, P1Function< real_t >& b ) const
{
   return dot< P1VariableOperator, forms::p1_linear_form_blending_q5 >( f, b );
}

// P2
template <>
void L2Space< 7, P2Function< real_t > >::dot( const std::function< real_t( const Point3D& ) >& f, P2Function< real_t >& b ) const
{
   return dot< P2ElementwiseOperator, forms::p2_linear_form_blending_q7 >( f, b );
}

} // namespace hyteg
