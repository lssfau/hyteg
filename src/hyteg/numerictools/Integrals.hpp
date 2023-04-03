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

// todo make L2 a class and rename file l2space
/// \file
/// This source code file contains free functions for evaluating integrals

#pragma once

#include "core/Abort.h"
#include "core/math/KahanSummation.h"

#include "hyteg/forms/form_hyteg_generated/p0/p0_linear_form_blending_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_linear_form_blending_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_linear_form_blending_q6.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"
#include "hyteg/volumedofspace/FaceDoFIndexing.hpp"

namespace hyteg {
namespace integrals {

using walberla::real_t;
using walberla::uint_t;

// =================================================================================================

/// \brief Integrate a form over a given macro element
/// \tparam PrimitiveType  type of T, i.e, either cell or face
/// \tparam P0Form         type of linear form
///
/// \param T            primitive over which f shall be integrated
/// \param form         P0form describing the local integral
/// \param lvl          grid level to work on
///
/// \return value of integral over T
template < class PrimitiveType, class P0Form >
real_t integrate( const PrimitiveType& T, P0Form& form, const uint_t& lvl )
{
   form.setGeometryMap( T.getGeometryMap() );

   walberla::math::KahanAccumulator< real_t > integral; // value of integral over T;

   Matrixr< 1, 1 > integral_micro; // value of integral over micro element

   if constexpr ( std::is_same< PrimitiveType, Cell >::value )
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
            form.integrateAll( coords, integral_micro );

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
            form.integrateAll( coords, integral_micro );

            integral += integral_micro( 0, 0 );
         }
      }
   }

   return integral.get();
}

/// \brief Compute an approximate L^2 dot product (u,v)_L^2(Ω) = ∫_Ω uv using a quadrature rule of order Q
/// \tparam Q           order of quadrature rule to compute the integrals on each micro element
///
/// \param storage      storage corresponding to the domain Ω
/// \param u            function in L^2(Ω)
/// \param v            function in L^2(Ω)
/// \param lvl          grid level to work on
///
/// \return (u,v)_L^2(Ω)
template < int Q = 6 >
real_t L2dot( const std::shared_ptr< PrimitiveStorage >&                           storage,
              const std::function< real_t( const Point3D&, const PrimitiveID& ) >& u,
              const std::function< real_t( const Point3D&, const PrimitiveID& ) >& v,
              const uint_t&                                                        lvl )
{
   // setup linear form with correct quadrature rule
   if constexpr ( Q != 6 )
   {
      WALBERLA_ABORT( "L2norm only implemented with Q=6" );
   }
   std::shared_ptr< forms::p0_linear_form_blending_q6 > form;

   real_t localsum = 0.0;

   if ( storage->hasGlobalCells() )
   {
      // 3D: integrate over all cells
      std::vector< PrimitiveID > cellIDs = storage->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : localsum )
#endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         auto id   = cellIDs[uint_c( i )];
         auto cell = storage->getCell( id );

         auto uv = [&]( const Point3D& x ) { return u( x, id ) * v( x, id ); };
         form    = std::make_shared< forms::p0_linear_form_blending_q6 >( uv, uv );

         localsum += integrate<>( *cell, *form, lvl );
      }
   }
   else
   {
      // 2D: integrate over all faces
      std::vector< PrimitiveID > faceIDs = storage->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : localsum )
#endif
      for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
      {
         auto id   = faceIDs[uint_c( i )];
         auto face = storage->getFace( id );

         auto uv = [&]( const Point3D& x ) { return u( x, id ) * v( x, id ); };
         form    = std::make_shared< forms::p0_linear_form_blending_q6 >( uv, uv );

         localsum += integrate<>( *face, *form, lvl );
      }
   }

   // compute global sum
   walberla::mpi::allReduceInplace( localsum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   return localsum;
}

/// \brief Compute the L^2 norm ||av + bf||_L2 using a quadrature rule of order Q
///        where a and b are scalars, v is a FE function and f is an L^2 function.
/// \tparam FE_t        type of FE function v
/// \tparam Q           order of quadrature rule to compute the integrals
///
/// \param storage      storage corresponding to the domain Ω
/// \param a            scalar factor for v
/// \param b            scalar factor for f
/// \param v            function in the chosen FE space
/// \param f            function in L^2(Ω)
/// \param lvl          grid level to work on
///
/// \return ||f||_L2
template < class FE_t, int Q = 6 >
real_t L2norm( const std::shared_ptr< PrimitiveStorage >&       storage,
               const real_t&                                    a,
               const real_t&                                    b,
               const FE_t&                                      v,
               const std::function< real_t( const Point3D& ) >& f,
               const uint_t&                                    lvl )
{
   // setup linear form
   if constexpr ( Q != 6 )
   {
      WALBERLA_ABORT( "L2norm only implemented with Q=6" );
   }
   std::shared_ptr< forms::p0_linear_form_blending_q6 > form;
   // integrand (av + bf)^2
   std::function< real_t( const Point3D& ) > av_p_bf_2 = [&]( const Point3D& x ) {
      real_t vx;
      v.evaluate( x, lvl, vx );
      auto av_p_bf = a * vx + b * f( x );
      return av_p_bf * av_p_bf;
   };
   form = std::make_shared< forms::p0_linear_form_blending_q6 >( av_p_bf_2, av_p_bf_2 );

   real_t localsum = 0.0;

   if ( storage->hasGlobalCells() )
   {
      // 3D: integrate over all cells
      std::vector< PrimitiveID > cellIDs = storage->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : localsum )
#endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         auto cell = storage->getCell( cellIDs[uint_c( i )] );

         // todo use P1Function::evaluate for given primitiveID
         // todo implement P2Function::evaluate for given primitiveID

         localsum += integrate<>( *cell, *form, lvl );
      }
   }
   else
   {
      // 2D: integrate over all faces
      std::vector< PrimitiveID > faceIDs = storage->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for reduction( + : localsum )
#endif
      for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
      {
         auto face = storage->getFace( faceIDs[uint_c( i )] );

         // todo use P1Function::evaluate for given primitiveID
         // todo implement P2Function::evaluate for given primitiveID

         localsum += integrate<>( *face, *form, lvl );
      }
   }

   // compute global sum
   walberla::mpi::allReduceInplace( localsum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   return sqrt( localsum );
}

} // namespace integrals
} // namespace hyteg
