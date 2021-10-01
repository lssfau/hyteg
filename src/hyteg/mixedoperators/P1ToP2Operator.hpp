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
#pragma once

#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P1ToP2FenicsForm.hpp"
#include "hyteg/forms/P1WrapperForm.hpp"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t;

template< class P1ToP2Form >
class P1ToP2ConstantOperator : public Operator<P1Function < real_t>, P2Function<real_t> > {
 public:

  P1ToP2ConstantOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel), vertexToVertex(storage, minLevel, maxLevel),
        vertexToEdge(storage, minLevel, maxLevel)
  {
  }


  P1ConstantOperator< P1WrapperForm<P1ToP2Form> > const & getVertexToVertexOpr() const {
    return vertexToVertex;
  }

  VertexDoFToEdgeDoFOperator< P1ToP2Form > const & getVertexToEdgeOpr() const {
    return vertexToEdge;
  }

   void apply( const P1Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      vertexToVertex.apply( src, dst.getVertexDoFFunction(), level, flag, updateType );
      vertexToEdge.apply( src, dst.getEdgeDoFFunction(), level, flag, updateType );
  }

  void smooth_gs(P1Function< real_t > & dst, P2Function< real_t > & rhs, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("not implemented");
  }

private:

  P1ConstantOperator< P1WrapperForm<P1ToP2Form> > vertexToVertex;
  VertexDoFToEdgeDoFOperator< P1ToP2Form > vertexToEdge;

};

typedef P1ToP2ConstantOperator< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > > P1ToP2ConstantDivTxOperator;
typedef P1ToP2ConstantOperator< P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > > P1ToP2ConstantDivTyOperator;
typedef P1ToP2ConstantOperator< P1ToP2FenicsForm< fenics::NoAssemble,                      p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > > P1ToP2ConstantDivTzOperator;

}
