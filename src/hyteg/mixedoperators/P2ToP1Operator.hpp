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
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
#include "hyteg/forms/P1WrapperForm.hpp"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t;

template< class P2ToP1Form >
class P2ToP1ConstantOperator : public Operator<P2Function < real_t>, P1Function<real_t> > {
 public:

  P2ToP1ConstantOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel),
        vertexToVertex(storage, minLevel, maxLevel),
        edgeToVertex(storage, minLevel, maxLevel)
  {
  }

  P1ConstantOperator< P1WrapperForm<P2ToP1Form> > const & getVertexToVertexOpr() const {
    return vertexToVertex;
  }

  EdgeDoFToVertexDoFOperator< P2ToP1Form > const & getEdgeToVertexOpr() const {
    return edgeToVertex;
  }

   void apply( const P2Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      vertexToVertex.apply( src.getVertexDoFFunction(), dst, level, flag, updateType );
      edgeToVertex.apply( src.getEdgeDoFFunction(), dst, level, flag, Add );
   }

  void smooth_gs(P2Function< real_t > & dst, P1Function< real_t > & rhs, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("not implemented");
  }

private:

  P1ConstantOperator< P1WrapperForm<P2ToP1Form> > vertexToVertex;
  EdgeDoFToVertexDoFOperator< P2ToP1Form > edgeToVertex;

};

typedef P2ToP1ConstantOperator< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > > P2ToP1ConstantDivxOperator;
typedef P2ToP1ConstantOperator< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > > P2ToP1ConstantDivyOperator;
typedef P2ToP1ConstantOperator< P2ToP1FenicsForm< fenics::NoAssemble,                     p2_to_p1_tet_div_tet_cell_integral_2_otherwise > > P2ToP1ConstantDivzOperator;

}
