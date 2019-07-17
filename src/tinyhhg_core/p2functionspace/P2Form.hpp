#pragma once

#include "tinyhhg_core/form/Form.hpp"

namespace hhg {

class P2Form : public Form
{
 public:
   virtual ~P2Form() {}

   // 2D P2 VertexDoF
   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   // 3D P2 VertexDoF
   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }
  
  /// Describe a degree of freedom in a P2 element on a tetrahedron
  /// by two vertex indices. There are two cases:
  ///
  /// (a) Both vertex indices are identical, then this is the
  ///     index of the dof associated with this vertex.
  ///
  /// (b) The two vertex indices are different, then this is the
  ///     index of the dof associated with the midpoint of the
  ///     tet's edge given by those two vertices.
  typedef std::array<uint_t,2> dofPosByVertexPair3D;

};

} // namespace hhg
