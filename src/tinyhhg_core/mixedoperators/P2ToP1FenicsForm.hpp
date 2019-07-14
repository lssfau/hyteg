#pragma once

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/fenics/ufc_traits.hpp"
#include "tinyhhg_core/p2functionspace/P2Form.hpp"

// P2 to P1
#include "tinyhhg_core/mixedoperators/generated/p2_to_p1_div.h"

namespace hhg {

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P2ToP1FenicsForm : public Form
{
public:
  void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
  {
     Matrixr<3, 6> localStiffnessMatrix;
     computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
     out[0] = localStiffnessMatrix( 0, 0 );
     out[1] = localStiffnessMatrix( 0, 1 );
     out[2] = localStiffnessMatrix( 0, 2 );
  }

  void integrateEdgeToVertex( const std::array< Point3D, 3 >& coords, Point3D& out ) const
  {
     Matrixr<3, 6> localStiffnessMatrix;
     computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
     out[0] = localStiffnessMatrix( 0, 3 );
     out[1] = localStiffnessMatrix( 0, 4 );
     out[2] = localStiffnessMatrix( 0, 5 );
  }

  void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

  real_t integrate( const std::array< Point3D, 4 >& coords, const P2Form::dofPosByVertexPair3D &cntrPos,
                    const P2Form::dofPosByVertexPair3D &leafPos ) const { WALBERLA_ABORT( "Not implemented." ); }

  bool assemble2D() const override { return !std::is_same< UFCOperator2D, hhg::fenics::NoAssemble >::value; }

  bool assemble3D() const override { return !std::is_same< UFCOperator3D, hhg::fenics::NoAssemble >::value; }

  bool assembly2DDefined() const override { return !std::is_same< UFCOperator2D, hhg::fenics::UndefinedAssembly >::value; }

  bool assembly3DDefined() const override { return !std::is_same< UFCOperator3D, hhg::fenics::UndefinedAssembly >::value; }

private:
  void computeLocalStiffnessMatrix( const std::array< Point3D, 3 >& coords, Matrixr<3, 6>& localStiffnessMatrix ) const
  {
     real_t fenicsCoords[6];
     fenicsCoords[0] = coords[0][0];
     fenicsCoords[1] = coords[0][1];
     fenicsCoords[2] = coords[1][0];
     fenicsCoords[3] = coords[1][1];
     fenicsCoords[4] = coords[2][0];
     fenicsCoords[5] = coords[2][1];
     UFCOperator2D gen;
     gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
  }
};

} // namespace hhg
