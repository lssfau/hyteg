#pragma once

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/fenics/ufc_traits.hpp"
#include "tinyhhg_core/p2functionspace/P2Form.hpp"

// P1 to P2
#include "tinyhhg_core/mixedoperators/generated/p1_to_p2_divt.h"

namespace hhg {

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P1ToP2FenicsForm : public Form
{
public:
  void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
  {
     Matrixr<6, 3> localStiffnessMatrix;
     computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
     out[0] = localStiffnessMatrix( 0, 0 );
     out[1] = localStiffnessMatrix( 0, 1 );
     out[2] = localStiffnessMatrix( 0, 2 );
  }

  void integrateVertexToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
  {
     Matrixr<6, 3> localStiffnessMatrix;
     computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
     out[0] = localStiffnessMatrix( 5, 0 );
     out[1] = localStiffnessMatrix( 5, 1 );
     out[2] = localStiffnessMatrix( 5, 2 );
  }

  void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

  real_t integrate( const std::array< Point3D, 4 >& coords, const P2Form::dofPosByVertexPair3D &cntrPos,
                    const P2Form::dofPosByVertexPair3D &leafPos ) const { WALBERLA_ABORT( "Not implemented." ); }

  bool assemble2D() const override { return !std::is_same< UFCOperator2D, hhg::fenics::NoAssemble >::value; }

  bool assemble3D() const override { return !std::is_same< UFCOperator3D, hhg::fenics::NoAssemble >::value; }

  bool assembly2DDefined() const override { return !std::is_same< UFCOperator2D, hhg::fenics::UndefinedAssembly >::value; }

  bool assembly3DDefined() const override { return !std::is_same< UFCOperator3D, hhg::fenics::UndefinedAssembly >::value; }

private:
  void computeLocalStiffnessMatrix( const std::array< Point3D, 3 >& coords, Matrixr<6, 3>& localStiffnessMatrix ) const
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
