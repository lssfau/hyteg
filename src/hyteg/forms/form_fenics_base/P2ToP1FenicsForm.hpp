#pragma once

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/fenics/ufc_traits.hpp"
#include "hyteg/forms/P2Form.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

// P2 to P1
#include "hyteg/forms/form_fenics_generated/p2_to_p1_div.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

namespace hyteg {

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

  // Method invoked by P1Elements3D (we return row of element matrix for vertex 0
  // at the vertex columns), so actually this should be better called
  // integrateVertexToVertex, but this would imply changing the P1 stuff, too.
  void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const {
    Matrixr<4,10> elMat;
    computeLocalStiffnessMatrix( coords, elMat );
    uint_t rowIdx = fenics::P2DoFMap[ 0 ][ 0 ];
    out[0] = elMat( rowIdx, fenics::P2DoFMap[ 0 ][ 0 ] );
    out[1] = elMat( rowIdx, fenics::P2DoFMap[ 1 ][ 1 ] );
    out[2] = elMat( rowIdx, fenics::P2DoFMap[ 2 ][ 2 ] );
    out[3] = elMat( rowIdx, fenics::P2DoFMap[ 3 ][ 3 ] );
  }

  real_t integrate( const std::array< Point3D, 4 >& coords, const P2Form::dofPosByVertexPair3D &cntrPos,
                    const P2Form::dofPosByVertexPair3D &leafPos ) const {

    Matrixr<4,10> elMat;
    computeLocalStiffnessMatrix( coords, elMat );
    WALBERLA_ASSERT_LESS( cntrPos[0], 4 );
    WALBERLA_ASSERT_LESS( cntrPos[1], 4 );
    uint_t rowIdx = fenics::P2DoFMap[ cntrPos[0] ][ cntrPos[1] ];
    uint_t colIdx = fenics::P2DoFMap[ leafPos[0] ][ leafPos[1] ];

    return real_c( elMat( rowIdx, colIdx ) );
  }

  template <long unsigned int size>
  std::array<real_t,size> integrate(  const std::array< Point3D, 4 >& coords,
                                      const P2Form::dofPosByVertexPair3D &cntrPos,
                                      const std::array<P2Form::dofPosByVertexPair3D,size> &leafPos ) const {
    WALBERLA_ABORT( "Missing implementation in P2ToP1FenicsForm" );
  }

  bool assemble2D() const override { return !std::is_same< UFCOperator2D, hyteg::fenics::NoAssemble >::value; }

  bool assemble3D() const override { return !std::is_same< UFCOperator3D, hyteg::fenics::NoAssemble >::value; }

  bool assembly2DDefined() const override { return !std::is_same< UFCOperator2D, hyteg::fenics::UndefinedAssembly >::value; }

  bool assembly3DDefined() const override { return !std::is_same< UFCOperator3D, hyteg::fenics::UndefinedAssembly >::value; }

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

  void computeLocalStiffnessMatrix( const std::array< Point3D, 4 >& coords, Matrixr<4, 10>& localStiffnessMatrix ) const
  {
    real_t fenicsCoords[12];
    for( int node = 0; node < 4; ++node ) {
      for( int dim = 0; dim < 3; ++dim ) {
        fenicsCoords[node*3+dim] = coords[node][dim];
      }
    }
    UFCOperator3D gen;
    gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
   }

};

} // namespace hyteg
