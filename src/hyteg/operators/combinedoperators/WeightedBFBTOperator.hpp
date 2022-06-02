
/* Implements the w-BFBT approximation for the Schur complement from [1].
   This operator is the result of combined application of:
    S_{w-BFBT}^{-1} = Z^{-1}*(BM(w)^{-1}AM(w)^{-1}B^t)*Z^{-1}   (eq. 1.8 in [1])
     
    with: 
        -  A, the (1,1) block in the discretization of the Stokes operator
        -  B, the divergence operator
        -  M(w), the lumped, viscosity weighted velocity mass matrix (see eq. 1.4 and 1.8 in [1], M(w) is short for \tilde(M)_u(w_l))
        -  Z = BM(w)^{-1}B^T, the BFBT suboperator, introduced for convenience 
        (and assuming the left viscosity weighting function w_l is equal to the right viscosity weighting function w_r)

    The operator is intended to be used in a block diagonal preconditioner for the Stokes equation
    with highly heterogeneous viscosity.

[1]: "WEIGHTED BFBT PRECONDITIONER FOR STOKES FLOW PROBLEMS WITH HIGHLY HETEROGENEOUS VISCOSITY"
      by JOHANN RUDI, GEORG STADLER, and OMAR GHATTAS
*/

#include <type_traits>
namespace hyteg {

using walberla::real_t;
template < class VelocityBlockOpType, class DivOpType, class DivOpTType >
class wBFBTOperator : public Operator< typename DivOpType::dstType, typename DivOpType::dstType >
{
   // sub operator consistency
   static_assert( std::is_same< typename VelocityBlockOpType::srcType, typename VelocityBlockOpType::dstType >::value == true );
   static_assert( std::is_same< typename VelocityBlockOpType::srcType, typename DivOpType::srcType >::value == true );
   static_assert( std::is_same< typename DivOpType::dstType, typename DivOpTType::srcType >::value == true );
   using PressureFunctionType = typename DivOpType::dstType;
   using VelocityFunctionType = typename DivOpType::srcType;
   using VectorVelocityMassOp = VectorToVectorOperator< walberla::real_t, P2VectorFunction, P2VectorFunction >;

 public:
   wBFBTOperator( const std::shared_ptr< PrimitiveStorage >&       storage,
                  size_t                                           minLevel,
                  size_t                                           maxLevel,
                  std::function< real_t( const hyteg::Point3D& ) > viscosity )
   : Operator< PressureFunctionType, PressureFunctionType >( storage, maxLevel, maxLevel )
   , A( storage, maxLevel, maxLevel, viscosity ) // A should be a discretization of the velocity block with non constant viscosity
   , Z( storage, maxLevel, viscosity )
   , ZSolver( std::make_shared< PETScLUSolver< wBFBTSubOperator > >( storage, maxLevel ) )
   {}

   void apply( const PressureFunctionType& src,
               const PressureFunctionType& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      //S_{w-BFBT}^{-1} = Z^{-1}*(BM(w)^{-1}AM(w)^{-1}B^t)*Z^{-1}
      PressureFunctionType pTmp( "pTmp", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
      VelocityFunctionType vTmp( "vTmp", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
      VelocityFunctionType vTmp2( "vTmp2", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
      ZSolver->solve( Z, pTmp, src, level ); // pTmp = Z^{-1}*src
      Z.BT.apply( pTmp, vTmp, level, flag ); // vTmp = B^t*pTmp
      Z.M.apply( vTmp, vTmp2, level, flag ); // vTmp2= M(w)^{-1}*vTmp
      A.apply( vTmp2, vTmp, level, flag );   // vtmp= A*vTmp2
      Z.M.apply( vTmp, vTmp2, level, flag ); // vTmp2= M(w)^{-1}*vTmp
      Z.B.apply( vTmp2, pTmp, level, flag ); // pTmp = B*vTmp2
      ZSolver->solve( Z, pTmp, dst, level ); // dst = Z^{-1}*pTmp
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                  const typename PressureFunctionType::template FunctionType< idx_t >& src,
                  const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                  size_t                                                               level,
                  DoFType                                                              flag ) const
   {
      WALBERLA_ABORT( "Requires exact inversion of A, should not be used!" );
   }

   // this implements the suboperator Z = BM(w)^{-1}B^T
   class wBFBTSubOperator : public Operator< PressureFunctionType, PressureFunctionType >
   {
    public:
      wBFBTSubOperator( const std::shared_ptr< PrimitiveStorage >&       storage,
                        size_t                                           level,
                        std::function< real_t( const hyteg::Point3D& ) > viscosity )
      : Operator< PressureFunctionType, PressureFunctionType >( storage, level, level )
      , B( storage, level, level )
      , BT( storage, level, level )
      , M( storage, level, level )
      {
         auto P2LumpedViscMassMatrix = std::make_shared< P2BlendingLumpedInverseDiagonalOperator >(
             storage,
             level,
             level,
             std::make_shared< P2RowSumForm >( std::make_shared< forms::p2_sqrtk_mass_affine_q4 >( viscosity, viscosity ) ) );
         M.setSubOperator( 0, 0, P2LumpedViscMassMatrix );
         M.setSubOperator( 1, 1, P2LumpedViscMassMatrix );
         // M.setSubOperator(2,2,P2LumpedViscMassMatrix);
      }

      void apply( const PressureFunctionType& src,
                  const PressureFunctionType& dst,
                  size_t                      level,
                  DoFType                     flag,
                  UpdateType                  updateType = Replace ) const
      {
         VelocityFunctionType tmp0( "tmp0", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
         VelocityFunctionType tmp1( "tmp1", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
         BT.apply( src, tmp0, level, flag );
         M.apply( tmp0, tmp1, level, flag );
         B.apply( tmp1, dst, level, flag );
      }

      void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                     const typename PressureFunctionType::template FunctionType< idx_t >& src,
                     const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                     size_t                                                               level,
                     DoFType                                                              flag ) const
      {
         // Z = BM(w)^{-1}B^T to matrix

         // Sparse matrix objects for sub operators
         PETScSparseMatrix< VectorVelocityMassOp > MMat;
         PETScSparseMatrix< DivOpType >            BMat;
         PETScSparseMatrix< DivOpTType >           BTMat;

         // numerators
         typename VelocityFunctionType::template FunctionType< idx_t > VelocityNum(
             "VelocityNum", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
         typename PressureFunctionType::template FunctionType< idx_t > PressureNum(
             "PressureNum", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
         VelocityNum.enumerate( level );
         PressureNum.enumerate( level );

         // assembly
         MMat.createMatrixFromOperator( M, level, VelocityNum );
         BMat.createMatrixFromOperator( B, level, VelocityNum, PressureNum );
         BTMat.createMatrixFromOperator( BT, level, PressureNum, VelocityNum );

         // Proxy objects for createFromMatrixProduct
         auto MMatProxy  = std::make_shared< PETScSparseMatrixProxy >( MMat.get() );
         auto BMatProxy  = std::make_shared< PETScSparseMatrixProxy >( BMat.get() );
         auto BTMatProxy = std::make_shared< PETScSparseMatrixProxy >( BTMat.get() );

         // final BMBT matrix creation
         mat->createFromMatrixProduct( { BTMatProxy, MMatProxy, BMatProxy } );
      }

      DivOpType  B;
      DivOpTType BT;
      //P2BlendingLumpedInverseDiagonalOperator M;
      VectorVelocityMassOp M;
   };

 private:
   VelocityBlockOpType                                  A;
   wBFBTSubOperator                                     Z;
   std::shared_ptr< PETScLUSolver< wBFBTSubOperator > > ZSolver;
};

using WBFBT_P2P1 =
    wBFBTOperator< P2ElementwiseAffineEpsilonOperator, P2ToP1ElementwiseDivOperator, P1ToP2ElementwiseDivTOperator >;

} // namespace hyteg