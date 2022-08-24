
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

#include "hyteg/forms/form_hyteg_generated/p2/p2_invk_mass_affine_q6.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_k_mass_affine_q4.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/p2functionspace/P2VectorApplyOperator.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.cpp"

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
   , Z( storage, maxLevel, viscosity)
   , ZSolver( std::make_shared< CGSolver< wBFBTSubOperator > >( storage, maxLevel, maxLevel ) )
   {
      A.computeInverseDiagonalOperatorValues();
      Z.Minv.diagonal_ = A.getInverseDiagonalValues();
      ZSolver->setPrintInfo( true );
 /*   ZSolver->assumeSymmetry( false );
    std::function< real_t( const hyteg::Point3D& ) > ones = []( const hyteg::Point3D& ) { return real_c( 1 ); };
    hyteg::P1Function< real_t > ones_( "ones", storage, maxLevel, maxLevel ); 
    ones_.interpolate(ones,maxLevel);
    ZSolver->setNullSpace(ones_, maxLevel);*/
 //   ZSolver->setConstantNullSpace();

   }

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
      pTmp.setBoundaryCondition(BoundaryCondition(DoFType::Inner));
      vTmp.setBoundaryCondition(BoundaryCondition(DoFType::Inner));
      vTmp2.setBoundaryCondition(BoundaryCondition(DoFType::Inner));
      std::function< real_t( const hyteg::Point3D& ) > zero = []( const hyteg::Point3D& ) { return real_c( 0 ); };

      
      VTKOutput vtkOutput_pTmp(
          "../../output", "BFBT_intres_BFBTApply_pTmp", Operator< PressureFunctionType, PressureFunctionType >::storage_ );
      VTKOutput vtkOutput_vTmp(
          "../../output", "BFBT_intres_BFBTApply_vTmp", Operator< PressureFunctionType, PressureFunctionType >::storage_ );
      VTKOutput vtkOutput_vTmp2(
          "../../output", "BFBT_intres_BFBTApply_vTmp2", Operator< PressureFunctionType, PressureFunctionType >::storage_ );
      vtkOutput_pTmp.add( pTmp );
      vtkOutput_pTmp.add( dst );

      vtkOutput_vTmp.add( vTmp );
      vtkOutput_vTmp2.add( vTmp2 );
      vtkOutput_pTmp.write( level, 0 );
      

       // WALBERLA_LOG_INFO_ON_ROOT("ZSolver pre");
    
      ZSolver->solve( Z, pTmp, src, level ); // pTmp = Z^{-1}*src
      vtkOutput_pTmp.write( level, 1 );

      Z.BT.apply( pTmp, vTmp, level, hyteg::Inner ); // vTmp = B^t*pTmp
      //vtkOutput_vTmp.write( level, 0 );

      //vTmp2.interpolate(zero, level);
      
      Z.Minv.apply(vTmp, vTmp2,  level,  hyteg::Inner);
      //Z.MSolver->solve(Z.M, vTmp2, vTmp, level ); // vTmp2= M(w)^{-1}*vTmp
      vtkOutput_vTmp2.write( level, 0 );

      A.apply( vTmp2, vTmp, level,  hyteg::Inner ); // vtmp= A*vTmp2
      vtkOutput_vTmp.write( level, 1 );

      vTmp2.interpolate(zero, level);
      Z.Minv.apply(vTmp, vTmp2,  level,  hyteg::Inner);
      //Z.MSolver->solve(Z.M, vTmp2, vTmp, level ); // vTmp2= M(w)^{-1}*vTmp
      vtkOutput_vTmp2.write( level, 1 );

      Z.B.apply( vTmp2, pTmp, level,  hyteg::Inner ); // pTmp = B*vTmp2
      vtkOutput_pTmp.write( level, 2 );

    //  WALBERLA_LOG_INFO_ON_ROOT("ZSolver post");
      ZSolver->solve( Z, dst, pTmp, level ); // dst = Z^{-1}*pTmp
      vtkOutput_pTmp.write( level, 3 );
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
                        std::function< real_t( const hyteg::Point3D& ) > viscosity
                        //std::shared_ptr<P2VectorFunction< real_t >> diagonal 
                        )
      : Operator< PressureFunctionType, PressureFunctionType >( storage, level, level )
      , B( storage, level, level )
      , BT( storage, level, level )
      , Minv( storage, level )
      
      //, M( storage, level, level) 
      //, MSolver( std::make_shared< PETScLUSolver< VectorVelocityMassOp > >( storage, level ) )
      {

          //MSolver->setPrintInfo( true );
         /* Becomes singular
         auto P2LumpedViscMassMatrix = std::make_shared< P2BlendingLumpedInverseDiagonalOperator >(
             storage,
             level,
             level,
             std::make_shared< P2RowSumForm >( std::make_shared< forms::p2_sqrtk_mass_affine_q6 >( viscosity, viscosity ) ) );
             */
        
        /*
         auto form = hyteg::forms::p2_sqrtk_mass_affine_q4( viscosity, viscosity );
         auto P2Mass = std::make_shared<P2ElementwiseOperator< hyteg::forms::p2_sqrtk_mass_affine_q4>>(storage, level, level, form);
         M.setSubOperator( 0, 0, P2Mass );
         M.setSubOperator( 1, 1, P2Mass );
         */
         //MSolver.flag_ =  hyteg::Inner;
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
         
        tmp0.setBoundaryCondition(BoundaryCondition(DoFType::Inner));
        tmp1.setBoundaryCondition(BoundaryCondition(DoFType::Inner));
        
         VTKOutput vtkOutput_pTmp(
             "../../output", "BFBT_intres_BFBTSubOpApply", Operator< PressureFunctionType, PressureFunctionType >::storage_ );
         vtkOutput_pTmp.add( tmp0 );
         vtkOutput_pTmp.add( tmp1 );
         vtkOutput_pTmp.add( dst );
         vtkOutput_pTmp.add( src ); 

         vtkOutput_pTmp.write( level, 0 );

         BT.apply( src, tmp0, level,  hyteg::Inner );
         vtkOutput_pTmp.write( level, 1 );

        //TODO: apply M solver here
        //WALBERLA_LOG_INFO_ON_ROOT("Inner MSolver");
         //MSolver->solve(M, tmp1, tmp0, level);
         Minv.apply( tmp0, tmp1, level, flag );
         vtkOutput_pTmp.write( level, 2 );

         B.apply( tmp1, dst, level,  hyteg::Inner );
         vtkOutput_pTmp.write( level, 3 );
      }

      void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                     const typename PressureFunctionType::template FunctionType< idx_t >& src,
                     const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                     size_t                                                               level,
                     DoFType                                                              flag ) const
      {
         // Z = BM(w)^{-1}B^T to matrix

        WALBERLA_ABORT( "Should not be used yet!" );

        /*
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

         
        MatView(MMat.get(), PETSC_VIEWER_STDOUT_WORLD);
        MatView(BMat.get(), PETSC_VIEWER_STDOUT_WORLD);
        MatView(BTMat.get(), PETSC_VIEWER_STDOUT_WORLD);
         //TODO assume symmetry false
         //MMat.applyDirichletBC(VelocityNum, level);
         //BMat.applyDirichletBC(PressureNum, level);
         //BTMat.applyDirichletBC(VelocityNum, level);

         // Proxy objects for createFromMatrixProduct
         auto MMatProxy  = std::make_shared< PETScSparseMatrixProxy >( MMat.get() );
         auto BMatProxy  = std::make_shared< PETScSparseMatrixProxy >( BMat.get() );
         auto BTMatProxy = std::make_shared< PETScSparseMatrixProxy >( BTMat.get() );

         // final BMBT matrix creation
         mat->createFromMatrixProduct( { BTMatProxy, MMatProxy, BMatProxy } );*/
      }

      DivOpType  B;
      DivOpTType BT;
      //P2BlendingLumpedInverseDiagonalOperator M;
      P2VectorApplyOperator Minv;
      //VectorVelocityMassOp M;
      //std::shared_ptr< PETScLUSolver< VectorVelocityMassOp > >      MSolver;
   };

 private:
   VelocityBlockOpType                                  A;
   wBFBTSubOperator                                     Z;
   std::shared_ptr< CGSolver< wBFBTSubOperator > >      ZSolver;
};

using WBFBT_P2P1 =
    wBFBTOperator< P2ElementwiseAffineEpsilonOperator, P2ToP1ElementwiseDivOperator, P1ToP2ElementwiseDivTOperator >;

} // namespace hyteg