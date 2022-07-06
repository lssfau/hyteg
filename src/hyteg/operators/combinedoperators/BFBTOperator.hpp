

#include <type_traits>

//#include "hyteg/forms/form_hyteg_generated/p2/p2_invk_mass_affine_q6.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.cpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_k_mass_affine_q4.hpp"
#include "hyteg/p2functionspace/P2FunctionApplyOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/solvers/CGSolver.hpp"

namespace hyteg {

using walberla::real_t;
template < class VelocityBlockOpType, class DivOpType, class DivOpTType >
class BFBTOperator : public Operator< typename DivOpType::dstType, typename DivOpType::dstType >
{
   // sub operator consistency
   static_assert( std::is_same< typename VelocityBlockOpType::srcType, typename VelocityBlockOpType::dstType >::value == true );
   static_assert( std::is_same< typename VelocityBlockOpType::srcType, typename DivOpType::srcType >::value == true );
   static_assert( std::is_same< typename DivOpType::dstType, typename DivOpTType::srcType >::value == true );
   using PressureFunctionType = typename DivOpType::dstType;
   using VelocityFunctionType = typename DivOpType::srcType;
   using VectorVelocityMassOp = VectorToVectorOperator< walberla::real_t, P2VectorFunction, P2VectorFunction >;

 public:
   BFBTOperator( const std::shared_ptr< PrimitiveStorage >&       storage,
                 size_t                                           minLevel,
                 size_t                                           maxLevel,
                 std::function< real_t( const hyteg::Point3D& ) > viscosity,
                 const BoundaryCondition& VelocitySpaceBC )
   : Operator< PressureFunctionType, PressureFunctionType >( storage, maxLevel, maxLevel )
   , A( storage, maxLevel, maxLevel, viscosity ) // A should be a discretization of the velocity block with non constant viscosity
   , Z( storage, maxLevel, viscosity, VelocitySpaceBC )
   , ZSolver( std::make_shared< CGSolver< BFBTSubOperator > >( storage, maxLevel, maxLevel ) )
   , storage_( storage )
   {
      A.computeInverseDiagonalOperatorValues();
      Z.Minv.diagonal_ = A.getInverseDiagonalValues();
      //ZSolver->setPrintInfo( true );
      //printComponentMatrices(maxLevel, storage);
   }

   void apply( const PressureFunctionType& src,
               const PressureFunctionType& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      //S_{w-BFBT}^{-1} = Z^{-1}*(BM(w)^{-1}AM(w)^{-1}B^t)*Z^{-1}

      // helper
      PressureFunctionType pTmp( "pTmp", storage_, level, level );
      VelocityFunctionType vTmp( "vTmp", storage_, level, level, Z.VelocitySpaceBC_ );
      VelocityFunctionType vTmp2( "vTmp2", storage_, level, level, Z.VelocitySpaceBC_ );
      
      
    
      pTmp.copyBoundaryConditionFromFunction(src);
      
      vTmp.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All^hyteg::DirichletBoundary );
      vTmp2.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All^hyteg::DirichletBoundary );
  /*    vTmp.setBoundaryCondition(BoundaryCondition(hyteg::All^hyteg::DirichletBoundary));
      vTmp2.setBoundaryCondition(BoundaryCondition(hyteg::All^hyteg::DirichletBoundary));

  
      WALBERLA_LOG_INFO_ON_ROOT( "pTmp bc type for 1: " << pTmp.getBoundaryCondition().getBoundaryType( 1 ) );
      WALBERLA_LOG_INFO_ON_ROOT( "pTmp bc type for 0: " << pTmp.getBoundaryCondition().getBoundaryType( 0 ) );
      WALBERLA_LOG_INFO_ON_ROOT( "vTmp bc type for 1: " << vTmp.getBoundaryCondition().getBoundaryType( 1 ) );
      WALBERLA_LOG_INFO_ON_ROOT( "vTmp bc type for 0: " << vTmp.getBoundaryCondition().getBoundaryType( 0 ) );
      */
      /*
      pTmp.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All );
      vTmp.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All );
      vTmp2.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All );
      */
   
      ZSolver->solve( Z, pTmp, src, level ); // pTmp = Z^{-1}*src

      Z.BT.apply( pTmp, vTmp, level, hyteg::Inner | hyteg::NeumannBoundary ); // vTmp = B^t*pTmp

      Z.Minv.apply( vTmp, vTmp2, level, hyteg::Inner | hyteg::NeumannBoundary ); // vTmp2 = M^{-1}*vTmp

      A.apply( vTmp2, vTmp, level, hyteg::Inner | hyteg::NeumannBoundary ); // vtmp= A*vTmp2

      Z.Minv.apply( vTmp, vTmp2, level, hyteg::Inner | hyteg::NeumannBoundary ); // vTmp2 = M^{-1}*vTmp

      Z.B.apply( vTmp2, pTmp, level, hyteg::Inner | hyteg::NeumannBoundary ); // pTmp = B*vTmp2

      ZSolver->solve( Z, dst, pTmp, level ); // dst = Z^{-1}*pTmp
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                  const typename PressureFunctionType::template FunctionType< idx_t >& src,
                  const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                  size_t                                                               level,
                  DoFType                                                              flag ) const
   {
      WALBERLA_ABORT( "Requires exact inversion of A, should not be used!" );
   }

   void printComponentMatrices( size_t level, std::shared_ptr< PrimitiveStorage > storage )
   {
      typename VelocityFunctionType::template FunctionType< idx_t > VelocityNum(
          "VelocityNum", storage /*Operator< PressureFunctionType, PressureFunctionType >::storage_*/, level, level );

      VelocityNum.enumerate( level );

      PETScSparseMatrix< P2FunctionApplyOperator > MinvMat;
      MinvMat.createMatrixFromOperator( Z.Minv, level, VelocityNum );

      PETScSparseMatrix< P2ElementwiseAffineEpsilonOperator > AMat;
      AMat.createMatrixFromOperator( A, level, VelocityNum );

      WALBERLA_LOG_INFO_ON_ROOT( "Minv:" );
      //MatView( MinvMat.get(), PETSC_VIEWER_STDOUT_WORLD );
      MinvMat.print( "MinvMat", false, PETSC_VIEWER_DEFAULT );
      WALBERLA_LOG_INFO_ON_ROOT( "A:" );
      //MatView( AMat.get(), PETSC_VIEWER_STDOUT_WORLD );
      AMat.print( "AMat", false, PETSC_VIEWER_DEFAULT );
   }

   // this implements the suboperator Z = BM(w)^{-1}B^T
   class BFBTSubOperator : public Operator< PressureFunctionType, PressureFunctionType >
   {
    public:
      BFBTSubOperator( const std::shared_ptr< PrimitiveStorage >&       storage,
                       size_t                                           level,
                       std::function< real_t( const hyteg::Point3D& ) > viscosity,
                       const BoundaryCondition& VelocitySpaceBC 
                       //std::shared_ptr<P2VectorFunction< real_t >> diagonal
                       )
      : Operator< PressureFunctionType, PressureFunctionType >( storage, level, level )
      , B( storage, level, level )
      , BT( storage, level, level )
      , Minv( storage, level )
      , storage_( storage )
      , VelocitySpaceBC_(VelocitySpaceBC)
      {}

      void apply( const PressureFunctionType& src,
                  const PressureFunctionType& dst,
                  size_t                      level,
                  DoFType                     flag,
                  UpdateType                  updateType = Replace ) const
      {
         VelocityFunctionType tmp0( "tmp0", storage_, level, level,VelocitySpaceBC_ );
         VelocityFunctionType tmp1( "tmp1", storage_, level, level,VelocitySpaceBC_ );

         tmp0.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All^hyteg::DirichletBoundary );
         tmp1.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All^hyteg::DirichletBoundary );

         
     // tmp0.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All );
    //  tmp1.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::All );
         //tmp0.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::DirichletBoundary );
         // tmp1.interpolate( []( const hyteg::Point3D& ) { return real_c( 0 ); }, level, hyteg::DirichletBoundary );

         //tmp0.setBoundaryCondition(BoundaryCondition(hyteg::Inner));
         //tmp1.setBoundaryCondition(BoundaryCondition(hyteg::Inner));

         //WALBERLA_LOG_INFO_ON_ROOT( "tmp0 bc type for 1: " << tmp0.getBoundaryCondition().getBoundaryType( 1 ) );

         //WALBERLA_LOG_INFO_ON_ROOT( "tmp0 bc type for 0: " << tmp0.getBoundaryCondition().getBoundaryType( 0 ) );
         BT.apply( src, tmp0, level, hyteg::Inner | hyteg::NeumannBoundary ); // tmp0 = B^T*src

         Minv.apply( tmp0, tmp1, level, hyteg::Inner | hyteg::NeumannBoundary ); // tmp1 = M^{-1}*tmp0

         B.apply( tmp1, dst, level, hyteg::Inner | hyteg::NeumannBoundary ); // dst = B*tmp1
      }

      void toMatrix( const std::shared_ptr< SparseMatrixProxy >&                          mat,
                     const typename PressureFunctionType::template FunctionType< idx_t >& src,
                     const typename PressureFunctionType::template FunctionType< idx_t >& dst,
                     size_t                                                               level,
                     DoFType                                                              flag ) const
      {
         // Z = BM(w)^{-1}B^T to matrix

         //WALBERLA_ABORT( "Should not be used yet!" );

         // Sparse matrix objects for sub operators
         PETScSparseMatrix< P2FunctionApplyOperator > MinvMat;
         PETScSparseMatrix< DivOpType >               BMat;
         PETScSparseMatrix< DivOpTType >              BTMat;

         // numerators
         typename VelocityFunctionType::template FunctionType< idx_t > VelocityNum(
             "VelocityNum", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
         typename PressureFunctionType::template FunctionType< idx_t > PressureNum(
             "PressureNum", Operator< PressureFunctionType, PressureFunctionType >::storage_, level, level );
         VelocityNum.enumerate( level );
         PressureNum.enumerate( level );

         // assembly
         MinvMat.createMatrixFromOperator( Minv, level, VelocityNum );
         BMat.createMatrixFromOperator( B, level, VelocityNum, PressureNum );
         BTMat.createMatrixFromOperator( BT, level, PressureNum, VelocityNum );

         //MatView(MinvMat.get(), PETSC_VIEWER_STDOUT_WORLD);
         //MatView(BMat.get(), PETSC_VIEWER_STDOUT_WORLD);
         //MatView(BTMat.get(), PETSC_VIEWER_STDOUT_WORLD);

         //TODO assume symmetry false
         //MMat.applyDirichletBC(VelocityNum, level);
         //BMat.applyDirichletBC(PressureNum, level);
         //BTMat.applyDirichletBC(VelocityNum, level);

         // Proxy objects for createFromMatrixProduct
         auto MinvMatProxy = std::make_shared< PETScSparseMatrixProxy >( MinvMat.get() );
         auto BMatProxy    = std::make_shared< PETScSparseMatrixProxy >( BMat.get() );
         auto BTMatProxy   = std::make_shared< PETScSparseMatrixProxy >( BTMat.get() );

         // final BMinvBT matrix creation
         mat->createFromMatrixProduct( { BTMatProxy, MinvMatProxy, BMatProxy } );
      }

      DivOpType                                 B;
      DivOpTType                                BT;
      P2FunctionApplyOperator                   Minv;
      const std::shared_ptr< PrimitiveStorage > storage_;
      const BoundaryCondition VelocitySpaceBC_;
   };

 private:
   VelocityBlockOpType                            A;
   BFBTSubOperator                                Z;
   std::shared_ptr< CGSolver< BFBTSubOperator > > ZSolver;
   const std::shared_ptr< PrimitiveStorage >      storage_;
};

using BFBT_P2P1 = BFBTOperator< P2ElementwiseAffineEpsilonOperator, P2ToP1ElementwiseDivOperator, P1ToP2ElementwiseDivTOperator >;

} // namespace hyteg
