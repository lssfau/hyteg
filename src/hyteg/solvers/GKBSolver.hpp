/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "core/Abort.h"
#include "core/timing/TimingTree.h"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/operators/AugmentedLargangianOperator.hpp"
#include "hyteg/operators/BlockOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t; 
using walberla::uint_t;

/// Golub-Kahan-Bidiagonalization Solver implementation from "GENERALIZED GOLUB–KAHAN BIDIAGONALIZATION AND STOPPING CRITERIA" by Mario Arioli 

template<
   class SaddlePointOp,// this op is not needed in GKB but added in order to fit the solver class hierarchy and for global residual calculation
   
   // these are necessary for the GKB
   class ConstraintOpT,  
   class ConstraintOp,   
   class AugmentedLagrangianOp,
   class AugmentedLagrangianSolver
>
class GKBSolver : public Solver< SaddlePointOp >
{
 public:

   // source and destination functions for the GKB with m,n and n+m DoFs
   // n are the pressure DoFs, m are the velocity DoFs
   typedef typename SaddlePointOp::srcType nmFunction;
   typedef typename ConstraintOp::srcType mFunction;
   typedef typename ConstraintOp::dstType nFunction;

   GKBSolver(
       const std::shared_ptr< PrimitiveStorage >& storage,
       uint_t                                     level,
       AugmentedLagrangianSolver                  aLS,
       uint_t                                     mI             = std::numeric_limits< uint_t >::max(),
       real_t                                     tol            = 1e-5,
       real_t                                     Gamma          = 0, // Augmented Lagrangian Parameter
       real_t                                     S              = 0.2, 
       uint_t                                     Delay          = 5
   )  :  flag( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   , printInfo( false )
   , tolerance( tol )
   , maxIter( mI )
   , timingTree( storage->getTimingTree())
   , ALSolver(aLS)  // make this default PETScLUSolver => petsc build check   
   , gamma(Gamma)
   , s(S)
   , delay(Delay)
   , A(ConstraintOpT(storage, level, level))
   , AT(ConstraintOp(storage, level, level))
   , M(AugmentedLagrangianOp(storage,level,gamma))
   , u("u", storage, level, level), p("p", storage, level, level)
   , q("q", storage, level, level), v("v", storage, level, level)
   , d("d", storage, level, level), tmp_v("tmp_v", storage, level, level)
   , tmp_q("tmp_q", storage, level, level),tmp_w("tmp_w", storage, level, level)
   , dualr("dualr", storage, level, level)
   , z(std::vector<real_t>(delay,0))
   {
   }

   void solve( const SaddlePointOp& A, const nmFunction& x, const nmFunction& b, const uint_t level ) override
   {
      if ( maxIter == 0 )
         return;

      if ( x.isDummy() || b.isDummy() )
         return;

      //TODO check for uvw and p attribute in x and b, check if functions match between SaddlePointOp and sub operators, need concepts

      // set up rhs side for p and u
      mFunction f = b.uvw;
      nFunction g = b.p;


      timingTree->start( "GKBSolver" );

      // Init: first q 
      ALSolver.solve(M,u,f,level);
      AT.apply(u,tmp_q,level,flag);
      dualr.assign({1,-1},{g,tmp_q},level,flag);
      q.assign({1/gamma},{dualr},level,flag);
      beta = sqrt(dualr.dotGlobal(q,level,flag));

       
      std::cout << "Init beta=" << beta << std::endl;
      
      uint k = 0;
      bool convergence = false;
      //while()



      timingTree->stop( "GKBSolver" );
   }

   void setPrintInfo( bool pI ) { printInfo = pI; }
   
   /*
   void ComputeNextQ() {
      AT.apply(v, tmp_q, level, flag);
      q.assign({1/gamma,-alpha/gamma}, {tmp_q,q}, level, flag);
      beta = sqrt(q.dotGlobal(q, level, flag)/gamma);
      q.assign({1/beta},{q}, level, flag);
   }

   void ComputeNextV() {
      A.apply(q, tmp_v, level, flag);
      M.apply(v, tmp_v2, level, flag);
      tmp_v.assign({1,-beta},{tmp_v,tmp_w}, level, flag);
      ALSolver.solve(M,tmp_w,tmp_v,level);
      M.apply(tmp_w,tmp_v,level,flag);
      alpha = sqrt(tmp_v.dotGlobal(tmp_w,level,flag));
      v.assign({1/alpha},{tmp_w},level,flag);

   }
   */
   void SolutionUpdate() {
      
   }

  
   void Mnorm() {

   }
   


 private:

   // general parameters
   hyteg::DoFType                            flag;
   bool                                      printInfo;
   real_t                                    tolerance;
   uint_t                                    maxIter;   
   std::shared_ptr< walberla::WcTimingTree > timingTree;

   // GKB related 
   // parameters
   real_t                                    gamma; // AL parameter
   real_t                                    s;   // for upper bound
   uint_t                                    delay; // delay for lower bound

   // operators
   ConstraintOp                              AT;
   ConstraintOpT                             A;
   AugmentedLagrangianOp                     M;
   AugmentedLagrangianSolver                 ALSolver;

   // scalars holding M and N norms
   real_t                                    alpha;
   real_t                                    beta;
   
   // functions
   mFunction                                 u; // velocity solution
   nFunction                                 p; // pressure solution
   mFunction                                 v; // basis vector for velocity solution in K(AA')
   nFunction                                 q; // basis vector for pressure solution in K(A'A)                                    
   nFunction                                 d; // auxiliary vector for update of p
   std::vector<real_t>                       z; // coefficients of velocity in v_i, must be recorded for lower bound stopping critera

   // temporary functions
   mFunction                                 tmp_v; 
   mFunction                                 tmp_w; 
   nFunction                                 tmp_q;  
   nFunction                                 dualr;                                
   
};


// specification of GKB solver for P2P1TaylorHood Finite Elements with Sparse LU as inner solver
using  GKBSolver_P2P1THOP = GKBSolver< 
   P2P1TaylorHoodStokesOperator, 
   P1ToP2ConstantDivTOperator,
   P2ToP1ConstantDivOperator, 
   ALOP_P2P1TH, 
   PETScLUSolver<ALOP_P2P1TH>
>; 
} // namespace hyteg
