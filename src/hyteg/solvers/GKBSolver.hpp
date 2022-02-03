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
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/operators/AugmentedLagrangianOperator.hpp"
#include "hyteg/operators/BlockOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

using walberla::real_t; 
using walberla::uint_t;

/// Golub-Kahan-Bidiagonalization Solver implementation from "GENERALIZED GOLUB–KAHAN BIDIAGONALIZATION AND STOPPING CRITERIA" by Mario Arioli, 2013 
/*  solves 
   K x = b, in block form:

   M   A  * u  = f 
   A^T 0    p    g

   with K, an (m + n) x (m + n) matrix
   and with M (m x m) the Augmented Lagrangian matrix: M = Laplace + 1/gamma*div^T*div
   and modified right hand sides f,g (m,n vector) 
   f= //TODO add real AL approach (not necessary yet for P2P1 Stokes)
*/

template<
   // this op is not needed in GKB(not monolithic) but added in order to fit the solver class hierarchy and for global residual calculation
   class SaddlePointOp,
   
   // these are necessary for the GKB
   class ConstraintOpT, 
   class ConstraintOp,  // divergence for stokes
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

   // check if given operators are consistent
   static_assert(std::is_same<typename ConstraintOp::srcType,typename  ConstraintOpT::dstType>::value == true);
   static_assert(std::is_same<typename ConstraintOp::dstType,typename  ConstraintOpT::srcType>::value == true);
   static_assert(std::is_same<typename AugmentedLagrangianOp::dstType,typename  ConstraintOp::srcType>::value == true);
   static_assert(std::is_same<typename AugmentedLagrangianOp::srcType,typename  ConstraintOp::srcType>::value == true);

   GKBSolver(
       const std::shared_ptr< PrimitiveStorage >& storage,
       uint_t                                     level,
       AugmentedLagrangianSolver                  innerSolver,
       real_t                                     Nu             = 0, 
       uint_t                                     maxIt          = 30,
       // GKB stops if one of the tolerances is reached for one of the quantities: residual or error in energy norm/M norm
       real_t                                     merror_tol     = 1e-10, 
       real_t                                     res_tol        = 1e-10,
       // Augmented Lagrangian Parameter
//     real_t                                     S              = 0.2, // upper bound in check 3
       uint_t                                     Delay          = 5
   )  :  flag(  hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   , printInfo( true )
   , res_tolerance( res_tol )
   , merror_tolerance( merror_tol )
   , maxIter( maxIt )
   , timingTree( storage->getTimingTree())
   , ALSolver(innerSolver)  
   , nu(Nu)
   //, s(S)
   , delay(Delay)
   , A(ConstraintOpT(storage, level, level))
   , AT(ConstraintOp(storage, level, level))
   , M(AugmentedLagrangianOp(storage,level,Nu))
   //TODO check for unnecessary auxiliary functions
   , u("u", storage, level, level), p("p", storage, level, level)
   , q("q", storage, level, level), v("v", storage, level, level)
   , d("d", storage, level, level), tmp_v("tmp_v", storage, level, level)
   , tmp_q("tmp_q", storage, level, level),tmp_w("tmp_w", storage, level, level)
   , dualr("dualr", storage, level, level), globalR("globalR", storage, level, level)
   , globalX("globalX", storage, level, level), globalTmp("globalTmp", storage, level, level)
   , z(std::vector<real_t>(maxIt,0))
   {  
   }

   void solve( const SaddlePointOp& K, const nmFunction& x, const nmFunction& b, const uint_t Level ) override
   {
      level = Level;
   

      //TODO get AL matrix, print, check vs matlab AL matrix

      // copy boundary conditions
      mFunction u0 = x.uvw;
      nFunction p0 = x.p;
      CopyBCs(u0,p0);

      // set up rhs side for p and u
      mFunction f = b.uvw;
      nFunction g = b.p;      
      
      if(nu < 1e-14) {
         WALBERLA_LOG_INFO_ON_ROOT("Augmented Lagrangian Approach switched off.");
         nu = 1;
      } else {
         WALBERLA_LOG_INFO_ON_ROOT("Applying augmented Lagrangian Approach in GKB with nu = " << nu);
         A.apply(g,tmp_v,level,flag);
         f.assign({1,nu},{f,tmp_v},level,flag);
      }



      

      timingTree->start( "GKBSolver" );

      // Init: first q 
      ALSolver.solve(M,u,f,level);                                 // u = M^-1*f
      AT.apply(u,tmp_q,level,flag);                                 // store A'*u
      dualr.assign({1,-1},{g,tmp_q},level,flag);                    // r = g - A'*u
      q.assign({1/nu},{dualr},level,flag);                       // q = N^-1*r
      beta = sqrt(dualr.dotGlobal(q,level,flag));                   // beta = ||r'*q||_2
      q.assign({1/beta},{q},level,flag);                            // q = q/beta
      

      // Init: first v
      A.apply(q,tmp_v,level,flag);                                  // store A*q
      ALSolver.solve(M,v,tmp_v,level);                             // v = M^-1*A*q
      M.apply(v,tmp_v,level,flag);                                  // for Mnorm of v
      alpha = sqrt(v.dotGlobal(tmp_v,level,flag));                  // alpha = ||v||_M
      v.assign({1/alpha},{v},level,flag);                           // v = v/alpha
      
      d.assign({1/alpha},{q},level,flag);                           // d = q/alpha
      z.at(0) = (beta/alpha);                                      // z0 = beta/alpha
      
      
      // initial iterates of u,p
      u.assign({1,z.at(0)},{u,v},level,flag);                       // u = u + z0*v
      p.assign({1,-z.at(0)},{p,d},level,flag);                      // p = p -z0*d

      
      ////// main loop /////////
      bool convergence = false;
      while(k < maxIter - 1 && !convergence) {
            k += 1;
            ComputeNextQ();
            ComputeNextV();
            SolutionUpdate();
            convergence = StoppingCrit(K,b,f,g);
      }
      timingTree->stop( "GKBSolver" );

      if(printInfo) {
         if(k == maxIter-1) {
            WALBERLA_LOG_INFO_ON_ROOT("GKB did not converge to MErrtol=" << merror_tolerance<< "or Rtol="<<res_tolerance<<" in "<<maxIter<<" iterations.");
         } else {
            WALBERLA_LOG_INFO_ON_ROOT("GKB converged in "<<k<<" iterations.");
         }
      }

      x.uvw.assign({1},{u},level,flag);
      x.p.assign({1},{p},level,flag);
   }


   void setPrintInfo( bool pI ) { printInfo = pI; }
   
   void ComputeNextQ() {
      AT.apply(v, tmp_q, level, flag, Replace);                             // Store A'*v
      q.assign({1/nu,-alpha}, {tmp_q,q}, level, flag);          // Store q = 1/gamma*(A'*v -alpha*gamma*q)
      beta = sqrt(q.dotGlobal(q, level,flag)*nu);              // beta = ||q||_N
      q.assign({1/beta},{q}, level, flag);                         // q = q/beta
   }

   void ComputeNextV( ) {
      A.apply(q, tmp_v, level, flag, Replace);                              // Store A*q
      M.apply(v, tmp_w, level,flag, Replace);                             // Store M*V
      tmp_v.assign({1,-beta},{tmp_v,tmp_w}, level, flag);          // Store (A*q - beta*M*v)
      ALSolver.solve(M,v,tmp_v,level);                             // v = M^-1*(A*q - beta*M*v) TODO: Optimize M^-1*M
      M.apply(v,tmp_v,level, flag, Replace); 
      alpha = sqrt(tmp_v.dotGlobal(v,level,flag));                 // alpha = ||v||_M
      v.assign({1/alpha},{v},level,flag);                          // v = v/alpha
   }
   
   void SolutionUpdate( ) {
      z.at(k) = -beta/alpha*z.at((k-1));                              // z_k+1 = -beta/alpha*z_k
      d.assign({1/alpha,-beta/alpha},{q,d},level,flag );            // d = 1/alpha*(q - beta*d)
      u.assign({1,z.at(k)},{u,v},level,flag);                       // u = u + z*v
      p.assign({1,-z.at(k)},{p,d},level,flag);                          // p = p -z0*d
   }

   void CopyBCs(const mFunction& u0, const nFunction& p0) {
      u.copyBoundaryConditionFromFunction(u0);
      v.copyBoundaryConditionFromFunction(u0);
      tmp_v.copyBoundaryConditionFromFunction(u0);
      tmp_w.copyBoundaryConditionFromFunction(u0);
   	p.copyBoundaryConditionFromFunction(p0);
      d.copyBoundaryConditionFromFunction(p0);
      tmp_q.copyBoundaryConditionFromFunction(p0);
      dualr.copyBoundaryConditionFromFunction(p0);
      q.copyBoundaryConditionFromFunction(p0);
   }

   real_t ResidualNorm(const SaddlePointOp& K, const nmFunction& b,const mFunction& f,const nFunction& g) {
      globalX.uvw.assign({1},{u},level,flag);
      globalX.p.assign({1},{p},level,flag);
      // TODO: non monolithic residual calculation with incorporation of M
      M.apply(u,tmp_v,level,flag);
      A.apply(p,tmp_w,level,flag);
      globalR.uvw.assign({1,-1,-1},{f,tmp_v,tmp_w},level,flag);
      AT.apply(u,tmp_q,level,flag);
      globalR.p.assign({1,-1},{g,tmp_q},level,flag);
      
      //K.apply(globalX,globalTmp,level,Inner|NeumannBoundary);
      //globalR.assign({1,-1},{b,globalTmp},level,flag);
      real_t resnorm = sqrt(globalR.dotGlobal(globalR,level,flag));
      if(printInfo)
         WALBERLA_LOG_INFO_ON_ROOT("Iteration " << k << " ||r||="<< resnorm);
      return resnorm;
   }

   bool StoppingCrit(const SaddlePointOp& K, const nmFunction& b,const mFunction& f,const nFunction& g) {
      if(ResidualNorm(K,b,f,g) < res_tolerance) return true;
      //TODO add lower and upper bound stopping criterium
      else return false;
      
   }
  
   


 private:

   // general parameters
   hyteg::DoFType                            flag;
   bool                                      printInfo;
   real_t                                    merror_tolerance;
   real_t                                    res_tolerance;
   uint_t                                    maxIter;   
   std::shared_ptr< walberla::WcTimingTree > timingTree;
   uint_t                                    level;
   uint_t                                    k = 0;

   // GKB related 
   // parameters
   real_t                                    nu; // AL parameter
//   real_t                                    s;   // for upper bound
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

   // functions for global residual computation
   nmFunction                                globalR; 
   nmFunction                                globalX;
   nmFunction                                globalTmp;

   // temporary functions
   mFunction                                 tmp_v; 
   mFunction                                 tmp_w; 
   nFunction                                 tmp_q;  
   nFunction                                 dualr;                                
   
};


// specification of GKB solver for P2P1TaylorHood Finite Elements with CG as inner solver
// without AL approach!
using  GKBSolver_P2P1TH_NO_AL = GKBSolver< 
   P2P1TaylorHoodStokesOperator, 
   P1ToP2ConstantDivTOperator,
   P2ToP1ConstantDivOperator, 
   P2ConstantVectorLaplaceOperator, 
   CGSolver<P2ConstantVectorLaplaceOperator> 
>; 

// specification of GKB solver for P2P1TaylorHood Finite Elements with CG as inner solver
// now with AL
using  GKBSolver_P2P1TH = GKBSolver< 
   P2P1TaylorHoodStokesOperator, 
   P1ToP2ConstantDivTOperator,
   P2ToP1ConstantDivOperator, 
   ALOP_P2P1TH, 
   CGSolver<ALOP_P2P1TH> 
>;
} // namespace hyteg
