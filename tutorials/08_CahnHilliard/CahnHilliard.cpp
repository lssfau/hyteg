/*
 * Copyright (c) 2021 Andreas Wagner.
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

/**
 * \page 08_CahnHilliard A full app for the Cahn-Hilliard equation.
 *
 * \dontinclude tutorials/08_CahnHilliard/CahnHilliard.cpp
 *
 * \brief In this tutorial, we will use the LinearCombinationForm to define specialized bilinear forms for the Cahn-Hilliard (CH) equations.
 *        In combination with block functions and block operators, we will discretize the CH-Equations.
 *        We will speed up our solvers by implementing a simple diagonal preconditioner.
 *
 * We already showed how to solve the Poisson- and the Stokes-Equation, which already are well supported by HyTeG.
 * In this tutorial, we will implement a solver for a completely unrelated problem, which is not a core application for HyTeG (yet ;)).
 * In addition, we will learn about block-functions, block-operators, and the composition of finite-element forms.
 *
 * Before we will dive into coding, we will first introduce the equations, motivate our time-integrator, talk about the space-discretization
 * and finally present a simple diagonal preconditioner.
 *
 * \section Problem
 *
 * As an example, we will try to solve the Cahn-Hilliard equation,
 * which simulates the phase separation of a mixture of two fluids.
 * For convenience, let us call these fluids A and B.
 * They are represented by an indicator function \f$\phi \in [-1, +1]\f$,
 * such that \f$\phi(x) = +1\f$ if only fluid A is present and \f$\phi(x) = -1\f$  if only fluid B is present.
 * All intermediate values correspond to mixtures of the fluid, where either A (\f$\phi(x) > 0\f$) or B dominates (\f$\phi(x) < 0\f$).
 *
 * We want to assign an energy functional to these equations given by
 * \f{align*}{
 *     F(\phi) = \int \boldsymbol \psi(\phi) dx + \frac{\epsilon^2}{2}\int |\nabla \phi |^2 dx
 * \f}
 * where \f$\psi(\phi) = \frac{1}{4}(1-\phi^2)^2\f$ is a double-well potential, which penalizes mixed states,
 * while the second term models a surface energy and encourages small interfaces.
 * A gradient-flow for this functional is given by the Cahn-Hilliard equations
 * \f{align*}{
 *     \partial_t \phi &= \Delta \mu, \\
 *     \mu &= \boldsymbol \psi'(\phi) - \epsilon^2 \Delta \phi \quad \text{ and } \\
 *     \partial_n \phi &= \partial_n \mu = 0 \quad \text{ on } \quad \partial \Omega
 *     ,
 * \f}
 * which minimizes \f$F(\phi(t))\f$ as time passes, i.e. \f$F(\phi(s)) \leq F(\phi(t))\f$ for \f$t \leq s\f$.
 * The equations are a mixed system of the fields \f$\phi\f$ and \f$\mu\f$.
 *
 * \section Time-discretization
 *
 * The energy minimization of $F$ is usually something we want to enforce by our time discretization.
 *
 * For this, we split up the potential $\boldsymbol \psi$ into a convex and a concave part, such that
 * \f{align*}{
 *     \boldsymbol \psi = \boldsymbol \psi_c - \boldsymbol \psi_e
 * \f}
 * where \f$\boldsymbol \psi_e\f$, \f$\boldsymbol \psi_c\f$ are convex.
 * It can be shown that if we treat \f$\psi_e\f$ explicitly and \f$\psi_c\f$ implicitly the resulting scheme will be unconditionally gradient stable.
 * That means that \f$F(\phi^{n+1}) \leq F(\phi^{n})\f$ for our time-discrete solutions \f$\phi^{(n)}\f$ independent of the time-step size \f$\tau\f$.
 * The resulting scheme is
 *
 * \f{align}{
 *     \phi^{(n+1)}-\phi^{(n)} &= \tau \Delta \mu \\
 *     \mu &= \boldsymbol \psi_c'(\phi^{(n+1)}) - \boldsymbol \psi_e'(\phi^{(n)}) - \epsilon^2 \Delta \phi^{(n+1)}
 *     .
 * \f}
 *
 * The decomposition of the potential in a convex and concave function is typically not unique.
 * A possible choice is
 * \f{align}{
 *     \boldsymbol \psi_c(\phi) = \frac{1}{4}
 *     \left( 4 \phi^2 + 1 \right)
 *     \quad \text{ and } \quad
 *     \boldsymbol \psi_e(\phi) = \frac{1}{4}
 *     \left( 6 \phi^2 - \phi^4 \right)
 * \f}
 * yielding the derivatives
 * \f{align}{
 *     \boldsymbol \psi_c'(\phi) =
 *     2 \phi
 *     \quad \text{ and } \quad
 *     \boldsymbol \psi_e'(\phi) =
 *     3 \phi - \phi^3
 *     .
 * \f}
 *
 * \section Space-discretization
 *
 * Applying a FEM discretization on the equations yields the following system of equations
 * \f{align}{
 *     \begin{pmatrix}
 *         \tau K & M\\
 *         M & - \epsilon^2 K - 2 M
 *     \end{pmatrix}
 *     \begin{pmatrix}
 *         \boldsymbol \mu \\
 *         \boldsymbol{\phi}^{(n+1)}
 *     \end{pmatrix}
 *     =
 *     \begin{pmatrix}
 *         \boldsymbol f^{(n)} \\ \boldsymbol g^{(n)}
 *     \end{pmatrix}
 *     ,
 * \f}
 * with the stiffness matrix \f$K\f$, the mass matrix \f$M\f$ and the right-hand side given by
 * \f{align}{
 *     \boldsymbol{f}^{(n)}_i = \left\langle\phi^{(n)}_h, \nu_{h,i}\right\rangle_{L_2}
 * \f}
 * and
 * \f{align}{
 *     \boldsymbol g^{(n)}_i = \left\langle-3 \phi^{(n)}_h + (\phi^{(n)}_h)^3, \psi_{h,i} \right\rangle_{L_2}
 * \f}
 * with test functions \f$\nu_{h,i}\f$ and \f$\psi_{h,i}\f$.
 * We consider a \f$\mathcal P^1\f$ (continuous, piecewise linear polynomials) discretization of FE-shape functions.
 *
 *
 * \section Preconditioner
 *
 * To solve our problem efficiently we will need a preconditioner.
 * Here, we will follow the approach in [BRENNER2018] and get simply by rescaling the equations
 * \f{align}{
 *     \begin{pmatrix}
 *         \tau K & M\\
 *         M & - \epsilon^2 K - 2 M
 *     \end{pmatrix}^{-1}
 *     =
 *     \begin{pmatrix}
 *         \tau^{-1/4} & \\
 *         & \tau^{1/4}
 *     \end{pmatrix}
 *     \begin{pmatrix}
 *         \sqrt{\tau} K & M\\
 *         M & - \sqrt \tau \epsilon^2 K - 2 \sqrt \tau M
 *     \end{pmatrix}^{-1}
 *     \begin{pmatrix}
 *         \tau^{-1/4} & \\
 *         & \tau^{1/4}
 *     \end{pmatrix}
 *     .
 * \f}
 * We denote the middle matrix on the right-hand side by
 * \f{align}{
 *     B:=
 *     \begin{pmatrix}
 *         \sqrt{\tau} K & M\\
 *         M & - \sqrt \tau \epsilon^2 K - 2 \sqrt \tau M
 *     \end{pmatrix}
 *     .
 * \f}
 * As shown in [Brenner2018] \f$P^{-1} B\f$ can be solved efficiently with MINRES, where
 * \f{align}{
 *     P^{-1}
 *     :=
 *     \begin{pmatrix}
 *         (\sqrt{\tau} K + M)^{-1} & \\
 *         &
 *         (\sqrt{\tau} \epsilon^2 K + M)^{-1}
 *     \end{pmatrix}
 * \f}
 * is a diagonal preconditioner.
 *
 * To get an optimal solver, we can use a geometric multigrid to invert the diagonal blocks.
 *
 * \section Implementation
 *
 * \subsection defining-the-block-function-space Defining the block function space
 *
 * We start with defining a block function representing the state of our system at a given time step by inheriting from the hyteg::BlockFunction.
 * The first component will correspond to the chemical potential, while the second is the indicator function for our mixtures.
 *
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction definition
 *
 * The two components for the \f$\phi \f$ and \f$\mu\f$ are defined inside the constructor by pushing a hyteg::P1Function onto the
 * `subFunc_` vector.
 *
 * For our convenience, we define two helper functions `getMu` and `getPhi`, which return the subcomponents of our vector.
 * In principle, these be accessed with hyteg::BlockFunction< value_t >::getSubFunction, which returns a hyteg::GenericFunction.
 * To get the underlying hyteg::P1Function we have to "cast" it to this type using the hyteg::GenericFunction::unwrap method,
 * which internally delegates this task to the hyteg::FunctionWrapper.
 *
 * The compiler needs to be able to infer some static information about each function.
 * For instance, its name, its value type, and it needs a tag for certain meta programming features.
 * This can be achieved by extending the given trait system:
 *
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction meta-information
 *
 * If one of HyTeG's internal algorithms needs to know the value type of an unknown function type `UnknownType`, it can do this by calling
 * `hyteg::FunctionTrait< UnknownType >::ValueType` at compile time.
 * This is implemented for hyteg::P1Function, hyteg::P2Function, hyteg::P2P1TaylorHoodFunction, ... and now also for our
 * P1CahnHilliardFunction.
 *
 * A more in-depth tutorial about the rationality behind traits can be found at [TRAITS].
 *
 * To keep our type names within "reasonable bounds" we will introduce the shorthand `chType` for our `P1CahnHilliardFunction`:
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction chType-definition
 *
 * \subsection rhs-assembly The right-hand-side assembly
 *
 * Next, we will implement the assembly of the right-hand-side vector
 * \f{align}{
 *     \boldsymbol{f}^{(n)}_i = \left\langle\phi^{(n)}_h, \nu_{h,i}\right\rangle_{L_2}
 *     \quad \text{ and } \quad
 *     \boldsymbol g^{(n)}_i = \left\langle-3 \phi^{(n)}_h + (\phi^{(n)}_h)^3, \psi_{h,i} \right\rangle_{L_2}
 * \f}
 * which has to happen after every time step using the updated value for our indicator function \f$\phi^{(n)}_h\f$.
 *
 * We will hide the assembly process in the RhsAssembler class, given by
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction ClassRhsAssembler
 *
 * The `RhsAssembler::assemble` method, will take the solution of the previous time step `u_prev` as an input value
 * and write the result into the vector `rhs`.
 *
 * Since HyTeG does not support yet the assembly of vectors this will have to happen with a matrix-vector multiplication.
 * This is certainly not the most efficient way to implement this, but we will see in the final simulation that the additional
 * costs are negligible compared to solving the system.
 *
 * For an arbitrary scalar function \f$\kappa\f$, we define the bilinear form \f$m_\kappa(\nu_j, \nu_i) := \int \kappa \nu_i \nu_j \text{d}x\f$
 * for test functions \f$\nu_i, \nu_j\f$ and the associated matrix \f$M_\kappa\f$ given by \f$(M_{\kappa})_{i,j} = m_\kappa(\nu_j, \nu_i)\f$.
 * The right-hand-side can now be written with these generalized mass-matrices as
 * \f{align}{
 *  \boldsymbol{f}^{(n)} = M_1 \boldsymbol{\phi^{(n)}}
 *  \quad \text{ and } \quad \boldsymbol{g}^{(n)} = M_{\tilde \kappa} \boldsymbol{\phi^{(n)}}
 *  \quad \text{ with } \quad \tilde \kappa = -3 + (\phi^{(n)}_h)^2
 *  .
 * \f}
 *
 * In HyTeG \f$m_\kappa\f$ is given by hyteg::forms::p1_k_mass_affine_q4.
 * Here
 * - `p1` stands for the used function space of piecewise linear continuous functions,
 * - `k_mass` is the common name of the bilinear form,
 * - `affine` stresses that it does not support geometry blending and
 * - `q4` denotes the degree of the quadrature rule.
 *
 * It has two scalar callback functions \f$\kappa\f$ as constructor arguments, which correspond to \f$\kappa\f$ in 2D and 3D which are allowed to be different.
 *
 * The operator \f$M_{\kappa}\f$ is realized in the hyteg::P1ElementwiseKMassOperator, and takes the given bilinear form in its constructor.
 * The assembly code becomes
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction assemble
 * To calculate \f$\boldsymbol{g}\f$, we first define a lambda function `kappa` which takes a quadrature point and returns the value of \f$\tilde \kappa\f$ at that point.
 * For this, we get the \f$\phi\f$-component of `u_prev` and evaluate it at `p` with hyteg::P1Function::evaluate.
 * The functional is inserted into the form for \f$m_{\tilde \kappa}\f$, which is used to initialize the operator `op_m_kappa_tilde` or \f$M_{\tilde \kappa}\f$ in our previous notation.
 * Since the \f$\tilde \kappa\f$ is dimension-independent, we use it both for 2D and 3D.
 * Finally, we apply the operator to the \f$\phi\f$-component of `u_prev` and store the result in the \f$\phi\f$-component of `rhs`.
 *
 * The calculation of \f$\boldsymbol{f}\f$ is simpler: First, we define the usual mass-matrix in `op_m_1`.
 * We apply it to the \f$\phi\f$-component of `u_prev` and store the result in the \f$\mu\f$-component of `rhs`.
 *
 * \subsection lhs-assembly The left-hand-side assembly
 *
 * Next, we want to discuss the assembly of the matrix we want to invert, namely
 * \f{align}{
 *     \begin{pmatrix}
 *         \sqrt{\tau} K & M\\
 *         M & - \sqrt \tau \epsilon^2 K - 2 \sqrt \tau M
 *     \end{pmatrix}
 * \f}
 * All these operators are not yet defined in HyTeG and have to be constructed from their respective bilinear forms.
 * But since they are all linear combinations of the mass- and stiffness matrix, we can use the hyteg::P1LinearCombinationForm to construct
 * them from the mass and stiffness forms already defined in HyTeG. The mass and stiffness forms are
 *
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction forms-definitions
 *
 * Here, they are given as `hyteg::P1FenicsForm`s, with the first type-argument for the 2D form and the second type-argument for the 3D form.
 * The form for the lower-right block can now be constructed by
 *
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction block-form-11
 *
 * We first define the hyteg::P1LinearCombinationForm, which defaults to the zero form.
 * We then add the stiffness form scaled by \f$-\epsilon^2\sqrt{\tau}\f$.
 * Finally, we add the mass scaled with \f$ - 2 \sqrt \tau \f$.
 *
 * The upper-left block is a bit simpler to construct since it only requires to scale the stiffness form:
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction block-form-00
 *
 * We now have the correct forms, but still, need to construct proper operators from them.
 * We could do this explicitly by
 * ```
 * hyteg::P1ConstantLinearCombinationOperator operator00( storage, minLevel, maxLevel, form00)
 * ```
 *
 * For the off-diagonal blocks no additional forms are needed and we can resort directly to the hyteg::P1ConstantMassOperator.
 *
 * We construct the complete block-operator in the following factory-function:
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction lhs-assembly
 *
 * The hyteg::BlockOperator takes the dimension of the block-system in the last two arguments.
 * By calling the hyteg::BlockOperator::createSubOperator method we can construct the blocks.
 * The type of the block operator is given in the angled brackets, followed by the two indices of the respective block and the arguments for the constructor of the operator.
 *
 * For the diagonal blocks, we construct a hyteg::P1ConstantLinearCombinationOperator as described before
 * and pass the correct linear combination of forms in its constructor.
 *
 * This is all we have to do to get a working operator in HyTeG.
 *
 * \subsection time-evolution The time-evolution operator
 *
 * We will postpone the discussion on how to invert the given operator for now and concentrate on the time-evolution, i.e.
 * given \f$ \boldsymbol u^{(n)} = (\boldsymbol \mu^{(n)}, \boldsymbol \phi^{(n)}) \f$ we want to calculate the values at the next time step
 * \f$ \boldsymbol u^{n+1} = (\boldsymbol \mu^{n+1}, \boldsymbol \phi^{n+1}) \f$.
 * This will be handled by the following operator
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction ClassCahnHilliardEvolutionOperator
 *
 * Note that `lhsOp_` is the operator we constructed in the previous step and ignore all the references to solvers for now.
 * The apply method takes function `u_prev` and "evolves" it to `u_now` on the given level.
 * This method thus executes the operation
 * \f{align}{
 *     \begin{pmatrix}
 *          \boldsymbol \mu^{(n+1)} \\
 *          \boldsymbol \phi^{(n+1)}
 *     \end{pmatrix}
 *     =
 *     \begin{pmatrix}
 *         \tau^{-1/4} & \\
 *         & \tau^{1/4}
 *     \end{pmatrix}
 *     \begin{pmatrix}
 *         \sqrt{\tau} K & M\\
 *         M & - \sqrt \tau \epsilon^2 K - 2 \sqrt \tau M
 *     \end{pmatrix}^{-1}
 *     \begin{pmatrix}
 *         \tau^{-1/4} & \\
 *         & \tau^{1/4}
 *     \end{pmatrix}
 *     \begin{pmatrix}
 *          \boldsymbol f^{(n)} \\
 *          \boldsymbol g^{(n)}
 *     \end{pmatrix}
 * \f}
 * where \f$(\boldsymbol f^{(n)}, \boldsymbol g^{(n)})\f$ depend on \f$(\boldsymbol \mu^{(n)}, \boldsymbol \phi^{(n)})\f$.
 * This means, that it has to
 * - assemble the right-hand side,
 * - scale the vector components with the diagonal matrix,
 * - invert the `lhsOp_` matrix,
 * - and scale again the result with the diagonal matrix.
 *
 * The multiplication with the diagonal matrix is done implicitly in the `scale` method by using the hyteg::P1Function::assign on both components
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction scale
 *
 * The apply method can thus be written as
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction time-evolution-apply
 * which is essentially the list of operations given before, wrapped by some calls to the timing-tree API.
 * This is especially useful to get a feeling of how long solving the system takes compared to assembling the right-hand side.
 *
 * \subsection minres-solver The solver
 *
 * We will now discuss how to invert the `lhsOp_` Operator from the previous step.
 * As an outer solver we will use a diagonally preconditioned MINRES, which is configured in the `create_solver` factory function,
 * which was called in the constructor of the CahnHilliardEvolutionOperator of the previous step.
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction minres-solver
 * First, we construct the CahnHilliardDiagonalPreconditioner, which we will implement in a second and will have the usual hyteg::Solver interface.
 * The preconditioner is passed to the constructor of the hyteg::MinResSolver, together with the maximum number of iteration and the absolute tolerance of the residual.
 * The default preconditioner is the hyteg::IdentityPreconditioner which is commented out.
 *
 * For the diagonal preconditioner, we will have to evaluate the matrix
 * \f{align}{
 *     P^{-1}
 *     :=
 *     \begin{pmatrix}
 *         (\sqrt{\tau} K + M)^{-1} & \\
 *         &
 *         (\sqrt{\tau} \epsilon^2 K + M)^{-1}
 *     \end{pmatrix}
 *     .
 * \f}
 * We will use a multigrid algorithm to invert the diagonal blocks.
 *
 * The CahnHilliardDiagonalPreconditioner inherits from the solver interface:
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardFunction CahnHilliardDiagonalPreconditioner
 * The CahnHilliardDiagonalPreconditioner::create_block_00_form and CahnHilliardDiagonalPreconditioner::create_block_11_form methods
 * will create the hyteg::P1LinearCombinationForm%s for \f$ \sqrt{\tau} K + M \f$
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardDiagonalPreconditioner block00form
 * and \f$\sqrt{\tau} \epsilon^2 K + M\f$
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardDiagonalPreconditioner block11form
 * which are used to construct the operators `block_00_operator` and `block_11_operator`.
 *
 * The CahnHilliardDiagonalPreconditioner::create_multigrid_solver initializes our multigrid algorithm in the same way as in \ref 05_FullAppP1GMG:
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardDiagonalPreconditioner create-multigrid-solver
 *
 * Now that we are finished with the setup phase, implementing the action of the preconditioner in the `solve` method is simple:
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp CahnHilliardDiagonalPreconditioner solve
 * We apply the multigrid solver first on the \f$\mu\f$-component of our vector `numVCyclesMu_`-times.
 * After that, we do the same for the lower-right block.
 *
 * \subsection domain The domain
 *
 * We will solve the Cahn-Hilliard equations on a rectangular or cubic domain, depending on whether we are in 2D or 3D.
 * A rectangular domain can be constructed with
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp create_storage rectangle
 * The only notable thing which changes compared to the previous tutorial is, that we need Neumann boundary conditions on the whole domain.
 * For this we use the hyteg::SetupPrimitiveStorage::setMeshBoundaryFlagsOnBoundary method.
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp create_storage neumann
 *
 * Since we want to seamlessly switch between 2D and 3D the whole construction of the primitive storage gets moved into an auxillary function,
 * which takes the boolean flag `use3D` as an argument:
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp create_storage fun
 * We stress, that this will be the only place in our code where different code paths are taken depending on the spatial dimension.
 *
 * \subsection main The main function
 *
 * We are still in need of the main function to start the simulation, which ties all the previously defined classes and functions together.
 * In principle, it should not contain anything surprising which was not encountered in the previous tutorials yet.
 *
 * We start with configuring walberla
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp main config-walberla
 *
 * We do not want to hardcode most of our parameters and hence have to read in a configuration file:
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp main config-file
 * As before `minLevel` and `maxLevel` refer to the coarsest and finest level of our multigrid algorithm.
 * Since we do not want to output our solution after every time step, we define an `outputInterval`, so that we only get every `outputInterval`th step.
 * The parameters `tau` and `epsilon` refer to the time step width and interface width.
 *
 * In practice choices for `maxLevel`, `tau` and `epsilon` will depend on one another,
 * since if we refine the mesh for a smaller spatial error we will also want to decrease `tau` since otherwise the temporal error would dominate.
 * And since we often want the interfacial width to be as small as possible, we will decrease this parameter too.
 *
 * Finally, the time when our simulation stops is given in `t_end`.
 *
 * We create the storage with our previously defined helper function
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp main storage
 *
 * Next we define functions for storing the state of our CH equation at the current and previous time step
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp main functions
 *
 * Next we setup the CahnHilliardEvolutionOperator as
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp main evolution-operator
 * We have added the flag `printMinresConvergence` to activate the verbose output mode of our MINRES solver.
 * This allows us to check how many iterations we need until convergence and how the iteration count changes for finer meshes.
 *
 * A common benchmark for a Cahn-Hilliard solver is to start with a homogeneous mixture of fluids that is disturbed by some random noise.
 * We will thus use walberla::math::realRandom to initialize `u_prev` with random values in the interval \f$[-0.01, +0.01]\f$
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp main initial-values
 * The seed for these random functions can be set with walberla::math::seedRandomGenerator.
 * We use the processor rank as an initial seed to avoid having the same random values on all processes.
 * This would lead to symmetric patterns in our solution, due to the homogeneous mesh.
 *
 * The actual time loop is fairly simple
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp main time-loop
 * We increment our current time `t_now` in every time step by `tau`, apply our CahnHilliardEvolutionOperator on `u_prev` to get
 * the updated state of our PDE-system in `u_now`.
 * This gets copied to u_prev with the hyteg::P1Function::assign method.
 * After `ouputInterval` steps, we add our current solution to the hyteg::VTKOutput and write the result to disk.
 *
 * When our simulation finally reaches `t_end` and finishes we might want to get additional information about the runtime of the application to identify bottlenecks.
 * For instance, whether it is worth optimizing our inefficient right-hand-side assembly.
 * This can be done by printing the walberla::timing::TimingTree of the hyteg::PrimitiveStorage to the console.
 * \snippet tutorials/08_CahnHilliard/CahnHilliard.cpp main timing-tree
 *
 * \section results Results in 2D
 *
 * \htmlonly
  <center>
  <table>
  <tr>
  <td><img src="cahn_hilliard_2d.00.png" width="100%"/><center>initial random state</center></td>
  <td><img src="cahn_hilliard_2d.01.png" width="100%"/><center>at t=0.05</center></td>
  </tr>
  <tr>
  <td><img src="cahn_hilliard_2d.02.png" width="100%"/><center>at t=0.25</center></td>
  <td><img src="cahn_hilliard_2d.04.png" width="100%"/><center>at t=0.95</center></td>
  </tr>
  </table>
  </center>
  \endhtmlonly
 *
 * \section code Code
 *
 * The full code is:
 *
 * \include tutorials/08_CahnHilliard/CahnHilliard.cpp
 *
 * \section references References
 *
 * [Brenner2018] BRENNER, Susanne C.; DIEGEL, Amanda E.; SUNG, Li-Yeng.
 *               <a href="https://arxiv.org/pdf/1811.03430.pdf">A robust solver for a mixed finite element method for the Cahn–Hilliard equation</a>.
 *               Journal of Scientific Computing, 2018, 77. Jg., Nr. 2, S. 1234-1249.
 *
 * [TRAITS] https://www.boost.org/doc/libs/1_53_0/libs/geometry/doc/html/geometry/design.html
 *
 */

#include <cmath>
#include <utility>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/operators/BlockOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

// forward declarations:
class CahnHilliardDiagonalPreconditioner;

/// [CahnHilliardFunction definition]
template < typename value_t >
class P1CahnHilliardFunction : public BlockFunction< value_t >
{
 public:
   template < typename VType >
   using FunctionType = P1CahnHilliardFunction< VType >;

   P1CahnHilliardFunction( const std::string&                         name,
                           const std::shared_ptr< PrimitiveStorage >& storage,
                           size_t                                     minLevel,
                           size_t                                     maxLevel )
   : BlockFunction< value_t >( name )
   {
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_mu", storage, minLevel, maxLevel ) );
      this->subFunc_.push_back(
          std::make_shared< FunctionWrapper< P1Function< value_t > > >( name + "_phi", storage, minLevel, maxLevel ) );
   };

   [[nodiscard]] const P1Function< value_t >& getMu() const
   {
      return this->getSubFunction( 0 ).template unwrap< P1Function< value_t > >();
   }

   [[nodiscard]] const P1Function< value_t >& getPhi() const
   {
      return this->getSubFunction( 1 ).template unwrap< P1Function< value_t > >();
   }
};
/// [CahnHilliardFunction definition]

/// [CahnHilliardFunction meta-information]
class P1CahnHilliardFunctionTag
{};

namespace hyteg {

template < typename VType >
struct FunctionTrait< P1CahnHilliardFunction< VType > >
{
   typedef VType                     ValueType;
   typedef P1CahnHilliardFunctionTag Tag;

   static std::string getTypeName() { return "CahnHilliardFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

} // namespace hyteg
/// [CahnHilliardFunction meta-information]

/// [CahnHilliardFunction chType-definition]
using chType = P1CahnHilliardFunction< real_t >;
/// [CahnHilliardFunction chType-definition]

/// [CahnHilliardFunction forms-definitions]
using MassFormType    = hyteg::P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >;
using LaplaceFormType = hyteg::P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise >;
/// [CahnHilliardFunction forms-definitions]

/// [CahnHilliardFunction minres-solver]
std::shared_ptr< MinResSolver< BlockOperator< chType, chType > > >
    create_solver( const std::shared_ptr< PrimitiveStorage >& storage,
                   uint_t                                     minLevel,
                   uint_t                                     maxLevel,
                   real_t                                     tau,
                   real_t                                     epsilon )
{
   auto preconditioner =
       std::make_shared< CahnHilliardDiagonalPreconditioner >( storage, minLevel, maxLevel, tau, epsilon, 2, 2 );

   // auto preconditioner = std::make_shared< IdentityPreconditioner< BlockOperator< chType, chType > > >();

   const uint_t maxIter   = std::numeric_limits< uint_t >::max();
   const real_t tolerance = 1e-8;

   auto solver = std::make_shared< MinResSolver< BlockOperator< chType, chType > > >(
       storage, minLevel, maxLevel, maxIter, tolerance, preconditioner );

   return solver;
}
/// [CahnHilliardFunction minres-solver]

/// [CahnHilliardFunction lhs-assembly]
std::shared_ptr< BlockOperator< chType, chType > > create_lhs_operator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                        uint_t                                     minLevel,
                                                                        uint_t                                     maxLevel,
                                                                        real_t                                     tau,
                                                                        real_t                                     epsilon )
{
   /// [CahnHilliardFunction block-form-00]
   P1LinearCombinationForm form00;
   form00.addOwnedForm< LaplaceFormType >( std::sqrt( tau ) );
   /// [CahnHilliardFunction block-form-00]

   /// [CahnHilliardFunction block-form-11]
   P1LinearCombinationForm form11;
   form11.addOwnedForm< LaplaceFormType >( -std::pow( epsilon, 2 ) * std::sqrt( tau ) );
   form11.addOwnedForm< MassFormType >( -2 * std::sqrt( tau ) );
   /// [CahnHilliardFunction block-form-11]

   auto op = std::make_shared< BlockOperator< chType, chType > >( storage, minLevel, maxLevel, 2, 2 );
   op->createSubOperator< P1ConstantLinearCombinationOperator >( 0, 0, storage, minLevel, maxLevel, form00 );
   op->createSubOperator< P1ConstantLinearCombinationOperator >( 1, 1, storage, minLevel, maxLevel, form11 );
   op->createSubOperator< P1ConstantMassOperator >( 0, 1, storage, minLevel, maxLevel );
   op->createSubOperator< P1ConstantMassOperator >( 1, 0, storage, minLevel, maxLevel );

   return op;
}
/// [CahnHilliardFunction lhs-assembly]

/// [CahnHilliardFunction ClassRhsAssembler]
class RHSAssembler
{
 public:
   RHSAssembler( std::shared_ptr< PrimitiveStorage > storage, uint_t minLevel, uint_t maxLevel, real_t epsilon )
   : storage_( std::move( storage ) )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , epsilon_( epsilon )
   {}

   void assemble( const chType& u_prev, const chType& rhs, uint_t level, DoFType flag );

 protected:
   std::shared_ptr< PrimitiveStorage > storage_;
   uint_t                              minLevel_;
   uint_t                              maxLevel_;

   real_t epsilon_;
};
/// [CahnHilliardFunction ClassRhsAssembler]

/// [CahnHilliardFunction assemble]
void RHSAssembler::assemble( const chType& u_prev, const chType& rhs, uint_t level, DoFType flag )
{
   auto kappa = [&]( const Point3D& p ) -> real_t {
      real_t phi_prev_value;
      WALBERLA_CHECK( u_prev.getPhi().evaluate( p, level, phi_prev_value ) );

      return -3. + std::pow( phi_prev_value, 2 );
   };

   forms::p1_k_mass_affine_q4 form( kappa, kappa );

   P1ElementwiseKMassOperator op_m_kappa_tilde( storage_, minLevel_, maxLevel_, form );

   op_m_kappa_tilde.apply( u_prev.getPhi(), rhs.getPhi(), level, flag );

   P1ElementwiseMassOperator op_m_1( storage_, minLevel_, maxLevel_ );

   op_m_1.apply( u_prev.getPhi(), rhs.getMu(), level, flag );
}
/// [CahnHilliardFunction assemble]

/// [CahnHilliardFunction CahnHilliardDiagonalPreconditioner]
class CahnHilliardDiagonalPreconditioner : public Solver< BlockOperator< chType, chType > >
{
 public:
   CahnHilliardDiagonalPreconditioner( const std::shared_ptr< PrimitiveStorage >& storage,
                                       uint_t                                     minLevel,
                                       uint_t                                     maxLevel,
                                       real_t                                     tau,
                                       real_t                                     epsilon,
                                       uint_t                                     numVCyclesMu,
                                       uint_t                                     numVCyclesPhi )
   : block_00_operator( storage, minLevel, maxLevel, create_block_00_form( tau ) )
   , block_11_operator( storage, minLevel, maxLevel, create_block_11_form( tau, epsilon ) )
   , solver( create_multigrid_solver( storage, minLevel, maxLevel ) )
   , numVCyclesMu_( numVCyclesMu )
   , numVCyclesPhi_( numVCyclesPhi )
   {}

   // y = M^{-1} * x
   void solve( const BlockOperator< chType, chType >&, const chType& x, const chType& b, const uint_t level ) override;

 protected:
   static P1LinearCombinationForm create_block_00_form( real_t tau );

   static P1LinearCombinationForm create_block_11_form( real_t tau, real_t epsilon );

   static std::shared_ptr< GeometricMultigridSolver< P1ConstantLinearCombinationOperator > >
       create_multigrid_solver( std::shared_ptr< PrimitiveStorage > storage, const uint_t minLevel, const uint_t maxLevel );

 protected:
   P1ConstantLinearCombinationOperator block_00_operator;
   P1ConstantLinearCombinationOperator block_11_operator;

   std::shared_ptr< GeometricMultigridSolver< P1ConstantLinearCombinationOperator > > solver;

   uint_t numVCyclesMu_;
   uint_t numVCyclesPhi_;
};
/// [CahnHilliardFunction CahnHilliardDiagonalPreconditioner]

/// [CahnHilliardDiagonalPreconditioner solve]
void CahnHilliardDiagonalPreconditioner::solve( const BlockOperator< chType, chType >&,
                                                const chType& x,
                                                const chType& b,
                                                const uint_t  level )
{
   for ( uint_t k = 0; k < numVCyclesMu_; k += 1 )
      solver->solve( block_00_operator, x.getMu(), b.getMu(), level );
   for ( uint_t k = 0; k < numVCyclesPhi_; k += 1 )
      solver->solve( block_11_operator, x.getPhi(), b.getPhi(), level );
}
/// [CahnHilliardDiagonalPreconditioner solve]

/// [CahnHilliardDiagonalPreconditioner block00form]
P1LinearCombinationForm CahnHilliardDiagonalPreconditioner::create_block_00_form( real_t tau )
{
   P1LinearCombinationForm form;
   form.addOwnedForm< LaplaceFormType >( std::sqrt( tau ) );
   form.addOwnedForm< MassFormType >( 1. );
   return form;
}
/// [CahnHilliardDiagonalPreconditioner block00form]

/// [CahnHilliardDiagonalPreconditioner block11form]
P1LinearCombinationForm CahnHilliardDiagonalPreconditioner::create_block_11_form( real_t tau, real_t epsilon )
{
   P1LinearCombinationForm form;
   form.addOwnedForm< LaplaceFormType >( std::sqrt( tau ) * std::pow( epsilon, 2 ) );
   form.addOwnedForm< MassFormType >( 1. );
   return form;
}
/// [CahnHilliardDiagonalPreconditioner block11form]

/// [CahnHilliardDiagonalPreconditioner create-multigrid-solver]
std::shared_ptr< GeometricMultigridSolver< P1ConstantLinearCombinationOperator > >
    CahnHilliardDiagonalPreconditioner::create_multigrid_solver( std::shared_ptr< PrimitiveStorage > storage,
                                                                 const uint_t                        minLevel,
                                                                 const uint_t                        maxLevel )
{
   const real_t coarse_tolerance = 1e-16;
   const uint_t max_coarse_iter  = 1000;

   auto smoother         = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1ConstantLinearCombinationOperator > >();
   auto coarseGridSolver = std::make_shared< hyteg::CGSolver< hyteg::P1ConstantLinearCombinationOperator > >(
       storage, minLevel, minLevel, max_coarse_iter, coarse_tolerance );
   auto restrictionOperator  = std::make_shared< hyteg::P1toP1LinearRestriction >();
   auto prolongationOperator = std::make_shared< hyteg::P1toP1LinearProlongation >();

   auto gmg = std::make_shared< GeometricMultigridSolver< P1ConstantLinearCombinationOperator > >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

   gmg->setSmoothingSteps( 3, 3 );

   return gmg;
}
/// [CahnHilliardDiagonalPreconditioner create-multigrid-solver]

/// [CahnHilliardFunction ClassCahnHilliardEvolutionOperator]
class CahnHilliardEvolutionOperator
{
 public:
   CahnHilliardEvolutionOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                  uint_t                                     minLevel,
                                  uint_t                                     maxLevel,
                                  real_t                                     tau,
                                  real_t                                     epsilon )
   : epsilon_( epsilon )
   , tau_( tau )
   , lhsOp_( create_lhs_operator( storage, minLevel, maxLevel, tau, epsilon ) )
   , solver_( create_solver( storage, minLevel, maxLevel, tau, epsilon ) )
   , rhsAssembler_( storage, minLevel, maxLevel, epsilon )
   , rhs_( "rhs", storage, minLevel, maxLevel )
   {}

   void apply( const chType& u_old, const chType& u_now, uint_t level );

   void setPrintInfo( bool printInfo ) { solver_->setPrintInfo( printInfo ); }

 protected:
   void scale( const chType& src, const chType& dst, uint_t level, DoFType flag ) const;

 protected:
   real_t epsilon_;
   real_t tau_;

   std::shared_ptr< BlockOperator< chType, chType > > lhsOp_;

   std::shared_ptr< MinResSolver< BlockOperator< chType, chType > > > solver_;

   RHSAssembler rhsAssembler_;

   chType rhs_;
};
/// [CahnHilliardFunction ClassCahnHilliardEvolutionOperator]

/// [CahnHilliardFunction time-evolution-apply]
void CahnHilliardEvolutionOperator::apply( const chType& u_old, const chType& u_now, uint_t level )
{
   lhsOp_->getStorage()->getTimingTree()->start( "apply CH evolution operator" );

   lhsOp_->getStorage()->getTimingTree()->start( "assemble rhs vector" );
   rhsAssembler_.assemble( u_old, rhs_, level, All );
   lhsOp_->getStorage()->getTimingTree()->stop( "assemble rhs vector" );

   scale( rhs_, rhs_, level, All );

   // solve the system
   u_now.assign( { 1 }, { u_old }, level, All );
   solver_->solve( *lhsOp_, u_now, rhs_, level );

   scale( u_now, u_now, level, All );

   lhsOp_->getStorage()->getTimingTree()->stop( "apply CH evolution operator" );
}
/// [CahnHilliardFunction time-evolution-apply]

/// [CahnHilliardFunction scale]
void CahnHilliardEvolutionOperator::scale( const chType& src, const chType& dst, uint_t level, DoFType flag ) const
{
   const real_t factor = std::pow( tau_, 0.25 );
   dst.getMu().assign( { 1. / factor }, { src.getMu() }, level, flag );
   dst.getPhi().assign( { factor }, { src.getPhi() }, level, flag );
}
/// [CahnHilliardFunction scale]

/// [create_storage fun]
std::shared_ptr< hyteg::PrimitiveStorage > create_storage( bool use3D )
{
   uint_t numProcesses = walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   std::shared_ptr< SetupPrimitiveStorage > setupStorage;

   if ( use3D )
   {
      const hyteg::Point3D lowerRight( { -1, -1, -1 } );
      const hyteg::Point3D upperLeft( { +1, +1, +1 } );

      auto meshInfo = hyteg::MeshInfo::meshCuboid( lowerRight, upperLeft, 2, 2, 2 );
      setupStorage  = std::make_shared< hyteg::SetupPrimitiveStorage >( meshInfo, numProcesses );
   }
   else
   {
      /// [create_storage rectangle]
      hyteg::Point2D lowerRight( { -1, -1 } );
      hyteg::Point2D upperLeft( { +1, +1 } );

      auto meshInfo = hyteg::MeshInfo::meshRectangle( lowerRight, upperLeft, hyteg::MeshInfo::CROSS, 2, 2 );
      setupStorage  = std::make_shared< hyteg::SetupPrimitiveStorage >( meshInfo, numProcesses );
      /// [create_storage rectangle]
   }

   /// [create_storage neumann]
   // Neumann BC everywhere
   setupStorage->setMeshBoundaryFlagsOnBoundary( 2, 0, true );
   /// [create_storage neumann]

   return std::make_shared< hyteg::PrimitiveStorage >( *setupStorage );
}
/// [create_storage fun]

/// [main config-walberla]
int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();
   /// [main config-walberla]

   /// [main config-file]
   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "./CahnHilliard.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );
   parameters.listParameters();

   const uint_t minLevel       = parameters.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel       = parameters.getParameter< uint_t >( "maxLevel" );
   const uint_t outputInterval = parameters.getParameter< uint_t >( "outputInterval" );

   const real_t tau     = parameters.getParameter< real_t >( "tau" );
   const real_t epsilon = parameters.getParameter< real_t >( "epsilon" );

   const real_t t_end = parameters.getParameter< real_t >( "t_end" );
   /// [main config-file]

   /// [main storage]
   const auto storage = create_storage( parameters.getParameter< bool >( "use3D" ) );
   /// [main storage]

   /// [main functions]
   chType u_now( "u_now", storage, minLevel, maxLevel );
   chType u_prev( "u_prev", storage, minLevel, maxLevel );
   /// [main functions]

   /// [main evolution-operator]
   CahnHilliardEvolutionOperator chOperator( storage, minLevel, maxLevel, tau, epsilon );

   if ( parameters.getParameter< bool >( "printMinresConvergence" ) )
      chOperator.setPrintInfo( true );
   /// [main evolution-operator]

   /// [main initial-values]
   walberla::math::seedRandomGenerator( static_cast< std::mt19937::result_type >( walberla::mpi::MPIManager::instance()->rank() ) );

   u_prev.getPhi().interpolate(
       []( const Point3D& ) { return walberla::math::realRandom( real_t( -0.01 ), real_t( 0.01 ) ); }, maxLevel, All );
   u_prev.getMu().interpolate( 0, maxLevel, All );
   /// [main initial-values]

   /// [main time-loop]
   real_t t_now = 0;

   VTKOutput vtkOutput( ".", "CahnHilliard", storage );

   // save initial state:
   u_now.assign( { 1. }, { u_prev }, maxLevel, All );
   vtkOutput.add( u_now );
   vtkOutput.write( maxLevel, 0 );

   for ( uint_t k = 1; k < std::numeric_limits< uint_t >::max(); ++k )
   {
      t_now += tau;

      if ( t_now > t_end )
         break;

      WALBERLA_LOG_INFO_ON_ROOT( "iter = " << k << " t = " << t_now );
      chOperator.apply( u_prev, u_now, maxLevel );

      u_prev.assign( { 1. }, { u_now }, maxLevel, All );

      if ( k % outputInterval == 0 )
      {
         vtkOutput.add( u_now );
         vtkOutput.write( maxLevel, k );
      }
   }
   /// [main time-loop]

   /// [main timing-tree]
   if ( parameters.getParameter< bool >( "printTimingTree" ) )
      WALBERLA_LOG_INFO_ON_ROOT( *storage->getTimingTree() );
}
/// [main timing-tree]
