
#pragma once

#include "core/DataTypes.h"

namespace terraneo {

using walberla::real_t;
using walberla::uint_t;

//! Class providing methods for evaluation of Spherical Harmonics

//! This class provides routines for the evaluation and calculation of
//! spherical harmonics. It is a C++ adaptation of old TERRA routines.
class SphericalHarmonicsTool
{
 private:
   struct
   {
      real_t *f1, *f2, *fac1, *fac2, *srt;
   } plm0_;

   //! Remember, if we have initialized plmbar
   bool initialized_;

   //! Remember largest degree used in intialisation
   uint_t maxDegree_;

   //! Dynamic array for evaluating associated Legendre functions
   // real_t* plm_;
   real_t* plm_;

   //! Dynamic array for evaluating associated Legendre functions derivatives
   real_t* dplm_;

   //! Dynamic array for evaluating sine/cosine parts of \f$Y_l^m\f$
   real_t** csm_;

   /// Compute derivatives of associated Legendre functions

   /// The method evaluates the derivative \f$ dP_l^m(\cos(\theta)) / d\theta \f$
   /// for all associated Legendre functions for given \f$\theta\f$ up to a
   /// specified degree. The computation of the derivatives is based on the
   /// identity
   ///
   /// \f[
   ///    2\frac{d P_l^m}{d \theta}(\cos(\theta)) =
   ///    (l+m)(l-m+1)P_l^{m-1}(\cos(\theta)) - P_l^{m+1}(\cos(\theta))
   /// \f]
   ///
   /// which holds for orders (m>0) and the un-normalised associated Legendre
   /// functions. For details see <em>W. Bosch, On the Computation of
   /// Derivatives of Legendre Functions, Phys. Chem. Earth (A), vol. 25,
   /// no. 9-11, pp. 655-659, 2000</em>. Note that (A8) in that paper is missing
   /// a \f$\sqrt{2}\f$ term.
   ///
   /// Implemented is a slightly modified version which holds for the associated
   /// Legendre functions normalised in the Geodesy fashion. This version is
   /// also descibed in Bosch (2000).
   ///
   /// \param   dplm    array for storing computed derivative values
   /// \param   plm     array containing values of \f$P_l^m\cos(\theta)\f$,
   ///                  implicitely defines co-latitude
   /// \param   lmax    maximal degree of spherical harmonics
   void dplmbar( real_t* dplm, real_t* plm, uint_t lmax );

   /// Determine position of value of Legendre function in internal array
   uint_t getArrayIndex( uint_t deg, uint_t ord )
   {
      int idx = 0;
      for ( uint_t i = 0; i <= deg; ++i )
      {
         for ( uint_t j = 0; j <= i; ++j )
         {
            idx++;
            if ( i == deg && j == ord )
               goto loopEnd;
         }
      }
   loopEnd:
      // correct for 0-based indexing
      idx--;
      return idx;
   }

 public:
   SphericalHarmonicsTool( uint_t lmax );
   ~SphericalHarmonicsTool();

   //! Evaluate scalar function described by spherical harmonics coefficients
   //! at given cartesian coordinates. The point has to be on the sphere.
   //!
   //!  \param  sphdata  Storing the data in sph format. Indices as follows:
   //!                   sphdata[kk][ir][ll][mm][cs]
   //!                      kk: For vector data, kk=0..nk-1
   //!                      ir: radial layer, counting outwards, 0..nr-1.
   //!                      ll: Spherical harmonic coefficients, 0..lmax
   //!                      mm: Spherical harmonic coefficients, 0..ll
   //!                      cs: We need two coefficients, for the sine and for
   //!                          the cosine, 0..1
   //!  \param  x      x-coordinate of evaluation node
   //!  \param  y      y-coordinate of evaluation node
   //!  \param  z      z-coordinate of evaluation node
   //!  \param  ir     radial index in array sphdata
   //!  \param  lmax   largest spherical harmonics degree
   //!  \param  kk     index kk in array sphdata
   real_t shconvert_eval( real_t***** sphdata, real_t x, real_t y, real_t z, uint_t ir, uint_t lmax, uint_t kk );

   //! Evaluate a single spherical harmonics function of given degree and order

   //! Evaluate a single spherical harmonics function of specified degree and
   //! order (l,m) at given cartesian coordinates (x,y,z). For positive argument
   //! ord we use the sine variant, otherwise the cosine variant.
   //!
   //! \param  deg    degree l
   //! \param  ord    order  m (also decides on sin/cos)
   //! \param  x      x-coordinate of evaluation node
   //! \param  y      y-coordinate of evaluation node
   //! \param  z      z-coordinate of evaluation node
   real_t shconvert_eval( uint_t deg, int ord, real_t x, real_t y, real_t z );

   //! Evaluate scalar function described by spherical harmonics coefficients
   //! at given cartesian coordinates. In this variant we perform linear
   //! interpolation in radial direction between the different radial layers,
   //! for which the spherical harmonics expansion was computed.
   //!
   //!  \param  sphdata  Storing the data in sph format. Indices as follows:
   //!                   sphdata[kk][ir][ll][mm][cs]
   //!                      kk: For vector data, kk=0..nk-1
   //!                      ir: radial layer, counting outwards, 0..nr.
   //!                      ll: Spherical harmonic coefficients, 0..lmax
   //!                      mm: Spherical harmonic coefficients, 0..ll
   //!                      cs: We need two coefficients, for the sine and for
   //!                          the cosine, 0..1
   //!  \param  x      x-coordinate of evaluation node
   //!  \param  y      y-coordinate of evaluation node
   //!  \param  z      z-coordinate of evaluation node
   //!  \param  rmin   radius of innermost shell of expansion
   //!  \param  rmax   radius of outermost shell of expansion
   //!  \param  nr     radial index in outermost shell in array sphdata
   //!  \param  lmax   largest spherical harmonics degree
   //!  \param  kk     index kk in array sphdata
   real_t shconvert_eval( real_t***** sphdata,
                          real_t      x,
                          real_t      y,
                          real_t      z,
                          real_t      rmin,
                          real_t      rmax,
                          uint_t      nr,
                          uint_t      lmax,
                          uint_t      kk );

   /// Evaluate associated Legendre functions
   ///
   /// The subroutine evaluates all associated Legendre functions \f$ P_l^m \f$
   /// of the 1st kind for a given argument up to a given degree l.
   ///
   /// #### Normalisation ####
   /// Normalisation of the associated Legendre functions is done in the
   /// following way:
   /// - We do not include the Condon-Shortley phase \f$ (-1)^m \f$.
   /// - We normalise such that the associated real-valued spherical harmonics
   ///   functions (trig being either sin or cos)
   /// \f[
   ///     Y_l^m(\theta,\phi) = P_l^m(\cos(\theta)) \,\mbox{trig}\,(m \phi)
   /// \f]
   ///   satisfy
   /// \f[
   ///     \int_S Y_l^m(\theta,\phi)\cdot Y_l^m(\theta,\phi) = 4 \pi
   /// \f]
   ///   Here \f$(\theta, \phi)\f$ represent colatitude and longitude. The
   ///   normalisation is the one preferred in Geodesy. The normalisation factor
   ///   in detail is given by
   /// \f[
   ///     \sqrt{2-\delta_{0,m}}\sqrt{(2l+1)\frac{(l-m)!}{(l+m)!}}
   /// \f]
   ///
   /// #### Algorithm ####
   /// The computation of the associated Legendre functions employs the so called
   /// <b>Forward Column Method</b>, see e.g. *Holmes & Featherstone, 2002, J. of
   /// Geodesy*. This is based on the three-term recurrence relation
   ///
   /// \f[
   ///    (l-m)P_l^m(z) = z(2l-1)P_{l-1}^m(z) - (l+m-1)P_{l-2}^m(z)
   /// \f]
   ///
   /// After computing \f$P_m^m\f$ and \f$P_{m+1}^m\f$ to initiate the recurrence
   /// we compute \f$P_l^m\f$ by increasing l while keeping m fixed.
   ///
   /// #### Storage ####
   /// - Dimension in the calling part of the program must be
   ///   p( (lmax+1)*(lmax+2)/2 ) and the same for array dp.
   /// - p(k) contains p(l,m) with k=(l+1)*l/2+m+1; i.e. m increments through
   ///   range 0 to l before incrementing l.
   ///
   /// #### Original comments in Terra ####
   /// > Routine is stable in single and real_t precision to l,m = 511 at least;
   /// > timing is proportional to lmax**2 R.J.O'Connell 7 Sept. 1989;
   /// > added dp(z) 10 Jan. 1990
   /// > Added precalculation and storage of square roots srl(k) 31 Dec 1992
   ///
   /// \param    p        array for storing values \f$ P_l^m(z) \f$
   /// \param    lmax     we evaluate all functions up to degree lmax
   /// \param    z        argument to \f$ P_l^m \f$, i.e. z = cos(colatitude)
   void plmbar( real_t* p, uint_t lmax, real_t z );

   //! Compute derivatives of associated Legendre functions

   //! This public method first evaluates the associated Legendre functions up
   //! to the specified degree for the given co-latitude and then their
   //! derivatives with respect to co-latitude. The actual computations are
   //! performed by the private internal methods
   //!
   //! \param  dplm   array to store derivative values
   //! \param  lmax   maximal degree up to which to compute
   //! \param  z      \f$\cos(\theta)\f$, i.e. cosine of co-latitude
   void dplmbar( real_t* dplm, uint_t lmax, real_t z );

   /// Method for evaluating vector spherical harmonics

   /** We use vector spherical harmonics given by
    *  \f[ \begin{split}
    *    \mathcal{Y}_{lm}^0(\theta,\phi) &= \nu_0 \vec{e}_r \\
    *   \mathcal{Y}_{lm}^1(\theta,\phi) &= \nu_1 \vec{e}_\theta +
    *                                       \nu_2 \vec{e}_\phi \\
    *   \mathcal{Y}_{lm}^2(\theta,\phi) &= - \nu_2 \vec{e}_\theta
    *                                       + \nu_1 \vec{e}_\phi
    *  \end{split} \f]
    *  where \f$\{\vec{e}_r,\vec{e}_\theta,\vec{e}_\phi\}\f$ is the local
    *  triad at the point \f$(\theta,\phi)\f$ on the unit sphere and the
    *  coefficient functions are given by
    *  \f[ \begin{split}
    *  \nu_0(\theta,\phi) &= Y_{lm}(\theta,\phi) \\
    *  \nu_1(\theta,\phi) &= \frac{1}{\sqrt{l(l+1)}}\cdot
    *           \frac{\partial}{\partial \theta} Y_{lm}(\theta,\phi) \\
    *  \nu_2(\theta,\phi) &= \frac{1}{\sqrt{l(l+1)}}\cdot
    *            \frac{1}{\sin\theta}\frac{\partial}{\partial \phi}
    *            Y_{lm}(\theta,\phi)
    *  \end{split} \f]
    *  Here \f$Y_{lm}\f$ is the scalar spherical harmonics function.
    *
    *  \param  deg    degree l
    *  \param  ord    order  m
    *  \param  x      x-coordinate of evaluation node
    *  \param  y      y-coordinate of evaluation node
    *  \param  z      z-coordinate of evaluation node
    *  \param  ind    index for choosing vector from \f$\{\mathcal{Y}_{lm}^0,
    *                 \mathcal{Y}_{lm}^1,\mathcal{Y}_{lm}^2\}\f$
    *  \param  comp   choosing 0th, 1st or 2nd vector component
    *
    *  \note
    *  - For positive order we use scalar spherical harmonics with a
    *    \f$\sin(m\theta)\f$ factor, otherwise with \f$\cos(m\theta)\f$.
    *  - If \f$p=(x,y,z)\f$ is a point on the unit sphere, then we can
    *    express the vectors of the local triad by
    *  - As all methods of this class, also this one employs Geodesy style
    *    normalisation.
    *  \f[
    *  \vec{e}_r = \begin{pmatrix} x \\ y \\ z \end{pmatrix}
    *  \enspace,\quad
    *  \vec{e}_\theta = \begin{pmatrix} zx/\varrho \\ zy/\varrho \\
    *  -\varrho \end{pmatrix}
    *  \enspace,\quad
    *  \vec{e}_\phi = \begin{pmatrix} -y/\varrho \\ x/\varrho \\ 0 \end{pmatrix}
    *  \enspace,\quad
    *  \varrho = \sqrt{x^2 + y^2}
    *  \f]
   **/
   real_t evalVSH( uint_t deg, int ord, real_t x, real_t y, real_t z, uint_t ind, uint_t comp );
};

} // namespace terraneo
