
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

#include "core/debug/Debug.h"
#include "core/logging/Logging.h"

namespace terraneo {

using walberla::real_c;
using walberla::uint_t;

// Auxiliary routines for spherical harmonics calculations

SphericalHarmonicsTool::SphericalHarmonicsTool( uint_t lmax )
{
   initialized_ = false;
   maxDegree_   = lmax;

   plm_    = new real_t[( lmax + 1 ) * ( lmax + 2 ) / 2];
   dplm_   = new real_t[( lmax + 1 ) * ( lmax + 2 ) / 2];
   csm_    = new real_t*[2];
   csm_[0] = new real_t[lmax + 1];
   csm_[1] = new real_t[lmax + 1];

   plm0_.f1   = new real_t[( lmax + 1 ) * ( lmax + 2 ) / 2 - 1];
   plm0_.f2   = new real_t[( lmax + 1 ) * ( lmax + 2 ) / 2 - 1];
   plm0_.fac1 = new real_t[( lmax + 1 ) * ( lmax + 2 ) / 2 - 1];
   plm0_.fac2 = new real_t[( lmax + 1 ) * ( lmax + 2 ) / 2 - 1];
   plm0_.srt  = new real_t[2 * lmax + 2];
}

SphericalHarmonicsTool::~SphericalHarmonicsTool()
{
   delete[] plm_;
   delete[] dplm_;
   delete[] csm_[0];
   delete[] csm_[1];
   delete[] csm_;
   delete[] plm0_.f1;
   delete[] plm0_.f2;
   delete[] plm0_.fac1;
   delete[] plm0_.fac2;
   delete[] plm0_.srt;
}

// ================
//  shconvert_eval
// ================
//
// Evaluate scalar function described by spherical harmonics coefficients
// at given cartesian coordinates. The point has to be on the sphere.
real_t
    SphericalHarmonicsTool::shconvert_eval( real_t***** sphdata, real_t x, real_t y, real_t z, uint_t ir, uint_t lmax, uint_t kk )
{
   // local variables
   real_t phi, r, val;
   uint_t a;

   // compute radius of point
   r = sqrt( x * x + y * y + z * z );

   // compute longitude phi. Poles are no problem, since c++-atan2 implements
   // atan2(0,0)
   phi = atan2( y, x );

   // Calculate associated Legendre functions
   plmbar( plm_, lmax, z / r );

   for ( uint_t mm = 0; mm <= lmax; mm++ )
   {
      csm_[0][mm] = cos( real_c( mm ) * phi );
      csm_[1][mm] = sin( real_c( mm ) * phi );
   }

   // Assemble the contributions
   val = 0.;
   a   = 0;
   for ( uint_t ll = 0; ll <= lmax; ll++ )
   {
      for ( uint_t mm = 0; mm <= ll; mm++ )
      {
         val = val + plm_[a] * ( csm_[0][mm] * sphdata[kk][ir][ll][mm][0] + csm_[1][mm] * sphdata[kk][ir][ll][mm][1] );
         a++;
      }
   }

   return val;
}

// ================
//  shconvert_eval
// ================
//
// Evaluate a single spherical harmonics function of specified degree and order
// (l,m) at given cartesian coordinates (x,y,z). For positive argument ord
// we use the sine variant otherwise the cosine variant.
real_t SphericalHarmonicsTool::shconvert_eval( uint_t deg, int ord, real_t x, real_t y, real_t z )
{
   // process and check input
   uint_t order = (uint_t) std::abs( ord );
   WALBERLA_ASSERT( deg <= maxDegree_ );
   WALBERLA_ASSERT( order <= deg );

   // compute radius of point
   real_t rad = sqrt( x * x + y * y + z * z );

   // compute longitude phi. Poles are no problem, since c++-atan2 implements
   // atan2(0,0)
   real_t phi = atan2( y, x );

   // evaluate associated Legendre functions for theta
   plmbar( plm_, maxDegree_, z / rad );

   // determine position of value of Legendre functions in array
   uint_t idx = getArrayIndex( deg, order );

   // evaluate our spherical harmonics function
   real_t retVal = 0.0;
   if ( ord <= 0 )
   {
      retVal = plm_[idx] * cos( real_c( order ) * phi );
   }
   else
   {
      retVal = plm_[idx] * sin( real_c( order ) * phi );
   }

   return retVal;
}

// ================
//  shconvert_eval
// ================
//
// Evaluate scalar function described by spherical harmonics coefficients
// at given cartesian coordinates. In this variant we perform linear
// interpolation in radial direction between the different radial layers,
// for which the spherical harmonics expansion was computed.
real_t SphericalHarmonicsTool::shconvert_eval( real_t***** sphdata,
                                               real_t      x,
                                               real_t      y,
                                               real_t      z,
                                               real_t      rmin,
                                               real_t      rmax,
                                               uint_t      nr,
                                               uint_t      lmax,
                                               uint_t      kk )
{
   // ----------------------------
   //  compute radial information
   // ----------------------------

   // radius
   real_t r = sqrt( x * x + y * y + z * z );

   // push everything to the boundary that is out of range
   if ( r < rmin )
   {
      r = rmin;
   }
   else if ( r > rmax )
   {
      r = rmax;
   }

   // layer thickness
   real_t h = ( rmax - rmin ) / (double) ( nr - 1 );

   // radial layer index
   uint_t ir = uint_t( ( r - rmin ) / h );

   // -----------------------
   //  prepare interpolation
   // -----------------------

   // evaluate associated Legendre functions via class method
   plmbar( plm_, lmax, z / r );

   // set interpolation weights
   real_t wLow = ( real_t )( ir + 1 ) + ( rmin - r ) / h;
   real_t wUpp = ( r - rmin ) / h - real_t( ir );

   // compute longitude phi; poles are no problem, since c++-atan2
   // implements atan2(0,0)
   real_t phi = atan2( y, x );

   // pre-compute sine and cosine factors
   for ( uint_t mm = 0; mm <= lmax; ++mm )
   {
      csm_[0][mm] = cos( real_c( mm ) * phi );
      csm_[1][mm] = sin( real_c( mm ) * phi );
   }

   // -------------------------------------
   //  perform summation and interpolation
   // -------------------------------------

   // assemble the contributions
   real_t vLow = 0.0;
   real_t vUpp = 0.0;
   uint_t idx  = 0;
   for ( uint_t ll = 0; ll <= lmax; ++ll )
   {
      for ( uint_t mm = 0; mm <= ll; ++mm )
      {
         vLow += plm_[idx] * ( csm_[0][mm] * sphdata[kk][ir][ll][mm][0] + csm_[1][mm] * sphdata[kk][ir][ll][mm][1] );

         vUpp += plm_[idx] * ( csm_[0][mm] * sphdata[kk][ir + 1][ll][mm][0] + csm_[1][mm] * sphdata[kk][ir + 1][ll][mm][1] );
         idx++;
      }
   }

#ifndef NDEBUG
   real_t value = vLow * wLow + vUpp * wUpp;
   if ( std::isnan( value ) || std::isinf( value ) )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Warning: shconvert_eval detect NaN or Inf. Replacing value by zero." );
      return 0.0;
   }
#endif

   return vLow * wLow + vUpp * wUpp;
}

// ========
//  plmbar
// ========
void SphericalHarmonicsTool::plmbar( real_t* p, uint_t lmax, real_t z )
{
   // local variables
   real_t fden, fnum, pm1, pm2, pmm, sintsq, plm;
   uint_t k, kstart;

   if ( std::abs( z ) > 1.0 || lmax > maxDegree_ )
   {
      WALBERLA_ABORT( "Parameter inconsistency in SphericalHarmonicsTool::plmbar!" );
   }

   // ---------------------------------------
   //  set up sqrt and factors on first pass
   // ---------------------------------------
   if ( !initialized_ )
   {
      for ( k = 1; k <= 2 * lmax + 2; k++ )
      {
         plm0_.srt[k - 1] = sqrt( real_c( k ) );
      }

      if ( lmax == 0 )
      {
         p[0] = 1.0;
         return;
      }

      // case for m > 0
      kstart = 1;

      for ( uint_t m = 1; m <= lmax; m++ )
      {
         // case for P(m,m)
         kstart = kstart + m + 1;

         if ( m != lmax )
         {
            // case for P(m+1,m)
            k = kstart + m + 1;

            // case for P(l,m) with l > m+1
            if ( m < lmax - 1 )
            {
               for ( uint_t l = m + 2; l <= lmax; l++ )
               {
                  k               = k + l;
                  plm0_.f1[k - 1] = ( plm0_.srt[2 * l] * plm0_.srt[2 * l - 2] ) / ( plm0_.srt[l + m - 1] * plm0_.srt[l - m - 1] );
                  plm0_.f2[k - 1] = ( plm0_.srt[2 * l] * plm0_.srt[l - m - 2] * plm0_.srt[l + m - 2] ) /
                                    ( plm0_.srt[2 * l - 4] * plm0_.srt[l + m - 1] * plm0_.srt[l - m - 1] );
               }
            }
         }
      }

      k = 3;

      for ( uint_t l = 2; l <= lmax; l++ )
      {
         k = k + 1;
         for ( uint_t m = 1; m <= l - 1; m++ )
         {
            k                 = k + 1;
            plm0_.fac1[k - 1] = plm0_.srt[l - m - 1] * plm0_.srt[l + m];
            plm0_.fac2[k - 1] = plm0_.srt[l + m - 1] * plm0_.srt[l - m];
            if ( m == 1 )
            {
               plm0_.fac2[k - 1] = plm0_.fac2[k - 1] * plm0_.srt[1];
            }
         }
         k = k + 1;
      }

      initialized_ = true;
   }

   // --------------------------------
   //  start calculation of Plm, etc.
   // --------------------------------

   // case for P(l,0)
   pm2  = 1.;
   p[0] = 1.;

   if ( lmax == 0 )
      return;

   pm1  = z;
   p[1] = plm0_.srt[2] * pm1;
   k    = 2;

   for ( uint_t l = 2; l <= lmax; l++ )
   {
      k        = k + l;
      plm      = ( real_c( 2 * l - 1 ) * z * pm1 - real_c( l - 1 ) * pm2 ) / real_c( l );
      p[k - 1] = plm0_.srt[2 * l] * plm;
      pm2      = pm1;
      pm1      = plm;
   }

   // case for m > 0
   pmm    = 1.;
   sintsq = ( 1. - z ) * ( 1. + z );
   fnum   = -1.;
   fden   = 0.;
   kstart = 1;

   for ( uint_t m = 1; m <= lmax; m++ )
   {
      // case for P(m,m)

      kstart        = kstart + m + 1;
      fnum          = fnum + 2.;
      fden          = fden + 2.0;
      pmm           = pmm * sintsq * fnum / fden;
      pm2           = sqrt( real_c( 4 * m + 2 ) * pmm );
      p[kstart - 1] = pm2;

      if ( m != lmax )
      {
         // case for P(m+1,m)
         pm1      = z * plm0_.srt[2 * m + 2] * pm2;
         k        = kstart + m + 1;
         p[k - 1] = pm1;
         // case for P(l,m) with l > m+1
         if ( m < lmax - 1 )
         {
            for ( uint_t l = m + 2; l <= lmax; l++ )
            {
               k        = k + l;
               plm      = z * plm0_.f1[k - 1] * pm1 - plm0_.f2[k - 1] * pm2;
               p[k - 1] = plm;
               pm2      = pm1;
               pm1      = plm;
            }
         }
      }
   }
}

// =========
//  dplmbar
// =========
void SphericalHarmonicsTool::dplmbar( real_t* dplm, real_t* plm, uint_t lmax )
{
   if ( lmax > maxDegree_ )
   {
      WALBERLA_ABORT( "Parameter inconsistency in SphericalHarmonicsTool::dplmbar!" );
   }

   // Holy obsfuscation! We fetch address of srt and shift it back so that
   // srt[idx] gives sqrt(idx)
   real_t* srt = ( plm0_.srt ) - 1;

   dplm[1] = -plm[2];
   dplm[2] = plm[1];

   uint_t k = 2;
   uint_t m;

   for ( uint_t l = 2; l <= lmax; ++l )
   {
      k++;

      // treat m=0 and m=l separately (note the indexing offset needed for srt!)
      dplm[k]     = -srt[l] * srt[l + 1] / srt[2] * plm[k + 1];
      dplm[k + l] = srt[l] / srt[2] * plm[k + l - 1];

      // if( m == 1) fac2 = fac2*srt(2)
      for ( m = 1; m <= l - 1; ++m )
      {
         k++;
         dplm[k] = 0.5 * ( plm0_.fac2[k] * plm[k - 1] - plm0_.fac1[k] * plm[k + 1] );
      }
      k++;
   }
}

// =========
//  dplmbar
// =========
//
// Public method with different interface than private dplmbar()
void SphericalHarmonicsTool::dplmbar( real_t* dplm, uint_t lmax, real_t z )
{
   WALBERLA_ASSERT( lmax > 0 );
   WALBERLA_ASSERT( lmax <= maxDegree_ );

   // run plmbar to get ALFs needed for computing derivatives
   plmbar( plm_, lmax, z );

   // now compute derivates
   dplmbar( dplm, plm_, lmax );
}

// =========
//  evalVSH
// =========
//
// Method for evaluating vector spherical harmonics
real_t SphericalHarmonicsTool::evalVSH( uint_t deg, int ord, real_t x, real_t y, real_t z, uint_t ind, uint_t comp )
{
   WALBERLA_ASSERT( ind <= 2 );
   WALBERLA_ASSERT( comp <= 2 );
   WALBERLA_ASSERT( deg <= maxDegree_ );
   WALBERLA_ASSERT( walberla::uint_c( std::abs( ord ) ) <= deg );

   // compute radius of point
   real_t rad = sqrt( x * x + y * y + z * z );

   // map point onto unit sphere
   real_t xSph = x / rad;
   real_t ySph = y / rad;
   real_t zSph = z / rad;

   // compute longitude phi (poles are no problem, since c++-atan2 implements
   // atan2(0,0))
   real_t phi = atan2( y, x );

   // compute co-latitude theta
   real_t theta = acos( zSph );

   // compute origin distance in x-y-plane
   // real_t rho = sqrt( xSph*xSph + ySph*ySph );

   // scaling factors
   real_t scale    = 1.0 / sqrt( ( real_t )( deg * ( deg + 1 ) ) );
   real_t sinTheta = sqrt( 1.0 - zSph * zSph );
   // real_t sinTheta = sin( theta );

   real_t nu1    = 0.0;
   real_t nu2    = 0.0;
   real_t retVal = 0.0;

   // ------------------------------
   //  Poles need special treatment
   // ------------------------------

   // For the two tangential vectors, if the absolute order is different from 1,
   // then these become zero at the poles. The only problematic case is the one
   // where we consider the tangential vectors and the order = +-1.
   //
   // The radial vector is basically only the value of Y_l^m at the pole which
   // should be taken care of in the plmbar computation already.
   // if ( ind != 0 && ( zSph == 1.0 || zSph == -1.0 ) )
   if ( ind != 0 && sinTheta <= 1e-14 )
   {
      if ( ord != 1 && ord != -1 )
      {
         return 0.0;
      }
      else
      {
         static bool haveWarned = false;
         if ( haveWarned == false )
         {
            std::cerr << " WARNING: Using tricks at pole for m = 1!\n";
            haveWarned = true;
         }

         // z-component will always be zero
         // x-component will be zero for vector index 1
         // y-component will be zero for vector index 2
         if ( comp == 2 || ( ind == 1 && comp == 0 ) || ( ind == 2 && comp == 1 ) )
         {
            return 0.0;
         }
         real_t scalFac = sqrt( 7.0 / 6.0 );

         if ( ind == 1 )
         {
            return sqrt( (double) deg ) * scalFac;
         }
         else
         {
            return -sqrt( (double) deg ) * zSph * scalFac;
         }
      }
   }

   // ---------------
   //  Non-Pole case
   // ---------------
   // determine real-valued trigonometric components, we use sine for positive
   // order and cosine otherwise; set neutral order value
   real_t trig  = 0.0;
   real_t dtrig = 0.0;
   real_t order = (real_t) abs( ord );
   if ( ord > 0 )
   {
      trig  = sin( order * phi );
      dtrig = order * cos( order * phi );
   }
   else if ( ord < 0 )
   {
      trig  = cos( order * phi );
      dtrig = -order * sin( order * phi );
   }
   else
   {
      trig  = 1.0;
      dtrig = 0.0;
   }

   // determine component of local triad vectors
   real_t e_r_comp     = 0.0;
   real_t e_theta_comp = 0.0;
   real_t e_phi_comp   = 0.0;

   switch ( comp )
   {
      // case 0:
      //   e_r_comp     = xSph;
      //   e_theta_comp = zSph * xSph / rho;
      //   e_phi_comp   = -ySph / rho;
      //   break;

      // case 1:
      //   e_r_comp     = ySph;
      //   e_theta_comp = zSph * ySph / rho;
      //   e_phi_comp   = xSph / rho;
      //   break;

      // case 2:
      //   e_r_comp     = zSph;
      //   e_theta_comp = - rho;
      //   e_phi_comp   = 0.0;
      //  break;

   case 0:
      e_r_comp     = xSph;
      e_theta_comp = cos( theta ) * cos( phi );
      e_phi_comp   = -sin( phi );
      break;

   case 1:
      e_r_comp     = ySph;
      e_theta_comp = cos( theta ) * sin( phi );
      e_phi_comp   = cos( phi );
      break;

   case 2:
      e_r_comp     = zSph;
      e_theta_comp = -sin( theta );
      e_phi_comp   = 0.0;
      break;
   }

   if ( std::isnan( e_theta_comp ) || std::isnan( e_phi_comp ) )
   {
      WALBERLA_ABORT( "Problem with components!" );
   }

   // need to distinguish which vector is needed
   uint_t idx = 0;
   switch ( ind )
   {
   case 0:
      // ------------------------------
      //  evaluation of radial vector
      // ------------------------------
      retVal = this->shconvert_eval( deg, ord, x, y, z ) * e_r_comp;
      break;

   case 1:
      // -------------------------------------
      //  evaluation of 1st tangential vector
      // -------------------------------------

      // compute ALF derivatives for zSph = cos(theta)
      dplmbar( dplm_, maxDegree_, zSph );
      idx    = getArrayIndex( deg, (uint_t) std::abs( ord ) );
      nu1    = dplm_[idx] * trig * scale;
      nu2    = scale / sinTheta * dtrig * plm_[idx];
      retVal = nu1 * e_theta_comp + nu2 * e_phi_comp;
      break;

   case 2:
      // -------------------------------------
      //  evaluation of 2nd tangential vector
      // -------------------------------------
      dplmbar( dplm_, maxDegree_, zSph );
      idx    = getArrayIndex( deg, (uint_t) std::abs( ord ) );
      nu1    = dplm_[idx] * trig * scale;
      nu2    = scale / sinTheta * dtrig * plm_[idx];
      retVal = -nu2 * e_theta_comp + nu1 * e_phi_comp;

      break;
   }

   // ensure we did not run into "pole" problems
   WALBERLA_ASSERT( !std::isnan( retVal ) );
   return retVal;
}

} // namespace terraneo
