
#pragma once

#include <sstream>
#include <fstream>
#include <iostream>
#include <string.h>
#include <math.h>

#include "core/DataTypes.h"

#include "tinyhhg_core/types/pointnd.hpp"
#include "tinyhhg_core/geophysics/SphericalHarmonics.hpp"

namespace hhg {

using walberla::uint_t;

template<typename TYPENAME>
class TomoVolumeFunction
{

public:

    //! Constructor
    TomoVolumeFunction( uint_t lvlMax, TYPENAME sx, TYPENAME sy, TYPENAME sz,
                        double setrmin, double setrmax, double scaling,
                        std::string tomofile );

    //! Destructor
    ~TomoVolumeFunction();

    TYPENAME operator() (TYPENAME x, TYPENAME y, TYPENAME z);
    TYPENAME operator() ( const Point3D & p );

    void read_line( std::ifstream &infile );

private:

    const TYPENAME sx_, sy_, sz_;
    double rho(double x, double y, double z);

    // Storing the data in sph format. Indices as follows:
    // Access: sphdata[kk][ir][ll][mm][cs]
    // kk: For vector data, kk=0..nk-1.
    // ir: radial layer, counting outwards, 0..nr.
    // ll: Spherical harmonic coefficients, 0..lmax
    // mm: Spherical harmonic coefficients, 0..ll
    // cs: We need two coefficients, for the sine and for the cosine, 0..1
    double *sphdata_buf;
    double *****sphdata;

    // Converter
    SphericalHarmonicsTool *shutil;

    // Resolution options for the tomography
    // mt: Number of nodes along a diamond edge
    // nk: Dimension of a data point. 1 for scalars, 3 for vectors, ...
    // nr: Number of radial layers - 1
    // nd: number of diamonds
    int lmax, nk, nr;

    // Radii of the spherical shell
    double rmin, rmax;

    // Scaling factor for the RHS.
    double scaling_factor;

    // Number of radial layers in the hhg grid
    int nr_hhg;

    // Structure to store a line of a file and convert it into a stream.
    // Declared here to avoid calling a constructor and destructor all the
    // time.
    std::stringstream instream;
    std::string line;

    // The tomography model to open
    std::string filename;

}; // class tnTomoVolumeFunction


// Read a line from a file, ignoring comments
template<typename TYPENAME>
void TomoVolumeFunction<TYPENAME>::read_line( std::ifstream &infile ) {
  line = "#";
  while((line.compare(0,1,"#") == 0)) {
    if(!getline(infile, line)) {
      break;
    }
  }
  instream.str("");
  instream.clear();
  instream << line;
}


// ==========================================================================
// Implementation of methods
// ==========================================================================

// Constructor
template<typename TYPENAME>
TomoVolumeFunction<TYPENAME>::TomoVolumeFunction( uint_t lvlMax,
                                                  TYPENAME sx,
                                                  TYPENAME sy,
                                                  TYPENAME sz,
                                                  double setrmin,
                                                  double setrmax,
                                                  double scaling,
                                                  std::string tomofile )
: sx_(sx), sy_(sy), sz_(sz) {

  filename = tomofile;
  scaling_factor = scaling;

  rmin = setrmin;
  rmax = setrmax;
  nr_hhg = pow(2.0, static_cast<int>(lvlMax)) + 1;

  // Read the names of the entries from the input file
  std::ifstream infile;
  infile.open(filename.c_str());
  if (!infile.is_open())  {
    std::cout << "Error opening tomography file" << std::endl;
    exit(-1);
  }
  read_line(infile);
  instream >> lmax;
  instream >> nr;
  instream >> nk;

  // Skip the names of the columns
  read_line(infile);

  // Allocate the memory for the tomography data
  sphdata_buf = new double[nk*(nr+1)*(lmax+1)*(lmax+2)];
  sphdata = new double****[nk];
  int index = 0;
  for(int kk=0; kk<nk; kk++) {
    sphdata[kk] = new double***[nr+1];
    for(int ir=0; ir<=nr; ir++) {
      sphdata[kk][ir] = new double**[lmax+1];
      for(int ll=0; ll<=lmax; ll++) {
        sphdata[kk][ir][ll] = new double*[ll+1];
        for(int mm=0; mm<=ll; mm++) {
          sphdata[kk][ir][ll][mm] = &sphdata_buf[index];
          index+=2;
        }
      }
    }
  }

  // Ignore depth information and average values
  for(int ir=0; ir<=nr; ir++) {
    read_line(infile);
  }
  for(int ir=0; ir<=nr; ir++) {
    read_line(infile);
  }

  //      double rhomax = -1.e12;
  //      double rhomin = 1.e12;

  // Read the tomography data
  for(int ir=0; ir<=nr; ir++) {
    for(int ll=0; ll<=lmax; ll++) {
      for(int mm=0; mm<=ll; mm++) {
        read_line(infile);
        for(int kk=0; kk<nk; kk++) {
          instream >> sphdata[kk][nr-ir][ll][mm][0];
          instream >> sphdata[kk][nr-ir][ll][mm][1];
        }
      }
    }
  }

  infile.close();

  shutil = new SphericalHarmonicsTool(lmax);

}


// Destructor
template<typename TYPENAME>
TomoVolumeFunction<TYPENAME>::~TomoVolumeFunction()
{

  for(int kk=0; kk<nk; kk++) {
    for(int ir=0; ir<=nr; ir++) {
      for(int ll=0; ll<=lmax; ll++) {
        delete[] sphdata[kk][ir][ll];
      }
      delete[] sphdata[kk][ir];
    }
    delete[] sphdata[kk];
  }
  delete[] sphdata;
  delete[] sphdata_buf;
  delete shutil;

}

template<typename TYPENAME> inline
TYPENAME TomoVolumeFunction<TYPENAME>::operator() ( TYPENAME x,
                                                    TYPENAME y,
                                                    TYPENAME z )
{
  return rho(x,y,z);
}

template<typename TYPENAME> inline
TYPENAME TomoVolumeFunction<TYPENAME>::operator() ( const Point3D & p)
{
  return rho( p[0], p[1], p[2] );
}


// Read the tomography
template<typename TYPENAME>
double TomoVolumeFunction<TYPENAME>::rho( double x, double y, double z ) {

  double rho;

  // delegate computation
  rho = shutil->shconvert_eval( sphdata, x, y, z, rmin, rmax, nr, lmax, 0 );

  // Non-dimensionalize rho with a previously given scaling factor.
  rho = rho * scaling_factor;

  return rho;
}

} // namespace hhg


