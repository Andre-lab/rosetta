// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/zernike/ZernikeDescriptor
/// @brief      Zernike transform class
/// @details    ... TO FILL
/// @author     Novotni and Klein (Ingemar Andre)

#include <core/scoring/shape/zernike/ZernikeDescriptor.hh>
#include <iostream>

namespace core {
namespace scoring {
namespace shape {
namespace zernike {

template<class T, class TIn>
ZernikeDescriptor<T, TIn>::ZernikeDescriptor() {}

template<class T, class TIn>
ZernikeDescriptor<T, TIn>::ZernikeDescriptor (T* _voxels, int _dim, int _order) :
        voxels_ (_voxels), dim_ (_dim), order_ (_order)
{
    ComputeNormalization ();
    NormalizeGrid ();

    ComputeMoments ();
    ComputeInvariants ();
}

template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeMoments ()
{
    gm_.Init (voxels_, dim_, dim_, dim_, xCOG_, yCOG_, zCOG_, scale_, order_);
    //gm_.SetTransform (xCOG_, yCOG_, zCOG_, scale_);
    //gm_.Compute ();

    // Zernike moments
    zm_.Init (order_, gm_);
    zm_.Compute ();
}

/**
 * Cuts off the function : the object is mapped into the unit ball according to
* the precomputed center of gravity and scaling factor. All the voxels remaining
* outside the unit ball are set to zero.
*/
template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::NormalizeGrid ()
{
    T point[3];

    // it is easier to work with squared radius -> no sqrt required
    T radius = (T)1 / scale_;
    T sqrRadius = radius * radius;

    for (int x=0; x<dim_; ++x)
    {
        for (int y=0; y<dim_; ++y)
        {
            for (int z=0; z<dim_; ++z)
            {
                if (voxels_[(z*dim_ + y)*dim_ + x] != (T)0)
                {
                    point[0] = (T)x - xCOG_;
                    point[1] = (T)y - yCOG_;
                    point[2] = (T)z - zCOG_;

                    T sqrLen = point[0]*point[0] + point[1]*point[1] + point[2]*point[2];
                    if (sqrLen > sqrRadius)
                    {
                        voxels_[(z*dim_ + y)*dim_ + x] = 0.0;
                    }
                }
            }
        }
    }
}

/**
* Center of gravity and a scaling factor is computed according to the geometrical
* moments and a bounding sphere around the cog.
*/
template<class T, class TIn>
void ZernikeDescriptor<T, TIn>::ComputeNormalization ()
{
    ScaledGeometricalMoments<T, T> gm (voxels_, dim_, 0.0, 0.0, 0.0, 1.0);

    // compute the geometrical transform for no translation and scaling, first
    // to get the 0'th and 1'st order properties of the function
    //gm.Compute ();

    // 0'th order moments -> normalization
    // 1'st order moments -> center of gravity
    zeroMoment_ = gm.GetMoment (0, 0, 0);
    xCOG_ = gm.GetMoment (1, 0, 0) / zeroMoment_;
    yCOG_ = gm.GetMoment (0, 1, 0) / zeroMoment_;
    zCOG_ = gm.GetMoment (0, 0, 1) / zeroMoment_;

    // scaling, so that the function gets mapped into the unit sphere

    //T recScale = ComputeScale_BoundingSphere (voxels_, dim_, xCOG_, yCOG_, zCOG_);
    T recScale = 2.0 * ComputeScale_RadiusVar (voxels_, dim_, xCOG_, yCOG_, zCOG_);
    if (recScale == 0.0)
    {
        std::cerr << "\nNo voxels in grid!\n";
        exit (-1);
    }
    scale_ = (T)1 / recScale;
}


template<class T, class TIn>
double ZernikeDescriptor<T, TIn>::ComputeScale_BoundingSphere (T* _voxels, int _dim, T _xCOG, T _yCOG, T _zCOG)
{
    T max = (T)0;

    // the edge length of the voxel grid in voxel units
    int d = _dim;

    for (int x=0; x<d; ++x)
    {
        for (int y=0; y<d; ++y)
        {
            for (int z=0; z<d; ++z)
            {
                if (_voxels[(z + d * y) * d + x] > 0.9)
                {
                    T mx = (T)x - _xCOG;
                    T my = (T)y - _yCOG;
                    T mz = (T)z - _zCOG;
                    T temp = mx*mx + my*my + mz*mz;

                    if (temp > max)
                    {
                        max = temp;
                    }
                }
            }
        }
    }

    return std::sqrt (max);
}

template<class T, class TIn>
double ZernikeDescriptor<T, TIn>::ComputeScale_RadiusVar (T* _voxels, int _dim, T _xCOG, T _yCOG, T _zCOG)
{
    // the edge length of the voxel grid in voxel units
    int d = _dim;

    int nVoxels = 0;

    T sum = 0.0;

    for (int x=0; x<d; ++x)
    {
        for (int y=0; y<d; ++y)
        {
            for (int z=0; z<d; ++z)
            {
                if (_voxels[(z + d * y) * d + x] > 0.9)
                {
                    T mx = (T)x - _xCOG;
                    T my = (T)y - _yCOG;
                    T mz = (T)z - _zCOG;
                    T temp = mx*mx + my*my + mz*mz;

                    sum += temp;

                    nVoxels++;
                }
            }
        }
    }

    T retval = sqrt(sum/nVoxels);

    return retval;
}


/**
* Computes the Zernike moment based invariants, i.e. the norms of vectors with
* components of Z_nl^m with m being the running index.
*/
template<class T, class TIn>
void ZernikeDescriptor<T,TIn>::ComputeInvariants ()
{
    //invariants_.resize (order_ + 1);
                    moments_ = std::vector<std::complex<double>>();
    invariants_.clear ();
    for (int n=0; n<order_+1; ++n)
    {
        //invariants_[n].resize (n/2 + 1);

        //T sum = (T)0; BUG!!!
        //int l0 = n % 2, li = 0;
        int li = 0;

        for (int l = n % 2; l<=n; ++li, l+=2)
        {
                T sum = (T)0; // Should be here instead!
            for (int m=-l; m<=l; ++m)
            {
                ComplexT moment = zm_.GetMoment (n, l, m);
                sum += std::norm (moment);
                                                moments_.push_back(moment);
            }

            invariants_.push_back (sqrt (sum));
            //invariants_[n][li] = std::sqrt (sum);
        }
    }
}

/** Implemented by Mads
 * Simple wrapper for the original reconstruction function by Novotni and Klein
 */
template<class T, class TIn>
std::vector<std::vector<std::vector<std::complex<double>>>>
ZernikeDescriptor<T,TIn>::ReconstructionWrapper () {
    const int gridSize = dim_;
    ComplexT3D grid(gridSize, vector<vector<ComplexT>>(gridSize, vector<ComplexT>(gridSize)));
    zm_.Reconstruct(grid, xCOG_,  yCOG_, zCOG_, scale_);
return grid;
}


//assert (invariants_.size () == (order_ + 1));

// write the invariants
//   for (int n=0; n<order_+1; ++n)
//   {
//       int l0 = n % 2, li = 0;
//       for (int l=l0; l<=n; l+=2, ++li)
//       {
//           temp = (float)invariants_[n][li];
//           outfile.write ((char*)(&temp), sizeof(float));
//       }
//   }
//
//    int dim = invariants_.size ();
//    outfile.write ((char*)(&dim), sizeof(int));
//
//    for (int i=0; i<dim; ++i)
//    {
//        temp = invariants_[i];
//	 std::cout << temp << std::endl;
//        outfile.write ((char*)(&temp), sizeof(float));
//    }

template<class T, class TIn>
void ZernikeDescriptor<T,TIn>::SaveInvariants (const char* _fName)
{
//    std::ofstream outfile (_fName, std::ios_base::binary | std::ios_base::out);
    std::ofstream outfile (_fName );

    float temp;

    int dim = invariants_.size ();
    outfile << "Calculated from voxel with dim=" << dim_ << " coefficients:" << dim << std::endl;
    for (int i=0; i<dim; ++i)
    {
        temp = invariants_[i];
        outfile << temp << std::endl;
    }
}

template class ZernikeDescriptor<double,double>;
}
}
}
}
