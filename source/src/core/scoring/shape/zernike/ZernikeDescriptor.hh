// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/zernike/ZernikeTransform
/// @brief      Zernike transform class
/// @details    ... TO FILL
/// @author     Novotni and Klein (Ingemar Andre)


#ifndef INCLUDED_core_scoring_shape_zernike_ZernikeDescriptor_hh
#define INCLUDED_core_scoring_shape_zernike_ZernikeDescriptor_hh

// Utility headers
#include <ObjexxFCL/FArray3D.hh>

// core
#include <core/types.hh>

// C++ headers
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>

// other
#include "core/scoring/shape/zernike/ScaledGeometricMoments.hh"
#include "core/scoring/shape/zernike/ZernikeMoments.hh"

namespace core {
namespace scoring {
namespace shape {
namespace zernike {

/**
* This class serves as a wrapper around the geometrical and
* Zernike moments. It provides also the implementation of invariant Zernike
* descriptor_out, means of reconstruction_out of orig. function, etc.
*/
template<class T, class TIn>
class ZernikeDescriptor
{
public:
    // ---- exported typedefs ----
    /// complex type
    typedef std::complex<T>                         ComplexT;
    /// 3D array of complex type
    typedef vector<vector<vector<ComplexT> > >      ComplexT3D;


    typedef vector<T>                               T1D;
    typedef vector<T1D>                             T2D;            // 2D array of T type

    //typedef CumulativeMoments<T, T>                 CumulativeMomentsT;
    typedef ScaledGeometricalMoments<T, T>          ScaledGeometricalMomentsT;
    typedef ZernikeMoments<T, T>                    ZernikeMomentsT;

public:

    ZernikeDescriptor (
            T* _voxels,                 /**< the cubic voxel get_grid */
            int _dim,                   /**< dimension is $_dim^3$ */
            int _order                  /**< maximal order of the Zernike moments (N in paper) */
    );
    ZernikeDescriptor ();

    std::vector<std::vector<std::vector<std::complex<double>>>>
    ReconstructionWrapper();

    /**
    * Saves the computed invariants into a binary file
    */
    void SaveInvariants (
            const char* _fName      /**< name of the output file */
    );

    /// Access to invariants
    /// MADS: I HAVE CHANGED  THIS TO RETURN A REFERENCE INSTEAD
    T1D & GetInvariants (){ return invariants_; };
    std::vector<std::complex<double>> GetMoments (){ return moments_; };

private:
    // ---- private helper functions ----
    void NormalizeGrid ();
    void ComputeNormalization ();
    void ComputeMoments ();
    void ComputeInvariants ();

    double ComputeScale_BoundingSphere (
            T* _voxels,
            int _dim,
            T _xCOG,
            T _yCOG,
            T _zCOG
    );
    double ComputeScale_RadiusVar (
            T* _voxels,
            int _dim,
            T _xCOG,
            T _yCOG,
            T _zCOG
    );


private:
    // ---- member variables ----
    T*      voxels_;                // 1D array containing the voxels
    int     dim_;                   // length of the edge of the voxel get_grid (which is a cube)
    int     order_;                 // maximal order of the moments to be computed (max{n})
    T       zeroMoment_,            // zero order moment
            xCOG_, yCOG_, zCOG_,    // center of gravity
            scale_;                 // scaling factor mapping the function into the unit sphere

    //T2D                 invariants_;        // 2D vector of invariants under SO(3)
    T1D                 invariants_;        // 2D vector of invariants under SO(3)
	  std::vector<std::complex<double>> moments_;
    ZernikeMomentsT     zm_;
    //CumulativeMomentsT  cm_;
    ScaledGeometricalMomentsT gm_;
};

} // namespace zernike
} // namespace shape
} // namespace scoring
} // namespace core

#endif //INCLUDED_core_scoring_shape_zernike_ZernikeDescriptor_hh

