// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/zernike/ScaledGeometricMoments.hh
/// @brief      Zernike transform class
/// @details    ... TO FILL
/// @author     Novotni and Klein (Ingemar Andre)


#ifndef INCLUDED_core_scoring_shape_zernike_ScaledGeometricMoments_hh
#define INCLUDED_core_scoring_shape_zernike_ScaledGeometricMoments_hh

// core
#include <core/types.hh>

// C++ headers
#include <vector>

namespace core {
    namespace scoring {
        namespace shape {
            namespace zernike {

                using std::vector;

                /**
                 Class for computing the scaled, pre-integrated geometrical moments.
                 These tricks are needed to make the computation numerically stable.
                 See the paper for more details.
                 \param VoxelT   type of the voxel values
                 \param MomentT  type of the moments -- recommended to be double
                 */
                template<class VoxelT, class MomentT>
                class ScaledGeometricalMoments
                {
                public:
                    // ---- public typedefs ----
                    /// the moment type
                    typedef MomentT             T;
                    /// vector scalar type
                    typedef vector<T>           T1D;
                    /// 2D array scalar type
                    typedef vector<T1D>         T2D;
                    /// 3D array scalar type
                    typedef vector<T2D>         T3D;
                    /// vector scalar type
                    typedef vector<double>      Double1D;
                    /// vector scalar type
                    typedef vector<Double1D>    Double2D;

                    typedef typename T1D::iterator       T1DIter;

                    // ----- public methods -----

                    // ---- construction / init ----
                    /// Contructor
                    ScaledGeometricalMoments (
                            const VoxelT* _voxels,  /**< input voxel grid */
                            int _xDim,              /**< x-dimension of the input voxel grid */
                            int _yDim,              /**< y-dimension of the input voxel grid */
                            int _zDim,              /**< z-dimension of the input voxel grid */
                            double _xCOG,           /**< x-coord of the center of gravity */
                            double _yCOG,           /**< y-coord of the center of gravity */
                            double _zCOG,           /**< z-coord of the center of gravity */
                            double _scale,          /**< scaling factor */
                            int _maxOrder = 1       /**< maximal order to compute moments for */
                    );

                    /// Constructor for equal dimensions for each axis
                    ScaledGeometricalMoments (
                            const VoxelT* _voxels,  /**< input voxel grid */
                            int _dim,               /**< the grid is _dim^3 */
                            double _xCOG,           /**< x-coord of the center of gravity */
                            double _yCOG,           /**< y-coord of the center of gravity */
                            double _zCOG,           /**< z-coord of the center of gravity */
                            double _scale,          /**< scaling factor */
                            int _maxOrder = 1       /**< maximal order to compute moments for */
                    );

                    /// Default constructor
                    ScaledGeometricalMoments ();

                    /// The init function used by the contructors
                    void Init (
                            const VoxelT* _voxels,  /**< input voxel grid */
                            int _xDim,              /**< x-dimension of the input voxel grid */
                            int _yDim,              /**< y-dimension of the input voxel grid */
                            int _zDim,              /**< z-dimension of the input voxel grid */
                            double _xCOG,           /**< x-coord of the center of gravity */
                            double _yCOG,           /**< y-coord of the center of gravity */
                            double _zCOG,           /**< z-coord of the center of gravity */
                            double _scale,          /**< scaling factor */
                            int _maxOrder = 1       /**< maximal order to compute moments for */
                    );


                    /// Access function
                    T GetMoment (
                            int _i,                 /**< order along x */
                            int _j,                 /**< order along y */
                            int _k                  /**< order along z */
                    );

                private:
                    int xDim_,              // dimensions
                            yDim_,
                            zDim_,
                            maxOrder_;          // maximal order of the moments

                    T2D         samples_;   // samples of the scaled and translated grid in x, y, z
                    T1D         voxels_;    // array containing the voxel grid
                    T3D         moments_;   // array containing the cumulative moments

                    // ---- private functions ----
                    void Compute ();
                    void ComputeSamples (double _xCOG, double _yCOG, double _zCOG, double _scale);
                    void ComputeDiffFunction (T1DIter _iter, T1DIter _diffIter, int _dim);

                    T Multiply (T1DIter _diffIter, T1DIter _sampleIter, int _dim);

                };

            } // namespace zernike
        } // namespace shape
    } // namespace scoring
} // namespace core

#endif //INCLUDED_core_scoring_shape_zernike_ScaledGeometricMoments_hh


