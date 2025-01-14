// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/zernike/Zernike2Ddescriptor.hh
/// @brief      2D Zernike transform class
/// @details    ... TO FILL
/// @author     Ingemar Andre


#ifndef INCLUDED_core_scoring_shape_zernike_Zernike2Ddescriptor_hh
#define INCLUDED_core_scoring_shape_zernike_Zernike2Ddescriptor_hh

// C++ headers
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <complex>

// core
#include <core/types.hh>

namespace core {
namespace scoring {
namespace shape {
namespace zernike {

class Zernike2Ddescriptor
{
public:

public:

    Zernike2Ddescriptor (
            std::vector<std::vector<int>> _grid,                 /*The voxel surface*/
            int _dim,                   /**< dimension of 2D grid */
            int _order                  /**< maximal order of the Zernike moments */
    );

    Zernike2Ddescriptor (
            int _dim,                   /**< dimension of 2D grid */
            int _order                  /**< maximal order of the Zernike moments */
    );

    Zernike2Ddescriptor ();
//    /**
//   * Saves the computed invariants into a binary file
//    */
//    void SaveInvariants (
//            const char* _fName      /**< name of the output file */
//												 ) {};

    std::vector<core::Real> & GetInvariants (){ return invariants_; };
    std::vector<std::complex<double>> & GetMoments (){ return moments_; };

    bool initialized();

    void Transform( std::vector<std::vector<int>> grid );

private:

void setup();

double get_r(int x,int y, int N);

unsigned long int factorial(int n);

int count_nonzero(std::vector<std::vector<double>> & rho);

std::vector< std::vector<int> >
get_nm(int n);

std::vector<std::vector<double>>
rad_poly(std::vector<std::vector<double>> &rho,int n,int m);


private:
    // ---- member variables ----
    std::vector<std::vector<int>> grid_;                // 1D array containing the voxels
    int     dim_;                   // length of the edge of the voxel grid (which is a cube)
    int     order_;                 // maximal order of the moments to be computed (max{n})

    std::vector<core::Real> invariants_;        // 2D vector of invariants under SO(3)
    std::vector<std::complex<double>> moments_;        // 2D vector of moments under SO(3)

   std::vector<std::vector<int>> nm_;
    std::vector<std::vector<double>> rho_;
    std::vector<std::vector<double>> theta_;
		std::vector< std::vector<std::vector<double>> > rad_poly_store_;
    bool initialized_;
};

} // namespace zernike
} // namespace shape
} // namespace scoring
} // namespace core

#endif //INCLUDED_core_scoring_shape_zernike_Zernike2Ddescriptor_hh

