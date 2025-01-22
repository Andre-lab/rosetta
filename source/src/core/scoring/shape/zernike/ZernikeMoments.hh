// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/zernike/ZernikeMoments
/// @brief      Zernike transform class
/// @details    ... TO FILL
/// @author     Novotni and Klein (Ingemar Andre)


#ifndef INCLUDED_core_scoring_shape_zernike_ZernikeMoments_hh
#define INCLUDED_core_scoring_shape_zernike_ZernikeMoments_hh

// core
#include <core/types.hh>

// C++ headers
#include <vector>
#include <complex>
#include <set>
#include <ios>

// ----- local program includes -----
#include "core/scoring/shape/zernike/Factorial.hh"
#include "core/scoring/shape/zernike/Binomial.hh"
#include "core/scoring/shape/zernike/ScaledGeometricMoments.hh"

#define PI 3.141592653589793

namespace core {
namespace scoring {
namespace shape {
namespace zernike {

/**
* Struct representing a complex coefficient of a moment
* of order (p_,q_,r_)
*/
template<class T>
struct ComplexCoeff
{
    typedef     std::complex<T>     ValueT;

    ComplexCoeff (int, int, int, const ValueT&);
    ComplexCoeff (const ComplexCoeff<T>& _cc);
    ComplexCoeff ();

    int                 p_, q_, r_;
    ValueT     value_;
};

/**
* Class representing the Zernike moments
*/
template<class VoxelT, class MomentT>
class ZernikeMoments
{
public:
    // ---- public typedefs ----
    typedef MomentT             T;
    typedef vector<T>           T1D;        // vector of scalar type
    typedef vector<T1D>         T2D;        // 2D array of scalar type
    typedef vector<T2D>         T3D;        // 3D array of scalar type
    typedef vector<T3D>         T4D;        // 3D array of scalar type

    typedef std::complex<T>                      ComplexT;       // complex type
    typedef vector<vector<vector<ComplexT> > >   ComplexT3D;     // 3D array of complex type

    typedef ComplexCoeff<T>                      ComplexCoeffT;
    typedef vector<vector<vector<vector<ComplexCoeffT> > > >    ComplexCoeffT4D;


public:
    // ---- public member functions ----
    ZernikeMoments (int _order, ScaledGeometricalMoments<VoxelT,MomentT>& _gm);
    ZernikeMoments ();

    void Init (int _order, ScaledGeometricalMoments<VoxelT,MomentT>& _gm);
    void Compute ();

    ComplexT GetMoment (int _n, int _l, int _m);

    // ---- debug functions/arguments ----
    void Reconstruct (ComplexT3D&   _grid,                // grid containing the reconstructed function
                      T             _xCOG,                // center of gravity
                      T             _yCOG,
                      T             _zCOG,
                      T             _scale,               // scaling factor to map into unit ball
                      int           _minN = 0,            // min value for n freq index
                      int           _maxN = 100,          // min value for n freq index
                      int           _minL = 0,            // min value for l freq index
                      int           _maxL = 100);         // max value for l freq index


    void SaveGrid (ComplexT3D& _grid, std::string filename);

    void NormalizeGridValues (ComplexT3D& _grid);
    void CheckOrthonormality (int _n1, int _l1, int _m1, int _n2, int _l2, int _m2);

private:
    // ---- private member functions ----
    void ComputeCs ();
    void ComputeQs ();
    void ComputeGCoefficients ();

    // ---- private attributes -----
    ComplexCoeffT4D     gCoeffs_;           // coefficients of the geometric moments
    ComplexT3D          zernikeMoments_;    // nomen est omen
    T3D                 qs_;                // q coefficients (radial polynomial normalization)
    T2D                 cs_;                // c coefficients (harmonic polynomial normalization)

    ScaledGeometricalMoments<VoxelT,MomentT> gm_;
    int                 order_;             // := max{n} according to indexing of Zernike polynomials

    // ---- debug functions/arguments ----
    void PrintGrid (ComplexT3D& _grid);
    T EvalMonomialIntegral (int _p, int _q, int _r, int _dim);
};

} // namespace zernike
} // namespace shape
} // namespace scoring
} // namespace core

#endif //INCLUDED_core_scoring_shape_zernike_ZernikeMoments_hh


