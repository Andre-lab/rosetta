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


#ifndef INCLUDED_core_scoring_shape_zernike_Binomial_hh
#define INCLUDED_core_scoring_shape_zernike_Binomial_hh

// core
#include <core/types.hh>

// C++ headers
#include <vector>
#include <assert.h>

namespace core {
    namespace scoring {
        namespace shape {
            namespace zernike {

                using std::vector;

                /**
                * A template class facilitating fast computation and retrieval of binomials.
                * The binomials are computed only once at first call of Get function according
                * to Pascal's Triangle.
                */
                template<class T>
                class Binomial
                {
                public:
                    typedef vector<T>           VectorT;
                    typedef vector<VectorT>     VVectorT;

                    /** Retrieves the binomial _i "over" _j */
                    static T Get (int _i, int _j);
                    /** Sets the maximal value of upper binomial param to _max */
                    static void SetMax (int _max);
                    /** Gets the maximal value of upper binomial param */
                    static int GetMax ();

                private:
                    /** Computed Pascal's Triangle */
                    static void ComputePascalsTriangle ();

                    static VVectorT pascalsTriangle_;
                    static int max_;
                };


                template<class T>
                inline void Binomial<T>::SetMax (int _max)
                {
                    max_ = _max;
                    ComputePascalsTriangle ();
                }


                template<class T>
                inline int Binomial<T>::GetMax ()
                {
                    return max_;
                }


                template<class T>
                typename Binomial<T>::VVectorT Binomial<T>::pascalsTriangle_;

                template<class T>
                int Binomial<T>::max_ = 60;


            } // namespace zernike
        } // namespace shape
    } // namespace scoring
} // namespace core

#endif //INCLUDED_core_scoring_shape_zernike_Binomial_hh



