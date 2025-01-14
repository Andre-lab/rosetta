
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/zernike/Factorial.hh
/// @brief      Zernike transform class
/// @details    ... TO FILL
/// @author     Novotni and Klein (Ingemar Andre)


#ifndef INCLUDED_core_scoring_shape_zernike_Factorial_hh
#define INCLUDED_core_scoring_shape_zernike_Factorial_hh

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
                * A template class for precomputation and subsequent retrieval of factorials
                * of an integer number.
                *
                * The maximal input parameter is set to 19 at first, since double and __int64
                * can represent exactly numbers of 18 or less digits. There is an option to change this
                * in case one uses a data type being capable of representing bigger numbers.
                */
                template<class T>
                class Factorial
                {
                public:
                    /** Gets the factorial of _i */
                    static T Get (int _i);
                    /** Gets _i*(_i+1)*...*(_j-1)*_j */
                    static T Get (int _i, int _j);
                    /** Sets the maximal stored factorial value to _max */
                    static void SetMax (int _max);
                    /** Gets the maximal stored factorial value */
                    static int GetMax ();

                private:
                    /** Computes factorials of numbers [1..max] */
                    static void ComputeFactorials ();

                    static int max_;
                    static vector<T>           factorials_;
                    static vector<vector<T> >  subFactorials_;
                };

                // The obligatory initialization of static attributes
                template<class T>
                int Factorial<T>::max_ = 19;

                template<class T>
                vector<T> Factorial<T>::factorials_;

                /**
                * Computes factorials up to max_ and stores them internally
                */
                template<class T>
                inline void Factorial<T>::ComputeFactorials ()
                {
                    factorials_.resize (max_);

                    factorials_[0] = (T)1;

                    for (int i=1; i<max_; ++i)
                    {
                        factorials_[i] = factorials_[i-1] * (T)(i+1);
                    }
                }

                /**
                * Retrieves the factorial of _i. All factorials are computed only the first
                * time this function is called, after this they are just read from the store.
                */
                template<class T>
                inline T Factorial<T>::Get (int _i)
                {
                    assert (_i >= 0 && _i <= max_);

                    if (!factorials_.size ())
                    {
                        ComputeFactorials ();
                    }

                    // 0! = 1
                    if (!_i)
                    {
                        return 1;
                    }

                    return factorials_[_i-1];
                }

                template<class T>
                inline T Factorial<T>::Get (int _i, int _j)
                {
                    T result = (T)1;

                    for (int i=_j; i>=_i; --i)
                    {
                        result *= i;
                    }

                    return result;
                }

                /*
                * Modifies the maximum factorial input parameter. All factorials are recomputed here
                */
                template<class T>
                inline void Factorial<T>::SetMax (int _max)
                {
                    assert (_max >= (T)0);

                    // In fact, the previously computed factorials could be reused here,
                    // however, since this is rarely performed and takes only a couple of
                    // multiplications, we don't care.
                    max_ = _max;
                    if (max_ <= _max)
                    {
                        ComputeFactorials ();
                    }
                }


                template<class T>
                inline int Factorial<T>::GetMax ()
                {
                    return max_;
                }


            } // namespace zernike
        } // namespace shape
    } // namespace scoring
} // namespace core

#endif //INCLUDED_core_scoring_shape_zernike_Factorial_hh



