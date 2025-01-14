// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/zernike/Binomial
/// @brief      Zernike transform class
/// @details    ... TO FILL
/// @author     Novotni and Klein (Ingemar Andre)

#include <core/scoring/shape/zernike/Binomial.hh>

namespace core {
    namespace scoring {
        namespace shape {
            namespace zernike {

                template<class T>
                void Binomial<T>::ComputePascalsTriangle ()
                {
                    // allocate storage for the pascal triangle of size determined by max_
                    pascalsTriangle_.resize (max_+1);

                    for (int i=0; i<max_+1; ++i)
                    {
                        pascalsTriangle_[i].resize (max_+1-i);
                        for (int j=0; j<max_+1-i; ++j)
                        {
                            // the values are ones on the edges of the triangle
                            if (!i || !j)
                            {
                                pascalsTriangle_[i][j] = (T)1;
                            }
                                // use the familiar addition to generate values on lower levels
                            else
                            {
                                pascalsTriangle_[i][j] = pascalsTriangle_[i][j-1] + pascalsTriangle_[i-1][j];
                            }
                        }
                    }
                }

                template<class T>
                T Binomial<T>::Get (int _i, int _j)
                {
                    // the values are computed only first time this function is called
                    if (!pascalsTriangle_.size ())
                    {
                        ComputePascalsTriangle ();
                    }

                    assert (_i>=0 && _j>=0 && _i>=_j);

                    return pascalsTriangle_[_j][_i-_j];
                }
                template class Binomial<double>;
            }
        }
    }
}

