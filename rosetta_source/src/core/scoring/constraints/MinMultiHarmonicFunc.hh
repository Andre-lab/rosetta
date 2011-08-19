// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/HarmonicFunc.hh
/// @brief Definition for functions used in definition of constraints.
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_core_scoring_constraints_MinMultiHarmonicFunc_hh
#define INCLUDED_core_scoring_constraints_MinMultiHarmonicFunc_hh

// AUTO-REMOVED #include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

//Auto Headers
#include <utility/vector1_bool.hh>


// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

class MinMultiHarmonicFunc : public Func {
public:
	MinMultiHarmonicFunc( utility::vector1<Real> const & x0_in, utility::vector1<Real> const & sd_in );

	FuncOP
	clone() const { return new MinMultiHarmonicFunc( *this ); }

	Real func( Real const x ) const;
	Real dfunc( Real const x ) const;

	void read_data( std::istream & in );

	void show_definition( std::ostream &out ) const;

	Real x0(int n) const {
		return x0_[n];
	}

	Real sd(int n) const {
		return sd_[n];
	}

	void x0( int n, Real x ) {
		x0_[n] = x;
	}

	void sd( int n, Real sd ) {
		sd_[n] = sd;
	}

	Size
	show_violations( std::ostream& out, Real x, Size verbose_level, core::Real threshold = 1 ) const;

private:
	utility::vector1<Real> x0_;
	utility::vector1<Real> sd_;
	Size n_;
	mutable Size which_component_;
};

} // constraints
} // scoring
} // core

#endif
