// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/TopOutFunc.hh
/// @brief Implementation of phenix "top-out" function
///   Similar to Geman-McClure: harmonic near 'x0_', flat past 'limit_'
/// @author Frank DiMaio


#include <core/scoring/constraints/TopOutFunc.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <sstream>
#include <cmath>

// C++ Headers
namespace core {
namespace scoring {
namespace constraints {

Real
TopOutFunc::func( Real const x ) const {
	Real xoff = x-x0_;
	Real top = weight_ * limit_ * limit_;
	Real score = top * (1.0 - exp(-weight_*xoff*xoff/top));
	return score;
}

Real
TopOutFunc::dfunc( Real const x ) const {
	Real xoff = x-x0_;
	Real top = weight_ * limit_ * limit_;
	Real grad = 2.0 * weight_ * xoff * exp(-(weight_*xoff*xoff)/top);
	return (grad);
}

void
TopOutFunc::read_data( std::istream& in ) {
	in >> weight_ >> x0_ >> limit_;
}

void
TopOutFunc::show_definition( std::ostream &out ) const {
	out << "TOPOUT " << weight_ << " " << x0_ << " " << limit_ << std::endl;
}

Size
TopOutFunc::show_violations( std::ostream& out, Real x, Size verbose_level, Real threshold) const {
	if (verbose_level > 100 ) {
		out << "TOPOUT " <<  func(x) << std::endl;
	}
	return Func::show_violations( out, x, verbose_level, threshold);
}

} // namespace constraints
} // namespace scoring
} // namespace core

