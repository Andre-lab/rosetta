// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/constraints/SoedingFunc.hh
/// @brief Definition for Johannes Soeding-style distance restraints.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_constraints_SoedingFunc_hh
#define INCLUDED_core_scoring_constraints_SoedingFunc_hh

#include <core/scoring/constraints/Func.hh>
#include <core/types.hh>

// C++ Headers
#include <ostream>

namespace core {
namespace scoring {
namespace constraints {

/// @brief Derived class of class Func representing a Soeding distribution with a user-specified
/// mean and standard deviation.
class SoedingFunc : public Func {
public:

	SoedingFunc()
		: w1_(0.0), mean1_(0.0), sdev1_(0.0), w2_(0.0), mean2_(0.0), sdev2_(0.0)
	{}

	/// @brief returns a clone of this SoedingFunc
	FuncOP clone() const { return new SoedingFunc( *this ); }

	Real
	compute_func( Real const x ) const;

	/// @brief Returns the value of this SoedingFunc evaluated at distance x.
	Real func( Real const x ) const;

	/// @brief Returns the value of the first derivative of this SoedingFunc at distance x.
	Real dfunc( Real const x ) const;

	/// @brief show the definition of this SoedingFunc to the specified output stream.
	virtual void show_definition( std::ostream &out ) const;

	/// @brief Calls show( out ) on this SoedingFunc.
	friend std::ostream& operator<<(std::ostream& out, const SoedingFunc& f ) {
		f.show( out );
		return out;
	} // operator<<

	/// @brief Initializes this SoedingFunc from the given istream.
	void read_data( std::istream & in );

private:
	Real w1_, mean1_, sdev1_, w2_, mean2_, sdev2_;
};

} // constraints
} // scoring
} // core

#endif
