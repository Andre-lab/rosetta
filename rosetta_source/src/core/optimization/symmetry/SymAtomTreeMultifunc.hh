// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/optimization/AtomTreeMultifunc.hh
/// @brief  Atom tree multifunction class for symmetrical minimization
/// @author Ingemar Andre

#ifndef INCLUDED_core_optimization_symmetry_SymAtomTreeMultifunc_hh
#define INCLUDED_core_optimization_symmetry_SymAtomTreeMultifunc_hh

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.hh>
// AUTO-REMOVED #include <core/optimization/symmetry/sym_atom_tree_minimize.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/kinematics/Jump.hh>
#include <core/optimization/symmetry/SymMinimizerMap.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace optimization {
namespace symmetry {

/// @brief Atom tree multifunction class
class SymAtomTreeMultifunc : public Multifunc {

public: // Creation

	// c-tor
	SymAtomTreeMultifunc(
		pose::Pose & pose_in,
		SymMinimizerMap & symm_min_map,
		scoring::ScoreFunction const & scorefxn_in,
		bool const deriv_check_in = false,
		bool const deriv_check_verbose_in = false
	):
		pose_( pose_in ),
		symm_min_map_( symm_min_map ),
		score_function_( scorefxn_in ),
		deriv_check_( deriv_check_in ),
		deriv_check_verbose_( deriv_check_verbose_in )
	{}



	/// @brief Destructor
	inline
	virtual
	~SymAtomTreeMultifunc()
	{}


public: // Methods


	// func
	virtual
	Real
	operator ()( Multivec const & vars ) const;

	// dfunc
	virtual
	void
	dfunc( Multivec const & vars, Multivec & dE_dvars ) const;

	/// @brief Error state reached; dump out current pdb.
	virtual
	void
	dump( Multivec const & vars, Multivec const & vars2 ) const;

private: // data

	/// non-const since pose_ is modified by calls to operator()
	pose::Pose & pose_;

	/// non-const since min_map_ is modified by calls to dfunc()
	SymMinimizerMap & symm_min_map_;

	scoring::ScoreFunction const & score_function_;

	bool deriv_check_;
	bool deriv_check_verbose_;

}; // AtomTreeMultifunc

} // symmetry
} // namespace optimization
} // namespace core

#endif // INCLUDED_core_optimization_symmetry_SymAtomTreeMultifunc_HH
