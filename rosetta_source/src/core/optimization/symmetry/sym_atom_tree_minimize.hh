// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/optimization/atom_tree_minimize.hh
/// @brief  Atom tree minimization functions
/// @author Ingemar Andre

#ifndef INCLUDED_core_optimization_symmetry_sym_atom_tree_minimize_hh
#define INCLUDED_core_optimization_symmetry_sym_atom_tree_minimize_hh


// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/symmetry/SymMinimizerMap.fwd.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.hh>
// AUTO-REMOVED #include <core/optimization/Multifunc.hh>
// AUTO-REMOVED #include <core/optimization/DOF_Node.hh>

// Symmetry
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetricConformation.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.fwd.hh>

//Auto Headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/optimization/MinimizerMap.fwd.hh>
#include <core/optimization/Multifunc.fwd.hh>



namespace core {
namespace optimization {
namespace symmetry {

//typedef conformation::symmetry::SymmetricConformation SymmetricConformation;
//typedef conformation::symmetry::SymmetryInfoCOP SymmetryInfoCOP;

void
atom_tree_dfunc(
	pose::Pose & pose,
	SymMinimizerMap & symm_min_map,
	scoring::ScoreFunction const & scorefxn,
	Multivec const & vars,
	Multivec & dE_dvars
);


void
atom_tree_get_atompairE_deriv(
	pose::Pose & pose,
	SymMinimizerMap & symm_min_map,
	scoring::ScoreFunction const & scorefxn
);


void
numerical_derivative_check(
	SymMinimizerMap const & min_map,
	Multifunc const & func,
	Multivec const & start_vars,
	Multivec const & dE_dvars,
	bool const verbose // = true
);


// Real
// calculate_direct_dof_derivatives(
// 	DOF_Node const & tor,
// 	pose::Pose const & pose,
// 	scoring::ScoreFunction const & scorefxn,
// 	ObjexxFCL::FArray2D< Real > const & dunbrack_deriv // currently this is pre-computed
// );

} // symmetry
} // namespace optimization
} // namespace core

#endif  // INCLUDED_core_optimization_symmetry_sym_atom_tree_minimize_hh
