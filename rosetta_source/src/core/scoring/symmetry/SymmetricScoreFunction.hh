// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   core/scoring/symmetry/SymmetricScoreFunction.hh
/// @brief  Symmetric Score function class
/// @author Ingemar Andre


#ifndef INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_hh
#define INCLUDED_core_scoring_symmetry_SymmetricScoreFunction_hh

// Unit headers
#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>

// Project headers
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

namespace core {
namespace scoring {
namespace symmetry {

class SymmetricScoreFunction : public ScoreFunction
{
public:
	typedef ScoreFunction parent;
	typedef conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef conformation::symmetry::SymmetryInfoCOP SymmetryInfoCOP;

public:

	/// ctor
	SymmetricScoreFunction();

	SymmetricScoreFunction &
	operator=( SymmetricScoreFunction const & );

	SymmetricScoreFunction( SymmetricScoreFunction const & );

	SymmetricScoreFunction( ScoreFunction const & src );

	SymmetricScoreFunction( ScoreFunctionOP src );

	SymmetricScoreFunction( ScoreFunctionCOP src );

	ScoreFunctionOP clone() const;

  /////////////////////////////////////////////////////////////////////////////
  // score
  /////////////////////////////////////////////////////////////////////////////

	virtual Real
	operator ()( pose::Pose & pose ) const;

	/// @brief Initialize a MinimizationGraph and cache it in the pose's Energies object
	/// for use during minimization -- only add edges to the asymmetric unit and within it
	/// are added to the MinimizationGraph.
	virtual
	void
	setup_for_minimizing(
		pose::Pose & pose,
		kinematics::MinimizerMapBase const & min_map
	) const;

	///
	void
	eval_twobody_neighbor_energies( pose::Pose & pose ) const;

	void
	eval_long_range_twobody_energies( pose::Pose & pose ) const;

	///
	void
	eval_onebody_energies( pose::Pose & pose ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose ) const;

	virtual
	void
	eval_npd_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		Vector & F1,
		Vector & F2
	) const;

	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose
	) const;

	//void create_intersubunit_hbonds( pose::Pose & pose, hbonds::HBondSetOP hbond_set_subunit ) const;

	void
	intersubunit_hbond_energy( pose::Pose & pose, EnergyMap & intersubunit_energy  ) const;

	void
	symmetrical_allow_hbonds( pose::Pose & pose ) const;

	void
	set_symmetric_residue_neighbors_hbonds( pose::Pose & pose ) const;

	void
	set_symmetric_cenlist( pose::Pose & pose ) const;

	void
	correct_arrays_for_symmetry( pose::Pose & pose ) const;

	void
	correct_finalize_score( pose::Pose & pose ) const;

};


} // symmetry
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_SymmetricScoreFunction_HH
