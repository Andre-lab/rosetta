// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_SugarCloseEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_SugarCloseEnergy_HH
#define INCLUDED_core_scoring_rna_RNA_SugarCloseEnergy_HH

// Unit headers
#include <core/scoring/rna/RNA_SugarCloseEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace rna {

///
class RNA_SugarCloseEnergy : public core::scoring::methods::ContextIndependentOneBodyEnergy  {
public:
	typedef core::scoring::methods::ContextIndependentOneBodyEnergy  parent;

public:

	/// @brief ctor
	RNA_SugarCloseEnergy();

	/// @brief dtor
	virtual ~RNA_SugarCloseEnergy();

	/// clone
	virtual
	core::scoring::methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	///
	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		EnergyMap & emap
	) const;


	///
	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &,
		EnergyMap & emap
	) const;


	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	///
///////////////////////////////////////////////////////////////////////////////
	virtual
	void
	eval_atom_derivative(
											 id::AtomID const & id,
											 pose::Pose const & pose,
											 kinematics::DomainMap const &, // domain_map,
											 ScoreFunction const & sfxn,
											 EnergyMap const & weights,
											 Vector & F1,
											 Vector & F2
											 ) const;


	/// @brief RNA_SugarCloseEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	void
	setup_sugar_ring_closure_constraints( pose::Pose & pose ) const;

	void
	add_sugar_ring_closure_constraints( conformation::Residue const & rsd, constraints::ConstraintSet & cst_set ) const;

	// data
private:

	// Ribose closure
	Real const scale_rna_torsion_tether_;
	Real const scale_rna_torsion_sd_;
	Distance const o4star_c1star_bond_length_;
	Distance const o4star_c1star_sd_;
	constraints::HarmonicFuncOP o4star_c1star_dist_harm_func_;

	Real const angle_sd_;
	Real const o4star_c1star_c2star_bond_angle_;
	constraints::HarmonicFuncOP o4star_c1star_c2star_angle_harm_func_;
	Real const o4star_c1star_first_base_bond_angle_;
	constraints::HarmonicFuncOP o4star_c1star_first_base_angle_harm_func_;
	Real const c4star_o4star_c1star_bond_angle_;
	constraints::HarmonicFuncOP c4star_o4star_c1star_angle_harm_func_;

	mutable constraints::ConstraintSetOP rna_sugar_close_constraints_;

	virtual
	core::Size version() const;


};

} // rna
} // scoring
} // core


#endif // INCLUDED_core_scoring_EtableEnergy_HH
