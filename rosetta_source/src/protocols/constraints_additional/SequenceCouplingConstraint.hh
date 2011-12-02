// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraints_additional/SequenceCouplingConstraint.hh
/// @brief  This is a constraint that refers to a core::sequence::SequenceCoupling? in order to influence the scoring of amino acid types based on multiple sequence alignments (i.e. for biasing amino acid choices during design).
/// @author ashworth

#ifndef INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_hh
#define INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_hh

#include <protocols/constraints_additional/SequenceCouplingConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>

// AUTO-REMOVED #include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/sequence/SequenceCoupling.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

#include <core/sequence/SequenceCoupling.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace constraints_additional {

///@brief
class SequenceCouplingConstraint : public core::scoring::constraints::Constraint {
public:
	typedef core::sequence::SequenceCoupling SequenceCoupling;
	typedef core::sequence::SequenceCouplingOP SequenceCouplingOP;
	typedef core::sequence::SequenceCouplingCOP SequenceCouplingCOP;
	typedef core::id::SequenceMapping SequenceMapping;
	typedef core::pose::Pose Pose;
	typedef core::conformation::Conformation Conformation;
	typedef core::scoring::EnergyMap EnergyMap;
	typedef core::scoring::constraints::FuncFactory FuncFactory;
	typedef core::scoring::constraints::XYZ_Func XYZ_Func;
	typedef core::scoring::constraints::ConstraintOP ConstraintOP;

public:

	SequenceCouplingConstraint();

	SequenceCouplingConstraint(
		Pose const &,
		core::Size,
		core::Size,
		SequenceCouplingOP profile = NULL
	);

	SequenceCouplingConstraint(
		core::Size,
		core::Size,
		utility::vector1< AtomID > const &,
		SequenceCouplingOP profile = NULL
	);

	virtual ~SequenceCouplingConstraint();

	virtual ConstraintOP clone() const;

	virtual std::string type() const { return "SequenceCoupling"; }

	///@brief used by ConstraintIO and ConstraintFactory to construct this constraint from a input file stream (constraint file)
	virtual void
	read_def(
		std::istream &,
		Pose const &,
		FuncFactory const &
	);

	virtual void show_def( std::ostream &, Pose const & ) const;

	virtual void show( std::ostream & out ) const;

	core::Size seqpos1() const { return seqpos1_; }
	core::Size seqpos2() const { return seqpos2_; }

	void set_sequence_coupling( SequenceCouplingOP );
	SequenceCouplingOP sequence_coupling();
	SequenceCouplingCOP sequence_coupling() const;

	virtual core::Size natoms() const;
	//virtual ConstraintOP remap_resid( SequenceMapping const & ) const;
	virtual AtomID const & atom( core::Size const ) const;
	utility::vector1< AtomID > const & atom_ids() const;

	virtual void
	score(
		XYZ_Func const &,
		EnergyMap const &,
		EnergyMap &
	) const;

	virtual void
	fill_f1_f2(
		AtomID const &,
		XYZ_Func const &,
		core::Vector &,
		core::Vector &,
		EnergyMap const &
	) const;

private:
	core::Size seqpos1_;
	core::Size seqpos2_;
	SequenceCouplingOP sequence_coupling_;
	utility::vector1< AtomID > atom_ids_;
};


} // namespace constraints_additional
} // namespace protocols

#endif // INCLUDED_protocols_constraints_additional_SequenceCouplingConstraint_HH
