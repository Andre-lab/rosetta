// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/PeakAssignment.hh>
#include <protocols/noesy_assign/CrossPeak.hh>

#include <protocols/noesy_assign/ResonanceList.hh>

// Package Headers
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Project Headers
#include <protocols/noesy_assign/util.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <basic/options/option_macros.hh>


//// C++ headers
#include <cstdlib>
#include <string>

static basic::Tracer tr("protocols.noesy_assign.assignment");

namespace protocols {
namespace noesy_assign {
using namespace core;

PeakAssignment::PeakAssignment( CrossPeakAP const& cp, core::Size assign_spin1, core::Size assign_spin2 )
  : crosspeak_( cp ),
    spin_assign_index1_( assign_spin1 ),
    spin_assign_index2_( assign_spin2 ),
		chemshift_overlap_( 1.0 ),
		symmetry_compliance_( false ),
		covalent_compliance_( false ),
		decoy_compatibility_( 1.0 ),
		network_anchoring_( 1.0 ),
		network_anchoring_per_residue_( 200 ) //so we are not eliminated by default
{
	if ( cp )  update_resonances_from_peak();
}

ResonanceList const& PeakAssignment::resonances() const {
	return crosspeak_->resonances();
}

void PeakAssignment::dump_weights( std::ostream& os ) const {
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	using namespace ObjexxFCL::fmt;
	os << F( 5, 3, chemshift_overlap_ ) << " "
		 << F( 5, 3, symmetry_compliance_ ? params.symmetry_compliance_weight_ : 1.0 ) << " "
		 << F( 5, 3, covalent_compliance_ ? params.covalent_compliance_weight_ : 1.0 ) << " "
		 << F( 5, 3, decoy_compatibility_ ) << " "
		 << F( 5, 3, network_anchoring_ ) << " "
		 << F( 5, 3, network_anchoring_per_residue_ ) << " ";
}

Real PeakAssignment::peak_volume() const {
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	return chemshift_overlap_
		* std::min( params.smax_,
			( symmetry_compliance_ ? params.symmetry_compliance_weight_ : 1.0 )
			* ( covalent_compliance_ ? params.covalent_compliance_weight_ : 1.0 )
			* network_anchoring_ )
		* decoy_compatibility_;
}

Real PeakAssignment::normalized_peak_volume() const {
	if ( crosspeak_->cumulative_peak_volume() > 0.0 ) return peak_volume() / crosspeak_->cumulative_peak_volume();
	return 0.0;
}

void PeakAssignment::update_resonances_from_peak() {
  resonance1_ = crosspeak_->proton( 1 ).assignment( spin_id( 1 ) );
  resonance2_ = crosspeak_->proton( 2 ).assignment( spin_id( 2 ) );
}

void PeakAssignment::update_chemshiftscore_from_peak() {
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	Real const& weight( params.chemshift_overlap_weight_ );
	Real sum( 0.0 );
	for ( Size d=1; d<=crosspeak_->dimension(); d++ ) {
		CrossPeak::Spin const& spin( crosspeak_->spin( d ) );
		Resonance const& assigned_resonance( resonances()[ spin.assignment( spin_id( d>2 ? d-2 : d ) ) ] );
		Real diff( spin.freq()-assigned_resonance.freq() );
		Real s( diff/weight/std::max( crosspeak_->tolerance( d ), assigned_resonance.tolerance() ) );
		sum += s*s;
	}
	chemshift_overlap_ = exp( -0.5*sum );
}

void PeakAssignment::update_upperdistance_score(/*dmax*/ ) {
	core::id::NamedAtomID const& atom1( atom( 1 ) );
	core::id::NamedAtomID const& atom2( atom( 2 ) );
	covalent_compliance_ = covalent_compliance( atom1, atom2 );
}

PeakAssignment::NmrConstraintOP PeakAssignment::create_constraint(
					pose::Pose const& pose,
					core::scoring::constraints::FuncOP func
) const {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_PA_GEN_CST );

	bool flip = atom( 1 ).rsd() > atom( 2 ).rsd();
	core::id::NamedAtomID const& atom1( atom( flip ? 2 : 1 ) );
	core::id::NamedAtomID const& atom2( atom( flip ? 1 : 2 ) );

	tr.Debug << "create constraint for atom: " << atom1 << " " << atom2 << " from res_id: "
					 << resonance_id( 1 ) << " " << resonance_id( 2 ) << std::endl;

	using namespace core::scoring::constraints;
	if ( !func )
		func = new BoundFunc( 1.5,
		5.5,
		1.0,
			"VC "+ObjexxFCL::string_of( normalized_peak_volume(), 3 )
	);
	return new NmrConstraint( atom1, atom2, pose, func ); //later figure out what the func should be...
}


	///@brief return resonance_id, i.e., pointer into Resonance list that will resolve in assigned atom
core::Size PeakAssignment::label_resonance_id( core::Size select ) const {
  runtime_assert( select == 1 || select == 2 );
  runtime_assert( crosspeak_ );
  runtime_assert( crosspeak_->has_label( select ) );
  return crosspeak_->label( select ).assignment( spin_id( select ) );
}

PeakAssignment const BOGUS_ASSIGNMENT( NULL, 0, 0 );


// void PeakAssignment::invalidate_assignment() {
//   if ( crosspeak_ ) crosspeak_->invalidate_assignment( assignment_index_ );
// }

}
}
