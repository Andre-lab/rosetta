// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// Unit headers
#include <core/scoring/constraints/Constraints.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/XYZ_Func.hh>

// AUTO-REMOVED #include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.fwd.hh>

// Utility Headers

// C++ Headers
#include <algorithm>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace core {
namespace scoring {
namespace constraints {

Constraints::Constraints()
{}

Constraints::Constraints( Constraints const & other ) :
	ReferenceCount()
{
	copy_from( other );
}



ConstraintsOP
Constraints::clone() const
{

	ConstraintsOP new_constraints = new Constraints();
	new_constraints->copy_from( *this );
	return new_constraints;
}

Constraints const &
Constraints::operator = ( Constraints const & rhs )
{
	copy_from( rhs );
	return *this;
}


void
Constraints::copy_from( Constraints const & other ) {
	for ( ConstraintCOPs::const_iterator
			iter = other.constraints_.begin(), iter_end = other.constraints_.end();
			iter != iter_end; ++iter ) {
		// All Constraints are now immutable, and so do not ever need to be cloned.
		//constraints_.push_back( (*iter)->clone() );
		constraints_.push_back( (*iter) );
	}
}

// call the setup_for_derivatives for each constraint
void
Constraints::setup_for_scoring( XYZ_Func const & xyz, ScoreFunction const &scfxn ) const {
	for ( ConstraintCOPs::const_iterator it= constraints_.begin(), ite=constraints_.end();
				it != ite; ++it ) {
		(*it)->setup_for_scoring( xyz, scfxn );
	}
}

// call the setup_for_derivatives for each constraint
void
Constraints::setup_for_derivatives( XYZ_Func const & xyz, ScoreFunction const &scfxn ) const {
	for ( ConstraintCOPs::const_iterator it= constraints_.begin(), ite=constraints_.end();
				it != ite; ++it ) {
		(*it)->setup_for_derivatives( xyz, scfxn );
	}
}


void
Constraints::eval_intrares_atom_derivative(
	id::AtomID const & atom_id,
	conformation::Residue const & residue,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	ResidueXYZ resxyz( residue );
	for ( ConstraintCOPs::const_iterator it= constraints_.begin(), ite=constraints_.end();
				it != ite; ++it ) {
		Constraint const & cst( **it );
		Vector f1(0.0), f2(0.0);
		cst.fill_f1_f2( atom_id, resxyz, f1, f2, weights );
		F1 += f1;
		F2 += f2;
	}
}

void
Constraints::eval_respair_atom_derivative(
	id::AtomID const & atom_id,
	conformation::Residue const & residue1,
	conformation::Residue const & residue2,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	ResiduePairXYZ respairxyz( residue1, residue2 );
	for ( ConstraintCOPs::const_iterator it= constraints_.begin(), ite=constraints_.end();
			it != ite; ++it ) {
		Constraint const & cst( **it );
		Vector f1(0.0), f2(0.0);
		cst.fill_f1_f2( atom_id, respairxyz, f1, f2, weights );
		F1 += f1;
		F2 += f2;
	}
}

void
Constraints::eval_ws_atom_derivative(
	id::AtomID const & atom_id,
	conformation::Conformation const & conformation,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	ConformationXYZ confxyz( conformation );
	for ( ConstraintCOPs::const_iterator it= constraints_.begin(), ite=constraints_.end();
			it != ite; ++it ) {
		Constraint const & cst( **it );
		Vector f1(0.0), f2(0.0);
		cst.fill_f1_f2( atom_id, confxyz, f1, f2, weights );
		F1 += f1;
		F2 += f2;
	}
}



/// private
/// does not zero the emap entries before accumulating
///
void
Constraints::energy(
	XYZ_Func const & xyz_func,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	for ( ConstraintCOPs::const_iterator it=constraints_.begin(), ite = constraints_.end(); it != ite; ++it ) {
		Constraint const & cst( **it );
		cst.score( xyz_func, weights, emap );
//		cst.show( std::cout );
	}
}

/// will fail if Residues dont contain all the necessary atoms
void
Constraints::residue_pair_energy(
	Residue const & rsd1,
	Residue const & rsd2,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	ResiduePairXYZ const xyz_func( rsd1, rsd2 );
	energy( xyz_func, weights, emap );
}

/// will fail if Residues dont contain all the necessary atoms
void
Constraints::intra_residue_energy(
	Residue const & rsd,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	ResidueXYZ const xyz_func( rsd );
	energy( xyz_func, weights, emap );
}

///
void
Constraints::conformation_energy(
	Conformation const & conformation,
	EnergyMap const & weights,
	EnergyMap & emap
) const
{
	ConformationXYZ const xyz_func( conformation );
	energy( xyz_func, weights, emap );
}

///
void
Constraints::add_constraint( ConstraintCOP cst )
{
	constraints_.push_back( cst );
}

Constraints::const_iterator Constraints::begin() const { return constraints_.begin(); }
Constraints::const_iterator Constraints::end() const { return constraints_.end(); }


/// @details If this list contains the same constraint multiple times,
/// only one copy of it is removed.
/// I can't imagine *why* that scenario would ever come up, but still...
/// flo jan '11 added object_comparison bool that triggers removal of
/// constraint if the actual constraints are identical, even though
/// they're different objects
bool
Constraints::remove_constraint(
	ConstraintCOP cst,
	bool object_comparison
)
{

	if( object_comparison ){
		for( ConstraintCOPs::iterator cst_it = constraints_.begin(), cst_end = constraints_.end(); cst_it != cst_end; ++cst_it ){

			if( *cst == **cst_it ){
				constraints_.erase( cst_it );
				return true;
			}
		}
		return false;
	}

	ConstraintCOPs::iterator where = std::find( constraints_.begin(), constraints_.end(), cst );
	if( where == constraints_.end() ) return false;
	constraints_.erase( where );
	return true;
}

void
Constraints::show( std::ostream & out ) {
	for ( ConstraintCOPs::const_iterator it=constraints_.begin(), ite = constraints_.end();
				it != ite;
				++it ) {
		Constraint const & cst( **it );
		cst.show( out );
	}
}

void
Constraints::show_definition(
	std::ostream & out,
	pose::Pose const & pose
) const {
	for ( ConstraintCOPs::const_iterator it=constraints_.begin(), ite = constraints_.end();
				it != ite;
				++it ) {
		Constraint const & cst( **it );
		cst.show_def( out, pose );
	}
}

Size
Constraints::show_violations(
	std::ostream& out,
	pose::Pose const & pose,
	Size verbose_level,
	Real threshold
) {
	Size total_viol( 0 );
	Size total_cst( 0 );
	for ( ConstraintCOPs::const_iterator it=constraints_.begin(), ite = constraints_.end();
				it != ite;
				++it ) {
		Constraint const & cst( **it );
		total_viol+=cst.show_violations( out, pose, verbose_level, threshold );
		total_cst++;
	}
 	if ( verbose_level > 60 ) out << " of total: " << total_cst << " ";
	return total_viol;
}

Size
Constraints::size() const
{
	return constraints_.size();
}

void
Constraints::clear()
{
	constraints_.clear();
}

} // constraints
} // scoring
} // core
