// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/NPDHBondEnergy.fwd.hh
/// @brief  Hydrogen bond energy method forward declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit Headers
#include <core/energy_methods/NPDHBondEnergy.hh>
#include <core/energy_methods/NPDHBondEnergyCreator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/NPDHBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>




#include <core/scoring/methods/EnergyMethodOptions.hh>





// Project headers
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers

#include <core/scoring/hbonds/HBondDatabase.hh>
#include <utility/vector1.hh>

#include <core/scoring/hbonds/constants.hh> // AUTO IWYU For MAX_R



#ifdef SERIALIZATION
// Project serialization headers
#include <core/scoring/trie/RotamerTrie.srlz.hh>

// Utility serialization headers
#include <utility/serialization/serialization.hh>

#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.energy_methods.NPDHBondEnergy" );

namespace core {
namespace energy_methods {


/// @details This must return a fresh instance of the NPDHBondEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
NPDHBondEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	return utility::pointer::make_shared< NPDHBondEnergy >( options.hbond_options() );
}

core::scoring::ScoreTypes
NPDHBondEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( npd_hbond_lr_bb );
	sts.push_back( npd_hbond_sr_bb );
	sts.push_back( npd_hbond_bb_sc );
	sts.push_back( npd_hbond_sc );
	sts.push_back( npd_hbond_intra );
	sts.push_back( npd_hbond );
	return sts;
}


/// ctor
NPDHBondEnergy::NPDHBondEnergy( core::scoring::hbonds::HBondOptions const & opts ):
	parent( utility::pointer::make_shared< NPDHBondEnergyCreator >() ),
	options_( utility::pointer::make_shared< core::scoring::hbonds::HBondOptions >( opts )),
	database_( core::scoring::hbonds::HBondDatabase::get_database(opts.params_database_tag()) )
{
}

/// copy ctor
NPDHBondEnergy::NPDHBondEnergy( NPDHBondEnergy const & src ):
	parent( src ),
	options_( utility::pointer::make_shared< core::scoring::hbonds::HBondOptions >( *src.options_ )) ,
	database_( src.database_)
{
}

NPDHBondEnergy::~NPDHBondEnergy() = default;

/// clone
core::scoring::methods::EnergyMethodOP
NPDHBondEnergy::clone() const
{
	return utility::pointer::make_shared< NPDHBondEnergy >( *this );
}

void
NPDHBondEnergy::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & sfxn ) const
{
	using core::scoring::EnergiesCacheableDataType::NPD_HBOND_SET;

	pose.update_residue_neighbors();
	core::scoring::hbonds::NPDHBondSetOP hbond_set( utility::pointer::make_shared< core::scoring::hbonds::NPDHBondSet >( *options_, pose, sfxn.weights() ) );
	pose.energies().data().set( NPD_HBOND_SET, hbond_set );
}


void
NPDHBondEnergy::setup_for_derivatives( pose::Pose & pose, core::scoring::ScoreFunction const & sfxn ) const
{
	setup_for_scoring( pose, sfxn );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
NPDHBondEnergy::residue_pair_energy(
	conformation::Residue const & ,
	conformation::Residue const & ,
	pose::Pose const & ,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const
{}

void
NPDHBondEnergy::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals
) const
{
	using core::scoring::EnergiesCacheableDataType::NPD_HBOND_SET;

	// iterate across the hbonds in the cached hbond set
	// and increment their weighted energies into the totals emap
	auto const & hbond_set = static_cast< core::scoring::hbonds::NPDHBondSet const & > ( pose.energies().data().get( NPD_HBOND_SET ) );
	for ( Size ii = 1; ii <= hbond_set.nhbonds(); ++ii ) {
		core::scoring::hbonds::HBond const & ii_hb( hbond_set.hbond( ii ) );
		if ( ii_hb.energy() > 0 ) continue;
		// std::cout << " scoring hb: " << ii_hb.don_res() << " " << ii_hb.acc_res() << " e: " << ii_hb.energy() << " wtdE: " << ii_hb.energy() * ii_hb.don_npd_weight() * ii_hb.acc_npd_weight() << std::endl;
		Real acc_don_wtd_energy = ii_hb.energy() * ii_hb.don_npd_weight() * ii_hb.acc_npd_weight();
		increment_npd_hbond_energy( ii_hb.eval_type(), totals, acc_don_wtd_energy, ii_hb.don_res() == ii_hb.acc_res() );
	}
}

void
NPDHBondEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	using core::scoring::EnergiesCacheableDataType::NPD_HBOND_SET;


	// ok -- iterate across all the hbonds for this residue and calculate
	// the per-atom derivatives
	auto const & hbond_set = static_cast< core::scoring::hbonds::NPDHBondSet const & > ( pose.energies().data().get( NPD_HBOND_SET ) );
	utility::vector1< core::scoring::hbonds::HBondCOP > res_hbonds = hbond_set.residue_hbonds( atom_id.rsd(), false );
	if ( res_hbonds.size() == 0 ) return;

	for ( Size ii = 1; ii <= res_hbonds.size(); ++ii ) {
		core::scoring::hbonds::HBond const & ii_hb( * res_hbonds[ ii ] );
		if ( ii_hb.energy() > 0 ) continue;

		conformation::Residue const & ii_don_rsd( pose.residue( ii_hb.don_res() ));
		conformation::Residue const & ii_acc_rsd( pose.residue( ii_hb.acc_res() ));
		Size const ii_hatm(ii_hb.don_hatm());
		Size const ii_datm(ii_don_rsd.atom_base(ii_hatm));
		Size const ii_aatm(ii_hb.acc_atm());

		core::scoring::hbonds::HBEvalTuple ii_hbe_type( ii_datm, ii_don_rsd, ii_aatm, ii_acc_rsd );

		core::scoring::hbonds::HBDerivAssigner assigner( *options_, ii_hbe_type, ii_don_rsd, ii_hatm, ii_acc_rsd, ii_aatm );

		core::scoring::hbonds::which_atom_in_hbond ii_which = core::scoring::hbonds::which_hb_unassigned;
		if ( atom_id.rsd() == ii_hb.don_res() ) {
			if ( atom_id.atomno() == assigner.h_ind() ) {
				ii_which = core::scoring::hbonds::which_hb_hatm;
			} else if ( atom_id.atomno() == assigner.d_ind() ) {
				ii_which = core::scoring::hbonds::which_hb_datm;
			}
		} else {
			if ( atom_id.atomno() == assigner.a_ind() ) {
				ii_which = core::scoring::hbonds::which_hb_aatm;
			} else if ( atom_id.atomno() == assigner.abase_ind() ) {
				ii_which = core::scoring::hbonds::which_hb_abase;
			} else if ( atom_id.atomno() == assigner.abase_prime_ind() ) {
				ii_which = core::scoring::hbonds::which_hb_abase_prime;
			} else if ( atom_id.atomno() == assigner.abase2_ind() ) {
				ii_which = core::scoring::hbonds::which_hb_abase2;
			}
		}

		if ( ii_which == core::scoring::hbonds::which_hb_unassigned ) continue;

		Real unweighted_energy( 0.0 ); // we already know the energy, actually.
		core::scoring::hbonds::HBondDerivs deriv;
		hb_energy_deriv( *database_, *options_, ii_hbe_type,
			ii_don_rsd.atom( ii_datm ).xyz(),
			ii_don_rsd.atom( ii_hatm ).xyz(),
			ii_acc_rsd.atom( ii_aatm ).xyz(),
			ii_acc_rsd.atom( assigner.abase_ind() ).xyz(),
			ii_acc_rsd.atom( assigner.abase2_ind() ).xyz(),
			unweighted_energy, true /*eval deriv*/, deriv );

		// ok; now we need the derivative for just this atom
		core::scoring::hbonds::AssignmentScaleAndDerivVectID ii_asadvi = assigner.assignment( ii_which );
		core::scoring::DerivVectorPair ii_raw_deriv = ii_asadvi.scale_ * deriv.deriv( ii_asadvi.dvect_id_ ); // aka dE/dxyz

		core::scoring::DerivVectorPair ii_wtd_deriv;

		Real sfxn_weight = npd_hb_eval_type_weight( ii_hbe_type.eval_type(), weights );

		// note that dEtot_dEhb is actually dEtot_dEhb divided by the sfxnwt
		// so that we have to multiply the score-function weight into it twice
		if ( atom_id.rsd() == ii_hb.don_res() ) {
			ii_wtd_deriv = sfxn_weight * (
				ii_hb.don_npd_weight() * ii_hb.acc_npd_weight() +
				sfxn_weight *  hbond_set.dEtot_dEhb( ii_hb.index() ) ) * ii_raw_deriv;
		} else {
			ii_wtd_deriv = sfxn_weight * (
				ii_hb.don_npd_weight() * ii_hb.acc_npd_weight() +
				sfxn_weight * hbond_set.dEtot_dEhb( ii_hb.index() ) ) * ii_raw_deriv;
		}

		F1 += ii_wtd_deriv.f1();
		F2 += ii_wtd_deriv.f2();

		//// the way I had calculated the energies when the sfxn weights were not part
		//// of the calculation for the donor and acceptor weights
		//if ( atom_id.rsd() == ii_hb.don_res() ) {
		// ii_wtd_deriv = (
		//  ii_hb.don_npd_weight() * ii_hb.acc_npd_weight() +
		//  hbond_set.dEtot_dEhb( ii_hb.index() ) ) * ii_raw_deriv;
		//} else {
		// ii_wtd_deriv = (
		//  ii_hb.don_npd_weight() * ii_hb.acc_npd_weight() +
		//  hbond_set.dEtot_dEhb( ii_hb.index() ) ) * ii_raw_deriv;
		//}
		//
		//Real sfxn_weight = npd_hb_eval_type_weight( ii_hbe_type.eval_type(), weights );
		//
		//F1 += sfxn_weight * ii_wtd_deriv.f1();
		//F2 += sfxn_weight * ii_wtd_deriv.f2();


	}
}



bool
NPDHBondEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const { return true; }


/// @brief HACK!  core::scoring::hbonds::MAX_R defines the maximum donorH to acceptor distance.
// The atomic_interaction_cutoff method is meant to return the maximum distance
// between two *heavy atoms* for them to have a zero interaction energy.
// I am currently assuming a 1.35 A maximum distance between a hydrogen and the
// heavy atom it is bound to, stealing this number from the CYS.params file since
// the HG in CYS is much further from it's SG than aliphatic hydrogens are from their carbons.
// This is a bad idea.  Someone come up with a way to fix this!
//
// At 4.35 A interaction cutoff, the hbond energy function is incredibly short ranged!
Distance
NPDHBondEnergy::atomic_interaction_cutoff() const
{
	return core::scoring::hbonds::MAX_R + 1.35; // MAGIC NUMBER
}

/// @brief the atomic interaction cutoff and the hydrogen interaction cutoff are the same.
Real
NPDHBondEnergy::hydrogen_interaction_cutoff2() const
{
	return (core::scoring::hbonds::MAX_R + 1.35) * ( core::scoring::hbonds::MAX_R + 1.35 );
}


/// @brief NPDHBondEnergy is context sensitive
void
NPDHBondEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & context_graphs_required
) const
{
	context_graphs_required[ core::scoring::ten_A_neighbor_graph ] = true;
}

bool
NPDHBondEnergy::defines_intrares_energy( core::scoring::EnergyMap const & /*weights*/ ) const
{
	return false;
}


void
NPDHBondEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const
{}

/// @details
core::Size
NPDHBondEnergy::version() const
{
	return 1; // Initial versioning
}

} // scoring
} // core

#ifdef    SERIALIZATION

//CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_hbonds_NPDHBondEnergy )

#endif
