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
/// @author


// Unit headers
#include <core/scoring/hbonds/hbonds.hh>

// Package headers
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>

// // Project headers
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

//pba membrane specific hbond
#include <core/scoring/Membrane_FAPotential.hh>
#include <core/scoring/Membrane_FAPotential.fwd.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>

//Auto Headers
#include <ObjexxFCL/FArray2D.hh>


//#include <core/scoring/Energies.hh>

// // Numeric headers
// #include <numeric/numeric.functions.hh>

using namespace ObjexxFCL;
//pba
using namespace basic::options;
using namespace OptionKeys;

namespace core {
namespace scoring {
namespace hbonds {

static basic::Tracer tr("core.scoring.hbonds.hbonds");

/**
	This routine fills an hbond-set with hbonds. All hbonds are included,
	even ones which might be excluded later based on the backbone-hbond
	exclusion.

	WARNING WARNING WARNING
	The pose must have an update energies object, eg it must be scored.
	WARNING WARNING WARNING
**/

void
fill_hbond_set(
	pose::Pose const & pose,
	bool const calculate_derivative,
	HBondSet & hbond_set,
	bool const exclude_bb  /* default false */,
	bool const exclude_bsc /* default false */,
	bool const exclude_scb /* default false */,
	bool const exclude_sc  /* default false */
)
{
	assert( pose.energies().residue_neighbors_updated() );

	// clear old data
	hbond_set.clear();
	HBondDatabase const & database( * HBondDatabase::get_database(hbond_set.hbond_options().params_database_tag()));

	// need to know which residues are neighbors
	// and what the neighbor-numbers are for each residue since some of the
	// weights are environment-dependent.
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

	// loop over all nbr-pairs
	for ( Size res1 = 1; res1 <= pose.total_residue(); ++res1 ) {
		int const nb1 = tenA_neighbor_graph.get_node( res1 )->num_neighbors_counting_self_static();
		conformation::Residue const & rsd1( pose.residue( res1 ) );

		for ( graph::Graph::EdgeListConstIter
				iru = energy_graph.get_node(res1)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(res1)->const_upper_edge_list_end();
				iru != irue; ++iru ) {

			int const res2( (*iru)->get_second_node_ind() );

			conformation::Residue const & rsd2( pose.residue( res2 ) );
			if ( hbond_set.hbond_options().exclude_DNA_DNA() && rsd1.is_DNA() && rsd2.is_DNA() ) continue;

			int const nb2 = tenA_neighbor_graph.get_node( res2 )->num_neighbors_counting_self_static();

			//pba membrane specific hbond
			if ( hbond_set.hbond_options().Mbhbond() ) {
				identify_hbonds_1way_membrane(
					database,
					rsd1, rsd2, nb1, nb2, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set, pose);

				identify_hbonds_1way_membrane(
					database,
					rsd2, rsd1, nb2, nb1, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set, pose);
			} else {
			   identify_hbonds_1way(
					database,
					rsd1, rsd2, nb1, nb2, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set);

			   identify_hbonds_1way(
					database,
					rsd2, rsd1, nb2, nb1, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set);

			}
		} // nbrs of res1
		if(!hbond_set.hbond_options().exclude_self_hbonds() && hbond_set.hbond_options().exclude_DNA_DNA() && rsd1.is_DNA() ){
			//pba membrane specific hbond
			if ( hbond_set.hbond_options().Mbhbond() ) {
				identify_hbonds_1way_membrane(
					database,
					rsd1, rsd1, nb1, nb1, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set, pose);
			} else {
			   identify_hbonds_1way(
					database,
					rsd1, rsd1, nb1, nb1, calculate_derivative,
					exclude_bb, exclude_bsc, exclude_scb, exclude_sc, hbond_set);
			}
		}
	} // res1
}

/// @brief  Get the f1 and f2 contributions from all hbonds involving this atom
/*void
get_atom_hbond_derivative(
	id::AtomID const & atom,
	HBondSet const & hbond_set,
	EnergyMap const & weights,
	Vector & f1,
	Vector & f2
){
	f1 = Vector(0.0);
	f2 = Vector(0.0);

	utility::vector1< HBondCOP > const & hbonds
		( hbond_set.atom_hbonds( atom ) );

	for ( Size i=1; i<= hbonds.size(); ++i ) {
		HBond const & hbond( *hbonds[ i ] );
		Real sign_factor( 0.0 );
		if ( hbond.atom_is_donorH( atom ) ) sign_factor = 1.0;
		else {
			assert( hbond.atom_is_acceptor( atom ) );
			sign_factor = -1;
		}
		// get the appropriate type of hbond weight
		Real const weight(sign_factor * hbond.weight() * hb_eval_type_weight(hbond.eval_type(), weights));
		f1 += weight * hbond.deriv().first;
		f2 += weight * hbond.deriv().second;

	}
}*/


/// @details identify_hbonds_1way is overloaded to either add HBond objects to
/// an HBondSet or to accumulate energy into a EnergyMap
/// object.  This is done for performance reasons.  The allocation of
/// the temporary HBondSet on the heap causes a substatial slow down.
void
identify_hbonds_1way(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	// output
	HBondSet & hbond_set
)
{
	assert( don_rsd.seqpos() != acc_rsd.seqpos() );

	// <f1,f2> -- derivative vectors
	//std::pair< Vector, Vector > deriv( Vector(0.0), Vector(0.0 ) );
	HBondDerivs derivs;

	for ( chemical::AtomIndices::const_iterator
			hnum  = don_rsd.Hpos_polar().begin(),	hnume = don_rsd.Hpos_polar().end();
			hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		if (datm_is_bb){
			if (exclude_bb && exclude_scb) continue;
		} else {
			if (exclude_sc && exclude_bsc) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( chemical::AtomIndices::const_iterator
				anum  = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end();
				anum != anume; ++anum ) {
			Size const aatm( *anum );
			if (acc_rsd.atom_is_backbone(aatm)){
				if (datm_is_bb){
					if (exclude_bb) continue;
				} else {
					if (exclude_bsc) continue;
				}
			} else {
				if (datm_is_bb){
					if (exclude_scb) continue;
				} else {
					if (exclude_sc) continue;
				}
			}

			// rough filter for existance of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalType hbe_type( hbond_evaluation_type( datm, don_rsd, aatm, acc_rsd));

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			assert( base2 > 0 && base != base2 );


			hb_energy_deriv( database, hbond_set.hbond_options(),
				hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if (unweighted_energy >= MAX_HB_ENERGY) continue;

			Real environmental_weight
				(!hbond_set.hbond_options().use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, hbond_set.hbond_options()));

			//////
			// now we have identified a hbond -> append it into the hbond_set
			hbond_set.append_hbond( hatm, don_rsd, aatm, acc_rsd,
				hbe_type, unweighted_energy, environmental_weight, derivs );

			//////

		} // loop over acceptors
	} // loop over donors
}

void
identify_hbonds_1way(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	HBondOptions const & options,
	// output
	EnergyMap & emap
)
{
	assert( don_rsd.seqpos() != acc_rsd.seqpos() );

	// <f1,f2> -- derivative vectors
	//std::pair< Vector, Vector > deriv( Vector(0.0), Vector(0.0 ) );
	HBondDerivs derivs;

	for ( chemical::AtomIndices::const_iterator
			hnum = don_rsd.Hpos_polar().begin(), hnume = don_rsd.Hpos_polar().end();
			hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		if (datm_is_bb){
			if (exclude_bb && exclude_scb) continue;
		} else {
			if (exclude_sc && exclude_bsc) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( chemical::AtomIndices::const_iterator
				anum = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end();
				anum != anume; ++anum ) {
			Size const aatm( *anum );
			if (acc_rsd.atom_is_backbone(aatm)){
				if (datm_is_bb){
					if (exclude_bb) continue;
				} else {
					if (exclude_bsc) continue;
				}
			} else {
				if (datm_is_bb){
					if (exclude_scb) continue;
				} else {
					if (exclude_sc) continue;
				}
			}

			// rough filter for existance of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalType hbe_type( hbond_evaluation_type( datm, don_rsd, aatm, acc_rsd));

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			assert( base2 > 0 && base != base2 );

			hb_energy_deriv( database, options, hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if (unweighted_energy >= MAX_HB_ENERGY) continue;

			Real environmental_weight
				(!options.use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, options));

			////////
			// now we have identified an hbond -> accumulate its energy

			Real hbE = unweighted_energy /*raw energy*/ * environmental_weight /*env-dep-wt*/;
			switch(get_hbond_weight_type(hbe_type)){
			case hbw_NONE:
			case hbw_SR_BB:
				emap[hbond_sr_bb] += hbE; break;
			case hbw_LR_BB:
				emap[hbond_lr_bb] += hbE; break;
			case hbw_SR_BB_SC:
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_sr_bb_sc] += hbE; break;
			case hbw_LR_BB_SC:
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_lr_bb_sc] += hbE; break;
			case hbw_SC:
				emap[hbond_sc] += hbE; break;
			default:
				tr << "Warning: energy from unexpected HB type ignored "
					<< hbe_type << std::endl;
				runtime_assert(false);
				break;
			}
			/////////

		} // loop over acceptors
	} // loop over donors
}

void
identify_hbonds_1way_membrane(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	// output
	HBondSet & hbond_set,
	pose::Pose const & pose
){
	assert( don_rsd.seqpos() != acc_rsd.seqpos() );

        //pbadebug
        //std::cout << "entered in identify_hbonds_1way_membrane() " << std::endl;

	// <f1,f2> -- derivative vectors
	//std::pair< Vector, Vector > deriv( Vector(0.0), Vector(0.0 ) );
        HBondDerivs derivs;

	for ( chemical::AtomIndices::const_iterator
			hnum = don_rsd.Hpos_polar().begin(), hnume = don_rsd.Hpos_polar().end();
			hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		if (datm_is_bb){
			if (exclude_bb && exclude_scb) continue;
		} else {
			if (exclude_sc && exclude_bsc) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( chemical::AtomIndices::const_iterator
			anum = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end();
			anum != anume; ++anum ) {
			Size const aatm( *anum );
			if (acc_rsd.atom_is_backbone(aatm)){
				if (datm_is_bb){
					if (exclude_bb) continue;
				} else {
					if (exclude_bsc) continue;
				}
			} else {
				if (datm_is_bb){
					if (exclude_scb) continue;
				} else {
					if (exclude_sc) continue;
				}
			}

			// rough filter for existance of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalType hbe_type( hbond_evaluation_type( datm, don_rsd, aatm, acc_rsd));

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			assert( base2 > 0 && base != base2 );


			hb_energy_deriv( database, hbond_set.hbond_options(),
				hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if (unweighted_energy >= MAX_HB_ENERGY) continue;

			Real environmental_weight
				(!hbond_set.hbond_options().use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, hbond_set.hbond_options()));

			//pba membrane depth dependent correction to the environmental_weight 

			Real membrane_depth_dependent_weight(
				get_membrane_depth_dependent_weight(pose, don_nb, acc_nb, hatm_xyz,
				acc_rsd.atom(aatm ).xyz()));

			/*if (acc_rsd.seqpos() == 3 && don_rsd.seqpos() == 6) {
				std::cout << "in 1: acc " << acc_rsd.seqpos() << " don " << don_rsd.seqpos() << " mbhb_weight "
				 << membrane_depth_dependent_weight << std::endl;
			}*/

			environmental_weight = membrane_depth_dependent_weight;
			//////
			// now we have identified a hbond -> put it into the hbond_set

			hbond_set.append_hbond( hatm, don_rsd, aatm, acc_rsd,
				hbe_type, unweighted_energy, environmental_weight, derivs );

			//////

		} // loop over acceptors
	} // loop over donors
}

void
identify_hbonds_1way_membrane(
	HBondDatabase const & database,
	conformation::Residue const & don_rsd,
	conformation::Residue const & acc_rsd,
	Size const don_nb,
	Size const acc_nb,
	bool const evaluate_derivative,
	bool const exclude_bb,  /* exclude if acc=bb and don=bb */
	bool const exclude_bsc, /* exclude if acc=bb and don=sc */
	bool const exclude_scb, /* exclude if acc=sc and don=bb */
	bool const exclude_sc,  /* exclude if acc=sc and don=sc */
	HBondOptions const & options,
	// output
	EnergyMap & emap,
	pose::Pose const & pose
)
{
	assert( don_rsd.seqpos() != acc_rsd.seqpos() );

	// <f1,f2> -- derivative vectors
	//std::pair< Vector, Vector > deriv( Vector(0.0), Vector(0.0 ) );
        HBondDerivs derivs;

	for ( chemical::AtomIndices::const_iterator
		hnum = don_rsd.Hpos_polar().begin(), hnume = don_rsd.Hpos_polar().end();
		hnum != hnume; ++hnum ) {
		Size const hatm( *hnum );
		Size const datm(don_rsd.atom_base(hatm));
		bool datm_is_bb = don_rsd.atom_is_backbone(datm);

		if (datm_is_bb){
			if (exclude_bb && exclude_scb) continue;
		} else {
			if (exclude_sc && exclude_bsc) continue;
		}
		Vector const & hatm_xyz(don_rsd.atom(hatm).xyz());
		Vector const & datm_xyz(don_rsd.atom(datm).xyz());

		for ( chemical::AtomIndices::const_iterator
			anum = acc_rsd.accpt_pos().begin(), anume = acc_rsd.accpt_pos().end();
			anum != anume; ++anum ) {
			Size const aatm( *anum );
			if (acc_rsd.atom_is_backbone(aatm)){
				if (datm_is_bb){
					if (exclude_bb) continue;
				} else {
					if (exclude_bsc) continue;
				}
			} else {
				if (datm_is_bb){
					if (exclude_scb) continue;
				} else {
					if (exclude_sc) continue;
				}
			}

			// rough filter for existance of hydrogen bond
			if ( hatm_xyz.distance_squared( acc_rsd.xyz( aatm ) ) > MAX_R2 ) continue;

			Real unweighted_energy( 0.0 );

			HBEvalType hbe_type( hbond_evaluation_type( datm, don_rsd, aatm, acc_rsd));

			int const base ( acc_rsd.atom_base( aatm ) );
			int const base2( acc_rsd.abase2( aatm ) );
			assert( base2 > 0 && base != base2 );

			hb_energy_deriv( database, options,
				hbe_type, datm_xyz, hatm_xyz,
				acc_rsd.atom(aatm ).xyz(),
				acc_rsd.atom(base ).xyz(),
				acc_rsd.atom(base2).xyz(),
				unweighted_energy, evaluate_derivative, derivs);

			if (unweighted_energy >= MAX_HB_ENERGY) continue;

			Real environmental_weight
				(!options.use_hb_env_dep() ? 1 :
				get_environment_dependent_weight(hbe_type, don_nb, acc_nb, options));

			//pba membrane depth dependent weight

			Real membrane_depth_dependent_weight(
				get_membrane_depth_dependent_weight(pose, don_nb, acc_nb, hatm_xyz,
				acc_rsd.atom(aatm ).xyz()));

			/*if ((acc_rsd.seqpos() == 7 && don_rsd.seqpos() == 113) || (acc_rsd.seqpos() == 26 
				&& don_rsd.seqpos() == 83) || (acc_rsd.seqpos() == 50 && don_rsd.seqpos() == 58)) {
				std::cout << "in 2: acc " << acc_rsd.seqpos() << " don " << don_rsd.seqpos() << " mbhb_weight " 
				 << membrane_depth_dependent_weight << std::endl;
			}*/

			environmental_weight = membrane_depth_dependent_weight;

			////////
			// now we have identified an hbond -> accumulate its energy

			Real hbE = unweighted_energy /*raw energy*/ * environmental_weight /*env-dep-wt*/;
			switch(get_hbond_weight_type(hbe_type)){
			case hbw_NONE:
			case hbw_SR_BB:
				emap[hbond_sr_bb] += hbE; break;
			case hbw_LR_BB:
				emap[hbond_lr_bb] += hbE; break;
			case hbw_SR_BB_SC:
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_sr_bb_sc] += hbE; break;
			case hbw_LR_BB_SC:
				//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
				emap[hbond_bb_sc] += hbE;
				emap[hbond_lr_bb_sc] += hbE; break;
			case hbw_SC:
				emap[hbond_sc] += hbE; break;
			default:
				tr << "Warning: energy from unexpected HB type ignored "
					<< hbe_type << std::endl;
				runtime_assert(false);
				break;
			}
			/////////

		} // loop over acceptors
	} // loop over donors
}



//mjo this should be the only way to assign hbond energies.  If you
//feel the need to collect the energies from some of the bonds,
//consider filtering the hbond set instead.  Don't assign energies for
//"special cases" with other functions, add cases to the HBEvalType
//instead.

void
get_hbond_energies(
  HBondSet const & hbond_set,
  EnergyMap & emap
)
{
	for (Size i = 1; i <= hbond_set.nhbonds(); ++i){
		if (!hbond_set.allow_hbond(i)) continue;

		HBond const & hbond(hbond_set.hbond(i));
		HBEvalType const hbe_type = hbond.eval_type();

		Real hbE = hbond.energy() /*raw energy*/ * hbond.weight() /*env-dep-wt*/;

		switch(get_hbond_weight_type(hbe_type)){
		case hbw_NONE:
		case hbw_SR_BB:
			emap[hbond_sr_bb] += hbE; break;
		case hbw_LR_BB:
			emap[hbond_lr_bb] += hbE; break;
		case hbw_SR_BB_SC:
			//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
			emap[hbond_bb_sc] += hbE;
			emap[hbond_sr_bb_sc] += hbE; break;
		case hbw_LR_BB_SC:
			//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
			emap[hbond_bb_sc] += hbE;
			emap[hbond_lr_bb_sc] += hbE; break;
		case hbw_SC:
			emap[hbond_sc] += hbE; break;
		default:
			tr << "Warning: energy from unexpected HB type ignored "
				<< hbond << std::endl;
			runtime_assert(false);
			break;
		}
	}
}

// this overloads get_hbond_energies to take an EnergyMap because
// EnergyMap and EnergyMap cannot be substituted efficiently
/*
void
get_hbond_energies(
	HBondSet const & hbond_set,
	EnergyMap & emap)
{
	for (Size i = 1; i <= hbond_set.nhbonds(); ++i){
		if (!hbond_set.allow_hbond(i)) continue;

		HBond const & hbond(hbond_set.hbond(i));
		HBEvalType const hbe_type = hbond.eval_type();
		// raw energy * env-dep-wt
		Real hbE = hbond.energy() * hbond.weight();

		switch(get_hbond_weight_type(hbe_type)){
		case hbw_NONE:
		case hbw_SR_BB:
			emap[hbond_sr_bb] += hbE; break;
		case hbw_LR_BB:
			emap[hbond_lr_bb] += hbE; break;
		case hbw_SR_BB_SC:
			//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
			emap[hbond_bb_sc] += hbE;
			emap[hbond_sr_bb_sc] += hbE; break;
		case hbw_LR_BB_SC:
			//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
			emap[hbond_bb_sc] += hbE;
			emap[hbond_lr_bb_sc] += hbE; break;
		case hbw_SC:
			emap[hbond_sc] += hbE; break;
		default:
			tr << "Warning: energy from unexpected HB type ignored "
				<< hbond << std::endl;
			runtime_assert(false);
			break;
		}
	}
}
*/


Real
hb_eval_type_weight(
	HBEvalType const & hbe_type,
	EnergyMap const & emap)
{
	Real weight(0.0);
	switch(get_hbond_weight_type(hbe_type)){
	case hbw_NONE:
	case hbw_SR_BB:
		weight += emap[hbond_sr_bb]; break;
	case hbw_LR_BB:
		weight += emap[hbond_lr_bb]; break;
	case hbw_SR_BB_SC:
		//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
		weight += emap[hbond_bb_sc];
		weight += emap[hbond_sr_bb_sc]; break;
	case hbw_LR_BB_SC:
		//Note this is double counting if both hbond_bb_sc and hbond_sr_bb_sc have nonzero weight!
		weight += emap[hbond_bb_sc];
		weight += emap[hbond_lr_bb_sc]; break;
	case hbw_SC:
		weight += emap[hbond_sc]; break;
	default:
		tr << "Warning: Unexpected HBondWeightType " << get_hbond_weight_type(hbe_type) << std::endl;
		runtime_assert(false);
		break;
	}
	return weight;
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// functions for environment-dependent weighting
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

////////////////////////
// lin jiang's approach:

//JSS note: this is half of the weight from one atom;
// the burial weight is the sum from donor and acceptor.
inline core::Real
burial_weight(int const nb)
{
	if ( nb < 7 ) return 0.1;
	if ( nb > 24 ) return 0.5;
	return (nb-2.75)*(0.5/21.25);
}

core::Real
hb_env_dep_burial_lin(int const nb1, int const nb2)
{
	return (burial_weight(nb1) + burial_weight(nb2));
}

////////////////////////////
// tanja kortemme's approach

void
burial3class_weight_initializer( FArray2D_double & burial )
{
	burial( 1, 1 ) = 0.2 ; burial( 1, 2 ) = 0.2 ; burial( 1, 3 ) = 0.55;
	burial( 2, 1 ) = 0.2 ; burial( 2, 2 ) = 0.55; burial( 2, 3 ) = 1.0 ;
	burial( 3, 1 ) = 0.55; burial( 3, 2 ) = 1.0 ; burial( 3, 3 ) = 1.0 ;
}

inline int
get_burial_3(
	int const neighbors,
	int const threshold_1,
	int const threshold_3
)
{
	//tk get burial measure, three possible classes:
	//tk 1: exposed, 2: intermediate, 3:buried
	if ( neighbors > threshold_1 ) {
		if ( neighbors >= threshold_3 ) return 3;
		else return 2;
	}
	else return 1;
}

core::Real
hb_env_dep_burial_tk(int const nb1, int const nb2)
{
	//tk assign weight based on CB neighbors of two interacting residues

	// local
	int const exposed_threshold = { 10 };
	int const buried_threshold = { 20 };

	static FArray2D_double const burial3class_weight
		( 3, 3, burial3class_weight_initializer );

	return burial3class_weight(
		get_burial_3( nb1, exposed_threshold, buried_threshold ),
		get_burial_3( nb2, exposed_threshold, buried_threshold ) );
}

///////////////////////////////////////////////////////////////////////////////
Real
get_environment_dependent_weight(
	HBEvalType const & hbe_type,
	int const don_nb,
	int const acc_nb,
	HBondOptions const & options
)
{
	Real weight( 1.0 );
	// mjo why is this only applied to side chains!!!!
	if ( hbe_is_SC_type(hbe_type) ) {
		if ( options.smooth_hb_env_dep() ) {
			weight = hb_env_dep_burial_lin( acc_nb, don_nb );
		} else {
			weight = hb_env_dep_burial_tk( acc_nb, don_nb );
		}
	}
	// std::cout << "HB_ENV_WEIGHT: " << weight << std::endl;
	return weight;
}

///////////////////////////////////////////////////////////////////////////////
Real
get_membrane_depth_dependent_weight(
	pose::Pose const & pose,
	int const don_nb,
	int const acc_nb,
	Vector const & Hxyz, // proton
	Vector const & Axyz  // acceptor
)
{
	Real wat_weight(1.0), memb_weight(1.0), total_weight(1.0);

	// water phase smooth_hb_env_dep
	wat_weight = hb_env_dep_burial_lin( acc_nb, don_nb );

	// membrane phase dependent weight
	Vector const normal(MembraneEmbed_from_pose( pose ).normal());
	Vector const center(MembraneEmbed_from_pose( pose ).center());
	Real const thickness(Membrane_FAEmbed_from_pose( pose ).thickness());
	Real const steepness(Membrane_FAEmbed_from_pose( pose ).steepness());

	// Hdonor depth
	Real fa_depth_H = dot(Hxyz-center, normal);
	Real internal_product = std::abs(fa_depth_H);
	Real z = internal_product;
	z /= thickness;
	Real zn = std::pow( z, steepness );
	Real fa_proj_H = zn/(1 + zn);

	// Acc depth
	Real fa_depth_A = dot(Axyz-center, normal);
	internal_product = std::abs(fa_depth_A);
	z = internal_product;
	z /= thickness;
	zn = std::pow( z, steepness );
	Real fa_proj_A = zn/(1 + zn);

	Real fa_proj_AH = 0.5*(fa_proj_H+fa_proj_A);
	total_weight = fa_proj_AH * wat_weight + (1-fa_proj_AH) * memb_weight;

	//pbadebug
	//std::cout << "tot " << total_weight << " wat " << wat_weight << " memb " << memb_weight << " proj " 
	//<< fa_proj_AH << std::endl;

	return total_weight;
}
///////////////////////////////////////////////////////////////////////////////
Real
get_membrane_depth_dependent_weight(
	Vector const & normal,
	Vector const & center,
	Real const & thickness,
	Real const & steepness,
	int const don_nb,
	int const acc_nb,
	Vector const & Hxyz, // proton
	Vector const & Axyz  // acceptor
)
{
	Real wat_weight(1.0), memb_weight(1.0), total_weight(1.0);

	// water phase smooth_hb_env_dep
	wat_weight = hb_env_dep_burial_lin( acc_nb, don_nb );

	// membrane phase dependent weight
	// Hdonor depth
	Real fa_depth_H = dot(Hxyz-center, normal);
	Real internal_product = std::abs(fa_depth_H);
	Real z = internal_product;
	z /= thickness;
	Real zn = std::pow( z, steepness );
	Real fa_proj_H = zn/(1 + zn);

	// Acc depth
	Real fa_depth_A = dot(Axyz-center, normal);
	internal_product = std::abs(fa_depth_A);
	z = internal_product;
	z /= thickness;
	zn = std::pow( z, steepness );
	Real fa_proj_A = zn/(1 + zn);

	Real fa_proj_AH = 0.5*(fa_proj_H+fa_proj_A);
	total_weight = fa_proj_AH * wat_weight + (1-fa_proj_AH) * memb_weight;

	//pbadebug
	//std::cout << "tot " << total_weight << " wat " << wat_weight << " memb " << memb_weight << " proj " 
				//<< fa_proj_AH << std::endl;

	return total_weight;
}
///////////////////////////////////////////////////////////////////////////////

bool
nonzero_hbond_weight( ScoreFunction const & scorefxn )
{
	return ( scorefxn.has_nonzero_weight( hbond_lr_bb ) ||
		scorefxn.has_nonzero_weight( hbond_sr_bb ) ||
		scorefxn.has_nonzero_weight( hbond_bb_sc ) ||
		scorefxn.has_nonzero_weight( hbond_sr_bb_sc ) ||
		scorefxn.has_nonzero_weight( hbond_lr_bb_sc ) ||
		scorefxn.has_nonzero_weight( hbond_sc ) );
}

} // hbonds
} // scoring
} // core
