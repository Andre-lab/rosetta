// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Andrew Leaver-Fay

// Unit Heaaders
#include <core/scoring/ScoreTypeManager.hh>


// Utility headers
#include <utility/exit.hh>

// C++ headers
#include <map>
#include <string>
#include <iostream>

#include <sstream>
#include <utility/vector1_bool.hh>

namespace core {
namespace scoring {

bool ScoreTypeManager::initialized_( false );
std::map< std::string, ScoreType > ScoreTypeManager::name2score_type_;
utility::vector1< std::string > ScoreTypeManager::score_type2name_;


void fill_score_range(std::map< std::string, ScoreType > & M, std::string prefix, int first, int last)
{
	M[ prefix + "_first" ] = ScoreType(first);
	M[ prefix + "_last" ] = ScoreType(last);
	for(int i=first+1; i<last; i++) {
		std::ostringstream s; s << prefix << '_' << i;
		M[ s.str() ] = ScoreType(i);
	}
}


/// @brief initialize the ScoreType name vector and map
///
/// @details initialize all the SCORETYPE string name into the vector then set up
/// the look-up map from string name to enum type
void
ScoreTypeManager::setup_score_type_names()
{
	if ( initialized_ ) return;
	initialized_ = true;

	name2score_type_[ "fa_atr" ] = fa_atr;
	name2score_type_[ "fa_rep" ] = fa_rep;
	name2score_type_[ "fa_sol" ] = fa_sol;
	name2score_type_[ "lk_hack" ] = lk_hack;
	name2score_type_[ "lk_costheta" ] = lk_costheta;
	name2score_type_[ "lk_polar" ] = lk_polar;
	name2score_type_[ "lk_nonpolar" ] = lk_nonpolar;
	name2score_type_[ "fa_intra_atr" ] = fa_intra_atr;
	name2score_type_[ "fa_intra_rep" ] = fa_intra_rep;
	name2score_type_[ "fa_intra_sol" ] = fa_intra_sol;
	name2score_type_[ "coarse_fa_atr" ] = coarse_fa_atr;
	name2score_type_[ "coarse_fa_rep" ] = coarse_fa_rep;
	name2score_type_[ "coarse_fa_sol" ] = coarse_fa_sol;
	name2score_type_[ "coarse_beadlj" ] = coarse_beadlj;
	name2score_type_[ "mm_lj_intra_rep" ] = mm_lj_intra_rep;
	name2score_type_[ "mm_lj_intra_atr" ] = mm_lj_intra_atr;
	name2score_type_[ "mm_lj_inter_rep" ] = mm_lj_inter_rep;
	name2score_type_[ "mm_lj_inter_atr" ] = mm_lj_inter_atr;
	name2score_type_[ "mm_twist" ] = mm_twist;
	name2score_type_[ "mm_bend" ] = mm_bend;
	name2score_type_[ "mm_stretch" ] = mm_stretch;
	name2score_type_[ "cart_bonded" ] = cart_bonded;
//	name2score_type_[ "csd_torsion" ] = csd_torsion; kwk commenting out csd atom type related code until I have implemented them fully
	name2score_type_[ "hack_elec" ] = hack_elec;
	name2score_type_[ "hack_elec_bb_bb" ] = hack_elec_bb_bb;
	name2score_type_[ "hack_elec_bb_sc" ] = hack_elec_bb_sc;
	name2score_type_[ "hack_elec_sc_sc" ] = hack_elec_sc_sc;
	name2score_type_[ "hack_elec_rna_phos_phos" ] = hack_elec_rna_phos_phos;
	name2score_type_[ "hack_elec_rna_phos_sugr" ] = hack_elec_rna_phos_sugr;
	name2score_type_[ "hack_elec_rna_phos_base" ] = hack_elec_rna_phos_base;
	name2score_type_[ "hack_elec_rna_sugr_sugr" ] = hack_elec_rna_sugr_sugr;
	name2score_type_[ "hack_elec_rna_sugr_base" ] = hack_elec_rna_sugr_base;
	name2score_type_[ "hack_elec_rna_base_base" ] = hack_elec_rna_base_base;
	name2score_type_[ "hack_elec_aro_aro" ] = hack_elec_aro_aro;
	name2score_type_[ "hack_elec_aro_all" ] = hack_elec_aro_all;
	name2score_type_[ "hack_aro" ] = hack_aro;
	name2score_type_[ "h2o_hbond" ] = h2o_hbond;
	name2score_type_[ "dna_dr" ] = dna_dr;
	name2score_type_[ "dna_bs" ] = dna_bs;
	name2score_type_[ "dna_bp" ] = dna_bp;
	name2score_type_[ "pro_close" ] = pro_close;
	name2score_type_[ "vdw" ] = vdw;
	name2score_type_[ "cenpack" ] = cenpack;
	name2score_type_[ "hybrid_vdw" ] = hybrid_vdw;
	name2score_type_[ "fa_cust_pair_dist" ] = fa_cust_pair_dist;

	// PyRosetta score types
	fill_score_range(name2score_type_, "PyRosettaTwoBodyContextIndepenedentEnergy", PyRosettaTwoBodyContextIndepenedentEnergy_first, PyRosettaTwoBodyContextIndepenedentEnergy_last);
	fill_score_range(name2score_type_, "PyRosettaTwoBodyContextDependentEnergy", PyRosettaTwoBodyContextDependentEnergy_first, PyRosettaTwoBodyContextDependentEnergy_last);
	fill_score_range(name2score_type_, "PyRosettaEnergy", PyRosettaEnergy_first, PyRosettaEnergy_last);

	name2score_type_[ "python" ] = python;

	name2score_type_[ "fastsaxs" ] = fastsaxs;
	name2score_type_[ "saxs_score" ] = saxs_score;
	name2score_type_[ "saxs_fa_score" ] = saxs_fa_score;
	name2score_type_[ "saxs_cen_score" ] = saxs_cen_score;
	name2score_type_[ "pddf_score" ] = pddf_score;

	name2score_type_[ "fa_pair" ] = fa_pair; // fa_pair == fa_pair_pol_pol
	name2score_type_[ "fa_pair_aro_aro" ] = fa_pair_aro_aro;
	name2score_type_[ "fa_pair_aro_pol" ] = fa_pair_aro_pol;
	name2score_type_[ "fa_pair_pol_pol" ] = fa_pair_pol_pol;
	name2score_type_[ "fa_plane" ] = fa_plane;
	name2score_type_[ "hbond_sr_bb" ] = hbond_sr_bb;
	name2score_type_[ "hbond_lr_bb" ] = hbond_lr_bb;
	name2score_type_[ "hbond_bb_sc" ] = hbond_bb_sc;
	name2score_type_[ "hbond_sr_bb_sc" ] = hbond_sr_bb_sc;
	name2score_type_[ "hbond_lr_bb_sc" ] = hbond_lr_bb_sc;
	name2score_type_[ "hbond_sc"    ] = hbond_sc;
	name2score_type_[ "interchain_pair" ] = interchain_pair;
	name2score_type_[ "interchain_vdw" ] = interchain_vdw;
	name2score_type_[ "interface_dd_pair" ] = interface_dd_pair;

	name2score_type_[ "ch_bond"    ] = ch_bond;
	name2score_type_[ "ch_bond_bb_bb" ] = ch_bond_bb_bb;
	name2score_type_[ "ch_bond_sc_sc" ] = ch_bond_sc_sc;
	name2score_type_[ "ch_bond_bb_sc" ] = ch_bond_bb_sc;

	name2score_type_[ "neigh_vect"  ] = neigh_vect;
	name2score_type_[ "neigh_count" ] = neigh_count;
	name2score_type_[ "neigh_vect_raw"] = neigh_vect_raw;
	name2score_type_[ "symE_bonus"  ] = symE_bonus;
	name2score_type_[ "sym_lig"  ] = sym_lig;

	name2score_type_[ "orbitals_hpol" ] = orbitals_hpol;
	name2score_type_[ "orbitals_haro" ] = orbitals_haro;
	name2score_type_["orbitals_hpol_bb"] = orbitals_hpol_bb;
	name2score_type_["orbitals_orbitals"] = orbitals_orbitals;

	name2score_type_[ "geom_sol"    ] = geom_sol;
	name2score_type_[ "occ_sol_fitted"    ] = occ_sol_fitted;
	name2score_type_[ "occ_sol_fitted_onebody"    ] = occ_sol_fitted_onebody;
	name2score_type_[ "occ_sol_exact"    ] = occ_sol_exact;

	name2score_type_[ "gb_elec" ] = gb_elec;
	name2score_type_[ "PB_elec" ] = PB_elec;
	name2score_type_[ "dslf_ss_dst" ] = dslf_ss_dst;
	name2score_type_[ "dslf_cs_ang" ] = dslf_cs_ang;
	name2score_type_[ "dslf_ss_dih" ] = dslf_ss_dih;
	name2score_type_[ "dslf_ca_dih" ] = dslf_ca_dih;
	name2score_type_[ "dslf_cbs_ds" ] = dslf_cbs_ds;
	name2score_type_[ "dslfc_cen_dst" ] = dslfc_cen_dst;
	name2score_type_[ "dslfc_cb_dst"  ] = dslfc_cb_dst;
	name2score_type_[ "dslfc_ang"     ] = dslfc_ang;
	name2score_type_[ "dslfc_cb_dih"  ] = dslfc_cb_dih;
	name2score_type_[ "dslfc_bb_dih"  ] = dslfc_bb_dih;

	name2score_type_[ "dslfc_rot"  ] = dslfc_rot;
	name2score_type_[ "dslfc_trans"  ] = dslfc_trans;
	name2score_type_[ "dslfc_RT"  ] = dslfc_RT;

	name2score_type_[ "custom_atom_pair" ] = custom_atom_pair;
	name2score_type_[ "atom_pair_constraint" ] = atom_pair_constraint;
	name2score_type_[ "dunbrack_constraint" ] = dunbrack_constraint;
	name2score_type_[ "angle_constraint" ] = angle_constraint;
	name2score_type_[ "dihedral_constraint" ] = dihedral_constraint;
	name2score_type_[ "big_bin_constraint" ] = big_bin_constraint;
	name2score_type_[ "constant_constraint" ] = constant_constraint;
	name2score_type_[ "coordinate_constraint" ] = coordinate_constraint;
	name2score_type_[ "site_constraint" ] = site_constraint;

	name2score_type_[ "rama"    ] = rama;
	name2score_type_[ "rama2b"  ] = rama2b;
	name2score_type_[ "omega"    ] = omega;
	name2score_type_[ "fa_dun" ] = fa_dun;
	name2score_type_[ "p_aa_pp" ] = p_aa_pp;
	name2score_type_[ "h2o_intra" ] =  h2o_intra;
	name2score_type_[ "ref" ] = ref;
	name2score_type_[ "seqdep_ref" ] = seqdep_ref;
	name2score_type_[ "envsmooth" ] = envsmooth;
	name2score_type_[ "e_pH" ] = e_pH;
	name2score_type_[ "rna_bulge"] = rna_bulge;
	// Variant type to flag rotamers for alternative scoring with varying weight
	name2score_type_[ "special_rot"] = special_rot;

	name2score_type_[ "env" ]    = env;
	name2score_type_[ "burial" ] = burial;
	name2score_type_[ "abego" ]  = abego;
	name2score_type_[ "pair" ]   = pair;
	name2score_type_[ "cbeta" ]  = cbeta;
	name2score_type_[ "DFIRE" ]  = DFIRE;

	//bw membrane scoring terms
	name2score_type_[ "Menv" ] = Menv;
	name2score_type_[ "Menv_non_helix" ] = Menv_non_helix;
	name2score_type_[ "Menv_termini" ] = Menv_termini;
	name2score_type_[ "Menv_tm_proj" ] = Menv_tm_proj;
	name2score_type_[ "Mcbeta" ] = Mcbeta;
	name2score_type_[ "Mpair" ] = Mpair;
	name2score_type_[ "Mlipo" ] = Mlipo;
	//pba membrane all atom terms
	name2score_type_[ "fa_mbenv" ] = fa_mbenv;
	name2score_type_[ "fa_mbsolv" ] = fa_mbsolv;
	name2score_type_[ "Menv_smooth" ] = Menv_smooth;

	name2score_type_[ "rg" ] = rg;
	name2score_type_[ "co" ] = co;
	name2score_type_[ "peptide_bond" ] = peptide_bond;
	name2score_type_[ "pcs" ] = pcs;
	name2score_type_[ "pcs2" ] = pcs2;
	name2score_type_[ "dock_ens_conf" ] = dock_ens_conf;

	name2score_type_[ "hs_pair" ] = hs_pair;
	name2score_type_[ "ss_pair" ] = ss_pair;
	name2score_type_[ "rsigma" ] = rsigma;
	name2score_type_[ "sheet" ] = sheet;
	name2score_type_[ "rdc" ] = rdc;
	name2score_type_[ "rdc_segments" ] = rdc_segments;
	name2score_type_[ "rdc_rohl" ] =rdc_rohl;
	name2score_type_[ "holes" ] = holes;
	name2score_type_[ "holes_resl" ] = holes_resl;
	name2score_type_[ "holes_decoy" ] = holes_decoy;
	name2score_type_[ "holes_min" ] = holes_min;
	name2score_type_[ "holes_min_mean" ] = holes_min_mean;
	name2score_type_[ "dab_sasa" ] = dab_sasa;
	name2score_type_[ "dab_sev" ] = dab_sev;
	name2score_type_[ "sa" ] = sa;

	name2score_type_[ "interchain_env"] = interchain_env;
	name2score_type_[ "interchain_contact"] = interchain_contact;

	name2score_type_[ "rna_rg"] = rna_rg;
	name2score_type_[ "rna_vdw"] = rna_vdw;
	name2score_type_[ "rna_base_backbone"] = rna_base_backbone;
	name2score_type_[ "rna_backbone_backbone"] = rna_backbone_backbone;
	name2score_type_[ "rna_repulsive"] = rna_repulsive;

	name2score_type_[ "rna_base_pair"] = rna_base_pair;
	name2score_type_[ "rna_base_axis"] = rna_base_axis;
	name2score_type_[ "rna_base_stagger"] = rna_base_stagger;
	name2score_type_[ "rna_base_stack"] = rna_base_stack;
	name2score_type_[ "rna_base_stack_axis"] = rna_base_stack_axis;
	name2score_type_[ "rna_data_base"] = rna_data_base;
	name2score_type_[ "rna_data_backbone"] = rna_data_backbone;

	//Will these ever really be used?
	name2score_type_[ "rna_base_pair_pairwise"] = rna_base_pair_pairwise;
	name2score_type_[ "rna_base_axis_pairwise"] = rna_base_axis_pairwise;
	name2score_type_[ "rna_base_stagger_pairwise"] = rna_base_stagger_pairwise;
	name2score_type_[ "rna_base_stack_pairwise"] = rna_base_stack_pairwise;
	name2score_type_[ "rna_base_stack_axis_pairwise"] = rna_base_stack_axis_pairwise;

	name2score_type_[ "fa_stack"] = fa_stack;
	//	name2score_type_[ "fa_stack_purine"] = fa_stack_purine;
	//	name2score_type_[ "fa_stack_pyrimidine"] = fa_stack_pyrimidine;
	name2score_type_[ "rna_torsion"] = rna_torsion;
	name2score_type_[ "rna_sugar_close"] = rna_sugar_close;
	name2score_type_[ "rna_bond_geometry"] = rna_bond_geometry;
	name2score_type_[ "rna_fa_atr_base"] = rna_fa_atr_base;
	name2score_type_[ "rna_fa_rep_base"] = rna_fa_rep_base;

	name2score_type_[ "chainbreak" ] = chainbreak;
	name2score_type_[ "linear_chainbreak" ] = linear_chainbreak;
	name2score_type_[ "overlap_chainbreak" ] = overlap_chainbreak;
	name2score_type_[ "distance_chainbreak" ] = distance_chainbreak;
	name2score_type_[ "dof_constraint" ] = dof_constraint;
	name2score_type_[ "rms_energy" ] = rms;
	name2score_type_[ "suck" ] = suck;
	name2score_type_[ "res_type_constraint" ] = res_type_constraint;
	name2score_type_[ "pocket_constraint" ] = pocket_constraint;
	name2score_type_[ "backbone_stub_constraint" ] = backbone_stub_constraint;
	name2score_type_[ "suck" ] = suck;
	name2score_type_[ "pack_stat" ] = pack_stat;

	name2score_type_[ "surface" ] = surface;
	name2score_type_[ "hpatch" ] = hpatch;
	name2score_type_[ "p_aa" ] = p_aa;
	name2score_type_[ "unfolded" ] = unfolded;

	name2score_type_[ "elec_dens_fast" ] = elec_dens_fast;
	name2score_type_[ "elec_dens_window" ] = elec_dens_window;
	name2score_type_[ "elec_dens_whole_structure_ca" ] = elec_dens_whole_structure_ca;
	name2score_type_[ "elec_dens_whole_structure_allatom" ] = elec_dens_whole_structure_allatom;
	name2score_type_[ "patterson_cc" ] = patterson_cc;

	name2score_type_[ "natbias_ss" ] = natbias_ss;
	name2score_type_[ "natbias_hs" ] = natbias_hs;
	name2score_type_[ "natbias_hh" ] = natbias_hh;
	name2score_type_[ "natbias_stwist" ] = natbias_stwist;

	name2score_type_[ "aa_cmp" ] = aa_cmp;

	name2score_type_[ "total_score" ] = total_score;


	assert( name2score_type_.size() == end_of_score_type_enumeration );

	score_type2name_.resize( end_of_score_type_enumeration );
	for ( std::map< std::string, ScoreType >::const_iterator iter = name2score_type_.begin(),
					iter_end = name2score_type_.end(); iter != iter_end; ++iter ) {
		score_type2name_[ iter->second ] = iter->first;
	}

}



//////////////////////////////////////////////////////////////////////////////
/// @brief give a ScoreType string name and return its enum type
ScoreType
ScoreTypeManager::score_type_from_name( std::string const & name )
{
	setup_score_type_names();
	std::map< std::string, ScoreType >::const_iterator iter( name2score_type_.find( name ) );
	if ( iter == name2score_type_.end() ) {
		utility_exit_with_message("unrecognized score_type type "+name);
	}
	return iter->second;
}

std::string
ScoreTypeManager::name_from_score_type( ScoreType score_type )
{
	setup_score_type_names();
	return score_type2name_[ score_type ];
}

///@brief
bool
ScoreTypeManager::is_score_type( std::string const & name )
{
	setup_score_type_names();
	std::map< std::string, ScoreType >::const_iterator iter( name2score_type_.find( name ) );
	return iter != name2score_type_.end();
}

}
}
