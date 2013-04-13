// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/score.bench.cc
///
/// @brief  Scoring Each benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_apps_benchmark_ScoreEach_bench_hh
#define INCLUDED_apps_benchmark_ScoreEach_bench_hh


#include <apps/benchmark/benchmark.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/Energies.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>


class ScoreEachBenchmark : public Benchmark
{

public:
	ScoreEachBenchmark(
		std::string name,
		core::scoring::ScoreType score_type,
		core::Size base_scale_factor
	) :
		Benchmark(name),
		pose_(),
		score_type_(score_type),
		base_scale_factor_(base_scale_factor),
		scorefxn_(),
		setup_successful_(false)
	{};


	virtual
	void
	set_scorefxn(
		core::scoring::ScoreFunctionOP scorefxn
	) {
		scorefxn_ = scorefxn;
	}

	virtual void setUp() {
		core::import_pose::pose_from_pdb(pose_, "test_in.pdb");
		if(!scorefxn_){
			scorefxn_ = new core::scoring::ScoreFunction();
		}
		try{
			// do this once in case there are one time setup requirements
			scorefxn_->set_weight(score_type_, 1);
			scorefxn_->score(pose_);
			pose_.energies().clear();
			scorefxn_->set_weight(score_type_, 0);
			setup_successful_ = true;
		} catch (utility::excn::EXCN_Base& excn){
			TR.Error
				<< "Unable to Test scoring with score type '"
				<< core::scoring::ScoreTypeManager::name_from_score_type(score_type_)
				<< "'" << std::endl
				<< "Fail with error:" << std::endl
				<< excn;
		}
	}

	virtual void run(core::Real scaleFactor) {
		if(!setup_successful_){
			TR.Error
				<< "Skipping running this test becuase the setup was not completed successfully." << std::endl;
			return;
		}

		core::scoring::ScoreFunction scorefxn;
		try{
			for(int i=0; i < base_scale_factor_*scaleFactor; i++) {
				scorefxn_->set_weight(score_type_, 1);
				scorefxn_->score(pose_);
				pose_.energies().clear();
				scorefxn_->set_weight(score_type_, 0);
			}
		} catch (utility::excn::EXCN_Base& excn){
			TR.Error
				<< "Unable to Test scoring with score type '"
				<< core::scoring::ScoreTypeManager::name_from_score_type(score_type_)
				<< "'" << std::endl
				<< "Fail with error:" << std::endl
				<< excn;
		}
	};

	virtual void tearDown() {};

private:
	core::pose::Pose pose_;

	core::scoring::ScoreType score_type_;
	core::Size base_scale_factor_;
	core::scoring::ScoreFunctionOP scorefxn_;
	bool setup_successful_;
};

using namespace core::scoring;

// fast
ScoreEachBenchmark Score_pro_close_("core.scoring.Score_10000x_pro_close",pro_close,10000);
ScoreEachBenchmark Score_rama2b_("core.scoring.Score_10000x_rama2b",rama2b,10000);
ScoreEachBenchmark Score_fa_stack_("core.scoring.Score_10000x_fa_stack",fa_stack,10000);
ScoreEachBenchmark Score_fa_stack_aro_("core.scoring.Score_10000x_fa_stack_aro",fa_stack_aro,10000);
ScoreEachBenchmark Score_fa_pair_("core.scoring.Score_10000x_fa_pair",fa_pair,10000);
ScoreEachBenchmark Score_fa_pair_aro_aro_("core.scoring.Score_10000x_fa_pair_aro_aro",fa_pair_aro_aro,10000);
ScoreEachBenchmark Score_fa_pair_aro_pol_("core.scoring.Score_10000x_fa_pair_aro_pol",fa_pair_aro_pol,1000);
ScoreEachBenchmark Score_fa_pair_pol_pol_("core.scoring.Score_10000x_fa_pair_pol_pol",fa_pair_pol_pol,10000);
ScoreEachBenchmark Score_dslf_ss_dst_("core.scoring.Score_10000x_dslf_ss_dst",dslf_ss_dst,10000);
ScoreEachBenchmark Score_dslf_cs_ang_("core.scoring.Score_10000x_dslf_cs_ang",dslf_cs_ang,10000);
ScoreEachBenchmark Score_dslf_ss_dih_("core.scoring.Score_10000x_dslf_ss_dih",dslf_ss_dih,10000);
ScoreEachBenchmark Score_dslf_ca_dih_("core.scoring.Score_10000x_dslf_ca_dih",dslf_ca_dih,10000);
ScoreEachBenchmark Score_dslf_cbs_ds_("core.scoring.Score_10000x_dslf_cbs_ds",dslf_cbs_ds,10000);
ScoreEachBenchmark Score_rama_("core.scoring.Score_10000x_rama",rama,10000);
ScoreEachBenchmark Score_omega_("core.scoring.Score_10000x_omega",omega,10000);
ScoreEachBenchmark Score_fa_dun_("core.scoring.Score_10000x_fa_dun",fa_dun,10000);
ScoreEachBenchmark Score_fa_dun_dev_("core.scoring.Score_10000x_fa_dun_dev",fa_dun_dev,10000);
ScoreEachBenchmark Score_fa_dun_rot_("core.scoring.Score_10000x_fa_dun_rot",fa_dun_rot,10000);
ScoreEachBenchmark Score_fa_dun_semi_("core.scoring.Score_10000x_fa_dun_semi",fa_dun_semi,10000);
ScoreEachBenchmark Score_dna_chi_("core.scoring.Score_10000x_dna_chi",dna_chi,10000);
ScoreEachBenchmark Score_p_aa_pp_("core.scoring.Score_10000x_p_aa_pp",p_aa_pp,10000);
ScoreEachBenchmark Score_yhh_planarity_("core.scoring.Score_10000x_yhh_planarity",yhh_planarity,10000);
ScoreEachBenchmark Score_h2o_intra_("core.scoring.Score_10000x_h2o_intra",h2o_intra,10000);
ScoreEachBenchmark Score_ref_("core.scoring.Score_10000x_ref",ref,10000);
ScoreEachBenchmark Score_ref_nc_("core.scoring.Score_10000x_ref_nc",ref_nc,10000);
ScoreEachBenchmark Score_nmer_ref_("core.scoring.Score_10000x_nmer_ref",nmer_ref,10000);
ScoreEachBenchmark Score_nmer_pssm_("core.scoring.Score_10000x_nmer_pssm",nmer_pssm,10000);
ScoreEachBenchmark Score_rna_bulge_("core.scoring.Score_10000x_rna_bulge",rna_bulge,10000);
ScoreEachBenchmark Score_special_rot_("core.scoring.Score_10000x_special_rot",special_rot,10000);
ScoreEachBenchmark Score_env_("core.scoring.Score_10000x_env",env,10000);
ScoreEachBenchmark Score_cbeta_("core.scoring.Score_10000x_cbeta",cbeta,10000);
ScoreEachBenchmark Score_rg_("core.scoring.Score_10000x_rg",rg,10000);
ScoreEachBenchmark Score_co_("core.scoring.Score_10000x_co",co,10000);
ScoreEachBenchmark Score_hs_pair_("core.scoring.Score_10000x_hs_pair",hs_pair,10000);
ScoreEachBenchmark Score_ss_pair_("core.scoring.Score_10000x_ss_pair",ss_pair,10000);
ScoreEachBenchmark Score_rsigma_("core.scoring.Score_10000x_rsigma",rsigma,10000);
ScoreEachBenchmark Score_sheet_("core.scoring.Score_10000x_sheet",sheet,10000);
ScoreEachBenchmark Score_chainbreak_("core.scoring.Score_10000x_chainbreak",chainbreak,10000);
ScoreEachBenchmark Score_linear_chainbreak_("core.scoring.Score_10000x_linear_chainbreak",linear_chainbreak,10000);
ScoreEachBenchmark Score_overlap_chainbreak_("core.scoring.Score_10000x_overlap_chainbreak",overlap_chainbreak,10000);
ScoreEachBenchmark Score_distance_chainbreak_("core.scoring.Score_10000x_distance_chainbreak",distance_chainbreak,10000);
ScoreEachBenchmark Score_dof_constraint_("core.scoring.Score_10000x_dof_constraint",dof_constraint,10000);
ScoreEachBenchmark Score_surface_("core.scoring.Score_10000x_surface",surface,10000);
ScoreEachBenchmark Score_p_aa_("core.scoring.Score_10000x_p_aa",p_aa,10000);
ScoreEachBenchmark Score_unfolded_("core.scoring.Score_10000x_unfolded",unfolded,10000);
ScoreEachBenchmark Score_e_pH_("core.scoring.Score_10000x_e_pH",e_pH,10000);

// average
ScoreEachBenchmark Score_fa_atr_("core.scoring.Score_1000x_fa_atr",fa_atr,1000);
ScoreEachBenchmark Score_fa_rep_("core.scoring.Score_1000x_fa_rep",fa_rep,1000);
ScoreEachBenchmark Score_fa_sol_("core.scoring.Score_1000x_fa_sol",fa_sol,1000);
ScoreEachBenchmark Score_fa_intra_atr_("core.scoring.Score_1000x_fa_intra_atr",fa_intra_atr,1000);
ScoreEachBenchmark Score_fa_intra_rep_("core.scoring.Score_1000x_fa_intra_rep",fa_intra_rep,1000);
ScoreEachBenchmark Score_fa_intra_sol_("core.scoring.Score_1000x_fa_intra_sol",fa_intra_sol,1000);
ScoreEachBenchmark Score_lk_hack_("core.scoring.Score_1000x_lk_hack",lk_hack,1000);
ScoreEachBenchmark Score_lk_ball_("core.scoring.Score_1000x_lk_ball",lk_ball,1000);
ScoreEachBenchmark Score_lk_ball_iso_("core.scoring.Score_1000x_lk_ball_iso)",lk_ball_iso,1000);
ScoreEachBenchmark Score_mm_lj_intra_rep_("core.scoring.Score_1000x_mm_lj_intra_rep",mm_lj_intra_rep,1000);
ScoreEachBenchmark Score_mm_lj_intra_atr_("core.scoring.Score_1000x_mm_lj_intra_atr",mm_lj_intra_atr,1000);
ScoreEachBenchmark Score_mm_lj_inter_rep_("core.scoring.Score_1000x_mm_lj_inter_rep",mm_lj_inter_rep,1000);
ScoreEachBenchmark Score_mm_lj_inter_atr_("core.scoring.Score_1000x_mm_lj_inter_atr",mm_lj_inter_atr,1000);
ScoreEachBenchmark Score_mm_bend_("core.scoring.Score_1000x_mm_bend",mm_bend,1000);
ScoreEachBenchmark Score_mm_stretch_("core.scoring.Score_1000x_mm_stretch",mm_stretch,1000);
ScoreEachBenchmark Score_lk_costheta_("core.scoring.Score_1000x_lk_costheta",lk_costheta,1000);
ScoreEachBenchmark Score_lk_polar_("core.scoring.Score_1000x_lk_polar",lk_polar,1000);
ScoreEachBenchmark Score_lk_nonpolar_("core.scoring.Score_1000x_lk_nonpolar",lk_nonpolar,1000);
ScoreEachBenchmark Score_h2o_hbond_("core.scoring.Score_1000x_h2o_hbond",h2o_hbond,1000);
ScoreEachBenchmark Score_ch_bond_("core.scoring.Score_1000x_ch_bond",ch_bond,1000);
ScoreEachBenchmark Score_gauss_("core.scoring.Score_1000x_gauss",gauss,1000);
ScoreEachBenchmark Score_hbond_sr_bb_("core.scoring.Score_1000x_hbond_sr_bb",hbond_sr_bb,1000);
ScoreEachBenchmark Score_hbond_lr_bb_("core.scoring.Score_1000x_hbond_lr_bb",hbond_lr_bb,1000);
ScoreEachBenchmark Score_hbond_bb_sc_("core.scoring.Score_1000x_hbond_bb_sc",hbond_bb_sc,1000);
ScoreEachBenchmark Score_hbond_sr_bb_sc_("core.scoring.Score_1000x_hbond_sr_bb_sc",hbond_sr_bb_sc,1000);
ScoreEachBenchmark Score_hbond_lr_bb_sc_("core.scoring.Score_1000x_hbond_lr_bb_sc",hbond_lr_bb_sc,1000);
ScoreEachBenchmark Score_hbond_sc_("core.scoring.Score_1000x_hbond_sc",hbond_sc,1000);
ScoreEachBenchmark Score_geom_sol_("core.scoring.Score_1000x_geom_sol",geom_sol,1000);
ScoreEachBenchmark Score_occ_sol_fitted_("core.scoring.Score_1000x_occ_sol_fitted",occ_sol_fitted,1000);
ScoreEachBenchmark Score_occ_sol_fitted_onebody_("core.scoring.Score_1000x_occ_sol_fitted_onebody",occ_sol_fitted_onebody,1000);
ScoreEachBenchmark Score_envsmooth_("core.scoring.Score_1000x_envsmooth",envsmooth,1000);
ScoreEachBenchmark Score_cart_bonded_("core.scoring.Score_1000x_cart_bonded",cart_bonded,1000);
ScoreEachBenchmark Score_neigh_vect_("core.scoring.Score_1000x_neigh_vect",neigh_vect,1000);
ScoreEachBenchmark Score_hack_elec_("core.scoring.Score_1000x_hack_elec",hack_elec,1000);

// slow 100x
ScoreEachBenchmark Score_mm_twist_("core.scoring.Score_100x_mm_twist",mm_twist,100);
ScoreEachBenchmark Score_sa_("core.scoring.Score_100x_sa",sa,100);
ScoreEachBenchmark Score_hpatch_("core.scoring.Score_100x_hpatch",hpatch,100);

// very slow 10x
ScoreEachBenchmark Score_gb_elec_("core.scoring.Score_10x_gb_elec",gb_elec,10);
ScoreEachBenchmark Score_pack_stat_("core.scoring.Score_10x_pack_stat",pack_stat,10);
ScoreEachBenchmark Score_fa_cust_pair_dist_("core.scoring.Score_10x_fa_cust_pair_dist",fa_cust_pair_dist,10);
ScoreEachBenchmark Score_occ_sol_exact_("core.scoring.Score_10x_occ_sol_exact",occ_sol_exact,10);

#endif // include guard
