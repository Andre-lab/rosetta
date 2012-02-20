// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_LoopCloseSampler
/// @brief Loop Close Sampler...
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/RNA_LoopCloseSampler.hh>
#include <protocols/swa/rna/RNA_AnalyticLoopCloser.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Conformation.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/exit.hh>
#include <time.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>


using namespace core;
using core::Real;


namespace protocols {
namespace swa {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
RNA_LoopCloseSampler::RNA_LoopCloseSampler ( Size const moving_suite, Size const chainbreak_suite ) :
	moving_suite_ ( moving_suite ),
	chainbreak_suite_ ( chainbreak_suite ),
	scorefxn_ ( core::scoring::getScoreFunction() ),
	bin_size_ ( 20 ),
	n_construct_ ( 0 ),
	epsilon_range_ ( 40.0 ),
	rep_cutoff_ ( 0.1 ),
	torsion_range_ ( 20.0 ),
	torsion_increment_ ( 5.0 ),
	just_output_score_ ( false ),
	sample_only_ ( false ),
	sample_native_torsion_ ( false ) {
	RNA_AnalyticLoopCloser rna_analytic_loop_closer ( moving_suite, chainbreak_suite );
	initialize_rep_scorefxn();
}

//////////////////////////////////////////////////////////////////////////
//destructor
RNA_LoopCloseSampler::~RNA_LoopCloseSampler()
{}

//////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::apply ( core::pose::Pose & pose ) {
	using namespace pose;
	using namespace scoring;
	using namespace io::silent;
	using namespace chemical;
	using namespace id;
	using namespace scoring::rna;
	using namespace protocols::swa;
	using namespace protocols::swa::rna;
	using namespace numeric::conversions;
	Real const fa_rep_score_baseline = initialize_fa_rep ( pose, utility::tools::make_vector1 ( moving_suite_ ), rep_scorefxn_ );
	// iterate over 4 dofs:
	//     epsilon1, zeta1, alpha1;    alpha2
	// 	solve for the other 6 by chain closure:
	//     beta1, gamma1;               epsilon2, zeta2, beta2, gamma2
	// This should probably be encapsulated into a PoseSampleGenerator or something similar
	int const bins1_ = 360 / bin_size_ ; //This is total bins, default is 18
	int const bins2_ = bins1_ / 2; //This is total bins divided by 2; default is 9
	int const bins3_ = bins1_ / 3; //This is total bins divided by 3; default is 6
	int const bins4_ = 1 + 40 / bin_size_; //This is the bin for chi and episilon, these two torsion angles vary from -20+mean to 20+mean
	PuckerState pucker_state1 = Get_residue_pucker_state ( pose, moving_suite_ );
	Real epsilon1, zeta1, alpha1 ( 0.0 ), alpha2 ( 0.0 ), perturb_epsilon1, perturb_zeta1, perturb_alpha1, perturb_alpha2;
	Real beta1, beta2, epsilon2;
	// following is only used if we are estimating jacobians (numerically).
	utility::vector1< utility::vector1< utility::vector1< Real > > > perturbed_solution_torsions;
	utility::vector1< utility::vector1< Real > > J; // 6 x 6 Jacobian.
	utility::vector1< Real > six_zeros;

	for ( int i = 1; i <= 6; i++ ) six_zeros.push_back ( 0.0 );

	for ( int i = 1; i <= 6; i++ ) J.push_back ( six_zeros );

	Real const perturbation_size ( 1.0e-5 );
	Real const detJ_cutoff_ ( 1.0 );
	Real epsilon1_center = ( pucker_state1 == NORTH ) ? -150.17 : -98.45;
	Real epsilon1_min = epsilon1_center - epsilon_range_;
	Real epsilon1_max = epsilon1_center + epsilon_range_;
	if (pucker_state1 == SOUTH && epsilon1_min > -178) epsilon1_min = -178.45;
	Real epsilon1_increment = bin_size_;
	Real alpha1_min = 20.0;
	Real alpha1_max = 340.0 - bin_size_;
	Real alpha1_increment = bin_size_;
	Real alpha2_min = 20.0;
	Real alpha2_max = 340.0 - bin_size_;
	Real alpha2_increment = bin_size_;
	Real zeta1_min = 20.0;
	Real zeta1_max = 340.0 - bin_size_;
	Real zeta1_increment = bin_size_;
	// move this to its own function?
	PoseCOP native_pose = get_native_pose();
	Real const epsilon1_native = pose.torsion ( TorsionID ( moving_suite_  , id::BB, EPSILON ) );
	Real const zeta1_native    = pose.torsion ( TorsionID ( moving_suite_  , id::BB, ZETA ) );
	Real const alpha1_native   = pose.torsion ( TorsionID ( moving_suite_ + 1, id::BB, ALPHA ) );
	Real const alpha2_native   = pose.torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, ALPHA ) );

	if ( sample_native_torsion_ ) {
		epsilon1_max = epsilon1_max + bin_size_;
		alpha1_max = alpha1_max + bin_size_;
		zeta1_max = zeta1_max + bin_size_;
		alpha2_max = alpha2_max + bin_size_;
	}

	RNA_AnalyticLoopCloser rna_analytic_loop_closer ( moving_suite_, chainbreak_suite_ );

	for ( Real epsilon1 = epsilon1_min; epsilon1 <= epsilon1_max; epsilon1 += epsilon1_increment ) {
		for ( Real alpha1 = alpha1_min; alpha1 <= alpha1_max; alpha1 += alpha1_increment ) {
			for ( Real alpha2 = alpha2_min; alpha2 <= alpha2_max; alpha2 += alpha2_increment ) {
				for ( Real zeta1 = zeta1_min; zeta1 <= zeta1_max; zeta1 += zeta1_increment ) {
					pose.set_torsion ( TorsionID ( moving_suite_,       id::BB, EPSILON ), epsilon1 );
					pose.set_torsion ( TorsionID ( moving_suite_,       id::BB, ZETA ),    zeta1 );
					pose.set_torsion ( TorsionID ( moving_suite_ + 1,     id::BB, ALPHA ),   alpha1 );
					pose.set_torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, ALPHA ),   alpha2 );

					if ( sample_native_torsion_ ) {
						if ( epsilon1 == epsilon1_max ) pose.set_torsion ( TorsionID ( moving_suite_, id::BB, EPSILON ), epsilon1_native );

						if ( zeta1 == zeta1_max ) pose.set_torsion ( TorsionID ( moving_suite_, id::BB, ZETA ), zeta1_native );

						if ( alpha1 == alpha1_max ) pose.set_torsion ( TorsionID ( moving_suite_ + 1, id::BB, ALPHA ), alpha1_native );

						if ( alpha2 == alpha2_max ) pose.set_torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, ALPHA ), alpha2_native );
					}

					//close loop.
					rna_analytic_loop_closer.apply ( pose );
					bool perturb_solutions_calculated = false;

					// iterate over solutions -- anything with OK repulsive?
					for ( Size n = 1; n <= rna_analytic_loop_closer.nsol(); n++ ) {
						rna_analytic_loop_closer.fill_solution ( pose, n );
						utility::vector1< Real > const & solution_torsions = rna_analytic_loop_closer.get_torsions ( n );

						if ( !torsion_angles_within_cutoffs ( pose, moving_suite_, chainbreak_suite_, bin_size_, bins2_ ) ) continue;

						if ( !sample_only_ ) {
							if ( !check_clash ( pose, fa_rep_score_baseline, rep_cutoff_, rep_scorefxn_ ) ) continue;
						}

						// save data
						n_construct_++;
						torsion_info_.clear();
						torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_  , id::BB, EPSILON ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_  , id::BB, ZETA ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_ + 1, id::BB, ALPHA ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_ + 1, id::BB, BETA ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_ + 1, id::BB, GAMMA ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_  , id::BB, EPSILON ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_  , id::BB, ZETA ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, ALPHA ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, BETA ) ) );
						torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, GAMMA ) ) );
						all_torsion_info_.push_back ( torsion_info_ );
					}
				}
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::fill_pose ( pose::Pose & pose, Size construct_number ) {
	using namespace pose;
	using namespace scoring;
	using namespace io::silent;
	using namespace chemical;
	using namespace id;
	using namespace scoring::rna;
	using namespace protocols::swa;
	using namespace protocols::swa::rna;
	using namespace numeric::conversions;
	utility::vector1< Real > & torsion_info = all_torsion_info_[construct_number];
	pose.set_torsion ( TorsionID ( moving_suite_  , id::BB, EPSILON ), torsion_info[1] );
	pose.set_torsion ( TorsionID ( moving_suite_  , id::BB, ZETA ), torsion_info[2] );
	pose.set_torsion ( TorsionID ( moving_suite_ + 1, id::BB, ALPHA ), torsion_info[3] );
	pose.set_torsion ( TorsionID ( moving_suite_ + 1, id::BB, BETA ), torsion_info[4] );
	pose.set_torsion ( TorsionID ( moving_suite_ + 1, id::BB, GAMMA ), torsion_info[5] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_  , id::BB, EPSILON ), torsion_info[6] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_  , id::BB, ZETA ), torsion_info[7] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, ALPHA ), torsion_info[8] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, BETA ), torsion_info[9] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, GAMMA ), torsion_info[10] );
	return;
}
///////////////////////////////
void
RNA_LoopCloseSampler::clear_all() {
	all_torsion_info_.clear();
	n_construct_ = 0;
	return;
}

//////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::initialize_rep_scorefxn() {
	using namespace core::scoring;
	rep_scorefxn_ = new ScoreFunction;
	rep_scorefxn_->set_weight ( fa_rep, 0.12 );
}

///////////////////////////////////////////////////////////////
Real
RNA_LoopCloseSampler::initialize_fa_rep ( pose::Pose const & pose,
    utility::vector1< Size > const & moving_suites,
    scoring::ScoreFunctionOP rep_scorefxn ) {
	using namespace pose;
	using namespace kinematics;
	using namespace scoring;
	using namespace protocols::swa;

	if ( sample_only_ ) {
		return 0;
	}

	Pose pose_expand = pose;

	for ( Size n = 1; n <= moving_suites.size(); n++ ) {
		Size const jump_at_moving_suite = make_cut_at_moving_suite ( pose_expand, moving_suites[n] );
		Jump j = pose_expand.jump ( jump_at_moving_suite );
		j.set_translation ( Vector ( 1.0e4 * n, 0.0, 0.0 ) );
		pose_expand.set_jump ( jump_at_moving_suite, j );
	}

	( *rep_scorefxn ) ( pose_expand );
	EnergyMap const & energy_map = pose_expand.energies().total_energies();
	return energy_map[ fa_rep ] * rep_scorefxn->get_weight ( fa_rep );
}


///////////////////////////////////////////////////////////////
bool
RNA_LoopCloseSampler::check_clash ( pose::Pose & pose,
                                    Real const & fa_rep_score_baseline,
                                    Real const & rep_cutoff_,
                                    scoring::ScoreFunctionOP rep_scorefxn ) {
	using namespace scoring;
	( *rep_scorefxn ) ( pose );
	EnergyMap const & energy_map = pose.energies().total_energies();
	Real const fa_rep_score = energy_map[ fa_rep ] * rep_scorefxn->get_weight ( fa_rep );
	//	std::cout << fa_rep_score << " " << fa_rep_score_baseline << std::endl;

	if ( ( fa_rep_score - fa_rep_score_baseline ) > rep_cutoff_ ) return false;

	static Real const tolerance ( 1.0e-3 );

	if ( ( fa_rep_score - fa_rep_score_baseline ) < -1.0 * tolerance ) {
		std::cout << fa_rep_score << " " << fa_rep_score_baseline << std::endl;
		//		utility_exit_with_message( "Weird fa_rep?" );
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::set_scorefxn ( core::scoring::ScoreFunctionOP const & scorefxn ) {
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
bool
RNA_LoopCloseSampler::torsion_angles_within_cutoffs ( pose::Pose const & pose,
    Size const moving_suite,
    Size const chainbreak_suite,
    int const bin_size_,
    int const bins2_ ) {
	using namespace id;
	using namespace core::scoring::rna;
	using namespace protocols::swa::rna;
	// Quick check on range of epsilon, beta. Other torsion angles span full 360.0 range.
	Real const beta1 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( moving_suite + 1, id::BB, BETA ) ) );

	if ( beta1 < 60 && beta1 > -60 ) {
		return false;
	}

	Real const beta2 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( chainbreak_suite + 1, id::BB, BETA ) ) );

	if ( beta2 < 60 && beta2 > -60 ) {
		return false;
	}

	PuckerState pucker_state1 = Get_residue_pucker_state ( pose, moving_suite );
	Real const epsilon1 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( moving_suite, id::BB, EPSILON ) ) );
	Real const epsilon1_ideal = ( pucker_state1 == NORTH ) ? -150.17 : -98.45;

	if ( ( epsilon1 < ( epsilon1_ideal - epsilon_range_ ) ) ||
	     ( epsilon1 > ( epsilon1_ideal + epsilon_range_ ) ) ) {
		return false;
	}

	PuckerState pucker_state2 = Get_residue_pucker_state ( pose, chainbreak_suite );
	Real const epsilon2 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( chainbreak_suite, id::BB, EPSILON ) ) );
	Real const epsilon2_ideal = ( pucker_state2 == NORTH ) ? -150.17 : -98.45;

	if ( ( epsilon2 < ( epsilon2_ideal - epsilon_range_ ) ) ||
	     ( epsilon2 > ( epsilon2_ideal + epsilon_range_ ) ) ) {
		return false;
	}

	return true;
}

}
}
}
