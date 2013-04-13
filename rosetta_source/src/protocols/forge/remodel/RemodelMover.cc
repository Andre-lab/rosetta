// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelMover.cc
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/remodel/RemodelGlobalFrame.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.fwd.hh>
#include <protocols/forge/remodel/RemodelMoverCreator.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
// AUTO-REMOVED #include <core/pose/symmetry/util.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>

// project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Remarks.hh>
#include <core/pose/util.hh> // for pdbinfo
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/remodel/RemodelMoverCreator.hh>
#include <protocols/forge/methods/chainbreak_eval.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/forge/methods/util.hh>
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.fwd.hh>
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/forge/remodel/RemodelAccumulator.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_KIC.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/PackRotamersMover.fwd.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh> // dihedral constraint
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToNeighborhoodOperation.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
// AUTO-REMOVED #include <protocols/toolbox/task_operations/PreventChainFromRepackingOperation.hh>
// AUTO-REMOVED #include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/forge/remodel/RemodelAccumulator.hh>
#include <protocols/forge/remodel/RemodelEnzdesCstModule.hh>


// Parser headers
#include <protocols/moves/DataMap.hh>
#include <utility/tag/Tag.hh>

#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

// C++ headers
#include <utility>


//#define FILE_DEBUG 1

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace ObjexxFCL::fmt;
using namespace core;
using namespace core::scoring;

namespace protocols {
namespace forge {
namespace remodel {


static basic::Tracer TR( "protocols.forge.remodel.RemodelMover" );

// parser
std::string
RemodelMoverCreator::keyname() const {
	return RemodelMoverCreator::mover_name();
}

protocols::moves::MoverOP
RemodelMoverCreator::create_mover() const {
	return new RemodelMover;
}

std::string
RemodelMoverCreator::mover_name() {
	return "RemodelMover";
}


///
/// @begin RemodelMover::RemodelMover
///
/// @brief
/// Default constructor. Checks values of options packing::soft_rep_design, and remodel::dr_cycles
/// Creates and sets the centroid and fullatom score functions.
///
RemodelMover::RemodelMover() :
	Super( "RemodelMover" ),
	// User specified options
	op_user_remodel_vdw_(option[OptionKeys::remodel::vdw].user()),
	op_user_remodel_rama_(option[OptionKeys::remodel::rama].user()),
	op_user_remodel_cbeta_(option[OptionKeys::remodel::cbeta].user()),
	op_user_remodel_cenpack_(option[OptionKeys::remodel::cenpack].user()),
	op_user_remodel_hb_lrbb_(option[OptionKeys::remodel::hb_lrbb].user()),
	op_user_remodel_hb_srbb_(option[OptionKeys::remodel::hb_srbb].user()),
	op_user_remodel_rsigma_(option[OptionKeys::remodel::rsigma].user()),
	op_user_remodel_ss_pair_(option[OptionKeys::remodel::ss_pair].user()),
	op_user_domainFusion_insert_segment_from_pdb_(option[OptionKeys::remodel::domainFusion::insert_segment_from_pdb].user()),
	op_user_remodel_repeat_structure_(option[OptionKeys::remodel::repeat_structure].user()),
	op_user_symmetry_symmetry_definition_(option[OptionKeys::symmetry::symmetry_definition].user()),
	op_user_remodel_quick_and_dirty_(option[OptionKeys::remodel::quick_and_dirty].user()),
	op_user_enzdes_cstfile_(option[OptionKeys::enzdes::cstfile].user()),
	op_user_constraints_cst_file_(option[OptionKeys::constraints::cst_file].user()),
	op_user_remodel_build_disulf_(option[OptionKeys::remodel::build_disulf].user()),
	op_user_remodel_rg_(option[OptionKeys::remodel::rg].user()),
	op_user_remodel_check_scored_centroid_(option[OptionKeys::remodel::check_scored_centroid].user()),
	op_user_remodel_run_confirmation_(option[OptionKeys::remodel::run_confirmation].user()),
	op_user_remodel_free_relax_(option[OptionKeys::remodel::free_relax].user()),
	op_user_RemodelLoopMover_cyclic_peptide_(option[OptionKeys::remodel::RemodelLoopMover::cyclic_peptide].user()),
	op_user_remodel_swap_refine_confirm_protocols_(option[OptionKeys::remodel::swap_refine_confirm_protocols].user()),
	op_user_design_allow_rare_aro_chi_(option[OptionKeys::remodel::design::allow_rare_aro_chi].user()),
	op_user_packing_soft_rep_design_(option[OptionKeys::packing::soft_rep_design].user()),
	op_user_remodel_blueprint_(option[OptionKeys::remodel::blueprint].user()),
	op_user_remodel_checkpoint_(option[OptionKeys::remodel::checkpoint].user()),
	op_user_remodel_use_pose_relax_(option[OptionKeys::remodel::use_pose_relax].user()),
	op_user_remodel_use_cart_relax_(option[OptionKeys::remodel::use_cart_relax].user()),
	op_user_remodel_use_blueprint_sequence_(option[OptionKeys::remodel::use_blueprint_sequence].user()),
	op_user_remodel_bypass_fragments_(option[OptionKeys::remodel::bypass_fragments].user()),
	op_RemodelLoopMover_max_linear_chainbreak_(option[OptionKeys::remodel::RemodelLoopMover::max_linear_chainbreak]),
	op_remodel_dr_cycles_(option[OptionKeys::remodel::dr_cycles]),
	op_remodel_cen_sfxn_(option[OptionKeys::remodel::cen_sfxn]),
	op_remodel_vdw_(option[OptionKeys::remodel::vdw]),
	op_remodel_rama_(option[OptionKeys::remodel::rama]),
	op_remodel_cbeta_(option[OptionKeys::remodel::cbeta]),
	op_remodel_cenpack_(option[OptionKeys::remodel::cenpack]),
	op_remodel_hb_lrbb_(option[OptionKeys::remodel::hb_lrbb]),
	op_remodel_hb_srbb_(option[OptionKeys::remodel::hb_srbb]),
	op_remodel_rsigma_(option[OptionKeys::remodel::rsigma]),
	op_remodel_ss_pair_(option[OptionKeys::remodel::ss_pair]),
	op_run_show_simulation_in_pymol_(option[OptionKeys::run::show_simulation_in_pymol]),
	op_remodel_num_trajectory_(option[OptionKeys::remodel::num_trajectory]),
	op_remodel_repeat_structure_(option[OptionKeys::remodel::repeat_structure]),
	op_remodel_cen_minimize_(option[OptionKeys::remodel::cen_minimize]),
	op_remodel_cstfilter_(option[OptionKeys::remodel::cstfilter]),
	op_remodel_collect_clustered_top_(option[OptionKeys::remodel::collect_clustered_top]),
	op_out_prefix_(option[OptionKeys::out::prefix]),
	op_remodel_use_same_length_fragments_(option[OptionKeys::remodel::use_same_length_fragments]),
	op_packing_soft_rep_design_(false),
	op_remodel_rg_(0.0),
	op_remodel_blueprint_(""),
	op_remodel_checkpoint_(false),
	op_remodel_use_pose_relax_(false),
	op_remodel_use_cart_relax_(false),
	op_remodel_use_blueprint_sequence_(false),
	op_remodel_bypass_fragments_(false),
	op_remodel_free_relax_(false)
{
	// Do this first
	set_param_from_options();

	//use_fullmer_ = false;
	//use_sequence_bias_ = false;
	//max_linear_chainbreak_ = 0.15;
	max_linear_chainbreak_ = op_RemodelLoopMover_max_linear_chainbreak_;
	//centroid_loop_mover_str_ = "quick_ccd";
	centroid_loop_mover_str_ = "RemodelLoopMover";
	redesign_loop_neighborhood_ = false;
	dr_cycles_ = op_remodel_dr_cycles_;
	centroid_sfx_ = core::scoring::ScoreFunctionFactory::create_score_function( op_remodel_cen_sfxn_ );
	fullatom_sfx_ = core::scoring::getScoreFunction();

	register_user_options();

	if ( op_packing_soft_rep_design_ ) {
		TR << "SWITCHING FULLATOM FUNCITON TO SOFT_REP_DESIGN" << std::endl;
		fullatom_sfx_ = ScoreFunctionFactory::create_score_function( SOFT_REP_DESIGN_WTS );
	}
}

///
/// @begin RemodelMover::RemodelMover( RemodelMover )
///
/// @brief
/// Copy constructor.
///
RemodelMover::RemodelMover( RemodelMover const & rval )
:	//utility::pointer::ReferenceCount(),
	Super( rval ),

	//manager_( rval.manager_ ),
	design_info_( rval.design_info_ ),
	//use_fullmer_( rval.use_fullmer_ ),
	//use_sequence_bias_( rval.use_sequence_bias_ ),
	max_linear_chainbreak_( rval.max_linear_chainbreak_ ),
	centroid_loop_mover_str_( rval.centroid_loop_mover_str_ ),
	redesign_loop_neighborhood_( rval.redesign_loop_neighborhood_ ),
	//resfile_( rval.resfile_ ),
	dr_cycles_( rval.dr_cycles_ ),
	centroid_sfx_( rval.centroid_sfx_->clone() ),
	fullatom_sfx_( rval.fullatom_sfx_->clone() ),
	blueprint_( rval.blueprint_ ),

	// User specified options
	op_user_remodel_vdw_(rval.op_user_remodel_vdw_),
	op_user_remodel_rama_(rval.op_user_remodel_rama_),
	op_user_remodel_cbeta_(rval.op_user_remodel_cbeta_),
	op_user_remodel_cenpack_(rval.op_user_remodel_cenpack_),
	op_user_remodel_hb_lrbb_(rval.op_user_remodel_hb_lrbb_),
	op_user_remodel_hb_srbb_(rval.op_user_remodel_hb_srbb_),
	op_user_remodel_rsigma_(rval.op_user_remodel_rsigma_),
	op_user_remodel_ss_pair_(rval.op_user_remodel_ss_pair_),
	op_user_domainFusion_insert_segment_from_pdb_(rval.op_user_domainFusion_insert_segment_from_pdb_),
	op_user_remodel_repeat_structure_(rval.op_user_remodel_repeat_structure_),
	op_user_symmetry_symmetry_definition_(rval.op_user_symmetry_symmetry_definition_),
	op_user_remodel_quick_and_dirty_(rval.op_user_remodel_quick_and_dirty_),
	op_user_enzdes_cstfile_(rval.op_user_enzdes_cstfile_),
	op_user_constraints_cst_file_(rval.op_user_constraints_cst_file_),
	op_user_remodel_build_disulf_(rval.op_user_remodel_build_disulf_),
	op_user_remodel_rg_(rval.op_user_remodel_rg_),
	op_user_remodel_check_scored_centroid_(rval.op_user_remodel_check_scored_centroid_),
	op_user_remodel_run_confirmation_(rval.op_user_remodel_run_confirmation_),
	op_user_remodel_free_relax_(rval.op_user_remodel_free_relax_),
	op_user_RemodelLoopMover_cyclic_peptide_(rval.op_user_RemodelLoopMover_cyclic_peptide_),
	op_user_remodel_swap_refine_confirm_protocols_(rval.op_user_remodel_swap_refine_confirm_protocols_),
	op_user_design_allow_rare_aro_chi_(rval.op_user_design_allow_rare_aro_chi_),
	op_RemodelLoopMover_max_linear_chainbreak_(rval.op_RemodelLoopMover_max_linear_chainbreak_),
	op_remodel_dr_cycles_(rval.op_remodel_dr_cycles_),
	op_remodel_cen_sfxn_(rval.op_remodel_cen_sfxn_),
	op_remodel_vdw_(rval.op_remodel_vdw_),
	op_remodel_rama_(rval.op_remodel_rama_),
	op_remodel_cbeta_(rval.op_remodel_cbeta_),
	op_remodel_cenpack_(rval.op_remodel_cenpack_),
	op_remodel_hb_lrbb_(rval.op_remodel_hb_lrbb_),
	op_remodel_hb_srbb_(rval.op_remodel_hb_srbb_),
	op_remodel_rsigma_(rval.op_remodel_rsigma_),
	op_remodel_ss_pair_(rval.op_remodel_ss_pair_),
	op_run_show_simulation_in_pymol_(rval.op_run_show_simulation_in_pymol_),
	op_remodel_num_trajectory_(rval.op_remodel_num_trajectory_),
	op_remodel_repeat_structure_(rval.op_remodel_repeat_structure_),
	op_remodel_cen_minimize_(rval.op_remodel_cen_minimize_),
	op_remodel_cstfilter_(rval.op_remodel_cstfilter_),
	op_remodel_collect_clustered_top_(rval.op_remodel_collect_clustered_top_),
	op_out_prefix_(rval.op_out_prefix_),
	op_remodel_use_same_length_fragments_(rval.op_remodel_use_same_length_fragments_),
	op_packing_soft_rep_design_(rval.op_packing_soft_rep_design_),
	op_remodel_rg_(rval.op_remodel_rg_),
	op_remodel_blueprint_(rval.op_remodel_blueprint_),  // File
	op_remodel_checkpoint_(rval.op_remodel_checkpoint_),
	op_remodel_use_pose_relax_(rval.op_remodel_use_pose_relax_),
	op_remodel_use_cart_relax_(rval.op_remodel_use_cart_relax_),
	op_remodel_use_blueprint_sequence_(rval.op_remodel_use_blueprint_sequence_),
	op_remodel_bypass_fragments_(rval.op_remodel_bypass_fragments_),
	op_remodel_free_relax_(rval.op_remodel_free_relax_)
{
	if ( rval.vlb_.get() ) {
		vlb_ = new VarLengthBuild( *rval.vlb_ );
	}
}

void RemodelMover::set_param_from_options()
{
	// Options that don't have default values.  Test before access.
	if( op_user_packing_soft_rep_design_ )
		op_packing_soft_rep_design_ = option[OptionKeys::packing::soft_rep_design];
	if( op_user_remodel_rg_ )
		op_remodel_rg_ = option[OptionKeys::remodel::rg];
	if( op_user_remodel_blueprint_ )
		op_remodel_blueprint_ = option[OptionKeys::remodel::blueprint]();  // File
	if( op_user_remodel_checkpoint_ )
		op_remodel_checkpoint_ = option[OptionKeys::remodel::checkpoint];
	if( op_user_remodel_use_pose_relax_ )
		op_remodel_use_pose_relax_ = option[OptionKeys::remodel::use_pose_relax];
	if( op_user_remodel_use_cart_relax_ )
		op_remodel_use_cart_relax_= option[OptionKeys::remodel::use_cart_relax];
	if( option.has(OptionKeys::remodel::use_blueprint_sequence) )
		op_remodel_use_blueprint_sequence_ = option[OptionKeys::remodel::use_blueprint_sequence];
	if( option.has(OptionKeys::remodel::bypass_fragments) )
		op_remodel_bypass_fragments_ = option[OptionKeys::remodel::bypass_fragments];
	if( option.has(OptionKeys::remodel::free_relax) )
		op_remodel_free_relax_ = option[OptionKeys::remodel::free_relax];
}
///
/// @begin RemodelMover::~RemodelMover
///
/// @brief
/// Default destructor. Does this need to free the VarLengthBuild memory?
///
RemodelMover::~RemodelMover() {}


///
/// @begin RemodelMover::register_user_options
///
/// @brief
/// Checks for presence of any score term weight override options and calls set_weight on the centroid scorefunction.
/// Note: options only get applied to centroid scorefunction - fullatom scorefunction left as is.
///
void RemodelMover::register_user_options() {

	// set optional weights
	if ( op_user_remodel_vdw_ ) {
		centroid_sfx_->set_weight( vdw, op_remodel_vdw_ );
		TR << "USER OVERWRITE VDW: " << op_remodel_vdw_ << std::endl;
	}

	if ( op_user_remodel_rama_ ){
		centroid_sfx_->set_weight( rama, op_remodel_rama_ );
		TR << "USER OVERWRITE RAMA: " << op_remodel_rama_ << std::endl;
	}

	if ( op_user_remodel_cbeta_ ){
		centroid_sfx_->set_weight( cbeta, op_remodel_cbeta_ );
		TR << "USER OVERWRITE CBETA: " << op_remodel_cbeta_ << std::endl;
	}

	if ( op_user_remodel_cenpack_ ){
		centroid_sfx_->set_weight( cenpack, op_remodel_cenpack_ );
		TR << "USER OVERWRITE CENPACK: " << op_remodel_cenpack_ << std::endl;
	}

	if ( op_user_remodel_hb_lrbb_ ){
		centroid_sfx_->set_weight( hbond_lr_bb, op_remodel_hb_lrbb_ );
		TR << "USER OVERWRITE HB_LRBB: " << op_remodel_hb_lrbb_ << std::endl;
	}

	if ( op_user_remodel_hb_srbb_ ){
		centroid_sfx_->set_weight( hbond_sr_bb, op_remodel_hb_srbb_ );
		TR << "USER OVERWRITE HB_SRBB: " << op_remodel_hb_srbb_ << std::endl;
	}
	if ( op_user_remodel_rg_ ){
		centroid_sfx_->set_weight( rg, op_remodel_rg_ );
		TR << "USER OVERWRITE RG: " << op_remodel_rg_ << std::endl;
	}

	if ( op_user_remodel_rsigma_ ){
		centroid_sfx_->set_weight( rsigma, op_remodel_rsigma_ );
		TR << "USER OVERWRITE RSIGMA: " << op_remodel_rsigma_ << std::endl;
	}
	if ( op_user_remodel_ss_pair_ ){
		centroid_sfx_->set_weight( ss_pair, op_remodel_ss_pair_ );
		TR << "USER OVERWRITE SSPAIR: " << op_remodel_ss_pair_ << std::endl;
	}

}


/// @brief clone for parser
RemodelMover::MoverOP RemodelMover::clone() const {
	return new RemodelMover( *this );
}


/// @brief fresh instance for parser
RemodelMover::MoverOP RemodelMover::fresh_instance() const {
	return new RemodelMover();
}


/// @brief clone this object
RemodelMover::MoverOP RemodelMover::clone() {
	return new RemodelMover( *this );
}


/// @brief create this type of object
RemodelMover::MoverOP RemodelMover::fresh_instance() {
	return new RemodelMover();
}


/// @brief the centroid level score function, default "remodel_cen"
ScoreFunction const & RemodelMover::centroid_scorefunction() const {
	return *centroid_sfx_;
}


/// @brief the full-atom level score function, default score12
ScoreFunction const & RemodelMover::fullatom_scorefunction() const {
	return *fullatom_sfx_;
}


/// @brief set the centroid level score function
void RemodelMover::centroid_scorefunction( ScoreFunction const & sfx ) {
	centroid_sfx_ = sfx.clone();
}


/// @brief set the centroid level score function
void RemodelMover::centroid_scorefunction( ScoreFunctionOP const & sfx ) {
	centroid_sfx_ = sfx->clone();
}


/// @brief set the full-atom level score function
void RemodelMover::fullatom_scorefunction( ScoreFunction const & sfx ) {
	fullatom_sfx_ = sfx.clone();
}


/// @brief set the full-atom level score function
void RemodelMover::fullatom_scorefunction( ScoreFunctionOP const & sfx ) {
	fullatom_sfx_ = sfx->clone();
}


///
/// @begin RemodelMover::apply
///
/// @brief
/// Apply method for Mover.
/// Checks the values of the following options
/// -remodel::checkpoint
/// -remodel::domainFusion::insert_segment_from_pdb
/// -remodel::bypass_fragments,
/// -remodel::num_trajectory,
/// -remodel::repeat_structure
/// -remodel::build_disulf
/// -remodel::quick_and_dirty
/// -remodel::use_pose_relax
/// -remodel::run_confirmation
/// -enzdes::cstfile
/// -symmetry::symmetry_definition
/// -run::show_simulation_in_pymol
///
///
void RemodelMover::apply( Pose & pose ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace protocols;
	using namespace chemical;

	using core::pose::metrics::CalculatorFactory;
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_RETRY;
	using protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator;
	using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;

#if defined GL_GRAPHICS
	protocols::viewer::add_conformation_viewer( pose.conformation(), "Remodel" );
#endif

	TR << "apply(): entered RemodelMover apply(). pose.total_residue(): " << pose.total_residue() << std::endl;

	// store the starting pose for KIC confirmation RMSD calculation
	native_pose_ = pose;

	// assign secondary structure
	scoring::dssp::Dssp dssp( pose );

	forge::remodel::RemodelData remodel_data;
	forge::remodel::RemodelWorkingSet working_model;

	// read blueprint file
	TR << "apply(): reading blueprint file " << std::endl;

	//TR << "blueprint value 0 is " << blueprint_ << std::endl;

	if (blueprint_ == "") {
		//TR << "blueprint value 1 is " << blueprint_ << std::endl;
		blueprint_ = op_remodel_blueprint_;
		//TR << "blueprint value 2 is " << blueprint_ << std::endl;
	}

	remodel_data.getLoopsToBuildFromFile(blueprint_);

	TR << pose.total_residue() << std::endl;

	ObjexxFCL::FArray1D_char dsspSS( pose.total_residue() );
	dssp.dssp_reduced(dsspSS);

	TR << "apply(): input PDB dssp assignment: (based on start structure)" << std::endl;
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		TR << dsspSS(i);
	}
	TR << std::endl;

	remodel_data.updateWithDsspAssignment( dsspSS );
	dssp.insert_ss_into_pose( pose );

	// process domain insertion option
	if ( op_user_domainFusion_insert_segment_from_pdb_ ) {
		TR << "apply(): INSERT SEGMENT FROM PDB" << std::endl;
		remodel_data.collectInsertionPose();
		// remodel_data will have several class member variables updated with the insertion information at this point
	}

	// create a scorefunction and score the pose first
	scoring::ScoreFunctionOP sfx = scoring::ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );
	(*sfx)( pose );

	working_model.workingSetGen( pose, remodel_data );

	remodel_data_ = remodel_data; // will use the movemap for natro definition later
	working_model_ = working_model;

	// test PyMol viewer
	if ( op_run_show_simulation_in_pymol_ ) {
		moves::AddPyMolObserver( pose, false, core::Real( 0.50 ) );
	}



	/*
	// DEBUG
	std::set<core::Size> up = working_model.manager.undefined_positions();
	for ( std::set<core::Size>::iterator i = up.begin(); i!=up.end(); i++) {
		TR << *i << std::endl;
	}
	std::set<core::Size> uup = working_model.manager.union_of_intervals_containing_undefined_positions();
	for ( std::set<core::Size>::iterator i = uup.begin(); i!=uup.end(); i++){
		TR << *i <<  " UUP" <<  std::endl;
	}
	*/

	if (op_user_remodel_repeat_structure_) {
		//for cases involve jxn, need to make pose longer so manager won't complain
		//about missing residues

		//this is pre modify, so simply extend to 2x blueprint length,
		//with extensions, the pose will go beyond the correct length.  need to fix
		//that after modify. Residues beyond first copy+ jxn doesn't really matter
		if (pose.total_residue() < 2*remodel_data.sequence.length()){ //just making sure it's shorter before grow, input pose can be longer
			Size len_diff = (2*remodel_data_.sequence.length()) - pose.total_residue();
			// append a tail of the same length
			for (Size i = 1; i<= len_diff; ++i){
				core::chemical::ResidueTypeSet const & rsd_set = (pose.residue(1).residue_type_set());
				core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd_set.name_map("ALA") ) );
				pose.conformation().safely_append_polymer_residue_after_seqpos(* new_rsd,pose.total_residue(), true);
				pose.conformation().insert_ideal_geometry_at_polymer_bond(pose.total_residue()-1);
				pose.set_omega(pose.total_residue()-1,180);
			}
		}
	}

	//
	//	Pose testArc;
	//	testArc = pose;
	if (working_model.manager.size()!= 0){
		if (!op_remodel_bypass_fragments_){
			working_model.manager.modify(pose);
			}else{
			working_model.manager.dummy_modify(pose.total_residue());
		}
		//	protocols::forge::methods::restore_residues(working_model.manager.original2modified(), testArc, pose);
		//	pose.dump_pdb("testArcRestore.pdb");
		//testArc=pose;
		manager_ = working_model.manager;
		core::pose::renumber_pdbinfo_based_on_conf_chains(
			pose,
			true ,  // fix chain
			true, // start_from_existing_numbering
			false, // keep_insertion_code
			false // rotate_chain_id
		);
	}

	//finally recheck length to ensure blueprint compliance
	if (op_user_remodel_repeat_structure_) {
		Size max_pdb_index = remodel_data_.blueprint.size()*2;
		while (pose.total_residue() >= max_pdb_index){
			pose.conformation().delete_residue_slow(pose.total_residue());
		}

		if ( pose.total_residue() < (remodel_data_.sequence.length()*2) ) {
			Size len_diff = (2*remodel_data_.sequence.length()) - pose.total_residue();
			// append a tail of the same length
			for (Size i = 1; i<= len_diff; ++i){
				core::chemical::ResidueTypeSet const & rsd_set = (pose.residue(1).residue_type_set());
				core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd_set.name_map("ALA") ) );
				pose.conformation().safely_append_polymer_residue_after_seqpos(* new_rsd,pose.total_residue(), true);
				pose.conformation().insert_ideal_geometry_at_polymer_bond(pose.total_residue()-1);
				pose.set_omega(pose.total_residue()-1,180);
			}
		}

	}


	/*
	up = working_model.manager.undefined_positions();
	for ( std::set<core::Size>::iterator i = up.begin(); i!=up.end(); i++ ) {
		TR << *i << std::endl;
	}
	uup = working_model.manager.union_of_intervals_containing_undefined_positions();
	for ( std::set<core::Size>::iterator i = uup.begin(); i!=uup.end(); i++){
		TR << *i <<  " UUP2" <<  std::endl;
	}
	*/

	//manager_.dummy_modify(testArc.n_residue());
	//core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, true);
	//core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD, true);
	//protocols::forge::methods::restore_residues(manager_.original2modified(), testArc, pose);
	//pose.dump_pdb("testArcRestore2.pdb");
	//protocols::forge::methods::restore_residues(manager_.original2modified(), testArc, pose);
	//pose.update_residue_neighbors();
	//pose.dump_pdb("testArcRestore3.pdb");
	//testArc.dump_pdb("testArcRestoreSrc3.pdb");
	//protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD);
	//protocols::simple_moves::ReturnSidechainMover recover_sidechains( testArc);
	//to_all_atom.apply(pose);
	//recover_sidechains.apply(pose);
	//pose.dump_pdb("MoverREstore.pdb");

	// initialize symmetry

	// only symmetrize here if not in the repeat structure mode. for repeats, stay monomer until repeat generation.
	if ( op_user_symmetry_symmetry_definition_ && !op_remodel_repeat_structure_ )  {
		simple_moves::symmetry::SetupForSymmetryMover pre_mover;
		pre_mover.apply( pose );
		// Remodel assumes chain ID is ' '
		//pose::PDBInfoOP pdb_info ( pose.pdb_info() );
		//for ( Size i=1; i<= pdb_info->nres(); ++i ){
		//	pdb_info->chain(i,' ');
		//}
		//pose.pdb_info( pdb_info );
	}

	Size i = op_remodel_num_trajectory_;
	Size num_traj = i; // need this for checkpointing math
	Size prev_checkpoint = 0;

	forge::remodel::RemodelAccumulator accumulator( working_model );

	if ( op_remodel_checkpoint_ ) {
		prev_checkpoint = accumulator.recover_checkpoint();
		if ( prev_checkpoint >= i ) {
			i = 0;
		} else {
			i = i - prev_checkpoint;
		}
	}

	if ( working_model.manager.size() != 0 ) {
		// setup calculators
		pose::metrics::CalculatorFactory::Instance().remove_calculator( neighborhood_calc_name() );
		pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc_name(),
			new toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( manager_.union_of_intervals_containing_undefined_positions() ) );
	}

	/*
	up = working_model.manager.undefined_positions();
	for ( std::set<core::Size>::iterator i = up.begin(); i!=up.end(); i++){
		TR << *i << std::endl;
	}
	uup = working_model.manager.union_of_intervals_containing_undefined_positions();
	for ( std::set<core::Size>::iterator i = uup.begin(); i!=uup.end(); i++){
		TR << *i <<  " UUP2" <<  std::endl;
	}
	*/

	if ( op_user_remodel_repeat_structure_ ) {

		// turning on the res_type_linking constraint weight for designs
		fullatom_sfx_->set_weight( scoring::atom_pair_constraint, 1.0 );
		fullatom_sfx_->set_weight( scoring::res_type_linking_constraint, 0.3 );
		fullatom_sfx_->set_weight( scoring::res_type_constraint, 1.0 );

		//fullatom_sfx_->set_weight( core::scoring::fa_dun, 0 );
		//fullatom_sfx_->set_weight( core::scoring::fa_sol, 0 );
		//fullatom_sfx_->set_weight( core::scoring::fa_pair, 0 );
		//fullatom_sfx_->set_weight( core::scoring::hbond_sc, 0 );
		//fullatom_sfx_->set_weight( core::scoring::hbond_sc, 0 );
		//fullatom_sfx_->set_weight( core::scoring::fa_intra_rep, 0 );
		//fullatom_sfx_->set_weight( core::scoring::rama, 0 );
		//fullatom_sfx_->set_weight( core::scoring::hbond_bb_sc, 0 );
		//fullatom_sfx_->set_weight( core::scoring::p_aa_pp, 0 );
		//fullatom_sfx_->set_weight( core::scoring::ref, 10 );
	}

	// initializes a RemodelDesignMover which will be used in the loop below
	forge::remodel::RemodelDesignMover designMover( remodel_data, working_model, fullatom_sfx_ );

	Size no_attempts_at_centroid_build = 0;
	//Size no_attempts_to_make_at_centroid_build = 10;

	//bool quick_mode = op_user_remodel_quick_and_dirty_;

	Size repeat_number = op_remodel_repeat_structure_;
	//rerooting tree
	if (op_user_remodel_repeat_structure_ && pose.total_residue() == remodel_data.blueprint.size()*repeat_number) {
		core::kinematics::FoldTree f = pose.fold_tree();
		f.reorder(working_model.safe_root_);
		pose.fold_tree(f);
		TR << "rerooting tree: " << pose.fold_tree() << std::endl;
	}

	while ( i > 0 ) {

		// cache the modified pose first for REPEAT
		Pose cached_modified_pose( pose );

		// do centroid build
		TR << std::endl << "apply(): BUILD CYCLE REMAINING " << i << std::endl;
		kinematics::FoldTree originalTree = pose.fold_tree();
		TR << "ORIGINAL TREE: " << pose.fold_tree() << std::endl;
		if ( working_model.manager.size() != 0 ) {
			if ( !centroid_build( pose, working_model.manager ) ) { // build failed
				no_attempts_at_centroid_build++;
				TR << "apply(): number of attempts at loop closure made: " << no_attempts_at_centroid_build << std::endl;
				set_last_move_status( FAIL_RETRY );

				/*if ( !quick_mode ) {
					continue;
				} else {
					if ( no_attempts_at_centroid_build >= no_attempts_to_make_at_centroid_build ) {
						TR << "apply(): number of attempts at loop closure exceeded limit of " << no_attempts_to_make_at_centroid_build << ". quitting." << std::endl;
						set_last_move_status( FAIL_DO_NOT_RETRY );
						return;
					} else {
						// try again. omitting this continue causes the protocol to give up after one failed iteration.
						continue;
					}
				}*/
				i--;
				continue;
				//return;
			}
		}
		if ( op_user_remodel_repeat_structure_ ) {
			// should fold this pose to match just the first segment of a repeat, and that will be used for next round of building
			for ( Size res = 1; res <= cached_modified_pose.n_residue(); res++ ) {
				cached_modified_pose.set_phi( res, pose.phi(res) );
				cached_modified_pose.set_psi( res, pose.psi(res) );
				cached_modified_pose.set_omega( res, pose.omega(res) );
				ResidueType const & rsd_type(pose.residue_type(res));
				replace_pose_residue_copying_existing_coordinates(cached_modified_pose,res,rsd_type);
			}
		}
		core::pose::renumber_pdbinfo_based_on_conf_chains(
				pose,
		    true ,  // fix chain
			  true, // start_from_existing_numbering
			  false, // keep_insertion_code
		    false // rotate_chain_id
		 );

		//test
		//pose.dump_pdb("check.pdb");
		/*
		//extract the constraint currently in Pose for later Recycling
		ConstraintSetOP cst_set_post_built;
		if ( op_user_remodel_repeat_structure_ ) {

			// at this stage it should hold generic cstfile and res_type_linking
			// constraints
			cst_set_post_built = new ConstraintSet( *pose.constraint_set());
		}
		*/

		designMover.set_state("stage");

		// handle constraints as soon as centroid is done.  If applying sidechain
		// constraints, replace residue to the right AA right away.
		if ( op_user_enzdes_cstfile_ ) {
			TR << "apply(): constraint file found on command line. updating score functions to include constraint terms." << std::endl;

			forge::remodel::RemodelEnzdesCstModuleOP cstOP = new forge::remodel::RemodelEnzdesCstModule(remodel_data);

			// RemodelEnzdesCstModule cst(remodel_data);
			//safety
			pose.remove_constraints();
			//wipe out cst_cache
			toolbox::match_enzdes_util::get_enzdes_observer( pose ) -> set_cst_cache( NULL );
			//wipe out observer too
			pose.observer_cache().set( pose::datacache::CacheableObserverType::ENZDES_OBSERVER, NULL , false);

			//cstOP->remove_constraints_from_pose(pose,true /*keep covalent*/, false /*fail if missing*/);

			cstOP->use_all_blocks();
			TR << "apply(): calling RemodelEnzdesCstModule apply function." << std::endl;
			cstOP->apply(pose);
			cstOP->enable_constraint_scoreterms(fullatom_sfx_);
			designMover.scorefunction(fullatom_sfx_);
		}
		if (op_user_constraints_cst_file_){
				//safety
				pose.remove_constraints();

				protocols::simple_moves::ConstraintSetMoverOP constraint = new protocols::simple_moves::ConstraintSetMover();
				constraint->apply( pose );

				fullatom_sfx_->set_weight(core::scoring::atom_pair_constraint, 1.0);
				fullatom_sfx_->set_weight(core::scoring::dihedral_constraint, 10.0);
				designMover.scorefunction(fullatom_sfx_);
		}

		if (op_user_remodel_build_disulf_){
			TR << "apply(): build_disfulf option set. finding disulfides." << std::endl;
			utility::vector1<std::pair <Size, Size> > disulf_partners;
			bool disulfPass = false;
			disulfPass = designMover.find_disulfides_in_the_neighborhood( pose, disulf_partners );
			if ( disulfPass != true ) {
				i--; //for now control disulf with num_trajectory flag, too.
				continue;
			}

			for ( utility::vector1< std::pair< Size, Size > >::iterator itr = disulf_partners.begin(); itr != disulf_partners.end(); itr++ ) {
				Pose disulf_copy_pose = pose;
				utility::vector1< std::pair< Size, Size > > single_disulf;
				single_disulf.push_back(*itr);

				kinematics::MoveMapOP combined_mm = new kinematics::MoveMap;

				combined_mm->import( remodel_data_.natro_movemap_ );
				combined_mm->import( manager_.movemap() );

				designMover.make_disulfide( disulf_copy_pose, single_disulf, combined_mm );
				designMover.apply( disulf_copy_pose );

				// for now, accept all disulf build, as it is hard enough to do already.  Accept instead of cst filter?
				// accumulator.apply(disulf_copy_pose);
				if ( op_user_enzdes_cstfile_ ) {
					simple_filters::ScoreTypeFilter const pose_constraint( fullatom_sfx_, atom_pair_constraint, 10 );
					bool CScore(pose_constraint.apply( disulf_copy_pose ));
					if (!CScore){  // if didn't pass, rebuild
						continue;
					} else {
						accumulator.apply(disulf_copy_pose);
					}

				} else {
					accumulator.apply(disulf_copy_pose);
				}
			}

		} else {
			// option build_disulf not specified...
			//if ( op_user_remodel_repeat_structure_ || op_remodel_cen_minimize_ ) {
			if (op_remodel_cen_minimize_) {

				//cache current foldTree;
				kinematics::FoldTree cenFT = pose.fold_tree();
				pose.fold_tree(originalTree);

				kinematics::MoveMapOP cmmop = new kinematics::MoveMap;
				//pose.dump_pdb("pretest.pdb");

				cmmop->import( remodel_data_.natro_movemap_ );
				cmmop->import( manager_.movemap() );

				for (Size i = 1; i<= pose.total_residue(); ++i){
					std::cout << "bb at " << i << " " << cmmop->get_bb(i) << std::endl;
				}

				//adding angles and bonds dof
				//	cmmop->set(core::id::THETA, true);
				//	cmmop->set(core::id::D, true);

				 for(Size i = 1; i <= pose.n_residue(); i++) {
						for(Size j = 1; j <= pose.residue(i).nheavyatoms(); j++) {
							if (cmmop->get_bb(i) == 1){
								cmmop->set(core::id::DOF_ID(core::id::AtomID(j,i),core::id::THETA),true);
								cmmop->set(core::id::DOF_ID(core::id::AtomID(j,i),core::id::D),true);
							}
						}
					}


				//scoring::ScoreFunctionOP cen_min_sfxn = scoring::ScoreFunctionFactory::create_score_function("score4_smooth");

				TR << "centroid minimizing" << std::endl;
				pose::Pose archived_pose = pose;

				// flip residue type set for centroid minimize
				util::switch_to_residue_type_set( pose, chemical::CENTROID, true );

				protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;

				if (op_user_remodel_repeat_structure_){
							//Dihedral (NCS) Constraints, need to be updated each mutation cycle for sidechain symmetry

							Size repeat_number = op_remodel_repeat_structure_;
							Size segment_length = (pose.n_residue())/repeat_number;


							for (Size rep = 1; rep < repeat_number-1; rep++){ // from 1 since first segment don't need self-linking
								std::stringstream templateRangeSS;
								templateRangeSS << "2-" << segment_length+1; // offset by one to work around the termini
								std::stringstream targetSS;
								targetSS << 1+(segment_length*rep)+1 << "-" << segment_length + (segment_length*rep)+1;
								TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
								setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
							}

						for (Size rep = 1; rep < repeat_number-1; rep++){ // from 1 since first segment don't need self-linking
								std::stringstream templateRangeSS;
								templateRangeSS << "3-" << segment_length+2; // offset by one to work around the termini
								std::stringstream targetSS;
								targetSS << 1+(segment_length*rep)+2 << "-" << segment_length + (segment_length*rep)+2;
								TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
								setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
							}


							std::stringstream templateRangeSS;
							//take care of the terminal repeat, since the numbers are offset.
							templateRangeSS << "2-" << segment_length-1; // offset by one to work around the termini
							std::stringstream targetSS;
							targetSS << 1+(segment_length*(repeat_number-1))+1 << "-" << segment_length + (segment_length*(repeat_number-1))-1;
							TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
							setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
							setup_ncs.apply(pose);

				}

				//sfx->show(TR, pose);
				//TR << std::endl;

				centroid_sfx_->set_weight( core::scoring::atom_pair_constraint, 1.0);
				centroid_sfx_->set_weight(core::scoring::dihedral_constraint, 10.0 );
				//enable cartesian bond terms
				centroid_sfx_->set_weight(core::scoring::cart_bonded_angle,  0.1 );
				centroid_sfx_->set_weight(core::scoring::cart_bonded_length,  0.1 );
				centroid_sfx_->set_weight(core::scoring::cart_bonded_torsion,  0.1 );
				centroid_sfx_->set_weight(core::scoring::omega, 0.2 );

				//only use smooth hb if either of the term is used in centroid build level
				if (centroid_sfx_->get_weight(core::scoring::hbond_lr_bb) > 0 || centroid_sfx_->get_weight(core::scoring::hbond_sr_bb) > 0 ){
					centroid_sfx_->set_weight( core::scoring::cen_hb, 2.0);
				}
				/*
				if (centroid_sfx_->get_weight(core::scoring::env) > 0 ){
					centroid_sfx_->set_weight( core::scoring::cen_env_smooth, centroid_sfx_->get_weight(core::scoring::env));
					centroid_sfx_->set_weight( core::scoring::env, 0.0);
				}*/

				//simple_moves::MinMoverOP minMover = new simple_moves::MinMover( cmmop , centroid_sfx_, "dfpmin_armijo", 0.01, true);
				simple_moves::MinMoverOP minMover = new simple_moves::MinMover( cmmop , centroid_sfx_, "lbfgs_armijo", 0.01, true);
				TR << "cen_minimize pose foldtree: " << pose.fold_tree() << std::endl;
				minMover->apply(pose);

				//reset cen_hb to 0
				centroid_sfx_->set_weight( core::scoring::cen_hb, 0.0);
				//switch back the foldtree
				//pose.fold_tree(cenFT);

				// flip residue type set back, for repeat builds, currently don't do
				// restore_sidechain, as they should all be redesigned.  MAY NEED TO
				// CHANGE
				util::switch_to_residue_type_set( pose, chemical::FA_STANDARD, true );
				//forge::methods::restore_residues( manager_.original2modified(), archived_pose , pose );
				//pose.dump_pdb("test.pdb");

			}
			TR << "apply(): calling RemodelDesignMover apply function." << std::endl;
			designMover.apply(pose);
			//****Previously this loop didn't maintain the pose correctly which resulted in a difficult to track down seg fault.  Fixed, but it's questionable weather you would want to filter poses based on constraints.
			if ( op_user_enzdes_cstfile_ ||
					 op_user_constraints_cst_file_
					 ){
				simple_filters::ScoreTypeFilter const  pose_constraint( fullatom_sfx_, atom_pair_constraint, op_remodel_cstfilter_ );
				bool CScore(pose_constraint.apply( pose ));
				if (!CScore){  // if didn't pass, rebuild
					TR << "built model did not pass constraints test." << std::endl;
					if ( op_user_remodel_repeat_structure_ ) {
						//reset the pose to monomer
						pose = cached_modified_pose;
					} else {
						pose.fold_tree(originalTree);
					}
					continue;
				}
				else {
					accumulator.apply(pose);
				}
			} else {
				accumulator.apply(pose);
			}
			//*****END horrible code.
		}
		if ( op_remodel_checkpoint_ ) {
			// debug:
			TR << "writing chkpnt at step " << num_traj-i+prev_checkpoint << std::endl;
			accumulator.write_checkpoint(num_traj-i-prev_checkpoint);
		}
		// restore foldtree
		if(i > 1) {//messes up rosetta_scripts if done on last loop.
			if ( op_user_remodel_repeat_structure_ ) {
				//reset the pose to monomer
				pose = cached_modified_pose;
			} else {
				pose.fold_tree(originalTree);
			}
		}
		i--; // 'i' is the number of remaining trajectories

	}
	/* DONT USE THIS FOR NOW...
	if (get_last_move_status() == FAIL_RETRY){
		return;
	}
*/
	// take the lowest member and the cluster center
	// accumulator.shrink_cluster();
	std::vector< pose::PoseOP > results;
	if ( accumulator.cluster_switch() ) {
		results = accumulator.clustered_best_poses();
		//results = accumulator.clustered_top_poses(op_remodel_collect_clustered_top_);
	} else {
		results = accumulator.contents_in_pose_store();
	}

	// seriously refine the poses
	Size filecount = 1;
	core::Real current_score = 100000;

	TR << "clustered poses count: " << results.size() << std::endl;
	for ( std::vector< pose::PoseOP >::iterator it = results.begin(), end= results.end(); it!= end; it++ ) {
		bool bypass_refinement = op_user_remodel_quick_and_dirty_;

		if ( working_model.manager.size() == 0 )
			bypass_refinement = true;

		if ( !bypass_refinement ) {

			//std::stringstream SS1;
			//SS1 << "pre-ref_" << filecount << ".pdb";
			//(*(*it)).dump_scored_pdb(SS1.str(), *fullatom_sfx_);

			TR << "aggressively refine" << std::endl;
			if ( op_remodel_use_pose_relax_ ) {
				if ( !design_refine_seq_relax( *(*it), designMover ) ) {
					TR << "WARNING: DESIGN REFINE SEQ RELAX FAILED!! (one should never see this)" << std::endl;
					continue;
				}
			}
			else if (op_remodel_use_cart_relax_){
				if (!design_refine_cart_relax(*(*it), designMover)){
					TR << "WARNING: CARTESIAN MIN FAILED!! (one should never see this)" << std::endl;
					continue;
				}
			} else {
				if (! design_refine(*(*it), designMover)){
					TR << "WARNING: DESIGN REFINE FAILED TO CLOSE STRUCTURE!!" << std::endl;
					continue;
				}
			}
		} else // simple design
		{
			if (op_user_remodel_check_scored_centroid_){
					std::stringstream SS;
					std::string prefix =  op_out_prefix_;
					if (!prefix.empty()){
						SS << prefix << "_" << filecount << "_cen.pdb";
					}
					else {
						SS << filecount << "_cen.pdb";
					}
					core::util::switch_to_residue_type_set( *(*it), core::chemical::CENTROID, true);
					(*(*it)).dump_scored_pdb(SS.str(), *centroid_sfx_);
			}

			designMover.set_state("finish");
			designMover.apply(*(*it));

		}

		if ( op_user_remodel_run_confirmation_ ) {
			if ( !confirm_sequence(*(*it)) ) {
				TR << "WARNING: STRUCTURE DID NOT PASS KIC CONFIRMATION!!" << std::endl;
				continue;
			}
		}

		std::stringstream SS;
		std::string prefix =  op_out_prefix_;
		if (!prefix.empty()){
			SS << prefix << "_" << filecount << ".pdb";
		}
		else {
			SS << filecount << ".pdb";
		}
		// this is to make sure that the final scoring is done with SCORE12
		scoring::ScoreFunctionOP scorefxn = scoring::ScoreFunctionFactory::create_score_function( scoring::STANDARD_WTS, scoring::SCORE12_PATCH );

	  if (op_user_remodel_repeat_structure_) {
						//Experiment with RemodelGlobalFrame
						RemodelGlobalFrame RGF(remodel_data, working_model, scorefxn);
						RGF.align_segment(*(*it));
						RGF.apply(*(*it));

		}

		(*(*it)).dump_scored_pdb(SS.str(), *scorefxn);

		simple_filters::ScoreTypeFilter const pose_total_score( scorefxn, total_score, 100 );
		Real score( pose_total_score.compute( *(*it) ) );
		if ( score <= current_score ) {
			current_score = score ;
			pose = *(*it) ;
		}

		filecount++;
	}

	// update PDBinfo
	pose.pdb_info( new core::pose::PDBInfo( pose ));

	// setup calculators
	pose::metrics::CalculatorFactory::Instance().remove_calculator( neighborhood_calc_name() );
	pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc_name(),
		new toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( manager_.union_of_intervals_containing_undefined_positions() ) );

	/*
	// do design-refine iteration
	if ( dr_cycles_ > 0 ) {
		if ( !design_refine( pose ) ) { // design-refine failed
			set_last_move_status( FAIL_RETRY );
			return;
		}
	}
	*/

	// if we've gotten to this point, then the structure has been built properly
	set_last_move_status( MS_SUCCESS );

	// setup the PoseMetricCalculators and add them to the evaluators in the JobOutputter
	pose::metrics::CalculatorFactory::Instance().remove_calculator( loops_buns_polar_calc_name() );
	pose::metrics::CalculatorFactory::Instance().remove_calculator( neighborhood_buns_polar_calc_name() );

	pose::metrics::CalculatorFactory::Instance().register_calculator(
		loops_buns_polar_calc_name(),
		new toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator( "default", "default", manager_.union_of_intervals_containing_undefined_positions() )
	);

	basic::MetricValue< std::set< Size > > loops_neighborhood;
	pose.metric( neighborhood_calc_name(), "neighbors", loops_neighborhood );
	pose::metrics::CalculatorFactory::Instance().register_calculator(
		neighborhood_buns_polar_calc_name(),
		new toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator( "default", "default", loops_neighborhood.value() )
	);

}


/// @brief get_name function for JobDistributor
std::string RemodelMover::get_name() const {
	return "RemodelMover";
}

///
/// @begin RemodelMover::centroid_build
///
/// @brief
/// Does the same as function below, but takes in a BuildManager object.
/// Also checks the value of the option -remodel::bypass_fragments.
///
bool RemodelMover::centroid_build( Pose & pose, protocols::forge::build::BuildManager & manager ) {

	manager_ = manager;
	if ( op_remodel_bypass_fragments_ ) {
		TR << "-=BYPASSING FRAGMENT BUILD (REFINE ONLY) =-" << std::endl;
		return true;
	}

	if ( centroid_build( pose ) ) {
		//update external manager
		//manager = manager_;
		return true;
	} else {
		return false;
	}

}


///
/// @begin RemodelMover::centroid_build
///
/// @brief
/// Runs the centroid level build state. Returns true if loop was closed, false if not.
/// Also checks the value of options -remodel:use_same_length_fragments, -remodel:use_blueprint_sequence
/// and -remodel:repeat_structure
///
bool RemodelMover::centroid_build( Pose & pose ) {

	using namespace basic::options;
	using namespace core::scoring;
	using protocols::moves::MS_SUCCESS;

	using core::util::switch_to_residue_type_set;
	using protocols::forge::methods::restore_residues;
	//using protocols::toolbox::pose_manipulation::construct_poly_uniq_restype_pose;
	using namespace protocols::forge::components;

	// safety, clear the energies object
	pose.energies().clear();

	// make backup Pose for transferring sidechains
	Pose archive_pose = pose;
	Pose modified_archive_pose = archive_pose;
	//manager_.modify( modified_archive_pose );

	// ensure modified_archive_pose is completely full-atom, otherwise mismatch
	// will occur when restoring sidechains at the end of the procedure
	bool mod_ap_is_full_atom = true;
	for ( Size i = 1, ie = modified_archive_pose.n_residue(); mod_ap_is_full_atom && i != ie; ++i ) {
		mod_ap_is_full_atom &= ( modified_archive_pose.residue( i ).residue_type_set().name() == core::chemical::FA_STANDARD );
	}

	if ( !mod_ap_is_full_atom ) {
		core::util::switch_to_residue_type_set( modified_archive_pose, core::chemical::FA_STANDARD );
	}
/*
	// flip to poly-ala-gly-pro-disulf pose, only in the rebuilt segment
	utility::vector1< Size > protein_residues;
	for (std::set<Size>::iterator it=rebuild.begin(), end=rebuild.end(); it != end; it++){
	//for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
		if ( pose.residue( *it ).is_protein() ) {
			protein_residues.push_back( *it );
			TR<< "turning these to ala: " << *it << std::endl;
		}
	}
	TR << "default building restype: " << "ALA" << std::endl;
	construct_poly_uniq_restype_pose( pose, protein_residues, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA"), true, true, true );
*/
	// run VLB to build the new section, if no segments have been added/deleted
	// we use the same VLB so that fragment caching works properly
	if ( !vlb_.get() ) {
		vlb_ = new VarLengthBuild( manager_ , remodel_data_ );
	}
	if (!working_model_.abego.empty()){
		//the following block simply packages the string to feed to vlb
		utility::vector1<std::string> abego_vec;
		for (Size i = 0; i < working_model_.abego.length(); i++){
			std::string buffer;
			buffer.push_back(working_model_.abego[i]);
			abego_vec.push_back(buffer);
		}
		vlb_->set_abego(abego_vec);
	}

	vlb_->scorefunction( centroid_sfx_ );
	vlb_->vall_memory_usage( VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS );
	vlb_->use_fullmer( op_remodel_use_same_length_fragments_ );
	vlb_->max_linear_chainbreak( max_linear_chainbreak_ );
	vlb_->loop_mover_str( centroid_loop_mover_str_ );
	vlb_->restart_mode(true);
	vlb_->new_secondary_structure_override(working_model_.ss);
	if ( op_remodel_use_blueprint_sequence_ ) {
	    if (op_user_remodel_repeat_structure_) {
				Size copies = op_remodel_repeat_structure_;
				String rep_seq = remodel_data_.sequence;
				while (copies > 1){
					rep_seq.append(remodel_data_.sequence);
					copies--;
				}
			  vlb_->new_sequence_override( rep_seq );
			}
			else vlb_->new_sequence_override( remodel_data_.sequence );
	}

	TR << "centroid_build(): calling VariableLengthBuild apply." << std::endl;
	vlb_->apply( pose );

	if ( vlb_->get_last_move_status() == MS_SUCCESS ) {

		// record the used manager w/ all mapping info
		//manager_ = vlb_->manager();

		// safety, clear all the energies before restoring full-atom residues and scoring
		pose.energies().clear();
		bool denovo = true;
		//have to loop to identify denovo case
		for (int i = 0, ie = (int)remodel_data_.blueprint.size(); i < ie; i++){
			if (remodel_data_.blueprint[i].sstype == "."){ //if anywhere hits this assignment, not de novo
				denovo = false;
			}
		}
		if(denovo)
			modified_archive_pose = pose;

		if (op_user_remodel_repeat_structure_) {
		//this part really needs work....  currently doesn't allow growing a loop
		//in regional repeat building.  This section is used in de novo rebuild
		//cases where the monomer pose is extended, so to restore the sidechain,
		//the source of the sidechains has to be extended too in de novo cases.
		//but with refining an existing repeat pose, no need to extend
		//	if (modified_archive_pose.total_residue() == pose.total_residue()){ //dangerous, assuming no further length change
				//do nothing.
		//	} else { // if there's mismatch, assuming restoration source need extension... dangerous.
				// because of code-change for always using 2x modified
				using namespace protocols::loops;
				using protocols::forge::methods::intervals_to_loops;
				std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
				LoopsOP loops = new Loops( intervals_to_loops( loop_intervals.begin(), loop_intervals.end() ) );
				RemodelLoopMover RLM(loops);
				RLM.set_repeat_tail_length(remodel_data_.sequence.length());
				Pose bufferPose(modified_archive_pose);
				//due to code change, modified_archive_pose would always be 2x length now
				RLM.repeat_generation_with_additional_residue( bufferPose, modified_archive_pose );
				if ( op_user_symmetry_symmetry_definition_){
				//symmetrize if both rep+sym are used
					simple_moves::symmetry::SetupForSymmetryMover pre_mover;
					pre_mover.apply(modified_archive_pose);
					modified_archive_pose.pdb_info()->obsolete(true);
				}

		//	}
		}

		// Swap back original sidechains.  At the moment this is a two step process
		// in case any sidechains from SegmentInsert and the like that aren't in the
		// original archive pose need to be transferred.
		restore_residues( modified_archive_pose, pose );

		// since pose is setup modified in RemodelMover, only one step will do
		//restore_residues( manager_.original2modified(), archive_pose, pose );
		// go ahead and score w/ full-atom here; we do this in case there are no
		// design-refine cycles -- it's useful to have e.g. rama in the output
		(*fullatom_sfx_)( pose );

		TR << "centroid_build(): variable length build succeeded." << std::endl;
		if (op_user_remodel_repeat_structure_) {
			//return the modified pose to original state, otherwise it keeps growing.
			modified_archive_pose = archive_pose;
		}

		return true; // loop closed
	} else {
		TR << "centroid_build(): variable length build failed. resetting pose to archived pose." << std::endl;
		// set passed in pose reference to the pose this function was called with?  I guess that means the rebuild tried
		// to happen, failed, and so we don't change anything.
		pose = archive_pose;

	}
	TR << "centroid_build(): centroid_build unable to close loop. retrying." << std::endl;
	return false; // false if loop not closed
}


///
/// @begin RemodelMover::design_refine_seq_relax
///
/// @brief
/// Sets up constraints and a modified scorefunction and run design/relax cycles.
/// Checks the value of -remodel:repeat_structure.
///
bool RemodelMover::design_refine_seq_relax( Pose & pose, RemodelDesignMover & designMover ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace protocols;

	// collect new regions/positions
	std::set< forge::build::Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	forge::build::BuildManager::Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	loops::Loops loops = forge::methods::intervals_to_loops( loop_intervals.begin(), loop_intervals.end() );

  if (op_user_remodel_repeat_structure_ || op_remodel_free_relax_ ){
		//do nothing
	} else {
		protocols::forge::methods::fill_non_loop_cst_set(pose, loops);
	}

	// safety, clear the energies object
	pose.energies().clear();

	// for refinement, always use standard repulsive
	ScoreFunctionOP sfx = core::scoring::getScoreFunction();

	// turning on weights
	sfx->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfx->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	sfx->set_weight( core::scoring::angle_constraint, 1.0 );
	sfx->set_weight( core::scoring::dihedral_constraint, 10.0 ); // 1.0 originally
	sfx->set_weight( core::scoring::res_type_constraint, 1.0);
	sfx->set_weight( core::scoring::res_type_linking_constraint, 1.0);
	protocols::relax::FastRelax relaxMover(sfx);


	scoring::constraints::ConstraintSetOP cst_set_post_built;
	if ( op_user_remodel_repeat_structure_ ) {
		// at this stage it should hold generic cstfile and res_type_linking constraints
		cst_set_post_built = new scoring::constraints::ConstraintSet( *pose.constraint_set() );
	}

	protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;
	if ( op_user_remodel_repeat_structure_ ) {

		// Dihedral (NCS) Constraints, need to be updated each mutation cycle for sidechain symmetry
		Size repeat_number = op_remodel_repeat_structure_;
		Size segment_length = (pose.n_residue())/repeat_number;

		for ( Size rep = 1; rep < repeat_number-1; rep++ ) { // from 1 since first segment don't need self-linking
			std::stringstream templateRangeSS;
			templateRangeSS << "2-" << segment_length+1; // offset by one to work around the termini
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep)+1 << "-" << segment_length + (segment_length*rep)+1;
			TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}

		for (Size rep = 1; rep < repeat_number-1; rep++){ // from 1 since first segment don't need self-linking
			std::stringstream templateRangeSS;
			templateRangeSS << "3-" << segment_length+2; // offset by one to work around the termini
			std::stringstream targetSS;
			targetSS << 1+(segment_length*rep)+2 << "-" << segment_length + (segment_length*rep)+2;
			TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
			setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
		}

		std::stringstream templateRangeSS;
		// take care of the terminal repeat, since the numbers are offset.
		templateRangeSS << "2-" << segment_length-1; // offset by one to work around the termini
		std::stringstream targetSS;
		targetSS << 1+(segment_length*(repeat_number-1))+1 << "-" << segment_length + (segment_length*(repeat_number-1))-1;
		TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
		setup_ncs.add_group(templateRangeSS.str(), targetSS.str());

	}

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {

		TR << "design_refine_seq_relax(): dr_cycle: " << i << std::endl;

		designMover.set_state("finish");
		designMover.apply(pose);

		//update dihedral constraint for repeat structures
		if ( op_user_remodel_repeat_structure_ ) {
			setup_ncs.apply(pose);

			//total hack (for now), see if the restypeset fails when initialized twice.
			if (op_user_RemodelLoopMover_cyclic_peptide_){
				protocols::forge::methods::cyclize_pose(pose);
			}

			sfx->show(TR, pose);
			TR << std::endl;
		}

		TR << "design_refine_seq_relax(): calling RelaxMover apply()." << std::endl;
		relaxMover.apply(pose);

		// reset constraints without NCS
		pose.constraint_set(cst_set_post_built);
		TR << "\n";
		sfx->show(TR, pose);
		TR << std::endl;
	}


	// turning off weights
	sfx->set_weight(core::scoring::coordinate_constraint, 0.0 );
	sfx->set_weight(core::scoring::atom_pair_constraint, 0.0 );
	sfx->set_weight(core::scoring::angle_constraint, 0.0 );
	sfx->set_weight(core::scoring::dihedral_constraint, 0.0 );
	sfx->set_weight(core::scoring::res_type_constraint, 0.0);
	sfx->set_weight(core::scoring::res_type_linking_constraint, 0.0);

	(*sfx)( pose );

	return true;
}

bool RemodelMover::design_refine_cart_relax(
	Pose & pose,
	RemodelDesignMover & designMover
)
{
	using core::kinematics::FoldTree;
	using core::pack::task::operation::RestrictResidueToRepacking;
	using core::pack::task::operation::RestrictResidueToRepackingOP;
	using core::pack::task::operation::RestrictToRepacking;
	using core::scoring::STANDARD_WTS;
	using core::scoring::SCORE12_PATCH;
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::forge::build::SegmentInsert;
	using namespace protocols::loops;
	using protocols::loops::Loops;
	using protocols::loops::loop_mover::refine::LoopMover_Refine_CCD;
	using protocols::simple_moves::PackRotamersMover;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;
	using namespace core::scoring::constraints;
	using namespace basic::options;


	using core::pose::annotated_to_oneletter_sequence;
	using protocols::forge::methods::intervals_to_loops;
	using protocols::forge::methods::linear_chainbreak;
	using protocols::loops::remove_cutpoint_variants;

	typedef protocols::forge::build::BuildManager::Positions Positions;

	// collect new regions/positions
	std::set< Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	Loops loops = intervals_to_loops( loop_intervals.begin(), loop_intervals.end() );

	if (op_user_remodel_repeat_structure_ || op_remodel_free_relax_ ){
		//do nothing
	} else {
		protocols::forge::methods::fill_non_loop_cst_set(pose, loops);
	}

	// safety, clear the energies object
	pose.energies().clear();

	//set simple tree
	FoldTree minFT;
	minFT.simple_tree(pose.total_residue());
	pose.fold_tree(minFT);

	// for refinement, always use standard repulsive
	ScoreFunctionOP sfx = core::scoring::getScoreFunction();
	//turning on weights
	sfx->set_weight(core::scoring::coordinate_constraint, 1.0 );
	sfx->set_weight(core::scoring::atom_pair_constraint, 1.0 );
	sfx->set_weight(core::scoring::angle_constraint, 1.0 );
	sfx->set_weight(core::scoring::dihedral_constraint, 10.0 ); // 1.0 originally
	sfx->set_weight(core::scoring::res_type_constraint, 1.0);
	sfx->set_weight(core::scoring::res_type_linking_constraint, 1.0);
	sfx->set_weight(core::scoring::cart_bonded, 0.5);


 	core::kinematics::MoveMapOP cmmop = new core::kinematics::MoveMap;
	//pose.dump_pdb("pretest.pdb");

	if (op_user_remodel_free_relax_) {
		for (Size i = 1; i<= pose.total_residue(); ++i){
			cmmop->set_bb(i, true);
			cmmop->set_chi(i, true);
		}
	}
	else {
		cmmop->import(remodel_data_.natro_movemap_);
		cmmop->import( manager_.movemap() );
	}

	for (Size i = 1; i<= pose.total_residue(); ++i){
		std::cout << "bb at " << i << " " << cmmop->get_bb(i) << std::endl;
		std::cout << "chi at " << i << " " << cmmop->get_chi(i) << std::endl;
		cmmop->set_chi(i,true);
		std::cout << "chi at " << i << " " << cmmop->get_chi(i) << std::endl;
	}

	for (Size i = 1; i<= pose.total_residue(); ++i){
		std::cout << "bbM at " << i << " " << manager_.movemap().get_bb(i) << std::endl;
		std::cout << "chiM at " << i << " " << manager_.movemap().get_chi(i) << std::endl;
	}

	simple_moves::MinMoverOP minMover = new simple_moves::MinMover( cmmop , sfx , "lbfgs_armijo", 0.01, true);
	minMover->cartesian(true);


	ConstraintSetOP cst_set_post_built;
	if (op_user_remodel_repeat_structure_){

	// at this stage it should hold generic cstfile and res_type_linking
	// constraints
		cst_set_post_built = new ConstraintSet( *pose.constraint_set() );
	}

	protocols::simple_moves::symmetry::SetupNCSMover setup_ncs;
	if (op_user_remodel_repeat_structure_){
				//Dihedral (NCS) Constraints, need to be updated each mutation cycle for sidechain symmetry

				Size repeat_number = op_remodel_repeat_structure_;
				Size segment_length = (pose.n_residue())/repeat_number;


				for (Size rep = 1; rep < repeat_number-1; rep++){ // from 1 since first segment don't need self-linking
					std::stringstream templateRangeSS;
					templateRangeSS << "2-" << segment_length+1; // offset by one to work around the termini
					std::stringstream targetSS;
					targetSS << 1+(segment_length*rep)+1 << "-" << segment_length + (segment_length*rep)+1;
				  TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
					setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
				}

			for (Size rep = 1; rep < repeat_number-1; rep++){ // from 1 since first segment don't need self-linking
					std::stringstream templateRangeSS;
					templateRangeSS << "3-" << segment_length+2; // offset by one to work around the termini
					std::stringstream targetSS;
					targetSS << 1+(segment_length*rep)+2 << "-" << segment_length + (segment_length*rep)+2;
				  TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
					setup_ncs.add_group(templateRangeSS.str(), targetSS.str());
				}


				std::stringstream templateRangeSS;
				//take care of the terminal repeat, since the numbers are offset.
				templateRangeSS << "2-" << segment_length-1; // offset by one to work around the termini
				std::stringstream targetSS;
				targetSS << 1+(segment_length*(repeat_number-1))+1 << "-" << segment_length + (segment_length*(repeat_number-1))-1;
				TR << "NCS " << templateRangeSS.str() << " " << targetSS.str() << std::endl;
				setup_ncs.add_group(templateRangeSS.str(), targetSS.str());

	}

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {


		designMover.set_state("finish");
		designMover.apply(pose);

		//update dihedral constraint for repeat structures
		if (op_user_remodel_repeat_structure_){
			setup_ncs.apply(pose);

			//total hack (for now), see if the restypeset fails when initialized
			//twice.
			if (op_user_RemodelLoopMover_cyclic_peptide_){
		     protocols::forge::methods::cyclize_pose(pose);
		  }

	//		sfx->show(TR, pose);
	//		TR << std::endl;
		}

		minMover->apply(pose);

		//reset constraints without NCS
		pose.constraint_set(cst_set_post_built);
			sfx->show(TR, pose);
			TR << std::endl;
	}


//turning off weights
  sfx->set_weight(core::scoring::coordinate_constraint, 0.0 );
  sfx->set_weight(core::scoring::atom_pair_constraint, 0.0 );
  sfx->set_weight(core::scoring::angle_constraint, 0.0 );
  sfx->set_weight(core::scoring::dihedral_constraint, 0.0 );
  sfx->set_weight(core::scoring::res_type_constraint, 0.0);
  sfx->set_weight(core::scoring::res_type_linking_constraint, 0.0);

	(*sfx)( pose );

	return true;
}


///
/// @begin RemodelMover::design_refine
///
/// @brief
/// Run the design-refine stage.
/// Checks the value of -remodel:repeat_structure and -remodel:swap_refine_confirm_protocols
/// NOTE: CURRENTLY ALWAYS RETURNS TRUE regardless of if chain breaks test passes or fails
///
bool RemodelMover::design_refine( Pose & pose, RemodelDesignMover & designMover ) {

	using namespace core;
	using namespace protocols;
	using namespace protocols::toolbox::task_operations;

	//using core::pack::task::operation::RestrictResidueToRepacking;
	//using core::pack::task::operation::RestrictResidueToRepackingOP;
	//using core::pack::task::operation::RestrictToRepacking;
	//using protocols::forge::build::SegmentInsert;
	//using protocols::loops::Loops;
	//using protocols::loops::LoopMover_Refine_CCD;
	//using protocols::simple_moves::PackRotamersMover;
	//using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;

	//using core::pose::annotated_to_oneletter_sequence;
	//using protocols::forge::methods::intervals_to_loops;
	//using protocols::forge::methods::linear_chainbreak;
	//using protocols::loops::remove_cutpoint_variants;

	//typedef protocols::forge::build::BuildManager::Positions Positions;

	// collect new regions/positions
	std::set< forge::build::Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	forge::build::BuildManager::Original2Modified original2modified_interval_endpoints = manager_.original2modified_interval_endpoints();

	// collect loops
	loops::LoopsOP loops = new loops::Loops( forge::methods::intervals_to_loops( loop_intervals.begin(), loop_intervals.end() ) );

	// refine Mover used doesn't setup a fold tree, so do it here
	//FoldTree loop_ft = protocols::forge::methods::fold_tree_from_loops( pose, loops );
	kinematics::FoldTree loop_ft;
	loops::fold_tree_from_loops( pose, *loops, loop_ft, true /*term cut*/);

	// save original fold tree
	kinematics::FoldTree original_ft = pose.fold_tree();

	// define the score function
	//ScoreFunctionOP sfx = fullatom_sfx_->clone();
	//for refinement always use hard repulsive
	scoring::ScoreFunctionOP sfx = core::scoring::ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );

//turning on weights, for paranoya
  sfx->set_weight(core::scoring::atom_pair_constraint, 1.0 );
  sfx->set_weight(core::scoring::dihedral_constraint, 10.0 ); // 1.0 originally

	// setup the refine TaskFactory
	pack::task::TaskFactoryOP refine_tf = generic_taskfactory();
	refine_tf->push_back( new toolbox::task_operations::RestrictToNeighborhoodOperation( neighborhood_calc_name() ) );
	refine_tf->push_back( new pack::task::operation::RestrictToRepacking() );

	// safety, clear the energies object
	pose.energies().clear();

	// run design-refine cycle
	for ( Size i = 0; i < dr_cycles_; ++i ) {

		// design the new section
		//PackRotamersMover design( sfx );
		//design.task_factory( design_tf );
		//design.apply( pose );
		designMover.set_state("finish");
		designMover.apply( pose );

		// set loop topology
		pose.fold_tree( loop_ft );

		if ( !op_user_remodel_swap_refine_confirm_protocols_ ) {
			// refine the new section
			loops::loop_mover::refine::LoopMover_Refine_CCD refine( loops, sfx );
			kinematics::MoveMapOP combined_mm = new kinematics::MoveMap();

			////// fix dna
			for ( Size i=1; i<=pose.total_residue() ; ++i ) {
				if ( pose.residue(i).is_DNA() ) {
					TR << "NATRO movemap setup: turning off DNA bb and chi move for refinement stage" << std::endl;
					remodel_data_.natro_movemap_.set_bb( i, false );
					remodel_data_.natro_movemap_.set_chi( i, false );
				}
			}
			////// end fix dna

			combined_mm->import( remodel_data_.natro_movemap_ );
			combined_mm->import( manager_.movemap() );

				//remodel_data_.natro_movemap_.show(pose.total_residue());
				//manager_.movemap().show(pose.total_residue());
			//combined_mm->show(pose.total_residue());
			//modify task to accept NATRO definition
			utility::vector1<core::Size> natroPositions;
			for (Size i = 1; i<= pose.total_residue(); i++){
				if (remodel_data_.natro_movemap_.get_chi(i) == 0){
					natroPositions.push_back(i);
				}
			}

			pack::task::operation::OperateOnCertainResiduesOP natroRes = new pack::task::operation::OperateOnCertainResidues;
			natroRes->residue_indices( natroPositions );
			natroRes->op( new pack::task::operation::PreventRepackingRLT );
			refine_tf->push_back( natroRes );

			refine.false_movemap( combined_mm );
			refine.set_task_factory( refine_tf );
			refine.apply( pose );

		} else {
			loops::loop_mover::refine::LoopMover_Refine_KIC KIC( loops );
			KIC.apply(pose);
		}

		// remove cutpoint variants -- shouldn't this happen at the end
		// of the refine Mover?
		loops::remove_cutpoint_variants( pose );


#ifdef FILE_DEBUG
		std::stringstream SS;
		SS << "RefineStage" << i << ".pdb";
		pose.dump_pdb( SS.str() );
#endif
	}

//turning off weights, for paranoya
  sfx->set_weight(core::scoring::atom_pair_constraint, 0.0 );
  sfx->set_weight(core::scoring::dihedral_constraint, 0.0 ); // 1.0 originally

	// must score one last time since we've removed variants and set
	// new topology, otherwise component energies not correct for
	// e.g. structure output
	(*sfx)( pose );

	// evaluate all chainbreaks using linear chainbreak
	bool cbreaks_pass = true;
	for ( loops::Loops::const_iterator l = loops->begin(), le = loops->end(); l != le && cbreaks_pass; ++l ) {
		if ( l->cut() > 0 ) {
			Real const c = forge::methods::linear_chainbreak( pose, l->cut() );
			TR << "design_refine: final chainbreak = " << c  << " at " << l->cut() << std::endl;
			cbreaks_pass = c <= max_linear_chainbreak_;
		}
	}
		// set original topology
		pose.fold_tree( original_ft );

	return cbreaks_pass;
	//return true; //FOR NOW!!  change me back!
}


///
/// @begin RemodelMover::confirm_sequence
///
/// @brief
/// As best as I can tell, does some loop closure and calculates RMSD to native. Returns true.
/// NOTE: CURRENTLY ALWAYS RETURNS TRUE regardless of rmsd value, because this stage is not being used as a filter
/// Checks the value of -remodel::swap_refine_confirm_protocols
///
bool RemodelMover::confirm_sequence( core::pose::Pose & pose ) {

	using namespace protocols::forge::methods;
 	using protocols::forge::build::Interval;
	using protocols::toolbox::task_operations::RestrictToNeighborhoodOperation;
	using pack::task::operation::RestrictToRepacking;

	pose::Pose archive_pose = pose;  //for rmsd

	std::set< forge::build::Interval > loop_intervals = manager_.intervals_containing_undefined_positions();
	//pose.dump_pdb("pre_KICpose.pdb");

	// collect loops
	loops::LoopsOP confirmation_loops = new loops::Loops( intervals_to_confirmation_loops( loop_intervals.begin(), loop_intervals.end(), pose.total_residue() ) );

	// refine Mover used doesn't setup a fold tree, so do it here
	kinematics::FoldTree loop_ft;
	loops::fold_tree_from_loops( pose, *confirmation_loops, loop_ft, true );
	TR << "confirmation loops tree" << loop_ft << std::endl;

	// save original fold tree
	kinematics::FoldTree original_ft = pose.fold_tree();

	// switch over to new tree
	pose.fold_tree(loop_ft);

	//LoopMover_Refine_KIC KIC(confirmation_loops);

	TR << "fold tree entering confirmation: " << pose.fold_tree() << std::endl;

	//KIC.apply(pose);

	if ( op_user_remodel_swap_refine_confirm_protocols_ ) {
		TR << "REFINE USING CCD" << std::endl;
		// refine the new section
		// setup the refine TaskFactory
		//
		//	protocols::forge::remodel::RemodelLoopMover scramble_mover(confirmation_loops);
		//	scramble_mover.randomize_stage(pose);

		TaskFactoryOP refine_tf = generic_taskfactory();
		refine_tf->push_back( new RestrictToNeighborhoodOperation( neighborhood_calc_name() ) );
		refine_tf->push_back( new RestrictToRepacking() );

		loops::loop_mover::refine::LoopMover_Refine_CCD refine( confirmation_loops, fullatom_sfx_ );
		kinematics::MoveMapOP combined_mm = new kinematics::MoveMap();

		////// fix dna
		for ( Size i=1; i<=pose.total_residue() ; ++i ) {
			if ( pose.residue( i ).is_DNA() ) {
				TR << "NATRO movemap setup: turning off DNA bb and chi move for refinement stage" << std::endl;
				remodel_data_.natro_movemap_.set_bb( i, false );
				remodel_data_.natro_movemap_.set_chi( i, false );
			}
		}
		////// end fix dna

		combined_mm->import(remodel_data_.natro_movemap_);
		combined_mm->import( manager_.movemap() );

		refine.false_movemap( combined_mm );
		refine.set_task_factory( refine_tf );
		refine.apply( pose );

	} else {
		TR << "REFINE USING KIC" << std::endl;
		loops::loop_mover::refine::LoopMover_Refine_KIC KIC( confirmation_loops );
		KIC.apply(pose);
	}

	// reset to original foldtree
	pose.fold_tree( original_ft );

	//pose.dump_pdb("post_KICpose.pdb");

	// rmsd_calculation:

	Real sum_sd = 0;
	Real sum_sd_native = 0;
	Real sum_sd_archive2native=0;
	Size atom_count = 0;

	for ( loops::Loops::iterator it = confirmation_loops->v_begin(), end = confirmation_loops->v_end(); it!=end; it++) {
		for ( Size i = it->start(); i <= it->stop(); ++i ) {
			Real dist_squared = ( pose.residue(i).xyz( "CA" ) - archive_pose.residue(i).xyz( "CA" ) ).length_squared();
			Real dist_squared_native = ( pose.residue(i).xyz( "CA" ) - native_pose_.residue(i).xyz( "CA" ) ).length_squared();
			Real dist_squared_archive2native = ( archive_pose.residue(i).xyz( "CA" ) - native_pose_.residue(i).xyz( "CA" ) ).length_squared();
			sum_sd = sum_sd + dist_squared;
			sum_sd_native = sum_sd_native + dist_squared_native;
			sum_sd_archive2native = sum_sd_archive2native + dist_squared_archive2native;
			atom_count++;
#ifdef FILE_DEBUG
				std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared= " << dist_squared << std::endl;
				std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared_native= " << dist_squared_native << std::endl;
				std::cout << " Backbone atom= " <<  "CA res: " << i  << " dist_squared_archive2native= " << dist_squared_archive2native << std::endl;
#endif
		}
	}

	sum_sd = sum_sd / atom_count;
	sum_sd_native = sum_sd_native / atom_count;
	sum_sd_archive2native = sum_sd_archive2native / atom_count;

	Real rmsd = sqrt(sum_sd);
	Real rmsd_native = sqrt(sum_sd_native);
	Real rmsd_archive2native = sqrt(sum_sd_archive2native);

	pose::PDBInfoOP temp_pdbinfo = pose.pdb_info();

	pose::RemarkInfo remark;
	remark.value = "KIC confirmation RMSD: " + utility::to_string( rmsd ) + " to native RMSD: " + utility::to_string( rmsd_native );
	temp_pdbinfo->remarks().push_back( remark );

	remark.value = " ARCHIVE2NATIVE RMSD: " + utility::to_string(rmsd_archive2native);

	temp_pdbinfo->remarks().push_back( remark );

	pose.pdb_info(temp_pdbinfo);

	TR << "RMSD of KIC conformation: " << rmsd << std::endl;
	TR << "RMSD of KIC conformation to native: " << rmsd_native << std::endl;
	TR << "RMSD of ARCHIVE to NATIVE: " << rmsd_archive2native << std::endl;

	// currently the confirmation stage is not setup as filter so always return true
	if ( rmsd <= 1 ) {
		return true;
	} else {
		return true; //for now CHANGE IT BACK!!
	}

}


///
/// @begin RemodelMover::generic_taskfactory
///
/// @brief
/// Returns a TaskFactory useable as a starting point for either design or refinement.
/// Only adds the NoRepackDisulfides, IncludeCurrent, and Init from command line ops. ReadResfile is not included.
///
RemodelMover::TaskFactoryOP RemodelMover::generic_taskfactory() {
	using protocols::toolbox::task_operations::LimitAromaChi2Operation;

	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	TaskFactoryOP tf = new TaskFactory();

	tf->push_back( new InitializeFromCommandline() ); // also inits -ex options
	tf->push_back( new IncludeCurrent() ); // enforce keeping of input sidechains
	tf->push_back( new NoRepackDisulfides() );
  if (!op_user_design_allow_rare_aro_chi_){
		tf->push_back( new LimitAromaChi2Operation() );
	}

	// load resfile op only if requested
	/*if ( !resfile_.empty() ) {
		ReadResfileOP rrf = new ReadResfile();
		rrf->filename( resfile_ );
		tf->push_back( rrf );
	}
	*/
	return tf;
}


///
/// @begin RemodelMover::process_continuous_design_string
///
/// @brief
/// process a continuous design string, adding appropriate operations to the TaskFactory
///
void RemodelMover::process_continuous_design_string( Interval const & original_interval, String const & design_str,
	Original2Modified const & original2modified_interval_endpoints, TaskFactoryOP design_tf ) {

	using namespace core;

	Size const offset = original2modified_interval_endpoints.find( original_interval.left )->second;
	for ( Size i = 0, ie = design_str.length(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		switch ( design_str.at( i ) ) {
			case 's': // surface case, no CFWY
				allowed_aa_types = allowed_surface_aa();
				break;
			case '.': // protocol default design
				continue;
			default: // regular case, single aa type
				allowed_aa_types[ chemical::aa_from_oneletter_code( design_str.at( i ) ) ] = true;
				break;
		}

		design_tf->push_back( new pack::task::operation::RestrictAbsentCanonicalAAS( i + offset, allowed_aa_types ) );
	}
}


///
/// @begin RemodelMover::process_insert_design_string
///
/// @brief
/// process a design string containing an insert, adding appropriate operations to the TaskFactory
///
void RemodelMover::process_insert_design_string( Interval const & original_interval, String const & design_str,
	Original2Modified const & original2modified_interval_endpoints, TaskFactoryOP design_tf ) {

	using namespace core;

	char const insert_char = forge::build::SegmentInsert::insertion_char();

	// Figure out the number of residues in each section.
	forge::build::Interval const interval(
		original2modified_interval_endpoints.find( original_interval.left )->second,
		original2modified_interval_endpoints.find( original_interval.right )->second
	);

	Size const insert_char_idx = design_str.find( insert_char );
	Size const left_nres = insert_char_idx;
	Size const right_nres = design_str.size() - left_nres - 1;
	Size const insert_nres = interval.length() - left_nres - right_nres;

	// Make setup easy by building a new design string to expand the
	// insertion character into a series of the insertion character
	// the size of the insert.
	String aa = design_str;
	aa.replace( insert_char_idx, 1, insert_nres, insert_char );

	// setup TaskOperations
	pack::task::operation::RestrictResidueToRepackingOP repack_op = new pack::task::operation::RestrictResidueToRepacking();

	Size const left_offset = interval.left;
	for ( Size i = 0, ie = aa.size(); i < ie; ++i ) {
		utility::vector1< bool > allowed_aa_types( 20, false );

		if ( aa.at( i ) == insert_char ) { // repack only
			repack_op->include_residue( i + left_offset );
			continue;

		} else if ( aa.at( i ) == 's' ) { // surface case, no CFWY
			allowed_aa_types = allowed_surface_aa();

		} else if ( aa.at( i ) == '.' ) { // protocol default design
			continue;

		} else { // regular case, single aa type
			allowed_aa_types[ chemical::aa_from_oneletter_code( aa.at( i ) ) ] = true;
		}

		design_tf->push_back( new pack::task::operation::RestrictAbsentCanonicalAAS( i + left_offset, allowed_aa_types ) );
	}

	design_tf->push_back( repack_op );
}


///
/// @begin RemodelMover::process_insert_design_string
///
/// @brief
/// return a boolean vector specifying allowed a.a. when designing on the surface
///
utility::vector1< bool > const & RemodelMover::allowed_surface_aa() {
	using core::chemical::aa_from_oneletter_code;

	static String surface_aa = "ADEGHIKLMNPQRSTV";
	static utility::vector1< bool > v( 20, false );

	for ( Size i = 0, ie = surface_aa.length(); i < ie; ++i ) {
		v[ aa_from_oneletter_code( surface_aa.at( i ) ) ] = true;
	}

	return v;
}

/// @brief parse xml
void
RemodelMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	DataMap & /*data*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/ )
{
	if( tag->hasOption("blueprint") ) {
		blueprint_ = tag->getOption<std::string>( "blueprint" );
	}
}



} // namespace remodel
} // namespace forge
} // namespace protocols
