// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SymDockingHiRes
/// @brief protocols that are specific to high resolution docking
/// @detailed
///		This contains the functions that create initial positions for docking
/// 	Also contains docking mcm protocol
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Sid Chaudhury
/// @author Modified by Jacob Corn
/// @author ingemar André. Based on the standard docking protocol

#include <protocols/symmetric_docking/SymDockingHiRes.hh>
#include <protocols/symmetric_docking/SymSidechainMinMover.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>

// AUTO-REMOVED #include <core/conformation/Interface.hh>

#include <basic/options/option.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
// AUTO-REMOVED #include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>
#include <core/conformation/Residue.hh> // for design() flag
#include <core/pack/task/operation/NoRepackDisulfides.hh>
// AUTO-REMOVED #include <core/pack/task/operation/OperateOnCertainResidues.hh>
// AUTO-REMOVED #include <core/pack/task/operation/ResLvlTaskOperations.hh> // PreventRepackingRLT
// AUTO-REMOVED #include <core/pack/task/operation/ResFilters.hh> // ResidueLacksProperty
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymDof.hh>

#include <protocols/moves/symmetry/SymMinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/symmetry/SymPackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/moves/OutputMovers.hh>
#include <protocols/moves/symmetry/SymRotamerTrialsMover.hh>
// AUTO-REMOVED #include <protocols/moves/RotamerTrialsMinMover.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/JumpOutMover.hh>
// AUTO-REMOVED #include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/moves/RepeatMover.hh>
//for resfile reading
#include <basic/options/keys/packing.OptionKeys.gen.hh>
// Auto-header: duplicate removed #include <basic/options/option.hh>

// AUTO-REMOVED #include <protocols/loops/loops_main.hh>
// AUTO-REMOVED #include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover.fwd.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_KIC.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_Backrub.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_CCD.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

// AUTO-REMOVED #include <utility/tag/Tag.hh> // REQUIRED FOR WINDOWS

#include <basic/Tracer.hh>
using basic::T;

// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.symmetric_docking.SymDockingHiRes");

using namespace core;

namespace protocols {
namespace symmetric_docking {

// default constructor
SymDockingHiRes::SymDockingHiRes() : Mover()
{
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" ) ;
	scorefxn_pack_ = core::scoring::ScoreFunctionFactory::create_score_function( "standard" ) ;
	moves::Mover::type( "SymDockingHiRes" );
	init_task_factory_=NULL;
	design_ = false;
}

// constructor with arguments
SymDockingHiRes::SymDockingHiRes(
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::scoring::ScoreFunctionOP scorefxn_pack_in
) : Mover(), scorefxn_(scorefxn_in), scorefxn_pack_(scorefxn_pack_in)
{
	moves::Mover::type( "SymDockingHRes" );
	init_task_factory_=NULL;
	design_ = false;
}

//destructor
SymDockingHiRes::~SymDockingHiRes() {}

//clone
protocols::moves::MoverOP SymDockingHiRes::clone() const {
	return new SymDockingHiRes(*this);
}
// what type of minimization is used?
void SymDockingHiRes::set_min_type( std::string min_type_in ) { min_type_ = min_type_in;}

// repack sidechains during docking?
void SymDockingHiRes::set_repack( bool repack_switch){ repack_switch_ = repack_switch;}

// pointer to mc object
moves::MonteCarloOP SymDockingHiRes::get_mc() { return mc_; }

// set the packer task
void
SymDockingHiRes::task_factory( core::pack::task::TaskFactoryOP task )
{
	init_task_factory_ = task;
}

core::pack::task::TaskFactoryOP &
SymDockingHiRes::task_factory(){
	return init_task_factory_;
}

// Use design during docking. Currently not tested...
void SymDockingHiRes::design( bool const des ) {
	design_ = des;
}

// Are we designing during the docking procedure
bool SymDockingHiRes::design() const { return design_; }

// Set default values for the docking protocol
void SymDockingHiRes::set_default( core::pose::Pose & pose ) {
		using namespace basic::options;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		// SETS UP the stuff in pose
		(*scorefxn_)( pose );

		// Get the docking translation and rotation magnitutes from command line. If
		// none is given the default values of 0.1 A and 5 degrees are used
		trans_magnitude_ = option[ OptionKeys::docking::dock_mcm_trans_magnitude ]();
		rot_magnitude_ = option[ OptionKeys::docking::dock_mcm_rot_magnitude ]();
		// use rtmin_ and sc_min. rtmin for symmetry is not implemented?! Turned off 		// for now. The flag is still here which adds to the confusion.

		rtmin_ = option[ OptionKeys::docking::dock_rtmin ]();
		scmin_ = option[ OptionKeys::docking::sc_min ]();

		temperature_ = 0.8;
		repack_switch_ = true;
		repack_period_ = 8;

		//sets up MC object
		mc_ = new moves::MonteCarlo( pose, *scorefxn_, temperature_ );

		//sets up default movemap
		bb_ = false;
		chi_ = false;
		movemap_ = new kinematics::MoveMap();
		movemap_->set_chi( chi_ );
		movemap_->set_bb( bb_ );
		core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );

		//sets up minimization parameters
		min_tolerance_ = 0.01;
		min_type_ = std::string( "dfpmin_armijo_nonmonotone" );
		nb_list_ = true;

		setup_packing( pose );

		//set up docking protocol based on options
		set_protocol( pose );
	}

// Set the movemap
void SymDockingHiRes::set_move_map(core::kinematics::MoveMapOP movemap_in){ movemap_ = movemap_in; }

// define which protocol is used...
void SymDockingHiRes::set_protocol( pose::Pose & pose ){
	using namespace basic::options;

	if ( option[ OptionKeys::docking::dock_min ]() ) {
		set_dock_min_protocol();
		}
		else if ( option[ OptionKeys::docking::dock_ppk ]() ){
			set_dock_ppk_protocol( pose );
		}

	else {
	  set_dock_mcm_protocol( pose );
		}
}

/*void SymDockingHiRes::define_loops( pose::Pose const & pose, loops::Loops & loop_set, Real & interface_dist ) {
	//runtime_assert( movable_jumps_.size() == 1 ); // CURRENTLY ONLY SUPPORTED WITH SIMPLE DOCKING
	//core::Size const rb_jump = movable_jumps_[1];

	loop_set.clear();

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	//RestrictTaskForDockingOP rtfd = new RestrictTaskForDocking( scorefxn_, rb_jump_, true, interface_dist );
	pack::task::TaskFactory tf;
	//tf.push_back( rtfd );
	//tf.push_back( new RestrictTaskForDocking( scorefxn_, rb_jump_, true, interface_dist ) );
	tf.push_back( new RestrictToInterface( movable_jumps_, interface_dist ) );
	pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( pose );

	// extend one residue beyond borders of repackable regions, don't allow 1-residue loops
	core::Size const nres = pose.total_residue();
	utility::vector1<bool> flexible_region( nres, false );
	for ( Size i=2; i < nres; ++i ) {
		int num_flexible(0);
		if ( pose.pdb_info()->chain(i-1) == pose.pdb_info()->chain(i+1) ) {
			if ( task->pack_residue(i-1) ) ++num_flexible;
			if ( task->pack_residue(i) ) ++num_flexible;
			if ( task->pack_residue(i+1) ) ++num_flexible;
		}
		if ( num_flexible > 1 ) {
			flexible_region.at(i-1) = true;
			flexible_region.at(i) = true;
			flexible_region.at(i+1) = true;
		}
	}

	// if we have a single fixed residue between two loops, make this flexible
	for ( Size i=2; i < nres; ++i ) {
		if ( flexible_region.at(i-1) && flexible_region.at(i+1) ) flexible_region.at(i) = true;
	}

	// jk For now, don't let the first or last two residues of a chain be flexible
	flexible_region.at(1) = false;
	flexible_region.at(2) = false;
	for ( Size i=3; i < (nres-1); ++i ) {
		if ( pose.pdb_info()->chain(i-1) != pose.pdb_info()->chain(i) ) {
			flexible_region.at(i-2) = false;
			flexible_region.at(i-1) = false;
			flexible_region.at(i) = false;
			flexible_region.at(i+1) = false;
		}
	}
	flexible_region.at(nres-1) = false;
	flexible_region.at(nres) = false;

	// disallow one-residue loops
	for ( Size i=2; i < nres; ++i ) {
		if ( ( ! flexible_region.at(i-1)) && ( ! flexible_region.at(i+1)) ) flexible_region.at(i) = false;
	}
	// disallow two-residue loops
	for ( Size i=3; i < nres; ++i ) {
		if ( ( ! flexible_region.at(i-2)) && ( ! flexible_region.at(i+1)) ) {
			flexible_region.at(i-1) = false;
			flexible_region.at(i) = false;
		}
	}

	// setup loops
	core::Size loop_start=0;
	core::Size loop_stop=0;
	for ( Size i=1; i < nres; ++i ) {
		if ( flexible_region.at(i) ) {
			loop_start = i;
			loop_stop = i;
			for ( Size j=i+1; j <= nres; ++j ) {
				// if j is on a different chain than i, break
				if ( pose.pdb_info()->chain(j) != pose.pdb_info()->chain(i) ) {
					break;
				}
				// if j is not flexible, break
				if ( ( ! flexible_region.at(j) ) ) {
					break;
				}
				loop_stop = j;
			}
			loop_set.add_loop( loop_start, loop_stop, 0 );
			i = loop_stop;
		}
	}

	loop_set.choose_cutpoints( pose );

	return;
} */


////////////////////////////////////////////////////////////////////////////////
/// @begin docking high resolution apply function
/// @brief
/// @detailed
///		decides what to call according to options
void SymDockingHiRes::apply( core::pose::Pose & pose )
{
	using namespace scoring;
	using namespace basic::options;

	TR << "in SymDockingHiRes.apply" << std::endl;

	// jec sanity check to avoid overwriting newly-set minimizers on every apply
	if ( !mc_ ) {
		set_default( pose );
	}

	mc_->reset( pose );

	docking_highres_protocol_mover_->apply( pose );
	mc_->recover_low( pose );
}

std::string
SymDockingHiRes::get_name() const {
	return "SymDockingHiRes";
}

///////////////////////////////////////////////////////////////////////////////////
/// @begin minimize_trial
///
/// @brief main entrance for normal rigid-body minimization
/// @detailed
///		retrieve the structure in the low array and do the normal minimization
///		by calling using a min_mover to optimize the score accourding to the
///		scorefunction that has been set
///
/// @remarks
///
/// @references docking_minimize_trial from docking_minimize.cc
///				pose_docking_minimize_trial from pose_docking.cc
///
/// @authors Monica Berrondo June 14 2007, modified for symmetric docking by
///					 Ingemar Andre
///
/// @last_modified October 15 2007
/////////////////////////////////////////////////////////////////////////////////
void SymDockingHiRes::set_dock_min_protocol() {
	using namespace moves;

	TR << "::::::::::::::::::DOCK_MIN:::::::::::::::::::" << std::endl;

	moves::MinMoverOP min_mover = new moves::symmetry::SymMinMover( movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_ );
	TrialMoverOP minimize_trial = new TrialMover( min_mover, mc_ );
	docking_highres_protocol_mover_ = new SequenceMover;
	docking_highres_protocol_mover_->add_mover( minimize_trial );
}

///////////////////////////////////////////////////////////////////////////////////
/// @begin dock_mcm_protocol
///
/// @brief main entrance to do monte carlo minimization
/// @detailed
///			a total of 50 cycles of monte-carlo minimization will be
///			carried out if the minimized structure can pass the filter
///			after the first and fifth cycle.  Then it is rigid-body minimized
///			to a stringent tolerance.
///
/// @remarks
///
/// @references docking_mcm_protocol from docking_minimize.cc
///				pose_docking_monte_carlo_minimize from pose_docking.cc
///
/// @authors Sid Chaudhury May 28 2009, modified for symmetric docking
///					 by Ingemar Andre
///
/// @last_modified April 30 2008
/////////////////////////////////////////////////////////////////////////////////
void SymDockingHiRes::set_dock_mcm_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::conformation::symmetry;
	using namespace protocols::toolbox::task_operations;

	assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

  std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	//set up rigid body movers
	RigidBodyDofSeqPerturbMoverOP rb_perturb = new RigidBodyDofSeqPerturbMover( dofs , rot_magnitude_, trans_magnitude_ );

	//set up minimizer movers
	moves::MinMoverOP min_mover = new moves::symmetry::SymMinMover( movemap_, scorefxn_, min_type_, min_tolerance_, nb_list_ );

	//set up sidechain movers for each movable jump
	tf_->push_back( new RestrictToInterface( 1 ) );

	RotamerTrialsMoverOP pack_rottrial = new moves::symmetry::SymRotamerTrialsMover( scorefxn_pack_, tf_ );

	SequenceMoverOP interface_repack_and_move_loops = new moves::SequenceMover;

	std::string const flex_bb_docking_type = option[ OptionKeys::docking::flexible_bb_docking ]();
	if ( flex_bb_docking_type == "fixedbb" ) {
		// Call pack_rotamers, no backbone movement
		PackRotamersMoverOP pack_interface_repack = new moves::symmetry::SymPackRotamersMover( scorefxn_pack_ );
		pack_interface_repack->task_factory(tf_);
		interface_repack_and_move_loops->add_mover( pack_interface_repack );
	} else {

		// Call pack_rotamer before and after loop movement
		PackRotamersMoverOP pack_interface_repack = new moves::symmetry::SymPackRotamersMover( scorefxn_pack_ );
		pack_interface_repack->task_factory(tf_);
		interface_repack_and_move_loops->add_mover( pack_interface_repack );

		// This does not work for symmetry yet. Anyone is free to implement it...
		if ( ( flex_bb_docking_type == "ccd" ) || ( flex_bb_docking_type == "kic" ) ||  ( flex_bb_docking_type == "backrub" ) ) {

/*			core::kinematics::FoldTree docking_fold_tree( pose.fold_tree() );
			core::kinematics::FoldTree loop_fold_tree = docking_fold_tree;

			// jk Create a copy because we can't modify input pose, and we need energies and the docking fold tree
			// Note: we need the docking fold tree so that we identify the correct "interface" when picking loops (ie. jump #1)
			// jk For now, define loops based on the input structure
			// (so that they're always the same, and we don't have to worry about relaxing a priori, akin to prepacking)
			pose::Pose pose_for_loop_defn = *get_input_pose();
			(*scorefxn_pack_)(pose_for_loop_defn);
			pose_for_loop_defn.fold_tree( docking_fold_tree );

			protocols::loops::Loops loop_set;
			Real interface_dist = option[ OptionKeys::docking::flexible_bb_docking_interface_dist ];
			define_loops( pose_for_loop_defn, loop_set, interface_dist );

			loops::LoopMoverOP loop_refine;
			if ( flex_bb_docking_type == "ccd" ) {
				// jk CCD loop refinement (fullatom only)
				TR << "Setting up for ccd loop modeling" << std::endl;
				protocols::loops::fold_tree_from_loops( pose, loop_set, loop_fold_tree );
				// need to pass a clone of the scorefxn because LoopMover requires a non-const scorefxn
				loop_refine = new loops::LoopMover_Refine_CCD( loop_set, scorefxn_pack_->clone() );
			} else if ( flex_bb_docking_type == "kic" ) {
				// jk KIC loop refinement (fullatom only)
				TR << "Setting up for kinematic (kic) loop modeling" << std::endl;
				protocols::loops::fold_tree_from_loops( pose, loop_set, loop_fold_tree );
				// need to pass a clone of the scorefxn because LoopMover requires a non-const scorefxn
				loop_refine = new loops::LoopMover_Refine_KIC( loop_set, scorefxn_pack_->clone() );
			} else if ( flex_bb_docking_type == "backrub" ) {
				// jk backrub loop refinement (fullatom only)
				TR << "Setting up for backrub loop modeling" << std::endl;
				// jk backrub can't use segments that span a jump residue, so use a simple fold tree here
				// note: this assumes that the termini are not allowed to move (and in define_loops they aren't)
				loop_fold_tree.simple_tree( pose.total_residue() );
				// need to pass a clone of the scorefxn because LoopMover requires a non-const scorefxn
				loop_refine = new loops::LoopMover_Refine_Backrub( loop_set, scorefxn_pack_->clone() );
			}

			moves::ChangeFoldTreeMoverOP get_loop_ft = new moves::ChangeFoldTreeMover( loop_fold_tree );
			moves::ChangeFoldTreeMoverOP get_docking_ft = new moves::ChangeFoldTreeMover( docking_fold_tree );

			interface_repack_and_move_loops->add_mover( get_loop_ft );
			interface_repack_and_move_loops->add_mover( loop_refine );
			interface_repack_and_move_loops->add_mover( get_docking_ft );
			*/
			TR << "[ ERROR ] flexible_bb_docking is not implemented for symmetric docking yet..." << std::endl;
			exit(1);

		} else {
			TR << "[ ERROR ] Unknown flexible_bb_docking type: " << flex_bb_docking_type << std::endl;
			exit(1);
		}

		interface_repack_and_move_loops->add_mover( pack_interface_repack );
	}


	TrialMoverOP pack_interface_and_move_loops_trial = new TrialMover( interface_repack_and_move_loops, mc_ );

//	RotamerTrialsMinMoverOP rtmin = new RotamerTrialsMinMover( scorefxn_pack_, tf_ );
//	TrialMoverOP rtmin_trial = new TrialMover( rtmin, mc_ );

	//InterfaceSidechainMinMoverOP scmin_mover = new InterfaceSidechainMinMover(rb_jump_, scorefxn_pack_ );
	SymSidechainMinMoverOP scmin_mover = new SymSidechainMinMover(scorefxn_pack_, tf_ );
	TrialMoverOP scmin_trial = new TrialMover( scmin_mover, mc_ );

	// the standard mcm cycle : rb perturbation->rotamer trials->minimization->MC accept
	SequenceMoverOP rb_mover = new SequenceMover;
	rb_mover->add_mover( rb_perturb );
	if ( repack_switch_ ) rb_mover->add_mover( pack_rottrial );

	core::Real minimization_threshold = 15.0;
	JumpOutMoverOP rb_mover_min = new JumpOutMover( rb_mover, min_mover, scorefxn_, minimization_threshold );
	TrialMoverOP rb_mover_min_trial = new TrialMover( rb_mover_min, mc_ );

	//every step (rb_mover_min_trial): the standard mcm cycle
	//every 8th step (repack_step): standard mcm cycle + repacking
	//try moving loops too, if desired

	SequenceMoverOP repack_step = new SequenceMover;
	repack_step->add_mover(rb_mover_min_trial);

	if ( repack_switch_ ){
		repack_step->add_mover( pack_interface_and_move_loops_trial );
	//	if (rtmin_) repack_step->add_mover(rtmin_trial);
		if (scmin_) repack_step->add_mover(scmin_trial);
		}

	CycleMoverOP rb_mover_min_trial_repack  = new CycleMover;
	for ( Size i=1; i<repack_period_; ++i ) rb_mover_min_trial_repack->add_mover( rb_mover_min_trial );
	rb_mover_min_trial_repack->add_mover( repack_step );

	//set up initial repack mover
	SequenceMoverOP initial_repack = new SequenceMover;
	initial_repack->add_mover(pack_interface_and_move_loops_trial);
	//if (rtmin_) initial_repack->add_mover(rtmin_trial);
	if (scmin_) initial_repack->add_mover(scmin_trial);

	//set up initial and final min_trial movers for docking
	TrialMoverOP minimize_trial = new TrialMover( min_mover, mc_ );

	//set up mcm cycles and mcm_repack cycles
	RepeatMoverOP mcm_four_cycles = new RepeatMover( rb_mover_min_trial, 4 );
	RepeatMoverOP mcm_fortyfive_cycles = new RepeatMover( rb_mover_min_trial_repack, 45 );
	//set up protocol mover
	TR << "::::::::::::::::::DOCK_MCM:::::::::::::::::::" << std::endl;

	docking_highres_protocol_mover_ = new SequenceMover;
	if (repack_switch_) docking_highres_protocol_mover_->add_mover( initial_repack );
	docking_highres_protocol_mover_->add_mover( minimize_trial );
	docking_highres_protocol_mover_->add_mover( mcm_four_cycles );
	docking_highres_protocol_mover_->add_mover( mcm_fortyfive_cycles );
	docking_highres_protocol_mover_->add_mover( minimize_trial );

}

// noncost pose for load_unboundrot csts
// detup the options for packing including setup up the packer task
void SymDockingHiRes::setup_packing( core::pose::Pose & pose ) {
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	//set upconstructor packer options
	tf_ = new TaskFactory;
	if( init_task_factory_ )
	{
		TR << "Using user-defined TaskFactory." << std::endl;
		tf_ = new TaskFactory( *init_task_factory_ );
	}
	if( design_ ) {
		TR << "Designing during docking" << std::endl;
	}
	else { // default case -- restrict everything to repacking.
		tf_->push_back( new RestrictToRepacking );
	}
//	tf_->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
	tf_->push_back( new InitializeFromCommandline );
	tf_->push_back( new IncludeCurrent );
	tf_->push_back( new NoRepackDisulfides );
	if( option[OptionKeys::packing::resfile].user() ) tf_->push_back( new ReadResfile );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot = new core::pack::rotamer_set::UnboundRotamersOperation();
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation = new operation::AppendRotamerSet( unboundrot );
	tf_->push_back( unboundrot_operation );
	core::pack::dunbrack::load_unboundrot(pose); // adds scoring bonuses for the "unbound" rotamers, if any

	// note that RestrictToInterfaceOperation is added during set_dock_mcm_protocol
}

/// define the prepacking protocol
void SymDockingHiRes::set_dock_ppk_protocol( core::pose::Pose & pose ) {
	using namespace moves;
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::conformation::symmetry;

	assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

  std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	TR << "::::::::::::::::::DOCK_PPK:::::::::::::::::::" << std::endl;

	//set up translate-by-axis movers
	RigidBodyDofSeqTransMoverOP translate_away( new RigidBodyDofSeqTransMover( dofs ) );
	RigidBodyDofSeqTransMoverOP translate_back( new RigidBodyDofSeqTransMover( dofs ) );
	translate_away->step_size( 1000 );
	translate_back->step_size( -1000 );


	PackerTaskOP task = tf_->create_task_and_apply_taskoperations( pose ); // does not include restrict to interface

	PackRotamersMoverOP prepack_full_repack = new symmetry::SymPackRotamersMover( scorefxn_pack_, task );
	//RotamerTrialsMinMoverOP rtmin_mover = new symmetry::SymRotamerTrialsMinMover( scorefxn_pack_, *task );
	SymSidechainMinMoverOP scmin_mover = new SymSidechainMinMover(scorefxn_pack_, task);

	// set up protocol
	docking_highres_protocol_mover_ = new SequenceMover;
	if (scmin_) docking_highres_protocol_mover_->add_mover( scmin_mover );
	docking_highres_protocol_mover_->add_mover( translate_away );
	docking_highres_protocol_mover_->add_mover( prepack_full_repack );
	//if (rtmin_) docking_highres_protocol_mover_->add_mover( rtmin_mover );
	if (scmin_) docking_highres_protocol_mover_->add_mover( scmin_mover );
	docking_highres_protocol_mover_->add_mover( translate_back );
}



} // namespace docking
} // namespace protocols
