// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file protocols/loops/LoopMover_Backrub.cc
/// @brief backrub loop refine
/// @author J. Karanicolas

//// Unit Headers
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/LoopMover_Backrub.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/BackrubMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverStatus.hh>
//
//// Rosetta Headers
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
// AUTO-REMOVED #include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>

#include <basic/options/keys/run.OptionKeys.gen.hh>

// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>

// AUTO-REMOVED #include <basic/prof.hh> // profiling
#include <basic/Tracer.hh>

//Utility Headers
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>

// C++ Headers
#include <iostream>
#include <map>
#include <string>

// option key includes

// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/keys/Key3Vector.hh>


namespace protocols {
namespace loops {

///////////////////////////////////////////////////////////////////////////////
using namespace core;

static numeric::random::RandomGenerator RG(42678);

extern basic::Tracer tr;


LoopMover_Refine_Backrub::LoopMover_Refine_Backrub(
	protocols::loops::Loops  loops_in
) : LoopMover( loops_in )
{
	scorefxn_ = core::scoring::getScoreFunction();
	protocols::moves::Mover::type("LoopMover_Refine_Backrub");
	set_default_settings();
}

LoopMover_Refine_Backrub::LoopMover_Refine_Backrub(
	protocols::loops::Loops loops_in,
	core::scoring::ScoreFunctionOP  scorefxn
) : LoopMover( loops_in )
{
	scorefxn_ = scorefxn;
	protocols::moves::Mover::type("LoopMover_Refine_Backrub");
	set_default_settings();
}

//destructor
LoopMover_Refine_Backrub::~LoopMover_Refine_Backrub(){}

void LoopMover_Refine_Backrub::set_task_factory( core::pack::task::TaskFactoryOP value ){ task_factory = value;}
bool LoopMover_Refine_Backrub::get_task_factory(){return task_factory;}


void LoopMover_Refine_Backrub::apply(
	core::pose::Pose & pose
){

	using namespace core;
	using namespace scoring;
	using namespace basic::options;

	// Note: based heavily on KIC loop refine

	// scheduler
	int const fast = option[OptionKeys::loops::fast];
	int outer_cycles(3);
	if ( option[ OptionKeys::loops::outer_cycles ].user() ) {
		outer_cycles = option[ OptionKeys::loops::outer_cycles ]();
	}
	if ( option[ OptionKeys::run::test_cycles ]() ) {
		outer_cycles = 2;
	}
	int max_inner_cycles( 200 );
	if ( option[ OptionKeys::loops::max_inner_cycles ].user() ) {
		max_inner_cycles = option[ OptionKeys::loops::max_inner_cycles ]();
	}
	if ( option[ OptionKeys::run::test_cycles ]() ) {
		max_inner_cycles = 2;
	}

	int const inner_cycles = std::min( Size(max_inner_cycles), fast ? (int)loops_.loop_size() : 10*loops_.loop_size() );
	int repack_period = 20; // should be an option
	if ( option[ OptionKeys::loops::repack_period ].user() ) {
		repack_period = option[ OptionKeys::loops::repack_period ]();
	}

	int ntrials_per_cycle = 10;
	if ( option[ OptionKeys::loops::backrub_trials ].user() ) {
		ntrials_per_cycle = option[ OptionKeys::loops::backrub_trials ]();
	}

	// scorefxn
	scoring::ScoreFunctionOP scorefxn;
	if ( scorefxn_ != 0 ) scorefxn = scorefxn_->clone();
	else {
		scorefxn = ScoreFunctionFactory::create_score_function( STANDARD_WTS, SCORE12_PATCH );
	}
	(*scorefxn)(pose);

	// set up BackrubMover and read from the database
	protocols::moves::BackrubMover backrubmover;
	backrubmover.branchopt().read_database();

	//clear segments and set the input pose
	backrubmover.clear_segments();
	pose::PoseOP input_poseOP ( new pose::Pose( pose ) );
	backrubmover.set_input_pose( input_poseOP );

	// set backrub segments
	Size const nres( pose.total_residue() );
	utility::vector1< bool > is_loop( nres, false );
	for ( Loops::const_iterator it=loops_.begin(), it_end=loops_.end(); it != it_end; ++it ) {
		for ( core::Size seg_start = it->start(); seg_start <= it->stop(); ++seg_start ) {
			is_loop[seg_start] = true;
			for ( core::Size seg_end = seg_start + 2; seg_end <= it->stop(); ++seg_end ) {
				id::AtomID start_atom_id = id::AtomID( pose.residue( seg_start ).atom_index("CA") , seg_start );
				id::AtomID end_atom_id = id::AtomID( pose.residue( seg_end ).atom_index("CA"), seg_end );
				// add segment to the mover
				backrubmover.add_segment( start_atom_id, end_atom_id, 0 );
			}
		}
	}

	// optimize branch angles and idealize side chains
	backrubmover.optimize_branch_angles( pose );
	(*scorefxn)(pose);

	// setup monte carlo
	float const init_temp( option[ OptionKeys::loops::refine_init_temp ]() );
	float const	final_temp( option[ OptionKeys::loops::refine_final_temp ]() );
	float const gamma = std::pow( (final_temp/init_temp), 1.0f/(outer_cycles*inner_cycles) );
	float temperature = init_temp;
	protocols::moves::MonteCarlo mc( pose, *scorefxn, temperature );
	mc.reset(pose);

	// setup PackerTask
	// default move map
	bool const fix_natsc = option[OptionKeys::loops::fix_natsc];
	// the following TaskFactory usage allows user-defined PackerTask creation on-demand
	using namespace pack::task;
	if ( task_factory == 0 ) {
		task_factory = new TaskFactory;
		// TaskOperations replace the following kind of code:
		// base_packer_task->initialize_from_command_line().or_include_current( true );
		task_factory->push_back( new operation::InitializeFromCommandline );
		task_factory->push_back( new operation::IncludeCurrent );
		if ( option[ OptionKeys::packing::resfile ].user() ) {
			// Note - resfile is obeyed, so use NATAA as default to maintain protocol behavior
			task_factory->push_back( new core::pack::task::operation::ReadResfile );
			tr << "Activating design" << std::endl;
		}
	}
	PackerTaskOP base_packer_task = task_factory->create_task_and_apply_taskoperations( pose );
	// this could also be handled/controlled by a [externally-defined] TaskOperation
	if ( redesign_loop ) {
		// allow design at loop positions
		for ( Size i=1; i<= nres; ++i ) {
			if ( !is_loop[i] ) base_packer_task->nonconst_residue_task( i ).restrict_to_repacking();
		}
	} else {
		// restrict to repacking at all positions
		base_packer_task->restrict_to_repacking();
	}
	base_packer_task->set_bump_check( true );

	// perform initial repack trial
	pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
	utility::vector1<bool> allow_repacked( nres, false );
	pose.update_residue_neighbors(); // to update 10A nbr graph
	select_loop_residues( pose, loops_, !fix_natsc, allow_repacked, 10.0 /* neighbor_cutoff */);
	this_packer_task->restrict_to_residues( allow_repacked );
	core::pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	std::string move_type = "repack";
	pose.update_residue_neighbors(); // to update 10A nbr graph
	mc.boltzmann( pose, move_type );

	std::string backrub_move_type = backrubmover.type();
	for (int i=1; i<=outer_cycles; ++i) {
		tr << "loopmover backrub outer refinement cycle " << i << std::endl;
		mc.score_function( *scorefxn );
		mc.recover_low( pose );

		for ( int j=1; j<=inner_cycles; ++j ) {
			temperature *= gamma;
			mc.set_temperature( temperature );

			for ( int trial = 1; trial <= ntrials_per_cycle; ++trial ) {
				backrubmover.apply(pose);
				mc.boltzmann(pose, backrub_move_type);
			}

			// main_repack_trial
			if ( (j%repack_period)==0 || j==inner_cycles ) {
				// repack trial
				//pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
				// DJM: try updating existing packer task
				utility::vector1<bool> allow_repacked( nres, false );
				select_loop_residues( pose, loops_, !fix_natsc, allow_repacked, 10.0 /* neighbor_cutoff */);
				this_packer_task->restrict_to_residues( allow_repacked );
				pack::pack_rotamers( pose, *scorefxn, this_packer_task );
				std::string move_type = "repack";
				mc.boltzmann( pose, move_type );
			}
		} //inner_cycle
	} //outer_cycle

	pose = mc.lowest_score_pose();

}

std::string
LoopMover_Refine_Backrub::get_name() const {
	return "LoopMover_Refine_Backrub";
}


} // namespace loops
} // namespace protocols
