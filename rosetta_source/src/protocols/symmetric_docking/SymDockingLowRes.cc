// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file DockingLowRes
/// @brief protocols that are specific to docking low resolution
/// @detailed This is to a very large extent a copy of the docking
/// @detailed protocol. Should derive out of that class instead.
/// @author Ingemar Andre

#include <protocols/symmetric_docking/SymDockingLowRes.hh>

// Rosetta Headers
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/kinematics/FoldTree.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/OutputMovers.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>

//Utility Headers
// AUTO-REMOVED #include <numeric/conversions.hh>

#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <basic/Tracer.hh>

#include <protocols/moves/MoverContainer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.symetric_docking.SymDockingLowRes");

//     originally from dock_structure.cc Jeff Gray April 2001
//
//

using namespace core;

namespace protocols {
namespace symmetric_docking {

	// constructor with arguments
	SymDockingLowRes::SymDockingLowRes(
		core::scoring::ScoreFunctionCOP scorefxn_in
	) : Mover(), scorefxn_(scorefxn_in)
	{
		moves::Mover::type( "SymDockingLowRes" );
	}

	SymDockingLowRes::~SymDockingLowRes(){}

	moves::MoverOP
	SymDockingLowRes::clone() const {
		return new SymDockingLowRes( *this );
	}

	void
	SymDockingLowRes::set_default( core::pose::Pose & pose ) {
		using namespace basic::options;

		// sets up the stuff in pose
		(*scorefxn_)( pose );

		// cycles
		inner_cycles_ = option[ OptionKeys::docking::docking_centroid_inner_cycles ]();
		outer_cycles_ = option[ OptionKeys::docking::docking_centroid_outer_cycles ]();

		if ( option[ OptionKeys::docking::dock_mcm_trans_magnitude ].user() ) {
			trans_magnitude_ = option[ OptionKeys::docking::dock_mcm_trans_magnitude ]();
		} else {
			trans_magnitude_ = 1.5;
		}

		if ( option[ OptionKeys::docking::dock_mcm_rot_magnitude ].user() ) {
			rot_magnitude_ = option[ OptionKeys::docking::dock_mcm_rot_magnitude ]();
		} else {
			rot_magnitude_ = 4;
		}

		chi_ = false;
		bb_ = false;

		temperature_ = 0.8;

		nb_list_ = true; /// not sure if this should be true or not
		accept_rate_ = 0.0;

		set_default_mc( pose );
		set_default_move_map( pose );
	  set_default_protocol( pose );
	}

	moves::MonteCarloOP
	SymDockingLowRes::get_mc() { return mc_; }

	void
	SymDockingLowRes::set_default_mc( pose::Pose & pose ) {
		// create the monte carlo object and movemap
		mc_ = new moves::MonteCarlo( pose, *scorefxn_, temperature_ );
	}

void SymDockingLowRes::set_default_move_map( pose::Pose & pose ) {
	using namespace core::conformation::symmetry;

	movemap_ = new kinematics::MoveMap();
	movemap_->set_bb( bb_ );
	movemap_->set_chi( chi_ );
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap_ );

}

void SymDockingLowRes::set_default_protocol( pose::Pose & pose ){
	using namespace moves;
	using namespace conformation::symmetry;

	 assert( core::pose::symmetry::is_symmetric( pose ));
  SymmetricConformation & symm_conf (
        dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

  std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );

	rb_mover_ = new RigidBodyDofSeqPerturbMover( dofs , rot_magnitude_, trans_magnitude_ );

	docking_lowres_protocol_ = new SequenceMover;
	docking_lowres_protocol_->add_mover( rb_mover_ );

}
////////////////////////////////////////////////////////////////////////////////
/// @begin DockingLowRes.apply
///
/// @brief Perform several cycles of rigid-body Monte Carlo moves
///       and adapt the step size.
/// @detailed
///
/// @remarks
///       currently used only in the low-resolution step (centroid mode)
///
/// @references pose_docking_centroid_rigid_body_adaptive from pose_docking.cc and
///				rigid_body_MC_cycle_adaptive from dock_structure.cc
///
/// @authors Monica Berrondo October 22 2007
///
/// @last_modified October 22 2007
/////////////////////////////////////////////////////////////////////////////////
void SymDockingLowRes::apply( core::pose::Pose & pose )
{
	using namespace scoring;

	TR << "in DockingLowRes.apply\n";

	set_default( pose );

	TR << "::::::::::::::::::Centroid Rigid Body Adaptive:::::::::::::::::::\n";

	for ( int i=1; i<=outer_cycles_; ++i) {
		rigid_body_trial( pose );
		if ( accept_rate_ < 0.5 ) {
			trans_magnitude_ *= 0.9;
			rot_magnitude_ *= 0.9;
		} else {
			trans_magnitude_ *= 1.1;
			rot_magnitude_ *= 1.1;
		}
		// if ( jump_out_check() ) return;
	}
	mc_->recover_low( pose );
	TR.flush();
	//pose.energies().show( std::cout );
}

std::string
SymDockingLowRes::get_name() const {
	return "SymDockingLowRes";
}

////////////////////////////////////////////////////////////////////////////////
/// @begin rigid_body_trial
///
/// @brief Perform a cycle of rigid-body Monte Carlo moves
///
/// @detailed  Performs a number (nattempts) of MC rigid-body moves
///       (of size trans_magnitude, rot_magnitude). The number of successful
///       attempts is stored in accept_rate_ and used in adaptive trials.
///
/// @remarks the success_rate defines
///       whether the translation/rotation size is increased or decreased for
///       the next cycle.
///       currently used only in the low-resolution step (centroid mode)
///
/// @references pose_docking_rigid_body_trial from pose_docking.cc and
///				rigid_body_MC_cycle from dock_structure.cc
///
/// @authors Monica Berrondo October 22 2007
///
/// @last_modified October 22 2007
/////////////////////////////////////////////////////////////////////////////////
void SymDockingLowRes::rigid_body_trial( core::pose::Pose & pose )
{
	using namespace moves;

	PDBDumpMoverOP dump = new PDBDumpMover("lowres_cycle_");
//	dump->apply( pose );
	MCShowMoverOP mc_show = new MCShowMover( mc_ );
//	mc_show->apply( pose );

	rb_mover_->rot_magnitude( rot_magnitude_ );
	rb_mover_->trans_magnitude( trans_magnitude_ );

	TrialMoverOP rb_trial = new TrialMover( docking_lowres_protocol_, mc_ );
	rb_trial->keep_stats_type( moves::accept_reject );

	RepeatMoverOP rb_cycle = new RepeatMover( rb_trial, inner_cycles_ );

	rb_cycle->apply( pose );

	pose = mc_->lowest_score_pose();
	//pose.energies().show( std::cout );
	mc_->reset( pose );

	accept_rate_ = rb_trial->acceptance_rate();
}

} // namespace docking
} // namespace protocols
