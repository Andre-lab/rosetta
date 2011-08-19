// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file FoldandDockClaimer
/// @brief Fold-and-dock
/// @author Ingemar Andre

// Unit Headers
#include <protocols/topology_broker/FoldandDockClaimer.hh>
#include <protocols/moves/symmetry/SymFoldandDockRbTrialMover.hh>
#include <protocols/moves/symmetry/SymFoldandDockSlideTrialMover.hh>
#include <protocols/symmetric_docking/SymDockingInitialPerturbation.hh>
#include <protocols/moves/symmetry/SymFoldandDockMoveRbJumpMover.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetricConformation.hh>

// Package Headers
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <basic/Tracer.hh>

// Utility header
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <numeric/random/random.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>


//Auto Headers
#include <core/conformation/Conformation.hh>
#include <protocols/moves/MoverContainer.hh>


// Project Headers

static basic::Tracer tr("protocols.topo_broker.fold_and_dock",basic::t_info);
static numeric::random::RandomGenerator RG(24278234);

namespace protocols {
namespace topology_broker {

using namespace core;


FoldandDockClaimer::FoldandDockClaimer() {}

FoldandDockClaimer::FoldandDockClaimer( pose::Pose const& input_pose ) :
	input_pose_(input_pose)
{}

	//clone
TopologyClaimerOP
FoldandDockClaimer::clone() const {
	return new FoldandDockClaimer( *this );
}

///@brief type() is specifying the output name of the TopologyClaimer
std::string
FoldandDockClaimer::type() const {
	return _static_type_name();
}

std::string
FoldandDockClaimer::_static_type_name() {
	return "FoldandDockClaimer";
}

void
FoldandDockClaimer::add_mover(
    moves::RandomMover& random_mover,
		core::pose::Pose const& /*pose*/,
		abinitio::StageID stageID,  /*abinitio sampler stage */
		core::scoring::ScoreFunction const& scorefxn,
		core::Real /*progress  progress within stage */
)
{
	using namespace basic::options;

	moves::MoverOP move_anchor_mover =	new moves::symmetry::SymFoldandDockMoveRbJumpMover;
	moves::MoverOP rb_trial_mover =	(stageID==abinitio::STAGE_4) ?
		new moves::symmetry::SymFoldandDockRbTrialMover( &scorefxn, true ) :
		new moves::symmetry::SymFoldandDockRbTrialMover( &scorefxn );  // smooth RB moves in stage 4
	moves::MoverOP slide_mover = new moves::symmetry::SymFoldandDockSlideTrialMover;
	core::Real move_anchor_weight(1.0),
	           rb_weight(option[ OptionKeys::fold_and_dock::rigid_body_frequency ]()),
	           slide_weight(option[ OptionKeys::fold_and_dock::slide_contact_frequency ]());

	random_mover.add_mover( move_anchor_mover, move_anchor_weight );
	random_mover.add_mover( rb_trial_mover, rb_weight );
	random_mover.add_mover( slide_mover, slide_weight );
}

void FoldandDockClaimer::initialize_dofs(
	core::pose::Pose& pose,
	DofClaims const& init_dofs,
	DofClaims& /*failed_to_init*/ ) {

	using namespace core::conformation::symmetry;

	// Setup symmetry if we have nit already done it
	// slide chains into contact
	protocols::moves::symmetry::SetupForSymmetryMoverOP setup_mover = new
		protocols::moves::symmetry::SetupForSymmetryMover;
	setup_mover->slide_into_contact(true);
	if ( !core::pose::symmetry::is_symmetric( pose ) ) {
		setup_mover->apply( pose ); // calls SymDockingInitialPerturbation
		assert( core::pose::symmetry::is_symmetric( pose ) );
		// Save the pose into input pose
		input_pose_ = pose;
	} else {
		input_pose_ = pose;
		// Randomize the rigid body
		protocols::symmetric_docking::SymDockingInitialPerturbation initial( true /*slide into contact*/ );
		initial.apply( pose );
	}

	// Setup the movemap
	//SymmetricConformation & symm_conf (dynamic_cast<SymmetricConformation & > ( pose.conformation()) );
	kinematics::MoveMapOP movemap = new kinematics::MoveMap();
	movemap->set_bb( true );
	movemap->set_jump( false );
	core::pose::symmetry::make_symmetric_movemap( pose, *movemap );

	for ( DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
          it != eit; ++it ) {
		if ( (*it)->owner()==this ) {
			(*it)->toggle( *movemap, true );
		}
	}
}

void FoldandDockClaimer::generate_claims( DofClaims& new_claims ) {
	// Set all cuts to real cuts. We don't want to close any of them...
	utility::vector1< int > cuts( input_pose_.conformation().fold_tree().cutpoints() );
	for ( Size i = 1; i <= cuts.size(); ++i ) {
		new_claims.push_back( new CutClaim( this, cuts[i], DofClaim::INIT /* for now... eventually CAN_INIT ? */ ) );
	}
}



} //topology_broker
} //protocols
