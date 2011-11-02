// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
// rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//		under license.
// (c) The Rosetta software is developed by the contributing members of the
//		Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org.
//		Questions about this can be
// (c) addressed to University of Washington UW TechTransfer,
//		email: license@u.washington.edu.

/// @file docking_min_protocol
/// @brief
/// @author Robin A Thottungal (raugust1@jhu.edu)

// Unit Headers
#include <protocols/docking/DockMinMover.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialMover.hh>

#include <utility/vector1.hh>


using namespace protocols::moves;
using namespace core;

namespace protocols{
namespace docking{

//@TODO create default values in empty constructor
DockMinMover::DockMinMover() : DockingHighRes()
{
	//need to set this up with default values;
	set_default();
	mc_= new MonteCarlo( *scorefxn(), 0.8 );
	MinMoverOP min_mover =new moves::MinMover
				( movemap_, scorefxn(), min_type_, min_tolerance_, nb_list_ );
	minimize_trial_ = new TrialMover( min_mover, mc_ );
}

DockMinMover::DockMinMover(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionCOP scorefxn
) : DockingHighRes(movable_jumps, scorefxn, scorefxn) {
	//need to set this up with default values;
	set_default();
	mc_= new MonteCarlo( *scorefxn(), 0.8 );
	MinMoverOP min_mover =new moves::MinMover
				( movemap_, scorefxn(), min_type_, min_tolerance_, nb_list_ );
	minimize_trial_ = new TrialMover( min_mover, mc_ );
}

	
	
//JQX: made a new constructor, which can take the mc_ object	
DockMinMover::DockMinMover(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionCOP scorefxn,
	moves::MonteCarloOP mc
) : DockingHighRes(movable_jumps, scorefxn, scorefxn) {
	//need to set this up with default values;
	set_default();
	mc_=mc;
	MinMoverOP min_mover =new moves::MinMover
			    ( movemap_, scorefxn(), min_type_, min_tolerance_, nb_list_ );
	minimize_trial_ = new TrialMover( min_mover, mc_ );
}
	
	
	
	
	
DockMinMover::DockMinMover(
	DockJumps const movable_jumps,
	core::kinematics::MoveMapOP movemap,
	core::scoring::ScoreFunctionCOP scorefxn,
	std::string min_type,
	core::Real min_tolerance,
	bool nb_list,
	moves::MonteCarloOP mc
): DockingHighRes(movable_jumps, scorefxn, scorefxn) {
	movemap_=movemap;
	min_type_=min_type;
	min_tolerance_=min_tolerance;
	nb_list_=nb_list;
	mc_=mc;
	MinMoverOP min_mover =new moves::MinMover
				( movemap_, scorefxn(), min_type_, min_tolerance_, nb_list_ );
	minimize_trial_ = new TrialMover( min_mover, mc_ );
}

DockMinMover::~DockMinMover(){}

void DockMinMover::set_default() {
	//sets up default movemap
	movemap_ = new kinematics::MoveMap();
	movemap_->set_chi( false );
	movemap_->set_bb( false );
	for( DockJumps::const_iterator it = movable_jumps().begin(); it != movable_jumps().end(); ++it ) {
		movemap_->set_jump( *it, true );
	}

	//sets up minimization parameters
//	min_tolerance_ = 1.0; /////was 0.01, in r++ docking, it is actually 1.0!! with 0.02 as the "tight" tolerance
//	JQX commented the line right above, the Legacy code use 0.01. Just want to match it
	min_tolerance_ = 0.01; //JQX added this line
	min_type_ = std::string( "dfpmin_armijo_nonmonotone" );
	nb_list_ = true;
}

void DockMinMover::apply(core::pose::Pose & pose){
	// since mc_ object is created without a pose, this reset is required!
	( *scorefxn() )( pose );
//	mc_->reset(pose);   // JQX: I think it is not necessary
	minimize_trial_->apply(pose);
}

std::string DockMinMover::get_name() const {
	return "DockMinMover";
}

}
}
