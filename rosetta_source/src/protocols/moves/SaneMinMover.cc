// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SaneMinMover.cc
/// @brief
/// @author James Thompson

#include <protocols/moves/SaneMinMover.hh>
#include <protocols/moves/SaneMinMoverCreator.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/prof.hh>
// AUTO-REMOVED #include <utility/string_util.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>


namespace protocols {
namespace moves {

static basic::Tracer TR("protocols.moves.SaneMinMover");

std::string
SaneMinMoverCreator::keyname() const {
	return SaneMinMoverCreator::mover_name();
}

protocols::moves::MoverOP
SaneMinMoverCreator::create_mover() const {
	return new SaneMinMover;
}

std::string
SaneMinMoverCreator::mover_name() {
	return "SaneMinMover";
}

// default constructor
SaneMinMover::SaneMinMover() : Mover("SaneMinMover") {
	set_defaults_();
}

SaneMinMover::SaneMinMover( std::string const & name ) :
Mover(name) {
	set_defaults_();
}

SaneMinMover::~SaneMinMover() {}

// constructor with arguments
SaneMinMover::SaneMinMover(
	core::kinematics::MoveMapOP movemap_in,
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::optimization::MinimizerOptionsOP min_options_in,
	bool cartesian_in
) : Mover("SaneMinMover"),
		movemap_(movemap_in),
		scorefxn_(scorefxn_in),
		min_options_(min_options_in),
		cartesian_(cartesian_in)
{}

bool SaneMinMover::cartesian() const {
	return cartesian_;
}

core::kinematics::MoveMapOP SaneMinMover::move_map() const {
	return movemap_;
}

core::scoring::ScoreFunctionOP SaneMinMover::score_function() const {
	return scorefxn_;
}

core::optimization::MinimizerOptionsOP SaneMinMover::min_options() const {
	return min_options_;
}

void
SaneMinMover::apply( core::pose::Pose & pose ) {
	using namespace core::optimization;

	PROF_START( basic::MINMOVER_APPLY );
	(*scorefxn_)(pose);
	if (cartesian()) {
		CartesianMinimizer minimizer;
		minimizer.run( pose, *movemap_, *scorefxn_, *min_options_ );
	} else {
		AtomTreeMinimizer minimizer;
		minimizer.run( pose, *movemap_, *scorefxn_, *min_options_ );
	}
	PROF_STOP( basic::MINMOVER_APPLY );

	scorefxn_->show(TR, pose);
	TR.flush();
}

std::string
SaneMinMover::get_name() const {
	return SaneMinMoverCreator::mover_name();
}

MoverOP SaneMinMover::clone() const { return new SaneMinMover( *this ); }

void SaneMinMover::set_defaults_() {
	cartesian_   = false;
	movemap_     = new core::kinematics::MoveMap;
	scorefxn_    = core::scoring::getScoreFunction();
	min_options_ = new core::optimization::MinimizerOptions( "dfpmin", 1e-2, true, false, false );
}

} // moves
} // protocols
