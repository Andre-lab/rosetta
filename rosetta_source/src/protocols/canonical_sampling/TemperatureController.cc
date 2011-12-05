// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/TemperatureControllerMover.cc
/// @brief TemperatureController methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/TemperatureController.hh>


// protocols headers
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

#include <protocols/rosetta_scripts/util.hh>

//#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

// core headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cmath>

static basic::Tracer tr( "protocols.canonical_sampling.TemperatureController" );

namespace protocols {
namespace canonical_sampling {
using namespace core;


TemperatureController::TemperatureController() :
	protocols::canonical_sampling::ThermodynamicObserver()
{}

TemperatureController::TemperatureController(	TemperatureController const & other ) :
	//utility::pointer::ReferenceCount(),
	protocols::canonical_sampling::ThermodynamicObserver(other)
{}

///@brief make a move in temperature space depending on the current score
void
TemperatureController::observe_after_metropolis( protocols::canonical_sampling::MetropolisHastingsMover const& ) {
	Real const score( monte_carlo_->last_accepted_score() );
	temperature_move( score );
}

protocols::moves::MonteCarloCOP
TemperatureController::monte_carlo() const {
	return monte_carlo_;
}

void
TemperatureController::set_monte_carlo(
	protocols::moves::MonteCarloOP monte_carlo
)
{
	monte_carlo_ = monte_carlo;
}

protocols::moves::MonteCarloOP
TemperatureController::monte_carlo() {
	return monte_carlo_;
}

std::string
TemperatureController::get_name() const
{
	return "TemperatureController";
}

} //moves
} //protocols

