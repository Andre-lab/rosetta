// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /devel/replica_docking/TempWeightedMetropolisHastingsMover.hh
/// @brief inherited from MetropolisHastingsMover to make each replica holds several rb_mover
/// (with different step size), then apply rb_mover according to the current temp_level
/// @author Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_TempWeightedMetropolisHastingsMover_hh
#define INCLUDED_devel_replica_docking_TempWeightedMetropolisHastingsMover_hh

// Unit Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/Mover.fwd.hh>

// zhe
#include <devel/replica_docking/TempInterpolator.fwd.hh>

// Project Headers
#include <numeric/random/WeightedSampler.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

#include <protocols/canonical_sampling/ThermodynamicMover.fwd.hh>

#ifdef WIN32
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#endif


namespace devel {
namespace replica_docking {

/// @details
class TempWeightedMetropolisHastingsMover : public protocols::canonical_sampling::MetropolisHastingsMover {
	typedef utility::vector1< devel::replica_docking::TempInterpolatorBaseOP > Interpolators;
	typedef utility::vector1< core::Real > Weights;
	typedef utility::vector1< core::Size > GridCoord;
	typedef protocols::canonical_sampling::MetropolisHastingsMover Parent;
public:

	TempWeightedMetropolisHastingsMover();

	TempWeightedMetropolisHastingsMover(
		TempWeightedMetropolisHastingsMover const & tw_metropolis_hastings_mover
	);

	protocols::moves::MoverOP
	clone() const override;


	protocols::moves::MoverOP
	fresh_instance() const override;


	std::string
	get_name() const override;


	protocols::canonical_sampling::ThermodynamicMoverOP
	random_mover() const override;


	void
	add_mover(
		protocols::canonical_sampling::ThermodynamicMoverOP mover,
		core::Real weight,
		utility::tag::TagCOP const& subtag
	) override ;

protected:

private:

	//helper
	Interpolators weight_contro_1_;
	Interpolators weight_contro_2_;
	Weights overall_weights_;

	mutable core::Size last_temp_level_in_random_mover_;
	mutable numeric::random::WeightedSampler current_weighted_sampler_; //overwrites base-class member
};

} //namespace
} //namespace

#endif //INCLUDED_
