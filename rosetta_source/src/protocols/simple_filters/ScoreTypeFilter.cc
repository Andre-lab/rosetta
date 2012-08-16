// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ScoreTypeFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)


//Unit Headers
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>
#include <utility/exit.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/format.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
//Auto Headers
#include <protocols/simple_filters/ScoreTypeFilter.hh>


namespace protocols{
namespace simple_filters {

using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::fmt;

static basic::Tracer score_type_filter_tracer( "protocols.simple_filters.ScoreTypeFilter" );

protocols::filters::FilterOP
ScoreTypeFilterCreator::create_filter() const { return new ScoreTypeFilter; }

std::string
ScoreTypeFilterCreator::keyname() const { return "ScoreType"; }

ScoreTypeFilter::ScoreTypeFilter( core::scoring::ScoreFunctionCOP scorefxn, core::scoring::ScoreType const score_type, core::Real const score_type_threshold ) : Filter( "ScoreType" ) {
	score_type_ = score_type;
	score_type_threshold_ = score_type_threshold;
	scorefxn_ = scorefxn->clone();
}

ScoreTypeFilter::~ScoreTypeFilter() {}

void
ScoreTypeFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn", "score12" ) );
	// scorefxn_ = new ScoreFunction( *(data.get< ScoreFunction * >( "scorefxns", scorefxn_name )) );
	scorefxn_ = data.get< ScoreFunction * >( "scorefxns", scorefxn_name )->clone();

	score_type_ = core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	if( ! tag->hasOption( "threshold" ) ) utility_exit_with_message("Must specify 'threshold' for ScoreTypeFilter.");
	score_type_threshold_ = tag->getOption<core::Real>( "threshold" );

	score_type_filter_tracer<<"ScoreType filter for score_type "<<score_type_<<" with threshold "<<score_type_threshold_<<std::endl;
}

bool
ScoreTypeFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const score( compute( pose ) );
	score_type_filter_tracer<<"score "<<core::scoring::ScoreTypeManager::name_from_score_type( score_type_ )<<" is "<<score<<". ";
	if( score <= score_type_threshold_ ) {
		score_type_filter_tracer<<"passing." << std::endl;
		return true;
	}
	else {
		score_type_filter_tracer<<"failing."<<std::endl;
		return false;
	}
}

void
ScoreTypeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out<<"Weighted score of "<<core::scoring::ScoreTypeManager::name_from_score_type( score_type_ )<<" "<<compute( pose )<<'\n';
}

core::Real
ScoreTypeFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
ScoreTypeFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::pose;
	using namespace core::scoring;

	PoseOP in_pose = new Pose( pose );

	// make sure that scoring weights are compatible with pose's residue type set
	// check centroid case
	if( ( (*scorefxn_)[fa_rep] == 0.0 && (*scorefxn_)[fa_atr] == 0.0 ) // full atom terms are off
				&& ( (*scorefxn_)[interchain_vdw] > 0.0 || (*scorefxn_)[vdw] > 0.0)  ) // a centroid term is on
		{
			if( in_pose->is_fullatom() ) { // but pose is full atom
			core::util::switch_to_residue_type_set( *in_pose, core::chemical::CENTROID );
		}
	}
	else { // full atom case
		if( in_pose->is_centroid() ) { // but pose is centroid
			core::util::switch_to_residue_type_set( *in_pose, core::chemical::FA_STANDARD );
		}
	}

	(*scorefxn_)( *in_pose );
	/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies( *in_pose );
	core::Real const weight( (*scorefxn_)[ ScoreType( score_type_ ) ] );
	core::Real const score( in_pose->energies().total_energies()[ ScoreType( score_type_ ) ]);
	if( score_type_ == total_score ) return( score );
	core::Real const weighted_score( weight * score );
	return( weighted_score );
}

}
}
