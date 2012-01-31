// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/EnergyPerResidueFilter.hh
/// @brief definition of filter class EnergyPerResidueFilter.
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_simple_filters_EnergyPerResidueFilter_hh
#define INCLUDED_protocols_simple_filters_EnergyPerResidueFilter_hh

#include <protocols/simple_filters/EnergyPerResidueFilter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>

namespace protocols {
namespace simple_filters {

class EnergyPerResidueFilter : public filters::Filter
{
public:
	EnergyPerResidueFilter() : filters::Filter( "EnergyPerResidue" ) {}

	EnergyPerResidueFilter( core::Size const resnum, core::scoring::ScoreFunctionCOP scorefxn,
	   core::scoring::ScoreType const score_type, core::Real const threshold,
	   bool const whole_interface = false, core::Size const rb_jump = 1,
	   core::Real const interface_distance_cutoff =  8.0 , bool const bb_bb = false );

	EnergyPerResidueFilter( EnergyPerResidueFilter const &init );
	bool apply( core::pose::Pose const & pose ) const;
	filters::FilterOP clone() const {
		return new EnergyPerResidueFilter( *this );
	}
	filters::FilterOP fresh_instance() const{
		return new EnergyPerResidueFilter();
	}

	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~EnergyPerResidueFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & );
private:
	core::Size resnum_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType score_type_;
	core::Real threshold_;
	bool whole_interface_;
	core::Size rb_jump_;
	core::Real interface_distance_cutoff_;
	bool bb_bb_;

};

}
}

#endif
