// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/RmsdFilter.hh
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_RmsdFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_RmsdFilter_hh


#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.fwd.hh>
#include <list>

#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

class RmsdFilter : public protocols::filters::Filter
{
public:
	RmsdFilter();
	RmsdFilter(
		std::list<core::Size> const selection,
		bool const superimpose,
		core::Real const threshold,
		core::pose::PoseOP reference_pose
	);
	bool apply( core::pose::Pose const & pose ) const;
	protocols::filters::FilterOP clone() const;
	protocols::filters::FilterOP fresh_instance() const{
		return new RmsdFilter();
	}
	void report( std::ostream & out, core::pose::Pose const & pose ) const;
	core::Real report_sm( core::pose::Pose const & pose ) const;
	core::Real compute( core::pose::Pose const & pose ) const;
	virtual ~RmsdFilter();
	void parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose );
private:
	std::list< core::Size > selection_;
	bool superimpose_, symmetry_;
	core::Real threshold_;
	core::pose::PoseOP reference_pose_;

	bool selection_from_segment_cache_;

};

} // filters
} // protein_interface_design
} // protocols

#endif //INCLUDED_protocols_protein_interface_design_filters_RmsdFilter_HH_

