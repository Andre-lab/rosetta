// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Lei Shi <shilei@u.washington.edu>
/// @date 3/19/2013

#ifndef INCLUDED_protocols_protein_interface_design_movers_DockWithHotspotMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_DockWithHotspotMover_hh

#include <protocols/protein_interface_design/movers/DockWithHotspotMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/AA.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

#include <utility/vector1.hh>


// Utility headers

// C++ headers

// Unit headers

namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief A mover to mutate a single residue
class DockWithHotspotMover : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	///@brief default ctor
	DockWithHotspotMover();
	virtual ~DockWithHotspotMover();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return (protocols::moves::MoverOP( new protocols::protein_interface_design::movers::DockWithHotspotMover( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new DockWithHotspotMover );
	}

	void parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

private:
  utility::vector1< std::string > hotspot_filenames_;
  utility::vector1< core::Real > hotspot_distcb_weight_;
	core::Real hotspot_score_weight_;
	core::Real centroidscore_filter_;
	core::Real hotspotcst_filter_;
};

} //movers
} //protein_interface_design
} //protocols

#endif //INCLUDED_protocols_protein_interface_design_movers_DockWithHotspotMover_HH_

