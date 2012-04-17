// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/HeavyAtomFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/filters/HeavyAtomFilter.hh>
#include <protocols/filters/HeavyAtomFilterCreator.hh>


#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <utility/exit.hh>

//Auto Headers


namespace protocols {
namespace filters {

static basic::Tracer heavy_atom_tracer( "protocols.filters.HeavyAtomFilter" );

bool
HeavyAtomFilter::apply( core::pose::Pose const & pose ) const {
	assert(chain_.size()==1 );
	assert(heavy_atom_limit_ >0 );
	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const start = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);

	if(	core::pose::num_heavy_atoms(start,end,pose) > heavy_atom_limit_ ){
		heavy_atom_tracer<< "Reached heavy atom limit"<< std::endl;
		return false;
	}
	return true;
}

void
HeavyAtomFilter::parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	if ( tag->getName() != "HeavyAtom" ) {
		heavy_atom_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( ! (tag->hasOption("chain") && tag->hasOption("heavy_atom_limit") ) ){
		utility_exit_with_message("HeavyAtom filter needs a 'chain' and a 'heavy_atom_limit' option");
	}
	chain_ = tag->getOption<std::string>("chain");

}

FilterOP
HeavyAtomFilterCreator::create_filter() const { return new HeavyAtomFilter; }

std::string
HeavyAtomFilterCreator::keyname() const { return "HeavyAtom"; }



} // filters
} // protocols
