// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/CompositionMover.cc
/// @brief
/// @author

// Unit Headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/CompositionMover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

CompositionMover::CompositionMover()
	: Mover("CompositionMover")
{}

void
CompositionMover::apply( core::pose::Pose & pose ) {
	typedef utility::vector1< MoverOP >::iterator iter;
	for ( iter it = movers_.begin(), end = movers_.end(); it != end; ++it ) {
		(*it)->apply( pose );
	}
}

std::string
CompositionMover::get_name() const {
	return "CompositionMover";
}

void CompositionMover::clear() {
	movers_.clear();
}

void CompositionMover::add_mover( MoverOP m ) {
	movers_.push_back( m );
}

utility::vector1< MoverOP > CompositionMover::get_movers() {
	return movers_;
}

} // moves
} // protocols

