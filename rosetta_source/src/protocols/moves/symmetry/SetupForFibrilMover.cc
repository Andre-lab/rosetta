// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author

// Unit headers
#include <protocols/moves/symmetry/SetupForFibrilMover.hh>

// Package headers
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>


// ObjexxFCL Headers

// C++ Headers

// Utility Headers
#include <basic/Tracer.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/symmetry/fibril_util.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {
namespace symmetry {

static basic::Tracer TR("protocols.moves.symmetry.SetupForFibrilMover");

SetupForFibrilMover::SetupForFibrilMover()
	: Mover("SetupForFibrilMover") {}

SetupForFibrilMover::~SetupForFibrilMover(){}

void
SetupForFibrilMover::apply( core::pose::Pose & pose )
{
	// If we are alredy symmetric do nothing
	if ( core::pose::symmetry::is_symmetric( pose ) ) return;
	protocols::moves::symmetry::make_symmetric_fibril( pose );
	assert( core::pose::symmetry::is_symmetric( pose ) );
}

std::string
SetupForFibrilMover::get_name() const {
	return "SetupForFibrilMover";
}

void
SetupForFibrilMover::align(
	core::pose::Pose & pose,
	core::pose::Pose & monomer_pose,
	protocols::loops::Loops core,
  protocols::loops::Loops ref_core
)
{
	// If we are alredy symmetric do nothing
	if ( core::pose::symmetry::is_symmetric( pose ) ) return;
	std::cout<<"align: from "<<core<<" to " <<ref_core<<std::endl;
	protocols::moves::symmetry::superimpose_pose_on_subset_bb( pose, monomer_pose, core, ref_core );
}

} // symmetry
} // moves
} // protocols
