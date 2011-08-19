// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// ThiS file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SaneDockingProtocol.cc
/// @author James Thompson

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/docking/stateless/SaneDockingProtocol.hh>

namespace protocols {
namespace docking {
namespace stateless {

void
SaneDockingProtocol::apply( core::pose::Pose & pose ) {
	using namespace core::pose;
	if ( !get_input_pose() ) {
		PoseCOP input_pose_op = new core::pose::Pose(pose);
		set_input_pose (input_pose_op);
	}
	if ( !get_native_pose() ) {
		PoseCOP native_pose_op = new core::pose::Pose(pose);
		set_native_pose(native_pose_op);
	}
	DockingProtocol::apply(pose);
}

std::string
SaneDockingProtocol::get_name() const {
	return "DockingProtocol";
}

} // stateless
} // docking
} // protocols
