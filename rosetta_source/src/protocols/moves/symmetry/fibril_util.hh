// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file  protocols/moves/symmetry/fibril_util.hh
/// @brief utility functions for handling with symmetric fibril
/// @author Lin Jiang

#ifndef INCLUDED_protocols_moves_symmetry_fibril_util_hh
#define INCLUDED_protocols_moves_symmetry_fibril_util_hh


// Unit headers
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.fwd.hh>

// AUTO-REMOVED #include <protocols/loops/Loops.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {
namespace symmetry {

void
reorient_extended_fibril(
	core::conformation::Conformation & src_conformation,
  core::conformation::symmetry::SymmData & symmdata
);

void
make_symmetric_fibril(
	core::pose::Pose & pose
);

void
superimpose_pose_on_subset_bb(
	core::pose::Pose& pose,
  core::pose::Pose& ref_pose,
  protocols::loops::Loops core,
  protocols::loops::Loops ref_core
);


} // symmetry
} // moves
} // protocols



#endif
