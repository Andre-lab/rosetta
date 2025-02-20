// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/pack_missing_sidechains.hh
/// @brief  header for subroutine to run rotamer trials on residues missing sidechain density in a PDB
/// @author Steven Lewis smlewi@gmail.com
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitue.org) to allow multi-threaded interaction graph setup.

#ifndef INCLUDED_core_pack_pack_missing_sidechains_hh
#define INCLUDED_core_pack_pack_missing_sidechains_hh

// Project headers
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>

#ifdef MULTI_THREADED
#include <core/types.hh>
#endif

namespace core {
namespace pack {

/// @brief This function runs rotamer trials on residues missing sidechain density (as described
/// by the AtomID_Mask).
/// @details this function will run rotamer trials on sidechains with missing density.  It
/// first sets up a PackerTask with repacking freedom for residues with sidechain missing atoms
/// in the missing AtomID_Mask, then runs rotamer_trials.  This function is smart enough to
/// ignore missing virtual atoms
/// @note In multi-threaded builds, this function takes an extra parameter for
/// the number of threads to request, for parallel interaction graph precomputation.
void
pack_missing_sidechains(
	pose::Pose & pose,
	id::AtomID_Mask const& missing
#ifdef MULTI_THREADED
	,
	core::Size const threads_to_request=0
#endif
);

/// @brief return vector of bools with true for each residue that has >=1 atom in to_repack that is not VIRTUAL, ORBS or LPbb
bool figure_out_repackable_residues(
	core::pose::Pose & pose,
	core::id::AtomID_Mask const& to_repack,
	utility::vector1_bool& repackable
);


} //namespace pack
} //namespace core

#endif //INCLUDED_core_pack_pack_missing_sidechains_HH
