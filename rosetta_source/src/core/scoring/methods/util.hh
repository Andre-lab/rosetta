// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/methods/util.hh
/// @brief utility methods for scoring.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_methods_util_hh
#define INCLUDED_core_scoring_methods_util_hh

#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

namespace core {
namespace scoring {
namespace methods {

core::Real get_residue_weight_by_ss(
	char ss
);

bool residues_interact(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::Real const interaction_cutoff
);

bool atoms_interact(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	core::id::AtomID const & id1,
	core::id::AtomID const & id2,
	core::Real const interaction_cutoff
);

} // namespace methods
} // namespace scoring
} // namespace core

#endif
