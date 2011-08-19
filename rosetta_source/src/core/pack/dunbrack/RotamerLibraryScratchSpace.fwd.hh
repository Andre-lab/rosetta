// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dunbrack/RotamerLibraryScratchSpace.fwd.hh
/// @brief  Forward declaration of scratch space class for Dunbrack rotamer library
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_RotamerLibraryScratchSpace_fwd_hh
#define INCLUDED_core_pack_dunbrack_RotamerLibraryScratchSpace_fwd_hh

// Package Headers
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/fixedsizearray1.hh>

//Auto Headers
#include <utility/fixedsizearray1.fwd.hh>
#include <iterator>


namespace core {
namespace pack {
namespace dunbrack {

Size const DUNBRACK_MAX_BBTOR = 3;
Size const DUNBRACK_MAX_SCTOR = 4;

typedef utility::fixedsizearray1< Real, DUNBRACK_MAX_SCTOR > Real4;
typedef utility::fixedsizearray1< Size, DUNBRACK_MAX_SCTOR > Size4;

typedef utility::fixedsizearray1< Real, DUNBRACK_MAX_BBTOR > Real3;
typedef utility::fixedsizearray1< Size, DUNBRACK_MAX_BBTOR > Size3;


class RotamerLibraryScratchSpace;

typedef utility::pointer::owning_ptr< RotamerLibraryScratchSpace > RotamerLibraryScratchSpaceOP;
typedef utility::pointer::owning_ptr< RotamerLibraryScratchSpace const > RotamerLibraryScratchSpaceCOP;


}
}
}

#endif

