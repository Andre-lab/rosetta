// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file ChangeFoldTreeMover.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_protocols_moves_ChangeFoldTreeMover_fwd_hh
#define INCLUDED_protocols_moves_ChangeFoldTreeMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

// Package headers

namespace protocols {
namespace moves {

class ChangeFoldTreeMover;
typedef utility::pointer::owning_ptr< ChangeFoldTreeMover > ChangeFoldTreeMoverOP;
typedef utility::pointer::owning_ptr< ChangeFoldTreeMover const > ChangeFoldTreeMoverCOP;

} // moves
} // rosetta


#endif
