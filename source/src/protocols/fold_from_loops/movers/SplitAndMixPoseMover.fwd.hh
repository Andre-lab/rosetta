// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SplitAndMixPoseMover.hh
/// @brief  Splits a Pose according to a ResidueSelector
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL


#ifndef INCLUDED_protocols_fold_from_loops_movers_SplitAndMixPoseMover_FWD_hh
#define INCLUDED_protocols_fold_from_loops_movers_SplitAndMixPoseMover_FWD_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace fold_from_loops {
namespace movers {

class SplitAndMixPoseMover;
typedef utility::pointer::shared_ptr< SplitAndMixPoseMover >        SplitAndMixPoseMoverOP;
typedef utility::pointer::shared_ptr< SplitAndMixPoseMover const >  SplitAndMixPoseMoverCOP;

}
}
}

#endif
