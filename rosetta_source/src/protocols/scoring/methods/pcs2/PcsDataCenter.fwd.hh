// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
 //////////////////////////////////////////////
 /// @file protocols/scoring/methods/pcs2/PcsDataCenter.fwd.hh
 ///
 /// @authorsv Christophe Schmitz
 ///
 /// @last_modified February 2010
 ////////////////////////////////////////////////

#include <utility/pointer/owning_ptr.fwd.hh>

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsDataCenter_fwd_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsDataCenter_fwd_hh


namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

class PcsDataCenter;

typedef utility::pointer::owning_ptr< PcsDataCenter > PcsDataCenterOP;
typedef utility::pointer::owning_ptr< PcsDataCenter const > PcsDataCenterCOP;

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
