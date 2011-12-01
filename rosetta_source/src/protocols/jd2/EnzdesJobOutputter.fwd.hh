// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/EnzdesJobOutputter.fwd.hh
/// @brief forward declaration for EnzdesJobOutputters
/// @author Florian Richter (floric@u.washington.edu), september 2010


#ifndef INCLUDED_protocols_jd2_EnzdesJobOutputter_fwd_hh
#define INCLUDED_protocols_jd2_EnzdesJobOutputter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace jd2 {

class EnzdesJobOutputter;
typedef utility::pointer::owning_ptr< EnzdesJobOutputter > EnzdesJobOutputterOP;

}//jd2
}//protocols



#endif //INCLUDED_protocols_jd2_EnzdesJobOutputter_FWD_HH
