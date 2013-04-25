// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/coupledSidechainProtocol.fwd.hh
/// @brief
/// @author


#ifndef INCLUDED_devel_coupled_sidechains_coupledSidechainProtocol_fwd_hh
#define INCLUDED_devel_coupled_sidechains_coupledSidechainProtocol_fwd_hh

#include <utility/pointer/owning_ptr.hh>

// Package headers

namespace devel {
namespace coupled_sidechains {

class coupledSidechainProtocol;
typedef utility::pointer::owning_ptr< coupledSidechainProtocol > coupledSidechainProtocolOP;
typedef utility::pointer::owning_ptr< coupledSidechainProtocol const > coupledSidechainProtocolCOP;

} // coupled_sidechains
} // devel


#endif  //INCLUDED_devel_coupled_sidechains_coupledSidechainProtocol_fwd_HH
