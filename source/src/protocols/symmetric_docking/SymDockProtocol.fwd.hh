// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/symmetric_docking/SymDockProtocol.fwd.hh
///
/// @brief
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_symmetric_docking_SymDockProtocol_fwd_hh
#define INCLUDED_protocols_symmetric_docking_SymDockProtocol_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace symmetric_docking {


class SymDockProtocol; // fwd declaration
typedef utility::pointer::owning_ptr< SymDockProtocol > SymDockProtocolOP;
typedef utility::pointer::owning_ptr< SymDockProtocol const > SymDockProtocolCOP;


} // namespace symmetric_docking
} // namespace protocols

#endif // INCLUDED_protocols_symmetric_docking_SymDockProtocol_FWD_HH

