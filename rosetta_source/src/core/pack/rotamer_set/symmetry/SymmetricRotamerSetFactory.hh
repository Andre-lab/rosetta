// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/pack/rotamer_set/RotamerSetFactory.hh
/// @brief  Residue Set Factory class for symmetric packing
/// @author  Ingemar Andre

#ifndef INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSetFactory_hh
#define INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSetFactory_hh

// Unit headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

class SymmetricRotamerSetFactory : public RotamerSetFactory
{
public:
	RotamerSetOP create_rotamer_set( conformation::Residue const & );

};

}
}
}
}

#endif
