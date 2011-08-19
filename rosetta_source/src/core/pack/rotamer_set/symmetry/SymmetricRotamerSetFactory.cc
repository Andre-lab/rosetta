// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.hh
/// @brief  Residue Set Factory class for symmetric packing
/// @author Ingemar Andre

// Unit header
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.hh>

// Package headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/exit.hh>

// STL Headers
#include <iostream>

namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

RotamerSetOP
SymmetricRotamerSetFactory::create_rotamer_set( conformation::Residue const & res )
{
	if ( res.is_protein() ) // This check will be removed when we get rotamers for NAs and Ligands online
	{
		return new SymmetricRotamerSet_();
	}
	else
	{
		//std::cout << "[ WARNING ] PB HACK -- SHOULD DIE HERE!" << std::endl; // seems OK?
		return new SymmetricRotamerSet_();
		//utility_exit_with_message( "Error in RotamerSetFactory, unsupported packing object" ); // get backtrace in gdb
		//exit(1); // add grace
		//return new AminoAcidRotamerSet(); // appease compiler
	}
}

}
}
}
}
