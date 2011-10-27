// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_residue_io_hh
#define INCLUDED_core_chemical_residue_io_hh


// Unit headers

// Project headers
#include <core/chemical/ResidueType.fwd.hh>

#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
//#include <core/chemical/CSDAtomTypeSet.fwd.hh>
// Auto-header: duplicate removed #include <core/chemical/ResidueTypeSet.fwd.hh>
#include <string>

// Utility headers

// C++ headers

namespace core {
namespace chemical {

/// @brief virtual constructor for ResidueType objects
ResidueTypeOP
read_topology_file(
	std::string const & filename,
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
//	chemical::CSDAtomTypeSetCAP csd_atom_types kwk commenting out csd_atom_types until I have a chance to fully implement them.
	chemical::ResidueTypeSetCAP rsd_type_set
);

/// @brief writes a .params file from a given ResidueType object
void
write_topology_file(
	ResidueType const & rsd
);




} // chemical
} // core



#endif
