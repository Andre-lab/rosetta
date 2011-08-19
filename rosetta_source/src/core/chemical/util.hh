// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/util.hh
///
/// @brief Utilities for modifying and utilizing Residues and other core::chemical classes.
/// @author

#ifndef INCLUDED_core_chemical_util_hh
#define INCLUDED_core_chemical_util_hh


// Unit headers
#include <core/chemical/ResidueTypeSet.fwd.hh>


namespace core {
namespace chemical {

core::chemical::ResidueTypeSetCAP rsd_set_from_cmd_line();

} // namespace chemical
} // namespace core

#endif
