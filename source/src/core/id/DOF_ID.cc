// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/id/DOF_ID.cc
/// @brief  Kinematics DOF identifier class
/// @author Phil Bradley


// Unit headers
#include <core/id/DOF_ID.hh>

// C++ headers
// AUTO-REMOVED #include <ostream>

#include <iostream>



namespace core {
namespace id {


/// @brief stream << DOF_ID
std::ostream &
operator <<( std::ostream & os, DOF_ID const & a )
{
	os << " atom_id= " << a.atom_id() << " type= " << a.type() << ' ';
	return os;
}


/// @brief Globals
DOF_ID const BOGUS_DOF_ID( BOGUS_ATOM_ID, PHI );


} // namespace id
} // namespace core
