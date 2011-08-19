// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSet.fwd.hh
/// @brief  Residue set class forward declarations
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerCouplings_hh
#define INCLUDED_core_pack_rotamer_set_RotamerCouplings_hh

#include <core/pack/rotamer_set/RotamerCouplings.fwd.hh>

#include <core/conformation/ResidueMatcher.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

// AUTO-REMOVED #include <utility>

//Auto Headers
#include <utility/vector1_bool.hh>

// utility headers

namespace core {
namespace pack {
namespace rotamer_set {

class RotamerCouplings : public utility::pointer::ReferenceCount {
public:
	typedef conformation::ResidueMatcherCOP ResidueMatcherCOP;

public:
	void
	resize( Size const size_in )
	{
		couplings_.resize( size_in );
	}

	std::pair< int, ResidueMatcherCOP > const &
	operator[]( Size const index ) const
	{
		return couplings_[ index ];
	}

	std::pair< int, ResidueMatcherCOP > &
	operator[]( Size const index )
	{
		return couplings_[ index ];
	}

private:

	utility::vector1< std::pair< int, ResidueMatcherCOP > > couplings_;
};


} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif //
