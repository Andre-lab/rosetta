// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ProtocolUtil.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_MainChainTorsionSet_hh
#define INCLUDED_protocols_swa_MainChainTorsionSet_hh

#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>
#include <utility/vector1.fwd.hh>
// AUTO-REMOVED #include <string>

namespace protocols {
namespace swa {

//////////////////////////////////
class MainChainTorsionSet: public utility::pointer::ReferenceCount {


public:

	MainChainTorsionSet( core::Real const & phi, core::Real const & psi, core::Real const & omega );

	MainChainTorsionSet( core::Real const & phi, core::Real const & psi );

	~MainChainTorsionSet();

	MainChainTorsionSet
	operator=( MainChainTorsionSet const & src );

	core::Real phi() const;
	core::Real psi() const;
	core::Real omega() const;

private:
	core::Real phi_;
	core::Real psi_;
	core::Real omega_;
};

	typedef utility::vector1< MainChainTorsionSet > MainChainTorsionSetList ;


} //swa
} // protocols

#endif
