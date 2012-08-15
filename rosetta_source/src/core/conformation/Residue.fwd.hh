// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/Residue.fwd.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_Residue_fwd_hh
#define INCLUDED_core_conformation_Residue_fwd_hh

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

// C++ headers

namespace core {
namespace conformation {

class Residue;

typedef  utility::pointer::access_ptr< Residue >  ResidueAP;
typedef  utility::pointer::access_ptr< Residue const >  ResidueCAP;
typedef  utility::pointer::owning_ptr< Residue >  ResidueOP;
typedef  utility::pointer::owning_ptr< Residue const >  ResidueCOP;

typedef  utility::vector1< ResidueOP >  ResidueOPs;
typedef  utility::vector1< ResidueCOP >  ResidueCOPs;
typedef  utility::vector1< ResidueCAP >  ResidueCAPs;

#ifdef USEBOOSTSERIALIZE
// why is this here? you don't want to know
template<class Archive> void save_construct_data( Archive & ar, const Residue * t, const unsigned int file_version);
template<class Archive> void load_construct_data( Archive & ar, Residue * t, const unsigned int file_version);
#endif

} // namespace conformation
} // namespace core

#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/map.hpp>
#endif

#endif // INCLUDED_core_conformation_Residue_FWD_HH
