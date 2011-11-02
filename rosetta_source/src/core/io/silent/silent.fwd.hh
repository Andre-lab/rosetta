// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/silent_in.FWD_HH.hh
///
/// @brief silent input file reader for mini.
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_silent_fwd_hh
#define INCLUDED_core_io_silent_silent_fwd_hh

// mini headers
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <map>
// AUTO-REMOVED #include <string>

#include <core/types.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace io {
namespace silent {

class SharedSilentData;
class SimpleSequenceData;
class EnergyNames;

// abstract base classes
class SilentStruct;
class SilentFileData;

// derived classes
template <class T> class ProteinSilentStruct_Template;
typedef ProteinSilentStruct_Template< float > ProteinSilentStruct_SinglePrec;
typedef ProteinSilentStruct_Template< core::Real > ProteinSilentStruct;
class BinaryProteinSilentStruct;
class RNA_SilentStruct;

// owning pointers
typedef utility::pointer::owning_ptr< SharedSilentData > SharedSilentDataOP;
typedef utility::pointer::owning_ptr< EnergyNames > EnergyNamesOP;
typedef utility::pointer::owning_ptr< SimpleSequenceData > SimpleSequenceDataOP;

typedef utility::pointer::owning_ptr< SilentStruct > SilentStructOP;
typedef utility::pointer::owning_ptr< ProteinSilentStruct > ProteinSilentStructOP;
typedef utility::pointer::owning_ptr< BinaryProteinSilentStruct > BinaryProteinSilentStructOP;
typedef utility::pointer::owning_ptr< RNA_SilentStruct > RNA_SilentStructOP;

	// data types
typedef std::map< std::string, SilentStructOP > Structure_Map;
typedef utility::vector1< SilentStructOP > SilentStructOPs;

} // namespace silent
} // namespace io
} // namespace core

#endif
