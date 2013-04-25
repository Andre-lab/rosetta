// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file NativeResidueReader.hh
///
/// @brief

/// @author Tim jacobs

#ifndef NATIVERESIDUEREADER_HH_
#define NATIVERESIDUEREADER_HH_

//Core
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>

//Utility
#include <utility/file/FileName.hh>

//C++
#include <string>
#include <map>
#include <stdio.h>

class NativeResidueReader{

public:

  std::map<core::Size, utility::vector1<core::conformation::ResidueOP> > generateResiduesFromFile(utility::file::FileName file);

private:

};

#endif /* NATIVERESIDUEREADER_HH_ */
