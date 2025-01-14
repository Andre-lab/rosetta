// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/shape/ZernikeDescriptorData.cc
/// @brief
/// @author


//Unit headers
#include <core/scoring/shape/ZernikeDescriptorData.hh>

//Package headers
#include <core/pose/Pose.hh>

//C++ headers
//#include <iostream>

/// Utility headers
//#include <utility/exit.hh>
//#include <utility/excn/Exceptions.hh>
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>

// term specific headers


static basic::Tracer TR( "core.scoring.shape.ZernikeDescriptorData"  );

using namespace core::scoring;

namespace core {
namespace scoring {
namespace shape {

ZernikeDescriptorData::~ZernikeDescriptorData() = default;

ZernikeDescriptorData::ZernikeDescriptorData( std::vector<double> const && invariants ) {
    invariants_ = invariants;
}

} // shape
} // scoring
} // core
