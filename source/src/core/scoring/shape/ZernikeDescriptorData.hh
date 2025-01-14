// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/shape/ZernikeDescriptorData.hh
/// @brief
/// @details
/// @author


#ifndef INCLUDED_core_scoring_shape_ZernikeDescriptorData_hh
#define INCLUDED_core_scoring_shape_ZernikeDescriptorData_hh

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <vector>
#include <basic/datacache/CacheableData.hh>

// Utility headers

namespace core {
namespace scoring {
namespace shape {

class ZernikeDescriptorData;
typedef utility::pointer::shared_ptr< ZernikeDescriptorData > ZernikeDescriptorDataOP;

class ZernikeDescriptorData : public basic::datacache::CacheableData {

public:

	ZernikeDescriptorData() = default;

	ZernikeDescriptorData( std::vector<double> const && invariants );

	~ZernikeDescriptorData() override;

	std::vector<double> const &
	invariants() const {
		return invariants_;
	};

	//this class lives in the PoseCache.... need to provide clone()
	basic::datacache::CacheableDataOP clone() const override {
		return utility::pointer::make_shared< ZernikeDescriptorData >(*this);
	}

//	void show(std::ostream&) const;


private:
	std::vector<double> invariants_;

};


} //shape
} //scoring
} //core

#endif // INCLUDED_core_scoring_shape_ZernikeDescriptorData_HH
