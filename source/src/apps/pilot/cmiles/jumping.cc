// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/cmiles/jumping.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/cmiles.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/jumping/util.hh>

int main(int argc, char* argv[]) {
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		devel::init(argc, argv);

		core::pose::PoseCOP pose = core::import_pose::pose_from_file(option[OptionKeys::in::file::native](), core::import_pose::PDB_file);

		core::Size i = option[OptionKeys::cmiles::jumping::resi]();
		core::Size j = option[OptionKeys::cmiles::jumping::resj]();

		core::Size orientation, pleating;
		protocols::jumping::get_pleating(*pose, i, j, orientation, pleating);

		std::cout << "orientation: " << orientation << std::endl;
		std::cout << "pleating: " << pleating << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
