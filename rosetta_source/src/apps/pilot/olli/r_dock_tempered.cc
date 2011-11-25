// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @ r_dock_tempered.cc
/// @ author Oliver Lange


// Rosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/docking/TemperedDocking.hh>

#include <devel/init.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// option key includes
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/option.hh>


int
main( int argc, char * argv [] )
{
  using namespace basic::options;
  using namespace protocols::docking;
  using namespace protocols::jd2;

	protocols::jd2::register_options();
	TemperedDocking::register_options();

  // initialize core
  devel::init(argc, argv);
  //	core::init_random_generators(3,numeric::random::_RND_TestRun_, "mt19937"); //JQX from Sergery

	JobDistributor::get_instance()->go( new TemperedDocking() );
}

