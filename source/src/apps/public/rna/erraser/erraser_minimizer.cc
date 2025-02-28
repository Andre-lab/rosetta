// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <protocols/rna/movers/ErraserMinimizerMover.hh>

//////////////////////////////////////////////////
///////////////////////////////////////////////////

#include <core/types.hh>

//////////////////////////////////////////////////////////

// C++ headers

#include <protocols/jd2/JobDistributor.hh>

#include <protocols/viewer/viewers.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/options/keys/OptionKeys.hh> // AUTO IWYU For OptionKeys,

//////////////////////////////////////////////////////////////////
//
// Rhiju-style comment of long-term objectives and work to do.
//
// 1. check_in_bonded_list and check_in_bond_angle_list aren't
// ideal but it's not a real performance concern. But you could
// easily imagine storing them sorted, thus permitting much better
// searches. This is obviously not the ERRASER bottleneck, but
// std::set provides much of this functionality, just needing a
// sorting function to do the rest. (Also uniqueness is then
// enforced.
//
//        -- AMW, 2016
//
//////////////////////////////////////////////////////////////////

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;
using namespace protocols::rna::movers;

///////////////////////////////////////////////////////////////
void*
my_main ( void* ) {

	ErraserMinimizerMoverOP emm( new ErraserMinimizerMover );
	emm->initialize_from_options( option );

	protocols::jd2::JobDistributor::get_instance()->go( emm );

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main ( int argc, char * argv [] ) {
	try {
		devel::init( argc, argv );
		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}
