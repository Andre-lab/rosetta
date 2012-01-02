// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// libRosetta headers
// AUTO-REMOVED #include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>


// C++ headers
#include <iostream>
#include <string>

#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/relax_main.hh>
#include <utility/vector1.hh>



using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;
using namespace protocols;

using utility::vector1;



///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	using namespace protocols::moves;
	using namespace scoring;
	using namespace basic::options;


	relax::ClassicRelax::register_options();
	jd2::register_options();
	option.add_relevant( OptionKeys::in::file::fullatom );
	option.add_relevant( OptionKeys::relax::fast );
	devel::init(argc, argv);

	return relax::Relax_main( false );
}
