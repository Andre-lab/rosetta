// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_multithreaded_pdb_load.cc
/// @brief A pilot app to test the creation of several poses at once (without first initializing anything by having a single thread initialize the pose).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/types.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>

// protocol headers

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/init/init.hh>
#include <protocols/relax/FastRelax.hh>

#include <boost/bind.hpp>

#ifdef MULTI_THREADED
#include <thread>
#include <mutex>
#include <atomic>
#endif

static basic::Tracer TR("test_multithreaded_pdb_load");

#define NUM_THREADS 30

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );
}

void
thread_fxn( core::Size const thread_index ) {
	TR << "Launching thread " << thread_index << std::endl;

	core::init::init_random_generators( 111111, "standard" );

	core::pose::PoseOP pose( core::import_pose::pose_from_file( basic::options::option[basic::options::OptionKeys::in::file::s]()[1], false, core::import_pose::PDB_file ) );

	TR << "Thread " << thread_index << " reports that the pose has " << pose->total_residue() << " residues." << std::endl;

	TR << "Thread " << thread_index << " attempting FastRelax." << std::endl;

	core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
	protocols::relax::FastRelax frlx( sfxn, 3 );
	frlx.apply(*pose);

	TR << "Thread " << thread_index << " reports post-relax energy of " << pose->energies().total_energy() << "." << std::endl;
	TR << "Thread " << thread_index << " terminated." << std::endl;
}


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );
		register_options();

		if ( ( ! option [ in::file::l ].user() ) && ( ! option [ in::file::s ].user() ) ) {
			utility_exit_with_message("Please specify either -s or -l to specify the input PDB.");
		}

		utility::vector1< std::thread > threads;
		for ( core::Size i(1); i<=NUM_THREADS; ++i ) {
			threads.push_back( std::thread( boost::bind( &thread_fxn, i ) ) );
		}
		for ( core::Size i(1); i<=NUM_THREADS; ++i ) {
			threads[i].join();
		}

		TR << "Execution completed." << std::endl;


	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
