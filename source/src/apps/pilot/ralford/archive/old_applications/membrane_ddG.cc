// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/ralford/membrane_ddG.cc
///
/// @brief   C++ Version of Computing ddGs in Membranes
/// @details    waiting on Pyrosetta
///
/// @author  Rebecca Alford (rflford12@gmail.com)
/// @note       Last Modified: 6/24/14

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/ddG/MembraneDDGMover.hh>
#include <protocols/membrane/ddG/Mutation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <core/chemical/AA.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR( "apps.pilot.ralford.membrane_ddG" );

using protocols::membrane::ddG::MutationOP;

/// @brief ompLA Task
utility::vector1< MutationOP >
ompLA_task() {

	using namespace protocols::membrane::ddG;

	// Create a new list of mutations
	utility::vector1< MutationOP > mutations;
	mutations.resize(20);

	// Append a list
	mutations[1] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_ala );
	mutations[2] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_cys );
	mutations[3] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_asp );
	mutations[4] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_glu );
	mutations[5] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_phe );
	mutations[6] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_gly );
	mutations[7] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_his );
	mutations[8] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_ile );
	mutations[9] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_lys );
	mutations[10] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_leu );
	mutations[11] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_met );
	mutations[12] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_asn );
	mutations[13] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_pro );
	mutations[14] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_gln );
	mutations[15] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_arg );
	mutations[16] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_ser );
	mutations[17] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_thr );
	mutations[18] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_val );
	mutations[19] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_trp );
	mutations[20] = utility::pointer::make_shared< Mutation >( 181, core::chemical::aa_tyr );

	return mutations;
}

/// @brief Lukas Tamm - Aromatic residues at the interface regions task
/// @remarks protein name = ompA
utility::vector1< MutationOP >
ompA_task() {

	using namespace protocols::membrane::ddG;

	// Create a new list of mutations
	utility::vector1< MutationOP > mutations;
	mutations.resize(12);

	// Append a list
	mutations[1] = utility::pointer::make_shared< Mutation >( 7, core::chemical::aa_ala );
	mutations[2] = utility::pointer::make_shared< Mutation >( 15, core::chemical::aa_ala );
	mutations[3] = utility::pointer::make_shared< Mutation >( 57, core::chemical::aa_ala );
	mutations[4] = utility::pointer::make_shared< Mutation >( 143, core::chemical::aa_ala );
	mutations[5] = utility::pointer::make_shared< Mutation >( 43, core::chemical::aa_ala );
	mutations[6] = utility::pointer::make_shared< Mutation >( 55, core::chemical::aa_ala );
	mutations[7] = utility::pointer::make_shared< Mutation >( 129, core::chemical::aa_ala );
	mutations[8] = utility::pointer::make_shared< Mutation >( 141, core::chemical::aa_ala );
	mutations[9] = utility::pointer::make_shared< Mutation >( 168, core::chemical::aa_ala );
	mutations[10] = utility::pointer::make_shared< Mutation >( 51, core::chemical::aa_ala );
	mutations[11] = utility::pointer::make_shared< Mutation >( 123, core::chemical::aa_ala );
	mutations[12] = utility::pointer::make_shared< Mutation >( 170, core::chemical::aa_ala );

	return mutations;
}

/// @brief Bowie - Bacteriorhodopsin B Helix Mutations
/// @remarks PDBID = 1PY6
utility::vector1< MutationOP >
bacteriorhodpsin_task() {

	using namespace protocols::membrane::ddG;

	// Create a new list of mutations
	utility::vector1< MutationOP > mutations;
	mutations.resize(12);

	// Append a list
	mutations[1] = utility::pointer::make_shared< Mutation >( 31, core::chemical::aa_ala );
	mutations[2] = utility::pointer::make_shared< Mutation >( 32, core::chemical::aa_ala );
	mutations[3] = utility::pointer::make_shared< Mutation >( 33, core::chemical::aa_ala );
	mutations[4] = utility::pointer::make_shared< Mutation >( 34, core::chemical::aa_ala );
	mutations[5] = utility::pointer::make_shared< Mutation >( 36, core::chemical::aa_ala );
	mutations[6] = utility::pointer::make_shared< Mutation >( 37, core::chemical::aa_ala );
	mutations[7] = utility::pointer::make_shared< Mutation >( 38, core::chemical::aa_ala );
	mutations[8] = utility::pointer::make_shared< Mutation >( 39, core::chemical::aa_ala );
	mutations[9] = utility::pointer::make_shared< Mutation >( 41, core::chemical::aa_ala );
	mutations[10] = utility::pointer::make_shared< Mutation >( 42, core::chemical::aa_ala );
	mutations[11] = utility::pointer::make_shared< Mutation >( 45, core::chemical::aa_ala );
	mutations[12] = utility::pointer::make_shared< Mutation >( 44, core::chemical::aa_ala );
	mutations[10] = utility::pointer::make_shared< Mutation >( 45, core::chemical::aa_ala );
	mutations[11] = utility::pointer::make_shared< Mutation >( 46, core::chemical::aa_ala );
	mutations[12] = utility::pointer::make_shared< Mutation >( 48, core::chemical::aa_ala );
	mutations[11] = utility::pointer::make_shared< Mutation >( 50, core::chemical::aa_ala );
	mutations[12] = utility::pointer::make_shared< Mutation >( 51, core::chemical::aa_ala );
	mutations[10] = utility::pointer::make_shared< Mutation >( 52, core::chemical::aa_ala );
	mutations[11] = utility::pointer::make_shared< Mutation >( 53, core::chemical::aa_ala );
	mutations[12] = utility::pointer::make_shared< Mutation >( 54, core::chemical::aa_ala );
	mutations[12] = utility::pointer::make_shared< Mutation >( 55, core::chemical::aa_ala );
	mutations[10] = utility::pointer::make_shared< Mutation >( 56, core::chemical::aa_ala );
	mutations[11] = utility::pointer::make_shared< Mutation >( 57, core::chemical::aa_ala );
	mutations[12] = utility::pointer::make_shared< Mutation >( 58, core::chemical::aa_ala );

	return mutations;
}

/// @brief Main method
int
main( int argc, char * argv [] )
{

	using namespace core::scoring;
	using namespace protocols::membrane::ddG;
	using namespace protocols::jd2;

	try {

		// Devel init factories
		devel::init(argc, argv);

		// Register JD2 options
		protocols::jd2::register_options();

		// Create a new energy function
		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "fa_menv_pHmode_2014" );

		// Load in mutants for the appropriate task
		utility::vector1< MutationOP > ompA_mutations = ompA_task();
		utility::vector1< MutationOP > ompLA_mutations = ompLA_task();
		utility::vector1< MutationOP > bacteriorhodopsin_mutations = bacteriorhodpsin_task();

		// Make ddG calculation
		MembraneDDGMoverOP ddG( new MembraneDDGMover( 0.0001, sfxn, ompLA_mutations ) );
		JobDistributor::get_instance()->go( ddG );

		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

