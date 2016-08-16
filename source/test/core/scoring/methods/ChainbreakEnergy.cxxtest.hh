// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/methods/ChainbreakEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::methods::ChainbreakEnergy
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>

// Utility headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class ChainbreakEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-include_sugars" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	/// This can be used to ensure norm matches norm-numeric
	void test_chainbreak_deriv_check()
	{
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::scoring;

		Pose pose = pdb1ubq5to13_pose();
		FoldTree ft;
		ft.add_edge( 1, 3, -1 );
		ft.add_edge( 3, 7,  1 );
		ft.add_edge( 7, 9, -1 );
		ft.add_edge( 7, 6, -1 );
		ft.add_edge( 3, 5, -1 );

		pose.fold_tree( ft );
		pose.set_phi( 4, 180 );
		pose.set_psi( 4, 180 );
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, 5 );
		core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, 6 );


		//pose.dump_pdb( "test_ft1.pdb" );

		ScoreFunction sfxn;
		sfxn.set_weight( chainbreak, 0.5 );
		TS_ASSERT_DELTA( sfxn( pose ), 126.855, 1e-3 );

		//std::cout << "chainbreak score: " << sfxn( pose ) << std::endl;
		MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

	void test_chainbreak_deriv_check_with_sugars()
	{
		using namespace core::chemical;
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::scoring;

		Pose pose;
		make_pose_from_saccharide_sequence( pose, "Glcp-Glcp-Glcp-Glcp-Glcp-Glcp", "fa_standard" );
		for ( core::uint i = 1; i < 6; ++i ) {
			// Make extended.
			pose.set_phi( i, 100.0 );
			pose.set_psi( i, 220.0 );
		}
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, 3 );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, 4 );

		ScoreFunction sfxn;
		sfxn.set_weight( chainbreak, 0.5 );
		TS_ASSERT_DELTA( sfxn( pose ), 0.0, 1e-3 );

		MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}
};


