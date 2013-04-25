// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.cxxtest.hh
/// @brief  unit testing for core/scoring/ScoreFunction.*
/// @author Rocco Moretti (rmoretti@u.washington.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/UTracer.hh>
#include <test/core/init_util.hh>

// Unit headers
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
#include <core/types.hh>

// Package headers
#include <basic/database/open.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
#include <utility/vector1.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

//Auto Headers


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.ScoreFunction.cxxtest");

// using declarations
using namespace core;
using namespace scoring;

///////////////////////////////////////////////////////////////////////////
/// @name ScoreFunctionUtilityTest
/// @brief: unified tests for associated score function utilities
///////////////////////////////////////////////////////////////////////////
class ScoreFunctionUtilityTest : public CxxTest::TestSuite {

public:

	void setUp() {
	    core_init();
	}

	void tearDown() {}

	void test_find_weights_file() {
		// Local weights
		TS_ASSERT_EQUALS( find_weights_file("core/scoring/test",".wts"), "core/scoring/test.wts" );
		TS_ASSERT_EQUALS( find_weights_file("core/scoring/test.wts",".wts"), "core/scoring/test.wts" );
		// Local patch
		TS_ASSERT_EQUALS( find_weights_file("core/scoring/test",".wts_patch"), "core/scoring/test.wts_patch" );
		TS_ASSERT_EQUALS( find_weights_file("core/scoring/test.wts_patch",".wts_patch"), "core/scoring/test.wts_patch" );
		// Database weights
		TS_ASSERT_EQUALS( find_weights_file("standard",".wts"), basic::database::full_name( "scoring/weights/standard.wts" ) );
		TS_ASSERT_EQUALS( find_weights_file("standard.wts",".wts"), basic::database::full_name( "scoring/weights/standard.wts" ) );
		// Database patch
		TS_ASSERT_EQUALS( find_weights_file("score12",".wts_patch"), basic::database::full_name( "scoring/weights/score12.wts_patch" ) );
		TS_ASSERT_EQUALS( find_weights_file("score12.wts_patch",".wts_patch"), basic::database::full_name( "scoring/weights/score12.wts_patch" ) );
	}



};

///////////////////////////////////////////////////////////////////////////
/// @name ScoreFunctionTest
/// @brief: unified tests for ScoreFunction object
///////////////////////////////////////////////////////////////////////////
class ScoreFunctionTest : public CxxTest::TestSuite {

public:

	void setUp() {
	    core_init();
	}

	void tearDown() {}

	void test_add_weights_from_file() {

		// Test local loading
		ScoreFunction sfxn;
		sfxn.add_weights_from_file("core/scoring/test");
		TS_ASSERT_DELTA( sfxn.get_weight(fa_rep) , 3.1416, 0.0001 );

		// Test database loading
		ScoreFunction sfxn2;
		sfxn2.add_weights_from_file("standard");
		TS_ASSERT_DELTA( sfxn2.get_weight(fa_atr) , 0.8, 0.0001 );
	}

	void test_apply_patch_from_file() {

		// Test local loading
		ScoreFunction sfxn;
		sfxn.add_weights_from_file("core/scoring/test");
		sfxn.apply_patch_from_file("core/scoring/test");
		TS_ASSERT_DELTA( sfxn.get_weight(fa_rep) , 3.1416, 0.0001 );
		TS_ASSERT_DELTA( sfxn.get_weight(fa_atr) , 2.7183, 0.0001 );

		// Test database loading
		ScoreFunction sfxn2;

		sfxn2.add_weights_from_file("standard");
		sfxn2.apply_patch_from_file("score12");
		TS_ASSERT_DELTA( sfxn2.get_weight(fa_atr) , 0.8, 0.0001 );
		TS_ASSERT_DELTA( sfxn2.get_weight(p_aa_pp) , 0.32, 0.0001 );
	}


	void test_get_sub_score_exclude_res() {

		ScoreFunction scfxn;
		core::pose::Pose pose = create_trpcage_ideal_pose();
		Real sc(scfxn(pose));

		utility::vector1< Size > exclude_list;
		Real sc_exc1(scfxn.get_sub_score_exclude_res(pose, exclude_list));
		TS_ASSERT_DELTA(sc_exc1, sc, .0000001);

		utility::vector1< bool > residue_mask(pose.total_residue(), true);
		Real sc_exc2(scfxn.get_sub_score(pose, residue_mask));
		TS_ASSERT_DELTA(sc_exc2, sc, .0000001);
	}

	///@detail the long range terms maintain their own neighbor maps, so
	///require iterating over the energy methods then over the residue
	///pairs rather then the other way around.
	void test_get_sub_score_exclude_res_long_range_terms() {

		ScoreFunction scfxn;

		//CartesianBondedEnergy is a ContextIndependentLRTwoBodyEnergy
		scfxn.set_weight(core::scoring::cart_bonded, 1);

		core::pose::Pose pose = create_trpcage_ideal_pose();
		Real sc(scfxn(pose));

		utility::vector1< Size > exclude_list;
		Real sc_exc1(scfxn.get_sub_score_exclude_res(pose, exclude_list));
		TS_ASSERT_DELTA(sc_exc1, sc, .0000001);

		utility::vector1< bool > residue_mask(pose.total_residue(), true);
		Real sc_exc2(scfxn.get_sub_score(pose, residue_mask));
		TS_ASSERT_DELTA(sc_exc2, sc, .0000001);
	}

};
