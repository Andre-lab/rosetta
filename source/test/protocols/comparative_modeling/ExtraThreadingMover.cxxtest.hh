// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/comparative_modeling/ExtraThreadingMover.cxxtest.hh
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <core/pose/Pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/conformation/Residue.hh>

#include <protocols/comparative_modeling/ExtraThreadingMover.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

//Auto Headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <utility/vector1.hh>
#include <string>


class ExtraThreadingMover_Tests : public CxxTest::TestSuite {

public:
	ExtraThreadingMover_Tests() {}

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_basic_threading() {
		using namespace core::sequence;
		using namespace protocols::comparative_modeling;
		using core::Size;
		using core::Real;
		using utility::vector1;
		using core::pose::Pose;
		using core::import_pose::pose_from_file;
		using core::pose::make_pose_from_sequence;

		SequenceOP query( new Sequence( "MKNGEQNGPTTCTNCFTQTTPLWRRNPEGQPLCNACGLFLKLHGVVRPLSLKTDVIKKRNRNSANS", "4gat_prot" ) );
		SequenceOP templ( new Sequence( "MKNGEQNGPTTCTNCFTQTTPLWRRNPEGQPLCNACGLFLKLHGVVRPLSLKTDVIKKRNRNSANS", "4gat_all" ) );

		SequenceAlignment align;
		align.add_sequence(templ);
		align.add_sequence(query);

		Pose query_pose, template_pose;
		core::import_pose::pose_from_file( query_pose, "protocols/comparative_modeling/4gat_protein.pdb" , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( template_pose, "protocols/comparative_modeling/4gat_all.pdb" , core::import_pose::PDB_file);

		utility::vector1< Size > residues_to_steal;

		residues_to_steal.push_back( 67 ); // Zinc
		for ( Size ii = 68; ii <= 93; ++ii ) {
			residues_to_steal.push_back(ii); // DNA
		}

		// test basic threading
		ExtraThreadingMover mover( align, template_pose, residues_to_steal );
		mover.apply(query_pose);

		for ( Size ii = 1; ii <= query_pose.size(); ++ii ) {
			TS_ASSERT_EQUALS( query_pose.residue(ii).natoms(), template_pose.residue(ii).natoms() );
			TS_ASSERT_EQUALS( query_pose.residue(ii).nheavyatoms(), template_pose.residue(ii).nheavyatoms() );
			TS_ASSERT_EQUALS( query_pose.residue_type(ii).name(), template_pose.residue_type(ii).name() );
		}

		// make sure that chains are set correctly
		TS_ASSERT( query_pose.residue(67).chain() == 2 );
		for ( Size ii = 68; ii <= 80; ++ii ) TS_ASSERT(query_pose.residue(ii).chain() == 3);
		for ( Size ii = 81; ii <= 93; ++ii ) TS_ASSERT(query_pose.residue(ii).chain() == 4);

	} // test_basic_threading

}; // class ExtraThreadingMover_Tests
