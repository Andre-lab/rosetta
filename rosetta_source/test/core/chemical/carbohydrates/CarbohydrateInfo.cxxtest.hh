// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	CarbohydrateInfo.cxxtest.hh
/// @brief   Test suite for CarbohydrateInfo
/// @author  labonte

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <sstream>


class CarbohydrateInfoTests : public CxxTest::TestSuite {
public:
	// Standard methods ////////////////////////////////////////////////////////
	// Initialization
	void
	setUp()
	{
		using namespace core::import_pose;
		using namespace basic::options;

		core_init();
		option[OptionKeys::in::include_sugars](true);
		pose_from_pdb(maltotriose_, "core/chemical/carbohydrates/maltotriose.pdb");
		pose_from_pdb(isomaltose_,  "core/chemical/carbohydrates/isomaltose.pdb");

		// TODO: Test several other cases besides glucose.
	}

	// Destruction
	void
	tearDown()
	{}


	// Tests ///////////////////////////////////////////////////////////////////
	// Confirm that carbohydrate properties are loaded correctly.
	void test_CarbohydrateInfo_show()
	{
		using namespace std;
		using namespace core::conformation;

		TS_TRACE("Testing show() method of CarbohydrateInfo.");

		stringstream output1(stringstream::out);
		stringstream output2(stringstream::out);
		stringstream output3(stringstream::out);
		stringstream expected_out1(stringstream::out);
		stringstream expected_out2(stringstream::out);
		stringstream expected_out3(stringstream::out);

		expected_out1 << "Carbohydrate Properties for this Residue:" << endl;
		expected_out1 << " Basic Name: glucose" << endl;
		expected_out1 << " IUPAC Name: ->4)-D-glucopyranose" << endl;
		expected_out1 << " Classification: aldohexose" << endl;
		expected_out1 << " Stereochemistry: D" << endl;
		expected_out1 << " Ring Form: pyranose" << endl;
		expected_out1 << " Anomeric Form: alpha" << endl;
		expected_out1 << " Modifications: " << endl;
		expected_out1 << "  none" << endl;
		expected_out1 << "" << endl;
		expected_out1 << " Polymeric Information:" << endl;
		expected_out1 << "  Main chain connection: (_->4)" << endl;
		expected_out1 << "  Branch connections: branches not yet implemented" << endl;

		expected_out2 << "Carbohydrate Properties for this Residue:" << endl;
		expected_out2 << " Basic Name: glucose" << endl;
		expected_out2 << " IUPAC Name: ->4)-alpha-D-glucopyranosyl" << endl;
		expected_out2 << " Classification: aldohexose" << endl;
		expected_out2 << " Stereochemistry: D" << endl;
		expected_out2 << " Ring Form: pyranose" << endl;
		expected_out2 << " Anomeric Form: alpha" << endl;
		expected_out2 << " Modifications: " << endl;
		expected_out2 << "  none" << endl;
		expected_out2 << "" << endl;
		expected_out2 << " Polymeric Information:" << endl;
		expected_out2 << "  Main chain connection: (_->4)" << endl;
		expected_out2 << "  Branch connections: branches not yet implemented" << endl;

		expected_out3 << "Carbohydrate Properties for this Residue:" << endl;
		expected_out3 << " Basic Name: glucose" << endl;
		expected_out3 << " IUPAC Name: alpha-D-glucopyranosyl" << endl;
		expected_out3 << " Classification: aldohexose" << endl;
		expected_out3 << " Stereochemistry: D" << endl;
		expected_out3 << " Ring Form: pyranose" << endl;
		expected_out3 << " Anomeric Form: alpha" << endl;
		expected_out3 << " Modifications: " << endl;
		expected_out3 << "  none" << endl;
		expected_out3 << "" << endl;
		expected_out3 << " Polymeric Information:" << endl;
		expected_out3 << "  Main chain connection: N/A" << endl;
		expected_out3 << "  Branch connections: branches not yet implemented" << endl;

		Residue res1 = maltotriose_.residue(1);
		res1.carbohydrate_info()->show(output1);
		TS_ASSERT_EQUALS(output1.str(), expected_out1.str());

		Residue res2 = maltotriose_.residue(2);
		res2.carbohydrate_info()->show(output2);
		TS_ASSERT_EQUALS(output2.str(), expected_out2.str());

		Residue res3 = maltotriose_.residue(3);
		res3.carbohydrate_info()->show(output3);
		TS_ASSERT_EQUALS(output3.str(), expected_out3.str());
	}

	// Confirm that CarbohydrateInfo.short_name_ is assigned correctly.
	void test_Pose_chain_sequence_w_polysaccharide()
	{
		TS_TRACE("Testing chain_sequence() method of Pose with polysaccharide chains.");
		TS_ASSERT_EQUALS(maltotriose_.chain_sequence(1), "alpha-D-Glcp-(1->4)-alpha-D-Glcp-(1->4)-D-Glcp");
		TS_ASSERT_EQUALS(isomaltose_.chain_sequence(1), "alpha-D-Glcp-(1->6)-D-Glcp");
	}

	// Confirm that backbone torsion angles are assigned correctly.
	void test_Pose_phi_psi_omega_w_polysaccharide()
	{
		TS_TRACE("Testing phi(), psi(), and omega() methods of Pose with polysaccharide chains.");
		TS_TRACE(isomaltose_.phi(2));
		TS_TRACE(isomaltose_.psi(2));
		TS_TRACE(isomaltose_.omega(2));
		//TS_ASSERT_DELTA(isomaltose_.phi(2), 44.3268, 0.001);
		//TS_ASSERT_DELTA(isomaltose_.psi(2), 170.869, 0.001);
		//TS_ASSERT_DELTA(isomaltose_.omega(2), 49.383, 0.001);
	}

	// Confirm that branches are handled properly.
	void test_CarbohydrateInfo_branch_point()
	{
		// TODO
		TS_TRACE("Testing branch_point() method of CarbohydrateInfo.");
	}

	// Confirm that nu angles are mapped to the correct CHI identifiers.
	void test_CarbohydrateInfo_nu_id()
	{
		using namespace core;
		using namespace conformation;

		TS_TRACE("Testing nu_id() method of CarbohydrateInfo.");
		for (Size resnum = 1; resnum <= 3; ++resnum) {
			for (Size i = 1; i <= 4; ++i) {
				// For glucose, the 4 CHI ids should be 7, 8, 9, and 10.
				Residue res = maltotriose_.residue(resnum);
				TS_ASSERT_EQUALS(res.carbohydrate_info()->nu_id(i).second, 6 + i);
			}
		}
	}

private:
	// Private data ////////////////////////////////////////////////////////////
	core::pose::Pose maltotriose_;  // a (1alpha->4) trisaccharide of D-glucose
	core::pose::Pose isomaltose_;  // a (1alpha->6) disaccharide of D-glucose

};  // class CarbohydrateInfoTests
