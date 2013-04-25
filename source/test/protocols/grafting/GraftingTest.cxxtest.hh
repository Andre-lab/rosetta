// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/grafting/GraftingTest.cxxtest.hh
/// @brief  tests for container GraftMovers and utility functions
/// @author Jared Adolf-Bryfogle, Brian Weitzner

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>

// Project Headers
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/util.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>

// Protocol Headers
#include <basic/Tracer.hh>


static basic::Tracer TR("protocols.grafting.GraftTest");
class GraftingTest : public CxxTest::TestSuite {
    core::pose::Pose scaffold_pose; //Full PDB
	core::pose::Pose framework_pose; //PDB Missing a cdr.
    core::pose::Pose piece; //CDR to graft.
    core::Size start;
    core::Size end;
    core::Size flex;
    core::Size nter_overhang;
    core::Size cter_overhang;
    core::Size starting_residues;
	core::Size insert_size;
    protocols::grafting::AnchoredGraftMoverOP anchored_grafter;

    
    
public:
	
	void setUp(){
		
		core_init();
		core::import_pose::pose_from_pdb(scaffold_pose, "protocols/grafting/2j88.pdb");
		core::import_pose::pose_from_pdb(framework_pose, "protocols/grafting/2j88_NoL1.pdb");
		core::import_pose::pose_from_pdb(piece, "protocols/grafting/2j88_L1_overhang3_rotated.pdb");
		
		starting_residues = scaffold_pose.total_residue();
		nter_overhang=3;
		cter_overhang=3;
		flex=0;
		start = 23;
		end = 35;
		insert_size = piece.total_residue()-nter_overhang-cter_overhang;
		anchored_grafter = new protocols::grafting::AnchoredGraftMover(start, end);
		TR <<"Setup"<<std::endl;
    }
	
	void tearDown(){
		scaffold_pose.clear();
		framework_pose.clear();
		piece.clear();
	}
    
    void test_utility_functions(){
		
		TR<<"Return region"<<std::endl;
		core::pose::Pose new_region = protocols::grafting::return_region(piece, 1+nter_overhang, piece.total_residue()-cter_overhang);
		TS_ASSERT_EQUALS(insert_size, new_region.total_residue());
		
		TR<<"Replace Region"<<std::endl;
		protocols::grafting::replace_region(scaffold_pose, new_region, 1, start, new_region.total_residue());
		TS_ASSERT_EQUALS(starting_residues, scaffold_pose.total_residue());

		TR<<"Delete Region"<<std::endl;
		core::pose::Pose scaffold_copy = core::pose::Pose(scaffold_pose);
		protocols::grafting::delete_region(scaffold_copy, start+1, end-1);
		core::Size deleted_residues  = starting_residues-scaffold_copy.total_residue();
		TS_ASSERT_EQUALS(deleted_residues, insert_size);
	
		TR<<"Insert Region"<<std::endl;
		core::pose::Pose new_pose = protocols::grafting::insert_pose_into_pose(scaffold_copy, new_region, start, start+1);
		TS_ASSERT_EQUALS(new_pose.total_residue(), starting_residues);
	}
    
    void test_anchoredgraft(){
		
		//test::UTracer UT("protocols/grafting/GraftingAnchoredGraftFunctions.u");
		core::pose::Pose scaffold_copy = core::pose::Pose(scaffold_pose);
		
		TR<<"Testing AnchoredGraft"<<std::endl;
		anchored_grafter->set_piece(piece, nter_overhang, cter_overhang);
		anchored_grafter->superimpose_overhangs_heavy(scaffold_pose, true, true);
		anchored_grafter->set_scaffold_flexibility(flex, flex);
		anchored_grafter->set_insert_flexibility(0, 0);
		anchored_grafter->set_cycles(1);
		anchored_grafter->set_use_double_loop_double_CCD_arms(true);
		anchored_grafter->apply(scaffold_pose);
        
		TS_ASSERT_EQUALS(scaffold_pose.total_residue(), starting_residues);
		core::Real rms = core::scoring::CA_rmsd(scaffold_copy, scaffold_pose);
		core::Real rms_assert = 0.0;
		TS_ASSERT_EQUALS(rms, rms_assert);
    }
    
};



