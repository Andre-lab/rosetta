// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   debug_labontes_current_work.cc
/// @brief  This is simply a generic pilot app for testing things.
/// @author Labonte


// Package headers
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/io/carbohydrates/pose_io.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/annotated_sequence.hh>
//#include <core/pose/PDBInfo.hh>
//#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/kinematics/FoldTree.hh>
//#include <core/id/TorsionID.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>

#include <protocols/minimization_packing/MinMover.fwd.hh>
//#include <protocols/simple_moves/BackboneMover.hh>
//#include <protocols/minimization_packing/PackRotamersMover.hh>

// Utility headers

// Numeric headers

// C++ headers
#include <iostream>

#include <utility/excn/Exceptions.hh> // AUTO IWYU For Exception

//#include <algorithm>


int
main( int argc, char *argv[] )
{
	try {
		using namespace std;
		using namespace utility;
		using namespace core;
		using namespace kinematics;
		using namespace scoring;
		using namespace import_pose;
		using namespace pose;
		using namespace protocols;
		using namespace minimization_packing;

		// initialize core
		devel::init( argc, argv );

		// declare variables
		Pose pose, ref;

		// Make a test pose.
		//make_pose_from_sequence( pose, "AAAAAAAAAA", core::chemical::FA_STANDARD );

		pose_from_file( pose, "/home/labonte/Workspace/Carbohydrates/4NCO_fixed3.pdb" , core::import_pose::PDB_file);
		core::io::carbohydrates::dump_gws( pose, "/home/labonte/Workspace/Carbohydrates/4NCO_fixed3.gws" );
		cout << "GWS file generated." << endl;

		/*
		pose_from_file( pose, "/home/labonte/Workspace/Carbohydrates/MBP-G4_ref.pdb" , core::import_pose::PDB_file);

		vector1< int > movable_jumps( 1, 1 );
		docking::setup_foldtree( pose, "A_B", movable_jumps );

		cout << pose << endl << endl;

		//ScoreFunctionOP sf( get_score_function() );
		vector1< string > const patches( 1, "docking" );
		ScoreFunctionOP sf( ScoreFunctionFactory::create_score_function( "talaris2013", patches ) );

		cout << "Initial Score: " << ( *sf )( pose ) << endl;

		MoveMapOP mm( utility::pointer::make_shared< MoveMap >() );
		mm->set_jump( 1, true );

		MinMover jump_minimizer;
		jump_minimizer.movemap( mm );
		jump_minimizer.score_function( sf );

		docking::FaDockingSlideIntoContact slider( 1 );
		rigid::RigidBodyPerturbMover perturber( 1, 2.0, 0.5 );
		rigid::RigidBodyRandomizeMover randomizer( pose, 1, rigid::partner_downstream, 360, 360, false );

		cout << "Score Before Any Moves: " << ( *sf )( pose ) << endl;

		cout << "Randomizing..." << endl;
		randomizer.apply( pose );
		cout << "Randomized" << endl;
		cout << "Sliding..." << endl;
		slider.apply( pose );
		cout << "Slid" << endl;

		cout << "Score After Initial Randomization: " << ( *sf )( pose ) << endl;

		for ( core::uint i = 1; i <= 1000; ++i ) {
		cout << "Perturbing..." << endl;
		perturber.apply( pose );
		cout << "Perturbed" << endl;
		cout << "Minimizing..." << endl;
		jump_minimizer.apply( pose );

		cout << "Minimization Score After Rigid Move " << i << ": " << ( *sf )( pose ) << endl;
		}*/

		//pose.dump_pdb( "" );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
