// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/mads/shapedocking.cc
/// @brief  shapedocking test
/// @author Mads Jeppesen

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/simple_filters/TotalSasaFilter.hh>
#include <protocols/score_filters/ScoreTypeFilter.hh>

#include <protocols/rigid/RigidBodyDofAdaptiveMover.hh>
#include <protocols/rigid/RigidBodyDofGridMover.hh>
#include <protocols/symmetric_docking/SymDockAdaptiveMover.hh>
#include <core/kinematics/Jump.hh>

int main(int argc, char ** argv) {

    std::cout << "Hello";
    std::cout << "Hello Wo~rld!" << std::endl;
    devel::init( argc, argv );
    utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
//    std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;

    // import pose
    core::pose::PoseOP pose = core::import_pose::pose_from_file( "/Users/mads/Rosetta/main/source/src/apps/pilot/mads/1stm.cif" );
    pose->dump_pdb("initial.pdb");

    // GIVES A PROBLEM WITH PYMOL!!!!!!
//    // score function
    core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function(); // talaris or ref2015?
    core::Real score = sfxn->score( *pose );
    std::cout << "this is the score " << score  << std::endl;

    // Try to first set up the symmetry before mc (I think it replaces it if not)
    protocols::symmetry::SetupForSymmetryMover setup_symmetry = protocols::symmetry::SetupForSymmetryMover("/Users/mads/Rosetta/main/source/src/apps/pilot/mads/ico.symm");
    setup_symmetry.apply(*pose);

//    pose->dump_pdb("initial_symmetric.pdb");

//    protocols::rigid::RigidBodyDofGridMoverOP gridmover = protocols::rigid::RigidBodyDofGridMoverOP( new protocols::rigid::RigidBodyDofGridMover());

//
    protocols::moves::PyMOLMover pmm = protocols::moves::PyMOLMover();
    pmm.keep_history(true);
//    pmm.apply(*pose);

    protocols::symmetric_docking::SymGridMover symgridmover = protocols::symmetric_docking::SymGridMover();
    symgridmover.add_jump(*pose, "JUMP5fold1", "z",-5,5,2, "uniform");
    symgridmover.add_jump(*pose, "JUMP5fold1", "angle_z",-5,5,2, "gauss");
    symgridmover.add_jump(*pose, "JUMP5fold111", "x",-5,5,2, "grid");
//    symgridmover.add_jump(*pose, "JUMP5fold1111", "x_angle",-10,10,2, "uniform");
//    symgridmover.add_jump(*pose, "JUMP5fold1111", "y_angle",-10,10,2, "grid");
//    symgridmover.add_jump(*pose, "JUMP5fold1111", "z_angle",-10,10,2, "uniform");
//    symgridmover.add_mover(gridmover);

    // says this ERROR: Must define core::pose::metrics::simple_calculators::SasaCalculatorLegacy with name 'sasa' in order to use TotalSasaFilter.
//    protocols::filters::FilterOP totalsasaop = protocols::filters::FilterOP(new protocols::simple_filters::TotalSasaFilter());
    //  core::scoring::ScoreFunctionCOP scorefxn, core::scoring::ScoreType const score_type, core::Real const score_type_threshold
    // 1 is fa_atr
//    protocols::filters::FilterOP totalscoreop = protocols::filters::FilterOP(new protocols::score_filters::ScoreTypeFilter(sfxn, core::scoring::ScoreType(1), 0));

//    symgridmover.add_filter(totalsasaop);
//    symgridmover.add_filter(totalscoreop);
    while (true) {
        symgridmover.apply(*pose);
//        pmm.apply(*pose);
    }

//    symgridmover.go_to_grid_point(*pose, 20);
//    symgridmover.go_to_grid_point(*pose, 25);


//    for (int i  = 0; i != 2*2*2*2*2*2; i++) {
//        gridmover.apply(*pose);
////        pmm.apply(*pose);
//    }



//    protocols::symmetric_docking::SymGridMover symgridmover = protocols::symmetric_docking::SymGridMover();
//    symgridmover.add_mover(gridmover); //  ,100, 5000);
//
//    gridmover.apply(*pose);
//    gridmover.apply(*pose);
//    gridmover.apply(*pose);
//    gridmover.go_to_grid_point(*pose, 10);
//    gridmover.go_to_grid_point(*pose, 100);
//    gridmover.go_to_grid_point(*pose, 1000);
//    gridmover.go_to_grid_point(*pose, 1010);
//    gridmover.go_to_grid_point(*pose, 7999);
//    gridmover.go_to_grid_point(*pose, 8000);
//    gridmover.go_to_grid_point(*pose, 7444);
//    gridmover.go_to_grid_point(*pose, 10);
//
//    std::cout << gridmover.get_n_total_grid_points() << std::endl;
//    std::cout << counts << std::endl;




    // MC
//    protocols::moves::MonteCarloOP mc = utility::pointer::make_shared< protocols::moves::MonteCarlo >(*pose, *sfxn, 1.0 );
//    pose->dump_pdb("initial_symmetric_MC.pdb");

    // rigidbodydofadaptive
    protocols::rigid::RigidBodyDofAdaptiveMoverOP rb_mover1 = utility::pointer::make_shared< protocols::rigid::RigidBodyDofAdaptiveMover>("test1");
//    rb_mover1->add_jump(*pose, "JUMP5fold111_subunit", "x_angle", "add", 10); // 0.5, 0.01, true, 5 ,-5);
    rb_mover1->add_jump(*pose, "JUMP5fold111", "x_angle", "add", 10); //, 0.0, 0.01, true, 5 ,-5);
//    rb_mover1->add_jump(*pose, "JUMP5fold111_subunit", "z_angle", "add", 10); // ,0.0, 0.01, true, 5, -5);
//    auto jump1_17 = pose->jump(17);
//    auto jump1_15 = pose->jump(15);
//    auto jump1_3 = pose->jump(3);
    rb_mover1->apply(*pose);
    rb_mover1->apply(*pose);
    rb_mover1->apply(*pose);
////    protocols::rigid::RigidBodyDofAdaptiveMoverOP rb_mover2 = utility::pointer::make_shared< protocols::rigid::RigidBodyDofAdaptiveMover>("test2");
////    rb_mover2->add_jump("JUMP5fold111", "x");
//
//    // symdockadaptive
//    protocols::symmetric_docking::SymDockAdaptiveMover symdock_mover = protocols::symmetric_docking::SymDockAdaptiveMover();
//    symdock_mover.set_docking_options();
//    symdock_mover.add_protocol("Random", 1000, sfxn);
//    symdock_mover.add_mover(rb_mover1, true);
//////    symdock_mover.add_mover(rb_mover2, true);
//////
//////
////////    pose->dump_pdb("before.pdb");
//    symdock_mover.apply(*pose);

    protocols::symmetry::ExtractAsymmetricPoseMover extractsym = protocols::symmetry::ExtractAsymmetricPoseMover();

    auto jump2_17 = pose->jump(17);
    auto jump2_15 = pose->jump(15);
    auto jump2_3 = pose->jump(3);

    extractsym.apply(*pose);
    pose->dump_pdb("after.pdb");



    return 0;
}
