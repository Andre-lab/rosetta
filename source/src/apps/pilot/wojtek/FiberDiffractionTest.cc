// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/wojtek/FiberDiffractionTest.cc
/// @brief  FiberDiffraction simple test for derivatives
/// @author Wojciech Potrzebowski and Ingemar Andre



#include <protocols/rigid/RigidBodyMover.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>

#include <core/kinematics/MoveMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <devel/init.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>

#include <core/import_pose/import_pose.hh>
// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rescore.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// C++ headers
//#include <cstdlib>
#include <string>

using namespace protocols::moves;
using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

class MyScoreMover : public Mover {
public:
	MyScoreMover();

	void apply( core::pose::Pose& pose ) override;
	std::string get_name() const override { return "MyScoreMover"; }

	MoverOP clone() const override { return utility::pointer::make_shared< MyScoreMover >( *this ); }

	MoverOP fresh_instance() const override { return utility::pointer::make_shared< MyScoreMover >(); }

	void set_keep_input_scores(){ keep_scores_flag_ = true; }
	void set_skip_scoring(){ skip_scoring_ = true; }
private:
	core::scoring::ScoreFunctionOP sfxn_;

	bool keep_scores_flag_{false};   // retains the previous scores from the input silent file or whatever
	bool skip_scoring_{false};       // skips the actual scoring call, calling evaluators only
};

MyScoreMover::MyScoreMover()
{
	//using namespace basic::options;
	//using namespace basic::options::OptionKeys;
	//using namespace core;

	// get scorefxn and add constraints if defined
	sfxn_ = core::scoring::get_score_function();
	if ( option[ in::file::keep_input_scores ]() ) {
		set_keep_input_scores();
	}
	if ( option[ rescore::skip ]() ) {
		set_skip_scoring();
	}

}

void MyScoreMover::apply( core::pose::Pose& pose ) {
	if ( !keep_scores_flag_ ) {
		pose.energies().clear();
		pose.data().clear();
	}
	sfxn_->set_weight( core::scoring::linear_chainbreak, 4.0/3.0 );
	sfxn_->set_weight( core::scoring::overlap_chainbreak, 1.0 );

	if ( ! skip_scoring_ ) {
		(*sfxn_)( pose );
		//std::cout <<  (*sfxn_)( pose ) << std::endl;
	}
}


// option key includes

using namespace core;

int
main( int argc, char * argv [] )
{

	try {

		devel::init(argc, argv);
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, option[ in::file::s ]()[1] , core::import_pose::PDB_file);

		MoverOP scoremover (new MyScoreMover);
		// set pose for symmetry
		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::symmetry::SetupForSymmetryMover symmover = protocols::symmetry::SetupForSymmetryMover();
			symmover.apply(pose);
		}

		scoremover->apply(pose);

		using SymmetricConformation = core::conformation::symmetry::SymmetricConformation;

		auto & symm_conf (
			dynamic_cast<SymmetricConformation & > ( pose.conformation()) );

		std::map< Size, core::conformation::symmetry::SymDof > dofs = symm_conf.Symmetry_Info()->get_dofs();
		std::map< Size, core::conformation::symmetry::SymDof >::iterator dof_iterator;
		std::map< Size, core::conformation::symmetry::SymDof >::iterator it;
		auto it_begin = dofs.begin();
		auto it_end = dofs.end();

		for ( it = it_begin; it != it_end; ++it ) {
			core::conformation::symmetry::SymDof dof ( (*it).second );
			if ( dof.allow_dof(1) ) {
				dof_iterator = it;
			}
		}

		protocols::rigid::RigidBodyDofTransMover dofmover( (*dof_iterator).second, (*dof_iterator).first, 0.5 );
		dofmover.apply(pose);
		pose.dump_pdb("trans.pdb");

		core::kinematics::MoveMapOP movemap ( new core::kinematics::MoveMap());
		movemap->set_jump( false );
		movemap->set_jump((*dof_iterator).first ,true);
		movemap->set_bb( false );
		movemap->set_chi( false );


		protocols::minimization_packing::MinMoverOP min_mover(
			new protocols::minimization_packing::symmetry::SymMinMover( movemap, core::scoring::get_score_function(), "lbfgs_armijo_nonmonotone",0.2, false, false, false) ); //Last two should be true to have deriv check
		min_mover->apply(pose);
		scoremover->apply(pose);
		pose.dump_pdb("min.pdb");

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}
