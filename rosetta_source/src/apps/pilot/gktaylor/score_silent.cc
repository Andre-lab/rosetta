// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

// libRosetta headers


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/Func.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>



using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



int
main( int argc, char* argv [] )
{
	// options, random initialization
	devel::init( argc, argv );

	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	if ( option[ in::file::fullatom ]() ) {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
	}

	// configure score function
  core::scoring::ScoreFunctionOP energy12scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "standard", "score12" ) );
	core::scoring::ScoreFunctionOP constraintscorefxn( new core::scoring::ScoreFunction );
	constraintscorefxn->set_weight( atom_pair_constraint, 1.0 );

	core::pose::Pose pose;

	ConstraintSetOP cstset_;
	std::string cstfile = option[ basic::options::OptionKeys::constraints::cst_file ]();
	if ( ! utility::file::file_exists( cstfile ) ) {
		utility_exit_with_message( "Error: can't read file " + cstfile );
	}

	// configure silent-file data object
	core::io::silent::SilentFileData sfd;
	if ( option[ in::file::silent ].user() ) {
		sfd.read_file( *(option[ in::file::silent ]().begin()) );
	}

	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native_pose, *rsd_set, option[ in::file::native ]() );
	}

	for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end(); iter != end; ++iter ) {
		iter->fill_pose( pose, *rsd_set );
		cstset_ = ConstraintIO::read( cstfile, new ConstraintSet, pose );

		core::io::silent::ProteinSilentStruct ss( pose, iter->decoy_tag(), option[ in::file::fullatom ]() );
		if ( option[ in::file::native ].user() ) {
			ss.add_energy( "CA_rmsd",   core::scoring::CA_rmsd  ( native_pose, pose ) );
			ss.add_energy( "CA_maxsub", core::scoring::CA_maxsub( native_pose, pose ) );
		}
		sfd.write_silent_struct( ss, option[ out::file::silent ]() );

		std::string pdbname = iter->decoy_tag();

		if ( !cstset_ ) {
			cstset_ = ConstraintIO::read( cstfile, new ConstraintSet, pose );
		} else {
			pose.constraint_set( cstset_ );
		}

		core::Real constraint_score = (*constraintscorefxn)(pose);
		core::Real score12_score = (*energy12scorefxn)(pose);
		core::Real gdtmm_score = CA_gdtmm(pose, native_pose);
		std::cout << pdbname << '\t' << gdtmm_score << '\t' << constraint_score << '\t' << score12_score << '\n';
	}

	// and finally, score the native
	core::Real constraint_score = (*constraintscorefxn)(native_pose);
	core::Real score12_score = (*energy12scorefxn)(native_pose);
	core::Real gdtmm_score = CA_gdtmm(native_pose, native_pose);
	std::cout << "NATIVE" << '\t' << gdtmm_score << '\t' << constraint_score << '\t' << score12_score << '\n';

	return 0;
}

