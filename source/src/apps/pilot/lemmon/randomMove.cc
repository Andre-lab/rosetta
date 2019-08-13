// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/** @page randomMove
Read a PDB with at least one jump, move a chain/ligand (this is a rigid body move; no packing),
print the PDB
Try:
"readPDB.cc -in::file::s <pdb file> -in::path::database <DB root dir>"
*/


/// @file   apps/pilot/lemmon/randomMove.cc
///
/// @brief This is to illustrate packing residues in a PDB
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of one PDB with at least 2 residues to pack>
/// 2) in::path::database <list of one database root directory>
/// 3) in::file::extra_res_fa <list of extra params files, one per residue type>
/// @author Gordon Lemmon (glemmon@gmail.com)

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/rigid/RigidBodyMover.hh>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	try {
		devel::init(argc, argv);

		utility::vector0<std::string> pdbs;
		{// process the options
			using namespace basic::options::OptionKeys;
			using basic::options::option;
			pdbs= option[in::file::s]();
		}
		core::pose::Pose pose; // starts NULL, coords *never* modified!
		{
			std::string pdb=pdbs[0];
			core::import_pose::pose_from_file(pose, pdb, core::import_pose::PDB_file);
		}
		///  Randomly Perturb with a default rotational magnitude of 3.0 and a translational
		///  magnitude of 8.0; You can change this with constructor args.
		protocols::rigid::RigidBodyPerturbMover mover;
		mover.apply(pose);
		{
			const std::string output("output.pdb");
			pose.dump_scored_pdb(output, *mover.scorefxn());
		}
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
