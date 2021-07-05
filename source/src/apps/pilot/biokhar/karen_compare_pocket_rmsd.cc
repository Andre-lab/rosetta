// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk + dj

#include <iostream>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <devel/init.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <core/scoring/rms_util.hh>

// Utility Headers
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <fstream> // AUTO IWYU For ifstream

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//Need -interface_list -first_structure -second_structure
OPT_KEY( String, interface_list )
OPT_KEY( String, first_structure )
OPT_KEY( String, second_structure )

static basic::Tracer TR( "apps.pilot.karen_compare_pocket_rmsd.main" );

//set to store pdb info keys
std::set <std::string> interface;

//stores resid of the ligand residue
core::Size lig_res_num;

bool is_interface_heavyatom(
	core::pose::Pose const & pose,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	// ws get residue "key" for set
	std::ostringstream residuestream;
	residuestream << pose.pdb_info()->chain(resno) << pose.pdb_info()->number(resno);
	std::string res_id = residuestream.str();

	core::conformation::Residue const & rsd = pose.residue(resno);

	if ( interface.count( res_id ) > 0 ) return rsd.is_protein() && !rsd.atom_is_hydrogen(atomno);

	return false;
}

/// General testing code
int
main( int argc, char * argv [] ){

	try{

		NEW_OPT( interface_list, "interface residues", "interface residues" );
		NEW_OPT( first_structure, "first structure", "-1" );
		NEW_OPT( second_structure, "comparison structure", "-1" );

		devel::init(argc, argv);

		std::string const ifilename = option[ interface_list ];
		if ( ifilename != "" ) {
			std::ifstream ifs(ifilename.c_str(), std::ifstream::in);
			if ( !ifs.is_open() ) {
				std::cout<< "Error opening contact list file "<<ifilename<<std::endl;
				return -100;
			}
			std::string intres;
			while ( ifs.good() ) {
				ifs >> intres;
				interface.insert(intres);
			}
		}

		std::string const structure1 = option[ first_structure ] ;
		std::string const structure2 = option[ second_structure ] ;

		TR << "Starting recomputing scores and rmsds" << std::endl;

		// create pose from pdb
		pose::Pose pose1;
		core::import_pose::pose_from_file( pose1, structure1 , core::import_pose::PDB_file);
		pose::Pose pose2;
		core::import_pose::pose_from_file( pose2, structure2 , core::import_pose::PDB_file);


		TR << "Defined interface" << std::endl;

		core::Real CA_rms = rmsd_with_super( pose1, pose2, is_protein_CA );
		core::Real heavyatom_rms = rmsd_with_super( pose1, pose2, is_interface_heavyatom );

		TR << "All residue rmsd: " << CA_rms << " Interface rmsd: " << heavyatom_rms <<std::endl;


		TR << "Done computing rmsds" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;

}
