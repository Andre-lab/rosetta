// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @detailed
///
///
/// @author Nikolas Sgourakis



// Unit Headers
#include <protocols/evaluation/PalesEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/evaluation.OptionKeys.gen.hh>


// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>

//Auto Headers
#include <utility/options/keys/BooleanOptionKey.hh>

// C++ headers
#ifdef NATCL 
#define system(a) 1
#endif

static basic::Tracer tr("protocols.evaluation.PalesEvaluator");

namespace protocols {
namespace evaluation {

using namespace core;



PalesEvaluator::PalesEvaluator( std::string tag, std::string pales_rdc_file )
  : ExternalEvaluator( tag )
{
	#ifdef WIN32
	utility_exit_with_message("don't use PalesEvaluator on a BillBox");
  #endif
	//	if (!utility::file::file_exists( scratch_dir()+"/SPARTA" ) ) {
	// std::string command( "cp -Rf $HOME/SPARTA "+scratch_dir());
	std::string command( "rsync -azvu $HOME/pales "+scratch_dir());
	int ret(system(command.c_str()));
	if( ret ){
		utility_exit_with_message("System command failed:'" + command + "'" );
	}


	//  std::string command2( "rsync -azvu $HOME/scripts/calculate_pales_rms.pl "+scratch_dir());
	// int ret2(system(command2.c_str()));
	//if( ret2 ){
	//utility_exit_with_message("System command failed:'" + command2 + "'" );
	// }

	//}
	set_command( scratch_dir()+"/pales -inD "+pales_rdc_file+" -pdb __POSE.pdb -stPAles -wv 0.05 | grep \"DATA RMS\" |awk '{print $3}' > __RESULT" );
}

bool PalesEvaluator::applicable( pose::Pose const& pose ) const {
  return pose.is_fullatom();
}


}
}
