// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit Headers
#include <devel/domain_assembly/DomainAssemblyReader.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pack/task/ResfileReader.hh>

// Project Headers

// Utility Headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

#include <devel/domain_assembly/domain_assembly_setup.hh> // AUTO IWYU For DomainInfo

using basic::Error;
using basic::Warning;

//STL headers
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>

namespace devel {
namespace domain_assembly {

using namespace core;

//////////////////////
void
PDB::domain_action(
	utility::vector1< std::string > const & tokens,
	core::Size & which_token,
	DomainInfo & domain
) const

{
	debug_assert( tokens[ which_token ] == name() );
	std::string const & pdb_path = tokens[ ++which_token ];
	core::pose::Pose temp_pose;
	std::cout << "Reading in PDB file " << pdb_path << std::endl;
	core::import_pose::pose_from_file( temp_pose, pdb_path.c_str() , false, core::import_pose::PDB_file);
	domain.set_input_pose(temp_pose);
	++which_token;
}
/////////////
void
NTermLinker::domain_action(
	utility::vector1< std::string > const & tokens,
	core::Size & which_token,
	DomainInfo & domain
) const
{
	debug_assert( tokens[ which_token ] == name() );
	std::string linker = tokens[ ++which_token ];
	domain.set_Nterm_linker( linker );
	++which_token;
}
////////////////
void
CTermLinker::domain_action(
	utility::vector1< std::string > const & tokens,
	core::Size & which_token,
	DomainInfo & domain
) const
{
	debug_assert( tokens[ which_token ] == name() );
	std::string linker = tokens[ ++which_token ];
	domain.set_Cterm_linker( linker );
	++which_token;
}
/////////////////
void
NTrim::domain_action(
	utility::vector1< std::string > const & tokens,
	core::Size & which_token,
	DomainInfo & domain
) const
{
	debug_assert( tokens[ which_token ] == name() );
	std::string trim = tokens[ ++which_token ];
	core::Size ntrim = atoi(trim.c_str());
	domain.set_trim_nterm( ntrim);
	++which_token;
}
/////////////////////
void
CTrim::domain_action(
	utility::vector1< std::string > const & tokens,
	core::Size & which_token,
	DomainInfo & domain
) const
{
	debug_assert( tokens[ which_token ] == name() );
	std::string trim = tokens[ ++which_token ];
	core::Size ctrim = atoi(trim.c_str());
	domain.set_trim_cterm( ctrim );
	++which_token;
}
//////////////////////////////

//@details this creates a map linking the parsed strings from the da_option_file
//to the command objects.  NEW COMMANDS MUST BE ADDED HERE, HARD CODED
std::map< std::string, DomainAssemblyCommandOP >
create_command_map()
{
	using namespace std;

	map< string, DomainAssemblyCommandOP > command_map;
	command_map[ PDB::name() ] = utility::pointer::make_shared< PDB >();
	command_map[ NTermLinker::name() ] = utility::pointer::make_shared< NTermLinker >();
	command_map[ CTermLinker::name() ] = utility::pointer::make_shared< CTermLinker >();
	command_map[ NTrim::name() ] = utility::pointer::make_shared< NTrim >();
	command_map[ CTrim::name() ] = utility::pointer::make_shared< CTrim >();

	return command_map;
}
//@brief utility function for DomainAssembly reader (checks for a leading # signaling a comment)
bool
comment_begin( utility::vector1< std::string > const & tokens, core::Size which_token )
{
	return tokens[ which_token ][ 0 ] == '#';
}

/// @details Auto-generated virtual destructor
DomainAssemblyCommand::~DomainAssemblyCommand() = default;

void
parse_da_option_file( utility::vector1< DomainInfo > & domains, std::string filename )
{
	using namespace std;
	map< string, DomainAssemblyCommandOP > command_map = create_command_map();

	// T("DomainAssemblyReader") << "Reading da_option-file: " << filename << std::endl;
	ifstream da_option_file( filename.c_str() );
	//if ( ! da_option_file ) { Error() << "Domain Assembly reader could not find a file named" << filename << std::endl; utility_exit(); }

	int lineno = 0;
	while ( da_option_file ) {
		utility::vector1< std::string > tokens( core::pack::task::tokenize_line( da_option_file ));
		++lineno;
		core::Size which_token = 1;
		core::Size ntokens( tokens.size() );

		//ignore blank lines
		if ( ntokens == 0 ) continue;

		if ( core::pack::task::comment_begin( tokens, which_token ) ) { continue; } // ignore the rest of this line
		DomainInfo domain;
		while ( which_token <= ntokens ) {
			DomainAssemblyCommandOP command = command_map[ tokens[ which_token ] ];
			if ( !command ) {
				Error() << "da_option_file ERROR: line: " << lineno << " command not found: " << tokens[ which_token ] << std::endl;
				utility_exit();
			}
			command->domain_action( tokens, which_token, domain );
		}
		domains.push_back(domain);
	}
}

}
}
