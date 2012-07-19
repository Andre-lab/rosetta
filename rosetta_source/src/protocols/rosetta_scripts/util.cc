// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/RosettaScripts/util.cc
/// @brief Utility functions useful in RosettaScripts.
/// @authors Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu),
///					Rocco Moretti (rmoretti@u.washington.edu), Eva-Maria Strauch (evas01@uw.edu)

// Unit Headers
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/task/PackerTask.hh>
// Project Headers

#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>


// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>


// Utility Headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/sql_database/types.hh>
#include <utility/vector0.hh>




// C++ headers

static basic::Tracer TR( "protocols.RosettaScripts.util" );

namespace protocols {
namespace rosetta_scripts {

using namespace core::scoring;
using namespace protocols::moves;
using namespace core;
using namespace std;
using utility::vector1;

/// @details This is essentially a shameless copy of Justin's PackRotamersMover::parse_task_operations. In truth
/// DesignRepackMover should disappear into Justin's better organized class, but this will wait... (SJF)
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

  TaskFactoryOP new_task_factory( new TaskFactory );
  if ( ! tag->hasOption("task_operations") ) return new_task_factory;
  std::string const t_o_val( tag->getOption<std::string>("task_operations") );
  typedef utility::vector1< std::string > StringVec;
  StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
	TR<<"Adding the following task operations to mover "<<tag->getName()<<" called "<<tag->getOption<std::string>( "name", "no_name" )<<":\n";
  for ( StringVec::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
        t_o_key != end; ++t_o_key ) {
    if ( data.has( "task_operations", *t_o_key ) ) {
      new_task_factory->push_back( data.get< TaskOperation * >( "task_operations", *t_o_key ) );
			TR<<*t_o_key<<' ';
    } else {
      utility_exit_with_message("TaskOperation " + *t_o_key + " not found in DataMap.");
    }
  }
	TR<<std::endl;
  return new_task_factory;
}


///option to add or refer to a Taskfactory through the datamap, similar to how to add/refer to movemap OPs (EMS)
core::pack::task::TaskFactoryOP
parse_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap & data, core::pack::task::TaskFactoryOP & task_factory /*, bool const reset_taskfactory */)
{
  using namespace core::pack::task;
  using namespace core::pack::task::operation;

  if ( tag->hasOption("task_factory" ) ){
    std::string const name( tag->getOption< std::string >("task_factory") );
    TR <<"taskfacotry name: " << name << std::endl;

    if( data.has( "TaskFactory", name ) ){
       task_factory = data.get< TaskFactory *>( "TaskFactory", name );
      TR<<"found helper task_factory: "<< name <<" for mover: "<<tag->getName()<< std::endl;
    }

  else{ // ( !data.has( "TaskFactory", name ) ){
      std::string tf_string = "TaskFactory";
      task_factory = new TaskFactory;
      data.add( tf_string , name , task_factory );
      TR<<"adding new TaskFactory to the datamap: "<< name  <<std::endl;
    }
  }

  if ( ! tag->hasOption("task_operations") ) return task_factory;

  std::string const t_o_val( tag->getOption<std::string>("task_operations") );
  typedef utility::vector1< std::string > StringVec;
  StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
  TR<<"Adding the following task operations to mover "<<tag->getName()<<" called "<<tag->getOption<std::string>( "name", "no_name" )<<":\n";

  for ( StringVec::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
     t_o_key != end; ++t_o_key ) {

    if ( data.has( "task_operations", *t_o_key ) ) {
      task_factory->push_back( data.get< TaskOperation * >( "task_operations", *t_o_key ) );
      TR<<*t_o_key<<' ';
    }
    else {
      utility_exit_with_message("TaskOperation " + *t_o_key + " not found in DataMap.");
    }
  }
  TR<<std::endl;
  return task_factory;
}

utility::vector1< core::pack::task::operation::TaskOperationOP >
get_task_operations( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data )
{
	using core::pack::task::operation::TaskOperationOP;
	using core::pack::task::operation::TaskOperation;
	typedef std::string String;

	utility::vector1< TaskOperationOP > task_operations;
	String const t_o_val( tag->getOption<String>("task_operations", "" ) );
	if( t_o_val != "" ){
		utility::vector1< String > const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( utility::vector1< String >::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
				t_o_key != end; ++t_o_key ) {
			if ( data.has( "task_operations", *t_o_key ) ) {
				task_operations.push_back( data.get< TaskOperation* >( "task_operations", *t_o_key ) );
			} else {
				utility_exit_with_message("TaskOperation " + *t_o_key + " not found in DataMap.");
			}
		}
	}
	return task_operations;
}



/// @details Utility function to find a scorefunction from parser-provided data. This is essentially a shameless
/// copy of Justin's PackRotamersMover::parse_score_function.
core::scoring::ScoreFunctionOP
parse_score_function( utility::tag::TagPtr const tag, protocols::moves::DataMap const & data, std::string const dflt_key/*="score12"*/ )
{
	std::string const scorefxn_key( tag->getOption<std::string>("scorefxn", dflt_key ) );
	if ( ! data.has( "scorefxns", scorefxn_key ) ) {
		utility_exit_with_message("ScoreFunction " + scorefxn_key + " not found in DataMap.");
	}
	return data.get< ScoreFunction* >( "scorefxns", scorefxn_key );
}

core::pose::PoseOP
saved_reference_pose( utility::tag::TagPtr const in_tag, protocols::moves::DataMap & data_map ){

	if( in_tag->hasOption("reference_name") ){
		core::pose::PoseOP refpose(NULL);
		std::string refpose_name(in_tag->getOption<std::string>( "reference_name") );

		if( !data_map.has("spm_ref_poses",refpose_name) ){
			refpose = new core::pose::Pose();
			data_map.add("spm_ref_poses",refpose_name,refpose );
		}
		else refpose = data_map.get<core::pose::Pose *>("spm_ref_poses",refpose_name );

		return refpose;
	}
	else std::cerr << "WARNING: saved_reference_pose function called even though tag has no 'reference_name' entry. something's unclean somewhere." << std::endl;
	return NULL;
}

/// @brief utility function for parse_movemap which goes over each of the tags in a movemap section
void
foreach_movemap_tag( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP mm ){
	using namespace core::kinematics;
	using namespace utility::tag;

	foreach( TagPtr const tag, in_tag->getTags() ){
		std::string const name( tag->getName() );
		runtime_assert( name == "Jump" || name == "Chain" || name == "Span" );
		if( name == "Jump" ){
			core::Size const num( tag->getOption< core::Size >( "number" ) );
			bool const setting( tag->getOption< bool >( "setting" ) );
			if( num == 0 ) mm->set_jump( setting ); // set all jumps if number==0
			else mm->set_jump( num, setting );
		}
		if( name == "Chain" ){
			core::Size const num( tag->getOption< core::Size >( "number" ) );
			bool const chi( tag->getOption< bool >( "chi" ) );
			bool const bb( tag->getOption< bool >( "bb" ) );
			core::Size const chain_begin( pose.conformation().chain_begin( num ) );
			core::Size const chain_end( pose.conformation().chain_end( num ) );
			for( core::Size i( chain_begin ); i <= chain_end; ++i ){
				mm->set_chi( i, chi );
				mm->set_bb( i, bb );
			}
		}
		if( name == "Span" ){
			core::Size const begin( tag->getOption< core::Size >( "begin" ) );
			core::Size const end( tag->getOption< core::Size >( "end" ) );
			runtime_assert( end >= begin );
			bool const chi( tag->getOption< bool >( "chi" ) );
			bool const bb( tag->getOption< bool >( "bb" ) );
			for( core::Size i( begin ); i <= end; ++i ){
				mm->set_chi( i, chi );
				mm->set_bb( i, bb );
			}
		}
	}//foreach tag
}

/// @brief variant of parse_movemap that takes in a datamap and searches it for already existing movemaps
void
parse_movemap( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP & mm, protocols::moves::DataMap & data, bool const reset_map ){
	using utility::tag::TagPtr;
	using namespace core::kinematics;

	if( in_tag() == NULL ) return;

	utility::vector1< TagPtr > const branch_tags( in_tag->getTags() );
	utility::vector1< TagPtr >::const_iterator tag_it;
	for( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ){
		if( (*tag_it)->getName() == "MoveMap" ){
			break;
		}
	}
	if( reset_map ){
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( true );
	}
	if( tag_it == branch_tags.end() ) return;

	if( (*tag_it)->hasOption( "name" ) ){
		std::string const name( (*tag_it)->getOption< std::string >( "name" ) );
		if( data.has( "movemaps", name ) ){
			mm = data.get< MoveMap * >( "movemaps", name );
			TR<<"Found movemap "<<name<<" on datamap"<<std::endl;
		}
		else{
			data.add( "movemaps", name, mm );
			TR<<"Adding movemap "<<name<<" to datamap"<<std::endl;
		}
	}
	foreach_movemap_tag( *tag_it, pose, mm );
}

///@details modifies an existing movemap according to tag
/// the movemap defaults to move all bb, chi, and jumps.
void
parse_movemap( utility::tag::TagPtr const in_tag, core::pose::Pose const & pose, core::kinematics::MoveMapOP mm ){
	using utility::tag::TagPtr;
	using namespace core::kinematics;

	if( in_tag() == NULL ) return;

	utility::vector1< TagPtr > const branch_tags( in_tag->getTags() );
	utility::vector1< TagPtr >::const_iterator tag_it;
	for( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ){
		if( (*tag_it)->getName() == "MoveMap" ){
			break;
		}
	}
	mm->set_bb( true );
	mm->set_chi( true );
	mm->set_jump( true );
	if( tag_it == branch_tags.end() ) return;

	if( (*tag_it)->hasOption( "name" ) ){
		TR<<"ERROR in "<<*tag_it<<'\n';
		utility_exit_with_message( "Tag called with option name but this option is not available to this mover. Note that FastRelax cannot work with a prespecified movemap b/c its movemap is parsed at apply time. Sorry." );
	}
	if( tag_it == branch_tags.end() ) return;

	foreach_movemap_tag( *tag_it, pose, mm );
}

protocols::filters::FilterOP
parse_filter( std::string const filter_name, protocols::filters::Filters_map const & filters ){
  protocols::filters::Filters_map::const_iterator filter_it( filters.find( filter_name ) );
  if( filter_it == filters.end() )
    utility_exit_with_message( "Filter "+filter_name+" not found" );
  return filter_it->second;
}

protocols::moves::MoverOP
parse_mover( std::string const mover_name, protocols::moves::Movers_map const & movers ){
  protocols::moves::Movers_map::const_iterator mover_it( movers.find( mover_name ) );
  if( mover_it == movers.end() )
    utility_exit_with_message( "Mover "+mover_name+" not found" );
  return mover_it->second;
}

/// @brief utility function for parsing xyzVector
numeric::xyzVector< core::Real >
parse_xyz_vector( utility::tag::TagPtr const xyz_vector_tag ){
	if ( ! xyz_vector_tag->hasOption("x") ) utility_exit_with_message("xyz_vector requires 'x' coordinates option");
	if ( ! xyz_vector_tag->hasOption("y") ) utility_exit_with_message("xyz_vector requires 'y' coordinates option");
	if ( ! xyz_vector_tag->hasOption("z") ) utility_exit_with_message("xyz_vector requires 'z' coordinates option");

	numeric::xyzVector< core::Real > xyz_v (
		xyz_vector_tag->getOption<core::Real>("x"),
		xyz_vector_tag->getOption<core::Real>("y"),
		xyz_vector_tag->getOption<core::Real>("z")
	);

	return xyz_v;

}

///@detail build database connection from options in a tag, this is useful make sure the fields for constructing a database connection are consistent across different tags.
utility::sql_database::sessionOP
parse_database_connection(
	utility::tag::TagPtr const tag
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys::inout;
	using utility::sql_database::DatabaseSessionManager;

	utility::sql_database::DatabaseMode::e database_mode;
	if(tag->hasOption("database_mode")){
		database_mode = utility::sql_database::database_mode_from_name(
			tag->getOption<string>("database_mode"));
	} else {
		database_mode = utility::sql_database::database_mode_from_name(
			option[dbms::mode]);
	}

	std::string database_name;
	if(tag->hasOption("database_name")){
		database_name = tag->getOption<string>("database_name");
	} else {
		database_name = option[dbms::database_name];
	}

	std::string database_pq_schema;
	if(tag->hasOption("database_pq_schema")){
		database_pq_schema = tag->getOption<string>("database_pq_schema");
	} else {
		database_pq_schema = option[dbms::pq_schema];
	}

	switch(database_mode){

	case utility::sql_database::DatabaseMode::mysql:
		if(tag->hasOption("database_pq_schema")){
			TR << "WARNING: You must specify 'database_mode=postgres' ";
			TR << "to use the 'database_pq_schema' tag." << endl;
		}
		break;
	case utility::sql_database::DatabaseMode::postgres:
		if(tag->hasOption("database_separate_db_per_mpi_process")){
			TR << "WARNING: You must specify 'database_mode=sqlite3' ";
			TR << "to use the 'database_separate_db_per_mpi_process' tag." << endl;
		}
		if(tag->hasOption("database_read_only")){
			TR << "WARNING: You must specify 'database_mode=sqlite3' ";
			TR << "to use the 'database_read_only' tag." << endl;
		}
		break;

	case utility::sql_database::DatabaseMode::sqlite3:
		if(tag->hasOption("database_host")){
			TR << "WARNING: You must specify either 'database_mode=mysql' ";
			TR << "or database_mode=postgres' to use the 'database_host' tag." << endl;
		}

		if(tag->hasOption("database_user")){
			TR << "WARNING: You must specify either 'database_mode=mysql' ";
			TR << "or database_mode=postgres' to use the 'database_user' tag." << endl;
		}

		if(tag->hasOption("database_password")){
			TR << "WARNING: You must specify either 'database_mode=mysql' ";
			TR << "or database_mode=postgres' to use the 'database_password' tag." << endl;
		}

		if(tag->hasOption("database_port")){
			TR << "WARNING: You must specify either 'database_mode=mysql' ";
			TR << "or database_mode=postgres' to use the 'database_port' tag." << endl;
		}
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(database_mode) + "'");
	}

	switch(database_mode){
	case utility::sql_database::DatabaseMode::sqlite3:
		return DatabaseSessionManager::get_instance()->get_db_session(
			database_mode, database_name, "", "", "", "", 0,
			tag->getOption("database_read_only", false),
			tag->getOption("database_separate_db_per_mpi_process", false));

	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres:

		if(!tag->hasOption("database_host") && !option[dbms::host].user()){
			TR << "WARNING: To connect to a postgres or mysql database you must set ";
			TR << "the database_host tag or specify -dbms:host on the command line." << endl;
		}

		if(tag->hasOption("database_user") && !option[dbms::user].user()){
			TR << "WARNING: To connect to a postgres or mysql database you must set ";
			TR << "the database_user tag or specify -dbms:user on the command line." << endl;
		}

		if(tag->hasOption("database_password") && !option[dbms::password].user()){
			TR << "WARNING: To connect to a postgres or mysql database you must set ";
			TR << "the database_password tag or specify -dbms:password on the command line." << endl;
		}

		if(tag->hasOption("database_port") && !option[dbms::port].user()){
			TR << "WARNING: To connect to a postgres or mysql database you must set ";
			TR << "the database_port tag or specify -dbms:port on the command line." << endl;
		}

		return DatabaseSessionManager::get_instance()->get_db_session(
			database_mode, database_name, database_pq_schema,
			tag->getOption<string>("database_host", option[dbms::host]),
			tag->getOption<string>("database_user", option[dbms::user]),
			tag->getOption<string>("database_password", option[dbms::password]),
			tag->getOption<Size>("database_port", option[dbms::port]));

	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(database_mode) + "'");
	}
	return 0;
}


/// @brief Return the number of the residue on source that is nearest to res on target. If the distance
/// is greater than 2.0 returns 0 to indicate error
core::Size
find_nearest_res( core::pose::Pose const & source, core::pose::Pose const & target, core::Size const res, core::Size const chain/*=0*/ ){
  core::Real min_dist( 100000 ); core::Size nearest_res( 0 );
  for( core::Size i = 1; i <= source.total_residue(); ++i ){
		if( source.residue( i ).is_ligand() ) continue;
		if( chain && source.residue( i ).chain() != chain ) continue;
    core::Real const dist( target.residue( res ).xyz( "CA" ).distance( source.residue( i ).xyz( "CA" ) ) );
    if( dist <= min_dist ){
      min_dist = dist;
      nearest_res = i;
    }
  }
  if( min_dist <= 2.0 ) return nearest_res;
  else return 0;
}

utility::vector1< core::Size >
residue_packer_states( core::pose::Pose const & pose, core::pack::task::TaskFactoryCOP tf, bool const designable, bool const packable/*but not designable*/) {
  utility::vector1< core::Size > designable_vec, packable_vec, both;
  designable_vec.clear(); packable_vec.clear(); both.clear();
  core::pack::task::PackerTaskOP packer_task( tf->create_task_and_apply_taskoperations( pose ) );
  for( core::Size resi=1; resi<=pose.total_residue(); ++resi ){
    if( packer_task->being_designed( resi ) )
      designable_vec.push_back( resi );
		else if( packer_task->being_packed( resi ) )
			packable_vec.push_back( resi );
	}
	if( designable && packable ){
		both.insert( both.begin(), designable_vec.begin(), designable_vec.end() );
		both.insert( both.end(), packable_vec.begin(), packable_vec.end() );
 		return both;
 	}
	if( designable )
		return designable_vec;
	return packable_vec;
}
} //RosettaScripts
} //protocols
