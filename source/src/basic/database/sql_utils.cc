// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/database/sql_utils.cc
/// @brief Database utility functions
/// @author Sam DeLuca
/// @author Matthew O'Meara

#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <numeric/random/random.hh>

#include <platform/types.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/string_util.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/Tracer.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/scoped_ptr.hpp>
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/predicate.hpp>
#define foreach BOOST_FOREACH

#ifdef WIN32
#include <windows.h>
#endif

using std::string;
using std::stringstream;
using utility::sql_database::sessionOP;
using platform::Size;
using cppdb::statement;
using cppdb::cppdb_error;
using cppdb::result;
using utility::vector1;
using namespace utility::sql_database;

namespace basic {
namespace database {

static basic::Tracer TR( "basic.database.sql_utils" );
static numeric::random::RandomGenerator RG(345264);




sessionOP
get_db_session() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return get_db_session(
		option[inout::dbms::database_name],
		utility::sql_database::TransactionMode::standard,
		0,
		option[inout::dbms::pq_schema]);
}


sessionOP
get_db_session(
	string const & db_name,
	string const & pq_schema
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return get_db_session(
		database_mode_from_name(option[inout::dbms::mode]),
		utility::sql_database::TransactionMode::standard,
		0,
		db_name,
		pq_schema);
}

sessionOP
get_db_session(
	string const & db_name,
	TransactionMode::e transaction_mode,
	Size chunk_size,
	string const & pq_schema
){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return get_db_session(
		database_mode_from_name(option[inout::dbms::mode]),
		transaction_mode,
		chunk_size,
		db_name,
		pq_schema);
}

utility::sql_database::sessionOP
get_db_session(
	utility::sql_database::DatabaseMode::e db_mode,
	std::string const & db_name,
	std::string const & pq_schema){

	return get_db_session(
		db_mode,
		utility::sql_database::TransactionMode::standard,
		0,
		db_name,
		pq_schema);
}

sessionOP
get_db_session(
	DatabaseMode::e db_mode,
	TransactionMode::e transaction_mode,
	Size chunk_size,
	string const & db_name,
	string const & pq_schema
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	switch(db_mode) {
	case DatabaseMode::sqlite3:

		if(
			option[inout::dbms::host].user() ||
			option[inout::dbms::user].user() ||
			option[inout::dbms::password].user() ||
			option[inout::dbms::port].user()) {
			utility_exit_with_message(
				"You have specified options for a client-server database "
				"but the database mode is sqlite3. "
				"Please specify -inout:dbms:mode <db_mode>.");
		}

		if(pq_schema.compare("")){
			TR.Warning
				<< "You have specified a postgres schema but using a sqlite3 database. "
				<< "To use postgres, please specify -inout:dbms:mode postgres"
				<< std::endl;
		}

		return DatabaseSessionManager::get_instance()->get_session_sqlite3(
			db_name,
			transaction_mode,
			chunk_size,
			option[inout::dbms::readonly],
			option[inout::dbms::separate_db_per_mpi_process]);

	case DatabaseMode::mysql:

		if(option[inout::dbms::readonly]){
			utility_exit_with_message(
				"Restricting access to a mysql database is done at the user level "
				"rather that the connection level. "
				"So requesting a readonly connection cannot fullfilled.");
		}

		if(option[inout::dbms::separate_db_per_mpi_process]){
			utility_exit_with_message(
				"The -inout:dbms:separate_db_per_mpi_process flag "
				"only applies to sqlite3 databases.");
		}

		if(pq_schema.compare("")){
			TR.Warning
				<< "You have specified a postgres schema but using a sqlite3 database. "
				<< "To use postgres, please specify -inout:dbms:mode postgres"
				<< std::endl;
		}

		if( !(
			option[inout::dbms::host].user() &&
			option[inout::dbms::user].user() &&
			option[inout::dbms::password].user() &&
			option[inout::dbms::port].user())) {
			utility_exit_with_message(
				"To connect to a mysql database you must specify "
				"-inout:dbms:host -inout:dbms:user -inout:dbms:password and "
				"-inout:dbms:port");
		}

		return DatabaseSessionManager::get_instance()->get_session_mysql(
			db_name,
			transaction_mode,
			chunk_size,
			option[inout::dbms::host],
			option[inout::dbms::user],
			option[inout::dbms::password],
			option[inout::dbms::port]);


	case DatabaseMode::postgres:

		if(option[inout::dbms::readonly]){
			utility_exit_with_message(
				"Restricting access to a postgres database is done at the user level "
				"rather that the connection level. So requesting a readonly connection "
				"cannot fullfilled.");
		}

		if(option[inout::dbms::separate_db_per_mpi_process]){
			utility_exit_with_message(
				"The -inout:dbms:separate_db_per_mpi_process flag only applies to "
				"sqlite3 databases.");
		}

		if( !(
			option[inout::dbms::host].user() &&
			option[inout::dbms::user].user() &&
			option[inout::dbms::password].user() &&
			option[inout::dbms::port].user())) {
			utility_exit_with_message(
				"To connect to a postgres database you must specify "
				"-inout:dbms:host -inout:dbms:user -inout:dbms:password "
				"and -inout:dbms:port");
		}

		return DatabaseSessionManager::get_instance()->get_session_postgres(
			db_name,
			transaction_mode,
			chunk_size,
			pq_schema,
			option[inout::dbms::host],
			option[inout::dbms::user],
			option[inout::dbms::password],
			option[inout::dbms::port]);
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" + name_from_database_mode(db_mode) + "'");
	}
	return 0;
}


statement safely_prepare_statement(
	string const & statement_string,
	sessionOP db_session)
{
	statement stmt;
	try
	{
		stmt = db_session->prepare(statement_string);
		return stmt;
	}catch(cppdb_error error)
	{
		TR.Error << " Failed to safely prepare the following statement: " << std::endl;
		TR.Error << statement_string << std::endl;
		TR.Error << error.what() << std::endl;

		utility_exit_with_message(error.what());
	}
	return stmt; //there's no way this should happen
}

void
safely_write_to_database(
	statement & statement
) {

	platform::Size retry_limit = 30;
	platform::Size cycle = 0;
	while(true)
	{
		try
		{
			statement.exec();
			return;
		}catch(cppdb::bad_value_cast & except)
		{

			utility_exit_with_message(except.what());
		}catch(cppdb::empty_row_access & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::invalid_column & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::invalid_placeholder & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::multiple_rows_query & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::not_supported_by_backend & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::null_value_fetch & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::backend_deadlock & except)
		{
			//retry until we hit the retry limit and then die
			//deadlocks sometimes happen with mysql.
			//see http://dev.mysql.com/doc/refman/5.0/en/innodb-deadlocks.html

			if(cycle < retry_limit)
			{
				TR << "Backend deadlock detected, retrying SQL statement";

#ifdef WIN32
				Sleep(1000);
#else
				//Sleep some amount between 100-2000 ms
				usleep(100+1900*RG.uniform());
#endif
			}else
			{
				utility_exit_with_message(except.what());
			}

		}catch(cppdb::cppdb_error & except)
		{
#ifdef USEMPI
			if(except.what() == "database is locked"){
				stringstream err_msg;
				err_msg
					<< "database is locked" << std::endl
					<< std::endl
					<< "PSA: If this is an sqlite3 session, running under MPI, and you are using the database primarily for writing output," << std::endl
					<< "consider using the 'separate_db_per_mpi_process' option if you're not already using it." << std::endl
					<< "To do this add separate_db_per_mpi_process=1 to a RosettaScripts/resource_definitions xml tag that takes database options" << std::endl
					<< "or add -inout:dbms:separate_db_per_mpi_process to the command line or flags file" << std::endl
					<< "This option will append '_<mpi_rank>' to the database filename for each mpi process. Once the run has finished," << std::endl
					<< "these databases can be merged together using " << std::endl
					<< std::endl
					<< "        bash /path/to/rosetta_tests/features/sample_sources/merge.sh <output_db> <input_db_part_1> [<input_db_part_2> [ ... ] ]" << std::endl
					<< std::endl
					<< "For more information see the database connection options page in the Rosetta Wiki: RosettaScripts_database_connection_options" << std::endl
					<< except.what() << std::endl;
				utility_exit_with_message(err_msg.str());
			}
#endif

			utility_exit_with_message(except.what());
		}

		cycle++;
	}
}

result
safely_read_from_database(
	statement & statement
) {
	platform::Size retry_limit = 30;
	platform::Size cycle = 0;
	while(true)
	{
		try
		{
			return statement.query();
		}catch(cppdb::bad_value_cast & except)
		{

			utility_exit_with_message(except.what());
		}catch(cppdb::empty_row_access & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::invalid_column & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::invalid_placeholder & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::multiple_rows_query & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::not_supported_by_backend & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::null_value_fetch & except)
		{
			utility_exit_with_message(except.what());
		}catch(cppdb::backend_deadlock & except)
		{
			//retry until we hit the retry limit and then die
			//deadlocks sometimes happen with mysql.
			//see http://dev.mysql.com/doc/refman/5.0/en/innodb-deadlocks.html

			if(cycle < retry_limit)
			{
				TR << "Backend deadlock detected, retrying SQL statement";

#ifdef WIN32
				Sleep(1000);
#else
				//Sleep some amount between 100-2000 ms
				usleep(100+1900*RG.uniform());
#endif
			}else
			{
				utility_exit_with_message(except.what());
			}

		}catch(cppdb::cppdb_error & except)
		{
			utility_exit_with_message(except.what());
		}
		cycle++;
	}
}

vector1<boost::uuids::uuid>
struct_ids_from_tag(
	sessionOP db_session,
	string const & tag
) {
	string statement_string = "SELECT struct_id FROM structures WHERE tag=?;";
	statement stmt = safely_prepare_statement(statement_string,db_session);
	stmt.bind(1,tag);
	result res = stmt.query();

	vector1<boost::uuids::uuid> uuids;
	while(res.next()){
		boost::uuids::uuid uuid;
		res >> uuid;
		uuids.push_back(uuid);
	}
	return uuids;
}

bool
table_exists(
	sessionOP db_session,
	string const & table_name
) {

	// TODO: handle when the current database is not the one from the
	// option system, can someone with mysql try this and see if it
	// works?
	// "SHOW TABLES IN database();"
	string statement_string;
	statement stmt;
	Size i(1);
	switch(db_session->get_db_mode()){
	case DatabaseMode::sqlite3:
		statement_string = "SELECT name FROM sqlite_master WHERE name=?;";
		stmt = safely_prepare_statement(statement_string,db_session);
		break;
	case DatabaseMode::mysql:
		statement_string = "SHOW TABLES WHERE Tables_in_"+db_session->get_db_name()+" = ?;";
		stmt = safely_prepare_statement(statement_string,db_session);
		break;
	case DatabaseMode::postgres:
		if(db_session->get_pq_schema() == ""){
			statement_string =
				"SELECT tablename \n"
				"FROM pg_catalog.pg_tables \n"
				"WHERE tablename = ?;";
			stmt = safely_prepare_statement(statement_string, db_session);
		} else {
			statement_string =
				"SELECT tablename FROM pg_catalog.pg_tables \n"
				"WHERE schemaname = ? AND tablename = ?;";
			stmt = safely_prepare_statement(statement_string, db_session);
			stmt.bind(i, db_session->get_pq_schema());
			i++;
		}
		break;
	default:
		utility_exit_with_message("unknown database mode");
	}

	stmt.bind(i, table_name);
	result res(stmt.query());

	return res.next();
}

//Simply (probably overly so) protection from SQL injection
void
check_statement_sanity(
	string sql
) {
	if (!boost::istarts_with(sql, "SELECT")){
		utility_exit_with_message("ERROR: Database select statement is not safe! Only SELECT statements are allowed.");
	}

	int semicolon_count=0;
	for (size_t i = 0; i < sql.size(); i++){
		if (sql[i] == ';'){
			semicolon_count++;
		}
	}
	if(semicolon_count > 1){
		utility_exit_with_message("Database select statement is not safe! Only 1 SQL statement is allowed");
	}
}

//This should ideally only be used for reference tables that have static data that needs to only be written once(ex: dssp_codes)
void
insert_or_ignore(
	string table_name,
	std::vector<string> column_names,
	std::vector<string> values,
	sessionOP db_session
){

	string statement_string="";
	switch(db_session->get_db_mode()){
		case utility::sql_database::DatabaseMode::mysql:{
			statement_string = "INSERT IGNORE into "+table_name+"(";
			for(size_t i=0; i<column_names.size(); i++){
				statement_string+=column_names[i];
				if(i != column_names.size()-1){
					statement_string+=",";
				}
			}

			statement_string+=") VALUES(";
			for(size_t i=0; i<values.size(); i++){
				statement_string+=values[i];
				if(i != column_names.size()-1){
					statement_string+=",";
				}
			}
			statement_string+=");";

			statement stmt = (*db_session) << statement_string;
			safely_write_to_database(stmt);
			break;
		}
		case utility::sql_database::DatabaseMode::postgres:{
			//This is a dirty postgres hack and seems to be the easiest workaround for lack of INSERT IGNORE support in postgres
			string select_statement_string = "SELECT * FROM "+table_name+" WHERE ";
			for(size_t i=0; i<column_names.size(); i++){
				select_statement_string+=column_names[i] + "=" + values[i];
				if(i != column_names.size()-1){
					select_statement_string+=" AND ";
				}
			}
			select_statement_string+=";";

			statement select_stmt = (*db_session) << select_statement_string;
			result res = safely_read_from_database(select_stmt);

			if(!res.next()){
				statement_string += "INSERT into "+table_name+"(";
				for(size_t i=0; i<column_names.size(); i++){
					statement_string+=column_names[i];
					if(i != column_names.size()-1){
						statement_string+=",";
					}
				}

				statement_string+=") VALUES(";
				for(size_t i=0; i<values.size(); i++){
					statement_string+=values[i];
					if(i != column_names.size()-1){
						statement_string+=",";
					}
				}
				statement_string+=");";
				statement stmt = (*db_session) << statement_string;
				safely_write_to_database(stmt);
			}
			break;
		}
		case utility::sql_database::DatabaseMode::sqlite3:{
			statement_string = "INSERT OR IGNORE into "+table_name+"(";
			for(size_t i=0; i<column_names.size(); i++){
				statement_string+=column_names[i];
				if(i != column_names.size()-1){
				 statement_string+=",";
				}
			}

			statement_string+=") VALUES(";
			for(size_t i=0; i<values.size(); i++){
				statement_string+=values[i];
				if(i != column_names.size()-1){
					statement_string+=",";
				}
			}
			statement_string+=");";

			statement stmt = (*db_session) << statement_string;
			safely_write_to_database(stmt);
			break;
		}
		default:
			utility_exit_with_message(
				"Unrecognized database mode: '" +
				name_from_database_mode(db_session->get_db_mode()) + "'");
	}
}

void write_schema_to_database(
	string schema_str,
	sessionOP db_session)
{
	boost::char_separator< char > sep(";");
	boost::tokenizer< boost::char_separator< char > > tokens( schema_str, sep );
	foreach( std::string const & stmt_str, tokens){
		string trimmed_stmt_str(utility::trim(stmt_str, " \n\t"));
		if(trimmed_stmt_str.size()){
			try{
				statement stmt = (*db_session) << trimmed_stmt_str + ";";
				safely_write_to_database(stmt);
			} catch (cppdb_error e) {
				TR.Error
					<< "ERROR reading schema \n"
					<< trimmed_stmt_str << std::endl;
				TR.Error << e.what() << std::endl;
				utility_exit();
			}
		}
	}
}


void
set_cache_size(
	sessionOP db_session,
	Size cache_size
) {

	if(db_session->get_db_mode() == DatabaseMode::sqlite3){
		stringstream stmt_ss;
		stmt_ss << "PRAGMA cache_size = " << cache_size << ";";
		statement stmt(safely_prepare_statement(stmt_ss.str(), db_session));
		safely_write_to_database(stmt);
	} else {
		TR
			<< "WARNING: Attempting to set database cache size "
			<< "for a database type for which this is currently not supported: "
			<< "'" << name_from_database_mode(db_session->get_db_mode()) << "'." << std::endl;
	}
}

std::string make_compound_statement(
	std::string const & table_name,
	std::vector<std::string> const & column_names,
	platform::Size const & row_count)
{
	std::string table_definition = table_name + " (" + utility::join(column_names,",") + ")";
	std::string value_list;
	platform::Size column_count = column_names.size();
	for(platform::Size i = 0; i < row_count;++i)
	{
		std::string row_block= "(?";
		for(platform::Size j = 1; j < column_count; ++j)
		{
			row_block += ",?";
		}
		row_block += ")";

		value_list += row_block;
		if(i != row_count-1)
		{
			value_list += ", ";
		} }

	return "INSERT INTO "+table_definition+ " VALUES " + value_list +";";
}

///@detail build database connection from options in a tag, this is useful make sure the fields for
///constructing a database connection are consistent across different tags.
utility::sql_database::sessionOP
parse_database_connection(
	utility::tag::TagPtr const tag
){
	using std::endl;
	using namespace basic::options;
	using namespace basic::options::OptionKeys::inout;
	using utility::sql_database::DatabaseSessionManager;
	using namespace basic::resource_manager;

	if(tag->hasOption("database_resource")){
		std::string database_resource = tag->getOption<string>("database_resource");
		if ( ! ResourceManager::get_instance()->has_resource_with_description( database_resource ) )
		{
			throw utility::excn::EXCN_Msg_Exception
				( "You specified a database_resource of '" + database_resource +
					"', but the ResourceManager doesn't have a resource with that description." );
		}
		return get_resource< utility::sql_database::session >( database_resource );
	}

	if(tag->hasOption("database_resource_tag")){
		std::string database_resource_tag = tag->getOption<string>(
			"database_resource_tag");
		if ( ! ResourceManager::get_instance()->has_resource(
				database_resource_tag ) )
		{
			throw utility::excn::EXCN_Msg_Exception
				( "You specified a database_resource_tag of '" + database_resource_tag +
					"', but the ResourceManager doesn't have a resource with that tag." );
		}
		utility::sql_database::session * db_session(dynamic_cast< utility::sql_database::session * > (
				ResourceManager::get_instance()->find_resource(database_resource_tag)()));
		if(!db_session){
			stringstream err_msg;
			err_msg
				<< "You specified a database_resource_tag of '" + database_resource_tag + "', while the ResourceManager does have a resource with that tag, it couldn't cast into a database session.";
			throw utility::excn::EXCN_Msg_Exception(err_msg.str());
		}
		return db_session;
	}

	utility::sql_database::TransactionMode::e transaction_mode;
	if(tag->hasOption("transaction_mode")){
		transaction_mode = utility::sql_database::transaction_mode_from_name(
			tag->getOption<string>("transaction_mode"));
	} else {
		transaction_mode = utility::sql_database::TransactionMode::standard;
	}

	Size chunk_size;
	switch(transaction_mode){
		case(utility::sql_database::TransactionMode::none):
			if(tag->hasOption("chunk_size")){
				TR << "WARNING: You must specify 'transaction_mode=chunk' ";
				TR << "to use the 'chunk size' tag." << endl;
			}
			chunk_size=0;
			break;
		case(utility::sql_database::TransactionMode::standard):
			if(tag->hasOption("chunk_size")){
				TR << "WARNING: You must specify 'transaction_mode=chunk' ";
				TR << "to use the 'chunk size' tag." << endl;
			}
			chunk_size=0;
			break;
		case(utility::sql_database::TransactionMode::chunk):
			if(!tag->hasOption("chunk_size")){
				utility_exit_with_message(
					"Must specify chunk_size if using the chunk transaction mode");
			}
			chunk_size=
				tag->getOption<Size>("chunk_size");
			break;
		default:
			utility_exit_with_message(
				"Unrecognized transaction mode: '" +
				name_from_transaction_mode(transaction_mode) + "'");
	}

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
				database_mode, transaction_mode, chunk_size,
				database_name, "", "", "", "", 0,
				tag->getOption("database_read_only", false),
				tag->getOption("database_separate_db_per_mpi_process", false));

		case utility::sql_database::DatabaseMode::mysql:
		case utility::sql_database::DatabaseMode::postgres:{

			std::string database_host;
			if(!tag->hasOption("database_host")){
				if(!option[dbms::host].user()){
					utility_exit_with_message(
						"WARNING: To connect to a postgres or mysql database you must set"
						" the database_host tag or specify -dbms:host on the command line.");
				}
				else{
					database_host=option[dbms::host];
				}
			}
			else{
				database_host=tag->getOption<string>("database_host");
			}

			std::string database_user;
			if(!tag->hasOption("database_user")){
				if(!option[dbms::user].user()){
					utility_exit_with_message(
						"WARNING: To connect to a postgres or mysql database you must set"
						"the database_user tag or specify -dbms:user on the command line.");
				}
				else{
					database_user=option[dbms::user];
				}
			}
			else{
				database_user=tag->getOption<string>("database_user");
			}

			std::string database_password;
			if(!tag->hasOption("database_password")){
				if(!option[dbms::password].user()){
					utility_exit_with_message(
						"WARNING: To connect to a postgres or mysql database you must set"
						"the database_password tag or specify -dbms:password on the command line.");
				}
				else{
					database_password=option[dbms::password];
				}
			}
			else{
				database_password=tag->getOption<string>("database_password");
			}

			Size database_port;
			if(!tag->hasOption("database_port")){
				if(!option[dbms::port].user()){
					utility_exit_with_message(
						"WARNING: To connect to a postgres or mysql database you must set"
						"the database_port tag or specify -dbms:port on the command line.");
				}
				else{
					database_port=option[dbms::port];
				}
			}
			else{
				database_port=tag->getOption<Size>("database_port");
			}

			return DatabaseSessionManager::get_instance()->get_db_session(
				database_mode, transaction_mode, chunk_size,
				database_name, database_pq_schema,
				database_host,
				database_user,
				database_password,
				database_port);
		}
		default:
			utility_exit_with_message(
				"Unrecognized database mode: '" +
				name_from_database_mode(database_mode) + "'");
	}
	return 0;
}

}
}
