// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ScoreTypeFeatures.cc
/// @brief  report protocol level features to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ScoreTypeFeatures.hh>
#include <protocols/features/util.hh>

//External
#include <boost/uuid/uuid.hpp>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>

// Platform Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>
#include <string>


namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using core::pose::Pose;
using core::graph::Graph;
using core::scoring::Energies;
using core::scoring::EnergyGraph;
using core::scoring::EnergyMap;
using core::scoring::EnergyEdge;
using core::scoring::getScoreFunction;
using core::scoring::hbond_sr_bb;
using core::scoring::hbond_lr_bb;
using core::scoring::n_score_types;
using core::scoring::ScoreTypeManager;
using core::scoring::ScoreType;
using core::scoring::ScoreFunctionOP;
using core::scoring::hbonds::HBondSet;
using core::scoring::hbonds::get_hbond_energies;
using core::scoring::EnergiesCacheableDataType::HBOND_SET;
using core::Size;
using core::Real;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

ScoreTypeFeatures::ScoreTypeFeatures() :
	scfxn_(getScoreFunction())
{}

ScoreTypeFeatures::ScoreTypeFeatures(
	ScoreFunctionOP scfxn
) :
	scfxn_(scfxn)
{
	if ( scfxn_ == 0 ) {
		utility_exit_with_message( "ScoreTypeFeatures may not be constructed with a null-pointer ScoreFunctionOP" );
	}
}

ScoreTypeFeatures::ScoreTypeFeatures(
	ScoreTypeFeatures const & src
) :
	FeaturesReporter(),
	scfxn_(src.scfxn_)
{}

ScoreTypeFeatures::~ScoreTypeFeatures() {}

string
ScoreTypeFeatures::type_name() const { return "ScoreTypeFeatures"; }

void
ScoreTypeFeatures::write_schema_to_db(
	sessionOP db_session
) const {

	using namespace basic::database::schema_generator;

	//******score_types******//
	Column batch_id("batch_id", new DbInteger(), false);
	Column score_type_id("score_type_id", new DbInteger(), false);
	Column score_type_name("score_type_name", new DbText(), false);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(batch_id);
	pkey_cols.push_back(score_type_id);

	Schema score_types("score_types", PrimaryKey(pkey_cols));
	score_types.add_column(batch_id);
	score_types.add_column(score_type_id);
	score_types.add_column(score_type_name);

	score_types.add_foreign_key(ForeignKey(batch_id, "batches", "batch_id", true));

	score_types.write(db_session);
}

utility::vector1<std::string>
ScoreTypeFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ProtocolFeatures");
	return dependencies;
}

Size
ScoreTypeFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	boost::uuids::uuid const struct_id,
	sessionOP db_session
){
	Size const batch_id(get_batch_id(struct_id, db_session));

	insert_score_type_rows(batch_id, db_session);
	return 0;
}


Size
ScoreTypeFeatures::report_features(
	Size const batch_id,
	sessionOP db_session
){
	insert_score_type_rows(batch_id, db_session);
	return 0;
}




void ScoreTypeFeatures::delete_record(
	boost::uuids::uuid,
	utility::sql_database::sessionOP ){}

void
ScoreTypeFeatures::insert_score_type_rows(
	Size protocol_id,
	sessionOP db_session
) {
	std::string statement_string;

	switch(db_session->get_db_mode()){
	case utility::sql_database::DatabaseMode::sqlite3:
		statement_string = "INSERT OR IGNORE INTO score_types (batch_id, score_type_id, score_type_name) VALUES (?,?,?);";
		break;
	case utility::sql_database::DatabaseMode::mysql:
	case utility::sql_database::DatabaseMode::postgres:
		statement_string = "INSERT IGNORE INTO score_types (batch_id, score_type_id, score_type_name) VALUES (?,?,?);";
		break;
	default:
		utility_exit_with_message(
			"Unrecognized database mode: '" +
			name_from_database_mode(db_session->get_db_mode()) + "'");
	}

	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));

	for(Size score_type_id=1; score_type_id <= n_score_types; ++score_type_id){
		ScoreType type(static_cast<ScoreType>(score_type_id));

		string const score_type( ScoreTypeManager::name_from_score_type(type) );
		stmt.bind(1,protocol_id);
		stmt.bind(2,score_type_id);
		stmt.bind(3,score_type);
		basic::database::safely_write_to_database(stmt);


	}
	//transact_guard.commit();
}


} // namesapce
} // namespace
