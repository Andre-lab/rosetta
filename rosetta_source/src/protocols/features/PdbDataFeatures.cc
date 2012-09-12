// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/features/PdbDataFeatures.cc
/// @author Sam DeLuca

#include <protocols/features/PdbDataFeatures.hh>

//project headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

//External
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>

//Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <algorithm>
#include <limits>

namespace protocols {
namespace features {

using std::string;
using std::max;
using std::min;
using std::numeric_limits;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pose::Pose;
using utility::sql_database::sessionOP;
using utility::vector1;
using basic::database::safely_read_from_database;
using basic::database::safely_prepare_statement;
using basic::database::safely_write_to_database;
using basic::database::table_exists;
using cppdb::result;
using cppdb::statement;
using core::pose::PDBInfo;
using core::pose::PDBInfoOP;
using core::pose::PDBInfoCOP;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowData;
using basic::database::insert_statement_generator::RowDataBaseOP;


PdbDataFeatures::PdbDataFeatures()
{

}

PdbDataFeatures::PdbDataFeatures(PdbDataFeatures const & )
{

}

PdbDataFeatures::~PdbDataFeatures()
{

}

string PdbDataFeatures::type_name() const
{
	return "PdbDataFeatures";
}

void
PdbDataFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", new DbUUID(), false);
	Column residue_number("residue_number", new DbInteger(), false);

	utility::vector1<Column> pkey_cols;
	pkey_cols.push_back(struct_id);
	pkey_cols.push_back(residue_number);

	//******residue_pdb_identification******//
	Column chain_id("chain_id", new DbText(), false);
	Column insertion_code("insertion_code", new DbText(), false);
	Column pdb_residue_number("pdb_residue_number", new DbInteger(), false);

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	Schema residue_pdb_identification("residue_pdb_identification", PrimaryKey(pkey_cols));
	residue_pdb_identification.add_column(struct_id);
	residue_pdb_identification.add_column(residue_number);
	residue_pdb_identification.add_column(chain_id);
	residue_pdb_identification.add_column(insertion_code);
	residue_pdb_identification.add_column(pdb_residue_number);

	residue_pdb_identification.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true));

	residue_pdb_identification.write(db_session);


	//******residue_pdb_confidence******//
	Column max_temperature("max_temperature", new DbReal(), false);
	Column max_bb_temperature("max_bb_temperature", new DbReal(), false);
	Column max_sc_temperature("max_sc_temperature", new DbReal(), false);
	Column min_occupancy("min_occupancy", new DbReal(), false);
	Column min_bb_occupancy("min_bb_occupancy", new DbReal(), false);
	Column min_sc_occupancy("min_sc_occupancy", new DbReal(), false);

	utility::vector1<Column> pdb_ident_pkeys;
	pdb_ident_pkeys.push_back(struct_id);
	pdb_ident_pkeys.push_back(residue_number);

	Schema residue_pdb_confidence("residue_pdb_confidence", PrimaryKey(pkey_cols));
	residue_pdb_confidence.add_column(struct_id);
	residue_pdb_confidence.add_column(residue_number);
	residue_pdb_confidence.add_column(max_temperature);
	residue_pdb_confidence.add_column(max_bb_temperature);
	residue_pdb_confidence.add_column(max_sc_temperature);
	residue_pdb_confidence.add_column(min_occupancy);
	residue_pdb_confidence.add_column(min_bb_occupancy);
	residue_pdb_confidence.add_column(min_sc_occupancy);

	residue_pdb_confidence.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true));

	residue_pdb_confidence.write(db_session);

}

utility::vector1<std::string>
PdbDataFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}


Size PdbDataFeatures::report_features(
	Pose const & pose,
	vector1<bool> const &,
	boost::uuids::uuid struct_id,
	sessionOP db_session )
{
	insert_residue_pdb_identification_rows(struct_id,db_session,pose);
	insert_residue_pdb_confidence_rows(struct_id, db_session, pose);
	return 0;
}

void PdbDataFeatures::delete_record(
	boost::uuids::uuid struct_id,
	sessionOP db_session)
{
	string id_statement_string = "DELETE FROM residue_pdb_identification WHERE struct_id = ?;\n";
	statement id_stmt(safely_prepare_statement(id_statement_string,db_session));
	id_stmt.bind(1,struct_id);
	safely_write_to_database(id_stmt);

	string confidence_statement_string = "DELETE FROM residue_pdb_confidence WHERE struct_id = ?;\n";
	statement confidence_stmt(
		safely_prepare_statement(confidence_statement_string, db_session));
	confidence_stmt.bind(1,struct_id);
	safely_write_to_database(confidence_stmt);
}

void PdbDataFeatures::load_into_pose(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose)
{
	load_residue_pdb_identification(db_session, struct_id, pose);
	load_residue_pdb_confidence(db_session, struct_id, pose);
}


void PdbDataFeatures::load_residue_pdb_identification(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose)
{
	if(!table_exists(db_session, "residue_pdb_identification")) return;

	vector1<int> pdb_numbers;
	vector1<char> pdb_chains;
	vector1<char> insertion_codes;
	string statement_string =
		"SELECT\n"
				"	r_id.struct_id,\n"
		"	r_id.residue_number,\n"
		"	r_id.chain_id,\n"
		"	r_id.insertion_code,\n"
		"	r_id.pdb_residue_number\n"
		"FROM\n"
		"	residue_pdb_identification AS r_id\n"
		"WHERE\n"
		"	r_id.struct_id=?;";

	statement stmt(safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	result res(safely_read_from_database(stmt));

	while(res.next()) {
				boost::uuids::uuid temp;
		Size residue_number;
		//cppdb doesn't do char's
		string chain_id;
		string insertion_code;
		int pdb_residue_number;

		res >> temp >> residue_number >> chain_id >> insertion_code >> pdb_residue_number;

		pdb_chains.push_back(chain_id[0]);
		insertion_codes.push_back(insertion_code[0]);
				pdb_numbers.push_back(pdb_residue_number);
	}

	if(!pose.pdb_info()){
		pose.pdb_info(new PDBInfo(pose.total_residue()));
	}

	pose.pdb_info()->set_numbering(pdb_numbers);
	pose.pdb_info()->set_chains(pdb_chains);
	pose.pdb_info()->set_icodes(insertion_codes);
}

void PdbDataFeatures::insert_residue_pdb_identification_rows(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Pose const & pose)
{

	InsertGenerator pdb_ident_insert("residue_pdb_identification");
	pdb_ident_insert.add_column("struct_id");
	pdb_ident_insert.add_column("residue_number");
	pdb_ident_insert.add_column("chain_id");
	pdb_ident_insert.add_column("insertion_code");
	pdb_ident_insert.add_column("pdb_residue_number");

	RowDataBaseOP struct_id_data = new RowData<boost::uuids::uuid>("struct_id",struct_id);

	Size res_num(pose.n_residue());
	for(Size index = 1; index <= res_num; ++index)
	{
		string chain_id(& pose.pdb_info()->chain(index),1);
		string insertion_code(&pose.pdb_info()->icode(index),1);
		int pdb_residue_number = pose.pdb_info()->number(index);

		RowDataBaseOP index_data = new RowData<Size>("residue_number",index);
		RowDataBaseOP chain_id_data = new RowData<string>("chain_id",chain_id);
		RowDataBaseOP insertion_code_data = new RowData<string>("insertion_code",insertion_code);
		RowDataBaseOP pdb_resnum_data = new RowData<int>("pdb_residue_number",pdb_residue_number	);
		pdb_ident_insert.add_row(
			utility::tools::make_vector(struct_id_data,index_data,chain_id_data,insertion_code_data,pdb_resnum_data));


	}
	pdb_ident_insert.write_to_database(db_session);
}


void PdbDataFeatures::load_residue_pdb_confidence(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose)
{
	if(!table_exists(db_session, "residue_pdb_confidence")) return;

	if(!pose.pdb_info()){
		pose.pdb_info(new PDBInfo(pose.total_residue()));
	}

	std::string statement_string =
		"SELECT\n"
		"	r_conf.residue_number,\n"
		"	r_conf.max_bb_temperature,\n"
		"	r_conf.max_sc_temperature,\n"
		"	r_conf.min_bb_occupancy,\n"
		"	r_conf.min_sc_occupancy\n"
		"FROM\n"
		"	residue_pdb_confidence AS r_conf\n"
		"WHERE\n"
		"	r_conf.struct_id=?;";

	statement stmt(safely_prepare_statement(statement_string,db_session));
	stmt.bind(1,struct_id);
	result res(safely_read_from_database(stmt));

	while(res.next()) {
		Size residue_number;
		Real max_bb_temperature;
		Real max_sc_temperature;
		Real min_bb_occupancy;
		Real min_sc_occupancy;

		res
			>> residue_number
			>> max_bb_temperature
			>> max_sc_temperature
			>> min_bb_occupancy
			>> min_sc_occupancy;

		Residue const & residue(pose.residue(residue_number));

		pose.pdb_info()->resize_atom_records(
			residue_number, residue.nheavyatoms(), false);

		for(
			Size atom_index=1;
			atom_index <= residue.last_backbone_atom();
			++atom_index){
			pose.pdb_info()->temperature(
				residue_number, atom_index, max_bb_temperature);
			pose.pdb_info()->occupancy(
				residue_number, atom_index, min_bb_occupancy);
		}
		for(
			Size atom_index = residue.first_sidechain_atom();
			atom_index <= residue.nheavyatoms();
			++atom_index){
			pose.pdb_info()->temperature(
				residue_number, atom_index, max_sc_temperature);
			pose.pdb_info()->occupancy(
				residue_number, atom_index, min_sc_occupancy);
		}
	}
}


void PdbDataFeatures::insert_residue_pdb_confidence_rows(
	boost::uuids::uuid struct_id,
	sessionOP db_session,
	Pose const & pose)
{
	PDBInfoCOP pdb_info(pose.pdb_info());
	if(!pdb_info) return;


	InsertGenerator confidence_insert("residue_pdb_confidence");
	confidence_insert.add_column("struct_id");
	confidence_insert.add_column("residue_number");
	confidence_insert.add_column("max_temperature");
	confidence_insert.add_column("max_bb_temperature");
	confidence_insert.add_column("max_sc_temperature");
	confidence_insert.add_column("min_occupancy");
	confidence_insert.add_column("min_bb_occupancy");
	confidence_insert.add_column("min_sc_occupancy");


	RowDataBaseOP struct_id_data = new RowData<boost::uuids::uuid>("struct_id",struct_id);

	for(Size ri=1; ri <= pose.n_residue(); ++ri) {
		Residue const & r(pose.residue(ri));
		Real max_bb_temperature(-1), max_sc_temperature(-1);
		Real min_bb_occupancy(9999999);
		Real min_sc_occupancy(9999999);
		Size const n_bb(r.n_mainchain_atoms());
		Size const n_sc(r.nheavyatoms() - r.n_mainchain_atoms());
		for(Size ai=1; ai <= n_bb; ++ai){
			max_bb_temperature = max(max_bb_temperature, pdb_info->temperature(ri, ai));
			min_bb_occupancy = min(min_bb_occupancy, pdb_info->occupancy(ri, ai));
		}
		for(Size ai=1; ai <= n_sc; ++ai){
			max_sc_temperature = max(max_sc_temperature, pdb_info->temperature(ri, ai));
			min_sc_occupancy = min(min_sc_occupancy, pdb_info->occupancy(ri, ai));
		}
		Real const max_temperature = max(max_bb_temperature, max_sc_temperature);
		Real const min_occupancy = min(min_bb_occupancy, min_sc_occupancy);

		RowDataBaseOP residue_number_data = new RowData<Size>("residue_number",ri);
		RowDataBaseOP max_temp_data = new RowData<Real>("max_temperature",max_temperature);
		RowDataBaseOP max_bb_temp_data = new RowData<Real>("max_bb_temperature",max_bb_temperature);
		RowDataBaseOP max_sc_temp_data = new RowData<Real>("max_sc_temperature",max_sc_temperature);
		RowDataBaseOP min_occupancy_data  = new RowData<Real>("min_occupancy",min_occupancy);
		RowDataBaseOP min_bb_occupancy_data = new RowData<Real>("min_bb_occupancy",min_bb_occupancy);
		RowDataBaseOP min_sc_occupancy_data = new RowData<Real>("min_sc_occupancy",min_sc_occupancy);

		confidence_insert.add_row(
			utility::tools::make_vector(struct_id_data,residue_number_data,max_temp_data,max_bb_temp_data,max_sc_temp_data,
				min_occupancy_data,min_bb_occupancy_data,min_sc_occupancy_data)
			);

	}

	confidence_insert.write_to_database(db_session);

}




}
}
