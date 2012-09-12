// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ProteinResidueConformationFeatures.cc
/// @brief  report idealized torsional DOFs Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/ProteinResidueConformationFeatures.hh>

//External
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

// Project Headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/string_util.hh>
#include <utility/tools/make_vector.hh>

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/insert_statement_generator/InsertGenerator.hh>
#include <basic/database/insert_statement_generator/RowData.hh>


// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <cmath>
#include <sstream>

namespace protocols{
namespace features{

static basic::Tracer TR("protocols.features.ProteinResidueConformationFeatures");

using std::string;
using core::Size;
using core::Real;
using core::Vector;
using core::chemical::num_canonical_aas;
using core::conformation::Residue;
using core::pose::Pose;
using core::id::AtomID;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using cppdb::result;
using basic::database::insert_statement_generator::InsertGenerator;
using basic::database::insert_statement_generator::RowDataBaseOP;
using basic::database::insert_statement_generator::RowData;

string
ProteinResidueConformationFeatures::type_name() const {
	return "ProteinResidueConformationFeatures";
}

void
ProteinResidueConformationFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{

	using namespace basic::database::schema_generator;

	//******protein_residue_conformation******//
	Column struct_id("struct_id", new DbUUID(), false);
	Column seqpos("seqpos", new DbInteger(), false);
	Column secstruct("secstruct", new DbText(), false);
	Column phi("phi", new DbDouble(), false);
	Column psi("psi", new DbDouble(), false);
	Column omega("omega", new DbDouble(), false);
	Column chi1("chi1", new DbDouble(), false);
	Column chi2("chi2", new DbDouble(), false);
	Column chi3("chi3", new DbDouble(), false);
	Column chi4("chi4", new DbDouble(), false);


	utility::vector1<Column> prot_res_pkeys;
	prot_res_pkeys.push_back(struct_id);
	prot_res_pkeys.push_back(seqpos);

	utility::vector1<Column> fkey_cols;
	fkey_cols.push_back(struct_id);
	fkey_cols.push_back(seqpos);

	utility::vector1<std::string> fkey_reference_cols;
	fkey_reference_cols.push_back("struct_id");
	fkey_reference_cols.push_back("resNum");

	Schema protein_residue_conformation("protein_residue_conformation", PrimaryKey(prot_res_pkeys));
	protein_residue_conformation.add_column(struct_id);
	protein_residue_conformation.add_column(seqpos);
	protein_residue_conformation.add_column(secstruct);
	protein_residue_conformation.add_column(phi);
	protein_residue_conformation.add_column(psi);
	protein_residue_conformation.add_column(omega);
	protein_residue_conformation.add_column(chi1);
	protein_residue_conformation.add_column(chi2);
	protein_residue_conformation.add_column(chi3);
	protein_residue_conformation.add_column(chi4);
	protein_residue_conformation.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

	protein_residue_conformation.write(db_session);

	//******residue_atom_coords******//
	Column atomno("atomno", new DbInteger(), false);
	Column x("x", new DbDouble(), false);
	Column y("y", new DbDouble(), false);
	Column z("z", new DbDouble(), false);

	utility::vector1<Column> res_atm_coords_pkeys;
	res_atm_coords_pkeys.push_back(struct_id);
	res_atm_coords_pkeys.push_back(seqpos);
	res_atm_coords_pkeys.push_back(atomno);

	Schema residue_atom_coords("residue_atom_coords", PrimaryKey(res_atm_coords_pkeys));
	residue_atom_coords.add_column(struct_id);
	residue_atom_coords.add_column(seqpos);
	residue_atom_coords.add_column(atomno);
	residue_atom_coords.add_column(x);
	residue_atom_coords.add_column(y);
	residue_atom_coords.add_column(z);
	residue_atom_coords.add_foreign_key(ForeignKey(fkey_cols, "residues", fkey_reference_cols, true));

	residue_atom_coords.write(db_session);

}

utility::vector1<std::string>
ProteinResidueConformationFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}


Size
ProteinResidueConformationFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	boost::uuids::uuid struct_id,
	sessionOP db_session
){
	bool fullatom(pose.is_fullatom());

	//check to see if this structure is ideal
	bool ideal = true;
	if(!basic::options::option[basic::options::OptionKeys::out::file::force_nonideal_structure]())
	{
		core::conformation::Conformation const & conformation(pose.conformation());
		for(core::Size resn=1; resn <= pose.n_residue();++resn){
			if(!relevant_residues[resn]) continue;
			bool residue_status(core::conformation::is_ideal_position(resn,conformation));
			if(!residue_status){
				ideal = false;
				break;
			}
		}
	}else
	{
		ideal = false;
	}

	InsertGenerator conformation_insert("protein_residue_conformation");
	conformation_insert.add_column("struct_id");
	conformation_insert.add_column("seqpos");
	conformation_insert.add_column("secstruct");
	conformation_insert.add_column("phi");
	conformation_insert.add_column("psi");
	conformation_insert.add_column("omega");
	conformation_insert.add_column("chi1");
	conformation_insert.add_column("chi2");
	conformation_insert.add_column("chi3");
	conformation_insert.add_column("chi4");

	InsertGenerator atom_insert("residue_atom_coords");
	atom_insert.add_column("struct_id");
	atom_insert.add_column("seqpos");
	atom_insert.add_column("atomno");
	atom_insert.add_column("x");
	atom_insert.add_column("y");
	atom_insert.add_column("z");

	RowDataBaseOP struct_id_data = new RowData<boost::uuids::uuid>("struct_id",struct_id);

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		if(!relevant_residues[i]) continue;
		Residue const & resi = pose.residue(i);
		if(resi.aa() > num_canonical_aas)
		{
			continue;
		}
		std::string secstruct = utility::to_string<char>(pose.secstruct(i));

		Real phi = 0.0;
		Real psi = 0.0;
		Real omega = 0.0;

		//If you have a non ideal structure, and you store both cartesian coordinates and backbone
		//chi angles and read them into a pose, most of the backbone oxygens will be placed in correctly
		//I currently have absolutely no idea why this is, but this fixes it.
		//It is worth noting that the current implementation of Binary protein silent files does the same thing

		if(ideal)
		{
			 phi = resi.mainchain_torsion(1);
			 psi = resi.mainchain_torsion(2);
			 omega = resi.mainchain_torsion(3);
		}
		Real chi1 = fullatom && resi.nchi() >= 1 ? resi.chi(1) : 0.0;
		Real chi2 = fullatom && resi.nchi() >= 2 ? resi.chi(2) : 0.0;
		Real chi3 = fullatom && resi.nchi() >= 3 ? resi.chi(3) : 0.0;
		Real chi4 = fullatom && resi.nchi() >= 4 ? resi.chi(4) : 0.0;

		RowDataBaseOP seqpos_data = new RowData<Size>("seqpos",i);
		RowDataBaseOP secstruct_data = new RowData<std::string>("secstruct",secstruct);
		RowDataBaseOP phi_data = new RowData<Real>("phi",phi);
		RowDataBaseOP psi_data = new RowData<Real>("psi",psi);
		RowDataBaseOP omega_data = new RowData<Real>("omega",omega);
		RowDataBaseOP chi1_data = new RowData<Real>("chi1",chi1);
		RowDataBaseOP chi2_data = new RowData<Real>("chi2",chi2);
		RowDataBaseOP chi3_data = new RowData<Real>("chi3",chi3);
		RowDataBaseOP chi4_data = new RowData<Real>("chi4",chi4);


		conformation_insert.add_row(utility::tools::make_vector(
			struct_id_data,seqpos_data,secstruct_data,phi_data,
			psi_data,omega_data,chi1_data,chi2_data,chi3_data,chi4_data));

		if(!ideal)
		{
			for(Size atom = 1; atom <= resi.natoms(); ++atom)
			{
				core::Vector coords = resi.xyz(atom);

				RowDataBaseOP atomno_data = new RowData<Size>("atomno",atom);
				RowDataBaseOP x_data = new RowData<Real>("x",coords.x());
				RowDataBaseOP y_data = new RowData<Real>("y",coords.y());
				RowDataBaseOP z_data = new RowData<Real>("z",coords.z());

				atom_insert.add_row(utility::tools::make_vector(
					struct_id_data,seqpos_data,atomno_data,x_data,y_data,z_data));

			}
		}
	}

	conformation_insert.write_to_database(db_session);
	atom_insert.write_to_database(db_session);

	return 0;
}

void
ProteinResidueConformationFeatures::delete_record(
	boost::uuids::uuid struct_id,
	utility::sql_database::sessionOP db_session
){

	//std::string struct_id_string(to_string(struct_id));
	statement conf_stmt(basic::database::safely_prepare_statement("DELETE FROM protein_residue_conformation WHERE struct_id = ?;\n",db_session));
	conf_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(conf_stmt);
	statement atom_stmt(basic::database::safely_prepare_statement("DELETE FROM residue_atom_coords WHERE struct_id = ?;\n",db_session));
	atom_stmt.bind(1,struct_id);
	basic::database::safely_write_to_database(atom_stmt);

}


void
ProteinResidueConformationFeatures::load_into_pose(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose
){
	load_conformation(db_session, struct_id, pose);
}

void
ProteinResidueConformationFeatures::load_conformation(
	sessionOP db_session,
	boost::uuids::uuid struct_id,
	Pose & pose
){

	if(!basic::database::table_exists(db_session, "protein_residue_conformation")){
		TR << "WARNING: protein_residue_conformation table does not exist and thus respective data will not be added to the pose!" << std::endl;
		return;
	}

	set_coords_for_residues(db_session,struct_id,pose);


	if(pose.is_fullatom()){
		std::string statement_string =
			"SELECT\n"
			"	seqpos,\n"
			"	secstruct,\n"
			"	phi,\n"
			"	psi,\n"
			"	omega,\n"
			"	chi1,\n"
			"	chi2,\n"
			"	chi3,\n"
			"	chi4\n"
			"FROM\n"
			"	protein_residue_conformation\n"
			"WHERE\n"
			"	protein_residue_conformation.struct_id=?;";
		statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
		//std::string struct_id_string(to_string(struct_id));
		stmt.bind(1,struct_id);

		result res(basic::database::safely_read_from_database(stmt));

		while(res.next()){
			Size seqpos;
			Real phi,psi,omega,chi1,chi2,chi3,chi4;
			std::string secstruct;
			res >> seqpos >> secstruct >> phi >> psi >> omega >> chi1 >> chi2 >> chi3 >> chi4 ;
			if (!pose.residue_type(seqpos).is_protein()){
				// WARNING why are you storing non-protein in the ProteinSilentReport?
				continue;
			}
			if(!(phi < 0.00001 && psi < 0.00001 && omega < 0.00001) )
			{
				pose.set_phi(seqpos,phi);
				pose.set_psi(seqpos,psi);
				pose.set_omega(seqpos,omega);
			}
			pose.set_secstruct(seqpos,static_cast<char>(secstruct[0]));
			Size nchi(pose.residue_type(seqpos).nchi());
			if(1 <= nchi) pose.set_chi(1, seqpos, chi1);
			if(2 <= nchi) pose.set_chi(2, seqpos, chi2);
			if(3 <= nchi) pose.set_chi(3, seqpos, chi3);
			if(4 <= nchi) pose.set_chi(4, seqpos, chi4);
		}

	}else{
	//	statement stmt = (*db_session) <<
	//		"SELECT\n"
	//		"	seqpos,\n"
	//		"	phi,\n"
	//		"	psi,\n"
	//		"	omega\n"
	//		"FROM\n"
	//		"	protein_residue_conformation\n"
	//		"WHERE\n"
	//		"	protein_residue_conformation.struct_id=?;" << struct_id;
	//
	//	result res(basic::database::safely_read_from_database(stmt));
	//	while(res.next()){
	//		Size seqpos;
	//		Real phi,psi,omega;
	//		res >> seqpos >> phi >> psi >> omega;
	//		if (!pose.residue_type(seqpos).is_protein()){
	//			// WARNING why are you storing non-protein in the ProteinSilentReport?
	//			continue;
	//		}
	//
	//		//pose.set_phi(seqpos,phi);
	//		//pose.set_psi(seqpos,psi);
	//		//pose.set_omega(seqpos,omega);
	//	}
	}
}


//This should be factored out into a non-member function.
void ProteinResidueConformationFeatures::set_coords_for_residues(
		sessionOP db_session,
		boost::uuids::uuid struct_id,
		Pose & pose
){
	// lookup and set all the atoms at once because each query is
	// roughly O(n*log(n) + k) where n is the size of the tables and k
	// is the number of rows returned. Doing it all at once means you
	// only have to pay the n*log(n) cost once.

	std::string statement_string =
		"SELECT\n"
		"	seqpos,\n"
		"	atomno,\n"
		"	x,\n"
		"	y,\n"
		"	z\n"
		"FROM\n"
		"	residue_atom_coords\n"
		"WHERE\n"
		"	residue_atom_coords.struct_id=?;";
	statement stmt(basic::database::safely_prepare_statement(statement_string,db_session));
	//std::string struct_id_string(to_string(struct_id));
	stmt.bind(1,struct_id);

	result res(basic::database::safely_read_from_database(stmt));

	vector1< AtomID > atom_ids;
	vector1< Vector > coords;
	while(res.next())
	{
		Size seqpos, atomno;
		Real x,y,z;
		res >> seqpos >> atomno >> x >> y >> z;

		atom_ids.push_back(AtomID(atomno, seqpos));
		coords.push_back(Vector(x,y,z));
	}
	// use the batch_set_xyz because it doesn't trigger a coordinate
	// update after setting each atom.
	pose.batch_set_xyz(atom_ids,coords);

}

} // namespace
} // namespace
