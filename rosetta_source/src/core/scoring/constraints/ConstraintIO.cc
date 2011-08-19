// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IO-functionality for Constraints
/// @brief
/// @author Oliver Lange olange@u.washington.edu

// Unit headers
#include <core/scoring/constraints/ConstraintIO.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
//#include <core/scoring/constraints/ConstraintForest.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
//#include <core/scoring/constraints/BindingSiteConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
// AUTO-REMOVED #include <core/scoring/constraints/MixtureFunc.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ScalarWeightedFunc.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <utility/excn/Exceptions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

static basic::Tracer tr("core.io.constraints");

namespace core {
namespace scoring {
namespace constraints {

ConstraintIO*
ConstraintIO::get_instance() {
	if ( instance_ == 0 ) {
		instance_ = new ConstraintIO();
	}
	return instance_;
}

FuncFactory & ConstraintIO::get_func_factory(void) {
	return func_factory_;
}

ConstraintFactory & ConstraintIO::get_cst_factory(void) {
	return * ConstraintFactory::get_instance();
}

ConstraintIO* ConstraintIO::instance_ = 0;
FuncFactory ConstraintIO::func_factory_;
//ConstraintFactory ConstraintIO::cst_factory_;


void
ConstraintIO::read_cst_atom_pairs(
	std::istream & data,
	std::string & next_section_name,
	ConstraintSet & cst_set,
	pose::Pose const & pose
) {
	tr.Debug << "ConstraintIO::read_cst_atom_pairs" << std::endl;
	std::string line;
	while( getline( data, line ) ) {
		Size res1, res2;
		std::string tempres1, tempres2;
		std::string name1, name2;
		std::string func_type;
		std::string type;
		std::istringstream line_stream( line );
		line_stream
			>> name1 >> tempres1
			>> name2 >> tempres2
			>> func_type;

		parse_residue( pose, tempres1, res1 );
		parse_residue( pose, tempres2, res2 );

		if ( name1.find("[" )!=std::string::npos ) { //end of this section
			tr.Debug << "section end detected in line " << line << std::endl;
			next_section_name = line;
			return;
		}
		tr.Debug 	<< "read: " << name1 << " " << name2 << " "
							<< res1 << " " << res2 << " func: " << func_type
							<< std::endl;

		if ( res1>pose.total_residue() || res2> pose.total_residue() ||
		     ( !pose.residue_type( res1 ).has_atom_name( name1 ) ) ||
		     ( !pose.residue_type( res2 ).has_atom_name( name2 ) )
		){

			tr.Error << "error in constraint (no such atom in pose!)"
									<< name1 << " " << name2 << " " << res1 << " " << res2 << " func: " << func_type << std::endl;
			utility_exit_with_message ("Constraint data referred to atom which is not present in pose");
		}


		id::AtomID atom1( pose.residue_type( res1 ).atom_index( name1 ), res1 );
		id::AtomID atom2( pose.residue_type( res2 ).atom_index( name2 ), res2 );

		FuncOP aFunc = func_factory_.func_types_[ func_type ]->clone();
		aFunc->read_data( line_stream );

		if ( tr.Debug.visible() ) {
			aFunc->show_definition( tr.Debug ); tr.Debug << std::endl;
		}

		cst_set.add_constraint( new AtomPairConstraint( atom1, atom2, aFunc ) );
	} // while getline
	tr.Debug << "end of file reached" << std::endl;
	next_section_name = "";
}

ConstraintOP ConstraintIO::parse_atom_pair_constraint(
	std::istream & data,
	core::pose::Pose pose
) {
	Size res1, res2;
	std::string tempres1, tempres2;
	std::string name1, name2;
	std::string func_type;
	std::string type;

	data
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> func_type;

	parse_residue( pose, tempres1, res1 );
	parse_residue( pose, tempres2, res2 );

	tr.Info 	<< "read: " << name1 << " " << name2 << " "
						<< res1 << " " << res2 << " func: " << func_type
						<< std::endl;
	if ( res1>pose.total_residue() || res2> pose.total_residue() ) {
		tr.Warning 	<< "ignored constraint (no such atom in pose!)"
								<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		return NULL;
	}

	id::AtomID atom1( pose.residue_type( res1 ).atom_index( name1 ), res1 );
	id::AtomID atom2( pose.residue_type( res2 ).atom_index( name2 ), res2 );

	FuncOP aFunc = func_factory_.func_types_[ func_type ]->clone();
	aFunc->read_data( data );

	if ( tr.Debug.visible() ) {
		aFunc->show_definition( tr.Debug ); tr.Debug<<std::endl;
	}

	ConstraintOP cst_op( new AtomPairConstraint( atom1, atom2, aFunc ) );
	return cst_op;
} // parse_atom_pair_constraint

ConstraintOP ConstraintIO::parse_coordinate_constraint(
	std::istream & data,
	core::pose::Pose pose
) {
	Real x, y, z;
	Size fixed_res, other_res;
	std::string tempfixed_res, tempother_res;
	std::string fixed_res_name, other_res_name;
	std::string func_type;
	std::string type;

	data
		>> fixed_res_name >> tempfixed_res
		>> other_res_name >> tempother_res
		>> x >> y >> z
		>> func_type;

	parse_residue( pose, tempfixed_res, fixed_res );
	parse_residue( pose, tempother_res, other_res );

	tr.Debug 	<< "read: " << fixed_res_name << " " << other_res_name << " "
						<< fixed_res << " " << other_res << " func: " << func_type << std::endl;
	if ( fixed_res > pose.total_residue() || other_res > pose.total_residue() ) {
		tr.Warning 	<< "ignored constraint (no such atom in pose!)"
								<< fixed_res_name << " " << other_res_name << " "
								<< fixed_res << " " << other_res << std::endl;
		return NULL;
	}

	id::AtomID atom1( pose.residue_type( fixed_res ).atom_index( fixed_res_name ), fixed_res );
	id::AtomID atom2( pose.residue_type( other_res ).atom_index( other_res_name ), other_res );

	FuncOP aFunc = func_factory_.func_types_[ func_type ]->clone();
	aFunc->read_data( data );

	if ( tr.Debug.visible() ) {
		aFunc->show_definition( tr.Debug ); tr.Debug << std::endl;
	}

	Vector transform( x, y, z );
	ConstraintOP cst_op( new CoordinateConstraint( atom1, atom2, transform, aFunc ) );
	return cst_op;
} // parse_coordinate_constraint


void
ConstraintIO::read_cst_angles(
	std::istream & data,
	std::string & next_section_name,
	ConstraintSet & cst_set,
	pose::Pose const & pose
) {
	tr.Debug << "ConstraintIO::read_cst_angles" << std::endl;
	std::string line;
	while( getline( data, line ) ) {
		Size res1, res2, res3;
		std::string tempres1, tempres2, tempres3;
		std::string name1, name2, name3;
		std::string func_type;
		std::string type;
		std::istringstream line_stream( line );
		line_stream
			>> name1 >> tempres1
			>> name2 >> tempres2
			>> name3 >> tempres3
			>> func_type;

		parse_residue( pose, tempres1, res1 );
		parse_residue( pose, tempres2, res2 );
		parse_residue( pose, tempres3, res3 );

		if ( name1.find("[" )!=std::string::npos ) { //end of this section
			tr.Debug << "section end detected in line " << line << std::endl;
			next_section_name = line;
			return;
		}
		id::AtomID atom1( pose.residue_type( res1 ).atom_index( name1 ), res1 );
		id::AtomID atom2( pose.residue_type( res2 ).atom_index( name2 ), res2 );
		id::AtomID atom3( pose.residue_type( res3 ).atom_index( name3 ), res3 );

		tr.Debug	<< "read: " << name1 << " " << res1 << " "  <<  name2 << " "
							<< res2 << " " << name3 << " " << res3 << std::endl;
		FuncOP aFunc = func_factory_.func_types_[ func_type ]->clone();
		aFunc->read_data( line_stream );

		if ( tr.Debug.visible() ) {
			aFunc->show_definition( tr.Debug ); tr.Debug<<std::endl;
		}

		if ( res1 > pose.total_residue() || res2 > pose.total_residue() || res3 > pose.total_residue() ) {
			tr.Warning << "ignored constraint" << name1 << " " << name2 << " "
									<< res1 << " " << res2 << std::endl;
			continue;
		}

		cst_set.add_constraint( new AngleConstraint( atom1, atom2, atom3, aFunc ) );
	}// while getline
	tr.Debug << "end of file reached" << std::endl;
	next_section_name = "";
}

/*
void
ConstraintIO::read_cst_bindingsites(
	std::istream & data,
	std::string & next_section_name,
	ConstraintSet & cst_set,
	pose::Pose const & pose
) {
	tr.Debug << "ConstraintIO::read_cst_angles" << std::endl;
	std::string line;
	while( getline( data, line ) ) {
		Size res;
		std::string tempres;
		std::string name;
		std::istringstream line_stream( line );

		if ( line.find("[" )!=std::string::npos ) { //end of this section
			tr.Debug << "section end detected in line " << line << std::endl;
			next_section_name = line;
			return;
		}

		// to do? atm name -> centroid conversion
		utility::vector1< id::AtomID > atms;
		tr.Debug << "read: ";
		while ( line_stream >> name >> tempres ) {

			parse_residue( pose, tempres, res );

			tr.Debug << "   " << name << " " << res ;

			atms.push_back( id::AtomID( pose.residue_type( res ).atom_index( name ), res ) );
			if ( res > pose.total_residue() ) {
				tr.Debug << "** ignored **";
				continue;
			}
		}
		tr.Debug << std::endl;

		cst_set.add_constraint( new BindingSiteConstraint( atms , pose ) );
	}// while getline
	tr.Debug << "end of file reached" << std::endl;
	next_section_name = "";
}
*/

std::string
get_section_name ( std::string line ) {
	if ( line.size() == 0 ) return line;
	std::istringstream line_stream( line );
	std::string tok;
	line_stream >> tok;
	int start = 0;
	if ( tok == "[" ) { //
		line_stream >> tok;
		start = 0;
	} else {
		std::string::size_type start = tok.find("[");
		if ( start != 0 ) return "NO_SECTION";
		//utility_exit_with_message (" reading constraints file, expected: [ section ] ");
		start = 1;
	}

	std::string::size_type loc = tok.find("]");
	if ( loc != std::string::npos ) {
		return tok;
	} else {
		return tok.substr(start,loc);
	}
}

void
ConstraintIO::read_cst_bindingsites(
   std::istream & data,
	 std::string & next_section_name,
	 ConstraintSet & cst_set,
	 pose::Pose const & pose
) {
	tr.Debug << "ConstraintIO::read_cst_angles" << std::endl;
	std::string line;
	while( getline( data, line ) ) {
		//mjo commenting out 'res' because it is unused and causes a warning
		//Size res;
		std::string name;
		std::istringstream line_stream( line );
		if ( line.find("[" )!=std::string::npos ) { //end of this section
			tr.Debug << "section end detected in line " << line << std::endl;
			next_section_name = line;
			return;
		}

		// to do? atm name -> centroid conversion
		/*utility::vector1< id::AtomID > atms;
			tr.Debug << "read: ";
			while ( line_stream >> name >> res ) {
			tr.Debug << "   " << name << " " << res ;
			atms.push_back( id::AtomID( pose.residue_type( res ).atom_index( name ), res ) );
			if ( res > pose.total_residue() ) {
			tr.Debug << "** ignored **";
			continue;
			}
			}
			tr.Debug << std::endl;
			cst_set.add_constraint( new BindingSiteConstraint( atms , pose ) );*/

		ConstraintOP bsc = ConstraintFactory::get_instance()->newConstraint( "BindingSite" );
		bsc->read_def( line_stream, pose, get_func_factory() );
		cst_set.add_constraint( bsc );
	}// while getline
	tr.Debug << "end of file reached" << std::endl;
	next_section_name = "";
}


ConstraintSetOP
ConstraintIO::read_constraints(
	std::string const & fname,
	ConstraintSetOP cset,
	pose::Pose const & pose
) {
	utility::io::izstream data( fname.c_str() );
	tr.Info << "read constraints from " << fname << std::endl;
	if ( !data ) {
		utility_exit_with_message( "[ERROR] Unable to open constraints file: "+ fname );
	}

	std::string line;
	getline(data,line); // header line
	std::string section = get_section_name ( line );
	std::string pre_read;
	while ( section.size() ) {
		tr.Info << "read constraints section --" << section << "---" << std::endl;
		if ( section ==  "atompairs" ) {
			read_cst_atom_pairs( data, pre_read, *cset, pose );
		} else if ( section == "angles" ) {
			read_cst_angles( data, pre_read, *cset, pose );
		} else if ( section == "bindingsites" ) {
			read_cst_bindingsites( data, pre_read, *cset, pose );
		} else if ( section == "NO_SECTION" ) {
			tr.Info << " no section header [ xxx ] found, try reading line-based format... DON'T MIX"
							<< std::endl;
			return read_constraints_new( fname, cset, pose );
		} else { //section header, but unknown name
			utility_exit_with_message(
				"constraint-file: " + fname + " section " + section + " not recognized!"
			);
		}
		tr.Trace << "pre_read: " << pre_read << std::endl;
		section = get_section_name ( pre_read );
	}
//	pose.constraint_set( cset );
	return cset;
} // read_constraints

/////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details read in constraints with new format from a file
ConstraintSetOP
ConstraintIO::read_constraints_new(
	std::string const & fname,
	ConstraintSetOP cset,
	pose::Pose const & pose
) {
	utility::io::izstream data( fname.c_str() );
	tr.Info << "read constraints from " << fname << std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open constraints file: " + fname );
	}

	while( data.good() ) { // check if we reach the end of file or not
		// read in each constraint and add it constraint_set
		ConstraintOP cst_op;
		cst_op = read_individual_constraint_new( data, pose, get_func_factory() );
		if ( cst_op ) {
			cset->add_constraint( cst_op );
		} else if ( ! data.eof() ) { // not end of line
			tr.Error << "ERROR: reading constraints from file" << fname << std::endl;
			break;
		}
	} // while
	return cset;
} // read_constraints_new

ConstraintOP
ConstraintIO::read_individual_constraint_new(
	std::istream & data,
	pose::Pose const& pose,
	FuncFactory const & func_factory,
	std::string tag
)
{
	ConstraintOP cst_op;
	cst_op = get_cst_factory().newConstraint( tag );

	std::string error_msg("");
	//	bool error_seen( false );

	if ( cst_op ) {
		//		try {
		cst_op->read_def( data, pose, func_factory );
		// } catch ( utility::excn::EXCN_Exception &excn  ) {
// 			tr.Error << "ERROR: reading of " + tag + " failed.\n" << excn << std::endl;
// 			cst_op = NULL;
			//}
		if ( !data.good() && !data.eof()) {
			error_msg += "ERROR: reading of " + tag + " failed.\n";
			tr.Error << error_msg << std::endl;
			cst_op = NULL;
		}
	} else {
		error_msg += "ERROR: constraint type " + tag + " not known.\n";
		tr.Error << error_msg << std::endl;
		cst_op = NULL;
	}
	if ( !cst_op ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if ( option[ OptionKeys::constraints::exit_on_bad_read ]() ) { //OL Aug 12 2010,
			//changed option default from 'true' to 'false'. This was way confusing!
			//hard exits can be rather annoying... after waiting for 24h in the queue of a cluster.
			utility_exit_with_message( error_msg );
		}
	}
	return cst_op;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details read in each individual constraint. The constraints could be single-line ones such as AtomPair,
/// Angle and Diehedral, or multi-line ones such as Multiconstraint and Ambiguous. Return owing pointer of the
/// constraint if read in successfully, otherwise return NULL. Identify the type of constraint and call
/// each constraint's read_def function to finish reading. Skip lines beginning with '#' or '\n'
/// May be called recursively in the case of Multiconstraint and AmbiguousConstraint.
/// @note the istream data should point to the beginning of a line when this function is called.
ConstraintOP
ConstraintIO::read_individual_constraint_new(
	std::istream & data,
	pose::Pose const& pose,
	FuncFactory const & func_factory
)
{
 	std::string tag, dummy_line;
	// get the ConstraintType tag
	while( true ) {
		char c = data.peek(); // get first char of the line but not moving istream pointer
		if ( c == '#' || c == '\n' ) { //ignore # comment line and empty line
			while( data.good() && (data.get() != '\n') ) {}
			continue;
		}
		if ( data.eof() ) return NULL;
		data >> tag;
		if ( data.fail() ) {
			tr.Error << "can't read constraint type" << std::endl;
			return NULL;
		}
		break;
	}
	if (( tag.substr(0,3) == "END" )||( tag.substr(0,3) == "End" )) return NULL; // stopper for MultiConstraint or AmbiguousConstraint
	return read_individual_constraint_new( data, pose, func_factory, tag );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void
ConstraintIO::write_constraints( std::ostream& out, ConstraintSet const& cst_set, pose::Pose const& pose ) {
	cst_set.show_definition( out, pose );
}

void
ConstraintIO::write_constraints( std::string const& filename, ConstraintSet const& cst_set, pose::Pose const& pose ) {
	utility::io::ozstream dump_cst( filename );
	write_constraints( dump_cst, cst_set, pose );
}


//ConstraintForestOP
//ConstraintIO::read_constraint_forest( std::string const & fname, pose::Pose & pose ) {
//	utility::io::izstream data( fname.c_str() );
//	tr.Info << "read ConstraintForest from " << fname << std::endl;
//	if ( !data ) {
//		utility_exit_with_message( "[ERROR] Unable to open ConstraintForest file: "+fname );
//	}
//
//	ConstraintForestOP cf = new ConstraintForest();
//	cf->read_forest(data, pose);
//	return cf;
//}


void
ConstraintIO::parse_residue( pose::Pose const& pose, std::string const residue_string, Size & residue_num )
{
	std::stringstream data;
	char chain;
	Size resnum;

	data.str( residue_string );

	data >> resnum;

	if ( !(data >> chain).fail() ) {
		residue_num = pose.pdb_info()->pdb2pose( chain, resnum );
	} else residue_num = resnum;
}

} //constraints
} //scoring
} //core

