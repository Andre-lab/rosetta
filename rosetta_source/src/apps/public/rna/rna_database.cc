// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED
#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/pack_rotamers.hh>
// AUTO-REMOVED #include <core/sequence/util.hh>
// AUTO-REMOVED #include <core/io/silent/RNA_SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <basic/basic.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>

#include <numeric/conversions.hh>

#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
// AUTO-REMOVED #include <protocols/rna/RNA_DeNovoProtocol.hh>
// AUTO-REMOVED #include <protocols/rna/RNA_StructureParameters.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end



using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, jump_database )
OPT_KEY( Boolean, vall_torsions )

///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions_test( ){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	utility::vector1 < std::string >  infiles  = option[ in::file::s ]();
	std::string outfile  = option[ out::file::o ];

	utility::io::ozstream torsions_out( outfile );

	for (Size n = 1; n <= infiles.size(); n++ )
  {
	  pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, *rsd_set, infiles[n] );
	  /////////////////////////////////////////
	  protocols::rna::ensure_phosphate_nomenclature_matches_mini( pose );
	  /////////////////////////////////////////

		protocols::rna::create_rna_vall_torsions( pose, torsions_out );

		std::cout << "***********************************************************" << std::endl;
		std::cout << "Put torsions from PDB file " <<  infiles[n] << " into " << outfile << std::endl;
	  std::cout << "***********************************************************" << std::endl;

	}


}

////////////////////////////////////////////////////////////////////////////////////
bool
check_for_contacts( pose::Pose & pose, Size const i,
										Vector const & atom_j, Vector const dir_vector,
										char & edge_i, char & orientation )
{
	static Real const CONTACT_CUTOFF2( 3.0 * 3.0 );

	scoring::rna::RNA_ScoringInfo  & rna_scoring_info( scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	scoring::rna::RNA_CentroidInfo & rna_centroid_info( rna_scoring_info.rna_centroid_info() );
	utility::vector1< Vector > const & base_centroids( rna_centroid_info.base_centroids() );
  utility::vector1< kinematics::Stub > const & base_stubs( rna_centroid_info.base_stubs() );

	//Figure out jump.
	conformation::Residue const & rsd_i( pose.residue( i ) );

	bool found_contact( false );
  for ( Size ii=rsd_i.first_sidechain_atom()+1 ; ii<= rsd_i.nheavyatoms(); ++ii ) {
		//		if ( rsd_i.atom_name( ii ) == " O2*" ) continue;
		Real const dist2( (rsd_i.xyz( ii ) - atom_j ).length_squared() ) ;
    if ( dist2 < CONTACT_CUTOFF2 ) {
			//			std::cout << dist2 << " " <<  i << " " << rsd_i.atom_name( ii ) << std::endl;
			found_contact = true;
			break;
		}
  }

	if (!found_contact) return false;

	Vector const & centroid_i( base_centroids[i] );
	kinematics::Stub const & stub_i( base_stubs[i] );
	Matrix const & M_i( stub_i.M );
	Vector const & x_i = M_i.col_x();
	Vector const & y_i = M_i.col_y();
	Vector const & z_i = M_i.col_z();

	Vector d_ij = atom_j - centroid_i;
	Real const dist_x = dot_product( d_ij, x_i );
	Real const dist_y = dot_product( d_ij, y_i );
	//	Real const dist_z = dot_product( d_ij, z_i );
	//	Real const rho2 = dist_x*dist_x + dist_y*dist_y;

	Real const zeta = numeric::conversions::degrees( std::atan2( dist_y, dist_x) );
	if (std::abs(zeta) < 60.0) edge_i = 'W';  //Watson-Crick edge
	else if ( zeta > +60.0 )   edge_i = 'H'; // Hoogsteen edge
	else                       edge_i = 'S'; // Sugar edge

	orientation = dot( dir_vector, z_i ) > 0 ?  'A' : 'P';

	return true;
}

/////////////////////////////////////////////////////////////////////////////////
void
check_for_contacts_and_output_jump_o2star( pose::Pose & pose, Size const i, Size const j,
																					 utility::io::ozstream & dataout ){

	char edge_i, orientation;

	conformation::Residue const & rsd_i( pose.residue( i ) );
	conformation::Residue const & rsd_j( pose.residue( j ) );

	std::string const atom_name = " O2*";
	Vector const atom_vector = rsd_j.xyz( atom_name );
	Vector const dir_vector  = rsd_j.xyz( " O2*" ) - rsd_j.xyz( " C2*" );
	if ( !check_for_contacts( pose, i, atom_vector, dir_vector, edge_i, orientation) ) return;

	char const edge_j = '2';

	kinematics::Stub const stub1( rsd_i.xyz( rsd_i.chi_atoms(1)[4] ),
																rsd_i.xyz( rsd_i.chi_atoms(1)[3] ),
																rsd_i.xyz( rsd_i.chi_atoms(1)[2] ) );

	kinematics::Stub const stub2( rsd_j.xyz( atom_name ),
																rsd_j.xyz( " C2*" ),
																rsd_j.xyz( " C3*" ) );

	dataout << "PAIR " <<
		I(5, i) << ' ' << edge_i << ' ' <<
		I(5, j) << ' ' << edge_j << "   " <<
		orientation << "   " <<
		pose.residue(i).name1() << ' ' << pose.residue(j).name1() << " " <<
		rsd_i.atom_name( rsd_i.chi_atoms(1)[4] ) <<  " " <<
		atom_name <<  " " <<
		kinematics::Jump( stub1, stub2) <<
		std::endl;

}


/////////////////////////////////////////////////////////////////////////////////
void
check_for_contacts_and_output_jump_phosphate( pose::Pose & pose, Size const i, Size const j,
																							utility::io::ozstream & dataout ){

	char edge_i, orientation;

	conformation::Residue const & rsd_i( pose.residue( i ) );
	conformation::Residue const & rsd_j( pose.residue( j ) );
	conformation::Residue const & prev_rsd( pose.residue( j-1 ) );


	Vector dir_vector  = cross( rsd_j.xyz( " O2P" ) - rsd_j.xyz( " P  " ),
															rsd_j.xyz( " O1P" ) - rsd_j.xyz( " P  " ) );

	std::string atom_name = " O1P";
	Vector atom_vector = rsd_j.xyz( atom_name );
	if ( !check_for_contacts( pose, i, atom_vector, dir_vector, edge_i, orientation) ) {
		atom_name = " O2P";
		atom_vector = rsd_j.xyz( atom_name );
		if (!check_for_contacts( pose, i, atom_vector, dir_vector, edge_i, orientation) ) return;
	}


	char const edge_j = 'P';

	kinematics::Stub const stub1( rsd_i.xyz( rsd_i.chi_atoms(1)[4] ),
																rsd_i.xyz( rsd_i.chi_atoms(1)[3] ),
																rsd_i.xyz( rsd_i.chi_atoms(1)[2] ) );

	kinematics::Stub const stub2_fwd( rsd_j.xyz( atom_name ),
																		rsd_j.xyz( " P  " ),
																		prev_rsd.xyz( " O3*" ) );

	kinematics::Stub const stub2_back( rsd_j.xyz( atom_name ),
																		 rsd_j.xyz( " P  " ),
																		 rsd_j.xyz( " O5*" ) );

	dataout << "PAIR " <<
		I(5, i) << ' ' << edge_i << ' ' <<
		I(5, j) << ' ' << edge_j << "   " <<
		orientation << "   " <<
		pose.residue(i).name1() << ' ' << pose.residue(j).name1() << " " <<
		rsd_i.atom_name( rsd_i.chi_atoms(1)[4] ) <<  " " <<
	  atom_name <<  " " <<
		kinematics::Jump( stub1, stub2_fwd ) <<
		kinematics::Jump( stub1, stub2_back ) <<
		std::endl;

}

//////////////////////////////////////////////////////////////////////////////////////
// JUMP extractor.
void
create_bp_jump_database_test( ){

	using namespace chemical;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	std::string infile  = option[ in::file::s ][1];
	std::string outfile  = option[ out::file::o ];

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, *rsd_set, infile );
	protocols::rna::ensure_phosphate_nomenclature_matches_mini( pose );

	//	figure_out_reasonable_rna_fold_tree( pose );

	// Fill base pairing information... these are
	// all functions used in scoring... see RNA_BaseBaseEnergy.cc
	ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	(*lores_scorefxn)( pose );
	lores_scorefxn->show( std::cout, pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	Energy_base_pair_list scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	utility::io::ozstream dataout( outfile );

	for ( Energy_base_pair_list::const_iterator it = scored_base_pair_list.begin();
				it != scored_base_pair_list.end(); ++it ){

		Base_pair const base_pair = it->second;

		int const i = base_pair.res1;
		int const k = base_pair.edge1;

		int const j = base_pair.res2;
		int const m = base_pair.edge2;

		char const orientation = ( base_pair.orientation == 1) ? 'A' : 'P';

		char const edge_i = get_edge_from_num( k );
		char const edge_j = get_edge_from_num( m );

		//Figure out jump.
		conformation::Residue const & rsd_i( pose.residue( i ) );
		conformation::Residue const & rsd_j( pose.residue( j ) );
		kinematics::Stub const stub_i( rsd_i.xyz( rsd_i.chi_atoms(1)[4] ),
																	 rsd_i.xyz( rsd_i.chi_atoms(1)[3] ),
																	 rsd_i.xyz( rsd_i.chi_atoms(1)[2] ) );
		kinematics::Stub const stub_j( rsd_j.xyz( rsd_j.chi_atoms(1)[4] ),
																	 rsd_j.xyz( rsd_j.chi_atoms(1)[3] ),
																	 rsd_j.xyz( rsd_j.chi_atoms(1)[2] ) );

		dataout << "PAIR " <<
			I(5, i) << ' ' << edge_i << ' ' <<
			I(5, j) << ' ' << edge_j << "   " <<
			orientation << "   " <<
			pose.residue(i).name1() << ' ' << pose.residue(j).name1() << " " <<
			rsd_i.atom_name( rsd_i.chi_atoms(1)[4] ) <<  " " <<
			rsd_j.atom_name( rsd_j.chi_atoms(1)[4] ) <<  " " <<
		  kinematics::Jump( stub_i, stub_j) <<
			std::endl;

	}

	//How about 2' and Phosphate jumps?
	// Look at each base.
	core::scoring::EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for (Size i = 1; i <= pose.total_residue(); i++ ){

		// Neighboring residues making base-phosphate or base-2'OH contacts?
		for ( graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node(i)->const_edge_list_begin(),
						irue = energy_graph.get_node(i)->const_edge_list_end();
					iru != irue; ++iru ) {
			EnergyEdge const * edge( static_cast< EnergyEdge const *> ( *iru ) );
			Size const j( edge->get_other_ind(i) );
			//			EnergyGraph const & energy_graph( pose.energies().energy_graph() );

			check_for_contacts_and_output_jump_o2star( pose, i, j, dataout );

			if ( j > 1  )	check_for_contacts_and_output_jump_phosphate( pose, i, j, dataout );

		}
	}


	dataout.close();

	std::cout << "***********************************************************" << std::endl;
	std::cout << "Put jumps from PDB file " <<  infile << " into " << outfile << std::endl;
	std::cout << "***********************************************************" << std::endl;

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options;

	if ( option[ jump_database ] ) {
		create_bp_jump_database_test();
	} else if ( option[ vall_torsions ] ) {
		create_rna_vall_torsions_test();
	} else {
		std::cout << std::endl;
		std::cout << "Please specify: " << std::endl;
		std::cout << "  -vall_torsions   for generating torsions library " << std::endl;
		std::cout << "  -jump_database   for generating database of rigid-body orientations " << std::endl;
		std::cout << std::endl;
	}

	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace basic::options;

	//Uh, options? MOVE THESE TO OPTIONS NAMESPACE INSIDE CORE/OPTIONS.
	NEW_OPT( vall_torsions, "Generate a torsions file from a big RNA file", false );
	NEW_OPT( jump_database, "Generate a database of jumps extracted from base pairings from a big RNA file", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

}
