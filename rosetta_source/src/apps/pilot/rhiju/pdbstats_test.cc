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
#include <core/scoring/constraints/HarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <basic/Tracer.hh>
using basic::T;

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <numeric/xyz.functions.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


using basic::Error;
using basic::Warning;

using namespace core;
//using namespace protocols;
using namespace basic::options::OptionKeys;

using utility::vector1;
using io::pdb::dump_pdb;

OPT_KEY( Boolean, ch_o_bonds )
OPT_KEY( Boolean, fa_cenpack )
OPT_KEY( Boolean, aro_pack )

typedef  numeric::xyzMatrix< Real > Matrix;


///////////////////////////////////////////////////////////////////////////////
// Look more than four atoms away.
bool
path_distance_OK(
			conformation::Residue const & rsd1,
			conformation::Residue const & rsd2,
			Size const ii,
			Size const jj )
{
	Size const & i( rsd1.seqpos() );
	Size const & j( rsd2.seqpos() );
	//	std::cout << i << " " << j << " " << ii << " " << jj << std::endl;
	if ( i == j && rsd1.path_distance(ii,jj) <= 4) return false;
	if ( i == j - 1 ) {
		Size const path_size =
			rsd1.path_distance( ii, rsd1.upper_connect_atom() ) +
			rsd2.path_distance( jj, rsd2.lower_connect_atom() ) + 1;
		if ( path_size <= 4 ) return false;
	}
	if ( i == j + 1 ) {
		Size const path_size =
			rsd1.path_distance( ii, rsd1.lower_connect_atom() ) +
			rsd2.path_distance( jj, rsd2.upper_connect_atom() ) + 1;
		if ( path_size <= 4 ) return false;
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////
void
ch_o_pdbstats_from_pose( utility::io::ozstream & out, pose::Pose & pose, Size const count, char const chain, Size & total_residues )
{

	using numeric::conversions::degrees;

	static Real const DIST_CUTOFF2 ( 4.5 * 4.5 );

	Size const nres = pose.total_residue();

	Size res_count( 0 );

	// Go over acceptors on each residue.
	for (Size i = 1; i <= nres; i++) {

		if ( chain !='?' && pose.pdb_info()->chain(i) != chain ) continue;
		//		std::cout << pose.pdb_resid(i).first << " " << pose.pdb_resid(i).second  << " " << chain <<  std::endl;
		res_count++;

		conformation::Residue const & rsd1( pose.residue( i ) ) ;

		for ( chemical::AtomIndices::const_iterator
						anum  = rsd1.accpt_pos().begin(),
						anume = rsd1.accpt_pos().end(); anum != anume; ++anum ) {

			Size const acc_atm( *anum );
			Size const base_atm ( rsd1.atom_base( acc_atm ) );
			Size const base_atm2( rsd1.abase2   ( acc_atm ) );

			Vector const & acc_atm_xyz( rsd1.atom( acc_atm ).xyz() );
			Vector const & base_atm_xyz( rsd1.atom( base_atm ).xyz() );
			Vector const & base_atm2_xyz( rsd1.atom( base_atm2 ).xyz() );
			//These define a coordinate system.
			Vector const z = ( acc_atm_xyz - base_atm_xyz ).normalized();
			Vector const x_temp = base_atm2_xyz - base_atm_xyz;
			Vector const y = (cross(z,x_temp) ).normalized();
			Vector const x = cross(y, z);
			//			std::cout << "NORMALIZED? " << x.length() << std::endl;

			// Go over hydrogens on all residues.
			for (Size j = 1; j <= nres; j++) {

				if ( chain !='?' && pose.pdb_info()->chain(j) != chain ) continue;

				conformation::Residue const & rsd2( pose.residue( j ) ) ;

				for (Size don_h_atm = rsd2.nheavyatoms() + 1; don_h_atm <= rsd2.natoms(); don_h_atm++ ){

					Vector const & don_h_atm_xyz( rsd2.atom( don_h_atm ).xyz() );
					Size const don_atm ( rsd2.atom_base( don_h_atm ) );
					Vector const & don_atm_xyz( rsd2.atom( don_atm ).xyz() );

					//					if ( i==j && don_atm == acc_atm) continue;
					if ( !path_distance_OK( rsd1, rsd2, acc_atm, don_atm ) ) continue; // Look more than four atoms away.

					if ( ( don_h_atm_xyz - acc_atm_xyz ).length_squared() > DIST_CUTOFF2 ) continue;

					//
					// Save following information --
					// which pdb,
					//   residue number, atom number, residue type, and atom type of acceptor.
					//   residue number, atom number, residue type, and atom type of putative "donor-heavy".
					//   x,y,z of donor-heavy relative to acceptor
					//   x,y,z of donor-h relative to acceptor
					//    r H_A
					//    r D_A
					// angle D--H--A
					// angle H--A--ABASE
					// dihedral D--H--A--ABASE
					//
					// Pretty similar to Kortemme et al.
					//

					//This could be a matrix, I guess.
					Vector const don_h_atm_xyz_rel = don_h_atm_xyz - acc_atm_xyz;
					Vector const don_h_atm_xyz_rel_rotate =
            Vector( dot(x,don_h_atm_xyz_rel), dot(y,don_h_atm_xyz_rel), dot(z,don_h_atm_xyz_rel) );

					Vector const don_atm_xyz_rel = don_atm_xyz - acc_atm_xyz;
					Vector const don_atm_xyz_rel_rotate =
            Vector( dot(x,don_atm_xyz_rel), dot(y,don_atm_xyz_rel), dot(z,don_atm_xyz_rel) );

					out <<
						count << "    " <<
						I(3, rsd1.seqpos()) << " " <<
						I(3, acc_atm) << " " <<
						I(3, Size( rsd1.aa() ) ) << " " <<
						I(3, rsd1.atom_type_index( acc_atm ) ) << "    " <<
						I(3, rsd2.seqpos()) << " " <<
						I(3, don_atm) << " " <<
						I(3, Size( rsd2.aa() ) ) << " " <<
						I(3, rsd2.atom_type_index( don_atm ) ) << "     " <<
						F( 8,3, don_h_atm_xyz_rel_rotate( 1 ) ) << " " <<
						F( 8,3, don_h_atm_xyz_rel_rotate( 2 ) ) << " " <<
						F( 8,3, don_h_atm_xyz_rel_rotate( 3 ) ) << "    " <<
						F( 8,3, don_atm_xyz_rel_rotate( 1 ) ) << " " <<
						F( 8,3, don_atm_xyz_rel_rotate( 2 ) ) << " " <<
						F( 8,3, don_atm_xyz_rel_rotate( 3 ) ) << "    " <<
						F(8,3,don_h_atm_xyz_rel.length()) << " " <<
						F(8,3,don_atm_xyz_rel.length()) << "    " <<
						F(8,3,degrees( angle_radians( don_atm_xyz, don_h_atm_xyz, acc_atm_xyz ) ) )<< " " <<
						F(8,3,degrees( angle_radians( don_h_atm_xyz, acc_atm_xyz, base_atm_xyz ) ) )<< "   " <<
						F(8,3,dihedral_degrees( don_atm_xyz, don_h_atm_xyz, acc_atm_xyz, base_atm_xyz ) )<< std::endl;

				} // don atms
			} // j
		} // acc atms
	} // i

	std::cout << "Processed " << res_count << " residues from chain " << chain << std::endl;

	total_residues += res_count;

}


///////////////////////////////////////////////////////////////////////////////
Vector
get_centroid( conformation::Residue const & rsd )
{

  Vector centroid( 0.0 );
  Size numatoms = 0;
  for ( Size i=rsd.first_sidechain_atom(); i<= rsd.nheavyatoms(); ++i ) {
    centroid += rsd.xyz(i);
    numatoms++;
  }
  if (numatoms > 0 ) {
		centroid /= static_cast< Real >( numatoms );
	} else { //Yo, is this a glycine?
		assert( rsd.aa() == chemical::aa_gly );
		centroid = rsd.xyz( "CA" );
	}

  return centroid;
}


namespace centroid_stats {
	Real const bin_width( 0.25 );
	Size const max_bins( 120 );
	ObjexxFCL::FArray2D< Size > centroid_counts( 20, max_bins );
}


///////////////////////////////////////////////////////////////////////////////
void
fa_cenpack_pdbstats_from_pose( utility::io::ozstream & out, pose::Pose & pose, Size const count, char const chain, Size & total_residues )
{

	static Real const DIST_CUTOFF ( 12.0 );

	Size const nres = pose.total_residue();

	Size res_count( 0 );

	// First, precalculate all centroids.
	utility::vector1 < Vector> centroids;
	for (Size i=1; i <= nres; i++ ) {
		conformation::Residue const & rsd1 = pose.residue(i);
		centroids.push_back( get_centroid( rsd1 ) );
	}
	pose.conformation().update_actcoords(); //Probably not necessary.


	for (Size i = 1; i <= nres; i++) {

		if ( chain !='?' && pose.pdb_info()->chain(i) != chain ) continue;
		res_count++;
		conformation::Residue const & rsd1( pose.residue( i ) ) ;
		Vector const & centroid_i( centroids[i] );
		Vector const & actcoord_i( rsd1.actcoord() );

		for (Size j = i+1; j <= nres; j++) {

			if ( chain !='?' && pose.pdb_info()->chain(j) != chain ) continue;
			conformation::Residue const & rsd2( pose.residue( j ) ) ;
			Vector const & centroid_j( centroids[j] );
			Vector const & actcoord_j( rsd2.actcoord() );

			// Save following information --
			// which pdb,
			//   residue number i
			//   residue type i
			//   residue number j
			//   residue type j
			//   distance[ centroid/centroid ]
			//   distance[ centroid/action-atom ]
			//   distance[ action-atom/centroid ]
			//   distance[ action-atom/action-atom ]
			//  pretty straightforward
			//

			Distance cendist = (centroid_i - centroid_j).length();

			if (cendist < DIST_CUTOFF  && false ) {
				out <<
					count << "    " <<
					I(3, i) << " " <<
					I(3, Size( rsd1.aa() ) ) << " " <<
					I(3, j) << " " <<
					I(3, Size( rsd2.aa() ) ) << "   " <<
					F( 8,3, cendist ) << " " <<
					F( 8,3, (centroid_i - actcoord_j).length() ) << " " <<
					F( 8,3, (actcoord_i - centroid_j).length() ) << " " <<
					F( 8,3, (actcoord_i - actcoord_j).length() ) << std::endl;
			}

			Size cenbin = static_cast<Size> ( cendist / centroid_stats::bin_width) + 1;
			if (cenbin <= centroid_stats::max_bins ) {
				centroid_stats::centroid_counts( Size(rsd1.aa() ), cenbin )++;
				centroid_stats::centroid_counts( Size(rsd2.aa() ), cenbin )++;
			}

		} // j
	} // i

	std::cout << "Processed " << res_count << " residues from chain " << chain << std::endl;

	total_residues += res_count;

}


/////////////////////////////////////////////////////////////////////////////////
kinematics::Stub
get_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid )
{
  using namespace chemical;
  Size res_type = rsd.aa();

  Vector x,y,z;

  // Make an axis pointing from base centroid to Watson-Crick edge.
  std::string WC_atom;
  if ( res_type == aa_phe ) WC_atom = " CZ ";
  if ( res_type == aa_tyr ) WC_atom = " CZ ";
  if ( res_type == aa_trp ) WC_atom = " CZ2";

  Vector const WC_coord (rsd.xyz( WC_atom ) );
  x = WC_coord - centroid;
  x.normalize();

  // Make a perpendicular axis pointing from centroid towards
  // Hoogstein edge (e.g., major groove in a double helix).
  std::string H_atom;
  if ( res_type == aa_phe ) H_atom = " CD1";
  if ( res_type == aa_tyr ) H_atom = " CD1";
  if ( res_type == aa_trp ) H_atom = " CD1";

  Vector const H_coord (rsd.xyz( H_atom ) );
  y = H_coord - centroid; //not orthonormal yet...
  z = cross(x, y);
  z.normalize(); // Should poSize roughly 5' to 3' if in a double helix.

  y = cross(z, x);
  y.normalize(); //not necessary but doesn't hurt.

  //  std::cout << "WC : " << WC_coord << "   H : " << H_coord << "    centroid: " << centroid << std::endl;

  return kinematics::Stub( Matrix::cols( x, y, z ), centroid );
}

////////////////////////////////////////////////////////////////////////////////
void
aro_pack_output(utility::io::ozstream & out, Size const & count,
								Size const & i, Size const & j,
								kinematics::Stub const & stub1,
								conformation::Residue const &  rsd1,
								conformation::Residue const &  rsd2,
								Vector const & centroid2,
								kinematics::Stub const & stub2 )
{
	using namespace chemical;

	out <<
		count << "    " <<
		I(3, i) << " " <<
		I( 3, Size(rsd1.aa()) ) << " " <<
		I(3, j) << " " <<
		I( 3, Size(rsd2.aa()) ) << "    ";


	// output centroid and stub for second base.
	Vector centroid = stub1.global2local( centroid2 );
	out << F( 8,3,centroid(1) ) << " "
			<< F( 8,3,centroid(2) ) << " "
			<< F( 8,3,centroid(3) ) << "    ";

	Vector z_vec2( stub2.M.col_z() );
	Vector centroid_plus_z_vec2 = centroid2 + z_vec2;

	Vector centroid_plus_z_vec = stub1.global2local( centroid_plus_z_vec2 );
	out << F( 8,3,centroid_plus_z_vec(1) ) << " "
			<< F( 8,3,centroid_plus_z_vec(2) ) << " "
			<< F( 8,3,centroid_plus_z_vec(3) ) << "    ";

	// output coordinates for second base.
  utility::vector1< std::string > output_atoms;
	// Doesn't really matter too much...
	if ( rsd2.aa() == aa_phe ) {
		output_atoms = utility::tools::make_vector1< std::string >( " CG ",	" CD1",	" CD2",	" CE1",	" CE2",	" CZ ",	" HD1",	" HE1",	" HZ ",	" HE2",	" HD2");
	} else if ( rsd2.aa() == aa_tyr ) {
		output_atoms = utility::tools::make_vector1< std::string >( " CG ",	" CD1",	" CD2",	" CE1",	" CE2",	" CZ ",	" HD1",	" HE1",	" OH ",	" HE2",	" HD2");
	} else if ( rsd2.aa() == aa_trp ) {
		output_atoms = utility::tools::make_vector1< std::string >( " CG ",	" CD1",	" CD2",	" NE1",	" CE2",	" CE3",	" CZ2",	" CZ3",	" CH2",	" HE1",	" HD1");
	}



	for ( Size n = 1; n <= output_atoms.size(); n++ ) {
		Vector const local_xyz = stub1.global2local( rsd2.xyz( output_atoms[ n ] ) );
		out << F( 8,3,local_xyz(1) ) << " "
				<< F( 8,3,local_xyz(2) ) << " "
				<< F( 8,3,local_xyz(3) ) << "    ";
	}
	out << std::endl;
}


///////////////////////////////////////////////////////////////////////////////
void
aro_pack_pdbstats_from_pose( utility::io::ozstream & out, pose::Pose & pose, Size const count, char const chain, Size & total_residues )
{

	using namespace core::conformation;
	using namespace core::chemical;

	static Real const DIST_CUTOFF ( 8.0 );

	Size const nres = pose.total_residue();

	Size res_count( 0 );

	static bool init( false );

	for (Size i = 1; i <= nres; i++) {

		if ( chain !='?' && pose.pdb_info()->chain(i) != chain ) continue;
		res_count++;
		conformation::Residue const & rsd1( pose.residue( i ) ) ;

		if ( rsd1.aa() != aa_phe && rsd1.aa() != aa_trp && rsd1.aa() != aa_tyr) continue;

		// set up coordinate system for first base.
		Vector centroid1 = get_centroid( rsd1 );
		kinematics::Stub stub1 = get_base_coordinate_system( rsd1, centroid1 );

		if (!init ) {
			aro_pack_output( out, count, i, i, stub1, rsd1, rsd1, centroid1, stub1 );
			init = true;
		}

		for (Size j = 1; j <= nres; j++) {

			if ( i == j ) continue;

			if ( chain !='?' && pose.pdb_info()->chain(j) != chain ) continue;
			conformation::Residue const & rsd2( pose.residue( j ) ) ;

			if ( rsd2.aa() != aa_phe && rsd2.aa() != aa_trp && rsd2.aa() != aa_tyr) continue;

			Vector centroid2 = get_centroid( rsd2 );
			kinematics::Stub stub2 = get_base_coordinate_system( rsd2, centroid2 );

			if ( ( centroid1 - centroid2 ).length() > DIST_CUTOFF ) continue;

			// output coordinates for first base, if this is the first phe encountered.

			aro_pack_output( out, count, i, j, stub1, rsd1, rsd2, centroid2, stub2 );

		} // j
	} // i

	std::cout << "Processed " << res_count << " residues from chain " << chain << std::endl;

	total_residues += res_count;

}

//////////////////////////////////////////
void output_centroid_stats( utility::io::ozstream & out )
{
	using namespace centroid_stats;
	for (Size i = 1; i <= centroid_counts.size1() ; i++ ) {
		for (Size j = 1; j <= centroid_counts.size2() ; j++ ) {
			out << ' ' << I(7,centroid_counts( i,j ));
		}
		out << std::endl;
	}
}

//////////////////////////////////////////
void lowercase(char string[])
{ int i = 0;
	while ( string[i] )
		{
			string[i] = tolower(string[i]);
			i++;
		}
}


///////////////////////////////////////////////////////////////////////////////
void
rhiju_pdbstats()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );
	std::string const file_path( option[ in::path::pdb ]( 1 ) );
	std::string const pdb_list( option[ in::file::l ](1) );

	utility::io::izstream instream( pdb_list );
	if (!instream){
		std::cerr  << "Can't find list file " << pdb_list << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		return;
	}

	std::string outfile  = option[ out::file::o ];
	utility::io::ozstream out( outfile );

	//	for (Size i = 1; i <= pdb_files.size(); i++) {
	//	std::string const pdb_file = pdb_files[i];
	std::string pdb_file;
	char chain(' ');
	Size count( 0 );
	pose::Pose pose;

	Size total_residues( 0 );

	std::string line;
	while ( 	getline( instream, line )  ) {
		std::istringstream line_stream( line );

		line_stream >> pdb_file;

		line_stream >> chain;

		if ( line_stream.fail() ) chain = '?';
		//		chain =  pdb_file.at(4) ;
		//		pdb_file = pdb_file.substr(0,4);
		//		lowercase( pdb_file );

		if (chain == '_' ) chain = ' ';

		core::import_pose::pose_from_pdb( pose, file_path + '/' + pdb_file );

		count++;
		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << " " << chain << std::endl;

		if ( option[ ch_o_bonds] ) {
			ch_o_pdbstats_from_pose( out, pose, count, chain, total_residues );
		}	else if ( option[ fa_cenpack ] ) {
			fa_cenpack_pdbstats_from_pose( out, pose, count, chain, total_residues );
		} else if ( option[ aro_pack ] ) {
			aro_pack_pdbstats_from_pose( out, pose, count, chain, total_residues );
		}
		std::cout << "TOTAL RESIDUES Processed: " << total_residues << std::endl;

	}

	std::cout << "FINISHED ==> TOTAL RESIDUES Processed: " << total_residues << std::endl;

	if (option[fa_cenpack]) {
		utility::io::ozstream out2( "histo_"+outfile );
		output_centroid_stats( out2 );
	}

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	using namespace basic::options;

	NEW_OPT( ch_o_bonds, "Testing existence of CH<-->O bonds", false );
	NEW_OPT( fa_cenpack, "Centroid-centroid pairwise correlations", false );
	NEW_OPT( aro_pack, "Aromatic geometry", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	rhiju_pdbstats();

	exit( 0 );

}
