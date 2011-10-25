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
/// @author Phil Bradley

// Unit header
#include <core/chemical/residue_io.hh>


// Rosetta headers
//#include <core/chemical/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueSupport.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/orbitals/OrbitalType.hh>

// Commented by inclean daemon #include <core/chemical/Adduct.hh>

#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.hh>
#include <core/kinematics/Stub.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
//#include <core/pack/dunbrack/SingleLigandRotamerLibrary.hh>

// Commented by inclean daemon #include <basic/basic.hh>
#include <basic/Tracer.hh>

// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/ObjexxFCL.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility headers
//#include <utility/io/izstream.hh>
// Commented by inclean daemon #include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

// Numeric headers
#include <numeric/conversions.hh>

//Auto Headers
#include <core/id/DOF_ID.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/io/izstream.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


//option
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end



// C++ headers
// Commented by inclean daemon #include <iostream>
// Commented by inclean daemon #include <sstream>
// Commented by inclean daemon #include <fstream>

namespace core {
namespace chemical {

static basic::Tracer tr("core.chemical");


// /// must be a better place for this, probably already exists!
// inline
// std::string
// strip_whitespace( std::string const & name )
// {
// 	std::string trimmed_name( name );
// 	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to dothis?
// 	return trimmed_name;
// }

///////////////////////////////////////////////////////////////////////////////
/// helper fxn

id::AtomID
atom_id_from_icoor_line(
	std::string const name,
	ResidueType const & rsd
)
{
	using id::AtomID;
	ICoorAtomID id( name, rsd );

	switch ( id.type() ) {
	case ICoorAtomID::INTERNAL:
		return AtomID( id.atomno(), 1 );
	case ICoorAtomID::CONNECT:
		return AtomID( id.atomno(), 2 );
	case ICoorAtomID::POLYMER_LOWER:
		return AtomID( 1, 3 );
	case ICoorAtomID::POLYMER_UPPER:
		return AtomID( 2, 3 );
	default:
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
	return id::BOGUS_ATOM_ID;
}

///////////////////////////////////////////////////////////////////////////////
/// @details	 construct a ResidueType from a file. Example files are currently in
///	 minirosetta_database/chemical/residue_type_sets/fa_standard/residue_types/???.params
///  These files contain information
///	 about each basic ResidueType which can be patched to created various
///	 variant types.

ResidueTypeOP
read_topology_file(
	std::string const & filename,
	chemical::AtomTypeSetCAP atom_types,
	chemical::ElementSetCAP elements,
	chemical::MMAtomTypeSetCAP mm_atom_types,
	chemical::orbitals::OrbitalTypeSetCAP orbital_atom_types,
//	chemical::CSDAtomTypeSetCAP csd_atom_types kwk commenting out until they have been fully implemented
	chemical::ResidueTypeSetCAP rsd_type_set
)
{
	using id::AtomID;
	using id::DOF_ID;
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	using namespace basic;

	// read the file
	utility::vector1< std::string > lines;
	{
		std::string line;
		if( ! utility::file::file_exists( filename ) ) {
			utility_exit_with_message("Cannot find file '"+filename+"'");
		}
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message("Cannot open file '"+filename+"'");
		}
// 		utility::io::izstream data( filename );
		while ( getline( data, line ) ) {
			std::istringstream l( line );
			//if ( line.size() < 1 || line[0] == '#' ) continue;
			if ( line.size() < 1 ) continue;
			std::string::size_type pound= line.find('#', 0);
			if( pound == std::string::npos ) lines.push_back( line );
			else{
				std::string no_comment_line= line.substr(0, pound);
				lines.push_back(no_comment_line);
			}
		}
		tr.Debug << "Read " << lines.size() << " lines from file: " <<
			filename << std::endl;
		data.close();
	}


	// decide what type of Residue to instantiate
	// would scan through for the TYPE line, to see if polymer or ligand...
	//
	// Residue needs a pointer to the AtomTypeSet object for setting up
	// atom-type dependent data
	//
	// You may be asking yourself, at least I was, why the hell do we scan this file more than once?
	// The problem is that while reading the params (topology file), we need to assign private member variable
	// data in ResidueType. We want to make sure that we get the correct number of atoms assigned
	// before we start adding bonds, icoor, etc. This allows to provide checks to make sure certain
	// things are being assigned correctly, ie adding bonds correctly, setting icoor values with correct placement
	// of stub atoms, etc etc.
	ResidueTypeOP rsd( new ResidueType( atom_types, elements, mm_atom_types, orbital_atom_types ) ); //kwk commenting out until atom types are fully implemented , csd_atom_types ) );
	rsd->residue_type_set( rsd_type_set ); // give this rsd_type a backpointer to its set

	// add the atoms
	Size const nlines( lines.size() );
	Size natoms(0), norbitals(0);
	for (Size i=1; i<= nlines; ++i ) {
		std::string line( lines[i] );
		if (line.size() > 0) {
			while (line.substr(line.size()-1) == " ") {
				line = line.substr(0,line.size()-1);
				if (line.size() == 0) break;
			}
		}

		std::istringstream l( line );
		std::string tag;

		l >> tag;
		//		if ( line.size() < 5 || line.substr(0,5) != "ATOM " ) continue;
		if ( tag != "ATOM" ) continue;

		// the atom name for this atom
		std::string const atom_name( line.substr(5,4) );
		l >> tag; // std::string const atom_name( tag );

		// the atom type name -- must match one of the atom types
		// for which force-field parameters exist
		// std::string const atom_type_name( line.substr(10,4) );
		l >> tag; std::string const atom_type_name( tag );

		// read in the Molecular mechanics atom type
		// std::string const mm_atom_type_name( line.substr(15,4) );
		l >> tag; std::string const mm_atom_type_name( tag );

		// the atomic charge
		float charge;
		// std::istringstream l( line.substr(20) );
		l >> charge;
		float parse_charge(charge);
		if (!l.eof()) {
			l >> parse_charge;
		}
		
    if ( ! basic::options::option[ basic::options::OptionKeys::corrections::chemical::parse_charge ]() ) {
			rsd->add_atom( atom_name, atom_type_name, mm_atom_type_name, charge );
		}
		else {
			rsd->add_atom( atom_name, atom_type_name, mm_atom_type_name, parse_charge );
		}

		++natoms;
	}

	// No ATOM lines probably means an invalid file. Perhaps someone made a mistake with an -extra_res_fa flag.
	// Fail gracefully now, versus a segfault later.
	if(natoms == 0) {
		utility_exit_with_message("Residue topology file '" + filename + "' does not contain valid ATOM records.");
	}

	// add the bonds, parse rest of file
	bool found_AA_record = false;
	bool found_PDB_ROTAMERS_record = false;
	std::string pdb_rotamers_filename = "";
	for (Size i=1; i<= nlines; ++i ) {
		std::string const & line( lines[i] );
		std::istringstream l( line );
		std::string tag,atom1,atom2,atom3,atom4, rotate, connect_type, orbitals_tag, orbital;
		l >> tag;
		if ( l.fail() ) continue;
		if ( tag == "CONNECT" ) {
			l >> atom1;
			l >> rotate; // not used here
			if ( l >> connect_type) rsd->add_residue_connection( atom1,  connect_type);
			else rsd->add_residue_connection( atom1);
			//std::cout << "CONNECT record depricated " << std::endl;
		} else if ( tag == "TYPE" ) {
			// will probably handle this differently later on
			l >> tag;
			if ( tag == "POLYMER" ) {
				rsd->add_property( tag );
			} else if ( tag == "LIGAND" ) {
				rsd->add_property( tag );
			}
			// want to know if we're a polymer for other setup decisions
//		} else if ( tag == "CSD_TYPE" ) { //kwk commenting out until CSD types have been fully implemented
//			l >> atom1 >> atom2;
//			rsd->set_csd_atom_type( atom1, atom2 );
//
		} else if ( tag == "BOND" ) {
			l >> atom1 >> atom2;
			rsd->add_bond( atom1, atom2 );

		} else if ( tag == "CUT_BOND" ) {
			l >> atom1 >> atom2;
			rsd->add_cut_bond( atom1, atom2 );

		} else if ( tag == "CHI" ) {
			Size chino;
			l >> chino >> atom1 >> atom2 >> atom3 >> atom4;
			rsd->add_chi( chino, atom1, atom2, atom3, atom4 );
		} else if ( tag == "PROTON_CHI") {
			Size chino, nsamples, nextra_samples;
			std::string dummy;
			l >> chino;
			l >> dummy; // should be "SAMPLES"
			l >> nsamples;
			utility::vector1< Real > samples( nsamples );
			for ( Size ii = 1; ii <= nsamples; ++ii ) {
				l >> samples[ ii ];
			}
			l >> dummy; // should be "EXTRA"
			l >> nextra_samples;
			utility::vector1< Real > extra_samples( nextra_samples );
			for ( Size ii = 1; ii <= nextra_samples; ++ii ) {
				l >> extra_samples[ ii ];
			}
			rsd->set_proton_chi( chino, samples, extra_samples );

		} else if ( tag == "NBR_ATOM" ) {
			l >> atom1;
			rsd->nbr_atom( atom1 );

		} else if ( tag == "NBR_RADIUS" ) {
			Real radius;
			l >> radius;
			rsd->nbr_radius( radius );

		} else if ( tag == "PROPERTIES" ) {
			l >> tag;
			while ( !l.fail() ) {
				rsd->add_property( tag );
				l >> tag;
			}

		} else if ( tag == "VARIANT" ) {
			l >> tag;
			while ( !l.fail() ) {
				rsd->add_variant_type( tag );
				l >> tag;
			}
		} else if ( tag == "FIRST_SIDECHAIN_ATOM" ) {
			// note-- atoms are sidechain by default
			l >> tag;
			if ( tag == "NONE" ) {
				// set all atoms to backbone
				for ( Size j=1; j<= rsd->natoms(); ++j ) {
					rsd->set_backbone_heavyatom( rsd->atom_name(j) );
				}
			} else if ( rsd->has( tag ) ) {
				for ( Size j=1; j< rsd->atom_index( tag ); ++j ) {
					rsd->set_backbone_heavyatom( rsd->atom_name(j) );
				}
			}

		} else if ( tag == "IO_STRING" ) {
			assert( line.size() >= 15 );
			std::string const three_letter_code( line.substr(10,3) ),
							one_letter_code( line.substr(14,1) );
			rsd->name3( three_letter_code );
			rsd->name1( one_letter_code[0] );

		} else if ( tag == "AA" ) {
			l >> tag;
			rsd->aa( tag );
			found_AA_record = true;

		} else if ( tag == "NAME" ) {
			l >> tag;
			rsd->name( tag );

		} else if ( tag == "CHI_ROTAMERS" ) {
			Size chino;
			Real mean, sdev;
			l >> chino;
			l >> mean >> sdev;
			while ( !l.fail() ) {
				rsd->add_chi_rotamer( chino, mean, sdev );
				l >> mean >> sdev;
			}
		} else if ( tag == "PDB_ROTAMERS" ) {
			found_PDB_ROTAMERS_record = true;
			l >> pdb_rotamers_filename;
		} else if ( tag == "ACT_COORD_ATOMS" )
		{
			while ( l ) {
				l >> atom1;
				if ( atom1 == "END") break;
				rsd->add_actcoord_atom( atom1 );
			}
		} else if ( tag == "LOWER_CONNECT" ) {
			l >> atom1;
			rsd->set_lower_connect_atom( atom1 );
		} else if ( tag == "UPPER_CONNECT" ) {
			l >> atom1;
			rsd->set_upper_connect_atom( atom1 );
		} else if ( tag == "ADDUCT" ) {
			std::string adduct_name, adduct_atom_name, adduct_atom_type, adduct_mm_type;
			Real adduct_q, adduct_d, adduct_theta, adduct_phi;
			l >> adduct_name >> adduct_atom_name;
			l >> adduct_atom_type >> adduct_mm_type;
			l >> adduct_q >> adduct_phi >> adduct_theta >> adduct_d;
			l >> atom1 >> atom2 >> atom3;
			ObjexxFCL::lowercase(adduct_name);
			Adduct new_adduct( adduct_name, adduct_atom_name,
				adduct_atom_type, adduct_mm_type, adduct_q,
				adduct_phi, adduct_theta, adduct_d,
				atom1, atom2, atom3 );
			rsd->add_adduct( new_adduct );
		} else if ( tag == "NCAA_ROTLIB_PATH" ) {
			std::string path;
			l >> path;
			rsd->set_ncaa_rotlib_path( path );
			rsd->set_use_ncaa_rotlib( true );
		} else if ( tag == "NCAA_ROTLIB_NUM_ROTAMER_BINS" ) {
			Size n_rots(0);
			utility::vector1<Size> n_bins_per_rot;
			l >> n_rots;
			rsd->set_ncaa_rotlib_n_rotameric_bins( n_rots );
			n_bins_per_rot.resize( n_rots );
			for( Size i = 1; i <= n_rots; ++i ) {
				Size bin_size(0);
				l >> bin_size;
				n_bins_per_rot[i] = bin_size;
			}
			rsd->set_ncaa_rotlib_n_bin_per_rot( n_bins_per_rot );
		} else if( tag== "ORBITALS" ){ //begin parsing orbital information
			if(basic::options::option[ basic::options::OptionKeys::in::add_orbitals]){


				l >> orbitals_tag; //looking at the second set of text of the params file
				if(orbitals_tag == "BOND"){ //add pseudo bonds that will be used for iccor
					l >> atom1 >> orbital;
					rsd->add_orbital_bond( atom1, orbital );
				} else if( orbitals_tag == "ICOOR_INTERNAL"){
					Real phi, theta, d;
					std::string child_atom(""),  parent_atom(""), angle_atom(""), torsion_atom("");
					l >> child_atom >> phi >> theta >> d >> parent_atom >> angle_atom >> torsion_atom;

					phi = radians(phi); theta = radians(theta);
					rsd->set_orbital_icoor_id(child_atom, phi, theta, d, parent_atom, angle_atom, torsion_atom);
				}
				else { //assign the name of the orbital and the orbital type. This actually happens first
					std::string orbital_type_name("");
					//std::string orbital_name(l);
					//

					l >> orbital_type_name; //the orbital type, which is defined in OrbitalType.hh

					rsd->add_orbital( orbitals_tag, orbital_type_name); //orbitals tag is the name of the orbital
					++norbitals;
				}
			}

		}

	} // i=1,nlines

	//rsd->print_bonded_orbitals(); check to see if orbitals are being bonded


	if ( !found_AA_record ) {
		basic::Warning() << "No AA record found for " << rsd->name() << "; assuming " << name_from_aa( rsd->aa() ) << std::endl;
	}


	// set icoor coordinates, store information about polymer links
	// also sets up base_atom
	{

		std::map< std::string, Vector > rsd_xyz;

		for ( Size i=1; i<= nlines; ++i ) {

			std::string const & line( lines[i] );
			std::istringstream l( line );
			std::string tag, child_atom, parent_atom, angle_atom, torsion_atom;

			Real phi, theta, d;
			l >> tag;

			if ( tag != "ICOOR_INTERNAL" ) continue;

			l >> child_atom >> phi >> theta >> d >> parent_atom >> angle_atom >> torsion_atom;

			phi = radians(phi); theta = radians(theta); // in degrees in the file for human readability


			if ( natoms > 1 ) {
				/// build the cartesian coords for the new atom:
				if ( child_atom == parent_atom ) {
					assert( rsd_xyz.empty() ); // first atom
					rsd_xyz[ child_atom ] = Vector( 0.0 );

				} else if ( child_atom == angle_atom ) {
					assert( rsd_xyz.size() == 1 && rsd_xyz.count( parent_atom ) ); // second atom
					rsd_xyz[ child_atom ] = Vector( d, 0.0, 0.0 );

				} else {
					Vector torsion_xyz;
					if ( child_atom == torsion_atom ) {
						assert( rsd_xyz.size() == 2 );
						assert( rsd_xyz.count( parent_atom ) );
						assert( rsd_xyz.count( angle_atom ) ); // third atom
						torsion_xyz = Vector( 1.0, 1.0, 0.0 );
					} else {
						assert( rsd_xyz.count( parent_atom ) && rsd_xyz.count( angle_atom ) && rsd_xyz.count( torsion_atom ) );
						torsion_xyz = rsd_xyz[ torsion_atom ];
					}
					kinematics::Stub const stub( rsd_xyz[ parent_atom ], rsd_xyz[ angle_atom ], torsion_xyz );
					rsd_xyz[ child_atom ] = stub.spherical( phi, theta, d );
				}
			}


			// set atom_base
			if ( child_atom != "UPPER" && child_atom != "LOWER" && child_atom.substr(0,4) != "CONN" ) {
				// atom base only valid for genuine atoms of this residue
				if ( child_atom == parent_atom ) {
					// root of the tree
					if ( natoms == 1 ) {
						rsd->set_atom_base( child_atom, child_atom ); // 1st child of root atom
					} else {
						rsd->set_atom_base( child_atom, angle_atom ); // 1st child of root atom
					}
				} else {
					rsd->set_atom_base( child_atom, parent_atom );
				}
			}

			// set icoor
			rsd->set_icoor( child_atom, phi, theta, d, parent_atom, angle_atom, torsion_atom );

		} // loop over file lines looking for ICOOR_INTERNAL lines


		// fill in the rsd-xyz values
		if ( natoms == 1 ) {
			std::string const name( rsd->atom_name(1) );
			rsd->set_xyz( name, Vector(0.0) );

		} else {
			// now fill in the icoor values -- in principle the rsd itself could be doing this...
			for ( Size i=1; i<= natoms; ++i ) {
				std::string name( rsd->atom_name(i) );
				strip_whitespace( name );
				assert( rsd_xyz.count( name ) );
				rsd->set_xyz( name, rsd_xyz[ name ] );
				//rsd->set_xyz( rsd->atom_name(i), atom_tree.xyz( id::AtomID(i,1) ) );
				//rsd->atom(i).xyz( atom_tree.xyz( id::AtomID(i,1) ) );
			}
		}


		// if polymer, fill in upper/lower connect and mainchain info
		if ( rsd->is_polymer() ) {
			uint upper_connect( rsd->upper_connect_atom() ), lower_connect( rsd->lower_connect_atom() );

			// fill in the mainchain info -- shortest path between upper connect and lower connect
			if ( upper_connect && lower_connect ) {
				AtomIndices mainchain;
				FArray2D_int D( get_residue_path_distances( *rsd ) );
				uint atom( lower_connect );
				while ( atom != upper_connect ) {
					mainchain.push_back( atom );
					AtomIndices const & nbrs( rsd->nbrs( atom ) );
					int min_d( D( atom, upper_connect ) );
					uint next_atom( atom );
					//std::cout << "setup mainchain: " << rsd->name() << ' ' << rsd->atom_name( atom ) << ' ' <<
					//	min_d << std::endl;

					for ( uint i=1; i<= nbrs.size(); ++i ) {
						uint const nbr( nbrs[i] );
						if ( D( nbr, upper_connect ) < min_d ) {
							min_d = D( nbr, upper_connect );
							next_atom = nbr;
						}
					}
					assert( next_atom != atom );

					atom = next_atom;
				}
				mainchain.push_back( upper_connect );
				rsd->set_mainchain_atoms( mainchain );
			}
		}

		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////


		// now also need to store the information about the geometry
		// at the links...

	} // scope

// 			AtomID
// 				child_atom_id  ( atom_id_from_icoor_line(   child_atom, *rsd ) ),
// 				parent_atom_id ( atom_id_from_icoor_line(  parent_atom, *rsd ) ),
// 				angle_atom_id  ( atom_id_from_icoor_line(   angle_atom, *rsd ) ),
// 				torsion_atom_id( atom_id_from_icoor_line( torsion_atom, *rsd ) );

// 			bool child_is_bonded_atom( true );
// 			if ( atom_tree.has( parent_atom_id ) ) {
// 				// not the root of the tree
// 				atom_tree.add_atom( child_atom_id, parent_atom_id, child_is_bonded_atom, false /* by xyz */ );
// 			} else {
// 				// root of the tree
// 				child_is_bonded_atom = false;
// 				atom_tree.add_atom( child_atom_id, id::BOGUS_ATOM_ID, child_is_bonded_atom, false /* by xyz */);
// 			}

// 			if ( child_is_bonded_atom ) {
// 				// set the internal coordinates in the atom_tree
// 				atom_tree.set_dof( DOF_ID( child_atom_id, id::THETA ), theta );
// 				atom_tree.set_dof( DOF_ID( child_atom_id, id::D     ),     d );

// 				if ( ( child_atom_id == angle_atom_id ) || ( child_atom_id == torsion_atom_id ) ) {
// 					// one of the first three atoms in the tree, some internal coords not defined
// 					assert( std::abs( phi ) < 1e-3 );
// 					atom_tree.set_dof( DOF_ID( child_atom_id, id::PHI ), phi );
// 				} else if ( ( parent_atom_id ==   angle_atom_id ) ||
// 										( parent_atom_id == torsion_atom_id ) ||
// 										( angle_atom_id == torsion_atom_id ) ) {
// 					// this would not be good for set_torsion_angle
// 					// should be special case for single-atom residue
// 					// we are not going to use the xyz's from the atomtree anyhow, so we dont need to do anything
// 					assert( natoms == 1 );
// 				} else {
// 					atom_tree.set_dof( DOF_ID( child_atom_id, id::PHI ), 0.0 /*temporary value*/ );
// 					atom_tree.set_torsion_angle( child_atom_id, parent_atom_id, angle_atom_id, torsion_atom_id, phi );
// 				}
// 			} // is the child atom a bonded atom, ie is it not the root of the tree?



// 	// set icoor coordinates, store information about polymer links
// 	// also sets up base_atom
// 	{
// 		Size root_atomno, anchor_atomno;

// 		bool first_line( true );
// 		for ( Size i=1; i<= nlines; ++i ) {
// 			std::string const & line( lines[i] );
// 			std::istringstream l( line );
// 			std::string tag, child_atom, parent_atom;
// 			Real phi, theta, d;
// 			l >> tag;
// 			if ( tag != "ICOOR_INTERNAL" ) continue;
// 			l >> child_atom >> parent_atom >> phi >> theta >> d;
// 			phi = radians(phi); theta = radians(theta);
// 			// 		if ( parent_atom[0] == '-' ) {
// 			// 			assert( rsd->is_polymer() );
// 			// 			// anchor atom in previous residue
// 			// 			continue;
// 			// 		}
// 			Size child_rsd(1), parent_rsd(1);
// 			if ( parent_atom[0] == '-' ) {
// 				continue; // do nothing, temporarily
// // 				parent_rsd = 0;
// // 				parent_atom.erase(0,1);
// 			} else if ( child_atom[0] == '+' ) {
// 				child_rsd = 2;
// 				child_atom.erase(0,1);
// 			} else {
// 				rsd->set_atom_base( child_atom, parent_atom );
// 				if ( first_line ) {
// 					rsd->set_atom_base( parent_atom, child_atom );
// 				}
// 			}

// 			AtomID
// 				child_id ( rsd->atom_index( child_atom ), child_rsd ),
// 				parent_id( rsd->atom_index( parent_atom ), parent_rsd );

// 			if ( first_line ) {
// 				atom_tree.add_atom( parent_id, id::BOGUS_ATOM_ID, false,
// 																	 false /* by xyz */ );
// 				root_atomno = parent_id.atomno();
// 				first_line = false;
// 			} else if ( child_rsd == 2 ) {
// 				anchor_atomno = parent_id.atomno();
// 			}

// 			atom_tree.add_atom( child_id, parent_id, true /* bonded */,
// 													false /* by xyz */ );
// 			atom_tree.set_dof( DOF_ID( child_id, id::PHI   ),   phi );
// 			atom_tree.set_dof( DOF_ID( child_id, id::THETA ), theta );
// 			atom_tree.set_dof( DOF_ID( child_id, id::D     ),     d );
// 		}

// 		if ( rsd->is_polymer() ) {
// 			// fill in the mainchain info
// 			AtomIndices mainchain;
// 			mainchain.push_back( anchor_atomno );
// 			while ( true ) {
// 				anchor_atomno = rsd->atom_base( anchor_atomno );
// 				if ( std::find( mainchain.begin(), mainchain.end(), anchor_atomno ) !=
// 						 mainchain.end() ) {
// 					utility_exit_with_message("failure deriving mainchain from atom_base");
// 				}
// 				mainchain.insert( mainchain.begin(), anchor_atomno );
// 				if ( anchor_atomno == root_atomno ) break;
// 			}
// 			rsd->set_mainchain_atoms( mainchain );
// 		}


// 		// fill in the icoor values
// 		for ( Size i=1; i<= natoms; ++i ) {
// 			rsd->atom(i).xyz( atom_tree.xyz( id::AtomID(i,1) ) );
// 		}

// 		// now also need to store the information about the geometry
// 		// at the links...
// 	}


	// calculate any remaining derived data
	rsd->finalize();

	// If we found a PDB_ROTAMERS library, load them now that the ResidueType
	// is totally initialized:
	if( found_PDB_ROTAMERS_record ) {
		using namespace utility::file;
		// Assume name of rotamers file has no path info, etc
		// and is in same directory as the residue parameters file.
		FileName this_file( filename ), rot_file( pdb_rotamers_filename );
		rot_file.vol( this_file.vol() );
		rot_file.path( this_file.path() );

		rsd->set_RotamerLibraryName( rot_file() );

		tr.Debug << "Setting up conformer library for " << rsd->name() << std::endl;
		/*using namespace core::pack::dunbrack;
		using namespace utility::file;
		SingleLigandRotamerLibraryOP pdb_rotamers = new SingleLigandRotamerLibrary();
		// Assume name of rotamers file has no path info, etc
		// and is in same directory as the residue parameters file.

		pdb_rotamers->init_from_file( rot_file.name(), rsd );
		rsd->set_RotamerLibrary( pdb_rotamers );*/
	}

	return rsd;
}


/// @brief function to write out a topology file given a residue type, can be used to
/// @brief debug on the fly generated residue types. Note: not perfect yet, the enums for
/// @brief the connection types are given in numbers instead of names
void
write_topology_file(
	ResidueType const & rsd
)
{

	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	std::string filename = rsd.name() + ".params";

	std::ofstream out( filename.c_str() );

	out << "#rosetta residue topology file \n";
	out << "#version ??? \n";
	out << "#This automatically generated file is not really formatted, but should work, excpet that enums for connection types \n# are given in numbers and not strings. \n";

	//first write out all the general tags
	out << "NAME " << rsd.name() << " \n";
	out << "IO_STRING " << rsd.name3() << " " << rsd.name1() << " \n";
	if( rsd.is_polymer() ) { out << "TYPE POLYMER \n"; }
	else if ( rsd.is_ligand() ) { out << "TYPE LIGAND \n"; }
  else if ( rsd.is_surface() ) { out << "TYPE SURFACE \n"; }
	out << "AA " << rsd.aa() << " \n";

	//then write out the atoms
	for(Size i=1; i <= rsd.natoms(); i++){

		std::string atom_out = "ATOM " + rsd.atom_name( i ) + " " + rsd.atom_type( i ).name() + "  ";
		atom_out = atom_out + rsd.mm_atom_type(i).name();
		out << atom_out << " " << rsd.atomic_charge(i) << " \n";

	} //atom write out

	if( rsd.is_polymer() ) {
		if( !rsd.is_lower_terminus() ) { out << "LOWER_CONNECT " << rsd.atom_name( rsd.lower_connect().atomno() ) << " \n"; }
		if( !rsd.is_upper_terminus() ) { out << "UPPER_CONNECT " << rsd.atom_name( rsd.upper_connect().atomno() ) << " \n";}
	}

	//then all the bonds
	for(Size i=1; i <= rsd.natoms(); i++){
		foreach(Size atom_index, rsd.nbrs(i)){// bond_this_atom
			if( atom_index > i ) {  //don't write out bonds more than once
				out << "BOND  " << rsd.atom_name( i ) << "    " << rsd.atom_name( atom_index ) << " \n";
			}
		}
	} // bond write out


	//now the chis
	for(Size i=1; i <= rsd.nchi(); i++){
		out << "CHI " << i ;
		AtomIndices atoms_this_chi = rsd.chi_atoms( i );
		for(AtomIndices::iterator at_it = atoms_this_chi.begin(); at_it != atoms_this_chi.end(); at_it++ ){
			out << "   " << rsd.atom_name( *at_it );
		}
		out << " \n";
	} //chi write out

	//and now the proton chis
	Size n_proton_chi(0);
	for(Size i=1; i <= rsd.nchi(); i++){
		if( rsd.is_proton_chi( i ) ){

			n_proton_chi++;
			out << "PROTON_CHI " << i << " SAMPLES ";
			utility::vector1< Real > pchi_samples = rsd.proton_chi_samples( n_proton_chi );
			utility::vector1< Real > pchi_extra = rsd.proton_chi_extra_samples( n_proton_chi );
			out << pchi_samples.size() ;
			for( Size j = 1; j <= pchi_samples.size(); j++){ out << " " << pchi_samples[j]; }
			out << " EXTRA " << pchi_extra.size();
			for( Size j = 1; j <= pchi_extra.size(); j++){ out << " " << pchi_extra[j]; }
			out << " \n";

		}
	}//proton chi write out

	//now all the properties
	out << "PROPERTIES";
	if(rsd.is_protein() ) { out << " PROTEIN"; }
	if(rsd.is_DNA() ) { out << " DNA"; }
	if(rsd.is_RNA() ) { out << " RNA"; }
	if(rsd.is_polar() ) { out << " POLAR"; }
	if(rsd.is_charged() ) { out << " CHARGED"; }
	if(rsd.is_aromatic() ) { out << " AROMATIC"; }
	if(rsd.is_lower_terminus() ) { out << " LOWER_TERMINUS"; }
	if(rsd.is_upper_terminus() ) { out << " UPPER_TERMINUS"; }
	if(rsd.is_terminus() ) { out << " TERMINUS"; }
	out << " \n";

	out << "NBR_ATOM " << rsd.atom_name( rsd.nbr_atom() ) << " \n";
	out << "NBR_RADIUS " << rsd.nbr_radius() << " \n";


	//actcoord atoms
	if( rsd.actcoord_atoms().size() > 0 ){
		out << "ACT_COORD_ATOMS ";
		AtomIndices act_atoms = rsd.actcoord_atoms();
		for(AtomIndices::iterator at_it = act_atoms.begin(); at_it != act_atoms.end(); at_it++ ){
			out << rsd.atom_name( *at_it ) << " ";
		}
		out << "END \n";
	}


	//last but not least the internal coordinates
	for(Size i=1; i <= rsd.natoms(); i++){
		AtomICoor cur_icoor = rsd.icoor( i );
		out << "ICOOR_INTERNAL   " << rsd.atom_name( i ) << "  " << degrees( cur_icoor.phi() ) << "  ";
		out << degrees( cur_icoor.theta() ) << "  " << cur_icoor.d();
		if( ( cur_icoor.stub_atom1().atomno() <= rsd.natoms() ) && ( cur_icoor.stub_atom1().atomno() > 0 ) ) {
			out << "   " << rsd.atom_name( cur_icoor.stub_atom1().atomno() );
		}
		else{ out << "   " << cur_icoor.stub_atom1().type(); }

		if( ( cur_icoor.stub_atom2().atomno() <= rsd.natoms()  ) && ( cur_icoor.stub_atom2().atomno() > 0 )){
			out << "   " << rsd.atom_name( cur_icoor.stub_atom2().atomno() );
		}
		else{ out << "  "  << cur_icoor.stub_atom2().type();}

		if( ( cur_icoor.stub_atom3().atomno() <= rsd.natoms()  ) && ( cur_icoor.stub_atom3().atomno() > 0 )){
			out << "   " << rsd.atom_name( cur_icoor.stub_atom3().atomno() );
		}
		else{ out << "  "  << cur_icoor.stub_atom3().type() ;}

		out << " \n";

	} //atom icoor write out

	//now write out icoors for connections (polymer, other)

	//TO DO

	out.close();

} //write_topology_file





} // chemical
} // core
