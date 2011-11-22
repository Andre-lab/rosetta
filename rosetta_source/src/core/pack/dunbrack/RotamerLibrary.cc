
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
/// @author


// Unit headers
#include <core/pack/dunbrack/RotamerLibrary.hh>

// Package headers
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.tmpl.hh>
#include <core/pack/dunbrack/SingleLigandRotamerLibrary.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>

// Project headers
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/graph/Graph.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>

#include <basic/basic.hh>
#include <basic/interpolate.hh>
//XRW_B_T1
//#include <core/coarse/Translator.hh>
//#include <core/coarse/TranslatorSet.hh>
//XRW_E_T1
#include <core/pose/Pose.hh>
#include <core/pack/task/PackerTask.hh>
//#include <core/scoring/TenANeighborGraph.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <numeric/constants.hh>
// AUTO-REMOVED #include <numeric/angle.functions.hh>

#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/format.hh>

// C++ Headers
#include <string>
#include <iostream>
#include <fstream>
#if defined(WIN32) || defined(__CYGWIN__)
// AUTO-REMOVED #include <ctime>
#endif

// Boost Headers
#include <boost/cstdint.hpp>

// C Headers
// AUTO-REMOVED #include <sys/stat.h>
// AUTO-REMOVED #include <stdio.h>
// AUTO-REMOVED #include <stdlib.h>

#include <basic/Tracer.hh>
using basic::T;

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/vector1.hh>



using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.pack.dunbrack");

/**
	 rearrange the dunbrack arrays:

	 mapping from aa to single-residue-rotamer-library

	 to preserve current interpolation scheme, need ability to (quickly)
	 lookup the probability/chi-means/chi-sdevs for a given rotamer in
	 several torsion-angle bins.

	 so how to store prob/mean/sd for alt rotamers at different pp bins

	 quick conversion from rotno-tuple to index?

	 rotamer_data: probability, chi-mean, chi-sd

**/
namespace core {
namespace pack {
namespace dunbrack {

/// @brief helper that'll determine the rotamer bins for a residue by asking its associated
/// rotamer library for that information.
void
rotamer_from_chi(
	conformation::Residue const & rsd,
	RotVector & rot
)
{
	assert( rsd.aa() <= chemical::num_canonical_aas );

	SingleResidueRotamerLibraryCAP rotlib = RotamerLibrary::get_instance().get_rsd_library( rsd.type() );
	if ( rotlib ) {
		SingleResidueDunbrackLibraryCAP dunlib( static_cast< SingleResidueDunbrackLibrary const * > ( rotlib() ));
		dunlib->get_rotamer_from_chi( rsd.chi(), rot );
	} else { rot.resize( 0 ); }
}

/// @brief helper that'll determine the rotamer bins for a residue_type by asking its associated
/// rotamer library for that information. Takes chi vector as input parameter.
void
rotamer_from_chi(
	chemical::ResidueType const & rsd_type,
	ChiVector const & chi,
	RotVector & rot
)
{
	SingleResidueRotamerLibraryCAP rotlib = RotamerLibrary::get_instance().get_rsd_library( rsd_type );
	if ( rotlib ) {
		SingleResidueDunbrackLibraryCAP dunlib( static_cast< SingleResidueDunbrackLibrary const * > ( rotlib() ));
		dunlib->get_rotamer_from_chi( chi, rot );
	} else { rot.resize( 0 ); }
}


void
rotamer_from_chi_02(
	ChiVector const & chi, // pass by reference for speed
	chemical::AA const res,
	RotVector & rot
)
{
	rot.resize( chi.size() );

	Size4 rot4; Real4 chi4;
	Size nchi( std::min( chi.size(), DUNBRACK_MAX_SCTOR ) );
	for ( Size ii = 1; ii <= nchi; ++ii ) chi4[ ii ] = chi[ ii ];

	rotamer_from_chi_02( chi4, res, nchi, rot4 );

	for ( Size ii = 1; ii <= nchi; ++ii ) rot[ ii ] = rot4[ ii ];
	for ( Size ii = nchi+1; ii <= rot.size(); ++ii ) rot[ ii ] = 0;
}

/// @brief DEPRICATED
/// convert between the real-valued chi dihedrals and the rotamer well indices.
///
/// @details This code comes directly from Roland Dunbrack.  In certain rare edge cases,
/// floating point comparison breaks down (e.g. ! x >= 120.0 && ! x < 120.0 ) and
/// Dunbrack's code leaves the rotamer well unassigned -- a value of zero. I'm modifying
/// the code to guarantee that the rotamer well is assigned, but so that
/// Dunbrack's original form is still recognizable.
/// e.g.
/// if (x) a; if (y) b; if (z) c; now reads as:
/// if (x) a; else if ( y ) b; else /*if (z)*/ c;
void
rotamer_from_chi_02(
	Real4 const & chi,
	chemical::AA const res,
	Size nchi,
	Size4 & rot
)
{
	using namespace chemical;

	// default to 0
	// right now this means 0 for rot numbers larger than dunbrack's nchi
	// eg for H-chi or extra chi's

	for ( Size i = 1; i <= nchi; ++i ) {
		rot[i] = 0;
		Real chi_i = basic::periodic_range( chi[ i ], 360.0 ); // this saves us a vector construction/destruction
		if ( i == 1 ) { // chi 1 //KMa phospho_ser
			if ( res != aa_pro && res != aa_gly && res != aa_ala ) {
				// chi1 of all residues except P,G,A,RNA
				if ( ( chi_i >= 0.0 ) && ( chi_i <= 120.0 ) ) { rot[i] = 1; }
				else if ( std::abs(chi_i) >= 120.0 ) { rot[i] = 2; }
				else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/{ rot[i] = 3; }
			} else if ( res == aa_pro ) {
				// chi1 of pro
				if ( chi_i >= 0.0 ) { rot[i] = 1; }
				else /*if ( chi_i <= 0.0 )*/ { rot[i] = 2; }
			}

		} else if ( i == 2 ) { // chi 2
			if ( res == aa_arg ) {
				if ( ( chi_i >= 0.0 ) && ( chi_i <= 120.0 ) ) { rot[i] = 1; }
				else if ( std::abs(chi_i) >= 120.0 ) { rot[i] = 2; }
				else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/ { rot[i] = 3; }

			} else if ( res == aa_glu || res == aa_his || res == aa_ile ||
									res == aa_lys || res == aa_leu || res == aa_met ||
									res == aa_gln ) {
				if ( ( chi_i >= 0.0 ) && ( chi_i <= 120.0 ) ) { rot[i] = 1; }
				else if ( std::abs(chi_i) >= 120.0 ) { rot[i] = 2; }
				else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/ { rot[i] = 3; }

			} else if ( res == aa_asp ) {
				if ( ( ( chi_i >= 30.0 ) && ( chi_i <= 90.0 ) ) ||
						( ( chi_i <= -90.0 ) && ( chi_i >= -150.0 ) ) ) {
					rot[i] = 1;
				} else if ( ( ( chi_i >= -30.0 ) && ( chi_i <= 30.0 ) ) ||
						( std::abs(chi_i) >= 150.0 ) ) {
					rot[i] = 2;
				} else {
				/*if ( ( ( chi_i >= -90.0 ) && ( chi_i <= -30.0 ) ) ||
						( ( chi_i >= 90.0 ) && ( chi_i <= 150.0 ) ) ) */
					rot[i] = 3;
				}

			} else if ( ( res == aa_phe ) || ( res == aa_tyr ) ) {
				if ( ( ( chi_i >= 30.0 ) && ( chi_i <= 150.0 ) ) ||
						( ( chi_i <= -30.0 ) && ( chi_i >= -150.0 ) ) ) {
					rot[i] = 1;
				} else { /*if ( ( ( chi_i >= -30.0 ) && ( chi_i <= 30.0 ) ) ||
						 ( std::abs(chi_i) >= 150.0 ) ) */
					rot[i] = 2;
				}

			} else if ( res == aa_trp ) {
				if ( ( chi_i >= -180.0 ) && ( chi_i <= -60.0 ) ) { rot[i] = 1; }
				else if ( ( chi_i >= -60.0 ) && ( chi_i <= 60.0 ) ) { rot[i] = 2; }
				else /*if ( ( chi_i >= 60.0 ) && ( chi_i <= 180.0 ) ) */ { rot[i] = 3; }

			} else if ( res == aa_asn ) { // chi2 of asn
				// ctsa - note this is a special case of the new dunbrack rotamer set
				if ( rot[1] == 1 ) {
					if ( chi_i >= -150.0 && chi_i < -90.0 ) { rot[i] = 1; }
					else if ( chi_i >= -90.0 && chi_i < -30.0 ) { rot[i] = 2; }
					else if ( chi_i >= -30.0 && chi_i < 30.0 ) { rot[i] = 3; }
					else if ( chi_i >=  30.0 && chi_i < 90.0 ) { rot[i] = 4; }
					else if ( chi_i >=  90.0 && chi_i < 150.0 ) { rot[i] = 5; }
					else /*if ( std::abs(chi_i) >= 150.0 )*/ { rot[i] = 6; }
				} else if ( rot[1] == 2 ) {
					if ( chi_i >= -180.0 && chi_i < -90.0 ) { rot[i] = 1; }
					else if ( chi_i >= -90.0 && chi_i < -45.0 ) { rot[i] = 2; }
					else if ( chi_i >= -45.0 && chi_i <  0.0 ) { rot[i] = 3; }
					else if ( chi_i >=   0.0 && chi_i < 45.0 ) { rot[i] = 4; }
					else if ( chi_i >=  45.0 && chi_i < 90.0 ) { rot[i] = 5; }
					else /*if ( chi_i >=  90.0 && chi_i <= 180.0 )*/ { rot[i] = 6; }
				} else {
					if ( chi_i >= -180.0 && chi_i < -105.0 ) { rot[i] = 1; }
					else if ( chi_i >= -105.0 && chi_i < -45.0 ) { rot[i] = 2; }
					else if ( chi_i >= -45.0 && chi_i <  15.0 ) { rot[i] = 3; }
					else if ( chi_i >=  15.0 && chi_i <  60.0 ) { rot[i] = 4; }
					else if ( chi_i >=  60.0 && chi_i < 120.0 ) { rot[i] = 5; }
					else /*if ( chi_i >= 120.0 && chi_i <= 180.0 )*/ { rot[i] = 6; }
				}
			} else if ( res == aa_pro ) {
				// apl proline has only two rotamers, distinguishable by chi1.  Chi2 is always 1.
				rot[i]=1;
			}


		} else if ( i == 3 ) { // chi 3
			if ( res == aa_glu ) {
				if ( ( ( chi_i >= 30.0 ) && ( chi_i <= 90.0 ) ) ||
					( ( chi_i <= -90.0 ) && ( chi_i >= -150.0 ) ) ) { rot[i] = 1; }
				else if ( ( ( chi_i >= -30.0 ) && ( chi_i <= 30.0 ) ) ||
					( std::abs(chi_i) >= 150.0 ) ) { rot[i] = 2; }
				else /*if ( ( ( chi_i >= -90.0 ) && ( chi_i <= -30.0 ) ) ||
					( ( chi_i >= 90.0 ) && ( chi_i <= 150.0 ) ) )*/ { rot[i] = 3; }

			} else if ( res == aa_arg || res == aa_lys || res == aa_met ) {
				if ( ( chi_i >= 0.0 ) && ( chi_i <= 120.0 ) ) { rot[i] = 1; }
				else if ( std::abs(chi_i) > 120.0 ) { rot[i] = 2; }
				else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/ { rot[i] = 3; }

			} else if ( res == aa_gln ) { // chi3 of gln
				// ctsa - note this is a special case of the new dunbrack rotamer set
				if ( rot[2] == 2 ) {
					if ( chi_i >= 135.0 || chi_i < -135.0 ) { rot[i] = 1; }
					else if ( chi_i >= -135.0 && chi_i < -45.0 ) { rot[i] = 2; }
					else if ( chi_i >= -45.0 && chi_i <  45.0 ) { rot[i] = 3; }
					else /*if ( chi_i >=  45.0 && chi_i < 135.0 )*/ { rot[i] = 4; }
				} else {
					if ( chi_i >= -180.0 && chi_i < -90.0 ) { rot[i] = 1; }
					else if ( chi_i >= -90.0 && chi_i <   0.0 ) { rot[i] = 2; }
					else if ( chi_i >=   0.0 && chi_i <  90.0 ) { rot[i] = 3; }
					else /*if ( chi_i >=  90.0 && chi_i <= 180.0 )*/ { rot[i] = 4; }
				}
			} else if ( res == aa_pro ) {
				// apl
				rot[i]=1;
			}
		} else if ( i == 4 ) {  // chi 4
			if ( res == aa_arg || res == aa_lys ) {
				if ( ( chi_i >= 0.0 ) && ( chi_i <= 120.0 ) ) { rot[i] = 1; }
				else if ( std::abs(chi_i) > 120.0 ) { rot[i] = 2; }
				else /*if ( ( chi_i >= -120.0 ) && ( chi_i <= 0.0 ) )*/ { rot[i] = 3; }
			}
		}
	} // i=1, i<= nchi

	//std::cout << "Rotamer from chi02: " << res << " ";
	//for ( Size ii = 1; ii <= nchi; ++ii ) std::cout << chi[ii] << " ";
	//std::cout << " : ";
	//for ( Size ii = 1; ii <= nchi; ++ii ) std::cout << rot[ii] << " ";
	//std::cout << std::endl;


}

///////////////////////////////////////////////////////////////////////////////
// taken from rotamer_functions -- total temporary hack
//
// removed pser and rna
//


SingleResidueRotamerLibrary::~SingleResidueRotamerLibrary()
{}


//// RotamerLibrary
RotamerLibrary * RotamerLibrary::rotamer_library_( 0 );

RotamerLibrary::RotamerLibrary():
	reslibraries_(),
	libraries_(),
	aa_libraries_( chemical::num_canonical_aas, 0 ), // for the canonical aa's, fast access
	libraries_ops_()
{}

RotamerLibrary::~RotamerLibrary()
{}

Real
RotamerLibrary::rotamer_energy(
	Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	SingleResidueRotamerLibraryCAP const library( get_rsd_library( rsd.type() ) );

	if ( library ) {
		return library->rotamer_energy( rsd, scratch );
	}
	return 0.0;
}

Real
RotamerLibrary::rotamer_energy_deriv(
	Residue const & rsd,
	RotamerLibraryScratchSpace & scratch
) const
{
	SingleResidueRotamerLibraryCAP const library( get_rsd_library( rsd.type() ) );

	if ( library ) {
		return library->rotamer_energy_deriv( rsd, scratch );
	}

	Real3 & dE_dbb( scratch.dE_dbb() );
	Real4 & dE_dchi( scratch.dE_dchi() );

	// ensure that these guys are dimensioned
	//dE_dbb.resize ( rsd.mainchain_torsions().size() );
	//dE_dchi.resize( rsd.nchi() );
	std::fill( dE_dbb.begin(), dE_dbb.end(), 0 );
	std::fill( dE_dchi.begin(), dE_dchi.end(), 0 );
	return 0.0;

}


//XRW_B_T1
/*
RotamerLibraryOP RotamerLibrary::coarsify(
	coarse::TranslatorSet const & map_set
) const
{
	RotamerLibraryOP coarse_lib = new RotamerLibrary();
	for (library_iterator it=libraries_.begin(), eit=libraries_.end(); it!=eit; ++it) {
		chemical::AA aa = it->first;
		coarse::TranslatorCOP map=map_set.default_for_aa(aa);
		coarse_lib->add_residue_library(aa,it->second->coarsify(*map));
	}
	return coarse_lib;
}
*/
//XRW_E_T1

///
void
RotamerLibrary::add_residue_library(
	AA const & aa,
	SingleResidueRotamerLibraryCOP rot_lib
) const
{
	if ( libraries_.find(aa) != libraries_.end() ) {
		utility_exit_with_message("cant add rsd library twice");
	}
	libraries_.insert( std::make_pair( aa, rot_lib() ) );
	if ( aa <= chemical::num_canonical_aas ) {
		aa_libraries_[ aa ] = rot_lib(); // ASSUMPTION -- only fullatom rotamer libraries are added; centroid restypes don't have rotlibs.
	}
	libraries_ops_.push_back( rot_lib );
}

void
RotamerLibrary::add_residue_library(
	ResidueType const & rsd_type,
	SingleResidueRotamerLibraryCOP rot_lib
) const
{
	// May want to change the library for a specific residue type during the run...
	// This will not be threadsafe, though.
	//if ( reslibraries_.find(aa) != reslibraries_.end() ) {
	//	utility_exit_with_message("cant add rsd library twice");
	//}
	std::string const key( rsd_type.name()+"@"+rsd_type.residue_type_set().name() );
	reslibraries_.insert( std::make_pair( key, rot_lib() ) );
	reslibraries_.insert( std::make_pair( rsd_type.name(), rot_lib() ) );
	libraries_ops_.push_back( rot_lib );
}

/// @brief Returns true if the singleton has already loaded a library for a given residue type.
bool
RotamerLibrary::rsd_library_already_loaded( chemical::ResidueType const & rsd_type ) const
{
	if ( ncaa_rotlibs_.find( rsd_type.name3() ) != ncaa_rotlibs_.end() ) {
		return true;
	}

	if ( rsd_type.residue_type_set().name() == chemical::FA_STANDARD &&
			rsd_type.aa() <= chemical::num_canonical_aas ) {
		return aa_libraries_[ rsd_type.aa() ] != 0;
	}

	std::string const key( rsd_type.name()+"@"+rsd_type.residue_type_set().name() );
	if ( reslibraries_.find( key ) != reslibraries_.end() ) {
		return true;
	} else if ( reslibraries_.find( rsd_type.name() ) != reslibraries_.end() ) {
		return true;
	// Dunbrack library only applies to fullatom residue types, right?
	} else if ( rsd_type.residue_type_set().name() != chemical::CENTROID
			&& libraries_.find( rsd_type.aa() ) != libraries_.end() ) {
		return true;
	}

	return false;
}

/// @details returns NULL if no library is available for the given residue type
/// and no pdb rotamers have been defined for the given residue; if PDB rotamers HAVE
/// been defined, then
SingleResidueRotamerLibraryCAP
RotamerLibrary::get_rsd_library( chemical::ResidueType const & rsd_type ) const
{
	/// intercept, in case the residue type would prefer using the NCAARotamerLibrary
	if ( rsd_type.get_use_ncaa_rotlib() ) {
		return get_NCAARotamerLibrary( rsd_type );
	}

	if ( rsd_type.residue_type_set().name() == chemical::FA_STANDARD &&
			rsd_type.aa() <= chemical::num_canonical_aas ) {
		return aa_libraries_[ rsd_type.aa() ];
	}

	std::string const key( rsd_type.name()+"@"+rsd_type.residue_type_set().name() );
	if ( reslibraries_.find( key ) != reslibraries_.end() ) {
		return reslibraries_.find( key )->second;
	} else if ( reslibraries_.find( rsd_type.name() ) != reslibraries_.end() ) {
		return reslibraries_.find( rsd_type.name() )->second;
	// Dunbrack library only applies to fullatom residue types, right?
	} else if ( rsd_type.residue_type_set().name() != chemical::CENTROID
			&& libraries_.find( rsd_type.aa() ) != libraries_.end() ) {
		return libraries_.find( rsd_type.aa() )->second;
	}

	if ( rsd_type.is_ligand() && rsd_type.get_RotamerLibraryName() != "" ) {
		TR.Debug << "Initializing conformer library for " << rsd_type.get_RotamerLibraryName() << std::endl;
		SingleLigandRotamerLibraryOP pdb_rotamers = new SingleLigandRotamerLibrary();
		pdb_rotamers->init_from_file( rsd_type.get_RotamerLibraryName(), &rsd_type );
		add_residue_library( rsd_type, pdb_rotamers );
		return pdb_rotamers();
	}
	return 0;
}

SingleResidueRotamerLibraryCAP
RotamerLibrary::get_library_by_aa( chemical::AA const & aa ) const
{
	if ( (Size) aa <= aa_libraries_.size() ) {
		return aa_libraries_[ aa ];
	}
	return 0;
}


/// @details Generic decision tree for both 02 and 08 libraries.  Requires
/// that low-level subroutines dispatch to 02 and 08 versions based on
/// the dun08 flag.
void
RotamerLibrary::create_fa_dunbrack_libraries()
{

	if ( decide_read_from_binary() ) {
		create_fa_dunbrack_libraries_from_binary();
	} else {
		create_fa_dunbrack_libraries_from_ASCII();
		if ( decide_write_binary() ) write_binary_fa_dunbrack_libraries();
	}
}

bool
RotamerLibrary::decide_read_from_binary() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ in::file::no_binary_dunlib ] ) return false;

	std::string binary_filename = get_binary_name();
	utility::io::izstream binlib( binary_filename.c_str(), std::ios::in | std::ios::binary );

	if ( ! binlib ) {
		return false;
	}

	if ( ! binary_is_up_to_date( binlib ) ) return false;

	/// Look for specific objections for the 02 or 08 libraries (if any exist!)
	if ( option[ corrections::score::dun10 ] ) {
		if ( ! decide_read_from_binary_10() ) return false;
	} else if ( option[ corrections::score::dun08 ] ) {
		if ( ! decide_read_from_binary_08() ) return false;
	} else {
		if ( ! decide_read_from_binary_02() ) return false;
	}

	return true;
}

/// @details Put specific objections to reading the 02 library from binary here.
/// This function alone will not answer whether the binary file should be read;
/// rather, it is called by the more generic "decide_read_from_binary" function.
bool
RotamerLibrary::decide_read_from_binary_02() const
{
	/// TO DO: check time stamp of ascii file vs time stamp of binary file
	return true;
}

/// @details Put specific objections to reading the 08 library from binary here.
/// This function alone will not answer whether the binary file should be read;
/// rather, it is called by the more generic "decide_read_from_binary" function.
bool
RotamerLibrary::decide_read_from_binary_08() const
{
	/// TO DO: check time stamp of ascii files vs time stamp of binary file
	return true;
}

/// @details Put specific objections to reading the 10 library from binary here.
/// This function alone will not answer whether the binary file should be read;
/// rather, it is called by the more generic "decide_read_from_binary" function.
bool
RotamerLibrary::decide_read_from_binary_10() const
{
	/// TO DO: check time stamp of ascii files vs time stamp of binary file
	return true;
}

/// @details Dispatches binary file is up-to-date to appropriate library given
/// flag status (dun10 or dun08 yes/no).
bool
RotamerLibrary::binary_is_up_to_date( utility::io::izstream & binlib ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ corrections::score::dun10 ] ) {
		return binary_is_up_to_date_10( binlib );
	} else if ( option[ corrections::score::dun08 ] ) {
		return binary_is_up_to_date_08( binlib );
	} else {
		return binary_is_up_to_date_02( binlib );
	}
}

bool
RotamerLibrary::binary_is_up_to_date_02( utility::io::izstream & binlib ) const
{
	boost::int32_t version( 0 );
	binlib.read( (char*) & version, sizeof( boost::int32_t ));
	return static_cast< Size > (version) == current_binary_format_version_id_02();
}

/// @details Binary file is up-to-date if all the following are true:
/// 1. Version # for binary file is equal to the current version #
/// 2. All the hard coded constants for the rotameric and semirotameric libraries
///    are the same from the previous hard coded values.
///
/// DANGER DANGER DANGER
/// When we reopen the binary file in create_fa_dunbrack_libraries_08_from_binary
/// we will want to skip past all this data; the code in that section assumes that
/// it will read 8 arrays: 2 representing the rotameric amino acids, and 6
/// representing the semirotameric amino acids.  The arithmetic it performs depends
/// on this basic structure not changing, therefore, if you change this structure,
/// you must also change the seekg arithmetic in the create...08_from_binary function.
bool
RotamerLibrary::binary_is_up_to_date_08( utility::io::izstream & binlib ) const
{
	/// 1.
	boost::int32_t version( 0 );
	binlib.read( (char*) & version, sizeof( boost::int32_t ));
	if ( static_cast< Size > (version) != current_binary_format_version_id_08() ) return false;

	boost::int32_t nrotameric( 0 );
	binlib.read( (char*) & nrotameric, sizeof( boost::int32_t ) );
	boost::int32_t nsemirotameric( 0 );
	binlib.read( (char*) & nsemirotameric, sizeof( boost::int32_t ) );

	utility::vector1< chemical::AA > rotameric_amino_acids;
	utility::vector1< Size > rotameric_n_chi;

	utility::vector1< chemical::AA > sraa;
	utility::vector1< Size > srnchi;
	utility::vector1< bool > scind;
	utility::vector1< bool > sampind;
	utility::vector1< bool > sym;
	utility::vector1< Real > astr;

	initialize_dun08_aa_parameters(
		rotameric_amino_acids, rotameric_n_chi,
		sraa, srnchi, scind, sampind, sym, astr );

	/// 2.
	if ( static_cast< Size > (nrotameric) != rotameric_amino_acids.size() ) return false;
	if ( static_cast< Size > (nsemirotameric) != sraa.size() ) return false;

	boost::int32_t * rotaa_bin   = new boost::int32_t[ nrotameric ];
	boost::int32_t * rot_nchi_bin= new boost::int32_t[ nrotameric ];

	boost::int32_t * sraa_bin    = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * srnchi_bin  = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * scind_bin   = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * sampind_bin = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * sym_bin     = new boost::int32_t[ nsemirotameric ];
	Real    * astr_bin   = new Real[    nsemirotameric ];

	binlib.read( (char*) rotaa_bin, nrotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) rot_nchi_bin, nrotameric * sizeof( boost::int32_t  ));

	binlib.read( (char*) sraa_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) srnchi_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) scind_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) sampind_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) sym_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) astr_bin, nsemirotameric * sizeof( Real  ));

	bool good( true );
	for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
		if ( rotaa_bin[ ii - 1 ] != static_cast< boost::int32_t > ( rotameric_amino_acids[ ii ] ) ) good = false;
		if ( rot_nchi_bin[ ii - 1 ] != static_cast< boost::int32_t > ( rotameric_n_chi[ ii ] ) ) good = false;
	}
	for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
		if ( sraa_bin[ ii - 1 ] != static_cast< boost::int32_t > ( sraa[ ii ] ) ) good = false;
		if ( srnchi_bin[ ii - 1 ] != static_cast< boost::int32_t > ( srnchi[ ii ] ) ) good = false;
		if ( scind_bin[ ii - 1 ] != static_cast< boost::int32_t > ( scind[ ii ] ) ) good = false;
		if ( sampind_bin[ ii - 1 ] != static_cast< boost::int32_t > ( sampind[ ii ] ) ) good = false;
		if ( sym_bin[ ii - 1 ] != static_cast< boost::int32_t > ( sym[ ii ] ) ) good = false;
		if ( astr_bin[ ii - 1 ] !=  astr[ ii ] ) good = false;
	}

	delete [] rotaa_bin;
	delete [] rot_nchi_bin;
	delete [] sraa_bin;
	delete [] srnchi_bin;
	delete [] scind_bin;
	delete [] sampind_bin;
	delete [] sym_bin;
	delete [] astr_bin;

	return good;
}

bool
RotamerLibrary::binary_is_up_to_date_10( utility::io::izstream & binlib ) const
{
	/// 1.
	boost::int32_t version( 0 );
	binlib.read( (char*) & version, sizeof( boost::int32_t ));
	if ( static_cast< Size > (version) != current_binary_format_version_id_10() ) return false;

	boost::int32_t nrotameric( 0 );
	binlib.read( (char*) & nrotameric, sizeof( boost::int32_t ) );
	boost::int32_t nsemirotameric( 0 );
	binlib.read( (char*) & nsemirotameric, sizeof( boost::int32_t ) );

	utility::vector1< chemical::AA > rotameric_amino_acids;
	utility::vector1< Size > rotameric_n_chi;

	utility::vector1< chemical::AA > sraa;
	utility::vector1< Size > srnchi;
	utility::vector1< bool > scind;
	utility::vector1< bool > sampind;
	utility::vector1< bool > sym;
	utility::vector1< Real > astr;

	initialize_dun08_aa_parameters(
		rotameric_amino_acids, rotameric_n_chi,
		sraa, srnchi, scind, sampind, sym, astr );

	/// 2.
	if ( static_cast< Size > (nrotameric) != rotameric_amino_acids.size() ) return false;
	if ( static_cast< Size > (nsemirotameric) != sraa.size() ) return false;

	boost::int32_t * rotaa_bin   = new boost::int32_t[ nrotameric ];
	boost::int32_t * rot_nchi_bin= new boost::int32_t[ nrotameric ];

	boost::int32_t * sraa_bin    = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * srnchi_bin  = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * scind_bin   = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * sampind_bin = new boost::int32_t[ nsemirotameric ];
	boost::int32_t * sym_bin     = new boost::int32_t[ nsemirotameric ];
	Real    * astr_bin   = new Real[    nsemirotameric ];

	binlib.read( (char*) rotaa_bin, nrotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) rot_nchi_bin, nrotameric * sizeof( boost::int32_t  ));

	binlib.read( (char*) sraa_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) srnchi_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) scind_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) sampind_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) sym_bin, nsemirotameric * sizeof( boost::int32_t  ));
	binlib.read( (char*) astr_bin, nsemirotameric * sizeof( Real  ));

	bool good( true );
	for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
		if ( rotaa_bin[ ii - 1 ] != static_cast< boost::int32_t > ( rotameric_amino_acids[ ii ] ) ) good = false;
		if ( rot_nchi_bin[ ii - 1 ] != static_cast< boost::int32_t > ( rotameric_n_chi[ ii ] ) ) good = false;
	}
	for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
		if ( sraa_bin[ ii - 1 ] != static_cast< boost::int32_t > ( sraa[ ii ] ) ) good = false;
		if ( srnchi_bin[ ii - 1 ] != static_cast< boost::int32_t > ( srnchi[ ii ] ) ) good = false;
		if ( scind_bin[ ii - 1 ] != static_cast< boost::int32_t > ( scind[ ii ] ) ) good = false;
		if ( sampind_bin[ ii - 1 ] != static_cast< boost::int32_t > ( sampind[ ii ] ) ) good = false;
		if ( sym_bin[ ii - 1 ] != static_cast< boost::int32_t > ( sym[ ii ] ) ) good = false;
		if ( astr_bin[ ii - 1 ] !=  astr[ ii ] ) good = false;
	}

	delete [] rotaa_bin;
	delete [] rot_nchi_bin;
	delete [] sraa_bin;
	delete [] srnchi_bin;
	delete [] scind_bin;
	delete [] sampind_bin;
	delete [] sym_bin;
	delete [] astr_bin;

	return good;
}


/// @details First checks general conditions under which a binary library file should
/// not be writen. It relies on get_binary_name() and binary_is_up_to_date() to
/// dispatch to the 02 and 08 data.  If the binary file doesn't exist or isn't up-to-date
/// it then checks if there are any specific 02 or 08 objections to writing the
/// binary file.  If no one objects, then this function gives a green light.
bool
RotamerLibrary::decide_write_binary() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ in::file::no_binary_dunlib ] ) return false;
	if ( option[ out::file::dont_rewrite_dunbrack_database ] ) return false;

	/// place other objections here, e.g. if we're running from a (genuine) database or on BOINC

#ifdef BOINC
	return false;
#endif

#ifdef WIN32 // binary file handling doesnt appear to work properly on windows 64 bit machines.
	return false;
#endif


	std::string binary_filename = get_binary_name();
	utility::io::izstream binlib( binary_filename.c_str(), std::ios::in | std::ios::binary );


	if ( binlib && binary_is_up_to_date( binlib ) ) {
		/// Do not overwrite a binary if it already exists and is up-to-date.
		return false;
	}
	/// if we get this far, then the binary file either doesn't exist or isn't up-to-date.

	/// check if there are any objections to writing the library
	if ( option[ corrections::score::dun10 ] ) {
		if ( ! decide_write_binary_10() ) return false;
	} else if ( option[ corrections::score::dun08 ] ) {
		if ( ! decide_write_binary_08() ) return false;
	} else {
		if ( ! decide_write_binary_02() ) return false;
	}

	// otherwise, write it.
	return true;

}

/// @details Add specific objections to the writing of the 02 library here.
bool
RotamerLibrary::decide_write_binary_02() const
{
	return true;
}

/// @details Add specific objections to the writing of the 08 library here.
bool
RotamerLibrary::decide_write_binary_08() const
{
	return true;
}

/// @details Add specific objections to the writing of the 10 library here.
bool
RotamerLibrary::decide_write_binary_10() const
{
	return true;
}


std::string
RotamerLibrary::get_library_name_02() const
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
	return basic::database::full_name( option[ corrections::score::dun02_file ]);
}

std::string
RotamerLibrary::get_binary_name() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ corrections::score::dun10 ] ) {
		return get_binary_name_10();
	} else if ( option[ corrections::score::dun08 ] ) {
		return get_binary_name_08();
	} else {
		return get_binary_name_02();
	}
}

std::string
RotamerLibrary::get_binary_name_02() const
{
	/// Avoid name collision with old dunbrack binary.
	///return get_library_name_02() + ".bin";
	return basic::database::full_name( "Dunbrack02.lib.bin" );
}

std::string
RotamerLibrary::get_binary_name_08() const
{
	return basic::database::full_name( "Dunbrack08.lib.bin" );
}

std::string
RotamerLibrary::get_binary_name_10() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	return basic::database::full_name( option[ corrections::score::dun10_dir ] + "/" + "Dunbrack10.lib.bin" );
}


/// @details The older binary formats did not have a version number, instead,
/// a somewhat stupid time-stamp was hard coded and compared against the time
/// stamp of the binary file that existed... e.g. if a new binary format was
/// committed to trunk on 1/1/08 12:00 pm, and the time stamp on a particular binary
/// was from 12/1/07, then the old binary file would be overwritten.  However,
/// if the time stamp on some other binary file was from 2/1/08, the code assumed
/// that the binary file was up to date.  This is a poor assumption, since the
/// code from 12/1/07 might have been run on 2/1/08 to generate a binary file.
/// The time stamp would have falsely indicated that the 2/1/08 generated file
/// was of the same binary format imposed on 1/1/08.
///
/// Version numbering avoids this issue.  If the 12/1/07 code wrote version # 3,
/// the 1/1/08 code could look for version # 4 in a file generated on 2/1/08
/// and use the mismatch as an indication of an out-of-date binary format.
///
/// The challenge is to come up with a versioning scheme that handles the fact that
/// until now, there haven't been version #'s in the binary files.  How can a
/// binary file generated by a pre-versioning piece of code be detected?  The trick
/// is to look at the data in the older binary file at the position which will in
/// the future be used to hold a version number.  The first piece of data in the
/// old binary files was a count of the number of libraries in the binary file:
/// 18.  Therefore, any starting version number beyond 18 would distinguish the
/// pre-versioned files from the versioned files.
///
/// When changing the binary format, document the date of the change in the comment
/// section below and increment the version number in the function.  Do not delete
/// old comments from the comment block below.
///
/// Version 19: start vesion number.  8/9/08. Andrew Leaver-Fay
/// Version 20: Limit zero-probability rotamers to the resolution of the library (1e-4).
///             6/16/2009.  Andrew Leaver-Fay
Size
RotamerLibrary::current_binary_format_version_id_02() const
{
	return 20;
}

/// @details Version number for binary format.  See comments for 02 version.
///
/// When changing the binary format, document the date of the change in the comment
/// section below and increment the version number in the function.  Do not delete
/// old comments from the comment block below.
///
/// Version 1: start vesion number.  8/9/08. Andrew Leaver-Fay
/// Version 2: Limit zero-probability rotamers to the resolution of the library (1e-4).
///            Also replacing the Jun8 library with the cross-validated Feb24 '09 library.
///            6/16/2009.  Andrew Leaver-Fay
Size
RotamerLibrary::current_binary_format_version_id_08() const
{
	return 2;
}

/// @details Version number for binary format.  See comments for 02 version.
///
/// When changing the binary format, document the date of the change in the comment
/// section below and increment the version number in the function.  Do not delete
/// old comments from the comment block below.
///
/// Version 1: start vesion number.  6/3/11. Andrew Leaver-Fay
Size
RotamerLibrary::current_binary_format_version_id_10() const
{
	return 1;
}

void
RotamerLibrary::create_fa_dunbrack_libraries_from_ASCII()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ corrections::score::dun10 ] ) {
		return create_fa_dunbrack_libraries_10_from_ASCII();
	} else if ( option[ corrections::score::dun08 ] ) {
		return create_fa_dunbrack_libraries_08_from_ASCII();
	} else {
		return create_fa_dunbrack_libraries_02_from_ASCII();
	}
}

void
RotamerLibrary::create_fa_dunbrack_libraries_from_binary()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ corrections::score::dun10 ] ) {
		return create_fa_dunbrack_libraries_10_from_binary();
	} else if ( option[ corrections::score::dun08 ] ) {
		return create_fa_dunbrack_libraries_08_from_binary();
	} else {
		return create_fa_dunbrack_libraries_02_from_binary();
	}
}

void
RotamerLibrary::create_fa_dunbrack_libraries_02_from_ASCII()
{
	using namespace chemical;

	////////////////////////////////////////////
	//// DATA FOR THE ROTAMERIC AMINO ACIDS
	////////////////////////////////////////////
	utility::vector1< Size > nchi_for_aa( num_canonical_aas, 0 );
	Size nrot_aas;

	utility::vector1< Size > memory_use_rotameric( num_canonical_aas, 0 );

	initialize_dun02_aa_parameters( nchi_for_aa, nrot_aas );


	/// Now read in the 02 library.
	clock_t starttime = clock();
	utility::io::izstream libstream( get_library_name_02() );

	chemical::AA aan = chemical::aa_unk;
	std::string nextaa;
	libstream >> nextaa;

	//std::cout << "first aa: " << nextaa << std::endl;

	SingleResidueDunbrackLibraryOP rot_lib(0);
	Size count_libraries_read( 0 );
	while ( nextaa != "" ) {

		aan = chemical::aa_from_name( nextaa );
		SingleResidueDunbrackLibraryOP newlib = create_rotameric_dunlib(
			aan, nchi_for_aa[ aan ], libstream, true /*dun02*/, nextaa, true /*first aa string already read*/ );
		++count_libraries_read;

		libraries_[ aan ] = newlib();
		aa_libraries_[ aan ] = newlib();
		libraries_ops_.push_back( newlib );
		memory_use_rotameric[ aan ] = newlib->memory_usage_in_bytes();

		//std::cout << "next aa: " << nextaa << " #" << count_libraries_read + 1 << std::endl;

	}
	libstream.close();

	if ( count_libraries_read != nrot_aas ) {
		std::cerr << "ERROR: expected to read " << nrot_aas << " libraries from Dun02, but read " << count_libraries_read << std::endl;
		utility_exit();
	}

	clock_t stoptime = clock();
	TR << "Dunbrack library took " << ((double)stoptime-starttime)/CLOCKS_PER_SEC << " seconds to load from ASCII" << std::endl;

	//TR << "Memory usage: " << std::endl;
	//Size total( 0 );
	//for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {
	//	total += memory_use_rotameric[ ii ];
	//if ( memory_use_rotameric[ii] != 0 ) {
	//TR << chemical::name_from_aa( AA( ii ) ) << " with " << memory_use_rotameric[ ii ] << " bytes" << std::endl;
	//}
	//}
	//TR << "Total memory on Dunbrack Libraries: " << total << " bytes." << std::endl;

}


void
RotamerLibrary::create_fa_dunbrack_libraries_02_from_binary()
{
	clock_t starttime = clock();

	std::string binary_filename = get_binary_name_02();
	utility::io::izstream binlib( binary_filename.c_str(), std::ios::in | std::ios::binary );

	if ( ! binlib ) {
		utility_exit_with_message( "Could not open binary Dunbrack02 file from database -- how did this happen?");
	}

	/// READ PREABMLE
	/// Even if binary file should change its structure in the future, version number should always be
	/// the first piece of data in the file.
	boost::int32_t version( 0 );
	binlib.read( (char*) & version, sizeof( boost::int32_t ));
	/// END PREABMLE


	assert( static_cast< Size > (version) == current_binary_format_version_id_02() );

	read_from_binary( binlib );

	clock_t stoptime = clock();
	TR << "Dunbrack library took " << ((double)stoptime-starttime)/CLOCKS_PER_SEC << " seconds to load from binary" << std::endl;

}

void
RotamerLibrary::create_fa_dunbrack_libraries_08_from_ASCII()
{
	using namespace chemical;
	using namespace utility;

	TR << "Reading Dunbrack Libraries" << std::endl;

	clock_t starttime = clock();

	/// Hard code some stuff.
	/// Rotameric library file names are of the form
	/// <3lc>.Jun9.bbdep.regular.lib
	/// Semi-Rotameric library files names are of the form
	/// 1. <3lc>.Jun9.bbdep.regular.lib
	/// 2. <3lc>.Jun9.bbdep.contMin.lib
	/// 3. <3lc>.Jun9.bbind.chi<#>.Definitions.lib
	/// 4. <3lc>.Jun9.bbind.chi<#>.Probabities.lib

	//std::string const db_subdir = "dun08/";
	std::string const regular_lib_suffix = ".Feb24.Optim.bbdep.regular.lib";
	std::string const bbdep_contmin = ".Feb24.Optim.bbdep.contMinim.lib";
	std::string const bbind_midfix  = ".Feb24.Optim.bbind.chi";
	std::string const bbind_defs    = ".Definitions.lib";
	//std::string const bbind_probs   = ".Probabilities.lib";

	////////////////////////////////////////////
	//// DATA FOR THE ROTAMERIC AMINO ACIDS
	////////////////////////////////////////////
	utility::vector1< chemical::AA > rotameric_amino_acids;
	utility::vector1< Size > rotameric_n_chi;
	utility::vector1< Size > memory_use_rotameric;

	////////////////////////////////////////////
	//// DATA FOR THE SEMI-ROTAMERIC AMINO ACIDS
	////////////////////////////////////////////

	/// AA code for entry i
	utility::vector1< chemical::AA > sraa;
	/// Number of rotameric chi for entry i; the index of the non-rotameric chi is +1 of the value stored.
	utility::vector1< Size > srnchi;
	/// Decision: SCore the nonrotameric chi in a backbone INDependent manner?
	utility::vector1< bool > scind;
	/// Decision: SAMPle the nonrotameric chi in a backbone INDependent manner?
	utility::vector1< bool > sampind;
	/// Is there symmetry about the rotameric chi?  If so, periodicity is 180 degrees, else, 360.
	/// Furthermore, if it's symmetric, then the bbdep sampling is at 5 degree intervals and
	/// the bbind is at 0.5 degree intervals; otherwise, bbdep is at 10 degrees, and bbind is at 1.
	utility::vector1< bool > sym;
	/// At what starting angle does the chi data come from?
	utility::vector1< Real > astr;

	utility::vector1< Size > memory_use_semirotameric;

	/// Initialize the data.

	initialize_dun08_aa_parameters(
		rotameric_amino_acids, rotameric_n_chi,
		sraa, srnchi, scind, sampind, sym, astr
	);

	for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
		std::string next_aa_in_library;

		chemical::AA ii_aa( rotameric_amino_acids[ ii ] );
		std::string  ii_lc_3lc( chemical::name_from_aa( ii_aa ) );
		std::transform(ii_lc_3lc.begin(), ii_lc_3lc.end(), ii_lc_3lc.begin(), tolower);
		std::string library_name = basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + regular_lib_suffix ;
		utility::io::izstream lib( library_name.c_str() );
		if ( !lib.good() ) {
			lib.close();
			library_name = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + regular_lib_suffix );
			lib.open( library_name );
		}
		if ( !lib.good() ) {
			utility_exit_with_message( "Unable to open database file for Dun10 rotamer library: " + library_name );
		}
		TR << "Reading " << library_name << std::endl;

		SingleResidueDunbrackLibraryOP newlib = create_rotameric_dunlib(
			ii_aa, rotameric_n_chi[ ii ], lib, false, next_aa_in_library, false );
		if ( next_aa_in_library != "" ) {
			std::cerr << "ERROR: Inappropriate read from dun08 input while reading " << ii_aa
				<< ": \"" << next_aa_in_library << "\"" << std::endl;
			utility_exit();
		}

		memory_use_rotameric.push_back( newlib->memory_usage_in_bytes() );

		libraries_[ ii_aa ] = newlib();
		aa_libraries_[ ii_aa ] = newlib();
		libraries_ops_.push_back( newlib );
	}

	for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
		chemical::AA ii_aa( sraa[ ii ] );
		std::string  ii_lc_3lc( chemical::name_from_aa( ii_aa ) );
		std::transform(ii_lc_3lc.begin(), ii_lc_3lc.end(), ii_lc_3lc.begin(), tolower);
		Size const nrchi = srnchi[ ii ] + 1;

		std::string reg_lib_name = basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + regular_lib_suffix ;
		std::string conmin_name  = basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + bbdep_contmin ;
		std::string rotdef_name  = basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_defs ;
		//std::string probs_name   = basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_probs ;

		utility::io::izstream lib(     reg_lib_name.c_str() );
		utility::io::izstream contmin( conmin_name.c_str()  );
		utility::io::izstream defs(    rotdef_name.c_str()  );
		//utility::io::izstream probs(   probs_name.c_str()   );

		if ( !lib.good() ) {
			defs.close();
			lib.close();
			contmin.close();
			//probs.close();

			reg_lib_name = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + regular_lib_suffix );
			conmin_name  = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + bbdep_contmin );
			rotdef_name  = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_defs );
			//probs_name   = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_probs );

			lib.open(reg_lib_name);
			contmin.open(conmin_name);
			defs.open(rotdef_name);
			//probs.open(probs_name);
		}
		if ( !lib.good()  ) {
			utility_exit_with_message( "Unable to open database file for Dun10 rotamer library: " + reg_lib_name );
		}
		if ( !contmin.good()  ) {
			utility_exit_with_message( "Unable to open database file for Dun10 rotamer library: " + conmin_name );
		}
		if ( !defs.good() ) {
			utility_exit_with_message( "Unable to open database file for Dun10 rotamer library: " + rotdef_name );
		}

		TR << "Reading " << reg_lib_name << std::endl;
		TR << "Reading " << conmin_name << std::endl;
		TR << "Reading " << rotdef_name << std::endl;
		//TR << "Reading " << probs_name << std::endl;

		SingleResidueDunbrackLibraryOP newlib = create_semi_rotameric_dunlib(
			ii_aa, srnchi[ ii ],
			scind[ ii ], sampind[ ii ],
			sym[ ii ], astr[ ii ],
			defs, lib, contmin );

		memory_use_semirotameric.push_back( newlib->memory_usage_in_bytes() );

		libraries_[ ii_aa ] = newlib();
		aa_libraries_[ ii_aa ] = newlib();
		libraries_ops_.push_back( newlib );
	}
	TR << "Finished reading Dunbrack Libraries" << std::endl;

	clock_t stoptime = clock();
	TR << "Dunbrack 2008 library took " << ((double)stoptime-starttime)/CLOCKS_PER_SEC << " seconds to load from ASCII" << std::endl;


	TR << "Memory usage: " << std::endl;
	Size total( 0 );
	for ( Size ii = 1; ii <= memory_use_rotameric.size(); ++ii ) {
		total += memory_use_rotameric[ ii ];
		TR << chemical::name_from_aa(rotameric_amino_acids[ ii ]) << " with " << memory_use_rotameric[ ii ] << " bytes" << std::endl;
	}

	for ( Size ii = 1; ii <= memory_use_semirotameric.size(); ++ii ) {
		total += memory_use_semirotameric[ ii ];
		TR << chemical::name_from_aa(sraa[ ii ]) << " with " << memory_use_semirotameric[ ii ] << " bytes" << std::endl;
	}
	TR << "Total memory on Dunbrack Libraries: " << total << " bytes." << std::endl;

}


void
RotamerLibrary::create_fa_dunbrack_libraries_08_from_binary()
{
	clock_t starttime = clock();

	std::string binary_filename = get_binary_name_08();
	utility::io::izstream binlib( binary_filename.c_str(), std::ios::in | std::ios::binary );

	if ( ! binlib ) {
		utility_exit_with_message( "Could not open binary Dunbrack08 file from database -- how did this happen?");
	}

	/// READ PREABMLE
	/// Even if binary file should change its structure in the future, version number should always be
	/// the first piece of data in the file.
	boost::int32_t version( 0 );
	binlib.read( (char*) & version, sizeof( boost::int32_t ));

	/// Version 1: write out hard coded parameters for 08 library.  Verify when reading in
	/// the binary file that the parameters have not changed.
	/// DANGER DANGER DANGER
	/// If you change the binary format and adjust the preamble to the binary file,
	/// change the arithmetic in this function so that the correct number of bytes
	/// are jumped past.

	boost::int32_t n_rotameric_aas( 0 );
	binlib.read( (char*) & n_rotameric_aas, sizeof( boost::int32_t ) );

	boost::int32_t n_semirotameric_aas( 0 );
	binlib.read( (char*) & n_semirotameric_aas, sizeof( boost::int32_t ) );

	Size const n_rotameric_arrays_read_of_int32s = 2;
	Size const n_semirot_arrays_read_of_int32s = 5;
	Size const n_semirot_arrays_read_of_Reals = 1;
	Size const bytes_to_skip =
		( n_rotameric_arrays_read_of_int32s * n_rotameric_aas +
		n_semirot_arrays_read_of_int32s * n_semirotameric_aas ) * sizeof( boost::int32_t ) +
		n_semirot_arrays_read_of_Reals * n_semirotameric_aas * sizeof( Real );

	//binlib.seekg( bytes_to_skip, std::ios::cur ); //seekg not available for izstream...

	char * temp_buffer = new char[ bytes_to_skip ];
	binlib.read( temp_buffer, bytes_to_skip );
	delete [] temp_buffer;
	/// END PREABMLE


	read_from_binary( binlib );

	clock_t stoptime = clock();
	TR << "Dunbrack 2008 library took " << ((double)stoptime-starttime)/CLOCKS_PER_SEC << " seconds to load from binary" << std::endl;

}

void
RotamerLibrary::create_fa_dunbrack_libraries_10_from_ASCII()
{
	using namespace chemical;
	using namespace utility;

	TR << "Reading Dunbrack Libraries" << std::endl;

	clock_t starttime = clock();

	/// Hard code some stuff.
	/// Rotameric library file names are of the form
	/// <3lc>.Jun9.bbdep.regular.lib
	/// Semi-Rotameric library files names are of the form
	/// 1. <3lc>.Jun9.bbdep.regular.lib
	/// 2. <3lc>.Jun9.bbdep.contMin.lib
	/// 3. <3lc>.Jun9.bbind.chi<#>.Definitions.lib
	/// 4. <3lc>.Jun9.bbind.chi<#>.Probabities.lib

	//std::string const db_subdir = "dun08/";
	std::string const regular_lib_suffix = ".bbdep.rotamers.lib";
	std::string const bbdep_contmin = ".bbdep.densities.lib";
	std::string const bbind_midfix  = ".bbind.chi";
	std::string const bbind_defs    = ".Definitions.lib";
	// Note: No Backbone-Independent NRCHI representation in '10 library
	// std::string const bbind_probs   = ".Probabilities.lib";

	////////////////////////////////////////////
	//// DATA FOR THE ROTAMERIC AMINO ACIDS
	////////////////////////////////////////////
	utility::vector1< chemical::AA > rotameric_amino_acids;
	utility::vector1< Size > rotameric_n_chi;
	utility::vector1< Size > memory_use_rotameric;

	////////////////////////////////////////////
	//// DATA FOR THE SEMI-ROTAMERIC AMINO ACIDS
	////////////////////////////////////////////

	/// AA code for entry i
	utility::vector1< chemical::AA > sraa;
	/// Number of rotameric chi for entry i; the index of the non-rotameric chi is +1 of the value stored.
	utility::vector1< Size > srnchi;
	/// Decision: SCore the nonrotameric chi in a backbone INDependent manner?
	utility::vector1< bool > scind;
	/// Decision: SAMPle the nonrotameric chi in a backbone INDependent manner?
	utility::vector1< bool > sampind;
	/// Is there symmetry about the rotameric chi?  If so, periodicity is 180 degrees, else, 360.
	/// Furthermore, if it's symmetric, then the bbdep sampling is at 5 degree intervals and
	/// the bbind is at 0.5 degree intervals; otherwise, bbdep is at 10 degrees, and bbind is at 1.
	utility::vector1< bool > sym;
	/// At what starting angle does the chi data come from?
	utility::vector1< Real > astr;

	utility::vector1< Size > memory_use_semirotameric;

	/// Initialize the data.

	initialize_dun08_aa_parameters(
		rotameric_amino_acids, rotameric_n_chi,
		sraa, srnchi, scind, sampind, sym, astr
	); // use the same parameters as dun08

	for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
		std::string next_aa_in_library;

		chemical::AA ii_aa( rotameric_amino_acids[ ii ] );
		std::string  ii_lc_3lc( chemical::name_from_aa( ii_aa ) );
		std::transform(ii_lc_3lc.begin(), ii_lc_3lc.end(), ii_lc_3lc.begin(), tolower);
		std::string library_name = basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + regular_lib_suffix ;
		utility::io::izstream lib( library_name.c_str() );
		if ( !lib.good() ) {
			lib.close();
			library_name = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + regular_lib_suffix );
			lib.open( library_name );
		}
		if ( !lib.good() ) {
			utility_exit_with_message( "Unable to open database file for Dun10 rotamer library: " + library_name );
		}
		TR << "Reading " << library_name << std::endl;

		SingleResidueDunbrackLibraryOP newlib = create_rotameric_dunlib(
			ii_aa, rotameric_n_chi[ ii ], lib, false, next_aa_in_library, false );
		if ( next_aa_in_library != "" ) {
			std::cerr << "ERROR: Inappropriate read from dun10 input while reading " << ii_aa
				<< ": \"" << next_aa_in_library << "\"" << std::endl;
			utility_exit();
		}

		memory_use_rotameric.push_back( newlib->memory_usage_in_bytes() );

		libraries_[ ii_aa ] = newlib();
		aa_libraries_[ ii_aa ] = newlib();
		libraries_ops_.push_back( newlib );
	}

	for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
		chemical::AA ii_aa( sraa[ ii ] );
		std::string  ii_lc_3lc( chemical::name_from_aa( ii_aa ) );
		std::transform(ii_lc_3lc.begin(), ii_lc_3lc.end(), ii_lc_3lc.begin(), tolower);
		Size const nrchi = srnchi[ ii ] + 1;

		std::string reg_lib_name = basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + regular_lib_suffix ;
		std::string conmin_name  = basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + bbdep_contmin ;
		std::string rotdef_name  = basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_defs ;
		//std::string probs_name   = basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_probs ;

		utility::io::izstream lib(     reg_lib_name.c_str() );
		utility::io::izstream contmin( conmin_name.c_str()  );
		utility::io::izstream defs(    rotdef_name.c_str()  );
		//utility::io::izstream probs(   probs_name.c_str()   );

		if ( !lib.good() ) {
			defs.close();
			lib.close();
			contmin.close();
			//probs.close();

			reg_lib_name = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + regular_lib_suffix );
			conmin_name  = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + bbdep_contmin );
			rotdef_name  = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun10_dir]() + "/" + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_defs );
			//probs_name   = basic::database::full_name( basic::options::option[basic::options::OptionKeys::corrections::score::dun08_dir]() + "/" + ii_lc_3lc + bbind_midfix + to_string(nrchi) + bbind_probs );

			lib.open(reg_lib_name);
			contmin.open(conmin_name);
			defs.open(rotdef_name);
			//probs.open(probs_name);
		}

		if ( !lib.good()  ) {
			utility_exit_with_message( "Unable to open database file for Dun10 rotamer library: " + reg_lib_name );
		}
		if ( !contmin.good()  ) {
			utility_exit_with_message( "Unable to open database file for Dun10 rotamer library: " + conmin_name );
		}
		if ( !defs.good() ) {
			utility_exit_with_message( "Unable to open database file for Dun10 rotamer library: " + rotdef_name );
		}

		TR << "Reading " << reg_lib_name << std::endl;
		TR << "Reading " << conmin_name << std::endl;
		TR << "Reading " << rotdef_name << std::endl;
		//TR << "Reading " << probs_name << std::endl;

		SingleResidueDunbrackLibraryOP newlib = create_semi_rotameric_dunlib(
			ii_aa, srnchi[ ii ],
			scind[ ii ], sampind[ ii ],
			sym[ ii ], astr[ ii ],
			defs, lib, contmin );

		memory_use_semirotameric.push_back( newlib->memory_usage_in_bytes() );

		libraries_[ ii_aa ] = newlib();
		aa_libraries_[ ii_aa ] = newlib();
		libraries_ops_.push_back( newlib );
	}
	TR << "Finished reading Dunbrack Libraries" << std::endl;

	clock_t stoptime = clock();
	TR << "Dunbrack 2010 library took " << ((double)stoptime-starttime)/CLOCKS_PER_SEC << " seconds to load from ASCII" << std::endl;


	TR << "Memory usage: " << std::endl;
	Size total( 0 );
	for ( Size ii = 1; ii <= memory_use_rotameric.size(); ++ii ) {
		total += memory_use_rotameric[ ii ];
		TR << chemical::name_from_aa(rotameric_amino_acids[ ii ]) << " with " << memory_use_rotameric[ ii ] << " bytes" << std::endl;
	}

	for ( Size ii = 1; ii <= memory_use_semirotameric.size(); ++ii ) {
		total += memory_use_semirotameric[ ii ];
		TR << chemical::name_from_aa(sraa[ ii ]) << " with " << memory_use_semirotameric[ ii ] << " bytes" << std::endl;
	}
	TR << "Total memory on Dunbrack Libraries: " << total << " bytes." << std::endl;

}

void
RotamerLibrary::create_fa_dunbrack_libraries_10_from_binary()
{
	clock_t starttime = clock();

	std::string binary_filename = get_binary_name_10();
	utility::io::izstream binlib( binary_filename.c_str(), std::ios::in | std::ios::binary );

	if ( ! binlib ) {
		utility_exit_with_message( "Could not open binary Dunbrack10 file from database -- how did this happen?");
	}

	/// READ PREABMLE
	/// Even if binary file should change its structure in the future, version number should always be
	/// the first piece of data in the file.
	boost::int32_t version( 0 );
	binlib.read( (char*) & version, sizeof( boost::int32_t ));

	/// Version 1: write out hard coded parameters for 08 library.  Verify when reading in
	/// the binary file that the parameters have not changed.
	/// DANGER DANGER DANGER
	/// If you change the binary format and adjust the preamble to the binary file,
	/// change the arithmetic in this function so that the correct number of bytes
	/// are jumped past.

	boost::int32_t n_rotameric_aas( 0 );
	binlib.read( (char*) & n_rotameric_aas, sizeof( boost::int32_t ) );

	boost::int32_t n_semirotameric_aas( 0 );
	binlib.read( (char*) & n_semirotameric_aas, sizeof( boost::int32_t ) );

	Size const n_rotameric_arrays_read_of_int32s = 2;
	Size const n_semirot_arrays_read_of_int32s = 5;
	Size const n_semirot_arrays_read_of_Reals = 1;
	Size const bytes_to_skip =
		( n_rotameric_arrays_read_of_int32s * n_rotameric_aas +
		n_semirot_arrays_read_of_int32s * n_semirotameric_aas ) * sizeof( boost::int32_t ) +
		n_semirot_arrays_read_of_Reals * n_semirotameric_aas * sizeof( Real );

	//binlib.seekg( bytes_to_skip, std::ios::cur ); //seekg not available for izstream...

	char * temp_buffer = new char[ bytes_to_skip ];
	binlib.read( temp_buffer, bytes_to_skip );
	delete [] temp_buffer;
	/// END PREABMLE


	read_from_binary( binlib );

	clock_t stoptime = clock();
	TR << "Dunbrack 2010 library took " << ((double)stoptime-starttime)/CLOCKS_PER_SEC << " seconds to load from binary" << std::endl;
}

void
RotamerLibrary::write_binary_fa_dunbrack_libraries() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ corrections::score::dun10 ] ) {
		return write_binary_fa_dunbrack_libraries_10();
	} else if ( option[ corrections::score::dun08 ] ) {
		return write_binary_fa_dunbrack_libraries_08();
	} else {
		return write_binary_fa_dunbrack_libraries_02();
	}
}

void
RotamerLibrary::write_binary_fa_dunbrack_libraries_02() const
{
	std::string binary_filename = get_binary_name_02();
	std::string tempfilename = random_tempname( "dun02_binary" );

	TR << "Opening file " << tempfilename << " for output." << std::endl;

	utility::io::ozstream binlib( tempfilename.c_str(), std::ios::out | std::ios::binary );
	if ( binlib ) {

		/// WRITE PREABMLE
		/// Even if binary file should change its structure in the future, version number should always be
		/// the first piece of data in the file.
		boost::int32_t version( current_binary_format_version_id_02() );
		binlib.write( (char*) & version, sizeof( boost::int32_t ));
		/// END PREABMLE

		write_to_binary( binlib );
		binlib.close();

		// Move the temporary file to its permanent location
		TR << "Moving temporary file to " << binary_filename.c_str() << std::endl;
		rename( tempfilename.c_str(), binary_filename.c_str() );
	} else {
		TR << "Unable to open temporary file in rosetta database for writing the binary version of the Dunbrack02 library." << std::endl;
	}

}

void
RotamerLibrary::write_binary_fa_dunbrack_libraries_08() const
{
	std::string binary_filename = get_binary_name_08();
	std::string tempfilename = random_tempname( "dun08_binary" );

	TR << "Opening file " << tempfilename << " for output." << std::endl;

	utility::io::ozstream binlib( tempfilename.c_str(), std::ios::out | std::ios::binary );
	if ( binlib ) {

		/// WRITE PREABMLE
		/// Even if binary file should change its structure in the future, version number should always be
		/// the first piece of data in the file.
		boost::int32_t version( current_binary_format_version_id_08() );
		binlib.write( (char*) & version, sizeof( boost::int32_t ));

		utility::vector1< chemical::AA > rotameric_amino_acids;
		utility::vector1< Size > rotameric_n_chi;

		utility::vector1< chemical::AA > sraa;
		utility::vector1< Size > srnchi;
		utility::vector1< bool > scind;
		utility::vector1< bool > sampind;
		utility::vector1< bool > sym;
		utility::vector1< Real > astr;

		initialize_dun08_aa_parameters(
			rotameric_amino_acids, rotameric_n_chi,
			sraa, srnchi, scind, sampind, sym, astr );

		boost::int32_t nrotameric( static_cast< boost::int32_t > (rotameric_amino_acids.size() ));
		boost::int32_t nsemirotameric( static_cast< boost::int32_t > (sraa.size()) );

		binlib.write( (char*) & nrotameric, sizeof( boost::int32_t ) );
		binlib.write( (char*) & nsemirotameric, sizeof( boost::int32_t ) );

		/// 2.
		boost::int32_t * rotaa_bin   = new boost::int32_t[ nrotameric ];
		boost::int32_t * rot_nchi_bin= new boost::int32_t[ nrotameric ];

		boost::int32_t * sraa_bin    = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * srnchi_bin  = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * scind_bin   = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * sampind_bin = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * sym_bin     = new boost::int32_t[ nsemirotameric ];
		Real    * astr_bin    = new Real[    nsemirotameric ];

		for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
			rotaa_bin[ ii - 1 ] = static_cast< boost::int32_t > ( rotameric_amino_acids[ ii ] );
			rot_nchi_bin[ ii - 1 ] = static_cast< boost::int32_t > ( rotameric_n_chi[ ii ] );
		}
		for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
			sraa_bin[ ii - 1 ]    = static_cast< boost::int32_t > ( sraa[ ii ] );
			srnchi_bin[ ii - 1 ]  = static_cast< boost::int32_t > ( srnchi[ ii ] );
			scind_bin[ ii - 1 ]   = static_cast< boost::int32_t > ( scind[ ii ] );
			sampind_bin[ ii - 1 ] = static_cast< boost::int32_t > ( sampind[ ii ] );
			sym_bin[ ii - 1 ]     = static_cast< boost::int32_t > ( sym[ ii ] );
			astr_bin[ ii - 1 ]  =  astr[ ii ];
		}


		binlib.write( (char*) rotaa_bin, nrotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) rot_nchi_bin, nrotameric * sizeof( boost::int32_t  ));

		binlib.write( (char*) sraa_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) srnchi_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) scind_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) sampind_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) sym_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) astr_bin, nsemirotameric * sizeof( Real  ));

		delete [] rotaa_bin;
		delete [] rot_nchi_bin;
		delete [] sraa_bin;
		delete [] srnchi_bin;
		delete [] scind_bin;
		delete [] sampind_bin;
		delete [] sym_bin;
		delete [] astr_bin;

		/// END PREABMLE

		write_to_binary( binlib );
		binlib.close();

		// Move the temporary file to its permanent location
		//std::cerr << "Moving temporary file to " << binary_filename.c_str() << std::endl;
		rename( tempfilename.c_str(), binary_filename.c_str() );
	} else {
		TR << "Unable to open temporary file in rosetta database for writing the binary version of the Dunbrack08 library." << std::endl;
	}
}

void
RotamerLibrary::write_binary_fa_dunbrack_libraries_10() const
{
	std::string binary_filename = get_binary_name_10();
	std::string tempfilename = random_tempname( "dun10_binary" );

	TR << "Opening file " << tempfilename << " for output." << std::endl;

	utility::io::ozstream binlib( tempfilename.c_str(), std::ios::out | std::ios::binary );
	if ( binlib ) {

		/// WRITE PREABMLE
		/// Even if binary file should change its structure in the future, version number should always be
		/// the first piece of data in the file.
		boost::int32_t version( current_binary_format_version_id_10() );
		binlib.write( (char*) & version, sizeof( boost::int32_t ));

		utility::vector1< chemical::AA > rotameric_amino_acids;
		utility::vector1< Size > rotameric_n_chi;

		utility::vector1< chemical::AA > sraa;
		utility::vector1< Size > srnchi;
		utility::vector1< bool > scind;
		utility::vector1< bool > sampind;
		utility::vector1< bool > sym;
		utility::vector1< Real > astr;

		initialize_dun08_aa_parameters(
			rotameric_amino_acids, rotameric_n_chi,
			sraa, srnchi, scind, sampind, sym, astr );

		boost::int32_t nrotameric( static_cast< boost::int32_t > (rotameric_amino_acids.size() ));
		boost::int32_t nsemirotameric( static_cast< boost::int32_t > (sraa.size()) );

		binlib.write( (char*) & nrotameric, sizeof( boost::int32_t ) );
		binlib.write( (char*) & nsemirotameric, sizeof( boost::int32_t ) );

		/// 2.
		boost::int32_t * rotaa_bin   = new boost::int32_t[ nrotameric ];
		boost::int32_t * rot_nchi_bin= new boost::int32_t[ nrotameric ];

		boost::int32_t * sraa_bin    = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * srnchi_bin  = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * scind_bin   = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * sampind_bin = new boost::int32_t[ nsemirotameric ];
		boost::int32_t * sym_bin     = new boost::int32_t[ nsemirotameric ];
		Real    * astr_bin    = new Real[    nsemirotameric ];

		for ( Size ii = 1; ii <= rotameric_amino_acids.size(); ++ii ) {
			rotaa_bin[ ii - 1 ] = static_cast< boost::int32_t > ( rotameric_amino_acids[ ii ] );
			rot_nchi_bin[ ii - 1 ] = static_cast< boost::int32_t > ( rotameric_n_chi[ ii ] );
		}
		for ( Size ii = 1; ii <= sraa.size(); ++ii ) {
			sraa_bin[ ii - 1 ]    = static_cast< boost::int32_t > ( sraa[ ii ] );
			srnchi_bin[ ii - 1 ]  = static_cast< boost::int32_t > ( srnchi[ ii ] );
			scind_bin[ ii - 1 ]   = static_cast< boost::int32_t > ( scind[ ii ] );
			sampind_bin[ ii - 1 ] = static_cast< boost::int32_t > ( sampind[ ii ] );
			sym_bin[ ii - 1 ]     = static_cast< boost::int32_t > ( sym[ ii ] );
			astr_bin[ ii - 1 ]  =  astr[ ii ];
		}


		binlib.write( (char*) rotaa_bin, nrotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) rot_nchi_bin, nrotameric * sizeof( boost::int32_t  ));

		binlib.write( (char*) sraa_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) srnchi_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) scind_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) sampind_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) sym_bin, nsemirotameric * sizeof( boost::int32_t  ));
		binlib.write( (char*) astr_bin, nsemirotameric * sizeof( Real  ));

		delete [] rotaa_bin;
		delete [] rot_nchi_bin;
		delete [] sraa_bin;
		delete [] srnchi_bin;
		delete [] scind_bin;
		delete [] sampind_bin;
		delete [] sym_bin;
		delete [] astr_bin;

		/// END PREABMLE

		write_to_binary( binlib );
		binlib.close();

		// Move the temporary file to its permanent location
		//std::cerr << "Moving temporary file to " << binary_filename.c_str() << std::endl;
		rename( tempfilename.c_str(), binary_filename.c_str() );
	} else {
		TR << "Unable to open temporary file in rosetta database for writing the binary version of the Dunbrack '10 library." << std::endl;
	}
}


std::string
RotamerLibrary::random_tempname( std::string const & prefix ) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// Ugly C stuff...
#ifdef WIN32
	char * tmpname_output = _tempnam( option[ in::path::database ](1).name().c_str(), prefix.c_str() );
#elif __CYGWIN__
	char * tmpname_output = tempnam( option[ in::path::database ](1).name().c_str(), prefix.c_str() );
#elif USEMPI
	// jk change this to mkstemp to prevent race condition
	// apl -- any reason mkstemp can't become the standard?
	std::string str_name = option[ in::path::database ](1).name() + "/" + prefix + "XXXXXX";
	char *tmpname_output = new char [str_name.size()+20]; // Why + 20?  Apparently mkstemp can overflow.
	strcpy(tmpname_output, str_name.c_str()); // use strcpy to get char* instead of const char*
	int fid = mkstemp( tmpname_output );
	if ( fid == -1 ) { // error condition of mkstemp.
		utility_exit_with_message("Failed to open temporary file with mkstemp" );
	}
	close( fid ); // have to close the file in order to open it again.
	/// On the web, it says the file descriptor alone returned by mkstemp and not the file name is hacker safe.
	/// http://lists.virus.org/vulnwatch-0212/msg00023.html
	/// I suppose this is a security flaw with mini... just as soon as someone figures out
	/// how to open a stream with a file descriptor, this code should be replaced.
#else
	char * tmpname_output = tempnam( option[ in::path::database ](1).name().c_str(), prefix.c_str() );
#endif

	/// Copy c-style string to std::string.
	std::string tempfilename = tmpname_output;

#ifdef USEMPI
	delete [] tmpname_output;
#else
	free( tmpname_output );
#endif

	return tempfilename;
}

SingleResidueDunbrackLibraryOP
RotamerLibrary::create_rotameric_dunlib(
	chemical::AA aa,
	Size const n_chi,
	utility::io::izstream & library,
	bool dun02,
	std::string & next_aa_in_library,
	bool first_three_letter_code_already_read
) const
{
	SingleResidueDunbrackLibraryOP rotlib;
	/// scope the case statements
	switch ( n_chi ) {
		case 1: {
			RotamericSingleResidueDunbrackLibrary< ONE > * r1 =
				new RotamericSingleResidueDunbrackLibrary< ONE >( aa, dun02 );
			next_aa_in_library = r1->read_from_file( library, first_three_letter_code_already_read );
			rotlib = r1;
			break;
		}
		case 2: {
			RotamericSingleResidueDunbrackLibrary< TWO > * r2 =
				new RotamericSingleResidueDunbrackLibrary< TWO >( aa, dun02 );
			next_aa_in_library = r2->read_from_file( library, first_three_letter_code_already_read );
			rotlib = r2;
			break;
		}
		case 3: {
			RotamericSingleResidueDunbrackLibrary< THREE > * r3 =
				new RotamericSingleResidueDunbrackLibrary< THREE >( aa, dun02 );
			next_aa_in_library = r3->read_from_file( library, first_three_letter_code_already_read );
			rotlib = r3;
			break;
		}
		case 4: {
			RotamericSingleResidueDunbrackLibrary< FOUR > * r4 =
				new RotamericSingleResidueDunbrackLibrary< FOUR >( aa, dun02 );
			next_aa_in_library = r4->read_from_file( library, first_three_letter_code_already_read );
			rotlib = r4;
			break;
		}
		default:
			utility_exit_with_message( "ERROR: too many chi angles desired for dunbrack library: " + n_chi );
		break;
	}
	return rotlib;

}

SingleResidueDunbrackLibraryOP
RotamerLibrary::create_rotameric_dunlib(
	chemical::AA aa,
	Size const n_chi,
	bool dun02
) const
{
	SingleResidueDunbrackLibraryOP rotlib;
	switch ( n_chi ) {
		case 1:
			rotlib = new RotamericSingleResidueDunbrackLibrary< ONE >( aa, dun02 );
			break;
		case 2:
			rotlib = new RotamericSingleResidueDunbrackLibrary< TWO >( aa, dun02 );
			break;
		case 3:
			rotlib = new RotamericSingleResidueDunbrackLibrary< THREE >( aa, dun02 );
			break;
		case 4:
			rotlib = new RotamericSingleResidueDunbrackLibrary< FOUR >( aa, dun02 );
			break;
		default:
			utility_exit_with_message( "ERROR: too many chi angles desired for dunbrack library: " + n_chi );
		break;
	}
	return rotlib;

}


SingleResidueDunbrackLibraryOP
RotamerLibrary::create_semi_rotameric_dunlib(
	chemical::AA aa,
	Size const nchi,
	bool const use_bbind_rnchi_scoring,
	bool const use_bbind_rnchi_sampling,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle,
	utility::io::izstream & rotamer_definitions,
	utility::io::izstream & regular_library,
	utility::io::izstream & continuous_minimization_bbdep
	// phasing out bbind sampling -- utility::io::izstream & rnchi_bbind_probabilities
) const
{
	SingleResidueDunbrackLibraryOP rotlib;
	/// scope the case statements
	switch ( nchi ) {
		case 1: {
			SemiRotamericSingleResidueDunbrackLibrary< ONE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = r1;
			break;
		}
		case 2: {
			SemiRotamericSingleResidueDunbrackLibrary< TWO > * r2 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_and_read_srsrdl( *r2, nrchi_is_symmetric, nrchi_start_angle,
				rotamer_definitions, regular_library, continuous_minimization_bbdep );
			rotlib = r2;
			break;
		}

		default:
			utility_exit_with_message( "ERROR: too many chi angles desired for semi-rotameric dunbrack library: " + nchi );
		break;
	}
	return rotlib;

}


SingleResidueDunbrackLibraryOP
RotamerLibrary::create_semi_rotameric_dunlib(
	chemical::AA aa,
	Size const nchi,
	bool const use_bbind_rnchi_scoring,
	bool const use_bbind_rnchi_sampling,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle
) const
{
	SingleResidueDunbrackLibraryOP rotlib;
	switch ( nchi ) {
		case 1: {
			SemiRotamericSingleResidueDunbrackLibrary< ONE > * r1 =
				new SemiRotamericSingleResidueDunbrackLibrary< ONE >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r1, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = r1;
			break;
		}
		case 2: {
			SemiRotamericSingleResidueDunbrackLibrary< TWO > * r2 =
				new SemiRotamericSingleResidueDunbrackLibrary< TWO >( aa, use_bbind_rnchi_scoring, use_bbind_rnchi_sampling );
			initialize_srsrdl( *r2, nrchi_is_symmetric, nrchi_start_angle );
			rotlib = r2;
			break;
		}

		default:
			utility_exit_with_message( "ERROR: too many chi angles desired for semi-rotameric dunbrack library: " + nchi );
		break;
	}
	return rotlib;

}

/// @details Build and initialize some of the data in an SRDL in preparation for
/// binary file input -- that is, the SRDL does not have "read_from_file(s)" called
/// before it is returned to the calling function.  ATM not threadsafe.
SingleResidueDunbrackLibraryOP
RotamerLibrary::create_srdl(
	chemical::AA aa_in
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static bool initialized( false );

	static utility::vector1< chemical::AA > rotameric_amino_acids;
	static utility::vector1< Size > rotameric_n_chi;

	static utility::vector1< chemical::AA > sraa;
	static utility::vector1< Size > srnchi;
	static utility::vector1< bool > scind;
	static utility::vector1< bool > sampind;
	static utility::vector1< bool > sym;
	static utility::vector1< Real > astr;

	if ( ! initialized ) {
		/// put a mutex here...
		if ( option[ corrections::score::dun10 ] || option[ corrections::score::dun08 ] ) {
			initialize_dun08_aa_parameters(
				rotameric_amino_acids, rotameric_n_chi,
				sraa, srnchi, scind, sampind, sym, astr );
		} else {
			initialize_dun02_aa_parameters(
				rotameric_amino_acids, rotameric_n_chi );
		}
		initialized = true;
	}

	Size find_rot_aa = 1;
	while ( find_rot_aa <= rotameric_amino_acids.size() ) {
		if ( rotameric_amino_acids[ find_rot_aa ] == aa_in ) {
			return create_rotameric_dunlib( aa_in, rotameric_n_chi[ find_rot_aa ], ! (option[ corrections::score::dun10 ] || option[ corrections::score::dun08 ] )  );
		}
		++find_rot_aa;
	}

	Size find_semirot_aa = 1;
	while ( find_semirot_aa <= sraa.size() ) {
		if ( sraa[ find_semirot_aa ] == aa_in ) {
			Size i = find_semirot_aa;
			return create_semi_rotameric_dunlib( aa_in, srnchi[i], scind[i], sampind[i], sym[i], astr[i] );
		}
		++find_semirot_aa;
	}
	return 0;
}



/// @details Hard code this data in a single place; adjust the booleans for the backbone
/// independent
void
RotamerLibrary::initialize_dun08_aa_parameters(
	utility::vector1< chemical::AA > & rotameric_amino_acids,
	utility::vector1< Size > & rotameric_n_chi,
	utility::vector1< chemical::AA > & sraa,
	utility::vector1< Size > & srnchi,
	utility::vector1< bool > & scind,
	utility::vector1< bool > & sampind,
	utility::vector1< bool > & sym,
	utility::vector1< Real > & astr
)
{
	using namespace chemical;

	/// Rotameric specifics
	rotameric_amino_acids.reserve( 10 );         rotameric_n_chi.reserve( 10 );
	rotameric_amino_acids.push_back( aa_cys );   rotameric_n_chi.push_back( 1 );
	rotameric_amino_acids.push_back( aa_ile );   rotameric_n_chi.push_back( 2 );
	rotameric_amino_acids.push_back( aa_lys );   rotameric_n_chi.push_back( 4 );
	rotameric_amino_acids.push_back( aa_leu );   rotameric_n_chi.push_back( 2 );
	rotameric_amino_acids.push_back( aa_met );   rotameric_n_chi.push_back( 3 );
	rotameric_amino_acids.push_back( aa_pro );   rotameric_n_chi.push_back( 3 );
	rotameric_amino_acids.push_back( aa_arg );   rotameric_n_chi.push_back( 4 );
	rotameric_amino_acids.push_back( aa_ser );   rotameric_n_chi.push_back( 1 );
	rotameric_amino_acids.push_back( aa_thr );   rotameric_n_chi.push_back( 1 );
	rotameric_amino_acids.push_back( aa_val );   rotameric_n_chi.push_back( 1 );


	/// Semi rotameric specifics
	Size const n_semi_rotameric_aas = 8;
	sraa.resize( n_semi_rotameric_aas );
	srnchi.resize( n_semi_rotameric_aas );
	scind.resize( n_semi_rotameric_aas );
	sampind.resize( n_semi_rotameric_aas );
	sym.resize( n_semi_rotameric_aas );
	astr.resize( n_semi_rotameric_aas );


	Size i = 0;
	sraa[++i] = aa_asp; srnchi[i] = 1; scind[i] = 0; sampind[i] = 0; sym[i]=1; astr[i] =  -90;
	sraa[++i] = aa_glu; srnchi[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=1; astr[i] =  -90;
	sraa[++i] = aa_phe; srnchi[i] = 1; scind[i] = 0; sampind[i] = 0; sym[i]=1; astr[i] =  -30;
	sraa[++i] = aa_his; srnchi[i] = 1; scind[i] = 0; sampind[i] = 0; sym[i]=0; astr[i] = -180;
	sraa[++i] = aa_asn; srnchi[i] = 1; scind[i] = 0; sampind[i] = 0; sym[i]=0; astr[i] = -180;
	sraa[++i] = aa_gln; srnchi[i] = 2; scind[i] = 0; sampind[i] = 0; sym[i]=0; astr[i] = -180;
	sraa[++i] = aa_trp; srnchi[i] = 1; scind[i] = 0; sampind[i] = 0; sym[i]=0; astr[i] = -180;
	sraa[++i] = aa_tyr; srnchi[i] = 1; scind[i] = 0; sampind[i] = 0; sym[i]=1; astr[i] =  -30;

}

void
RotamerLibrary::initialize_dun02_aa_parameters(
	utility::vector1< chemical::AA > & rotameric_amino_acids,
	utility::vector1< Size > & rotameric_n_chi
)
{
	using namespace chemical;

	utility::vector1< Size > nchi_for_aa;
	Size n_rotameric_aas;
	initialize_dun02_aa_parameters( nchi_for_aa, n_rotameric_aas );
	rotameric_amino_acids.resize( 0 ); rotameric_amino_acids.reserve( n_rotameric_aas );
	rotameric_n_chi.resize( 0 );       rotameric_n_chi.reserve( n_rotameric_aas );

	for ( Size ii = 1; ii <= num_canonical_aas; ++ii ) {
		if ( nchi_for_aa[ ii ] != 0 ) {
			rotameric_amino_acids.push_back( AA(ii) );
			rotameric_n_chi.push_back( nchi_for_aa[ ii ] );
		}
	}

}

void
RotamerLibrary::initialize_dun02_aa_parameters(
	utility::vector1< Size > & nchi_for_aa,
	Size & n_rotameric_aas
)
{
	using namespace chemical;
	nchi_for_aa.resize( num_canonical_aas );
	std::fill( nchi_for_aa.begin(), nchi_for_aa.end(), 0 ); // overkill! ala, gly = 0...
	n_rotameric_aas = 0;

	nchi_for_aa[ aa_cys ] = 1; ++n_rotameric_aas;
	nchi_for_aa[ aa_asp ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_glu ] = 3; ++n_rotameric_aas;
	nchi_for_aa[ aa_phe ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_his ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_ile ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_lys ] = 4; ++n_rotameric_aas;
	nchi_for_aa[ aa_leu ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_met ] = 3; ++n_rotameric_aas;
	nchi_for_aa[ aa_asn ] = 2; ++n_rotameric_aas;
	nchi_for_aa[ aa_pro ] = 3; ++n_rotameric_aas;
	nchi_for_aa[ aa_gln ] = 3; ++n_rotameric_aas;
	nchi_for_aa[ aa_arg ] = 4; ++n_rotameric_aas;
	nchi_for_aa[ aa_ser ] = 1; ++n_rotameric_aas;
	nchi_for_aa[ aa_thr ] = 1; ++n_rotameric_aas;
	nchi_for_aa[ aa_val ] = 1; ++n_rotameric_aas;
	nchi_for_aa[ aa_trp ] = 2; ++n_rotameric_aas;
 	nchi_for_aa[ aa_tyr ] = 2; ++n_rotameric_aas;

}


template < Size T >
void
RotamerLibrary::initialize_and_read_srsrdl(
	SemiRotamericSingleResidueDunbrackLibrary< T > & srsrdl,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle,
	utility::io::izstream & rotamer_definitions,
	utility::io::izstream & regular_library,
	utility::io::izstream & continuous_minimization_bbdep
	// phasing out bbind sampling -- utility::io::izstream & rnchi_bbind_probabilities
) const
{
	initialize_srsrdl( srsrdl, nrchi_is_symmetric, nrchi_start_angle );

	srsrdl.read_from_files( rotamer_definitions, regular_library, continuous_minimization_bbdep );
}

template < Size T >
void
RotamerLibrary::initialize_srsrdl(
	SemiRotamericSingleResidueDunbrackLibrary< T > & srsrdl,
	bool const nrchi_is_symmetric,
	Real const nrchi_start_angle
) const
{
	Real const nrchi_periodicity     = nrchi_is_symmetric ? 180.0 : 360.0;
	Real const nrchi_bbdep_step_size = nrchi_is_symmetric ?   5.0 :  10.0;
	Real const nrchi_bbind_step_size = nrchi_is_symmetric ?   0.5 :   1.0;

	srsrdl.set_nrchi_periodicity( nrchi_periodicity );
	srsrdl.set_nonrotameric_chi_start_angle( nrchi_start_angle );
	srsrdl.set_nonrotameric_chi_bbdep_scoring_step_size( nrchi_bbdep_step_size );
	srsrdl.set_nonrotameric_chi_bbind_scoring_step_size( nrchi_bbind_step_size );
}

///////////////////////////////////////////////////////////////////////////////
//
// this is going to be wicked slow!!!
// apl -- I love Phil's comment above, and though it no longer makes sense to
// keep around, I have trouble deleting it.  Input is pretty fast now.
// This function seems unneccessary...
void
read_dunbrack_library(
	RotamerLibrary & rotamer_library
)
{
	rotamer_library.create_fa_dunbrack_libraries();
}


void
RotamerLibrary::write_to_file( std::string filename ) const
{
	utility::io::ozstream out( filename );
	for ( library_iterator it=libraries_.begin(), eit=libraries_.end(); it!=eit; ++it ) {
		//it->second->write_to_file( out );
	}
}

/// @brief Generic "write out a library to binary" which has a straightforward
/// nlibraries (aa, aa_data)* format.  Works for both 02 and 08 libraries.  The
/// kind of aa_data written is left to the discression of the SinResDunLib
/// derived classes.
void
RotamerLibrary::write_to_binary( utility::io::ozstream & out ) const
{

	using namespace boost;

	Size count_libraries( 0 );
	for ( library_iterator it=libraries_.begin(), eit=libraries_.end(); it!=eit; ++it ) {
		SingleResidueDunbrackLibrary const * srdl =
			dynamic_cast< SingleResidueDunbrackLibrary const * > ( it->second() );

		if ( srdl ) ++count_libraries; /// write out only the dunbrack libraries.

	}

	/// 1. How many libraries?
	boost::int32_t const nlibraries = count_libraries;
	out.write( (char*) & nlibraries, sizeof( boost::int32_t )  );

	for ( library_iterator it=libraries_.begin(), eit=libraries_.end(); it!=eit; ++it ) {
		SingleResidueDunbrackLibrary const * srdl =
			dynamic_cast< SingleResidueDunbrackLibrary const * > ( it->second() );

		if ( ! srdl ) continue; /// write out only the dunbrack libraries.

		/// 2. Amino acid type of next library
		boost::int32_t const which_aa( it->first );
		out.write( (char*) &which_aa, sizeof( boost::int32_t ) );

		/// 3. Data for this amino acid type.
		srdl->write_to_binary( out );
	}

}

/// @details Generic binary file reader.  Intended to work for both 02 and 08 Libraries.
/// If a "preamble" starts the binary format, the input stream to this function must have already
/// been walked past it.  Somewhat circularly, the preamble is anything in the binary file that
/// leads the "nlibraries (aa, aa_data)*" format perscribed in this function.
void
RotamerLibrary::read_from_binary( utility::io::izstream & in )
{

	using namespace boost;

	/// 1. How many libraries?
	boost::int32_t nlibraries( 0 );
	in.read( (char*) &nlibraries, sizeof( boost::int32_t )  );
	for ( Size ii = 1; ii <= Size( nlibraries ); ++ii ) {

		/// 2.  Which amino acid is next in the binary file?
		boost::int32_t which_aa32( chemical::aa_unk );
		in.read( (char*) &which_aa32, sizeof( boost::int32_t ) );
		AA which_aa( static_cast< AA >( which_aa32 ) );

		/// 3. Read the data associated with that amino acid.
		SingleResidueDunbrackLibraryOP single_lib = create_srdl( which_aa );
		if( !single_lib ){
			utility_exit_with_message( "Error reading single residue rotamer library for " + name_from_aa( which_aa ) );
		}
		single_lib->read_from_binary( in );
		add_residue_library( which_aa, single_lib );
	}
}

RotamerLibrary &
RotamerLibrary::get_instance()
{
	if ( rotamer_library_ == 0 ) {
		rotamer_library_ = new RotamerLibrary();
		read_dunbrack_library( *rotamer_library_ );
	}
	return *rotamer_library_;
}


SingleResidueRotamerLibraryCAP
RotamerLibrary::get_NCAARotamerLibrary( chemical::ResidueType const & rsd_type ) const
{
	// get some info about amino acid type
	std::string aa_name3( rsd_type.name3() );
	Size n_rotlib_chi( rsd_type.nchi() - rsd_type.n_proton_chi() );
	chemical::AA aan( rsd_type.aa() );
	bool dun02( true );

	if ( ncaa_rotlibs_.find( aa_name3 ) == ncaa_rotlibs_.end() ) {

		// create izstream from path
		std::string dir_name = basic::database::full_name( "/ncaa_rotlibs/" );
		std::string file_name = rsd_type.get_ncaa_rotlib_path();
		utility::io::izstream rotlib_in( dir_name + file_name );
		std::cout << "Reading in rot lib " << dir_name + file_name << "...";

		// get an instance of RotamericSingleResidueDunbrackLibrary, but need a RotamerLibrary to do it
		// this means that when ever you read in the NCAA libraries you will also read in the Dunbrack libraries
		// this may need to be a pointer to the full type and not just a SRRLOP
		std::string empty_string("");

		// this comes almost directally from RotmerLibrary.cc::create_rotameric_dunlib()
		SingleResidueRotamerLibraryOP ncaa_rotlib;
		RotamericSingleResidueDunbrackLibrary< ONE > * r1;
 		RotamericSingleResidueDunbrackLibrary< TWO > * r2;
 		RotamericSingleResidueDunbrackLibrary< THREE > * r3;
 		RotamericSingleResidueDunbrackLibrary< FOUR > * r4;

		switch ( n_rotlib_chi ) {
		case 1:
			r1 = new RotamericSingleResidueDunbrackLibrary< ONE >( aan, dun02 );
			r1->set_n_chi_bins( rsd_type.get_ncaa_rotlib_n_bin_per_rot() );
			r1->read_from_file( rotlib_in, false );
			ncaa_rotlib = r1; break;
		case 2:
			r2 = new RotamericSingleResidueDunbrackLibrary< TWO >( aan, dun02 );
			r2->set_n_chi_bins( rsd_type.get_ncaa_rotlib_n_bin_per_rot() );
			r2->read_from_file( rotlib_in, false );
			ncaa_rotlib = r2; break;
		case 3:
			r3 = new RotamericSingleResidueDunbrackLibrary< THREE >( aan, dun02 );
			r3->set_n_chi_bins( rsd_type.get_ncaa_rotlib_n_bin_per_rot() );
			r3->read_from_file( rotlib_in, false );
			ncaa_rotlib = r3; break;
		case 4:
			r4 = new RotamericSingleResidueDunbrackLibrary< FOUR >( aan, dun02 );
			r4->set_n_chi_bins( rsd_type.get_ncaa_rotlib_n_bin_per_rot() );
			r4->read_from_file( rotlib_in, false );
			ncaa_rotlib = r4; break;
		default:
			utility_exit_with_message( "ERROR: too many chi angles desired for ncaa library: " + n_rotlib_chi );
			break;
		}

		// add new rotamer library to map
		ncaa_rotlibs_[ aa_name3 ] = ncaa_rotlib;
		std::cout << "done!" << std::endl;
	}
	return ( ncaa_rotlibs_.find( aa_name3 )->second)();
}


} // dunbrack
} // namespace scoring
} // namespace core
