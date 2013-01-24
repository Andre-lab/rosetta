// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Ingemar Andre, Phil Bradley

// Unit Headers
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Package headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/id/types.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>

#include <core/conformation/Conformation.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <cassert>
#include <iostream>

// core utilities
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.conformation.SymmetryInfo");

namespace core {
namespace conformation {
namespace symmetry {

SymmetryInfo::SymmetryInfo() { use_symmetry_ = false; score_multiply_factor_ = 1; }
SymmetryInfo::~SymmetryInfo() {}


bool SymmetryInfo::operator==( SymmetryInfo const & s )
{
	return ( npseudo_ == s.npseudo_ &&
           bb_clones_ ==s.bb_clones_ &&
           chi_clones_ == s.chi_clones_ &&
           jump_clones_ == s.jump_clones_ );
}

bool SymmetryInfo::operator!=( SymmetryInfo const & s )
{
	return !( *this == s );
}

SymmetryInfo::SymmetryInfo( SymmData const & symm_data, Size const nres_subunit, Size const njump_subunit )
{
	Size joff = njump_subunit*symm_data.get_subunits();
	std::map<std::string,Size> const & name2num = symm_data.get_jump_string_to_jump_num();
	for(std::map<std::string,Size>::const_iterator i = name2num.begin(); i != name2num.end(); ++i) {
		dofname2jnum_[i->first] = i->second+joff;
		jnum2dofname_[i->second+joff] = i->first;
	}
	if (  symm_data.get_jump_clones().size() > 0 ) {
		initialize( nres_subunit, njump_subunit,
			symm_data.get_subunits(), symm_data.get_num_virtual(),
			symm_data.get_jump_clones(), symm_data.get_dofs(),
			symm_data.get_score_subunit(), symm_data.get_score_multiply_subunit(),
			symm_data.get_slide_info(), symm_data.get_interfaces() );
	} else {
		initialize( nres_subunit, njump_subunit,
			symm_data.get_subunits(), symm_data.get_num_virtual(),
			symm_data.get_dofs(),symm_data.get_score_subunit(),
			symm_data.get_score_multiply_subunit(), symm_data.get_slide_info(),
			symm_data.get_interfaces() );
	}
	if(symm_data.get_num_components()==1) TR       << *this << std::endl;
	else                                  TR.Debug << *this << std::endl;
}

SymmetryInfo::SymmetryInfo(
	Size const nres_monomer,
	Size const njump_monomer,
	Size const N,
	std::map< Size, SymDof > dofs,
	Size const score_subunit,
	utility::vector1< Size > score_multiply_subunit,
	SymSlideInfo slide_info,
	Size const num_interfaces,
	std::string const & type
)
{
	initialize( nres_monomer, njump_monomer, N, N, dofs, score_subunit,
		score_multiply_subunit, slide_info, num_interfaces, type );
}

///@details make a copy of this SymmetryInfo ( allocate actual memory for it )
SymmetryInfoOP
SymmetryInfo::clone() const
{
  return new SymmetryInfo( *this );
}

// This is an old style constructor. Should change soon...
void
SymmetryInfo::initialize(
	Size const nres_monomer,
	Size const njump_monomer,
	Size const n_subunits,
	Size const num_virtual,
	std::map< Size, SymDof > dofs,
	Size const score_subunit,
	utility::vector1< Size > score_multiply_subunit,
	SymSlideInfo slide_info,
	Size const num_interfaces,
	std::string const & type
)
{
	nres_monomer_ = nres_monomer;

	// set number of interfaces
	interfaces_ = num_interfaces;
	// store the score multiplication factors
	set_score_multiply_from_subunit_factors(score_multiply_subunit, nres_monomer, n_subunits);
	// store the number of monomer jumps
	njump_monomer_ = njump_monomer;
	// store the allowed dofs
	dofs_ = dofs;
	// set use symmetry
	use_symmetry_ = true;
	// set cp_weighting_during_minimization_ to false
	cp_weighting_during_minimization_ = false;
	// slide info
	slide_info_ = slide_info;
	// store type
	type_ = type;
	//scoring subunit
	scoring_subunit_ = score_subunit;
	// setup bb,chi clones
	bb_clones_.clear();
	chi_clones_.clear();
	jump_clones_.clear();

	//check that score_monomer makes sense...
	if ( score_subunit > n_subunits || score_subunit < 1 ) {
      utility_exit_with_message("score_subunit must be in the range 1-N");
   }

	//special case of no symmetry
	if ( type == "c1" ) {
		npseudo_ = num_virtual;

		//we need to map to an empty array in order for
		for ( Size i=1; i<= nres_monomer; ++i ) {
			Clones clones;
			clones.clear();
			bb_clones_.insert( std::make_pair( i, clones ) );
			chi_clones_.insert( std::make_pair( i, clones ) );
		}
		for ( Size i=1; i<= njump_monomer; ++i ) {
			Clones clones;
			clones.clear();
			jump_clones_.insert( std::make_pair( i, clones ) );
		}
		return;
	}//end c1 symmetry

	for ( Size i=1; i<= nres_monomer; ++i ) {
		Clones clones;
		int base ( i + ( score_subunit - 1 ) * nres_monomer );
		for ( Size k=0; k<n_subunits; ++k ) {
			if ( k+1 != score_subunit ) {
				clones.push_back( i + k * nres_monomer );
				add_bb_clone( base, i + k * nres_monomer );
				add_chi_clone( base, i + k * nres_monomer );
			}
		}
		bb_clones_.insert( std::make_pair( base, clones ) );
		chi_clones_.insert( std::make_pair( base, clones ) );
	}

	// the N*njump_monomer internal jumps
	for ( Size i=1; i<= njump_monomer; ++i ) {
		//Clones clones;
		for ( Size k=1; k<n_subunits; ++k ) {
			//clones.push_back( i + k * njump_monomer );
			add_jump_clone( i, i + k * njump_monomer, 0.0 );
		}
		//jump_clones_.insert( std::make_pair( i, clones ) );
	}

	if ( type == "no_pseudo" ) {
		npseudo_ = num_virtual;

	} else if ( type == "simple" ) {
		// 1                 --> N*njump_monomer  : the internal jumps
		// N*njump_monomer+1 --> N*njump_monomer+N: the pseudo-rsd--monomer jumps
		// last N-1 jumps                         : jumps between pseudo-rsds

		npseudo_ = num_virtual;

		// the N jumps from pseudo-residues to monomers
		{
			Size const base_jump( n_subunits*njump_monomer + 1 );
			//Clones clones;
			for ( Size k=1; k<n_subunits; ++k ) {
				//clones.push_back( base_jump + k );
				add_jump_clone( base_jump, base_jump + k, 0.0 );
			}
			//jump_clones_.insert( std::make_pair( base_jump, clones ) );
		}


		// the N-1 jumps between pseudo-residues
		{
/*			Size const base_jump( N*njump_monomer + N + 1 );
			//Clones clones;
			for ( Size k=1; k<N-1; ++k ) {
				//clones.push_back( base_jump + k );
				add_jump_clone( base_jump, base_jump + k, 0.0 );
			}
			//jump_clones_.insert( std::make_pair( base_jump, clones ) ); */
		}
	} else {
		std::cerr << "unrecognized type: " << type << std::endl;
		utility_exit();
	}

	update_score_multiply_factor();
}

/// @details  This is a helper function for some of the DOF_ID routines below. Really just a best guess...
/// this is a little tricky: the mapping from a DOF_ID to a TorsionID is not straightforward to
/// construct (see kinematics/util.cc:setup_dof_to_torsion_map)
/// So we don't really know whether a dof_id is a bb degree of freedom or a chi degree of freedom...
/// or even which residue it should be attached to, eg the phi of residue i might be dof_id with rsd i-1
/// So we take a guess based on whether id.atomno is a backbone or sidechain atom
id::TorsionType
guess_torsion_type_of_dof_id( id::DOF_ID const & id, Conformation const & conf )
{
	if ( id::RB1 <= id.type() && id.type() <= id::RB6 ) {
		return id::JUMP;
	} else {
		if ( conf.atom_is_backbone_norefold( id.rsd(), id.atomno() ) ) {
			return id::BB;
		} else {
			return id::CHI;
		}
	}
}

//
void
SymmetryInfo::update_score_multiply_factor()
{
	// compute the score_multiply_factor
	utility::vector1< bool > indep_res = independent_residues();
	for (int i=1; i<=(int)indep_res.size(); ++i) {
		if (indep_res[i]) {
			score_multiply_factor_ = score_multiply_[ i ];
			//std::cout<< "score_multiply_factor "<< score_multiply_factor_ <<std::endl;
			break;
		}
	}
}

// This is an old style constructor. Should change soon...
void
SymmetryInfo::initialize(
	Size const nres_monomer,
	Size const njump_monomer,
	Size const n_subunits,
	Size const num_virtual,
	std::map< Size, WtedClones > jump_clones,
	std::map< Size, SymDof > dofs,
	Size const score_subunit,
	utility::vector1< Size > score_multiply_subunit,
	SymSlideInfo slide_info,
	Size const num_interfaces,
	std::string const & type
)
{
	nres_monomer_ = nres_monomer;

	// set number of interfaces
	interfaces_ = num_interfaces;
	// store the score multiplication factors
	set_score_multiply_from_subunit_factors(score_multiply_subunit, nres_monomer, n_subunits);

	// store the number of monomer jumps
	njump_monomer_ = njump_monomer;
	// store the allowed dofs
	dofs_ = dofs;
	// set use symmetry
	use_symmetry_ = true;
	// set cp_weighting_during_minimization_ to false
	cp_weighting_during_minimization_ = false;
	// slide info
	slide_info_ = slide_info;
	// store type
	type_ = type;
	//scoring subunit
	scoring_subunit_ = score_subunit;
	// setup bb,chi clones
	bb_clones_.clear();
	chi_clones_.clear();
	jump_clones_.clear();

	// 1                 --> N*njump_monomer  : the internal jumps
	// N*njump_monomer+1 --> N*njump_monomer+N: the pseudo-rsd--monomer jumps
	// last N-1 jumps                         : jumps between pseudo-rsds

	npseudo_ = num_virtual;

	for ( Size i=1; i<= nres_monomer; ++i ) {
		Clones clones;
      int base ( i + ( score_subunit - 1 ) * nres_monomer );
      for ( Size k=0; k<n_subunits; ++k ) {
        if ( k+1 != score_subunit ) {
          clones.push_back( i + k * nres_monomer );
					add_bb_clone( base, i + k * nres_monomer );
		      add_chi_clone( base, i + k * nres_monomer );
        }
      }
      bb_clones_.insert( std::make_pair( base, clones ) );
      chi_clones_.insert( std::make_pair( base, clones ) );
	}

	// the N*njump_monomer internal jumps
	for ( Size i=1; i<= njump_monomer; ++i ) {
		for ( Size k=0; k<n_subunits; ++k ) {
			if (k != ( score_subunit - 1 ) )
				add_jump_clone( i + (score_subunit-1)*njump_monomer, i + k*njump_monomer, 0.0 );
		}
		//jump_clones_.insert( std::make_pair( i, clones ) );
	}

	std::map< Size,WtedClones >::const_iterator it, it_start=jump_clones.begin(), it_end=jump_clones.end();
	for ( it=it_start; it != it_end; ++it ) {
		//Clones clones;
		for ( Size i = 1; i<= it->second.size(); ++i ) {
			//clones.push_back( it->second[i] + N*njump_monomer );
			add_jump_clone( it->first + n_subunits*njump_monomer, it->second[i].first + n_subunits*njump_monomer, it->second[i].second );
		}
		//jump_clones_.insert( std::make_pair( it->first + N*njump_monomer, clones ) );
	}

	// compute the score_multiply_factor
	update_score_multiply_factor();
}

/////////////////////////////////////////////////////////////////////////////
template< class T >
void
comma_strings_to_vector_map(
	std::istream & is,
	Size const nbase,
	std::map< Size, utility::vector1< T > > & clones,
	std::string tag=""
)
{
	bool fail( false );
	std::string tag0;
	if( tag != "" ) {
		is >> tag0;
		if( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}
	for ( Size i=1; !fail && i<= nbase; ++i ) {
		std::string jump_string;
		is >> jump_string;
		if ( is.fail() ) {
			fail = true;
			break;
		}
		std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
		std::istringstream l( jump_string );
		Size base_jump;
		l >> base_jump;
		if ( l.fail() ) {
			fail = true;
			break;
		}
		while ( true ) {
			T j;
			l >> j;
			if ( l.fail() ) break;
			clones[ base_jump ].push_back( j );
		}
		if ( clones[ base_jump ].size() < 1 ) {
			fail = true;
			break;
		}
	}
	if ( fail ) {
		is.setstate( std::ios_base::failbit );
	}
}

/////////////////////////////////////////////////////////////////////////////
template <class T>
void
comma_strings_to_map(
	std::istream & is,
	Size const nbase,
	std::map< Size, T > & clones,
	std::string tag=""
) {
	bool fail( false );
	std::string tag0;
	if( tag != "" ) {
		is >> tag0;
		if( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}
	for ( Size i=1; !fail && i<= nbase; ++i ) {
		std::string jump_string;
		is >> jump_string;
		if ( is.fail() ) {
			fail = true;
			break;
		}
		std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
		std::istringstream l( jump_string );
		Size base_jump;
		l >> base_jump;
		if ( l.fail() ) {
			fail = true;
			break;
		}
		l >> clones[ base_jump ];
	}
	if ( fail ) {
		is.setstate( std::ios_base::failbit );
	}
}

/////////////////////////////////////////////////////////////////////////////
void
comma_strings_to_map(
	std::istream & is,
	Size const nbase,
	std::map< Size, SymDof > & clones,
	std::string tag=""
)
{
	bool fail( false );
	std::string tag0;
	if( tag != "" ) {
		is >> tag0;
		if( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}
	for ( Size i=1; !fail && i<= nbase; ++i ) {
		std::string jump_string;
		is >> jump_string;
		if ( is.fail() ) {
			fail = true;
			break;
		}
		std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
		std::istringstream l( jump_string );
		Size base_jump;
		l >> base_jump;
		std::string dof_line;
		l >> dof_line;
		clones[base_jump].read(dof_line);
		if ( l.fail() ) {
			fail = true;
			break;
		}
	}
	if ( fail ) {
		is.setstate( std::ios_base::failbit );
	}
}

/////////////////////////////////////////////////////////////////////////////
void
comma_strings_to_vector(
	std::istream & is,
	Size const nbase,
	utility::vector1< Size > & clones,
	std::string tag=""
)
{
	bool fail( false );
	std::string tag0;
	if( tag != "" ) {
		is >> tag0;
		if( tag0 != tag ) {
			TR << "Input failed: tag mismatch " << tag << " " << tag0 << std::endl;
			return;
		}
	}

	std::string jump_string;
	is >> jump_string;
	if ( is.fail() ) fail = true;
	std::replace( jump_string.begin(), jump_string.end(), ',', ' ' );
	std::istringstream l( jump_string );
	while ( true ) {
		Size j;
		l >> j;
		if ( l.fail() ) break;
		clones.push_back( j );
	}
	if ( clones.size() != nbase ) {
		fail = true;
	}

	if ( fail ) {
		is.setstate( std::ios_base::failbit );
	}
}

/////////////////////////////////////////////////////////////////////////////
template<class CloneType>
void vector_map_to_comma_strings(
	std::ostream & out,
	std::map< Size, utility::vector1< CloneType > > clones,
	std::string tag=""
) {
	typename std::map< Size,utility::vector1<CloneType> >::const_iterator it;
	if( tag != "" ) out << ' ' << tag ;
	for ( it = clones.begin(); it != clones.end(); ++it ) {
		out << ' ' << it->first;
		utility::vector1< CloneType > const & l( it->second );
		for ( Size i=1; i<= l.size(); ++i ) {
			out << ',' << l[i];
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
template< class CloneType >
void map_to_comma_strings(
	std::ostream & out,
	std::map< Size, CloneType > clones,
	std::string tag=""
) {
	typename std::map< Size , CloneType >::const_iterator it;
	if( tag != "" ) out << ' ' << tag ;
	for ( it = clones.begin(); it != clones.end(); ++it ) {
		out << ' ' << it->first << ',' << it->second ;
	}
}

/////////////////////////////////////////////////////////////////////////////
void
map_to_comma_strings(
	std::ostream & out,
	std::map< Size, SymDof > clones,
	std::string tag=""
)
{
	if( tag != "" ) out << ' ' << tag ;
	for ( std::map< Size , SymDof >::const_iterator
					it = clones.begin(); it != clones.end(); ++it ) {
		//Dof const & dof (it->second);
		out << " " << it->first << "," << it->second ;
	}
}

/////////////////////////////////////////////////////////////////////////////
void
vector_to_comma_strings(
	std::ostream & out,
	utility::vector1 < Size > clones,
	std::string tag=""
)
{
	if( tag != "" ) out << ' ' << tag ;
	for ( Size i=1; i<= clones.size(); ++i ) {
		if( i == 1 ) {
			out << " " << clones[i] ;
		} else {
			out << ',' << clones[i] ;
		}
	}

}

/////////////////////////////////////////////////////////////////////////////
std::istream& operator>> ( std::istream & s, SymmetryInfo & symminfo )
{
	bool fail( false );

	std::string tag;
	Size num_bb_indep, num_chi_indep, num_jump_indep;
	Size num_bb_dep, num_chi_dep, num_jump_dep;
	Size num_dof, num_score_multiply;

	symminfo.set_use_symmetry(true);
	symminfo.set_cp_weighting_during_minimization(false);

	s >> tag ;
	if ( tag != "SYMMETRY_INFO" || s.fail() ) {
		fail = true;
	} else {
		s >> tag >> tag >> tag;

		//fpd try to ensure backwards compatibility
		bool old_stream = false;
		if (tag == "N_RES_MONOMER") {
			s >> symminfo.nres_monomer_
			  >> tag >> symminfo.scoring_subunit_
			  >> tag >> symminfo.njump_monomer_
			  >> tag;
		} else {
			old_stream = true;
		}
		s >> symminfo.npseudo_
			>> tag >> symminfo.interfaces_
		  >> tag >> symminfo.type_
		  >> tag >> num_bb_indep
		  >> tag >> num_chi_indep
		  >> tag >> num_jump_indep
		  >> tag >> num_bb_dep
		  >> tag >> num_chi_dep
		  >> tag >> num_jump_dep
		  >> tag >> num_dof
		  >> tag >> num_score_multiply;

		if ( s.fail() ) fail = true;

		// clones
		comma_strings_to_vector_map( s,   num_bb_indep,  symminfo.bb_clones_, "BB_CLONES" );
		comma_strings_to_vector_map( s,  num_chi_indep,  symminfo.chi_clones_, "CHI_CLONES" );
		comma_strings_to_vector_map( s, num_jump_indep,  symminfo.jump_clones_, "JUMP_CLONES" );

		if (old_stream) {
			TR << "Warning: Symmetric input stream is out of date! Trying to recover." << std::endl;
			// set master jumps to 1; clones to 0
			for (std::map<Size,SymmetryInfo::Clones>::const_iterator map_it=symminfo.jump_clones_.begin(),
			     map_end=symminfo.jump_clones_.end();
			     map_it != map_end; ++map_it) {
				for (Size i=1; i<=map_it->second.size(); ++i) {
					symminfo.jump_clone_wts_[ map_it->second[i] ] = 0;
				}
			}

			// guess at missing parameters
			symminfo.nres_monomer_ = num_bb_indep;
			std::map<Size,SymmetryInfo::Clones>::const_iterator first_bb_clone=symminfo.bb_clones_.begin();
			symminfo.scoring_subunit_ = 1 + ((first_bb_clone->first-1) / symminfo.nres_monomer_);
			symminfo.njump_monomer_ = 0;
		} else {
			comma_strings_to_map( s, num_jump_dep-num_jump_indep,  symminfo.jump_clone_wts_, "JUMP_CLONE_WEIGHTS" );
			old_stream |= (symminfo.jump_clone_wts_.size() == 0);
		}

		// follows
		comma_strings_to_map( s,   num_bb_dep,  symminfo.bb_follows_, "BB_FOLLOWS" );
		comma_strings_to_map( s,  num_chi_dep,  symminfo.chi_follows_, "CHI_FOLLOWS" );
		comma_strings_to_map( s, num_jump_dep,  symminfo.jump_follows_, "JUMP_FOLLOWS" );

		// dof_
		comma_strings_to_map( s, num_dof,  symminfo.dofs_, "DOFS" );

		// score_multiply_
		comma_strings_to_vector( s, num_score_multiply,  symminfo.score_multiply_, "SCORE_MULTIPLY" );
		symminfo.update_score_multiply_factor();

		//
		symminfo.set_use_symmetry( true );

		if ( fail ) {
			std::cout << "Symmetry_info operator>>: Input failed" << std::endl;
			s.setstate( std::ios_base::failbit );
			return s;
		}
	}

	return s;
}

/////////////////////////////////////////////////////////////////////////////
std::ostream& operator<< ( std::ostream & s, const SymmetryInfo & symminfo )
{
	s << "SYMMETRY_INFO " <<
		"N " << symminfo.subunits() << ' ' <<
		"N_RES_MONOMER " << symminfo.nres_monomer_ << ' ' <<
		"SCORING_SUBUNIT " << symminfo.scoring_subunit_ << ' ' <<
		"N_JUMP_MONOMER " << symminfo.njump_monomer_ << ' ' <<
		"N_VIRT " << symminfo.npseudo_ << ' ' <<
		"N_INTERFACE " << symminfo.num_interfaces() << ' ' <<
		"TYPE " << symminfo.type_ << ' ' <<
		"BB_CLONES_SIZE " << symminfo.bb_clones_.size() << ' ' <<
		"CHI_CLONES_SIZE " << symminfo.chi_clones_.size() << ' ' <<
		"JUMP_CLONES_SIZE " << symminfo.jump_clones_.size() << ' ' <<
		"BB_FOLLOWS_SIZE " << symminfo.bb_follows_.size() << ' ' <<
		"CHI_FOLLOWS_SIZE " << symminfo.chi_follows_.size() << ' ' <<
		"JUMP_FOLLOWS_SIZE " << symminfo.jump_follows_.size() << ' ' <<
		"DOFS_SIZE " << symminfo.dofs_.size() << ' ' <<
		"SCORE_MULTIPLY_SIZE " << symminfo.score_multiply_.size() ;

	// clones
	vector_map_to_comma_strings( s, symminfo.bb_clones_, "BB_CLONES" );
	vector_map_to_comma_strings( s, symminfo.chi_clones_, "CHI_CLONES" );
	vector_map_to_comma_strings( s, symminfo.jump_clones_, "JUMP_CLONES" );
	map_to_comma_strings( s, symminfo.jump_clone_wts_, "JUMP_CLONE_WEIGHTS" );

	// follows
	map_to_comma_strings( s, symminfo.bb_follows_, "BB_FOLLOWS" );
	map_to_comma_strings( s, symminfo.chi_follows_, "CHI_FOLLOWS" );
	map_to_comma_strings( s, symminfo.jump_follows_, "JUMP_FOLLOWS" );

	//dof
	map_to_comma_strings( s, symminfo.dofs_, "DOFS" );

	//score_multiply_
	vector_to_comma_strings( s, symminfo.score_multiply_, "SCORE_MULTIPLY" );

	return s;
}

bool
SymmetryInfo::write_silent_struct(
std::string const & filename
)
{
	bool success = false;

	utility::io::ozstream output;
	if ( !utility::file::file_exists( filename ) ) {
		output.open( filename );
	} else {
		output.open_append( filename );
	}

	output << *this << '\n';

	output.close();

	success = true;
	return success;
}


bool
SymmetryInfo::read_silent_struct(
std::string const & filename
)
{
	bool success = false;

	utility::io::izstream input ( filename.c_str() );
	std::istringstream line_stream;
	std::string line("");
	if ( !input ) {
		std::cerr << "ERROR:: Unable to open symmetry info file: "
					<< filename << std::endl;
		return success;
	}

	while( !input.eof() ) {
		getline(input,line);
		line_stream.clear();
		line_stream.str(line);
		line_stream >> *this;
	}

	input.close();

	success = true;
	return success;
}

bool
SymmetryInfo::is_virtual( Size const seqpos ) const {
	return ( seqpos > num_total_residues_without_pseudo() );
}

Size
SymmetryInfo::bb_follows( Size const seqpos ) const
{
	std::map< Size, Size >::const_iterator it( bb_follows_.find( seqpos ) );
	return ( it == bb_follows_.end() ? 0 : it->second );
}

Size
SymmetryInfo::chi_follows( Size const seqpos ) const
{
	std::map< Size, Size >::const_iterator it( chi_follows_.find( seqpos ) );
	return ( it == chi_follows_.end() ? 0 : it->second );
}

Size
SymmetryInfo::jump_follows( Size const seqpos ) const
{
	std::map< Size, Size >::const_iterator it( jump_follows_.find( seqpos ) );
	return ( it == jump_follows_.end() ? 0 : it->second );
}

std::vector < std::pair < Size, Size > >
SymmetryInfo::map_symmetric_res_pairs( Size res1, Size res2 )
{
	std::vector < std::pair < Size, Size > > map;
	int delta ( res2 - res1 );
	int mapped_res;
	for ( std::vector< Size>::const_iterator
        clone     = bb_clones( res1 ).begin(),
        clone_end = bb_clones( res1 ).end();
        clone != clone_end; ++clone ){
		if ( *clone + delta > num_total_residues() ) {
			mapped_res = (*clone + delta)%num_total_residues();
		} else {
			mapped_res = *clone + delta;
		}
		if ( mapped_res < 0 )
			mapped_res += num_independent_residues();
		map.push_back( std::make_pair( *clone, mapped_res ) );
	}
	return map;
}

bool
SymmetryInfo::bb_is_independent( Size const seqpos ) const
{
	return bb_follows(seqpos) == 0;
}

bool
SymmetryInfo::chi_is_independent( Size const seqpos ) const
{
  return chi_follows(seqpos) == 0;
}

bool
SymmetryInfo::fa_is_independent( Size const seqpos ) const
{
	return ( bb_is_independent(seqpos) && chi_is_independent(seqpos) );
}

bool
SymmetryInfo::jump_is_independent( Size const seqpos ) const
{
  return jump_follows(seqpos) == 0;
}

Size
SymmetryInfo::subunits() const
{
	return num_bb_clones() + 1;
}

utility::vector1< bool >
SymmetryInfo::independent_residues() const
{
	utility::vector1 < bool > residues;
	for ( Size i=1; i <=num_total_residues_with_pseudo(); ++i ){
		if ( bb_is_independent(i) )
			residues.push_back(true);
		else
			residues.push_back(false);
  }
  return residues;
}

Size
SymmetryInfo::num_bb_clones() const
{
	// all these lists have the same size
	if ( bb_clones_.empty() ) {
		return 0;
	}

	return bb_clones_.begin()->second.size();
}

Size
SymmetryInfo::num_chi_clones() const
{
  // all these lists have the same size
  return chi_clones_.begin()->second.size();
}


Size
SymmetryInfo::num_jump_clones() const
{
  // all these lists have the same size
  return jump_clones_.begin()->second.size();
}

Size
SymmetryInfo::num_independent_residues() const
{
	return bb_clones_.size();
}

Size
SymmetryInfo::num_total_residues() const
{
	return num_independent_residues()*( num_bb_clones() + 1 );
}

Size
SymmetryInfo::num_total_residues_with_pseudo() const
{
	return num_independent_residues()*( num_bb_clones() + 1 ) + npseudo_;
}

Size
SymmetryInfo::num_total_residues_without_pseudo() const
{
  return num_independent_residues()*( num_bb_clones() + 1 );
}

Size
SymmetryInfo::num_interfaces() const
{
	return interfaces_;
}

Size
SymmetryInfo::score_multiply_factor() const
{
	return score_multiply_factor_;
}

Size
SymmetryInfo::num_virtuals() const
{
		return npseudo_;
}

// bool
// SymmetryInfo::scoring_residue(Size residue ) const
// {
// 	if ( score_multiply( residue ) == 0 ) return false;
// 	return true;
// }

SymmetryInfo::Clones const &
SymmetryInfo::bb_clones( Size const seqpos ) const
{
	std::map< Size, Clones >::const_iterator it( bb_clones_.find( seqpos ) );
  if ( it == bb_clones_.end() ) {
  	return empty_list;
  }
  return it->second;
}

SymmetryInfo::Clones const &
SymmetryInfo::chi_clones( Size const seqpos ) const
{
  std::map< Size, Clones >::const_iterator it( chi_clones_.find( seqpos ) );
  if ( it == chi_clones_.end() ) {
    return empty_list;
  }
  return it->second;
}

SymmetryInfo::Clones const &
SymmetryInfo::jump_clones( Size const seqpos ) const
{
  std::map< Size, Clones >::const_iterator it( jump_clones_.find( seqpos ) );
  if ( it == jump_clones_.end() ) {
    return empty_list;
  }
  return it->second;
}

//fpd remap bb_clones/chi_clones when the ASU size changes
//fpd this recreates the arrays from scratch so it may be somewhat inefficient
void
SymmetryInfo::resize_asu( Size nres_new ) {
	if (nres_new == nres_monomer_) return; // nothing to do

	Size N = subunits();

	nres_monomer_ = nres_new;
	bb_clones_.clear();
	bb_follows_.clear();
	chi_clones_.clear();
	chi_follows_.clear();

	// make empty clones array
	for ( Size i=1; i<= nres_monomer_; ++i ) {
		Clones clones;
		clones.clear();
		bb_clones_.insert( std::make_pair( i, clones ) );
		chi_clones_.insert( std::make_pair( i, clones ) );
	}

	for ( Size i=1; i<= nres_monomer_; ++i ) {
		Clones clones;
		int base ( i + ( scoring_subunit_ - 1 ) * nres_monomer_ );
		for ( Size k=0; k<N; ++k ) {
			if ( k+1 != scoring_subunit_ ) {
				clones.push_back( i + k * nres_monomer_ );
				add_bb_clone( base, i + k * nres_monomer_ );
				add_chi_clone( base, i + k * nres_monomer_ );
			}
		}
		bb_clones_.insert( std::make_pair( base, clones ) );
		chi_clones_.insert( std::make_pair( base, clones ) );
	}
}


//fpd remap jump_clones when the number of monomer jumps changes
void
SymmetryInfo::update_nmonomer_jumps( Size njump_monomer ) {
	if (njump_monomer == njump_monomer_) return; // nothing to do

	//std::cerr << "SymmetryInfo::update_nmonomer_jumps(" << njump_monomer << ")  [old=" << njump_monomer_ << "]\n";
	Size N = subunits();

	// remember previous
	std::map< Size, Clones > old_jump_clones = jump_clones_;
	std::map< Size, Size > old_jump_follows = jump_follows_;
	std::map< Size, Real > old_jump_clone_weights = jump_clone_wts_;
	Size old_njump_monomer = njump_monomer_;

	njump_monomer_ = njump_monomer;
	jump_clones_.clear();
	jump_follows_.clear();
	jump_clone_wts_.clear();

	// make new monomer jumps from scratch
	// the N*njump_monomer internal jumps
	for ( Size i=1; i<= njump_monomer; ++i ) {
		for ( Size k=0; k<N; ++k ) {
			if (k != ( scoring_subunit_ - 1 ) )
				add_jump_clone( i + (scoring_subunit_-1)*njump_monomer_, i + k*njump_monomer_, 0.0 );
		}
		//jump_clones_.insert( std::make_pair( i, clones ) );
	}

	// 1                 --> N*njump_monomer  : the internal jumps
	// N*njump_monomer+1 --> N*njump_monomer+N: the pseudo-rsd--monomer jumps
	// last N-1 jumps                         : jumps between pseudo-rsds
	for ( std::map<Size,Clones>::const_iterator it=old_jump_clones.begin(), it_end = old_jump_clones.end();
	      it != it_end; it++) {
		Size source = it->first;
		Clones target = it->second;

		if (source > N*old_njump_monomer) {
			// a symm jump
			Size new_source = source + N*( njump_monomer - old_njump_monomer );
			for (Size i=1; i<=target.size(); ++i) {
				add_jump_clone( new_source, target[i] + N*(njump_monomer-old_njump_monomer), old_jump_clone_weights[target[i]] );
				//std::cerr << "Map (" << source << " , " << target[i] << ") to (" << new_source << " , " << target[i] + N*(njump_monomer-old_njump_monomer) << ")\n";
			}
		}
	}

	// dofs
	std::map< Size, SymDof > dofs_new;
	for ( std::map< Size, SymDof >::iterator it = dofs_.begin(), it_end = dofs_.end(); it!=it_end; ++it ) {
		dofs_new.insert( std::make_pair( it->first + N*( njump_monomer - old_njump_monomer ), it->second ) );
	}
	set_dofs( dofs_new );
}

void
SymmetryInfo::add_bb_clone( Size const base_pos, Size const clone_pos )
{
	if ( bb_follows_[ base_pos ] != 0 ) {
		std::cerr << "Error: add_bb_clone: base_pos is not independent: " <<
			base_pos << ' ' << bb_follows_[ base_pos ] << std::endl;
		utility_exit();
	}
	if ( bb_follows_[ clone_pos ] != 0 &&
			 bb_follows_[ clone_pos ] != base_pos ) {
		std::cerr << "Error: add_bb_clone: clone_pos already a follower: " <<
			clone_pos << ' ' << bb_follows_[ clone_pos ] << ' ' << base_pos <<
			std::endl;
		utility_exit();
	}

	bb_follows_[ clone_pos ] = base_pos;
	bb_clones_[ base_pos ].push_back( clone_pos );
}

void
SymmetryInfo::add_chi_clone( Size const base_pos, Size const clone_pos )
{
	if ( chi_follows_[ base_pos ] != 0 ) {
		std::cerr << "Error: add_chi_clone: base_pos is not independent: " <<
			base_pos << ' ' << chi_follows_[ base_pos ] << std::endl;
		utility_exit();
	}
	if ( chi_follows_[ clone_pos ] != 0 &&
			 chi_follows_[ clone_pos ] != base_pos ) {
		std::cerr << "Error: add_chi_clone: clone_pos already a follower: " <<
			clone_pos << ' ' << chi_follows_[ clone_pos ] << ' ' << base_pos <<
			std::endl;
		utility_exit();
	}

	chi_follows_[ clone_pos ] = base_pos;
	chi_clones_[ base_pos ].push_back( clone_pos );
}

void
SymmetryInfo::add_jump_clone( Size const base_pos, Size const clone_pos, Real const jump_wt )
{
	if ( jump_follows_[ base_pos ] != 0 ) {
		std::cerr << "Error: add_jump_clone: base_pos is not independent: " <<
			base_pos << ' ' << bb_follows_[ base_pos ] << std::endl;
		utility_exit();
	}
	if ( jump_follows_[ clone_pos ] != 0 &&
			 jump_follows_[ clone_pos ] != base_pos ) {
		std::cerr << "Error: add_jump_clone: clone_pos already a follower: " <<
			clone_pos << ' ' << jump_follows_[ clone_pos ] << ' ' << base_pos <<
			std::endl;
		utility_exit();
	}

	jump_follows_[ clone_pos ] = base_pos;
	jump_clones_[ base_pos ].push_back( clone_pos );
	jump_clone_wts_[ clone_pos ] = jump_wt;
}

std::map< Size, SymDof > const &
SymmetryInfo::get_dofs() const
{
	return dofs_;
}

void
SymmetryInfo::set_dofs( std::map< Size, SymDof > const & dofs )
{
	dofs_ = dofs;
}

Size
SymmetryInfo::score_multiply( Size const res1, Size const res2 ) const
{
	assert ( res1 <= score_multiply_.size() && res2 <= score_multiply_.size() );
	assert ( res1 <= score_multiply_.size() && res2 <= score_multiply_.size() );

	//fpd  if one of the reses is a symm vrt
	//fpd     return score_multiply_[i] if the other is in scoring subunit
	//fpd     return 0 otherwise (since these terms can't be properly minimized)
	//fpd  this fixes problems with coordinate constraints
	if ( res1 > num_total_residues_without_pseudo() ) {
		return ( bb_is_independent(res2) ) ? score_multiply_[res2] : 0;
	} else if  (res2 > num_total_residues_without_pseudo() ) {
		return ( bb_is_independent(res1) ) ? score_multiply_[res1] : 0;
	} else {
		return score_multiply_[ ( bb_is_independent(res1) ) ? res2 : res1 ];
	}
}


Size
SymmetryInfo::interface_number( Size const res1, Size const res2 ) const
{
	return subunit_index(bb_is_independent(res1) ? res2 : res1);
}

void
SymmetryInfo::set_score_multiply_from_subunit_factors( utility::vector1< Size > const & score_multiply_vector_subunit, Size const nres_subunit, Size const n_subunits )
{
	score_multiply_.clear();
	for ( Size i = 1; i<= n_subunits; ++i ) {
		for ( Size j = 1; j<= nres_subunit; ++j ) {
			score_multiply_.push_back( score_multiply_vector_subunit[i] );
		}
	}
	for ( Size i = n_subunits + 1 ; i <= score_multiply_vector_subunit.size() ; ++i ) {
	    score_multiply_.push_back( score_multiply_vector_subunit[i] );
	}
}

void
SymmetryInfo::set_score_multiply( Size const res, Size const factor )
{
	assert ( res <= score_multiply_.size() );
	score_multiply_[ res ] = factor;
}


Size
SymmetryInfo::get_nres_subunit() const
{
	return nres_monomer_;
}

Size
SymmetryInfo::get_njumps_subunit() const
{
	return njump_monomer_;
}


bool
SymmetryInfo::get_use_symmetry() const
{
	return use_symmetry_;
}

bool
SymmetryInfo::cp_weighting_during_minimization() const
{
	return cp_weighting_during_minimization_;
}

void
SymmetryInfo::set_cp_weighting_during_minimization( bool setting )
{
  cp_weighting_during_minimization_ = setting;
}

SymSlideInfo
SymmetryInfo::get_slide_info() const
{
	return slide_info_;
}

void
SymmetryInfo::set_use_symmetry( bool setting )
{
	use_symmetry_ = setting;
}

/// @details  Returns set of DOF_IDs that follow a given one. Inefficient: it creates a list each time
SymmetryInfo::DOF_IDs
SymmetryInfo::dependent_dofs( DOF_ID const & id, Conformation const & conf ) const
{
	if ( !dof_is_independent( id, conf ) ) {
		utility_exit_with_message( "SymmetryInfo::dependent_dofs: dof is not independent!" );
	}

	Size const seqpos( id.rsd() );
	Size const atomno( id.atomno() );
	id::TorsionType const type( guess_torsion_type_of_dof_id( id, conf ) );
	Clones const & clones( type == id::JUMP ? jump_clones( conf.fold_tree().get_jump_that_builds_residue( id.rsd() ) ):
												 ( type == id::BB ? bb_clones( seqpos ) : chi_clones( seqpos ) ) );

	DOF_IDs dofs;
	for ( Clones::const_iterator pos= clones.begin(), epos=clones.end(); pos != epos; ++pos ) {
		if ( type == id::JUMP ) {
			dofs.push_back( DOF_ID( id::AtomID( atomno, conf.fold_tree().downstream_jump_residue( *pos )  ), id.type() ) );
		} else {
			dofs.push_back( DOF_ID( id::AtomID( atomno, *pos ), id.type() ) );
		}
	}
	return dofs;
}

///
bool
SymmetryInfo::dof_is_independent( DOF_ID const & id, Conformation const & conf ) const
{
	id::TorsionType const type( guess_torsion_type_of_dof_id( id, conf ) );

	switch ( type ) {
	case id::BB:
		return bb_is_independent( id.rsd() );
	case id::CHI:
		return chi_is_independent( id.rsd() );
	case id::JUMP:
		return jump_is_independent( conf.fold_tree().get_jump_that_builds_residue( id.rsd() ) );
	}

	utility_exit_with_message("dof_is_independent: unrecognized TorsionType!");
	return false;
}

// get a weight for derivative calculations
// weights are 1 for indep DOFs, 0 for dependent NON-JUMP DOFs
//    and may be any real for dependent jump dofs
core::Real
SymmetryInfo::get_dof_derivative_weight( DOF_ID const & id, Conformation const & conf ) const {
	id::TorsionType const type( guess_torsion_type_of_dof_id( id, conf ) );

	if ( type == id::BB ) {
		return bb_is_independent( id.rsd() ) ? 1. : 0.;
	} else if ( type == id::CHI ) {
		return chi_is_independent( id.rsd() ) ? 1. : 0.;
	} else if ( type == id::JUMP ) {
		int jumpnum = conf.fold_tree().get_jump_that_builds_residue( id.rsd() );
		std::map< Size, Real >::const_iterator it( jump_clone_wts_.find( jumpnum ) );
		return ( it == jump_clone_wts_.end() ? 1. : it->second );
	}

	utility_exit_with_message("get_dof_derivative_weight: unrecognized TorsionType!");
	return false;
}



///
bool
SymmetryInfo::torsion_is_independent( TorsionID const & id ) const
{
	return ( ( id.type() == id::BB   &&   bb_is_independent( id.rsd() ) ) ||
					 ( id.type() == id::CHI  &&  chi_is_independent( id.rsd() ) ) ||
					 ( id.type() == id::JUMP && jump_is_independent( id.rsd() ) ) );
}

bool
SymmetryInfo::atom_is_independent( AtomID const & id ) const
{
	return fa_is_independent( id.rsd() );
}

/// @details  Returns set of TorsionIDs that follow a given one. Inefficient: it creates a list each time
SymmetryInfo::TorsionIDs
SymmetryInfo::dependent_torsions( TorsionID const & id ) const
{
	if ( !torsion_is_independent( id ) ) {
		utility_exit_with_message( "SymmetryInfo::dependent_torsions: torsion is not independent!" );
	}

	Size const seqpos( id.rsd() );
	Clones const & seqpos_clones( id.type() == id::BB ? bb_clones( seqpos ) : chi_clones( seqpos ) );

//std::cerr << " dependent_torsions( TorsionID const & )  bb=" << (id.type() == id::BB) << "   chi=" << (id.type() == id::CHI) << "   jump=" << (id.type() == id::JUMP) << std::endl;

	TorsionIDs tors;
	for ( Clones::const_iterator pos= seqpos_clones.begin(), epos=seqpos_clones.end(); pos != epos; ++pos ) {
		tors.push_back( TorsionID( *pos, id.type(), id.torsion() ) );
	}
	return tors;
}

/// @details  Returns set of AtomIDs that follow a given one. Inefficient: it creates a list each time
SymmetryInfo::AtomIDs
SymmetryInfo::dependent_atoms( AtomID const & id ) const
{
  if ( !atom_is_independent( id ) ) {
    utility_exit_with_message( "SymmetryInfo::dependent_atoms: atom is not independent!" );
  }
	Size const seqpos( id.rsd() );
  Clones const & seqpos_clones( bb_clones( seqpos ) );

	AtomIDs atoms;
	for ( Clones::const_iterator pos= seqpos_clones.begin(), epos=seqpos_clones.end(); pos != epos; ++pos ) {
		atoms.push_back( AtomID( id.atomno(), *pos ) );
	}
	return atoms;
}

bool
SymmetryInfo::is_asymmetric_seqpos( Size const res ) const
{
	Size nres_monomer = num_independent_residues(), num_monomers = subunits() ;
	return ( !get_use_symmetry() || res > nres_monomer*num_monomers || res <= nres_monomer );
}

Size
SymmetryInfo::get_asymmetric_seqpos( Size const res ) const
{
	if( res > 0 && get_use_symmetry() ) {
		Size nres_monomer = num_independent_residues(), num_monomers = subunits() ;
		return ( res > nres_monomer*num_monomers ? res - nres_monomer*(num_monomers-1) : (res-1)%nres_monomer + 1 );
	} else {
		return res;
	}
}

Size
SymmetryInfo::subunit_index( Size const seqpos ) const {
	return ( (seqpos-1) / num_independent_residues() + 1 );
}

std::string
SymmetryInfo::get_jump_name(Size jnum) const {
	if( 0 == jnum2dofname_.count(jnum) ) utility_exit_with_message("bad jump num");
	return jnum2dofname_.find(jnum)->second;
}

Size
SymmetryInfo::get_jump_num(std::string jname) const {
	if( 0 == dofname2jnum_.count(jname) ) utility_exit_with_message("bad jump name");
	return dofname2jnum_.find(jname)->second;
}

void
SymmetryInfo::set_jump_name(Size jnum, std::string jname) {
	jnum2dofname_[jnum] = jname;
	dofname2jnum_[jname] = jnum;
}

Size
SymmetryInfo::num_slidablejumps() const {
	Size retval = 0;
	for(std::map<Size,SymDof>::const_iterator i = dofs_.begin(); i != dofs_.end(); i++) {
		if (i->second.allow_dof(1) || i->second.allow_dof(2) || i->second.allow_dof(3)) retval++;
	}
	return retval;
}

utility::vector1<char> const &
SymmetryInfo::get_components() const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	return components_;
}

std::map<char,std::pair<Size,Size> > const &
SymmetryInfo::get_component_bounds() const { 
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	return component_bounds_;
}

std::map<std::string,char> const &
SymmetryInfo::get_subunit_name_to_component() const { 
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	return name2component_;
}

std::map<std::string,utility::vector1<char> > const &
SymmetryInfo::get_jump_name_to_components() const { 
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	return jname2components_;
}

std::map<std::string,utility::vector1<Size> > const & 
SymmetryInfo::get_jump_name_to_subunits() const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	return jname2subunits_;
}
std::pair<Size,Size> const & 
SymmetryInfo::get_component_bounds(char c) const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	if( component_bounds_.find(c) == component_bounds_.end() ){
		utility_exit_with_message(std::string("no symmetry component ")+c);
	}
	return component_bounds_.find(c)->second;
}
Size
SymmetryInfo::get_component_lower_bound(char c) const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	if( component_bounds_.find(c) == component_bounds_.end() ){
		utility_exit_with_message(std::string("no symmetry component ")+c);
	}
	return component_bounds_.find(c)->second.first;
}
Size
SymmetryInfo::get_component_upper_bound(char c) const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	if( component_bounds_.find(c) == component_bounds_.end() ){
		utility_exit_with_message(std::string("no symmetry component ")+c);
	}
	return component_bounds_.find(c)->second.second;
}
char
SymmetryInfo::get_component_of_residue(Size ir) const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	if( ir > num_total_residues_without_pseudo() || ir < 1 ){
		utility_exit_with_message(std::string("no symmetry component for residue "));
	}
	Size irindep = (ir-1)%num_independent_residues()+1;
	for(std::map<char,std::pair<Size,Size> >::const_iterator i = component_bounds_.begin(); i != component_bounds_.end(); ++i){
		char component = i->first;
		Size lower = i->second.first;
		Size upper = i->second.second;
		// std::cerr << component << " " << lower << " " << upper << " " << irindep << std::endl;
		if( lower <= irindep && irindep <= upper ) return component;
	}
	utility_exit_with_message(std::string("no symmetry component for residue "));
	return ' ';
}
char
SymmetryInfo::get_subunit_name_to_component(std::string const & vname) const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	if( name2component_.find(vname) == name2component_.end() ){
		utility_exit_with_message(std::string("no symmetry component for ")+vname);
	}
	return name2component_.find(vname)->second;
}
utility::vector1<char> const & 
SymmetryInfo::get_jump_name_to_components(std::string const & jname) const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	if( jname2components_.find(jname) == jname2components_.end() ){
		utility_exit_with_message(std::string("no symmetry component for ")+jname);
	}
	return jname2components_.find(jname)->second;
}
utility::vector1<Size> const & 
SymmetryInfo::get_jump_name_to_subunits(std::string const & jname) const {
	if(components_.size()==0) utility_exit_with_message("function not for use in single component symmetry");
	if( jname2subunits_.find(jname) == jname2subunits_.end() ){
		utility_exit_with_message(std::string("no symmetry component for ")+jname);
	}
	return jname2subunits_.find(jname)->second;
}


void 
SymmetryInfo::set_multicomponent_info(
	utility::vector1<char> const & components,
	std::map<char,std::pair<Size,Size> > const & component_bounds,
	std::map<std::string,char> const & name2component,
	std::map<std::string,utility::vector1<char> > const & jname2component,
		std::map<std::string,utility::vector1<Size> > const & jname2subunits
){
	components_ = components;
	component_bounds_ = component_bounds;
	name2component_ = name2component;
	jname2components_ = jname2component;
	jname2subunits_ = jname2subunits;
}


} // symmetry
} // conformation
} // core
