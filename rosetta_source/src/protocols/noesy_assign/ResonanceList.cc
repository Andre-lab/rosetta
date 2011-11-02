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
/// @detailed
///
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/Resonance.hh>
#include <protocols/noesy_assign/ResonanceList.hh>



// Package Headers
#include <protocols/noesy_assign/Exceptions.hh>
#include <core/id/NamedAtomID.hh>
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Project Headers
#include <core/chemical/AA.hh>

// Utility headers
#include <ObjexxFCL/format.hh>

// #include <utility/exit.hh>
// #include <utility/excn/Exceptions.hh>
// #include <utility/vector1.fwd.hh>
// #include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <deque>

#include <utility/vector1.hh>



static basic::Tracer tr("protocols.noesy_assign.resonances");

using core::Real;
using namespace core;
using namespace basic;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

namespace protocols {
namespace noesy_assign {


ResonanceList::ResonanceList( std::string const& sequence ) : sequence_( sequence )
{}

/// translate sequence information into AA
core::chemical::AA ResonanceList::aa_from_resid( core::Size resi ) const {
  runtime_assert(  resi <= sequence_.size() );
  return core::chemical::aa_from_oneletter_code( sequence_[ resi-1 ] );
}


/* helper function for 'read_from_stream'
	 each element is first put into a DEQUE such that atoms like methyl-protons can be combined if they are assigned the same frequency.
	 -HB1, HB2 -- > QB

 example:
	 input:
	  HB1  1.00
		HB2  1.00
		CD  13.00

	 output:
	  QB   1.00
		CD  13.00
*/
void process_last_resonances( std::deque< Resonance >& last_resonances, bool drain=false ) {

  std::string combine_name;
  Resonance first_res = last_resonances.front();
  last_resonances.pop_front();

	//is this a HB1, HG1 or a 1HG2 type of name ?
  bool single_char ( first_res.name().size() == 3 );

	//do we have enough protons for a QQX type of combined atom. (2 methyl groups... )
  bool double_methyl( last_resonances.size() == 6 );

	//
  Size str_cmp_length( 1 );

  std::string pseudo_letter( "Q" ); //default, single methyl group, proton
  if ( first_res.name().substr(0,1)=="C" ) pseudo_letter = "C"; //not default, we don't have a proton
	if ( first_res.name().substr(0,1)=="Q" ) pseudo_letter = "QQ"; //not default, this is a double-methyl...
	//if the two protons stuffed into the ambiguity queue are, .e.g, QG1 + QG2 --> QQG

	// construct combined atom name, QX, QQX
	if ( single_char ) { 	//things like HB1+HB2 or QG1+QG2
    combine_name=pseudo_letter+first_res.name().substr(1,1);
  } else if ( double_methyl ) { //-->QQX
    combine_name=pseudo_letter+pseudo_letter+first_res.name().substr(1,1);
  } else { // HG11+HG12+HG13 --> QG1
    combine_name=pseudo_letter+first_res.name().substr(1,2);
    str_cmp_length=2;
  }

	//now figure out, how many atoms can be combined...
	Size limit( drain ? 0 : 1 ); //should we go to the very end, or leave last atom behind...
  while( last_resonances.size() > limit ) {//check others
    if ( first_res.name().substr( 1, str_cmp_length ) == last_resonances.front().name().substr( 1, str_cmp_length ) ) {
			last_resonances.pop_front();
		} else { //could not combine...
			combine_name = first_res.name();
			break;
		}
  } //now replace front of deque with the combined atom
  last_resonances.push_front( Resonance( first_res.label(), first_res.freq(), first_res.error(), core::id::NamedAtomID( combine_name, first_res.resid() ) ) );
}

void ResonanceList::read_from_stream( std::istream& is ) {
  using namespace core::chemical; //for AA
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
  std::string line;
  std::deque< Resonance > last_resonances;
  while( getline( is, line ) ) {
    core::Size label;
    core::Real freq;
    core::Real error;
    std::string name;
    core::Size resn;
    AA aa;
    std::istringstream line_stream( line );

		//read from stream...
		tr.Trace << "read line: " << line << std::endl;
    line_stream >> label >> freq >> error >> name;
		if ( params.ignore_resonancefile_tolerances_ ) error = 0.0;

		//check for stream-error
		if ( !line_stream.good() ) {
			tr.Info << "ignore line : " << line << std::endl;
			continue; //ignore weird lines
		}

		if ( name=="" ) break;//got an empty line at end of file ?

		//read optional fields
		line_stream >> resn;
    std::string aa_name;
    line_stream >> aa_name;

		//process optional fields...
    if ( line_stream.good() ) { //optional field present?
      if ( aa_name.size() == 1 ) {
				aa = aa_from_oneletter_code( aa_name[ 0 ] );
      } else if ( aa_name.size() == 3 ) {
				aa = aa_from_name( aa_name );
      } else {
				throw utility::excn::EXCN_BadInput( "did not recognize aminoacid: " + aa_name);
      }
      if ( sequence_.size() < resn ) {
				while ( sequence_.size() < resn-1 ) sequence_.push_back('X');
				sequence_[ resn-1 ]=oneletter_code_from_aa( aa );
      } else if ( sequence_[ resn-1 ]!=oneletter_code_from_aa( aa ) ) {
				tr.Warning << "sequence mismatch!! check your data: found " << name_from_aa( aa ) << " in line " << line
									 << " which does not match " << sequence_[ resn-1 ] << "\n" << sequence_ << std::endl;
      }
    } else { //optional field was not present... use the pre-supplied fasta-sequence
		  aa = aa_from_resid( resn );
		}

		// replace some names...
		if ( name == "HN" ) name ="H";
		if ( aa == aa_leu && name == "HD1" ) name = "QD1";
		if ( aa == aa_leu && name == "HD2" ) name = "QD2";
		if ( aa == aa_val && name == "HG1" ) name = "QG1";
		if ( aa == aa_val && name == "HG2" ) name = "QG2";
		if ( aa == aa_thr && name == "HG2" ) name = "QG2";
		if ( aa == aa_ala && name == "HB"  ) name = "QB";
		if ( aa == aa_gly && name == "HA"  ) name = "QA";
		if ( aa == aa_ile && name == "HG1" ) name = "QG1";
		if ( aa == aa_ile && name == "HG2" ) name = "QG2";
		if ( aa == aa_ile && name == "HD1" ) name = "QD1";

		// before assigning to Resonance List put it in DEQUE (push_back), maybe it will get combined with next resonance...
    last_resonances.push_back( Resonance( label, freq, error, core::id::NamedAtomID( name, resn ) ) );
    if ( freq != last_resonances.front().freq() ) { ///if we have just read a new frequency we need to finalize what is in the DEQUE
      if ( last_resonances.size() > 2 ) process_last_resonances( last_resonances ); //just 2 --> 1 old freq and 1 new freq
      map_[ last_resonances.front().label() ] = last_resonances.front(); ///assign most forward value in our DEQUE
      last_resonances.pop_front(); ///remove the just-assigned value
    }
  }

	//still unprocessed data in DEQUE ?
  if ( last_resonances.size() > 1 ) process_last_resonances( last_resonances, true /*drain*/ );
  map_[ last_resonances.front().label() ]= last_resonances.front();

	//ASSERT: we have captured everything...
  runtime_assert( last_resonances.size() == 1 );

	///no problems ?
  if ( is.fail() && is.eof() && is.bad() ) {
    tr.Error << "[ERROR WHILE READING]" << std::endl;
  }

	//post-process input...
	update_residue_map();
}

///@brief write ResonanceList to stream
void ResonanceList::write_to_stream( std::ostream& os   ) const {
  for ( ResonanceIDs::const_iterator it = map_.begin(); it != map_.end(); ++it ) {
    runtime_assert( it->first == it->second.label() );
    it->second.write_to_stream( os );
    if ( sequence_.size() ) {
      using namespace core::chemical; //for AA
			AA aa( aa_from_resid( it->second.resid() ) );
      os << " " << name_from_aa( aa ) << " " << oneletter_code_from_aa( aa );
    }
    os << std::endl;
  }
}

///@brief write ResonanceList in TALOS format
void ResonanceList::write_talos_format( std::ostream& os, bool backbone_only ) const {
	if ( sequence_.size() == 0 ) {
		utility_exit_with_message( "sequence information required to write TALOS format -- use -in:file:fasta" );
	}

	//write sequence section
	Size const TALOS_SEQ_LINE_LEN( 50 );
	Size const TALOS_SEQ_BLCK_SIZE( 10 );
	for ( Size ct=0; ct<sequence_.size(); ct++ ) {
		if ( ct % TALOS_SEQ_LINE_LEN == 0 && ct > 1 ) os <<"\n";
		if ( ct % TALOS_SEQ_LINE_LEN == 0 ) os << "DATA SEQUENCE";
		if ( ct % TALOS_SEQ_BLCK_SIZE == 0 ) os << " ";
		os << sequence_[ ct ];
	}

	///write format section
	os << "\n\n";
	os << "VARS RESID RESNAME ATOMNAME SHIFT\n";
	os << "FORMAT %4d %1s %4s %8.3f\n";
	os << "\n";


	///write resonances
  for ( ResonanceIDs::const_iterator it = map_.begin(); it != map_.end(); ++it ) {
    runtime_assert( it->first == it->second.label() );
		if ( sequence_.size() < it->second.resid() ) {
			tr.Error << " no sequence information for residue " << it->second.resid() << std::endl;
			utility_exit_with_message( "sequence information required for all residues to write TALOS format -- use -in:file:fasta" );
		}
		using namespace core::chemical; //for AA
		std::string atom( it->second.name() );

		//are we backbone
		bool const is_backbone( //backbone im Sinne von SPARTA
													 atom=="H" || atom=="HN" ||
													 atom=="CA" || atom=="CB" ||
													 atom=="N" ||
													 atom=="1HA" || atom=="2HA" || atom=="3HA" );
		//write individual atom
		if ( is_backbone || !backbone_only ) {
			if ( atom=="1HA" ) atom="HA3";
			if ( atom=="2HA" ) atom="HA2";
			if ( atom=="3HA" ) atom="HA1";
			AA aa( aa_from_resid( it->second.resid() ) );
			os << ObjexxFCL::fmt::RJ( 5, it->second.resid() ) << " ";
			os << oneletter_code_from_aa( aa ) << " ";
			os << ObjexxFCL::fmt::RJ( 3, atom=="H" ? "HN" : atom ) << " ";
			os << ObjexxFCL::fmt::F( 7, 3, it->second.freq() ) << std::endl;
		}
	}
	os << std::endl;
}


///retrieve Resonance --- throws EXCN_UnknonwResonance if atom not found
Resonance const& ResonanceList::operator[] ( core::id::NamedAtomID const& atom ) const {
	ResidueMap::const_iterator it_res( by_resid_.find( atom.rsd() ) );
	if ( it_res != by_resid_.end() ) {
		Resonances const& reso_list( it_res->second );
		for ( Resonances::const_iterator it = reso_list.begin(); it != reso_list.end(); ++it ) {
			if ( it->atom() == atom ) return *it;
		}
	}
	throw EXCN_UnknownResonance( atom, "can't find atom ");
  return map_.begin()->second; //to make compiler happy
}

///retrieve Resonance ---  throws EXCN_UnknonwResonance if atom not found
Resonance const& ResonanceList::operator[] ( core::Size key ) const {
  ResonanceIDs::const_iterator iter = map_.find( key );
  if ( iter == map_.end() ) {
    throw EXCN_UnknownResonance( id::BOGUS_NAMED_ATOM_ID, "can't find resonance " + ObjexxFCL::string_of( key ) );
  }
  return iter->second;
}

///create map with all resonances sorted by residue number
void ResonanceList::update_residue_map() {
	by_resid_.clear();
	for ( ResonanceIDs::const_iterator it = map_.begin(); it != map_.end(); ++it ) {
    runtime_assert( it->first == it->second.label() );
		by_resid_[ it->second.resid() ].push_back( it->second );
	}
}

///retrieve list of Resonance at certain residue ---  throws EXCN_UnknonwResonance if residue number not found
ResonanceList::Resonances const& ResonanceList::resonances_at_residue( core::Size resid ) const {
	ResidueMap::const_iterator it_res( by_resid_.find( resid ) );
	if ( it_res != by_resid_.end() ) {
		return it_res->second;
	}
	throw EXCN_UnknownResonance( id::BOGUS_NAMED_ATOM_ID, "can't find resonance with residue " + ObjexxFCL::string_of( resid ) );
	return by_resid_.begin()->second; //to make compile happy
}


bool ResonanceList::has_residue( core::Size resid ) const {
	return by_resid_.find( resid ) != by_resid_.end();
}

}
}
