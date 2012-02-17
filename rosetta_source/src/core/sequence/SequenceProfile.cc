// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SequenceProfile.cc
/// @brief class definition for a given scoring scheme for an alignment.
/// @detailed Simply based on comparing single characters from two protein
/// sequences, along with affine gap penalties of the form penalty = A + Bk, where
/// A represents the penalty for starting a gap, and B represents the penalty for
/// extending a previously opened gap by k characters.
/// @author James Thompson

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/sequence/MatrixScoringScheme.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/pointer/owning_ptr.hh>

#include <core/chemical/AA.hh>

// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <iostream>
#include <string>

#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>


namespace core {
namespace sequence {

static basic::Tracer tr( "core.sequence.SequenceProfile" );

void SequenceProfile::read_from_checkpoint(
		utility::file::FileName const & fn,
		bool negative_better
) {
	utility::io::izstream input( fn );
	negative_better_ = negative_better;
	// order of amino acids in the .checkpoint file
	static utility::vector1< core::chemical::AA > order;
	order.resize( 20 );
	order[ 1] = core::chemical::aa_from_oneletter_code( 'A' );
	order[ 2] = core::chemical::aa_from_oneletter_code( 'C' );
	order[ 3] = core::chemical::aa_from_oneletter_code( 'D' );
	order[ 4] = core::chemical::aa_from_oneletter_code( 'E' );
	order[ 5] = core::chemical::aa_from_oneletter_code( 'F' );
	order[ 6] = core::chemical::aa_from_oneletter_code( 'G' );
	order[ 7] = core::chemical::aa_from_oneletter_code( 'H' );
	order[ 8] = core::chemical::aa_from_oneletter_code( 'I' );
	order[ 9] = core::chemical::aa_from_oneletter_code( 'K' );
	order[10] = core::chemical::aa_from_oneletter_code( 'L' );
	order[11] = core::chemical::aa_from_oneletter_code( 'M' );
	order[12] = core::chemical::aa_from_oneletter_code( 'N' );
	order[13] = core::chemical::aa_from_oneletter_code( 'P' );
	order[14] = core::chemical::aa_from_oneletter_code( 'Q' );
	order[15] = core::chemical::aa_from_oneletter_code( 'R' );
	order[16] = core::chemical::aa_from_oneletter_code( 'S' );
	order[17] = core::chemical::aa_from_oneletter_code( 'T' );
	order[18] = core::chemical::aa_from_oneletter_code( 'V' );
	order[19] = core::chemical::aa_from_oneletter_code( 'W' );
	order[20] = core::chemical::aa_from_oneletter_code( 'Y' );

	std::string aa_seq;
	utility::vector1< utility::vector1< core::Real > > new_prof;
	// profile is indexed by the order of amino acids in core::chemical::AA

	if ( !input ) {
		utility_exit_with_message( "ERROR: Unable to open file!" );
	}
	std::string line;
	getline( input, line ); // SKIPS FIRST LINE
	while( getline( input, line ) ) {
		if ( line.substr(0,3) == "END" ) break; // end of profile
		std::string aa;
		utility::vector1< core::Real > prof_row;
		prof_row.resize( order.size() );

	 	std::istringstream ls( line );
		ls >> aa;

		core::Real aa_prob;
		ls >> aa_prob;
		Size index(1);
		while ( !ls.fail() ) {
			//prof_row.push_back( aa_prob );
			prof_row[ order[index] ] = aa_prob;

			ls >> aa_prob;
			++index;
		}

		aa_seq += aa;
		new_prof.push_back( prof_row );
	}

	sequence( aa_seq );
	profile( new_prof );
}

void SequenceProfile::read_from_file(
	utility::file::FileName const & fn
) {
	negative_better_ = false;
	// order of amino acids in the .pssm file
	static utility::vector1< core::chemical::AA > order;
	order.resize( 20 );
	order[ 1] = core::chemical::aa_from_oneletter_code( 'A' );
	order[ 2] = core::chemical::aa_from_oneletter_code( 'R' );
	order[ 3] = core::chemical::aa_from_oneletter_code( 'N' );
	order[ 4] = core::chemical::aa_from_oneletter_code( 'D' );
	order[ 5] = core::chemical::aa_from_oneletter_code( 'C' );
	order[ 6] = core::chemical::aa_from_oneletter_code( 'Q' );
	order[ 7] = core::chemical::aa_from_oneletter_code( 'E' );
	order[ 8] = core::chemical::aa_from_oneletter_code( 'G' );
	order[ 9] = core::chemical::aa_from_oneletter_code( 'H' );
	order[10] = core::chemical::aa_from_oneletter_code( 'I' );
	order[11] = core::chemical::aa_from_oneletter_code( 'L' );
	order[12] = core::chemical::aa_from_oneletter_code( 'K' );
	order[13] = core::chemical::aa_from_oneletter_code( 'M' );
	order[14] = core::chemical::aa_from_oneletter_code( 'F' );
	order[15] = core::chemical::aa_from_oneletter_code( 'P' );
	order[16] = core::chemical::aa_from_oneletter_code( 'S' );
	order[17] = core::chemical::aa_from_oneletter_code( 'T' );
	order[18] = core::chemical::aa_from_oneletter_code( 'W' );
	order[19] = core::chemical::aa_from_oneletter_code( 'Y' );
	order[20] = core::chemical::aa_from_oneletter_code( 'V' );

	utility::io::izstream input( fn );

	if ( !input ) {
		std::string msg(
			"ERROR: Unable to open file " +
			static_cast< std::string > (fn) +
			"!"
		);
		utility_exit_with_message( msg );
	}
	std::string line;

	tr.Debug << "reading from " << fn << std::endl;

	// read in two header lines
	getline( input, line );
	getline( input, line );

	// initialize headers
	getline( input, line );
	std::istringstream line_stream( line );
	while ( !line_stream.fail() ) {
		std::string aa;
		line_stream >> aa;
		if ( line_stream.fail() ) continue;

		alphabet_.push_back( aa );
		if ( alphabet_.size() >= 20 ) break; // super-hack for pssm file format! bad james, bad!
	}

	std::string seq;
	while( getline( input, line ) ) {
		std::istringstream line_stream( line );
		//std::cout << "line = " << line << std::endl;
		Size pos;
		string aa;

		line_stream >> pos >> aa;
		if ( line_stream.fail() ) continue;

		utility::vector1< Real > prof_row;
		prof_row.resize( order.size() );

		Real score;
		line_stream >> score;
		Size index(1);
		while ( !line_stream.fail() && index <= order.size() ) {
			prof_row[ order[index] ] = score;
			line_stream >> score;
			++index;
		}

		profile_.push_back( prof_row );

		seq += aa;
	} // while( getline( input, line ) )

	tr.Debug << "Read sequence " << seq << " from  " << fn << std::endl;
	tr.Debug << "profile dimensions are " << profile_.size() << "x"
		<< profile_.front().size() << "." << std::endl;
	sequence( seq );
	id( std::string(fn) );
} // read_from_file

void SequenceProfile::generate_from_sequence( Sequence const & seq, std::string matrix ) {
	MatrixScoringScheme mat;
	mat.read_from_database(matrix);

	utility::vector1< utility::vector1< core::Real > > new_prof;
	for (core::Size ii(1), end(seq.length()); ii <= end; ++ii) {
		utility::vector1< Real > prof_row( mat.values_for_aa( seq[ii] ) );
		prof_row.resize( core::chemical::num_canonical_aas, -10000 );
		new_prof.push_back(prof_row);
	}

	sequence( seq.sequence() );
	profile( new_prof );
	negative_better_ = false;
}

/// @brief Returns the 2D vector1 of Real-values representing this profile.
utility::vector1< utility::vector1< Real > > const &
SequenceProfile::profile() const {
	return profile_;
}

/// @brief Multiply all profile weights by factor
void SequenceProfile::rescale(core::Real factor) {
	for ( Size ii = 1; ii <= profile().size(); ++ii ) {
		for ( utility::vector1< core::Real >::iterator
						it = profile_[ii].begin(), end = profile_[ii].end(); it != end; ++it ) {
			*it *= factor;
		}
	}
	if ( factor < 0 ) {
		tr << "Flipping sense of negative_better." << std::endl;
		negative_better_ = ! negative_better_;
	}
}

void SequenceProfile::convert_profile_to_probs( core::Real temp /*= 1.0*/ ) {
	temp_ = temp;
	utility::vector1< utility::vector1< Real > > new_prof;
	for ( Size ii = 1; ii <= profile().size(); ++ii ) {
		utility::vector1< core::Real > new_prof_row = prof_row(ii);
		scores_to_probs_( new_prof_row, temp, negative_better_ );
		new_prof.push_back( new_prof_row );
	}
	profile( new_prof );
	negative_better_ = false;
}

void SequenceProfile::global_auto_rescale() {
	core::Real maxval(0.0);
	for ( Size ii = 1; ii <= profile_.size(); ++ii ) {
		for ( Size jj = 1; jj <= profile_[ii].size(); ++jj ) {
			if ( maxval < std::abs(profile_[ii][jj]) ) {
				maxval = std::abs(profile_[ii][jj]);
			}
		}
	}
	rescale( 1/maxval );
}


void SequenceProfile::profile(
	utility::vector1< utility::vector1< core::Real > > const & new_prof
) {
	profile_ = new_prof;
}

void SequenceProfile::insert_char(
	core::Size pos,
	char new_char
) {
	using core::Real;
	using utility::vector1;
	// add in a profile column of zeroes.
	vector1< Real > zero_col( width(), 0.0 );
	vector1< vector1< Real > > new_prof, old_prof( profile() );
	for ( Size i = 0; i <= length() + 1; ++i ) {
		if ( i == pos )                new_prof.push_back( zero_col );
		if ( i >= 1 && i <= length() ) new_prof.push_back( old_prof[i] );
	}

	profile( new_prof );
	Sequence::insert_char( pos, new_char );
}

void SequenceProfile::delete_position(
	core::Size pos
) {
	//std::cout << "sequence is " << sequence() << std::endl;
	using core::Real;
	using utility::vector1;

	runtime_assert( pos <= length() );

	vector1< vector1< Real > > new_prof( profile() );
	vector1< vector1< Real > >::iterator it
		= new_prof.begin() + pos - start();

	runtime_assert( it != new_prof.end() );

	new_prof.erase( it );
	profile( new_prof );

	Sequence::delete_position( pos );
	//std::cout << "sequence is " << sequence() << std::endl;
}

/// @brief Returns the number of distinct values at each position in this profile.
Size SequenceProfile::width() const {
	assert( check_internals_() );
	return profile_[1].size();
}

/// @brief Returns the vector1 of values at this position.
utility::vector1< Real > const &
SequenceProfile::prof_row( Size pos ) const {
	runtime_assert( pos <= profile_.size() );
	return profile_[ pos ];
}

void SequenceProfile::scores_to_probs_(
	utility::vector1< core::Real > & scores,
	core::Real kT,
	bool negative_better /* = false */
) const {
	using std::exp;
	using utility::vector1;

	// calculate partition (aka Z), with this definition:
	// Z = sum( exp( score / kT ) )
	core::Real partition( 0.0 );
	for ( vector1< core::Real >::iterator
				it = scores.begin(), end = scores.end(); it != end; ++it
	) {
		if ( negative_better ) {
			*it = exp( -1 * *it / kT );
		} else {
			*it = exp( *it / kT );
		}
		partition += *it;
	}

	// transform scores using the partition calculated above:
	// P(s) = exp( -1  * score / kT ) ) / Z
	for ( utility::vector1< core::Real >::iterator
				it = scores.begin(), end = scores.end(); it != end; ++it
	) {
		*it = *it / partition;
	}
} // scores_to_probs_

std::ostream & operator<<( std::ostream & out, const SequenceProfile & p ) {
	Size width = 8;
	Size precision = 3;

	out << p.to_string() << std::endl;
	for ( Size i = 1; i <= p.length(); ++i ) {
		for ( Size j = 1; j <= p.width(); ++j ) {
			out << ObjexxFCL::fmt::F( width, precision, p.prof_row(i)[j] );
		}
		out << std::endl;
	}

	return out;
}

bool SequenceProfile::check_internals_() const {
	using core::Real;
	using utility::vector1;

	runtime_assert( profile_.size() == length() );

	for ( Size i = 1; i <= profile_.size(); ++i ) {
		runtime_assert( profile_[i].size() == alphabet_.size() );
	}

	return true;
}

} // sequence
} // core
