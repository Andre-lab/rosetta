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

#include <core/types.hh>
#include <core/id/SequenceMapping.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

#include <algorithm>

#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>


namespace core {
namespace id {

using namespace ObjexxFCL;

/// @brief ctor
SequenceMapping::SequenceMapping( utility::vector1< Size > const & mapping ):
	size2_( mapping.size() ),
	mapping_( mapping )
{}

SequenceMapping::SequenceMapping(
	conformation::signals::LengthEvent const & event
)
	: size2_( 0 )
{

	//int direction(0);
	bool longer( event.length_change > 0 );
	Size upstream_res(0), old_num_res(0);

	mapping_.clear();

	if ( event.tag == conformation::signals::LengthEvent::INVALIDATE ) {
		//utility_exit_with_message("invalidated LengthEvent passed");
		//i guess in this case it's better to return an empty sequence mapping and let
		//the observers deal with it
		return;
	}

	else if ( event.tag == conformation::signals::LengthEvent::RESIDUE_APPEND ) {
		//direction = 1;
		upstream_res = event.position;
	}
	else if( event.tag == conformation::signals::LengthEvent::RESIDUE_PREPEND ) {
		//direction = 1;
		upstream_res = event.position - 1;
	}
	else if( event.tag == conformation::signals::LengthEvent::RESIDUE_DELETE ) {
		//direction = -1;
		upstream_res = event.position - 1;
	}
	else {
		utility_exit_with_message(
			"unknown signal triggered by conformation length change. please update this file"
		);
	}

	old_num_res = event.conformation->size() - event.length_change;

	for ( Size i = 1; i <= upstream_res; ++i ) mapping_.push_back( i );

	//if ( direction == 1 ) mapping_.push_back( upstream_res + 1 + event.length_change );
	if ( longer ) {
		if (upstream_res < old_num_res)
			mapping_.push_back( upstream_res + 1 + event.length_change );
	} else {
		for( int i = event.length_change; i < 0; ++i ) mapping_.push_back( 0 );
	}
	//for ( Size i = upstream_res + 2; i <= old_num_res; ++i ) mapping_.push_back( i + direction );
	Size downstream_res( mapping_.size() + 1 );
	runtime_assert( downstream_res <= old_num_res + 1);
	for( Size i = downstream_res; i <= old_num_res; ++i ) mapping_.push_back( i + event.length_change );
}

SequenceMapping::~SequenceMapping() {}


SequenceMapping::SequenceMapping( SequenceMapping const & src )
	: ReferenceCount(src)
{
	*this = src;
}

SequenceMapping &
SequenceMapping::operator = ( SequenceMapping const & src ) {
	size2_      = src.size2();
	mapping_    = src.mapping();

	return *this;
}

void
SequenceMapping::resize( Size const s1, Size const s2 )
{
	size2_ = s2;
	mapping_.clear();
	mapping_.resize(s1,0);
}

///
void
SequenceMapping::reverse()
{
	// update size2!
	size2_ = *max_element( mapping_.begin(), mapping_.end() );

	utility::vector1< Size > new_mapping( size2_, 0 );
	for ( Size i = 1; i <= size1(); ++i ) {
		if ( mapping_[i] != 0 ) {
			new_mapping[ mapping_[i] ] = i;
		}
	}

	size2_ = mapping_.size();
	mapping_.swap( new_mapping );
}

/// access
Size
SequenceMapping::size1() const
{
	return mapping_.size();
}

Size
SequenceMapping::size2() const
{
	return size2_;
}

//
bool
SequenceMapping::all_aligned() const
{
	bool aligned( true );
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[i] == 0 ) {
			aligned = false;
			break;
		}
	}
	return aligned;
}

//
bool
SequenceMapping::is_identity() const {
	bool identity( true );
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[i] != i ) {
			identity = false;
			break;
		}
	}
	return identity;
}

bool
SequenceMapping::is_identity_ignore_gaps() const {
	bool identity( true );
	for ( Size i = 1; i <= size1(); ++i ) {
		if ( mapping_[i] != 0 && mapping_[i] != i ) {
			identity = false;
			break;
		}
	}
	return identity;

}

void
SequenceMapping::size2( Size const s2 )
{
	size2_ = s2;
}

void
SequenceMapping::push_back( Size const al )
{
	mapping_.push_back( al );
}

void
SequenceMapping::delete_source_residue( Size const pos1 )
{
	mapping_.erase( mapping_.begin() + pos1-1 );
}

void
SequenceMapping::show() const
{
	show( basic::T("id.SequenceMapping") );
}

void
SequenceMapping::show( std::ostream & output ) const {
   for ( Size i=1; i<= size1(); ++i ) {
      output << ("id.SequenceMapping ") << i << " --> ";
      if ( mapping_[i] ) output << mapping_[i] << std::endl;
      else output << "----" << std::endl;
   }
}


void
SequenceMapping::set_offset( int setting ) {
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[ i ] ) {
			int pos =  static_cast< int >(mapping_[ i ]) - setting;
			mapping_[ i ] = ( pos > 0 ) ? pos : 0;
		}
	}
}


void
SequenceMapping::insert_source_residue( Size const pos1 )
{
	mapping_.insert( mapping_.begin() + pos1-1, 0 );
}

void
SequenceMapping::insert_aligned_residue( Size const pos1, Size const pos2 )
{
	mapping_.insert( mapping_.begin() + pos1-1, pos2 );
}

// same as insert_aligned_residue, but a couple of extra checks on size1 and size2.
void
SequenceMapping::insert_aligned_residue_safe(
	Size const pos1,
	Size const pos2
) {
	if ( pos1 == 0 ) return;
	size2( std::max( pos2, size2() ) );
	if ( pos1 >= size1() ) {
		mapping_.resize( pos1, 0 );
	}

	mapping_[ pos1 ] = pos2;
}

void
SequenceMapping::insert_target_residue( Size const pos )
{
	++size2_;
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[i] >= pos ) ++mapping_[i];
	}
}

void
SequenceMapping::delete_target_residue( Size const pos )
{
	--size2_;
	for ( Size i=1; i<= size1(); ++i ) {
		if ( mapping_[i] == pos ) mapping_[i] = 0;
		else if ( mapping_[i] > pos ) --mapping_[i];
	}
}

void
SequenceMapping::clear()
{
	mapping_.clear();
	size2_ = 0;
}

Size
SequenceMapping::operator[]( Size const pos1 ) const {
	if ( pos1 > mapping_.size() ) return 0;
	return mapping_[ pos1 ];
}

Size &
SequenceMapping::operator[]( Size const pos1 ) {
	if ( pos1 > mapping_.size() ) mapping_.resize( pos1, 0 );
	return mapping_[ pos1 ];
}

bool SequenceMapping::operator==( SequenceMapping const & other ) const {
	if ( other.size1() != size1() ) return false;

	for ( Size ii = 1; ii <= size1(); ++ii ) {
		if ( mapping_[ii] != other[ii] ) return false;
	}

	return true;
}

std::string SequenceMapping::to_string() const {
	std::string retval( "" );
	for ( Size i = 1; i <= size1(); ++i ) {
		retval += string_of(i) + " --> ";
		if ( mapping_[i] ) {
			retval += string_of(mapping_[i]) + "\n";
		} else {
			retval += static_cast< std::string > ("----") + "\n";
		}
	}
	return retval;
}

/// @details combine all input sequence mappings into one.
/// sequentially, that is
core::id::SequenceMappingOP
combine_sequence_mappings(
	utility::vector1< core::id::SequenceMapping > const & smaps
){

	using namespace core::id;

	//gigo :)
	if( smaps.size() == 0 ) return new SequenceMapping() ;

	SequenceMappingOP composite_smap = new SequenceMapping();
	*composite_smap = smaps[1];

	for( core::Size i = 2; i <= smaps.size(); ++i ){
		combine_sequence_mappings( *composite_smap, smaps[i] );
	}

	return composite_smap;

} //combine_sequence_mappings



/// @details combine smap_to_add into smap,
/// i.e. smap[j] becomes smap_to_add[ smap[j] ]
void
combine_sequence_mappings(
	core::id::SequenceMapping & smap,
	core::id::SequenceMapping const & smap_to_add )
{

	for( core::Size i = 1; i <= smap.size1(); ++i){

		if( smap[i] != 0 ){

			if( smap[i] <= smap_to_add.size1() ) smap[i] = smap_to_add[ smap[i] ];

			else smap[i] = 0;
		}
	}
}

} // id
} // core

