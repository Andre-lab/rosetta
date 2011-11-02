// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragSet.cc
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson (tex@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///

// Unit Headers
#include <core/fragment/OrderedFragSet.hh>

// Package Headers
//#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/OrderedFragSetIterator_.hh>
#include <core/fragment/Frame.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh>

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/types.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
#include <ostream>

#include <core/fragment/FrameIterator.hh>
#include <utility/vector1.hh>

namespace core {
namespace fragment {

using namespace kinematics;

static basic::Tracer tr("core.fragments");
// preliminary reader method --- reads classic rosetta++ frag files


OrderedFragSet::OrderedFragSet() {}
OrderedFragSet::~OrderedFragSet() {}

FragSetOP OrderedFragSet::clone() const {
	return new OrderedFragSet( *this );
}
FragSetOP OrderedFragSet::empty_clone() const {
	return new OrderedFragSet();
}

///@brief get fragments that start somewhere between start and end
Size OrderedFragSet::region(
  MoveMap const&,
  core::Size start,
  core::Size end, //not used
  core::Size, //min_overlap not used
  core::Size, //min_length not used
  FrameList &frame_list
) const {
  Size count( 0 );
  for ( Size pos=start; pos<=end; pos++ ) {
    count += frames( pos, frame_list );
  }
  return count;
}


/// @brief Accessor for the Frame at the specified insertion position. Returns false if
/// there is no frame at the specified position.
Size OrderedFragSet::frames( Size pos, FrameList &out_frames ) const
{
	FrameMap::const_iterator it = frames_.find(pos);
	if ( it == frames_.end() ) return 0;
	if ( it->second.begin() != it->second.end() ) {
		copy( it->second.begin(), it->second.end(), back_inserter( out_frames ) ); // should append frames
		return it->second.size();
	}

	return 0;
}

FrameIterator OrderedFragSet::begin() const {
 return FrameIterator( new OrderedFragSetIterator_( frames_.begin(), frames_.end() ) );
}

FrameIterator OrderedFragSet::end() const {
 return FrameIterator( new OrderedFragSetIterator_( frames_.end(), frames_.end() ) );
}

bool OrderedFragSet::empty() const {
	return frames_.size()==0;
}



void OrderedFragSet::add_( FrameOP aframe )
{
  Size seqpos( aframe->start() );
  frames_[ seqpos ].push_back( aframe );
}

}//fragment
}// core
