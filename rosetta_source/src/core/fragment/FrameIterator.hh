// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FrameIterator.hh
/// @brief
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date   Wed Aug 22 12:08:31 2007
///

#ifndef core_fragments_FrameIterators_HH
#define core_fragments_FrameIterators_HH

// Unit Headers
#include <core/fragment/FrameIterator.fwd.hh>

// AUTO-REMOVED #include <core/fragment/FragSet.fwd.hh>

// Package Headers
// no fwd here. The FrameIteratorWorker_ is really a HIDDEN detail. since this isn't a NODE header there is no danger of polluting the
// header graph.
// can't expect users of these classes to figure out that they have to include it, too.
// AUTO-REMOVED #include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/Frame.hh>
// AUTO-REMOVED #include <core/fragment/FragID_Iterator.fwd.hh>

// Project Headers
#include <core/types.hh>

// ObjexxFCL Headers

// Utility header
// AUTO-REMOVED #include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/pointer/ReferenceCount.hh>

// std Headers
// AUTO-REMOVED #include <iterator>

//Auto Headers
#include <core/fragment/FrameIteratorWorker_.hh>
#include <iostream>


/*
Might change the FrameIterator to return FrameOP instead of Frame&...
This has the advantage of bein closer to the std:: implementations...
for instance one could do a
FrameList my_list;
FrameIterator it=fragset.begin(), eit=fragset.end();
copy( it, eit,back_inserter( my_list ));

the it-> wouldn't change, it still returns a Frame*
the *it would change and return a Frame* ...

*/

namespace core {
namespace fragment {


class FrameIterator : public std::iterator< std::forward_iterator_tag, Frame > {
	friend class FragID_Iterator;

public:
	FrameIterator( FrameIteratorWorker_OP it );
	FrameIterator();
	~FrameIterator();

	bool operator != ( FrameIterator const& fi) const;
	bool operator == ( FrameIterator const& fi) const;
	FrameIterator & operator++ ();
	FrameIterator & operator+ ( Size offset );

	FrameIterator const & operator = ( FrameIterator const& itr );

	Frame* operator* ();

	Frame const* operator* () const;

	Frame* operator-> ();

	Frame const* operator-> () const;

	Frame* frame_ptr();

	Frame const* frame_ptr() const;

protected:
	FrameIteratorWorker_OP it_;

};

}
}

#endif
