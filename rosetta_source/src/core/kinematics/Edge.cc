// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/Edge.cc
/// @brief  Fold tree edge class
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/Edge.hh>

// AUTO-REMOVED #include <utility/stream_util.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>
#include <sstream>


// will remove the following in a couple of weeks OL/ 8/4/2008
// just to have some backward compatibiliy for a transient time
#include <basic/options/option.hh>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>


namespace core {
namespace kinematics {

static basic::Tracer tr("core.kinematics");

// PHIL -- need to update these routines to include atom info?
// ANDREW -- made a stab at adding atom info to comparison operators and ostream operators; no guarantee

/////////////////////////////////////////////////////////////////////////////
// these two should be inverses:

std::ostream &
operator <<( std::ostream & os, const Edge & e )
{
	std::string tag ( "EDGE" );
	if ( e.is_jump() && e.has_atom_info() ) tag = "JEDGE";
	os << " " << tag << " " << e.start() << ' ' << e.stop() << ' ' << e.label() << ' ';
	if ( e.label() == Edge::CHEMICAL ) os << e.start_atom() << ' ' << e.stop_atom() << ' ';
	if ( e.is_jump() ) {
		if ( e.start_atom().size() ) {
			os << e.start_atom() << ' ' << e.stop_atom() << ' ';
		} else {
			//		os << " X X "; not-necessary with JEDGE tag //otherwise reading becomes difficult
		}
	}
	if ( e.is_jump() && e.has_atom_info() ) {
		if ( e.keep_stub_in_residue() ) {
			os << " INTRA_RES_STUB ";
			assert( e.start_atom().size() );
		} else {
			os << " END ";
		}
	}
	return os;
}

/////////////////////////////////////////////////////////////////////////////
std::istream &
operator >>( std::istream & is, Edge & e )
{
	std::string tag;

//============================================
// this is a temporary hack to allow read of the fold-tree format that existed for the weeks from revision 21447 30/3 -> 8/4 2008
// in a couple of weeks this support shall be removed
 	if ( basic::options::option[ basic::options::OptionKeys::in::use_stupid_foldtree_format ] ) {

		is >> tag;
		if ( ! (tag == "EDGE" ) ) {
			tr.Trace << "failed reading EDGE tag --- found instead: " << tag << std::endl;
			is.setstate( std::ios_base::failbit );
			return is;
		}
		is >> e.start_ >> e.stop_ >> e.label_;
		if ( e.label() == Edge::CHEMICAL ) is >> e.start_atom_ >> e.stop_atom_;

		e.start_atom_ = ""; e.stop_atom_ = "";
		if ( e.is_jump() ) {
			is >> e.start_atom_;
			is >> e.stop_atom_;
		} else return is;
	} else { // normal ( new ) format
//=================================================
// END OF TEMPORARY HACK

		is >> tag;
		if ( ! (tag == "EDGE" || tag == "JEDGE") ) {
			tr.Trace << "failed reading EDGE tag --- found instead: " << tag << std::endl;
			is.setstate( std::ios_base::failbit );
			return is;
		}
		is >> e.start_ >> e.stop_ >> e.label_;
		if ( e.label() == Edge::CHEMICAL ) is >> e.start_atom_ >> e.stop_atom_;

		e.start_atom_ = ""; e.stop_atom_ = "";
		if ( e.is_jump() && tag == "JEDGE" ) {
			is >> e.start_atom_;
			is >> e.stop_atom_;
		} else return is;
	} // from here both pathways are the same.
	if ( e.start_atom_ == "X" ) e.start_atom_ = "";
	if ( e.stop_atom_ == "X" ) e.stop_atom_ = "";

	// allow either both atoms set or both unset
	assert( ( e.start_atom_.size() && e.stop_atom_.size() )
		|| ( e.start_atom_.size() == 0) && (e.stop_atom_.size() ==0 ));

	is >> tag;
	e.bKeepStubInResidue_ = false;
	if ( tag == "END" ) return is;
	assert( tag == "INTRA_RES_STUB" ); //only allowed if also atoms are specified;
	e.bKeepStubInResidue_ = true;
	return is;
}


/////////////////////////////////////////////////////////////////////////////
/// @details compare start residue number first, then stop residue number, then label
/// index number, then start_atom, then stop_atom
bool
operator <( Edge const & a, Edge const & b )
{
	//return ( a.start() <  b.start() ||
	//	a.start() == b.start() && a.stop() <  b.stop() ||
	//	a.start() == b.start() && a.stop() == b.stop() && a.label() < b.label() ||
	//	a.start() == b.start() && a.stop() == b.stop() && a.label() == b.label() && a.start_atom() < b.start_atom() ||
	//	a.start() == b.start() && a.stop() == b.stop() && a.label() == b.label() && a.start_atom() == b.start_atom() && a.stop_atom() < b.stop_atom() );
	//);
	return ( a.start() == b.start() ? ( a.stop() == b.stop() ? ( a.label() == b.label() ?
		( a.start_atom() == b.start_atom() ? a.stop_atom() < b.stop_atom() : a.start_atom() < b.start_atom() ) :
		a.label() < b.label() ) : a.stop() < b.stop() ) : a.start() < b.start() );
}


/////////////////////////////////////////////////////////////////////////////
/// @details when start residue number, stop residue number and label index number are all equal
bool
operator ==( Edge const & a, Edge const & b )
{
	return ( a.start() == b.start() && a.stop() == b.stop() && a.label() == b.label()
		&& a.start_atom() == b.start_atom() && a.stop_atom() == b.stop_atom() );
}


/////////////////////////////////////////////////////////////////////////////
/// @details when any of start residue number, stop residue number and label index number is not equal
bool
operator !=( Edge const & a, Edge const & b )
{
	return ( a.start() != b.start() || a.stop() != b.stop() || a.label() != b.label() );
}


} // namespace kinematics
} // namespace core
