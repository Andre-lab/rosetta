// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loops_definers/LoopsExplicitDefiner.hh
/// @brief  A loops definer is creates a serialized loops list
/// @author Matthew O'Meara (mattjomear@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_LoopsExplicitDefiner_HH
#define INCLUDED_protocols_loops_loops_definers_LoopsExplicitDefiner_HH

// Unit Headers
#include <protocols/loops/loops_definers/LoopsDefiner.hh>
#include <protocols/loops/loops_definers/LoopsExplicitDefiner.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>

// Platform Headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>


namespace protocols {
namespace loops {
namespace loops_definers {

class LoopsExplicitDefiner : public LoopsDefiner {
public:

	LoopsExplicitDefiner();

	virtual
	~LoopsExplicitDefiner();

	LoopsExplicitDefiner(
		LoopsExplicitDefiner const & src);

	/// @brief Create another loops definer of the type matching the most-derived
	/// version of the class.
	virtual
	LoopsDefinerOP
	clone() const;


	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		moves::DataMap const & data,
		core::pose::Pose const &);


	SerializedLoopList
	apply(
		core::pose::Pose const &);

private:
	SerializedLoop
	parse_loop_tag(
		utility::tag::TagPtr const tag,
		std::string const & loops_name);


private:
	SerializedLoopList loop_list_;

};

// do not add any derived classes to this file, unless they are
// generalized abstract base classes and do not actually 'do any work'

} //namespace
} //namespace
} //namespace

#endif


