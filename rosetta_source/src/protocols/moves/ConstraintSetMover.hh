// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Assigns a ConstraintSet to a pose. Reads and creats ConstraintSet from file via command line option -constraints::cst_file, unless a ConstraintSet is supplied via the constructor or the constraint_set() method.
/// @author ashworth

#ifndef INCLUDED_protocols_moves_ConstraintSetMover_hh
#define INCLUDED_protocols_moves_ConstraintSetMover_hh

#include <protocols/moves/ConstraintSetMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>


namespace protocols {
namespace moves {

class ConstraintSetMover : public Mover {

public:
	typedef core::scoring::constraints::ConstraintSet ConstraintSet;
	typedef core::scoring::constraints::ConstraintSetOP ConstraintSetOP;
	typedef core::scoring::constraints::ConstraintSetCOP ConstraintSetCOP;

public:
	ConstraintSetMover();
	virtual ~ConstraintSetMover();
	ConstraintSetMover( std::string const & );

	void read_options();
	void constraint_file( std::string const & );

	void constraint_set( ConstraintSetCOP );
	ConstraintSetOP constraint_set();
	ConstraintSetCOP constraint_set() const;

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual MoverOP clone() const;
	virtual MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagPtr const, DataMap &, Filters_map const &, Movers_map const &, Pose const & );

private:
	ConstraintSetOP constraint_set_low_res_;
	ConstraintSetOP constraint_set_high_res_;
	std::string cst_file_;
	std::string cst_fa_file_;

};

} // moves
} // protocols

#endif
