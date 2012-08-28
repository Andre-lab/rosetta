// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief  SymDof data container
/// @file   core/conformation/symmetry/SymDof.hh
/// @author Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_SymDof_hh
#define INCLUDED_core_conformation_symmetry_SymDof_hh

// Utility headers
#include <core/conformation/symmetry/SymDof.fwd.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <core/types.hh>

// C++ headers
// AUTO-REMOVED #include <string>

#include <utility/vector1_bool.hh>


namespace core {
namespace conformation {
namespace symmetry {

// Available dofs
enum dof_type {
	X_DOF = 1,
	Y_DOF,
	Z_DOF,
	X_ANGLE_DOF,
	Y_ANGLE_DOF,
	Z_ANGLE_DOF
};

class SymDof {

	public:

	/// @brief constructor
	SymDof();

	/// @brief copy constructor
	SymDof( SymDof const & src );

	SymDof &
  operator=( SymDof const & src );

	~SymDof();

	void
	add_dof_from_string( utility::vector1< std::string > dof_string );

	// io
	void read( std::string dof_line);
	friend std::ostream & operator<< ( std::ostream & s, const SymDof & dof );

	// @details is df allowed to move?
	bool
	allow_dof( int df ) const;

	bool has_dof();

	// @details the lower boundary of range1
	core::Real
	range1_lower( int df ) const;

	// @details the upper boundary of range1
	core::Real
  range1_upper( int df ) const;

	// @details the lower boundary of range2
	core::Real
	range2_lower( int df ) const;

	// @details the upper boundary of range1
	core::Real
  range2_upper( int df ) const;

	// details Have a range1 been specified?
	bool
	has_range1( int df ) const;

	// details Have a range2 been specified?
	bool
	has_range2( int df ) const;

	// @details has a lower boundary of range1 been specified?
	bool
	has_range1_lower( int df ) const;

	// @details has a upper boundary of range1 been specified?
	bool
  has_range1_upper( int df ) const;

	// @details has a lower boundary of range2 been specified?
	bool
  has_range2_lower( int df ) const;

	// @details has a upper boundary of range2 been specified?
	bool
  has_range2_upper( int df ) const;

	// @detail return the direction( upstream or downstream )
	// of the jump for a dof
	int
	jump_direction( int df ) const;

	friend
	bool
	operator==(SymDof const & a, SymDof const & b);

	friend
	bool
	operator!=(SymDof const & a, SymDof const & b);

	private:

	utility::vector1< bool > allowed_dof_jumps_; // is a particular dof allowed to move?
	utility::vector1< Real > lower_range_dof_jumps1_; // store the lower boundary of range1
	utility::vector1< Real > upper_range_dof_jumps1_; // store the upper boundary of range1
	utility::vector1< Real > lower_range_dof_jumps2_; // store the lower boundary of range2
	utility::vector1< Real > upper_range_dof_jumps2_; // store the upper boundary of range2
	utility::vector1< bool > has_range1_lower_;	// Is there a lower bound on range1?
	utility::vector1< bool > has_range1_upper_; // Is there a upper bound on range1?
	utility::vector1< bool > has_range2_lower_; // Is there a lower bound on range2?
	utility::vector1< bool > has_range2_upper_; // Is there a upper bound on range2?
	utility::vector1< int > jump_dir_;	// store jump dir for each dof

};

} // symmetry
} // conformation
} // core
#endif
