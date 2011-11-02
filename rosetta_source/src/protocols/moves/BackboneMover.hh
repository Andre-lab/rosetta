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


#ifndef INCLUDED_protocols_moves_BackboneMover_hh
#define INCLUDED_protocols_moves_BackboneMover_hh

// Unit headers
#include <protocols/moves/BackboneMover.fwd.hh>

// Package headers
#include <protocols/moves/ThermodynamicMover.hh>

#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/id/DOF_ID_Range.fwd.hh>

#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// ObjexxFCL Headers

// C++ Headers
#include <map>
#include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

///////////////////////////////////////////////////////////////////////////////
/// BackboneMover class has elements of the MC temperature to do repetitions
/// of bb moves (small, shear, wobble, etc.).
/// @todo change this to some kind of 'protocol' so the MC is managed separately from
/// conformational moves
class BackboneMover : public ThermodynamicMover {

public:
	typedef core::Real Real;

public:

	// empty constructor fills the mover with default values
	// default values from smallmoves.cc of Rosetta++ (small_move_param)
	BackboneMover();

	BackboneMover(
		core::kinematics::MoveMapOP movemap_in,
		core::Real temperature_in,
		core::Size nmoves_in
	);

	//destructor
	~BackboneMover();

	/// virtual functions that get overloaded or called from the inheriting classes
	virtual void apply( core::pose::Pose & );
	virtual std::string get_name() const;

	virtual void setup_list( core::pose::Pose & ) = 0;

	virtual void set_angles( core::Real ) = 0;

	virtual bool make_move( core::pose::Pose & ) = 0;

	void clear();

	bool check_rama();

/// Properties set/get functions
	void temperature( core::Real const temperature_in );
	void nmoves( core::Size const nmoves_in );
	core::kinematics::MoveMapCOP movemap();

	// Because this function is in a .hh and not a .cc, we must #include <MoveMap.hh> and
	// are not able to #include <MoveMap.fwd.hh> alone.  As a general rule, do not put
	// function definitions in .hh files.  Break that rule only if you're able to demonstrate a
	// genuine inlining performance boost.
	void movemap(core::kinematics::MoveMapOP new_movemap);
	/// @brief Sets the maximum angle of perturbation,
	/// independent of secondary structure
	///
	/// example:
	///     bbmover.angle_max(25)
	/// See also:
	///     ShearMover
	///     SmallMover
	void angle_max( core::Real const angle );
	/// @brief Sets the max angle of perturbation, for residues with  <type>
	/// secondary structure (<type>  must be 'H', 'E', and 'L')
	///
	/// example:
	///     bbmover.angle_max("H",25)
	///
	/// See also:
	///     ShearMover
	///     SmallMover
	void angle_max( char const type, core::Real const angle );
	// note pass in by value for one direction assignment.
/// @brief sets the max angle of perturbation, for secondary structure 'H', 'E', and 'L'
	void angle_max( std::map< char, core::Real > angle_max_in );

	core::Real new_phi();
	core::Real new_psi();

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	/// @brief get whether detailed balance is preserved (i.e. no Ramachandran biasing)
	bool
	preserve_detailed_balance() const;

	/// @brief set whether detailed balance is preserved (i.e. no Ramachandran biasing)
	void
	set_preserve_detailed_balance(
		bool preserve_detailed_balance
	);

	/// @brief get the DOF_IDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::DOF_ID_Range>
	dof_id_ranges(
		core::pose::Pose & pose
	) = 0;

protected:

	core::kinematics::MoveMapOP movemap_;

	/// controls bias w/which uphill moves are accepted
	core::Real temperature_;

	/// number of positions at which to make moves
	Size nmoves_;

	/// max allowed angle-change as a function of ss type
	std::map< char, core::Real > angle_max_;

	// variables for the apply
	int num_, resnum_, tries_;
	core::Real big_angle_, small_angle_;
	utility::vector1< std::pair< int, core::Real > > pos_list_;
	utility::vector1< int > already_moved_;

	core::Real old_phi_, new_phi_, old_psi_, new_psi_;
	core::Real old_rama_score_, new_rama_score_;

	bool preserve_detailed_balance_;
};

///////////////////////////////////////////////////////////////////////////////


/// @brief A mover that makes independent random perturbations of the phi and
/// psi torsion angles of residue i. It selects residue i at random among
/// movable residues (set by its MoveMap), and the final torsion angle
/// is subject to a metropolis criterion using the rama score to ensure that
/// only favorable backbone torsion angles are being selected. The number of
/// perturbations, and the magnitude of perturbations, and the temperature
/// in the rama check, can all be modified.
///
/// Common Methods:
///     SmallMover.apply
///     SmallMover.angle_max
class SmallMover : public BackboneMover {

public:

	// default constructor
	SmallMover();
//	SmallMover() : BackboneMover() { Mover::type( "SmallMover" ); }

	/// @brief Constructs a SmallMover
	/// smallmover = SmallMover( movemap , kT , n_moves )
	///
	/// MoveMap        movemap   /object storing BB torsion movability
	/// Real (float)   kT        /used in rama Metropolis Criterion
	/// Size (int)     n_moves   /the number of perturbations for one move
	SmallMover(
		core::kinematics::MoveMapOP movemap_in,
		core::Real temperature_in,
		core::Size nmoves_in
	);

	//destructor
	~SmallMover();
	virtual std::string get_name() const;

	protocols::moves::MoverOP clone() const;

	virtual void setup_list( core::pose::Pose & pose );
	virtual void set_angles( core::Real angle_in );
	virtual bool make_move( core::pose::Pose & pose );

	virtual void test_move( core::pose::Pose & );

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges(
		core::pose::Pose & pose
	);

	/// @brief get the DOF_IDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::DOF_ID_Range>
	dof_id_ranges(
		core::pose::Pose & pose
	);
};


///////////////////////////////////////////////////////////////////////////////


/// @brief A mover that perturbs the phi of residue i and the psi of residue
/// i-1 such that they create a 'shearing' effect, minimizing the downstream
/// consequences of this torsional perturbation. The final torsion angle
/// is subject to a metropolis criterion using the rama score to ensure that
/// only favorable backbone torsion angles are being selected. The number of
/// perturbations, and the magnitude of perturbations, and the temperature
/// in the rama check, can all be modified.
///
/// Common Methods:
///     ShearMover.apply
class ShearMover : public BackboneMover {

public:

	// default constructor
	ShearMover();
	//ShearMover() : BackboneMover() { Mover::type( "ShearMover" );	}

	/// @brief Constructs a ShearMover
	/// shearmover = ShearMover( movemap , kT , n_moves )
	///
	/// MoveMap        movemap   /object storing BB torsion movability
	/// Real (float)   kT        /used in rama Metropolis Criterion
	/// Size (int)     n_moves   /the number of perturbations for one move
	ShearMover(
		core::kinematics::MoveMapOP movemap_in,
		core::Real temperature_in,
		core::Size nmoves_in
	);

	//destructor
	~ShearMover();
	virtual std::string get_name() const;

	protocols::moves::MoverOP clone() const;

	virtual void setup_list( core::pose::Pose & pose );
	virtual void set_angles( core::Real angle_in );
	virtual bool make_move( core::pose::Pose & pose );

	virtual void test_move( core::pose::Pose & );

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges(
		core::pose::Pose & pose
	);

	/// @brief get the DOF_IDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::DOF_ID_Range>
	dof_id_ranges(
		core::pose::Pose & pose
	);
};
} // moves
} // protocols


#endif //INCLUDED_protocols_moves_BackboneMover_HH

