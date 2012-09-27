// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    RingConformationMover.hh
/// @brief   Declarations and simple accessor/mutator definitions for RingConformationMover.
/// @author  labonte

#ifndef INCLUDED_protocols_simple_moves_carbohydrates_RingConformationMover_HH
#define INCLUDED_protocols_simple_moves_carbohydrates_RingConformationMover_HH

// Unit headers
#include <protocols/simple_moves/carbohydrates/RingConformationMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>
#include <iostream>

namespace protocols {
namespace simple_moves {
namespace carbohydrates {

/// @details  Based on a given MoveMap, this mover selects movable carbohydrate residues and flips their rings to an
/// idealized energy minimum ring conformer.  Six-membered rings sample chair and twist-boat conformations; five-
/// membered rings sample envelope and half-chair.  The sampling is biased towards ideally lower-energy conformers.
/// @remarks  This class is a work in progress....
class RingConformationMover: public moves::Mover {
public:
	// Standard methods ////////////////////////////////////////////////////////
	/// @brief  Default constructor
	RingConformationMover();

	/// @brief  Copy constructor
	RingConformationMover(RingConformationMover const & object_to_copy);

	/// @brief  Constructor with MoveMap input option
	RingConformationMover(core::kinematics::MoveMapOP input_movemap);

	// Assignment operator
	RingConformationMover & operator=(RingConformationMover const & object_to_copy);

	// Destructor
	virtual ~RingConformationMover();


	// Standard Rosetta methods ////////////////////////////////////////////////
	// General methods
	/// @brief  Register options with the option system.
	static void register_options();

	/// @brief  Generate string representation of RingConformationMover for debugging purposes.
	void show(std::ostream & output=std::cout) const;

	// Insertion operator (overloaded so that RingConformationMover can be "printed" in PyRosetta).
	friend std::ostream & operator<<(std::ostream & output, RingConformationMover const & object_to_output);

	
	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief  Apply the corresponding move to <input_pose>.
	virtual void apply(core::pose::Pose & input_pose);


	// Accessors/Mutators
	/// @brief  Get the current MoveMap.
	core::kinematics::MoveMapCOP movemap() const;

	/// @brief  Set the MoveMap.
	void movemap(core::kinematics::MoveMapOP new_movemap);

private:
	// Private methods /////////////////////////////////////////////////////////
	// Initialize data members from arguments.
	void init(core::kinematics::MoveMapOP movemap);

	// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
	void copy_data(RingConformationMover object_to_copy_to, RingConformationMover object_to_copy_from);

	// Setup list of movable carbohydrate residues from MoveMap.
	void setup_residue_list(core::pose::Pose & pose);


	// Private data ////////////////////////////////////////////////////////////
	core::kinematics::MoveMapOP movemap_;
	utility::vector1<core::Size> residue_list_;  // list of movable carbohydrate residues

	// Lists of conformers as a pair:
	// first: the type of conformer (e.g., chair, twist-boat, etc.)
	// second: a list of torsion angle values for that conformation
	utility::vector1<std::pair<std::string, utility::vector1<core::Real> > > five_membered_ring_conformers_;
	utility::vector1<std::pair<std::string, utility::vector1<core::Real> > > six_membered_ring_conformers_;

};  // class RingConformationMover

}  // namespace carbohydrates
}  // namespace simple_moves
}  // namespace protocols

#endif  // INCLUDED_simple_moves_protocols_carbohydrates_RingConformationMover_HH
