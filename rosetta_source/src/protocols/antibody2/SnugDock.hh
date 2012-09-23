// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/SnugDock.hh
/// @brief Dock and antigen to an antibody while optimizing the rigid body orientation of the VH and VL chains and performing CDR loop minimization.
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_antibody2_SnugDock_HH
#define INCLUDED_protocols_antibody2_SnugDock_HH

// Unit headers
#include <protocols/antibody2/SnugDock.fwd.hh>
#include <protocols/docking/DockingHighRes.hh>

// Package headers
#include <protocols/antibody2/AntibodyInfo.fwd.hh>
#include <protocols/antibody2/RefineOneCDRLoop.fwd.hh>
#include <protocols/antibody2/CDRsMinPackMin.fwd.hh>


// Project headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.fwd.hh>

// C++ headers
#include <iostream>

using namespace core;
namespace protocols {
namespace antibody2 {

class SnugDock: public docking::DockingHighRes {
public: // boiler plate / virtuals
	// default constructor
	SnugDock();

	// copy constructor
	SnugDock( SnugDock const & rhs );

	// assignment operator
	SnugDock & operator=( SnugDock const & rhs );

	// destructor
	virtual ~SnugDock();
	
	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	///@brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	virtual bool reinitialize_for_new_input() const;

	/// @brief Associates relevant options with the SnugDock class
	static void register_options();

public:
	Size number_of_high_resolution_cycles() const;
	void number_of_high_resolution_cycles( Size const number_of_high_resolution_cycles );
	void set_antibody_info( AntibodyInfoOP antibody_info );

	void show( std::ostream & out=std::cout );
	friend std::ostream & operator<<(std::ostream& out, SnugDock const & snugdock );

private: // methods
	void setup_objects( Pose const & pose );
	void init();
	void init_for_equal_operator_and_copy_constructor( SnugDock & lhs, SnugDock const & rhs);
	void init_options();

private: // data
	AntibodyInfoOP antibody_info_;
	moves::MonteCarloOP mc_;

	// Movers
	moves::RandomMoverOP high_resolution_step_;
	std::string loop_refinement_method_;
	antibody2::CDRsMinPackMinOP pre_minimization_;

private:
	Size number_of_high_resolution_cycles_;

}; // class SnugDock

} // namespace antibody2
} // namespace protocols

#endif // INCLUDED_protocols_antibody2_SnugDock_HH
