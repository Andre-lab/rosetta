// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/SnugDockProtocol.hh
/// @brief Dock and antigen to an antibody while optimizing the rigid body orientation of the VH and VL chains and performing CDR loop minimization.
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_antibody2_SnugDockProtocol_HH
#define INCLUDED_protocols_antibody2_SnugDockProtocol_HH

// Unit headers
#include <protocols/antibody2/SnugDockProtocol.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/antibody2/AntibodyInfo.fwd.hh>
#include <protocols/antibody2/RefineOneCDRLoop.fwd.hh>

// Project headers
#include <protocols/docking/DockingProtocol.fwd.hh>

// C++ headers
#include <iostream>

using namespace core;
namespace protocols {
namespace antibody2 {

class SnugDockProtocol: public moves::Mover {
public: // boiler plate / virtuals
	// default constructor
	SnugDockProtocol();

	// copy constructor
	SnugDockProtocol( SnugDockProtocol const & rhs );

	// assignment operator
	SnugDockProtocol & operator=( SnugDockProtocol const & rhs );

	// destructor
	virtual ~SnugDockProtocol();
	
	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	///@brief This mover retains state such that a fresh version is needed if the input Pose is about to change
	virtual bool reinitialize_for_new_input() const;

	/// @brief Associates relevant options with the SnugDockProtocol class
	static void register_options();

public:
	void show( std::ostream & out=std::cout );
	friend std::ostream & operator<<(std::ostream& out, SnugDockProtocol const & snugdockprotocol );

private: // methods
	void setup_objects( Pose const & pose );
	void setup_loop_refinement_movers( Pose const & pose );
	void init();
	void init_for_equal_operator_and_copy_constructor( SnugDockProtocol & lhs, SnugDockProtocol const & rhs);
	void init_options();

	docking::DockingProtocolOP docking() const;

private: // data
	AntibodyInfoOP antibody_info_;

	// Movers
	RefineOneCDRLoopOP low_res_refine_cdr_h2_;
	RefineOneCDRLoopOP low_res_refine_cdr_h3_;
	mutable docking::DockingProtocolOP docking_;
	
	std::string loop_refinement_method_;

}; // class SnugDockProtocol

} // namespace antibody2
} // namespace protocols

#endif // INCLUDED_protocols_antibody2_SnugDockProtocol_HH
