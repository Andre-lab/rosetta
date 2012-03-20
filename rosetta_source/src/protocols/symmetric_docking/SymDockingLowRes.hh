// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file SymDockingLowRes
/// @brief low resolution mode for symmetrical docking
/// @author Ingemar Andre


#ifndef INCLUDED_protocols_symmetric_docking_SymDockingLowRes_hh
#define INCLUDED_protocols_symmetric_docking_SymDockingLowRes_hh

// Package headers

#include <protocols/symmetric_docking/SymDockingLowRes.fwd.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>

#include <protocols/moves/MoverContainer.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace symmetric_docking {

class SymDockingLowRes : public moves::Mover
{
	typedef core::Real Real;

public:

	// default constructor
	//SymDockingLowRes();

	// constructor with arguments
	SymDockingLowRes( core::scoring::ScoreFunctionCOP scorefxn_in );
	// default constructor

		moves::MoverOP clone() const;

	~SymDockingLowRes();

	void set_default( core::pose::Pose & pose );
	void set_default_mc( core::pose::Pose & pose );
	void set_default_protocol( core::pose::Pose & pose );
	void set_default_move_map( core::pose::Pose & pose );

	moves::MonteCarloOP get_mc();

	// protocol functions
	virtual void apply( core::pose::Pose & pose );
	void rigid_body_trial( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	// protocol stuff
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::kinematics::MoveMapOP movemap_;
	moves::MonteCarloOP mc_;
	moves::SequenceMoverOP docking_lowres_protocol_;
	rigid::RigidBodyDofSeqPerturbMoverOP rb_mover_;

	// docking
	int inner_cycles_, outer_cycles_;
	Real trans_magnitude_, rot_magnitude_, accept_rate_;
	bool chi_, bb_, nb_list_;

	Real temperature_;

};

} // symmetric_docking
} // protocols

#endif
