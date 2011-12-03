// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author J. Karanicolas

#ifndef INCLUDED_protocols_loops_LoopMover_Backrub_hh
#define INCLUDED_protocols_loops_LoopMover_Backrub_hh


#include <protocols/loops/LoopMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {


class LoopMover_Refine_Backrub: public LoopMover {
public:
	LoopMover_Refine_Backrub();

	LoopMover_Refine_Backrub(
		protocols::loops::Loops  loops_in
	);

	LoopMover_Refine_Backrub(
		protocols::loops::Loops  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn
	);

	//destructor
	~LoopMover_Refine_Backrub();

	virtual std::string get_name() const;

	void set_default_settings(){
		redesign_loop = false;
	}

	void set_redesign_loop( bool value = true ){ redesign_loop = value; }
	bool get_redesign_loop(){ return redesign_loop; }

	void set_task_factory( core::pack::task::TaskFactoryOP value );
	bool get_task_factory();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const {
		return new LoopMover_Refine_Backrub(*this); // <--- TaskFactory.hh has to be #included here because this class's copy constructor is undefined, and this function is being invoked in the header
	}

	void apply( core::pose::Pose & pose );

protected:

	core::pack::task::TaskFactoryOP task_factory;
	bool redesign_loop;
};



} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_LoopMover_Backrub_HH
