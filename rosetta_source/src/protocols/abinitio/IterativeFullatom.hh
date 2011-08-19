// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file AbrelaxMover
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @detailed responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_IterativeFullatom_hh
#define INCLUDED_protocols_abinitio_IterativeFullatom_hh

// Unit Headers
//#include <protocols/abinitio/IterativeFullatom.fwd.hh>

// Package Headers
#include <protocols/abinitio/IterativeBase.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

//// C++ headers
#include <string>


namespace protocols {
namespace abinitio {

class IterativeFullatom : public IterativeBase {
public:
	static void register_options();

  IterativeFullatom( jd2::archive::ArchiveManagerAP ptr );
	virtual bool ready_for_batch() const;

	virtual void generate_batch();

protected:
	void gen_resample_core( jd2::archive::Batch& batch, bool flex );

private:
	static bool options_registered_;
	core::Real perturb_start_structures_;

};


}
}

#endif
