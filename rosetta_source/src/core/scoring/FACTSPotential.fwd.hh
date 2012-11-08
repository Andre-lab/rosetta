// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file:   core/scoring/FACTSPotential.fwd.hh
// @author: Hahnbeom Park


#ifndef INCLUDED_core_scoring_FACTSPotential_fwd_hh
#define INCLUDED_core_scoring_FACTSPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

  //Declaring a class of type FACTSResidueInfo
  class FACTSResidueInfo;
	//Creating an alias for a pointer of type FACTSResidueInfo
  typedef utility::pointer::owning_ptr< FACTSResidueInfo > FACTSResidueInfoOP;

	//Declaring a class of type FACTSPoseInfo
	class FACTSPoseInfo;
	//Creating an alias for a pointer of type FACTSPoseInfo
	typedef utility::pointer::owning_ptr< FACTSPoseInfo > FACTSPoseInfoOP;

	class FACTSRotamerSetInfo;
	typedef utility::pointer::owning_ptr< FACTSRotamerSetInfo > FACTSRotamerSetInfoOP;

	//Declaring a class of type FACTSPotential
	class FACTSPotential;

  typedef utility::pointer::owning_ptr< FACTSPotential > FACTSPotentialOP;
	typedef utility::pointer::owning_ptr< FACTSPotential const > FACTSPotentialCOP;

} // scoring
} // core

#endif
