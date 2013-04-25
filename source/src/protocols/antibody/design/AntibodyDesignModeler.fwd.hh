// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDesignModeler.fwd.hh
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignModeler_fwd_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignModeler_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// C++ Headers
namespace protocols {
namespace antibody {
namespace design {

// Forward
class AntibodyDesignModeler;

typedef utility::pointer::owning_ptr< AntibodyDesignModeler > AntibodyDesignModelerOP;
typedef utility::pointer::owning_ptr< AntibodyDesignModeler const > AntibodyDesignModelerCOP;


} //namespace design
} //namespace antibody
} //namespace protocols



#endif	//INCLUDED_protocols_antibody_design_ANTIBODYDESIGNMODELER_FWD_HH

