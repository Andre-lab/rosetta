// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/filters/DesignBySecondaryStructureCreator.cc
/// @brief Design residues that don't match the predicted secondary structure.
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_flxbb_filters_designbysecondarystructurecreator_hh
#define INCLUDED_protocols_flxbb_filters_designbysecondarystructurecreator_hh

// unit headers

// package headers

// project headers
#include <core/pack/task/operation/TaskOperationCreator.hh>

namespace protocols {
namespace flxbb {
namespace filters {

class DesignBySecondaryStructureOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
  virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
  virtual std::string keyname() const;
};

}
}
}

#endif
