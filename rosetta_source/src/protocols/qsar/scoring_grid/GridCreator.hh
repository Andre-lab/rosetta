// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/qsar/scoring_grid/GridCreator.hh
/// @brief Base class for GridCreators for the Grid load time factory registration scheme
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_qsar_scoring_grid_GridCreator_HH_
#define INCLUDED_protocols_qsar_scoring_grid_GridCreator_HH_

#include <core/types.hh>
#include <protocols/qsar/scoring_grid/GridBase.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {

/// @brief Abstract class fora  mover factory.  The creator class is responsible
/// for creating a particular mover class
class GridCreator : public utility::pointer::ReferenceCount
{
public:
	GridCreator();
	virtual ~GridCreator();

	virtual GridBaseOP create_grid(utility::tag::TagPtr const tag) const = 0;
	virtual std::string keyname() const = 0;
private:
	core::Real weight_;
};

typedef utility::pointer::owning_ptr<GridCreator> GridCreatorOP;
typedef utility::pointer::owning_ptr<GridCreator const> GridCreatorCOP;

}
}
}

#endif /* GRIDCREATOR_HH_ */
