// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file SymRestrictTaskForDocking.cc
/// @brief When passed to a PackerTask, pack/design is limited to the protein interface
/// @author Ingemar Andre

#include <protocols/symmetric_docking/SymRestrictTaskForDocking.hh>

#include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
// Symmetry
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetricConformation.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>


#include <protocols/scoring/Interface.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>


namespace protocols {
namespace symmetric_docking {

using namespace core;
using namespace scoring;
using namespace pack;

SymRestrictTaskForDocking::SymRestrictTaskForDocking()
	: scorefxn_( 0 ),
		include_current_( true ),
		distance_( 0 )
{}

SymRestrictTaskForDocking::SymRestrictTaskForDocking(
	core::scoring::ScoreFunctionCOP scorefxn,
	bool include_current,
	core::Real distance
) : scorefxn_( scorefxn ),
		include_current_( include_current ),
		distance_( distance )
{}

SymRestrictTaskForDocking::~SymRestrictTaskForDocking(){}


task::operation::TaskOperationOP SymRestrictTaskForDocking::clone() const
{
	return new SymRestrictTaskForDocking( *this );
}

void
SymRestrictTaskForDocking::apply(
	core::pose::Pose const & pose,
	core::pack::task::PackerTask & task
) const
{
	task.initialize_from_command_line().restrict_to_repacking().or_include_current( include_current_ );

	assert( scorefxn_ != 0 );
	// (existing comment) /// why is this still necessary???
//	(*scorefxn_)(pose);
//	scorefxn_->accumulate_residue_total_energies( pose );

	runtime_assert( scorefxn_ != 0 );
	runtime_assert( distance_ );

	protocols::scoring::Interface interface( 1 );
	interface.distance( distance_ );
	interface.calculate( pose );
	interface.set_symmetric_pack( pose, &task );
}

} // namespace symmetric_docking
} // namespace protocols

