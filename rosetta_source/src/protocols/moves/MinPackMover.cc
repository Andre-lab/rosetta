// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/moves/MinPackMover.cc
/// @brief  Implementation of the MinPackMover class; a wrapper class for invoking core::pack::min_pack
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/moves/MinPackMover.hh>
#include <protocols/moves/MinPackMoverCreator.hh>

// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/pack/min_pack.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
// AUTO-REMOVED #include <utility/string_util.hh>

// option key includes
// AUTO-REMOVED #include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


namespace protocols {
namespace moves {

using namespace core;
	using namespace basic::options;
	using namespace pack;
		using namespace task;
			using namespace operation;
	using namespace scoring;

using basic::Warning;
using basic::t_warning;
static basic::Tracer TR("protocols.moves.MinPackMover");


std::string
MinPackMoverCreator::keyname() const
{
	return MinPackMoverCreator::mover_name();
}

protocols::moves::MoverOP
MinPackMoverCreator::create_mover() const {
	return new MinPackMover;
}

std::string
MinPackMoverCreator::mover_name()
{
	return "MinPackMover";
}

MinPackMover::MinPackMover() :
	Mover("MinPackMover"),
	scorefxn_(0),
	task_(0),
	task_factory_(0),
	stochastic_pack_( false )
{}

MinPackMover::MinPackMover( std::string const & type_name ) :
	Mover( type_name ),
	scorefxn_(0),
	task_(0),
	task_factory_(0),
	stochastic_pack_( false )
{}

	// constructors with arguments
MinPackMover::MinPackMover(
	ScoreFunctionCOP scorefxn
) :
	Mover("MinPackMover"),
	scorefxn_( scorefxn ),
	task_( 0 ),
	task_factory_(0),
	stochastic_pack_( false )
{}


MinPackMover::MinPackMover(
	ScoreFunctionCOP scorefxn,
	PackerTaskCOP task
) :
	Mover("MinPackMover"),
	scorefxn_( scorefxn ),
	task_( task ),
	task_factory_(0),
	stochastic_pack_( false )
{}

MinPackMover::~MinPackMover(){}

MinPackMover::MinPackMover( MinPackMover const & other ) :
	//utility::pointer::ReferenceCount(),
	Mover( other ),
	stochastic_pack_( other.stochastic_pack_ )
{
	scorefxn_ = other.score_function();
	task_ = other.task();
	task_factory_ = other.task_factory();
}

void
MinPackMover::apply( Pose & pose )
{
	if ( scorefxn_ == 0 ) {
		Warning() << "undefined ScoreFunction -- creating a default one" << std::endl;
		scorefxn_ = ScoreFunctionFactory::create_score_function( STANDARD_WTS );
	}

	core::pack::task::PackerTaskCOP task;
	if ( task_factory_ ) {
		task = task_factory_->create_task_and_apply_taskoperations( pose );
	} else {
		runtime_assert( task_ );
		runtime_assert( task_is_valid( pose ) );
		task = task_;
		//core::pack::min_pack( pose, *scorefxn_, task_ );
		core::pack::stochastic_pack( pose, *scorefxn_, task_ );
	}

	if ( stochastic_pack_ ) {
		core::pack::stochastic_pack( pose, *scorefxn_, task );
	} else {
		core::pack::min_pack( pose, *scorefxn_, task );
	}

}

std::string
MinPackMover::get_name() const {
	return MinPackMoverCreator::mover_name();
}

///@brief when the PackerTask was not generated locally, verify compatibility with pose
///@details the pose residue types must be equivalent to the ones used to generate the ResidueLevelTasks, because of the way that prevent_repacking and its associated flags work
bool
MinPackMover::task_is_valid( Pose const & pose ) const
{
	if ( task_->total_residue() != pose.total_residue() ) return false;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		if ( ! task_->residue_task(i).is_original_type( &pose.residue_type(i) ) ) return false;
	}
	return true;
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
MinPackMover::parse_my_tag(
	TagPtr const tag,
	DataMap & datamap,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->getName() != "MinPackMover" ) {
		TR(t_warning) << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

///@brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
MinPackMover::parse_score_function(
	TagPtr const tag,
	DataMap const & datamap,
	Filters_map const &,
	Movers_map const &,
	Pose const &
)
{
	ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == 0 ) return;
	score_function( new_score_function );
}

///@brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
MinPackMover::parse_task_operations(
	TagPtr const tag,
	DataMap const & datamap,
	Filters_map const &,
	Movers_map const &,
	Pose const &
)
{
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0) return;
	task_factory( new_task_factory );
}

///@brief required in the context of the parser/scripting scheme
MoverOP
MinPackMover::fresh_instance() const
{
	return new MinPackMover;
}

///@brief required in the context of the parser/scripting scheme
MoverOP
MinPackMover::clone() const
{
	return new MinPackMover( *this );
}


// setters
void MinPackMover::score_function( ScoreFunctionCOP sf )
{
	runtime_assert( sf );
	scorefxn_ = sf;
}

void MinPackMover::task( task::PackerTaskCOP t ) { task_ = t; }

void MinPackMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf );
	task_factory_ = tf;
}

// accessors
ScoreFunctionCOP MinPackMover::score_function() const { return scorefxn_; }
PackerTaskCOP MinPackMover::task() const { return task_; }
TaskFactoryCOP MinPackMover::task_factory() const { return task_factory_; }

void MinPackMover::stochastic_pack( bool setting ) { stochastic_pack_ = setting; }
bool MinPackMover::stochastic_pack() const { return stochastic_pack_; }

} // moves
} // protocols

