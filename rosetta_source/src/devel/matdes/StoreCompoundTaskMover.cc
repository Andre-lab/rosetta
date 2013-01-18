// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/matdes/StoreCompoundTaskMover.cc
/// @brief  Combine tasks using boolean logic for residues that are packable or designable,
/// assign new packing behavior to residues that match or do not match the specified criteria, 
/// and store the resulting task in the pose's cacheable data.
/// TODO (Jacob): Add option to combine allowed amino acid sets.
/// @author Jacob Bale (balej@uw.edu) (Much of this code was adapted from the CompoundStatementFilter, StoreTaskMover, and rosetta_scripts/util.cc)

// Unit Headers
#include <devel/matdes/StoreCompoundTaskMover.hh>
#include <devel/matdes/StoreCompoundTaskMoverCreator.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/make_symmetric_task.hh>
#include <devel/matdes/STMStoredTask.hh>
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/format.hh>

// Boost Headers
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

static basic::Tracer TR("devel.matdes.StoreCompoundTaskMover");

namespace devel {
namespace matdes {

// @brief default constructor
StoreCompoundTaskMover::StoreCompoundTaskMover() {}

// @brief destructor
StoreCompoundTaskMover::~StoreCompoundTaskMover() {}

void
StoreCompoundTaskMover::clear()
{
	compound_task_.clear();
}

StoreCompoundTaskMover::iterator
StoreCompoundTaskMover::begin()
{
	return( compound_task_.begin() );
}
StoreCompoundTaskMover::const_iterator
StoreCompoundTaskMover::begin() const
{
	return( compound_task_.begin() );
}

StoreCompoundTaskMover::iterator
StoreCompoundTaskMover::end()
{
	return( compound_task_.end() );
}

StoreCompoundTaskMover::const_iterator
StoreCompoundTaskMover::end() const
{
	return( compound_task_.end() );
}

void
StoreCompoundTaskMover::invert( bool const inv )
{
	invert_ = inv;
}

void
StoreCompoundTaskMover::verbose( bool const verb )
{
	verbose_ = verb;
}

void
StoreCompoundTaskMover::overwrite( bool const ow )
{
	overwrite_ = ow;
}

void
StoreCompoundTaskMover::task_name( std::string const tn )
{
	task_name_ = tn;
}

void
StoreCompoundTaskMover::mode( std::string const md )
{
	mode_ = md;
}

void
StoreCompoundTaskMover::true_behavior( std::string const tb )
{
	true_behavior_ = tb;
}

void
StoreCompoundTaskMover::false_behavior( std::string const fb )
{
	false_behavior_ = fb;
}

void
StoreCompoundTaskMover::CompoundPackableTask( core::Size & total_residue, core::pack::task::PackerTaskOP & task)
{
	std::string select_true_pos(task_name_+": select true_positions, resi ");
	for( core::Size resi=1; resi<=total_residue; ++resi ){

		bool value( true );

		for( StoreCompoundTaskMover::const_iterator it=compound_task_.begin(); it!=compound_task_.end(); ++it ) {
			if( it - compound_task_.begin() == 0 ){
				// first logical op may only be NOT
				// ANDNOT and ORNOT are also treated as NOT (with a warning)
				value = it->first->being_packed( resi );
				if (it->second == NOT) value = !value;
				if (it->second == ORNOT) {
					TR << "WARNING: StoreCompoundTaskMover treating operator ORNOT as NOT" << std::endl;
					value = !value;
				}
				if (it->second == ANDNOT) {
					TR << "WARNING: StoreCompoundTaskMover treating operator ANDNOT as NOT" << std::endl;
					value = !value;
				}
			} else {
				switch( it->second  ) {
					case ( AND ) : value = value && it->first->being_packed( resi ); break;
					case ( OR  ) : value = value || it->first->being_packed( resi ); break;
					case ( XOR ) : value = value ^ it->first->being_packed( resi ); break;
					case ( ORNOT ) : value = value || !it->first->being_packed( resi ); break;
					case ( ANDNOT ) : value = value && !it->first->being_packed( resi ); break;
					case ( NOR ) : value = !( value || it->first->being_packed( resi ) ); break;
					case (NAND ) : value = !( value && it->first->being_packed( resi ) ); break;
					case (NOT ) :
						TR << "WARNING: StoreCompoundTaskMover treating operator NOT as ANDNOT" << std::endl;
						value = value && !it->first->being_packed( resi );
						break;
				}
			}
		}
		if( invert_ ) value = !value;
		if( !value ) { 
			if( false_behavior_ == "prevent_repacking" ) {
				task->nonconst_residue_task(resi).prevent_repacking();
			} else if( false_behavior_ == "restrict_to_repacking" ) { 
				task->nonconst_residue_task(resi).restrict_to_repacking();
			}
		} else { 
			if( true_behavior_ == "prevent_repacking" ) {
				task->nonconst_residue_task(resi).prevent_repacking();
			} else if( true_behavior_ == "restrict_to_repacking" ) { 
				task->nonconst_residue_task(resi).restrict_to_repacking();
			}
			select_true_pos.append(ObjexxFCL::string_of(resi) + "+");   
		}
	}
	if( verbose_ ) {
		TR << select_true_pos << std::endl;
	}	
}

void
StoreCompoundTaskMover::CompoundDesignableTask( core::Size & total_residue, core::pack::task::PackerTaskOP & task)
{
	std::string select_true_pos(task_name_+": select true_positions, resi ");
	for( core::Size resi=1; resi<=total_residue; ++resi ){

		bool value( true );

		for( StoreCompoundTaskMover::const_iterator it=compound_task_.begin(); it!=compound_task_.end(); ++it ) {
			if( it - compound_task_.begin() == 0 ){
				// first logical op may only be NOT
				// ANDNOT and ORNOT are also treated as NOT (with a warning)
				value = it->first->being_designed( resi );
				if (it->second == NOT) value = !value;
				if (it->second == ORNOT) {
					TR << "WARNING: StoreCompoundTaskMover treating operator ORNOT as NOT" << std::endl;
					value = !value;
				}
				if (it->second == ANDNOT) {
					TR << "WARNING: StoreCompoundTaskMover treating operator ANDNOT as NOT" << std::endl;
					value = !value;
				}
			} else {
				switch( it->second  ) {
					case ( AND ) : value = value && it->first->being_designed( resi ); break;
					case ( OR  ) : value = value || it->first->being_designed( resi ); break;
					case ( XOR ) : value = value ^ it->first->being_designed( resi ); break;
					case ( ORNOT ) : value = value || !it->first->being_designed( resi ); break;
					case ( ANDNOT ) : value = value && !it->first->being_designed( resi ); break;
					case ( NOR ) : value = !( value || it->first->being_designed( resi ) ); break;
					case (NAND ) : value = !( value && it->first->being_designed( resi ) ); break;
					case (NOT ) :
						TR << "WARNING: StoreCompoundTaskMover treating operator NOT as ANDNOT" << std::endl;
						value = value && !it->first->being_designed( resi );
						break;
				}
			}
		}
		if( invert_ ) value = !value;
		if( !value ) { 
			if( false_behavior_ == "prevent_repacking" ) {
				task->nonconst_residue_task(resi).prevent_repacking();
			} else if( false_behavior_ == "restrict_to_repacking" ) { 
				task->nonconst_residue_task(resi).restrict_to_repacking();
			}
		} else { 
			if( true_behavior_ == "prevent_repacking" ) {
				task->nonconst_residue_task(resi).prevent_repacking();
			} else if( true_behavior_ == "restrict_to_repacking" ) { 
				task->nonconst_residue_task(resi).restrict_to_repacking();
			}
			select_true_pos.append(ObjexxFCL::string_of(resi) + "+");   
		}
	}
	if( verbose_ ) {
		TR << select_true_pos << std::endl;
	}	
}

void
StoreCompoundTaskMover::apply( core::pose::Pose & pose )
{

	core::Size total_residue;
	if(core::pose::symmetry::is_symmetric( pose )) { 
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		total_residue = symm_info->num_independent_residues();		
	} else {
		total_residue = pose.total_residue(); 
	}

	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );

	////////////////////////////////////////////////////////////////////////////////////////////////
	/// For each of the possible modes: 
	/// 1) loop over all of the independent residues
	/// 2) for each residue
	///		1) loop through the task/boolean_operation pairs
	///		2) determine the final state of each residue (ie, packable, designable, designable with amino acids TWER... etc)
	///		3) apply this state to the new compound task
	////////////////////////////////////////////////////////////////////////////////////////////////

	if( mode_ == "packable" ) {
		CompoundPackableTask( total_residue, task );
	} else if( mode_ == "designable" ) {
		CompoundDesignableTask( total_residue, task );
	}

/* TODO (Jacob): Add option to combine aa_sets. 
	Does a getter function exist to get a length-20 vector of bools for allowed_aas for a given residue?
	Does a function exist to check if an amino acid (given by either one or three letter aa code) is allowed at a given position?
	Add addition loop over allowed_aa vector for each residue, storing value in a vector of bools, which is then used to set
	the new set of allowed amino acids via restrict_absent_canonical_aas( final vector of bools )
	} else if( mode_ == "allowed_aas" ) {
	// utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false )
	// for(ResidueLevelTask::ResidueTypeCOPListConstIter aa_iter(storage_task->residue_task(i).allowed_residue_types_begin()),
	// aa_end(storage_task->residue_task(i).allowed_residue_types_end());
	// aa_iter != aa_end; ++aa_iter){
	// residues_to_mutate[i]=((*aa_iter)->aa());

	}
*/
	//////////////////////////////////////////////////////////////////////////////////////////////////
	/// Store the new compound task in the pose's cacheable data.
	/////////////////////////////////////////////////////////////////////////////////////////////////
	if (core::pose::symmetry::is_symmetric(pose))
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // Does this need to be fixed or omitted?
	// If the pose doesn't have STM_STORED_TASK data, put a blank STMStoredTask in there.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) {
		devel::matdes::STMStoredTaskOP blank_tasks = new devel::matdes::STMStoredTask();
		pose.data().set( core::pose::datacache::CacheableDataType::STM_STORED_TASKS, blank_tasks );
	}
	// Grab a reference to the data
	devel::matdes::STMStoredTask & stored_tasks = *( static_cast< devel::matdes::STMStoredTask* >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS )() ) );
	// If you haven't set overwrite to true and your task name already exists, fail. Otherwise, put the task you've made into the data cache.
	if ( overwrite_ || !stored_tasks.has_task(task_name_) ) {
		stored_tasks.set_task( task, task_name_ );
	} else {
		utility_exit_with_message("A stored task with the name " + task_name_ + " already exists; you must set overwrite flag to true to overwrite." );
	}
}

void
StoreCompoundTaskMover::parse_my_tag( TagPtr const tag, protocols::moves::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{

  typedef utility::vector1< std::string > StringVec;

	TR<<"StoreCompoundTask: "<<tag->getName()<<std::endl;
	task_name_ = tag->getOption< std::string >( "task_name", "" );
	true_behavior_ = tag->getOption<std::string>( "true_behavior", "" );
	if( !( (true_behavior_ == "prevent_repacking") || (true_behavior_ == "restrict_to_repacking") || (true_behavior_ == "") ) ) { 
		throw utility::excn::EXCN_RosettaScriptsOption( "Error: true_behavior in tag is undefined." );
	}
	false_behavior_ = tag->getOption<std::string>( "false_behavior", "prevent_repacking" );
	if( !( (false_behavior_ == "prevent_repacking") || (false_behavior_ == "restrict_to_repacking") || (false_behavior_ == "") ) ) { 
		throw utility::excn::EXCN_RosettaScriptsOption( "Error: false_behavior in tag is undefined." );
	}
	mode_ = tag->getOption<std::string>( "mode", "packable" );
	if( !( (mode_ == "packable") || (mode_ == "designable") ) ) { 
		throw utility::excn::EXCN_RosettaScriptsOption( "Error: mode in tag is undefined." );
	}
	invert_ = tag->getOption<bool>( "invert", false );
	verbose_ = tag->getOption<bool>( "verbose", false );
	overwrite_ = tag->getOption< bool >( "overwrite", false );

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// Loop through all user-provided subtags (ex. < AND task_name="bbi" />) and put these into a 
	/// vector of (PackerTaskOP, boolean_operation) pairs.
	/////////////////////////////////////////////////////////////////////////////////////////////////
	foreach(TagPtr cmp_tag_ptr, tag->getTags() ){
		std::string const operation( cmp_tag_ptr->getName() );
		std::pair< core::pack::task::PackerTaskOP, boolean_operations > task_pair;
		if( operation == "AND" ) task_pair.second = AND;
		else if( operation == "OR" ) task_pair.second = OR;
		else if( operation == "XOR" ) task_pair.second = XOR;
		else if( operation == "NOR" ) task_pair.second = NOR;
		else if( operation == "NAND" ) task_pair.second = NAND;
		else if( operation == "ORNOT" ) task_pair.second = ORNOT;
		else if( operation == "ANDNOT" ) task_pair.second = ANDNOT;
		else if( operation == "NOT" ) task_pair.second = NOT;
		else {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: Boolean operation in tag is undefined." );
		}

  	core::pack::task::TaskFactoryOP new_task_factory = new core::pack::task::TaskFactory;
  	std::string const t_o_val( cmp_tag_ptr->getOption<std::string>("task_operations") );
 
	 	TR<<"Defined with operator: "<<operation<<" and tasks: "<<t_o_val<<std::endl;

  	StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
  	for ( StringVec::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
  	      t_o_key != end; ++t_o_key ) {
    	if ( data_map.has( "task_operations", *t_o_key ) ) {
      	new_task_factory->push_back( data_map.get< core::pack::task::operation::TaskOperation * >( "task_operations", *t_o_key ) );
    	} else {
      	utility_exit_with_message("TaskOperation " + *t_o_key + " not found in DataMap.");
    	}
  	}
		core::pack::task::PackerTaskOP new_packer_task = new_task_factory->create_task_and_apply_taskoperations( pose );
		task_pair.first = new_packer_task->clone(); //clone?
		runtime_assert( new_packer_task );
		compound_task_.push_back( task_pair );
	}
}

// @brief Identification
std::string StoreCompoundTaskMoverCreator::keyname() const { return StoreCompoundTaskMoverCreator::mover_name(); }
std::string StoreCompoundTaskMoverCreator::mover_name() { return "StoreCompoundTaskMover"; }
std::string StoreCompoundTaskMover::get_name() const { return "StoreCompoundTaskMover"; }

protocols::moves::MoverOP
StoreCompoundTaskMoverCreator::create_mover() const {
	return new StoreCompoundTaskMover;
}

protocols::moves::MoverOP
StoreCompoundTaskMover::clone() const {
	return new StoreCompoundTaskMover( *this );
}

protocols::moves::MoverOP
StoreCompoundTaskMover::fresh_instance() const {
	return new StoreCompoundTaskMover;
}

} // matdes
} // devel

