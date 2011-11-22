// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/protein_interface_design/movers/TryRotamers.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/TryRotamers.hh>
#include <protocols/protein_interface_design/movers/TryRotamersCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <protocols/protein_interface_design/util.hh>
#include <utility/string_util.hh>

#include <core/graph/Graph.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <basic/Tracer.hh>
#include <protocols/moves/DataMap.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.TryRotamers" );

std::string
TryRotamersCreator::keyname() const
{
	return TryRotamersCreator::mover_name();
}

protocols::moves::MoverOP
TryRotamersCreator::create_mover() const {
	return new TryRotamers;
}

std::string
TryRotamersCreator::mover_name()
{
	return "TryRotamers";
}


TryRotamers::TryRotamers() :
	protocols::moves::Mover( TryRotamersCreator::mover_name() )
{ }

TryRotamers::TryRotamers( Size resnum,
	ScoreFunction const & scorefxn,
	protocols::filters::Filter const & final_filter,
	Size explosion, // rotamer explosion
	Size jump_num,
	bool clash_check
) :
	protocols::moves::Mover( TryRotamersCreator::mover_name() ),
	scorefxn_(new ScoreFunction(scorefxn)),
	resnum_(resnum),
	jump_num_(jump_num),
	clash_check_(clash_check),
	explosion_(explosion),
	final_filter_(final_filter.clone())
{}

	/// @note Pass everything through the final filter (True Filter)
TryRotamers::TryRotamers( core::Size resnum,
	core::scoring::ScoreFunction const& scorefxn,
	Size explosion, // rotamer explosion
	Size jump_num,
	bool clash_check
	) :
	protocols::moves::Mover(),
	scorefxn_(new ScoreFunction(scorefxn)),
	resnum_(resnum),
	jump_num_(jump_num),
	clash_check_(clash_check),
	explosion_(explosion),
	final_filter_(new protocols::filters::TrueFilter)
{}

TryRotamers::~TryRotamers() {}

void
TryRotamers::apply ( pose::Pose & pose )
{
	//using namespace rotamer_set;
	using namespace core::scoring;
	using namespace core::pack::task;
	using namespace core::pack::rotamer_set;

	core::kinematics::FoldTree const saved_ft( pose.fold_tree() );

	Residue const & res_ = pose.residue( resnum_);

	if( automatic_connection_ ){
		core::kinematics::FoldTree const new_ft( make_hotspot_foldtree( pose ) );
		TR<<"New foldtree:\n"<<new_ft<<std::endl;
		pose.fold_tree( new_ft );
	}

	foreach( core::Size const resid, shove_residues_ )
		core::pose::add_variant_type_to_pose_residue( pose, "SHOVE_BB", resid );

	pose.update_residue_neighbors();

 if( !rotset_ || rotset_->num_rotamers() == 0 )
  {
	PackerTaskOP ptask( TaskFactory::create_packer_task( pose ) );
	ptask->set_bump_check( clash_check_ );

	ResidueLevelTask & restask( ptask->nonconst_residue_task( resnum_ ) );
	graph::GraphOP packer_graph = new graph::Graph( pose.total_residue() );
		  restask.or_ex1( true );
      restask.or_ex2( true );
      restask.or_ex3( true );
      restask.or_ex4( true );
			if( explosion_ > 0 ) restask.or_ex1_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
			if( explosion_ > 1	) restask.or_ex2_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
			if( explosion_ > 2	) restask.or_ex3_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
			if( explosion_ > 3	) restask.or_ex4_sample_level(core::pack::task::EX_FOUR_HALF_STEP_STDDEVS);
			restask.or_include_current( true );

	 restask.restrict_to_repacking();
	RotamerSetFactory rsf;
	rotset_ = rsf.create_rotamer_set( res_ );

	rotset_->set_resid( resnum_ );
	rotset_->build_rotamers( pose, *scorefxn_, *ptask, packer_graph, false );
	rotamer_it_ = rotset_->begin();
	TR<<"building rotamer set of " <<rotset_->num_rotamers()<< " different rotamers ...\n";

	}//end building rotamer set

	// job distributor iterates ...
	if ( rotamer_it_ == rotset_->end() )
	{
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
		TR<<"reached the end of the rotamer ensemble "  <<std::endl;
		pose.fold_tree( saved_ft );
		return;
	}

	if ( final_filter_->apply( pose ) )
	{
		if( jump_num_ > 0 ) { // Let 0 indicate no jumps
		core::kinematics::Jump const saved_jump( pose.jump( jump_num_ ) );
		pose.replace_residue ( resnum_, **rotamer_it_, false/*orient bb*/ );
			pose.set_jump_now( jump_num_, saved_jump );
		}
		else { // no jumps to save
			pose.replace_residue ( resnum_, **rotamer_it_, false/*orient bb*/ );
		}
		TR << "TryRotamers passed the final_filter" <<std::endl;
		set_last_move_status ( protocols::moves::MS_SUCCESS );
		++rotamer_it_;
		pose.fold_tree( saved_ft );
		return;
	}
	else
	{
		set_last_move_status( protocols::moves::FAIL_RETRY );
		++rotamer_it_;
		TR << std::endl;
	}
	pose.fold_tree( saved_ft );
}

std::string
TryRotamers::get_name() const {
	return TryRotamersCreator::mover_name();
}

void
TryRotamers::parse_my_tag( TagPtr const tag,
		DataMap & data,
		protocols::filters::Filters_map const &filters,
		Movers_map const &,
		Pose const & pose)
{
	resnum_ = protocols::rosetta_scripts::get_resnum( tag, pose );
	automatic_connection_ = tag->getOption< bool >( "automatic_connection", 1 );
	string const scorefxn_name( tag->getOption<string>( "scorefxn", "score12" ) );
	scorefxn_ = new ScoreFunction( *(data.get< ScoreFunction * >( "scorefxns", scorefxn_name) ));
	jump_num_ = tag->getOption<core::Size>( "jump_num", 1);
	std::string const final_filter_name( tag->getOption<std::string>( "final_filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator find_filter( filters.find( final_filter_name ));
	clash_check_ = tag->getOption<bool>("clash_check", 0 );
	explosion_ = tag->getOption<core::Size>( "explosion", 0);
	bool const filter_found( find_filter != filters.end() );
	if( filter_found )
		final_filter_ = find_filter->second->clone();
	else {
		if( final_filter_name != "true_filter" ){
			TR<<"***WARNING WARNING! Filter defined for TryRotamers not found in filter_list!!!! Defaulting to truefilter***"<<std::endl;
			runtime_assert( filter_found );
		}
		else
			final_filter_ = new protocols::filters::TrueFilter;
	}

	if( tag->hasOption( "shove" ) ){
  	std::string const shove_val( tag->getOption< std::string >( "shove" ) );
  	utility::vector1< std::string > const shove_keys( utility::string_split( shove_val, ',' ) );
  	foreach( std::string const key, shove_keys ){
			core::Size const resnum( protocols::rosetta_scripts::parse_resnum( key, pose ) );
			shove_residues_.push_back( resnum);
			TR<<"Using shove atomtype for "<<key<<'\n';
		}
	}

	TR<<"TryRotamers was instantiated using scorefxn="<<scorefxn_name<<", jump_number="<<jump_num_<< ", clash_check=" << clash_check_ << ", and explosion=" << explosion_ << std::endl;

}

} //movers
} //protein_interface_design
} //protocols

