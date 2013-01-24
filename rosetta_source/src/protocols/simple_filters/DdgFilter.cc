// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/DdgFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/DdgFilterCreator.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/elscripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/format.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.DdgFilter" );

const core::Real DdgFilter::DEFAULT_TRANSLATION_DISTANCE = 100.0;

protocols::filters::FilterOP
DdgFilterCreator::create_filter() const { return new DdgFilter; }

std::string
DdgFilterCreator::keyname() const { return "Ddg"; }


DdgFilter::DdgFilter() :
	filters::Filter( "Ddg" ),
	ddg_threshold_( -15.0 ),
	scorefxn_( NULL ),
	rb_jump_( 1 ),
	use_custom_task_(false),
	repack_bound_(true),
	relax_bound_(false),
	repeats_( 1 ),
	symmetry_(false),
	repack_( true ),
	relax_mover_( NULL ),
	pb_enabled_(false),
	translate_by_(DEFAULT_TRANSLATION_DISTANCE)
{
	scorename_ = "ddg";
}


DdgFilter::DdgFilter( core::Real const ddg_threshold, 
											core::scoring::ScoreFunctionCOP scorefxn, 
											core::Size const rb_jump/*=1*/, 
											core::Size const repeats/*=1*/, 
											bool const symmetry /*=false*/ ) :
	Filter("Ddg" ),
	ddg_threshold_(ddg_threshold),
	scorefxn_(scorefxn->clone()),
	rb_jump_(rb_jump),
	use_custom_task_( false ),
	repack_bound_( true ),
	relax_bound_( false ),
	repeats_(repeats),
	symmetry_(symmetry),
	repack_( true ),
	relax_mover_( NULL ),
	pb_enabled_(false),
	translate_by_(DEFAULT_TRANSLATION_DISTANCE)
{
	// Determine if this PB enabled.
	if( scorefxn_->get_weight(core::scoring::PB_elec) != 0.) {
		// Set this to PB enabled
		pb_enabled_ = true;
		TR << "PB enabled" << std::endl;
	}
	else{
		pb_enabled_ = false;
	}
}

DdgFilter::~DdgFilter() {}

filters::FilterOP DdgFilter::clone() const {
	return new DdgFilter( *this );
}
filters::FilterOP DdgFilter::fresh_instance() const{
	return new DdgFilter();
}

void
DdgFilter::repack( bool const repack )
{
	repack_ = repack;
}

bool
DdgFilter::repack() const
{
	return repack_;
}

void
DdgFilter::parse_my_tag( utility::tag::TagPtr const tag, 
												 moves::DataMap & data, 
												 filters::Filters_map const & , 
												 moves::Movers_map const & movers, 
												 core::pose::Pose const & )
{
	using namespace core::scoring;

	std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn", "score12" ) );
	scorefxn_ = new ScoreFunction( *(data.get< ScoreFunction * >( "scorefxns", scorefxn_name )) );
	ddg_threshold_ = tag->getOption<core::Real>( "threshold", -15 );
	rb_jump_ = tag->getOption< core::Size >( "jump", 1 );
	repeats( tag->getOption< core::Size >( "repeats", 1 ) );
	repack( tag->getOption< bool >( "repack", 1 ) );
	symmetry_ = tag->getOption<bool>( "symmetry", 0 );
	use_custom_task( tag->hasOption("task_operations") );
  task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	repack_bound( tag->getOption<bool>( "repack_bound", 1 ) );
	relax_bound( tag->getOption<bool>( "relax_bound", 0 ) );
	translate_by_ = tag->getOption<int>( "translate_by", DEFAULT_TRANSLATION_DISTANCE );

	if( tag->hasOption( "relax_mover" ) )
		relax_mover( protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "relax_mover" ), movers ) );

	if(tag->hasOption("chain_num"))
	{
		chain_ids_ = utility::string_split(tag->getOption<std::string>("chain_num"),',',core::Size());
	}

	if( repeats() > 1 && !repack() )
		throw utility::excn::EXCN_RosettaScriptsOption( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );

	if ( symmetry_ )
		TR<<"ddg filter with threshold "<< ddg_threshold_<<" repeats="<<repeats()<<" and scorefxn "<<scorefxn_name<<" with symmetry " <<std::endl;
	else
		TR<<"ddg filter with threshold "<< ddg_threshold_<<" repeats="<<repeats()<<" and scorefxn "<<scorefxn_name<<" over jump "<<rb_jump_<<" and repack "<<repack()<<std::endl;

	// Determine if this PB enabled.
	if( scorefxn_->get_weight(core::scoring::PB_elec) != 0.) {
		// Set this to PB enabled
		pb_enabled_ = true;
		TR << "PB enabled.  Translation distance = " << translate_by_ << " A" << std::endl;
		if( translate_by_ > DEFAULT_TRANSLATION_DISTANCE ) {
			TR.Warning << "Translation distance may be too large for PB-enabled scoring.  Consider 100 (default) if you run out of memory: " << translate_by_ << std::endl;
			TR.Warning.flush();
		}
	}
	else{
		pb_enabled_ = false;
	}
	TR.flush();
}

void DdgFilter::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & score_fxns,
				utility::lua::LuaObject const & /*tasks*/ ) {
	using namespace core::scoring;
	if( def["scorename"] )
		scorename_ = def["scorename"].to<std::string>();
	if( def["scorefxn"] ) {
		scorefxn_ = protocols::elscripts::parse_scoredef( def["scorefxn"], score_fxns );
	} else {
		scorefxn_ = score_fxns["score12"].to<ScoreFunctionSP>()->clone();
	}
	ddg_threshold_ = def["threshold"] ? def["threshold"].to<core::Real>() : -15;
	rb_jump_ = def["jump"] ? def["jump"].to<core::Size>() : 1;
	repeats( def["repeats"] ? def["repeats"].to<core::Size>() : 1 );
	repack( def["repack"] ? def["repack"].to<bool>() : true );
	symmetry_ = def["symmetry"] ? def["symmetry"].to<bool>() : false;
	repack_bound_ = def["repack_bound"] ? def["repack_bound"].to<bool>() : true;
	relax_bound_ = def["relax_bound"] ? def["relax_bound"].to<bool>() : false;
	translate_by_ = def["translate_by"] ? def["translate_by"].to<int>() : DEFAULT_TRANSLATION_DISTANCE;
	// ignoring relax_mover option
	if( def["chain_num"] ) {
		chain_ids_.clear();
		for (utility::lua::LuaIterator i=def["chain_num"].begin(), end; i != end; ++i) {
			chain_ids_.push_back( (*i).to<core::Size>() );
		}
	}

	if( repeats() > 1 && !repack() )
		utility_exit_with_message( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );

	if ( symmetry_ )
		TR<<"ddg filter with threshold "<< ddg_threshold_<<" repeats="<<repeats()<<" with symmetry " <<std::endl;
	else
		TR<<"ddg filter with threshold "<< ddg_threshold_<<" repeats="<<repeats()<<" over jump "<<rb_jump_<<" and repack "<<repack()<<std::endl;
}

bool
DdgFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const pose_ddg( compute( pose ) );
	TR<<"ddg is "<<pose_ddg<<" ";
	if( pose_ddg <= ddg_threshold_ ) {
		TR<<"passing"<<std::endl;
		return true;
	}
	TR<<"failing"<<std::endl;
	TR.flush();
	return false;
}

void
DdgFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const pose_ddg( compute( pose ) );
	out<<"ddg "<<pose_ddg<<'\n';
}

core::Real
DdgFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const pose_ddg( compute( pose ) );
	return( pose_ddg );
}

core::Size
DdgFilter::repeats() const
{
	return( repeats_ );
}

void
DdgFilter::repeats( core::Size const repeats )
{
	repeats_ = repeats;
}

core::Real
DdgFilter::compute( core::pose::Pose const & pose_in ) const {
	core::pose::Pose pose(pose_in);
	if( repack() ){
		protocols::simple_moves::ddG ddg( scorefxn_, rb_jump_, chain_ids_, symmetry_ );
		if ( use_custom_task() ) {
			ddg.use_custom_task( use_custom_task() );
			ddg.task_factory( task_factory() );
		}
		if ( repack_bound() ) {
			ddg.repack_bound( repack_bound() );
		}
		if ( relax_bound() ) {
			ddg.relax_bound( relax_bound() );
		}
		ddg.relax_mover( relax_mover() );
		core::Real average( 0.0 );
		for( core::Size i = 1; i<=repeats_; ++i ){
			ddg.calculate( pose );
			average += ddg.sum_ddG();
			ddg.report_ddG( TR );
		}
		return average / (core::Real)repeats_;
	} else {
		if( repeats() > 1 && !repack() )
			utility_exit_with_message( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );
		using namespace protocols::moves;

		simple_filters::ScoreTypeFilter const stf( scorefxn_, core::scoring::total_score, 10000/*threshold*/ );
		core::pose::Pose split_pose( pose );
		if(chain_ids_.size() > 0)
		{
			//We want to translate each chain the same direction, though it doesnt matter much which one
			core::Vector translation_axis(1,0,0); 
			for(utility::vector1<core::Size>::const_iterator chain_it = chain_ids_.begin(); chain_it != chain_ids_.end();++chain_it)
			{
				core::Size current_chain_id = *chain_it;
				core::Size current_jump_id = core::pose::get_jump_id_from_chain_id(current_chain_id,split_pose);
				rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( split_pose, current_jump_id) );
				translate->step_size( translate_by_ );
				translate->trans_axis(translation_axis);
				translate->apply( split_pose );
			}
		}else
		{
			rigid::RigidBodyTransMoverOP translate( new rigid::RigidBodyTransMover( split_pose, rb_jump_ ) );
			translate->step_size( translate_by_ );
			translate->apply( split_pose );
		}

		core::Real const bound_energy( stf.compute( pose ));
		core::Real const unbound_energy( stf.compute( split_pose ));
		core::Real const dG( bound_energy - unbound_energy );
		return( dG );
	}
}

void
DdgFilter::relax_mover( protocols::moves::MoverOP m ){
	relax_mover_ = m;
}

protocols::moves::MoverOP
DdgFilter::relax_mover() const{ return relax_mover_; }


}
}
