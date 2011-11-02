// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlaceOnLoop.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PlaceOnLoop.hh>
#include <protocols/protein_interface_design/movers/PlaceOnLoopCreator.hh>


#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <utility/tag/Tag.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
// Auto-header: duplicate removed #include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/toolbox/task_operations/RestrictChainToRepackingOperation.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/DataMap.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/KinematicMover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/LoopMover_KIC.hh>
#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
// AUTO-REMOVED #include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/loops/KinematicWrapper.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/loops/loops_main.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <numeric/random/random.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>


static numeric::random::RandomGenerator RG( 14071789 );

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.PlaceOnLoop" );

std::string
PlaceOnLoopCreator::keyname() const
{
	return PlaceOnLoopCreator::mover_name();
}

protocols::moves::MoverOP
PlaceOnLoopCreator::create_mover() const {
	return new PlaceOnLoop;
}

std::string
PlaceOnLoopCreator::mover_name()
{
	return "PlaceOnLoop";
}

PlaceOnLoop::PlaceOnLoop() :
	DesignRepackMover( PlaceOnLoopCreator::mover_name() ),
	loop_begin_( 0 ), loop_end_( 0 ),
	hires_scorefxn_( NULL ), lores_scorefxn_( NULL ),
	chain_closing_attempts_( 100 ), host_chain_( 2 ), stub_set_( NULL ), minimize_toward_stub_( true )
{
	delta_length_.clear();
	delta_length_.push_back( 0 );
	kinematic_mover_ = new protocols::moves::KinematicMover;
	set_kinematic_defaults();
}

PlaceOnLoop::~PlaceOnLoop() {}

MoverOP PlaceOnLoop::fresh_instance() const {
	return MoverOP( new PlaceOnLoop );
}

protocols::moves::MoverOP
PlaceOnLoop::clone() const {
  return( MoverOP( new PlaceOnLoop( *this )));
}

bool
PlaceOnLoop::minimize_toward_stub( core::pose::Pose & pose ) const
{
	using namespace protocols::loops;

	core::kinematics::FoldTree const saved_ft( pose.fold_tree() );
	core::kinematics::FoldTree f_new;
	Loop loop( loop_begin_, curr_loop_end_ );
	loop.choose_cutpoint( pose );
	Loops loops;
	loops.push_back( loop );
	loops::fold_tree_from_loops( pose, loops, f_new, true /* include terminal cutpoints */);
	pose.fold_tree( f_new );
	pose.update_residue_neighbors();
	protocols::loops::KinematicWrapper kinwrap( kinematic_mover_, loop, chain_closing_attempts_ );
	kinwrap.apply( pose );
	if( kinwrap.get_last_move_status() != MS_SUCCESS ){
		TR<<"Kinematic loop closure failed to close loop."<<std::endl;
		pose.fold_tree( saved_ft );
		return( false );
	}
	LoopMover_Refine_KIC refine( loops, hires_scorefxn_ );
	using namespace core::pack::task;
	using namespace protocols::toolbox::task_operations;
	core::pack::task::TaskFactoryOP tf = new TaskFactory;
	tf->push_back( new ProteinInterfaceDesignOperation );
	tf->push_back( new RestrictChainToRepackingOperation( 1 ) );
	tf->push_back( new RestrictChainToRepackingOperation( 2 ) );
	refine.set_task_factory( tf );
	refine.apply( pose );

	pose.fold_tree( saved_ft );
	pose.update_residue_neighbors();
	return( true );
}

/// @details bb_csts only affect loop. Transforms loop to poly-ala
void
PlaceOnLoop::add_bb_csts_to_loop( core::pose::Pose & pose ) const
{
	using namespace core::pack::task;
	using namespace core::chemical;

	if( stub_set_() == NULL ) return;
	PackerTaskOP ptask = TaskFactory::create_packer_task( pose );
	for( core::Size i( 1 ); i<=pose.total_residue(); ++i ){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		if( i<=loop_begin_ || i>=curr_loop_end_ )
			ptask->nonconst_residue_task( i ).prevent_repacking();
		if( (pose.residue( i ).aa() == aa_gly || pose.residue( i ).aa() == aa_pro ) )
			ptask->nonconst_residue_task( i ).prevent_repacking();
	}
	stub_set_->pair_with_scaffold( pose, host_chain_, new protocols::filters::TrueFilter );
  core::Size fixed_res(1);
  if( host_chain_ == 1 ) fixed_res = pose.total_residue();
	core::id::AtomID const fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );
	stub_set_->add_hotspot_constraints_to_pose( pose, fixed_atom_id, ptask, stub_set_, 0.7/*cbforce*/, 0/*worst allowed stub bonus*/, false/*apply self energies*/, 10.0/*bump cutoff*/, true/*apply ambiguous constraints*/ );
}

void
PlaceOnLoop::ala_pose_loop( core::pose::Pose & pose ) const
{
	using namespace core::pack::task;
	using namespace core::chemical;

	utility::vector1< bool > ala_only( num_canonical_aas, false );
	ala_only[ aa_ala ] = true;
	PackerTaskOP ptask = TaskFactory::create_packer_task( pose );
	for( core::Size i( 1 ); i<=pose.total_residue(); ++i ){
		if( i>=loop_begin_ && i<=curr_loop_end_ )
			ptask->nonconst_residue_task( i ).restrict_absent_canonical_aas( ala_only );
		else
			ptask->nonconst_residue_task( i ).prevent_repacking();
	}
	core::pack::pack_rotamers( pose, *hires_scorefxn_, ptask );
}

/// @details Chooses one user-defined loop length at random, changes the loop accordingly,
/// and tries to close the loop with that change.
bool
PlaceOnLoop::loop_length( core::pose::Pose & pose )
{
	core::scoring::constraints::ConstraintCOPs saved_bb_constraints = protocols::hotspot_hashing::remove_hotspot_constraints_from_pose( pose );
	Pose const saved_pose( pose );
	curr_loop_end_ = loop_end_;
	std::random_shuffle( delta_length_.begin(), delta_length_.end() );
	bool loop_closed( false );
	int const delta( *delta_length_.begin() );
	TR<<"changing loop length by "<<delta<<std::endl;
	if( delta < 0 ){
		for( int del(-1); del>=delta; --del ){
			pose.delete_polymer_residue( loop_begin_ + 1 );
			curr_loop_end_--;
		}
	}
	else if( delta > 0 ){
		using namespace core::chemical;
		using namespace core::conformation;

		ResidueTypeSet const & residue_set( pose.residue( 1 ).residue_type_set() ); // residuetypeset is noncopyable
		ResidueCOP new_res = ResidueFactory::create_residue( residue_set.name_map( name_from_aa( aa_from_oneletter_code( 'A' ) ) ) );
		for( core::Size leng(1); leng<=(core::Size) delta; ++leng ){
			pose.conformation().safely_append_polymer_residue_after_seqpos( *new_res, loop_begin_ + 1, true/*build_ideal_geometry*/ );
			curr_loop_end_++;
		}
	}

	pose.update_residue_neighbors();
	ala_pose_loop( pose );
	add_bb_csts_to_loop( pose );
	if( minimize_toward_stub_ )
		loop_closed = minimize_toward_stub( pose );
	else
		loop_closed = position_stub( pose );

	return( loop_closed );
}

bool
PlaceOnLoop::position_stub( core::pose::Pose & ) const
{
	utility_exit_with_message( "position_stub is not yet implemented." );
	return( false );
}

/// @details set the defaults for the kinmover. Taken from Jacob's LoopRemodel
void
PlaceOnLoop::set_kinematic_defaults()
{
	runtime_assert( kinematic_mover_ );
	kinematic_mover_->set_idealize_loop_first( true );
	kinematic_mover_->set_temperature( 0.6/*mc_kt*/ );
	kinematic_mover_->set_sfxn(hires_scorefxn_);
	kinematic_mover_->set_vary_bondangles( true );
	kinematic_mover_->set_sample_nonpivot_torsions( true );
	kinematic_mover_->set_rama_check( true );
}

void
PlaceOnLoop::apply( Pose & pose )
{
	bool const success( loop_length( pose ) );
	if( success ) set_last_move_status( MS_SUCCESS );
	else set_last_move_status( FAIL_RETRY );
}

std::string
PlaceOnLoop::get_name() const {
	return PlaceOnLoopCreator::mover_name();
}

void
PlaceOnLoop::parse_my_tag( TagPtr const tag, DataMap &data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	host_chain_ = tag->getOption< core::Size >( "host_chain", 2 );
	loop_begin_ = tag->getOption< core::Size >( "loop_begin" );
	loop_end_ = tag->getOption< core::Size >( "loop_end" );
	runtime_assert( loop_begin_ < loop_end_ );
	runtime_assert( loop_begin_ > 1 );
	core::Size const chain_begin( pose.conformation().chain_begin( host_chain_ ) );
	core::Size const chain_end( pose.conformation().chain_end( host_chain_ ) );
	runtime_assert( chain_begin < loop_begin_ && loop_begin_ < chain_end );
	runtime_assert( chain_begin < loop_end_ && loop_end_ < chain_end );
	minimize_toward_stub_ = tag->getOption< bool >( "minimize_toward_stub", 1 );
	if( tag->hasOption( "stubfile" ) ){
		std::string const stub_fname = tag->getOption< std::string >( "stubfile" );
		stub_set_ = new protocols::hotspot_hashing::HotspotStubSet;
		stub_set_->read_data( stub_fname );
	}

	std::string const score_hi_name = tag->getOption< std::string >( "score_high", "score12" );
	std::string const score_lo_name = tag->getOption< std::string >( "score_low", "score4L" );
	hires_scorefxn_ = data.get< core::scoring::ScoreFunction * >( "scorefxns", score_hi_name );
	lores_scorefxn_ = data.get< core::scoring::ScoreFunction * >( "scorefxns", score_lo_name );
	chain_closing_attempts_ = tag->getOption< core::Size >( "closing_attempts", 100 );

  typedef utility::vector1< std::string > StringVec;
	std::string shorten_by(""), lengthen_by("");
	if( tag->hasOption( "shorten_by" ) )
		shorten_by = tag->getOption< std::string >( "shorten_by" );
	if( tag->hasOption( "lengthen_by" ) )
		lengthen_by = tag->getOption< std::string >( "lengthen_by" );
  StringVec const shorten_by_keys( utility::string_split( shorten_by, ',' ) );
  StringVec const lengthen_by_keys( utility::string_split( lengthen_by, ',' ) );
	core::Size const loop_length( loop_end_ - loop_begin_ + 1 );
	foreach( std::string const shorten, shorten_by_keys ){
		if( shorten == "" ) continue;
		int const shorten_i( -1 * atoi( shorten.c_str() ) );
		runtime_assert( ( core::Size ) abs( shorten_i ) < loop_length );
		delta_length_.push_back( shorten_i );
	}

	foreach( std::string const lengthen, lengthen_by_keys ){
		if( lengthen == "" ) continue;
		delta_length_.push_back( atoi( lengthen.c_str() ) );
	}

	sort( delta_length_.begin(), delta_length_.end() );
	runtime_assert( loop_end_ - loop_begin_ + 1 + *delta_length_.begin() >= 3 );
	unique( delta_length_.begin(), delta_length_.end() );
	TR<<"PlaceOnLoop mover defined with kinematic mover ";
	TR<<" will change the loop by these values: ";
	foreach( int const i, delta_length_ ) TR<<i<<", ";
	TR<<std::endl;
}

} //movers
} //protein_interface_design
} //protocols
