// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/GenericMonteCarloMover.cc
/// @brief perform a given mover and sample structures by MonteCarlo
/// @detailed The score evaluation of pose during MC after applying mover is done by
/// either FilterOP that can do report_sm() or ScoreFunctionOP.
/// By setting sample_type_ to high, you can also sample the pose that have higher score.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


// Unit Headers
#include <protocols/moves/GenericMonteCarloMover.hh>
#include <protocols/moves/GenericMonteCarloMoverCreator.hh>
#include <protocols/moves/DataMapObj.hh>
#include <protocols/moves/DataMap.hh>

// C/C++ headers
#include <iostream>
#include <iterator>
#include <string>

// External headers
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#define foreach BOOST_FOREACH
#include <algorithm>

// Utility headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/filters/Filter.hh>

// Package Headers
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <fstream>
#include <utility/io/izstream.hh>
#include <sstream>
#include <core/pose/util.hh>
#include <protocols/simple_filters/OperatorFilter.hh>
#include <protocols/filters/BasicFilters.hh>

static basic::Tracer TR("protocols.moves.GenericMonteCarloMover");
static basic::Tracer TR_energies("protocols.moves.GenericMonteCarloMover.individual_energies");
static numeric::random::RandomGenerator mc_RG(61452); // <- Magic number, do not change it!!!

using namespace core;

namespace protocols {
namespace moves {

using namespace ObjexxFCL::fmt;

std::string
GenericMonteCarloMoverCreator::keyname() const
{
  return GenericMonteCarloMoverCreator::mover_name();
}

protocols::moves::MoverOP
GenericMonteCarloMoverCreator::create_mover() const {
  return new GenericMonteCarloMover;
}

std::string
GenericMonteCarloMoverCreator::mover_name()
{
  return "GenericMonteCarlo";
}

/// @brief default constructor
GenericMonteCarloMover::GenericMonteCarloMover():
	Mover("GenericMonteCarlo"),
	maxtrials_( 10 ),
	task_scaling_( 5 ),
	mover_( NULL ),
	scorefxn_( NULL ),
	temperature_( 0.0 ),
	sample_type_( "low" ),
	drift_( true ),
	preapply_( true ),
	recover_low_( true ),
	rank_by_filter_( 1 ),
	boltz_rank_( false ),
	last_accepted_pose_( NULL ),
	lowest_score_pose_( NULL ),
	stopping_condition_( NULL ),
	mover_stopping_condition_( NULL ),
	adaptive_movers_( false ),
	adaptation_period_( 0 ),
	saved_accept_file_name_( "" ),
	saved_trial_number_file_( "" ),
	mover_tag_( NULL )
{
  initialize();
}


/// @brief value constructor without a score function

GenericMonteCarloMover::GenericMonteCarloMover(
  Size const maxtrials,
  Size const task_scaling,
  MoverOP const & mover,
  Real const temperature,
  String const sample_type,
  bool const drift ) :
	Super("GenericMonteCarlo"),
	maxtrials_( maxtrials ),
	task_scaling_( task_scaling ),
	mover_( mover ),
	temperature_( temperature ),
	sample_type_( sample_type ),
	drift_( drift ),
	preapply_( true ),
	recover_low_( true ),
	rank_by_filter_(1),
	boltz_rank_( false ),
	last_accepted_pose_( NULL ),
	lowest_score_pose_( NULL ),
	saved_accept_file_name_( "" ),
	mover_tag_( NULL )
{
  initialize();
}


/// @brief value constructor with a TaskFactory
GenericMonteCarloMover::GenericMonteCarloMover(
  Size const maxtrials,
  Size const task_scaling,
  MoverOP const & mover,
	TaskFactoryOP factory_in,
  Real const temperature,
  String const sample_type,
  bool const drift ) :
  Super("GenericMonteCarlo"),
  maxtrials_( maxtrials ),
	task_scaling_( task_scaling ),
  mover_( mover ),
	task_( NULL ),
	factory_ (factory_in),
  temperature_( temperature ),
  sample_type_( sample_type ),
  drift_( drift ),
	preapply_( true ),
  recover_low_( true ),
	rank_by_filter_(1),
	boltz_rank_( false ),
	saved_accept_file_name_( "" ),
	mover_tag_( NULL )
{
  initialize();
}


/// @brief destructor
GenericMonteCarloMover::~GenericMonteCarloMover(){}

/// @brief clone this object
GenericMonteCarloMover::MoverOP
GenericMonteCarloMover::clone() const
{
  return new GenericMonteCarloMover( *this );
}

void GenericMonteCarloMover::task_factory( core::pack::task::TaskFactoryOP tf ) { factory_ = tf; }

/// @brief create this type of object
GenericMonteCarloMover::MoverOP
GenericMonteCarloMover::fresh_instance() const
{
  return new GenericMonteCarloMover();
}

/// @brief initialize
void
GenericMonteCarloMover::initialize()
{
  last_accepted_scores_.clear();
  if( sample_type_ == "high" ){
    flip_sign_ = -1;
  }else if( sample_type_ == "low" ){
    flip_sign_ = 1;
  }else{
    TR << "WARNING: the sample type, " << sample_type_ << ", is not defined." << std::endl;
    runtime_assert( false );
  }
  trial_counter_ = 0;
  accept_counter_ = 0;
  energy_gap_counter_ = 0.0;

  // trigger initialization
  next_trigger_id_ = 1;
}

/// @brief return the last accepted pose
GenericMonteCarloMover::PoseOP
GenericMonteCarloMover::last_accepted_pose() const
{
  return last_accepted_pose_;
}

/// @brief return the last accepted score
GenericMonteCarloMover::Real
GenericMonteCarloMover::last_accepted_score() const
{
  return last_accepted_score_;
}

/// @brief return the lowest score pose
GenericMonteCarloMover::PoseOP
GenericMonteCarloMover::lowest_score_pose() const
{
  return lowest_score_pose_;
}

/// @brief return the lowest score
GenericMonteCarloMover::Real
GenericMonteCarloMover::lowest_score() const
{
  return lowest_score_;
}

/// @brief return the lowest score
GenericMonteCarloMover::Real
GenericMonteCarloMover::current_score() const
{
  return current_score_;
}

/// @brief return mc_accepted
MCA
GenericMonteCarloMover::mc_accpeted() const
{
  return mc_accepted_;
}

/// @brief set max trials of monte carlo iterations
void
GenericMonteCarloMover::set_maxtrials( Size const ntrial )
{
  maxtrials_ = ntrial;
}

/// @brief set task multiplier to calculate trials from task
void
GenericMonteCarloMover::set_task_scaling( Size const scaling )
{
  task_scaling_ = scaling;
}

/// @brief set mover
void
GenericMonteCarloMover::set_mover( MoverOP const & mover )
{
  mover_ = mover;
}

/// @brief set filter
/// Pose is evaluated by FilterOP which can do report_sm() or ScoreFunctionOP during MC trials
/// You can choose either way FilterOP or ScoreFunction.
void
GenericMonteCarloMover::add_filter( FilterOP filter, bool const adaptive, Real const temp, String const sample_type, bool rank_by)
{
  filters_.push_back( filter );
	if(rank_by) {
		rank_by_filter_ = filters_.size();
	}
  adaptive_.push_back( adaptive );
  temperatures_.push_back( temp );
  sample_types_.push_back( sample_type );
  last_accepted_scores_.assign( filters_.size(), 100000 );
	num_rejections_.push_back(0);
  scorefxn_ = NULL;
}

/// @brief set scorefxn
/// Pose is evaluated by FilterOP which can do report_sm() or ScoreFunctionOP during MC trials
/// You can choose either way FilterOP or ScoreFunction.
void
GenericMonteCarloMover::set_scorefxn( ScoreFunctionOP const & sfxn )
{
  scorefxn_ = sfxn;
  filters_.clear();
}

/// @brief set temperatrue
void
GenericMonteCarloMover::set_temperature( Real const temp )
{
  temperature_ = temp;
}

/// @brief set sample type, high or low
/// when sample_type == high, sample pose which have higher value of scorey
/// when sample_type == low, sample pose which have lower value of score
void
GenericMonteCarloMover::set_sampletype( String const & type )
{
  if( sample_type_ != "high" && sample_type_ != "low" ){
    TR << "WARNING !! the sample type, " << type << ", is not defined." << std::endl;
    runtime_assert( false );
  }
  sample_type_ = type;
}

/// @brief if drift=false, the pose is set back to the initial pose
/// Of course, this is not MC sampling.
void
GenericMonteCarloMover::set_drift( bool const drift ){
  drift_ = drift;
}

/// @brief if preapply=true, auto-accept the first application of the submover,
/// ignoring boltzman criteria.
void
GenericMonteCarloMover::set_preapply( bool const preapply ) {
	preapply_ = preapply;
}

/// @brief if recover_low=true, after apply() the returned
/// is the lowest energy structure, rather than the last accepted structure.
void
GenericMonteCarloMover::set_recover_low( bool const recover_low ){
  recover_low_ = recover_low;
}

/// @brief if boltz_rank=true, rank structures by the temperature-weighted
/// sum of scores, rather than a single filter
void
GenericMonteCarloMover::set_boltz_rank( bool const boltz_rank ){
  boltz_rank_ = boltz_rank;
}

/// @brief show scores of last_accepted_score and lowest_score
void
GenericMonteCarloMover::show_scores( std::ostream & out ) const
{
  if( sample_type_ == "high" ){
    out << "Higher score sampled: trial=" << I( 5, trial_counter_ ) << ", score(current/last_accepted/best)=";
  }else{
    out << "Lower score sampled: trial=" << I( 5, trial_counter_ ) << ", score(current/last_accepted/best)=";
  }
  out << F( 9, 3, flip_sign_*current_score() ) << '/'
      << F( 9, 3, flip_sign_*last_accepted_score() )<< '/'
      << F( 9, 3, flip_sign_*lowest_score() ) << std::endl;
}

/// @brief show counters
void
GenericMonteCarloMover::show_counters( std::ostream & out ) const
{
  String const & mover( mover_->get_name() );
  String evaluation;
  if( scorefxn_ ){
    evaluation = "ScoreXXX"; // should we have name for ScoreFunction ?
  }
  int  const ntrials( trial_counter_ );
  int  const accepts( accept_counter_ );
  Real const energy_gap( energy_gap_counter_ );
  if( accepts > 0 ){
    out << "mover=" << LJ( 16, mover ) << " Score_eval=" << LJ( 16, evaluation ) <<
      " trials= " << I( 5, ntrials ) << "; " <<
      " accepts= " << F( 6, 3, Real( accepts )/ntrials ) << "; " <<
      " energy_gap/trial= " << F( 8, 3, Real( flip_sign_ * energy_gap ) / ntrials ) << std::endl;
  }else{
    out << "mover=" << A( 16, mover ) << " Score_eval=" << A( 16, evaluation )
        << " trials= " << I( 6, ntrials ) <<  " NO ACCEPTS." << std::endl;
  }
	if(num_rejections_.size()) {
		out << "Number of rejections per filter: ";
		for( core::Size ii(1); ii <= num_rejections_.size(); ++ii) {
			out << num_rejections_[ii] << "   ";
		}
		out << std::endl;
	}
}

/// @brief return the simulation state to the lowest energy structure we've seen
void
GenericMonteCarloMover::recover_low( Pose & pose )
{
  if( lowest_score_pose_ ) {
    pose = *lowest_score_pose_;
    *last_accepted_pose_ = *lowest_score_pose_;
  }else{
    // Case of that lowest_score_pose_ was never updated in MonteCarlo
    TR << " No lowest_score_pose_ exists, probably because all movers or filters failed. " << std::endl;
  }
}


/// @brief reset this GenericMonteCarloMover
void
GenericMonteCarloMover::reset( Pose & pose )
{
  if( filters_.size() == 0 ) {
    lowest_score_ = scoring( pose );
	} else {
    last_accepted_scores_.clear();
    for( Size index = 1; index <= filters_.size(); ++index ){
      protocols::filters::FilterCOP filter( filters_[ index ] );
      Real const flip( sample_types_[ index ] == "high" ? -1 : 1 );
      last_accepted_scores_.push_back( flip * filter->report_sm( pose ) );
    }
  }// fi filters_.size()

  lowest_score_pose_ = new Pose( pose );
  last_accepted_pose_ = new Pose( pose );
  last_accepted_score_ = lowest_score_;

  trial_counter_ = 0;
  accept_counter_ = 0;
  energy_gap_counter_ = 0;
	if ( num_rejections_.size() ) {
		num_rejections_.assign( num_rejections_.size(), 0 );
	}
  TR << "Initialization done " << std::endl;
}

/// @brief score pose based on filter or scorefxn
GenericMonteCarloMover::Real
GenericMonteCarloMover::scoring( Pose & pose )
{
  Real score( 0.0 );
  if( scorefxn_ ){
    score = flip_sign_ * (*scorefxn_)( pose );
  }
  return score;
}
GenericMonteCarloMover::Size
GenericMonteCarloMover::num_designable( Pose & pose, PackerTaskOP & task )
{
	Size number_designable = 0;
	TR << "Designable Residues: ";
	for(Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i){
		if(task->design_residue(i)) {
			++number_designable;
			TR << i << ", ";
		}
	}
	TR << std::endl << "Calculated number designable residues" << ": ";
	TR << number_designable << std::endl;
	return number_designable;
}

// @brief
bool
GenericMonteCarloMover::boltzmann( Pose & pose )
{
  ++trial_counter_;
  TR.Debug <<"filters.size() "<<filters_.size()<<std::endl;
  if( filters_.size() ){
    runtime_assert( filters_.size() == adaptive_.size() && filters_.size() == temperatures_.size() && filters_.size() == sample_types_.size() &&  filters_.size() == num_rejections_.size());
    bool accept( false );
    utility::vector1< Real > provisional_scores;
    provisional_scores.clear();
    Real ranking_score( 0.0 );
    for( core::Size index( 1 ); index <= filters_.size(); ++index ){
      TR.Debug <<"Filter #"<<index<<std::endl;
      protocols::filters::FilterCOP filter( filters_[ index ] );
      bool const adaptive( adaptive_[ index ] );
      Real const temp( temperatures_[ index ] );
      Real const flip( sample_types_[ index ] == "high" ? -1 : 1 );
			core::Real const filter_val( filter->report_sm( pose ));
			TR<<"Filter "<<index<<" reports "<<filter_val<<std::endl;

      provisional_scores.push_back( flip * filter_val );
      if( index == rank_by_filter_ ) {
        ranking_score = provisional_scores[ rank_by_filter_ ];
			}
      Real const boltz_factor = ( last_accepted_scores_[ index ] - provisional_scores[ index ] ) / temp;
      TR_energies.Debug <<"energy index, last_accepted_score, current_score "<<index<<" "<< last_accepted_scores_[ index ]<<" "<<provisional_scores[ index ]<<std::endl;
      TR.Debug <<"Current, best, boltz "<<provisional_scores[ index ]<<" "<<last_accepted_scores_[ index ]<<" "<<boltz_factor<<std::endl;
      if( !adaptive ) { // return the starting score
        provisional_scores[ index ] = last_accepted_scores_[ index ];
			}
      Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
      Real const random_num( mc_RG.uniform() );
      bool const reject_filter( provisional_scores[ index ] > last_accepted_scores_[ index ] && random_num >= probability );
      if( reject_filter ){
        accept = false;
				++num_rejections_[index];
        break;
      }
      if( !reject_filter ) {
        accept = true;
			}
    }//for index
    if( accept ){
      TR<<"Accept"<<std::endl;
      mc_accepted_ = MCA_accepted_thermally;
      copy( provisional_scores.begin(), provisional_scores.end(), last_accepted_scores_.begin() );
      ++accept_counter_;
			if(boltz_rank_) {
				ranking_score = 0.0;
				for(core::Size ii(1); ii <= filters_.size(); ++ii ){
					ranking_score += provisional_scores[ii] / temperatures_[ii];
				}
			}
      energy_gap_counter_ += ranking_score - last_accepted_score();
      last_accepted_score_ = ranking_score;
      *last_accepted_pose_ = pose;
      if( ranking_score <= lowest_score() ){
        *lowest_score_pose_ = pose;
        lowest_score_ = ranking_score;
        mc_accepted_ = MCA_accepted_score_beat_low; //3;
				if( saved_accept_file_name_ != "" ){
					TR<<"Dumping accepted file to disk as: "<<saved_accept_file_name_<<std::endl;
					pose.dump_pdb( saved_accept_file_name_ );
					if( mover_tag_() != NULL ){
						TR<<"Dumping accepted mover tag to disk"<<std::endl;
						std::ofstream f;
						std::string const fname( saved_accept_file_name_ + ".mover_tag" );
						f.open( fname.c_str(), std::ios::out );
						if( !f.good() )
							utility_exit_with_message( "Unable to open MC mover_tag file " + fname );
						f<<mover_tag_->obj;
						f.close();
						core::pose::add_comment( pose, user_defined_mover_name_, mover_tag_->obj ); /// adding comment to the pose to save the mover's tag since it's accepted
					}
				}
      }
      return( true );
    }// fi accept
    else{
      TR.Debug <<"Reject"<<std::endl;
      mc_accepted_ = MCA_rejected;
      return( false );
    }
  }//fi filters_.size()
  else{
    Real score = scoring( pose );
    current_score_ = score; // for debugging
    show_scores( TR.Debug );
    if ( score > last_accepted_score() ) {
      if( temperature_ < 1e-8 ){
        mc_accepted_ = MCA_rejected; // rejected
      }else{
        Real const boltz_factor = ( last_accepted_score() - score ) / temperature_;
        Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );
        if ( mc_RG.uniform() >= probability ) {
          mc_accepted_ = MCA_rejected; // rejected
        }else{
          mc_accepted_ = MCA_accepted_thermally; // accepted thermally
        }
      }
    }else{
      mc_accepted_ = MCA_accepted_score_beat_last; // accepted: energy is lower than last_accepted
    }

    if( mc_accepted_ >= 1 ){ // accepted
      ++accept_counter_;
      energy_gap_counter_ += score - last_accepted_score();
      *last_accepted_pose_ = pose;
      last_accepted_score_ = score;
      if ( score < lowest_score() ){  // energy is lower than last_accepted
        lowest_score_ = score;
        *lowest_score_pose_ = pose;
        mc_accepted_ = MCA_accepted_score_beat_low; //3;
      }
      return true;
    }else{ //rejected
      return false;
    }
  }
} // boltzmann

core::Size
GenericMonteCarloMover::load_trial_number_from_checkpoint() const{
	using namespace std;
	if( saved_trial_number_file_ == "" )
		return 1;
	ifstream f( saved_trial_number_file_.c_str(), ios::in );
	if( !f.good() )
		return 1;

	core::Size const begin = f.tellg();
	f.seekg( 0, ios::end );
	core::Size const end = f.tellg();
	if( end - begin == 0 )//file size == 0
		return 1;
	f.seekg( 0, ios::beg );// return to the beginning

	TR<<"Loading trial number from checkpoint"<<std::endl;
	std::string line;
	getline( f, line );
	std::istringstream line_stream( line );
	core::Size trial;
	line_stream >> trial;
	TR<<"Loaded trial number: "<<trial<<std::endl;
	if( mover_tag_() != NULL ){
		std::string const fname( saved_trial_number_file_ + ".mover_tag" );
		ifstream f_mover_tag( fname.c_str(), ios::in );
		if( f_mover_tag.good() ){
			f_mover_tag >> mover_tag_->obj;
			TR<<"Loaded mover_tag from checkpointing file: "<<mover_tag_->obj<<std::endl;
		}
		else
			TR<<"File containing mover_tag "<<fname<<" not found. Not loading movertag from checkpoint"<<std::endl;
	}

	return trial;
}

void
GenericMonteCarloMover::save_trial_number_to_checkpoint( core::Size const i ) const{
	if( saved_trial_number_file_ == "" )
		return;
	std::ofstream f;
	f.open( saved_trial_number_file_.c_str(), std::ios::out );
	if( !f.good() )
		utility_exit_with_message( "Unable to open MC checkpointing file " + saved_trial_number_file_ );

	f<<i;
	f.close();
}

/// @Brief
///comment
void
GenericMonteCarloMover::apply( Pose & pose )
{
  using protocols::moves::FAIL_DO_NOT_RETRY;
  using protocols::moves::FAIL_BAD_INPUT;
  using protocols::moves::FAIL_RETRY;

  if( !mover_ ){
    TR.Warning << "Mover is empty ! " << std::endl;
    return;
  }
  if( !filters_.size() && !scorefxn_ ){
    TR.Warning << "Both ScorefunctionOP and FilterOP are empty ! " << std::endl;
    return;
  }

	core::pack::task::PackerTaskOP task;

	if ( factory_ ) {
		task = factory_->create_task_and_apply_taskoperations( pose );
	} else {
		TR << "No task inputted" << std::endl; //LGN
	}

	if (factory_ ){
		number_designable_ = num_designable (pose, task);
		TR << "Input number of trials is " << maxtrials_ << std::endl;
		TR << "The task_scaling is: " << task_scaling_ << std::endl;
		maxtrials_ = task_scaling_ * number_designable_;
		TR << "Resetting number of trials based on your input task to: " << maxtrials_ << std::endl;
	}
	TR << "The number of trials for this run is: " << maxtrials_ << std::endl;

	bool const stop_at_start( ( mover_stopping_condition_() != NULL && mover_stopping_condition_->obj ) || stopping_condition()->apply( pose ) );
	if( stop_at_start ){
		TR<<"MC stopping condition met at the start, so failing without retrying "<<std::endl;
		set_last_move_status( FAIL_DO_NOT_RETRY );
		return;
	}

  //fpd
  if (mover_->get_additional_output())
      utility_exit_with_message("Movers returning multiple poses are unsupported by GenericMontoCarloMover.");
	if( saved_accept_file_name_ != "" ){
		TR<<"Saving initial pose entering the MC trajectory, for use in checkpoint recovery. Checkpointing filename: "<<saved_accept_file_name_<<std::endl;
		pose.dump_pdb( saved_accept_file_name_ );
	}

  PoseOP initial_pose = new Pose( pose );
	reset( pose ); //(re)initialize MC statistics
  MoverStatus ms( FAIL_RETRY );
	core::Size accept( 0 ), reject( 0 );
	using namespace protocols::rosetta_scripts;
	ParsedProtocolOP mover_pp( dynamic_cast< ParsedProtocol * >( mover_() ) );
	if( adaptive_movers() ){
		bool is_single_random( mover_pp->mode() == "single_random" );
		if( mover_pp && !is_single_random ){ // dig in one level (at most) to find the correct ParsedProtocol; if this becomes more generally useful then it would make sense to generatlize this to look for all parsedprotocols of type single_random that are being called by the MC mover. A simple recursion could do it, but I'm not sure how useful this would be
			foreach( ParsedProtocol::mover_filter_pair const mfp, *mover_pp ){
				ParsedProtocolOP tmp( dynamic_cast< ParsedProtocol * >( mfp.first.first() ) );
				if( tmp && tmp->mode() == "single_random" ){/// the parsedprotocol mover must be run in mode single_random for the apply_probabilities to be modified
					mover_pp = tmp;
					is_single_random = true;
				}
			}//foreach mfp
		}//fi mover_pp && is_single_random
		runtime_assert( is_single_random ); /// yes, we found a single-random parsedprotocol mover somewhere; notice that this part boils down to finding the first parsed protocol of single-random mode up to a depth level of 2
	}//fi adaptive_movers()
	utility::vector1< core::Size > mover_accepts;// count how many accepts each mover in the parsedprotocol made
	mover_accepts.clear();
	if( adaptive_movers() ){
		runtime_assert( mover_pp );
		runtime_assert( mover_pp->mode() == "single_random" );
		mover_accepts = utility::vector1< core::Size >( mover_pp->size(), 1 ); /// each mover gets a pseudocount of 1. This ensures that the running probability of each mover never goes to 0
	}
  for( Size i=load_trial_number_from_checkpoint(); i<=maxtrials_; i++ ){
    TR<<"Trial number: "<<i<<std::endl;
		if( i == 1 ){
			foreach( protocols::filters::FilterOP comp_statement_filt, filters_ ){ /// User defined filters in RosettaScripts are all CompoundFilter, so poke inside...
				protocols::filters::CompoundFilterOP comp_filt_op( dynamic_cast< protocols::filters::CompoundFilter * >( comp_statement_filt() ) );
				runtime_assert( comp_filt_op );
				for( protocols::filters::CompoundFilter::CompoundStatement::iterator cs_it = comp_filt_op->begin(); cs_it != comp_filt_op->end(); ++cs_it ){
					protocols::filters::FilterOP filt( cs_it->first );
					if( filt->get_type() == "Operator" ){
						TR<<"Resetting Operator filter's baseline"<<std::endl;
						protocols::simple_filters::OperatorOP operator_filter( dynamic_cast< protocols::simple_filters::Operator * >( filt() ) );
						operator_filter->reset_baseline( pose );
					}// fi Operator
				}// for cs_it
			} //foreach
		}//fi i == 1
		if( i > 1 && adaptive_movers() && i % adaptation_period() == 0 ){
/// The probability for each mover within a single-random parsedprotocol is determined by the number of accepts it had during the previous adaptation period:
/// each mover is assigned a pseducount of 1, and then any additional accept favors it over others. At the adaptation stage, the total number of accepts (including pseudocounts) is used to normalize the individual movers' number of accepts and the probability is the mover's accepts / by the total accepts
			core::Size sum_prev_accepts( 0 );
			foreach( core::Size const a, mover_accepts )
				sum_prev_accepts += a;
			utility::vector1< core::Real > new_probabilities;
			new_probabilities = utility::vector1< core::Real >( mover_accepts.size(), 0.0 );
			TR<<"Adapting running probabilities: old/new probabilities per mover\n";
			for( core::Size ma = 1; ma<= mover_accepts.size(); ma++ ){
				new_probabilities[ ma ] = ( core::Real )mover_accepts[ ma ] / ( core::Real )sum_prev_accepts;
				TR<<mover_pp->apply_probability()[ ma ]<<"/"<<new_probabilities[ ma ]<<'\n';
			}
			TR.flush();
			mover_pp->apply_probability( new_probabilities );
			mover_accepts = utility::vector1< core::Size >( mover_accepts.size(), 1 );
		}
		bool const stop( ( mover_stopping_condition_() != NULL && mover_stopping_condition_->obj ) || stopping_condition()->apply( pose ) );
		if( stop ){
			TR<<"MC stopping condition met at trial "<<i<<std::endl;
			break;
		}
    Pose store_pose( pose );
    // Mover apply
    mover_->apply( pose );
    ms = mover_->get_last_move_status();
    if( ms == FAIL_RETRY ){
      TR.Warning << "Mover failed. The mover, " << mover_->get_name() << ", is performed again. " << std::endl;
      pose = store_pose;
      continue;
    }else if( ms == FAIL_DO_NOT_RETRY || ms == FAIL_BAD_INPUT ){
      TR.Error << "Mover failed. Exit from GenericMonteCarloMover." << std::endl;
      break;
    }
    // MonteCarlo
    if( preapply_ && i==1 ){ // Auto-accept first application in order to deal with movers that e.g. change the length of the pose.
      reset( pose );
    }else{
      if( ! boltzmann( pose ) ){ // evaluate pose by scorefxn_ or filter_.report_sm()
        pose = *last_accepted_pose();
				reject++;
      }
			else{
				accept++;
				if( adaptive_movers() )
					mover_accepts[ mover_pp->last_attempted_mover_idx() ]++;
			}
    }
    if( !drift_ ){
      pose = (*initial_pose); // set back pose to initial one, of course this is not way of Monte Carlo
    }

    // Iterate over the list of triggers, executing each of them in turn
    fire_all_triggers(i, maxtrials_, pose, score_function());
		save_trial_number_to_checkpoint( i );
  } // i<=maxtrials_

	// Output final diagnositics, for potential tuning
	show_scores(TR);
	show_counters(TR);
	TR<<"Finished MC. Out of "<<maxtrials_<<" "<<accept<<" accepted "<<" and "<<reject<<" rejected."<<std::endl;
  // Recover pose that have the lowest score, or the last accepted pose, as appropriate
	if(recover_low_) {
		recover_low( pose );
	} else {
		if (last_accepted_pose()) {
			pose = *last_accepted_pose();
		}
	}

}// apply

std::string
GenericMonteCarloMover::get_name() const {
  return GenericMonteCarloMoverCreator::mover_name();
}

/// @brief parse xml file
void
GenericMonteCarloMover::parse_my_tag( TagPtr const tag, DataMap & data, Filters_map const &filters, Movers_map const &movers, Pose const & )
{
	//using core::pack::task::operation::TaskOperation;
	//using core::pack::task::TaskFactoryOP;
	//using core::pack::task::TaskFactory;

	maxtrials_ = tag->getOption< core::Size >( "trials", 10 );
	temperature_ = tag->getOption< Real >( "temperature", 0.0 );
	task_scaling_ = tag->getOption< core::Size >( "task_scaling", 5 );
	adaptive_movers( tag->getOption< bool >( "adaptive_movers", false ) );
	if( adaptive_movers() )
		adaptation_period( tag->getOption< core::Size >( "adaptation_period", std::max( (int) maxtrials_ / 10, 10 ) ) );

	String const  user_defined_mover_name_( tag->getOption< String >( "mover_name" ,""));
	if( data.has( "stopping_condition", user_defined_mover_name_ ) ){
		TR<<user_defined_mover_name_<<" defines its own stopping condition, and GenericMC will respect this stopping condition"<<std::endl;
		mover_stopping_condition_ = data.get< DataMapObj< bool > * >( "stopping_condition", user_defined_mover_name_ );
	}

	String const filter_name( tag->getOption< String >( "filter_name", "true_filter" ) );
	Movers_map::const_iterator  find_mover ( movers.find( user_defined_mover_name_ ));
	Filters_map::const_iterator find_filter( filters.find( filter_name ));
	if( find_mover == movers.end() && user_defined_mover_name_ != "" ) {
		TR.Error << "ERROR !! mover not found in map: \n" << tag << std::endl;
		runtime_assert( find_mover != movers.end() );
	}
	if( find_filter == filters.end() ) {
		TR.Error << "ERROR !! filter not found in map: \n" << tag << std::endl;
		runtime_assert( find_filter != filters.end() );
	}
	if( user_defined_mover_name_ != "" )
		mover_ = find_mover->second;

	if( adaptive_movers() ) /// adaptive movers only works if the mover being called is of type parsedprotocol
		runtime_assert( dynamic_cast< protocols::rosetta_scripts::ParsedProtocol * >( mover_() ));

	bool const adaptive( tag->getOption< bool >( "adaptive", true ) );
	add_filter( find_filter->second->clone(), adaptive, temperature_, sample_type_ );
	String const sfxn ( tag->getOption< String >( "scorefxn_name", "" ) );
	if( sfxn != "" ){
		//scorefxn_ = new ScoreFunction( *data.get< ScoreFunction * >( "scorefxns", sfxn ));
		scorefxn_ = data.get< ScoreFunction * >( "scorefxns", sfxn )->clone();   //fpd use clone
		TR << "Score evaluation during MC is done by" << sfxn << ", ";
		TR << filter_name << " is ignored." << std::endl;
		filters_.clear();
	}else{
		scorefxn_ = NULL;
	}

	parse_task_operations( tag, data, filters, movers );

	if( filter_name == "true_filter" && !scorefxn_ ){
		TR.Error << "You need to set filter_name or scorefxn_name for MC criteria." << std::endl;
		runtime_assert( false );
	}

	if( filters_.size() == 0 ){
		TR << "Apply mover of " << user_defined_mover_name_ << ", and evaluate score by " << sfxn
		<< " at Temperature=" << temperature_ << ", ntrails= " << maxtrials_ << std::endl;
	}

	stopping_condition( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "stopping_condition", "false_filter" ), filters ) );
	if( tag->hasOption( "stopping_condition" ) )
		TR<<"Generic MC using stopping condition "<< stopping_condition()->get_user_defined_name()<<std::endl;
	drift_ = tag->getOption< bool >( "drift", 1 );
	preapply_ = tag->getOption< bool >( "preapply", 1 ); // Default true for historical reasons
	recover_low_ = tag->getOption< bool >( "recover_low", 1 );
	boltz_rank_ = tag->getOption< bool >( "bolz_rank", 0 );
	sample_type_ = tag->getOption< String >( "sample_type", "low" );

	utility::vector1< TagPtr > const branch_tags( tag->getTags() );
	foreach( TagPtr const btag, branch_tags ){
		if( btag->getName() == "Filters" ){
			utility::vector1< TagPtr > const filters_tags( btag->getTags() );
			foreach( TagPtr const ftag, filters_tags ){
				String const filter_name( ftag->getOption< String >( "filter_name" ) );
				Filters_map::const_iterator find_filt( filters.find( filter_name ));
				if( find_filt == filters.end() ) {
					TR.Error << "Error !! filter not found in map: \n" << tag << std::endl;
					runtime_assert( find_filt != filters.end() );
				}
				Real const temp( ftag->getOption< Real >( "temperature", 1 ) );
				bool const adap( ftag->getOption< bool >( "adaptive", true ));
				String const samp_type( ftag->getOption< String >( "sample_type", "low" ));
				bool const rank( ftag->getOption< bool >( "rank", false ) );
				if( rank ) {
					if( rank_by_filter_ != 1 ) {
						TR.Warning << "WARNING Multiple filters set to rank! Using most recent (currently "<< filter_name << ")." << std::endl;
					}
					if( boltz_rank_ ) {
						TR.Warning << "WARNING Setting of rank on sub-filter "<< filter_name << " will be ignored as parent is set to Boltzmann rank all." << std::endl;
					}
				}
				add_filter( find_filt->second, adap, temp, samp_type, rank );
			} //foreach ftag
		}// fi Filters
		else
			utility_exit_with_message( "tag name " + btag->getName() + " unrecognized." );
	}//foreach btag

  if( !drift_ ){
    TR << "Pose is set back to initial pose every after applying mover and score evaluation." << std::endl;
  }

  if( sample_type_ == "high" ){
    TR << "Pose that have higher score is sampled." << sfxn << std::endl;
  }

	saved_accept_file_name_ = tag->getOption< std::string >( "saved_accept_file_name", "" );
	saved_trial_number_file_ = tag->getOption< std::string >( "saved_trial_number_file", "" );
	if( tag->hasOption( "mover_tag" ) )
		mover_tag_ = protocols::moves::get_set_from_datamap< DataMapObj< std::string > >( "tags", tag->getOption< std::string >( "mover_tag" ), data );
  initialize();
}

///@brief parse "task_operations" XML option
void GenericMonteCarloMover::parse_task_operations(
	TagPtr const tag,
	DataMap const & datamap,
	Filters_map const &,
	Movers_map const &
)
{
	if ( ( tag->hasOption("task_operations") )){
			TR << "Found a task operation" << std::endl;
			TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
			if ( new_task_factory == 0) return;
			task_factory( new_task_factory );
		}
}


void GenericMonteCarloMover::fire_all_triggers(
	Size cycle,
	Size num_cycles,
	const Pose& pose,
	ScoreFunctionOP scoring)
{
	boost::unordered_map<Size, GenericMonteCarloMoverTrigger>::iterator i;
	for (i = triggers_.begin(); i != triggers_.end(); ++i) {
		GenericMonteCarloMoverTrigger& t = i->second;
		bool rescore = t(cycle, num_cycles, pose, scoring);

		if (rescore) {
			last_accepted_score_ = scoring->score(*last_accepted_pose_);
			lowest_score_ = scoring->score(*lowest_score_pose_);
		}
	}
}

Size GenericMonteCarloMover::add_trigger(const GenericMonteCarloMoverTrigger& trigger) {
  Size tid = next_trigger_id_++;
  triggers_[tid] = trigger;
  return tid;
}

void GenericMonteCarloMover::remove_trigger(Size trigger_id) {
  boost::unordered_map<Size, GenericMonteCarloMoverTrigger>::iterator i = triggers_.find(trigger_id);
  if (i == triggers_.end()) {
    TR.Warning << "Attempt to remove invalid trigger_id => " << trigger_id << std::endl;
    return;
  }
  triggers_.erase(i);
}

Size GenericMonteCarloMover::num_triggers() const {
  return triggers_.size();
}

void
GenericMonteCarloMover::stopping_condition( protocols::filters::FilterOP f ){
	stopping_condition_ = f;
}

protocols::filters::FilterOP
GenericMonteCarloMover::stopping_condition() const{
	return stopping_condition_;
}

std::string
GenericMonteCarloMover::saved_accept_file_name() const{
	return saved_accept_file_name_;
}

void
GenericMonteCarloMover::saved_accept_file_name( std::string const s ){
	saved_accept_file_name_ = s;
}

std::string
GenericMonteCarloMover::saved_trial_number_file() const{
	return saved_trial_number_file_;
}

void
GenericMonteCarloMover::saved_trial_number_file( std::string const s ){
	saved_trial_number_file_ = s;
}
} // ns moves
} // ns protocols
