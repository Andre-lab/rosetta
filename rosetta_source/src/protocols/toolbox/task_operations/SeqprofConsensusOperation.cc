// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/SeqprofConsensusOperation.cc
/// @brief set designable residues to those observed in sequence profile
/// @author Florian Richter, floric@u.washington.edu, april 2011


// Unit Headers
#include <protocols/toolbox/task_operations/SeqprofConsensusOperation.hh>
#include <protocols/toolbox/task_operations/SeqprofConsensusOperationCreator.hh>

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/sequence/SequenceProfile.hh>

#include <core/io/ddg/PositionDdGInfo.hh>

// Utility Headers
#include <utility/string_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>

#include <utility/vector0.hh>


static basic::Tracer tr("protocols.toolbox.tas_operations.SeqprofConsensusOperation");

namespace protocols{
namespace toolbox{
namespace task_operations{

core::pack::task::operation::TaskOperationOP
SeqprofConsensusOperationCreator::create_task_operation() const
{
	return new SeqprofConsensusOperation;
}


/// @brief default constructor
SeqprofConsensusOperation::SeqprofConsensusOperation():
	TaskOperation(),
	seqprof_filename_( basic::options::option[ basic::options::OptionKeys::in::file::pssm ][1] ),
	seqprof_(NULL),
	min_aa_probability_(0.0),
	prob_larger_current_(true),
	ignore_pose_profile_length_mismatch_(false)
{
	if( utility::file::file_exists( seqprof_filename_ ) ){
		core::sequence::SequenceProfileOP seqprof = new core::sequence::SequenceProfile( seqprof_filename_ );
		seqprof->convert_profile_to_probs(); // was previously implicit in from-filename constructor
		seqprof_ = seqprof;
	}
}


/// @brief destructor
SeqprofConsensusOperation::~SeqprofConsensusOperation() {}

/// @brief clone
core::pack::task::operation::TaskOperationOP
SeqprofConsensusOperation::clone() const {
	return new SeqprofConsensusOperation( *this );
}

/// @brief all AA that have a higher probability in the seqprofile
/// than the native residue are allowed. probability also
/// needs to be higher than min_aa_probability_
/// @details NOTE ON SYMMETRIC POSE BEHAVIOR:
/// pssm files are usually for one chain only, therefore
/// this task operation will only set the residue behavior for
/// the first chain/asymetric unit.
/// it could be possible to handle the symmetry setup here, i.e.
/// set up the residue level task for every symmetric copy, but
/// it's prolly better to let the symmetry machinery deal with that
/// mode of packer task symmetrization should be intersection
void
SeqprofConsensusOperation::apply( Pose const & pose, PackerTask & task ) const
{
	if( !seqprof_) utility_exit_with_message("No sequence profile set. option -in:file:pssm not specified? no filename in tag specified?");

	core::Size asymmetric_unit_res( pose.total_residue() );
	if ( core::pose::symmetry::is_symmetric(pose) ) {
    core::conformation::symmetry::SymmetricConformation const & SymmConf (
      dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
    asymmetric_unit_res = SymmConf.Symmetry_Info()->num_independent_residues();
		task.request_symmetrize_by_intersection();
  }
	core::Size last_res (asymmetric_unit_res <= seqprof_->profile().size() ? pose.total_residue() : seqprof_->profile().size() );
	for( core::Size i = 1; i <= last_res; ++i){

		if( !pose.residue_type( i ).is_protein() ) continue;
		//std::cout << "SCO at pos " << i << " allows the following residues: ";
		utility::vector1< Real > const & pos_profile( (seqprof_->profile())[ i ] );
		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		core::Real current_prob( pos_profile[ pose.residue_type(i).aa() ] );

		for( core::Size aa = core::chemical::aa_ala; aa <= core::chemical::num_canonical_aas; ++aa){
			core::Real prob( pos_profile[ aa ] );
			if( prob >= min_aa_probability_ ){
				if( prob_larger_current_) {
					if( prob >= current_prob ) keep_aas[ aa ] = true;
				}
				else	keep_aas[ aa ] = true;
				//std::cout << " " << static_cast<core::chemical::AA>(aa) << " prob=" << prob << ", ";
			}
		}
		keep_aas[  pose.residue_type(i).aa() ] = true; //current always allowed
		//std::cout << " native " << pose.residue_type(i).aa() << " prob=" << native_prob << "." << std::endl;

		task.nonconst_residue_task(i).restrict_absent_canonical_aas( keep_aas );

	} //loop over all residues for which profile information exists

	bool prot_res_without_profile_information_exist(false);
	for( core::Size i = last_res + 1; i <= asymmetric_unit_res; ++i){
		task.nonconst_residue_task(i).restrict_to_repacking();
		if( pose.residue_type( i ).is_protein() ) prot_res_without_profile_information_exist = true;
	}

	if( prot_res_without_profile_information_exist ){
		if( ignore_pose_profile_length_mismatch_ ) tr << "WARNING WARNING: the passed in pose is longer than the sequence profile specified. Double check whether the used sequence profile is correct. Setting every excess pose residue to repacking.";

		else utility_exit_with_message("The passed in pose is longer than the sequence profile specified. Double check whether the used sequence profile is correct.");
	}
} // apply

void
SeqprofConsensusOperation::parse_tag( TagPtr tag )
{
	if( tag->hasOption("filename") ){
		seqprof_filename_ = tag->getOption< String >( "filename" );
		tr<<"Loading seqprof from a file named: "<<seqprof_filename_<<std::endl;
		core::sequence::SequenceProfileOP seqprof = new core::sequence::SequenceProfile( seqprof_filename_ );
		seqprof->convert_profile_to_probs(); // was previously implicit in from-filename constructor
		seqprof_ = seqprof;
	}
	else{
		tr<<"Seqprof not loaded. Expecting another mover/filter to provide a sequence profile..."<<std::endl;
	}
	if( tag->hasOption("min_aa_probability") ) min_aa_probability_ = tag->getOption< Real >("min_aa_probability" );
	if( tag->hasOption("probability_larger_than_current") ) prob_larger_current_ = tag->getOption< bool >("probability_larger_than_current");

	if( tag->hasOption("ignore_pose_profile_length_mismatch") ) ignore_pose_profile_length_mismatch_ = tag->getOption< bool >("ignore_pose_profile_length_mismatch");
}

core::sequence::SequenceProfileCOP
SeqprofConsensusOperation::seqprof() const
{
	return seqprof_;
}

void
SeqprofConsensusOperation::set_seqprof( core::sequence::SequenceProfileOP seqprof, bool reweight )
{
	if ( reweight ) {
		core::sequence::SequenceProfileOP reweightedprof = new core::sequence::SequenceProfile( *seqprof );
		reweightedprof->convert_profile_to_probs(); // was previously implicit in from-filename constructor
		seqprof_ = reweightedprof;
	} else {
		seqprof_ = seqprof;
	}
}

core::pack::task::operation::TaskOperationOP
RestrictConservedLowDdgOperationCreator::create_task_operation() const
{
	return new RestrictConservedLowDdgOperation;
}

RestrictConservedLowDdgOperation::RestrictConservedLowDdgOperation()
: Parent(),
	ddG_predictions_filename_(basic::options::option[ basic::options::OptionKeys::in::file::ddg_predictions_file ].value()),
	conservation_cutoff_(0.6),
	ddG_cutoff_(1.5),
	verbose_(false)
{
	position_ddGs_.clear();
	if( utility::file::file_exists( ddG_predictions_filename_ ) ){
		position_ddGs_ = core::io::PositionDdGInfo::read_ddg_predictions_file( ddG_predictions_filename_ );
	}
}

RestrictConservedLowDdgOperation::~RestrictConservedLowDdgOperation()
{}

core::pack::task::operation::TaskOperationOP
RestrictConservedLowDdgOperation::clone() const
{
	return new RestrictConservedLowDdgOperation( *this );
}

void
RestrictConservedLowDdgOperation::parse_tag( TagPtr tag )
{
	Parent::parse_tag( tag );
	if( tag->hasOption("ddG_filename") ){
		ddG_predictions_filename_ = tag->getOption< std::string >("ddG_filename" );
		position_ddGs_ = core::io::PositionDdGInfo::read_ddg_predictions_file( ddG_predictions_filename_ );
	}

	if( tag->hasOption("conservation_cutoff") ) conservation_cutoff_ = tag->getOption< Real >("conservation_cutoff" );
	if( tag->hasOption("ddG_cutoff") ) ddG_cutoff_ = tag->getOption< Real >("ddG_cutoff" );
	if( tag->hasOption("verbose") ) verbose_ = tag->getOption< bool >("verbose" );
}

void
RestrictConservedLowDdgOperation::apply(
	Pose const & pose,
	PackerTask & task
) const
{
	if( !this->seqprof()) utility_exit_with_message("No sequence profile set. option -in:file:pssm not specified? no filename in tag specified?");

	if( position_ddGs_.size() == 0 ) utility_exit_with_message("No ddG infos were read in. option -in:file:ddg_predictions_file not specified? no filename in tag specified?");

	for( core::Size i = 1; i <= pose.total_residue(); ++i){

		if( !pose.residue_type( i ).is_protein() ) continue;
		core::chemical::AA seqprof_wt_aa( this->seqprof_wt_aa( i ) );

		if( position_untouchable( i, seqprof_wt_aa ) ){
			if( seqprof_wt_aa == pose.residue_type(i).aa() ) task.nonconst_residue_task( i ).restrict_to_repacking();
			else{
				utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
				keep_aas[ seqprof_wt_aa ] = true;
				keep_aas[ pose.residue_type(i).aa() ] = true;
				task.nonconst_residue_task(i).restrict_absent_canonical_aas( keep_aas );
			}
		} // if untouchable
	} //loop over all residues
}

bool
RestrictConservedLowDdgOperation::position_untouchable(
	core::Size seqpos,
	core::chemical::AA seqprof_wt
) const
{

	//note: first we deal with the alanine special case
	//obviousluy there is no ddG associated with mutating Ala to Ala,
	//so for alanine residues we return true through the conservation
	//criterion alone, the rationale being that conserved alanines
	//are probably important structurally
	if( seqprof_wt == core::chemical::aa_ala ){
		if( (seqprof()->profile())[ seqpos ][ seqprof_wt ] > conservation_cutoff_) return true;
		else return false;
	}

	std::map< core::Size, core::io::PositionDdGInfo::PositionDdGInfoOP >::const_iterator posddg_it( position_ddGs_.find( seqpos ) );

	//note: if the position wasn't found, this could mean that the predictions file
	//was incomplete or that the ddG protocol couldn't calculate a proper ddG,
	//which is the case for disulfides and some modified residue types such as
	//phospho-ser etc. so let's spit out a warning and just apply the conservation_cutoff_
	if( posddg_it == position_ddGs_.end() ){

		tr << "Warning: no ddG information read for sequence position " << seqpos << ". This could either mean that the ddG predictions input file is incomplete or that the original PDB had a disulfide cys or other modified residue at this position. Decision whether residue is untouchable will be made based on sequence conservation alone." << std::endl;
		if( (seqprof()->profile())[ seqpos ][ seqprof_wt ] > conservation_cutoff_) return true;
		else return false;
	}

	core::io::PositionDdGInfo::PositionDdGInfo const & pos_ddg( *(posddg_it->second) );
	if( seqprof_wt != pos_ddg.wt_aa() ) utility_exit_with_message ("The wildtype aa for position "+utility::to_string( seqpos ) + " is different in the ddG file and the pssm file. Something's unclean somewhere." );
	std::map< core::chemical::AA, core::Real >::const_iterator ala_it( pos_ddg.mutation_ddGs().find( core::chemical::aa_ala ) );
	if( ala_it ==  pos_ddg.mutation_ddGs().end() ) utility_exit_with_message("The ddG of mutating to Ala was not found for position "+utility::to_string( seqpos )+" in file "+ddG_predictions_filename_ + ".");

	if( (ala_it->second > ddG_cutoff_) && ( (seqprof()->profile())[ seqpos ][ seqprof_wt ] > conservation_cutoff_) ){
		if( verbose_ ) tr << "Pos " << seqprof_wt << seqpos << " has ddG_cutoff of " << ala_it->second << " and profile frequency of " << (seqprof()->profile())[ seqpos ][ seqprof_wt ] << ", considered untouchable." << std::endl;
		return true;
	}
	return false;
}

core::chemical::AA
RestrictConservedLowDdgOperation::seqprof_wt_aa( core::Size seqpos ) const
{
	return core::chemical::aa_from_oneletter_code( (*(this->seqprof()))[seqpos] );
}

core::Real
RestrictConservedLowDdgOperation::position_ala_ddG( core::Size seqpos ) const
{
	if( seqprof_wt_aa( seqpos ) == core::chemical::aa_ala ) return 0.0;

	std::map< core::Size, core::io::PositionDdGInfo::PositionDdGInfoOP >::const_iterator posddg_it( position_ddGs_.find( seqpos ) );
	if( posddg_it == position_ddGs_.end() ) utility_exit_with_message("no ddg information read for sequence position "+ utility::to_string( seqpos ) );
	core::io::PositionDdGInfo::PositionDdGInfo const & pos_ddg( *(posddg_it->second) );
	std::map< core::chemical::AA, core::Real >::const_iterator ala_it( pos_ddg.mutation_ddGs().find( core::chemical::aa_ala ) );
	if( ala_it ==  pos_ddg.mutation_ddGs().end() ) utility_exit_with_message("The ddG of mutating to Ala was not found for position "+utility::to_string( seqpos )+" in file "+ddG_predictions_filename_ + ".");
	return ala_it->second;
}

} // TaskOperations
} // toolbox
} // protocols

