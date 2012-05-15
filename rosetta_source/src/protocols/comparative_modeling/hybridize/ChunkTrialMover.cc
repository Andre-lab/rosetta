// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Wrapper for InsertChunkMover. It can take a random template and steal coordinates of all chunks or a random one
/// @detailed
/// @author Yifan Song

#include <protocols/comparative_modeling/hybridize/ChunkTrialMover.hh>

#include <core/pose/PDBInfo.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static numeric::random::RandomGenerator RG(57029435);
static basic::Tracer TR( "protocols.comparative_modeling.hybridize.ChunkTrialMover" );

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace core::kinematics;
using namespace ObjexxFCL;
using namespace protocols::moves;
using namespace protocols::loops;
using namespace numeric::model_quality;
using namespace id;
using namespace basic::options;
using namespace basic::options::OptionKeys;
	
ChunkTrialMover::ChunkTrialMover(
					utility::vector1 < core::pose::PoseCOP > const & template_poses,
					utility::vector1 < protocols::loops::Loops > const & template_chunks,
					Loops ss_chunks_pose,
					bool random_template,
					AlignOption align_option,
					Size max_registry_shift ) :
    template_poses_(template_poses),
    template_chunks_(template_chunks),
    random_template_(random_template),
    align_option_(align_option),
    align_chunk_(),
    max_registry_shift_input_(max_registry_shift)
{
	bool alignment_from_template = option[cm::hybridize::alignment_from_template_seqpos]();
	
	Size count = 0;
	for (core::Size i_template=1; i_template<=template_poses_.size(); ++i_template) {
		if (template_chunks_[i_template].size() != 0) ++count;
	}
	if (count == 0) {
		utility_exit_with_message("Template structures need at least one secondary structure for this protocol");
	}
	
	sequence_alignments_.clear();
	for (core::Size i_template=1; i_template<=template_poses_.size(); ++i_template) {
		std::map <core::Size, core::Size> sequence_alignment;
        //TR << "template: " << i_template << std::endl;
        sequence_alignment.clear();
		if (alignment_from_template) {
			get_alignment_from_template(template_poses_[i_template], sequence_alignment);
		}
		else {
			std::map <core::Size, core::Size> chunk_mapping;
			if (option[cm::hybridize::alignment_from_chunk_mapping].user()) {
				for (Size i=1;i<=option[cm::hybridize::alignment_from_chunk_mapping]().size();++i) {
					chunk_mapping[i] = option[cm::hybridize::alignment_from_chunk_mapping]()[i];
				}
			}
			
			get_alignment_from_chunk_mapping(chunk_mapping, template_chunks_[i_template], ss_chunks_pose, sequence_alignment);
		}
		sequence_alignments_.push_back(sequence_alignment);
	}
}

void
ChunkTrialMover::get_alignment_from_template(
                                             core::pose::PoseCOP template_pose,
                                             std::map <core::Size, core::Size> & seqpos_alignment
                                             ) {
	// specific to this case, alignment comes from residue number
	for (core::Size ires=1; ires<=template_pose->total_residue(); ++ires) {
		//TR << "Sequence aln: " << template_pose->pdb_info()->number(ires) << " " << ires << std::endl;
		seqpos_alignment[template_pose->pdb_info()->number(ires)] = ires;
	}
}

void
ChunkTrialMover::get_alignment_from_chunk_mapping(std::map <core::Size, core::Size> const & chunk_mapping,
								 Loops const template_ss_chunks,
								 Loops const target_ss_chunks,
								 std::map <core::Size, core::Size> & sequence_alignment)
{
	max_registry_shift_.resize(target_ss_chunks.size());
	for (Size i_chunk_pose = 1; i_chunk_pose <= target_ss_chunks.size(); ++i_chunk_pose) {
		if (chunk_mapping.find(i_chunk_pose) == chunk_mapping.end()) continue;
		core::Size j_chunk_template = chunk_mapping.find(i_chunk_pose)->second;
		
		Size respos_mid_pose = (target_ss_chunks[i_chunk_pose].start() + target_ss_chunks[i_chunk_pose].stop()) / 2;
		Size respos_mid_template = (template_ss_chunks[j_chunk_template].start() + template_ss_chunks[j_chunk_template].stop()) / 2;
		int offset = respos_mid_template - respos_mid_pose;
		
		using namespace ObjexxFCL::fmt;
		if (target_ss_chunks[i_chunk_pose].length() <= template_ss_chunks[j_chunk_template].length()) {
			max_registry_shift_[i_chunk_pose] = max_registry_shift_input_ + template_ss_chunks[j_chunk_template].length() - target_ss_chunks[i_chunk_pose].length();
			for (Size ires=target_ss_chunks[i_chunk_pose].start(); ires<=target_ss_chunks[i_chunk_pose].stop(); ++ires) {
				sequence_alignment[ires] = ires+offset;
				//std::cout << I(4, ires) << I(4, ires+offset) << std::endl;
			}
		}
		else {
			max_registry_shift_[i_chunk_pose] = max_registry_shift_input_ + target_ss_chunks[i_chunk_pose].length() - template_ss_chunks[j_chunk_template].length();
			for (Size ires_templ=template_ss_chunks[j_chunk_template].start(); ires_templ<=template_ss_chunks[j_chunk_template].stop(); ++ires_templ) {
				sequence_alignment[ires_templ-offset] = ires_templ;
				//std::cout << I(4, ires_templ) << I(4, ires_templ-offset) << std::endl;
			}
		}
	}
}
	
void ChunkTrialMover::pick_random_template()
{
	assert(template_poses_.size() != 0);

    set_template(0);
	while (!template_number()) {
		set_template( RG.random_range(1, template_poses_.size()) );
		if (template_chunks_[template_number()].size() == 0) set_template(0);
	}
}

void ChunkTrialMover::set_template(core::Size const template_number)
{
    template_number_ = template_number;
}

core::Size ChunkTrialMover::template_number()
{
    return template_number_;
}
    
void ChunkTrialMover::pick_random_chunk(core::pose::Pose & pose) {
	int ntrials=500;
	jump_number_ = RG.random_range(1, pose.num_jump());
	core::Size jump_residue_pose = pose.fold_tree().downstream_jump_residue(jump_number_);
	while ( pose.residue(jump_residue_pose).aa() == core::chemical::aa_vrt && --ntrials>0) {
		jump_number_ = RG.random_range(1, pose.num_jump());
		jump_residue_pose = pose.fold_tree().downstream_jump_residue(jump_number_);	
	}
	if (ntrials == 0) {
		utility_exit_with_message( "Fatal error in ChunkTrialMover::pick_random_chunk()");
	}
}

Size ChunkTrialMover::trial_counter(Size ires) {
	return align_chunk_.trial_counter(ires);	
}
	
void
ChunkTrialMover::apply(core::pose::Pose & pose) {
	max_registry_shift_.resize(pose.num_jump(), max_registry_shift_input_);

	// pick a random template
	if (random_template_) {
		pick_random_template();
  }
  //TR << "templ number: " << template_number() << std::endl;
	align_chunk_.set_template(template_poses_[template_number()], template_number(), sequence_alignments_[template_number()]);

	// random chunk or loop all chunks
	if (align_option_ == random_chunk) {
		// pick a random jump
		pick_random_chunk(pose);
		align_chunk_.set_aligned_chunk(pose, jump_number_, false);
		align_chunk_.set_reset_torsion_unaligned(false);

		// apply alignment
		int registry_shift = RG.random_range(-max_registry_shift_[jump_number_], max_registry_shift_[jump_number_]);
		align_chunk_.set_registry_shift(registry_shift);
		align_chunk_.apply(pose);
	} else {
		// loop over all jumps (we're initializing)
		align_chunk_.reset_torsion(true);
		for (core::Size jump_number=1; jump_number<=pose.num_jump(); ++jump_number) {
			align_chunk_.set_aligned_chunk(pose, jump_number, true);

			// apply alignment
			int registry_shift = RG.random_range(-max_registry_shift_[jump_number], max_registry_shift_[jump_number]);
			align_chunk_.set_registry_shift(registry_shift);
			align_chunk_.apply(pose);
			if (!align_chunk_.success()) {
				TR.Debug << "Warning! This chunk might not be aligned, jump number " << jump_number << std::endl;
			}
		}
	}
}

std::string ChunkTrialMover::get_name() const
{
	return "ChunkTrialMover";
}

} // hybridize 
} // comparative_modeling 
} // protocols

