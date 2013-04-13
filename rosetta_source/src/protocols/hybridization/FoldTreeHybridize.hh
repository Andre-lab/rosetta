// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Align a random jump to template
/// @detailed
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_FoldTreeHybridize_hh
#define INCLUDED_protocols_hybridization_FoldTreeHybridize_hh

#include <protocols/hybridization/InsertChunkMover.hh>
#include <protocols/hybridization/FoldTreeHybridize.fwd.hh>
#include <protocols/hybridization/HybridizeFoldtreeDynamic.hh>
#include <protocols/hybridization/WeightedFragmentTrialMover.hh>
#include <protocols/hybridization/WeightedFragmentSmoothTrialMover.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// pairings
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/JumpSample.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <set>

namespace protocols {
namespace hybridization {

using namespace core;
using namespace protocols::moves;
using namespace protocols::loops;

class FoldTreeHybridize: public protocols::moves::Mover {

public:
	FoldTreeHybridize();

	FoldTreeHybridize(
		core::Size const initial_template_index,
		utility::vector1 < core::pose::PoseOP > const & template_poses,
		utility::vector1 < core::Real > const & template_weights,
		utility::vector1 < protocols::loops::Loops > const & template_chunks,
		utility::vector1 < protocols::loops::Loops > const & template_contigs,
		core::fragment::FragSetOP fragments_small_in,
		core::fragment::FragSetOP fragments_big_in );

	// initialize options to defaults
	void init();

	void revert_loops_to_original(core::pose::Pose & pose, Loops loops);

	void set_loops_to_virt_ala(core::pose::Pose & pose, Loops loops);

	Real gap_distance(Size Seq_gap);

	void add_gap_constraints_to_pose(core::pose::Pose & pose, Loops const & chunks, int gap_edge_shift=0, Real stdev=0.1);

	void restore_original_foldtree(core::pose::Pose & pose);

	void setup_foldtree(core::pose::Pose & pose);

	numeric::xyzVector<Real> center_of_mass(core::pose::Pose const & pose);

	void translate_virt_to_CoM(core::pose::Pose & pose);

	utility::vector1< core::Real > get_residue_weights_for_1mers(core::pose::Pose & pose);
	utility::vector1< core::Real > get_residue_weights_for_small_frags(core::pose::Pose & pose);
	utility::vector1< core::Real > get_residue_weights_for_big_frags(core::pose::Pose & pose);

	void set_constraint_file(std::string cst_file_in) { cst_file_=cst_file_in; }
	void set_increase_cycles(core::Real increase_cycles_in) { increase_cycles_=increase_cycles_in; }
	void set_stage1_1_cycles(core::Size stage1_1_cycles_in) { stage1_1_cycles_=stage1_1_cycles_in; }
	void set_stage1_2_cycles(core::Size stage1_2_cycles_in) { stage1_2_cycles_=stage1_2_cycles_in; }
	void set_stage1_3_cycles(core::Size stage1_3_cycles_in) { stage1_3_cycles_=stage1_3_cycles_in; }
	void set_stage1_4_cycles(core::Size stage1_4_cycles_in) { stage1_4_cycles_=stage1_4_cycles_in; }
	void set_add_non_init_chunks(bool add_non_init_chunks_in) { add_non_init_chunks_=add_non_init_chunks_in; }
	void set_domain_assembly(bool domain_assembly_in) { domain_assembly_=domain_assembly_in; }
	void set_add_hetatm(bool add_hetatm_in, core::Real hetatm_self_cst_weight_in, core::Real hetatm_prot_cst_weight_in) { add_hetatm_=add_hetatm_in; hetatm_self_cst_weight_=hetatm_self_cst_weight_in; hetatm_prot_cst_weight_=hetatm_prot_cst_weight_in;}
	void set_frag_1mer_insertion_weight(core::Real frag_1mer_insertion_weight_in) { frag_1mer_insertion_weight_=frag_1mer_insertion_weight_in; }
	void set_small_frag_insertion_weight(core::Real small_frag_insertion_weight_in) { small_frag_insertion_weight_=small_frag_insertion_weight_in; }
	void set_big_frag_insertion_weight(core::Real big_frag_insertion_weight_in) { big_frag_insertion_weight_=big_frag_insertion_weight_in; }
	void set_frag_weight_aligned(core::Real frag_weight_aligned_in) { frag_weight_aligned_=frag_weight_aligned_in; }
	void set_auto_frag_insertion_weight(bool auto_frag_insertion_weight_in) { auto_frag_insertion_weight_ = auto_frag_insertion_weight_in; }
	void set_max_registry_shift(core::Size max_registry_shift_in) { max_registry_shift_=max_registry_shift_in; }
	void set_top_n_big_frag(core::Size top_n_big_frag_in) { top_n_big_frag_=top_n_big_frag_in; }
	void set_top_n_small_frag(core::Size top_n_small_frag_in) { top_n_small_frag_=top_n_small_frag_in; }

	// strand pairings
	void set_pairings_file( std::string pairings_file_in ) { pairings_file_ = pairings_file_in; }
	void set_sheets( utility::vector1< core::Size > sheets_in ) { sheets_ = sheets_in; }
	void set_random_sheets( utility::vector1< core::Size > random_sheets_in ) { random_sheets_ = random_sheets_in; }
	void set_filter_templates( bool filter_templates_in ) { filter_templates_ = filter_templates_in; }

	inline void set_scorefunction(core::scoring::ScoreFunctionOP const scorefxn) { scorefxn_ = scorefxn; }

	void set_movable_region( utility::vector1< bool > allowed_to_move_in ) { allowed_to_move_ = allowed_to_move_in; }

	void setup_scorefunctions(
		core::scoring::ScoreFunctionOP score0,
		core::scoring::ScoreFunctionOP score1,
		core::scoring::ScoreFunctionOP score2,
		core::scoring::ScoreFunctionOP score5,
		core::scoring::ScoreFunctionOP score3);

	void apply(core::pose::Pose & pose);

	std::string	get_name() const;

	// strand pairings
	utility::vector1< std::pair< core::Size, core::Size > > get_strand_pairs() { return strand_pairs_; };
	std::set< core::Size > get_pairings_residues();

private:

	void normalize_template_wts();

	// automatically set the fragment insertion weight
	void auto_frag_insertion_weight(
		WeightedFragmentTrialMoverOP & frag_1mer_trial_mover,
		WeightedFragmentTrialMoverOP & small_frag_trial_mover,
		WeightedFragmentTrialMoverOP & big_frag_trial_mover
	);

	utility::vector1< core::Size > get_jump_anchors();

	// strand pairings
	void add_strand_pairings();
	void add_strand_pairing( core::scoring::dssp::Pairing const & pairing );
	void filter_templates( std::set< core::Size > const & templates_to_remove );
	void superimpose_strand_pairings_to_templates(core::pose::Pose & pose);
	core::Size map_pdb_info_number( const core::pose::Pose & pose, core::Size pdb_res );
	protocols::simple_moves::ClassicFragmentMoverOP get_pairings_jump_mover();


private:
	core::Real increase_cycles_;
	core::Size stage1_1_cycles_;
	core::Size stage1_2_cycles_;
	core::Size stage1_3_cycles_;
	core::Size stage1_4_cycles_;

	// 1mer fragment insertion weight where large and small fragments are not allowed (across anchor points) vs. chunk insertion + big and small fragments
	core::Real frag_1mer_insertion_weight_;
	// small fragment insertion weight where large fragments are not allowed (across anchor points) vs. chunk insertion + big fragments
	core::Real small_frag_insertion_weight_;
	// fragment insertion weight, vs. chunk insertion + small gap fragments
	core::Real big_frag_insertion_weight_;
	bool add_non_init_chunks_;
	bool domain_assembly_;
	bool add_hetatm_;
	core::Real hetatm_self_cst_weight_, hetatm_prot_cst_weight_;
	core::Real frag_weight_aligned_; // fragment insertion to the aligned region, vs. unaligned region
	bool auto_frag_insertion_weight_; // automatically set the fragment insertion weight
	core::Size max_registry_shift_;
	std::string cst_file_;

	core::Size initial_template_index_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1 < core::Real > template_wts_;
	utility::vector1 < core::pose::PoseOP > template_poses_;
	utility::vector1 < protocols::loops::Loops > template_chunks_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;
	utility::vector1 < core::fragment::FragSetOP > frag_libs_1mer_;
	utility::vector1 < core::fragment::FragSetOP > frag_libs_small_;
	utility::vector1 < core::fragment::FragSetOP > frag_libs_big_;

	Loops ss_chunks_pose_;

	core::Size top_n_big_frag_;
	core::Size top_n_small_frag_;

	// strand pairings
	std::string pairings_file_;
	utility::vector1< core::Size > sheets_;
	utility::vector1< core::Size > random_sheets_;
	bool filter_templates_;
	std::string target_sequence_;
	utility::vector1< std::pair< core::Size, core::Size > > strand_pairs_;
	std::set< core::Size > strand_pairings_template_indices_;
	std::set< core::Size > templates_with_incorrect_strand_pairings_;
	protocols::jumping::JumpSample jump_sample_;
	core::fragment::FragSetOP jump_frags_;
	std::set< core::Size > floating_pairs_;

	bool overlap_chainbreaks_;

	HybridizeFoldtreeDynamic foldtree_mover_;

	core::pose::PoseOP native_;

	utility::vector1<bool> allowed_to_move_;

}; //class FoldTreeHybridize

} // hybridization
} // protocols

#endif
