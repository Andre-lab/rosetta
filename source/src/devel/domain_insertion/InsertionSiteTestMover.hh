// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/domain_insertion/InsertionSiteTestMover.hh
/// @brief  hh file for InsertionSiteTestMover
/// @author Florian Richter, flosopher@gmail.com, february 2013

#ifndef INCLUDED_devel_domain_insertion_InsertionSiteTestMover_hh
#define INCLUDED_devel_domain_insertion_InsertionSiteTestMover_hh

// Unit headers
#include <devel/domain_insertion/InsertionSiteTestMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// package headers

// Project headers
#include <core/types.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <devel/enzdes/EnzdesRemodelProtocol.hh>

// Utility headers
//#include <utility/VirtualBase.hh>

// C++ headers

//#include <utility/vector1.fwd.hh>


namespace devel {
namespace domain_insertion {


/// @brief  mover to test whether a pose wouldd tolerate
///         insertion of another domain at certain positions.
///         strategy: at all test positions, insert a short
///         stretch of glycines (~8) using rosetta remodel
///         make sure the pose doesn't have horrific energy,
///         minimize/relax, and check rmsd to prerelaxed pose.
///         if the energy doesn't get worse and the pose
///         doesn't move too much, the site might be able to
///         tolerate larger insertions.
///         I guess this could also be a filter, but it changes
///         the pose, so let's have it be a mover for now
class InsertionSiteTestMover : public protocols::moves::Mover {

public:

	typedef protocols::moves::Mover parent;

	typedef core::Size Size;
	typedef core::Real Real;

	InsertionSiteTestMover();

	InsertionSiteTestMover( InsertionSiteTestMover const & other );

	~InsertionSiteTestMover() override;

	protocols::moves::MoverOP
	clone() const override;




	void
	apply( core::pose::Pose & pose ) override;

	/// @brief Adds the given ResidueSelector to the selection. (As in all residues in both selections.)
	void add_insert_test_pos_selector( core::select::residue_selector::ResidueSelectorOP const & sele  );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		core::pose::Pose const & pose ) override;

	core::pack::task::PackerTaskCOP
	make_insert_task(
		core::pose::Pose const & pose,
		core::Size insert_pos )
	const;

	devel::enzdes::EnzdesRemodelMoverOP
	make_enzremodel_mover(
		core::pose::Pose & pose,
		core::Size insert_pos
	) const;

	/// @brief returns false if remodel mover didn't work
	bool
	create_raw_insert_pose(
		core::pose::Pose & pose,
		core::Size const insert_pos
	);

	void
	relax_raw_insert_pose(
		core::pose::Pose & pose,
		core::pose::Pose const & raw_pose,
		core::Size const insert_pos
	);

	void
	evaluate_insert_pose(
		core::pose::Pose const & start_pose,
		core::pose::Pose const & rawinsert_pose,
		core::pose::Pose const & relax_pose,
		core::Size const insert_pos
	);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	/*
	void
	output_eval_results(
	core::Size const insert_pos,
	core::Real const insert_pos_sasa,
	core::Real const diffsco_raw_start,
	core::Real const diffsco_relax_start,
	core::Real const diffsco_anchor,
	core::Real const rms_neighbors,
	core::Real const rms_anchors,
	core::Real const insert_sasa,
	core::Real const insert_burial
	);
	*/

private:

	core::scoring::ScoreFunctionOP sfxn_;
	core::select::residue_selector::ResidueSelectorOP insert_test_pos_;
	core::Size flex_window_; //how many res up- and downstream of the insert pos to move

	std::string test_insert_ss_;
	Real insert_allowed_score_increase_, insert_attempt_sasa_cutoff_;
	core::Size length_of_insert_, num_repeats_;
	bool pdb_numbering_;

	//actual insertion handled by enzdes machinery
	protocols::enzdes::EnzdesFlexBBProtocolOP enz_flexbb_prot_;

	core::id::SequenceMappingCOP insert_seqmap_;

};


}
}

#endif
