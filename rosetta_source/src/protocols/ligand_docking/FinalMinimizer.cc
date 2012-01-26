	// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com), adapted from the ResfileReader code
/// by Steven Lewis (smlewi@unc.edu) and Andrew Leaver-Fay

// Unit Headers
#include <protocols/ligand_docking/FinalMinimizer.hh>
#include <protocols/ligand_docking/FinalMinimizerCreator.hh>
#include <protocols/ligand_docking/MoveMapBuilder.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/MinimizeBackbone.hh>

//Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/ScoreFunction.hh>

// Scripter Headers
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>

// Utility Headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

//STL headers

namespace protocols {
namespace ligand_docking {

static basic::Tracer FinalMinimizer_tracer("protocols.ligand_docking.ligand_options.FinalMinimizer", basic::t_debug);

std::string
FinalMinimizerCreator::keyname() const
{
	return FinalMinimizerCreator::mover_name();
}

protocols::moves::MoverOP
FinalMinimizerCreator::create_mover() const {
	return new FinalMinimizer;
}

std::string
FinalMinimizerCreator::mover_name()
{
	return "FinalMinimizer";
}

FinalMinimizer::FinalMinimizer():
		Mover("FinalMinimizer"),
		score_fxn_(NULL),
		movemap_builder_(NULL)
		{}

FinalMinimizer::FinalMinimizer(FinalMinimizer const & that):
	    //utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		score_fxn_(that.score_fxn_),
		movemap_builder_(that.movemap_builder_)
{}

FinalMinimizer::~FinalMinimizer() {}

protocols::moves::MoverOP FinalMinimizer::clone() const {
	return new FinalMinimizer( *this );
}

protocols::moves::MoverOP FinalMinimizer::fresh_instance() const {
	return new FinalMinimizer;
}

std::string FinalMinimizer::get_name() const{
	return "FinalMinimizer";
}

//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
FinalMinimizer::parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & datamap,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "FinalMinimizer" ) utility_exit_with_message("This should be impossible");

	/// Score Function ///
	if ( ! tag->hasOption("scorefxn") ) utility_exit_with_message("'FinalMinimizer' requires 'scorefxn' tag");
	std::string scorefxn_name= tag->getOption<std::string>("scorefxn");
	score_fxn_= datamap.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name);

	/// MoveMapBuilder///
	if ( ! tag->hasOption("movemap_builder") ) utility_exit_with_message("'FinalMinimizer' requires 'movemap_builder' tag");
	std::string movemap_builder_name= tag->getOption<std::string>("movemap_builder");
	movemap_builder_= datamap.get< protocols::ligand_docking::MoveMapBuilder * >( "movemap_builders", movemap_builder_name);

}

void
FinalMinimizer::apply( core::pose::Pose & pose ){
	(*score_fxn_)(pose); // Debug Line, Remove Later...
	assert(movemap_builder_);

	if(movemap_builder_->minimize_backbone()){
		core::kinematics::FoldTree fold_tree_copy;
		FinalMinimizer_tracer<< "setting up fold_tree for min_bb"<< std::endl;
		fold_tree_copy= pose.fold_tree();
		InterfaceBuilderOP bb_interface_builder= movemap_builder_->get_bb_interface_builder();
		MinimizeBackbone backbone_foldtree_setup(bb_interface_builder);
		backbone_foldtree_setup.apply(pose);

		protocols::simple_moves::MinMoverOP const dfpMinTightTol = get_final_min_mover(pose);
		dfpMinTightTol->min_options()->nblist_auto_update(true);
		dfpMinTightTol->apply(pose);
		pose.fold_tree(fold_tree_copy);
	}
	else{
		protocols::simple_moves::MinMoverOP const dfpMinTightTol = get_final_min_mover(pose);
		dfpMinTightTol->min_options()->nblist_auto_update(true);
		dfpMinTightTol->apply(pose);
	}

}

protocols::simple_moves::MinMoverOP const
FinalMinimizer::get_final_min_mover(core::pose::Pose const & pose) const{
	std::string min_type= "dfpmin_armijo_nonmonotone_atol";
	core::Real tolerance= 0.02;
	bool use_nb_list=true;
	core::kinematics::MoveMapOP movemap= movemap_builder_->build(pose);
	movemap->show(FinalMinimizer_tracer, pose.n_residue());
	FinalMinimizer_tracer<< std::endl;
	return new protocols::simple_moves::MinMover(movemap, score_fxn_, min_type, tolerance, use_nb_list);
}

} //namespace ligand_docking
} //namespace protocols
