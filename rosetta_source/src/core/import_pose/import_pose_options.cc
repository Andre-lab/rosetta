// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/import_pose/import_pose_options.cc
///
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Unit headers
#include <core/import_pose/import_pose_options.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// C++ headers

namespace core {
namespace import_pose {

basic::Tracer tr("core.import_pose.import_pose_options");

ImportPoseOptions::ImportPoseOptions() { init_from_options(); }

ImportPoseOptions::~ImportPoseOptions() {}

std::string ImportPoseOptions::type() const { return "ImportPoseOptions"; }

void ImportPoseOptions::parse_my_tag( utility::tag::TagPtr tag )
{
	PDB_DReaderOptions::parse_my_tag( tag );
	
	set_centroid( tag->getOption< bool >( "centroid", 0 ));
	set_fold_tree_io( tag->getOption< bool >( "fold_tree_io", 0 ));
	set_no_optH( tag->getOption< bool >( "no_optH", 0 ));
	set_pack_missing_sidechains( tag->getOption< bool >( "pack_missing_sidechains", 1 ));
	set_read_fold_tree( tag->getOption< bool >( "read_fold_tree", 0 ));
	set_rna( tag->getOption< bool >( "rna", 0 ));
	set_skip_set_reasonable_fold_tree_( tag->getOption< bool >( "skip_set_reasonable_fold_tree", 0 ));
	set_residue_type_set( tag->getOption< std::string >( "residue_type_set", "fa_standard" ));
}

// accessors
bool ImportPoseOptions::centroid() const { return centroid_; }
bool ImportPoseOptions::fold_tree_io() const { return fold_tree_io_; }
bool ImportPoseOptions::no_optH() const { return no_optH_; }
bool ImportPoseOptions::pack_missing_sidechains() const { return pack_missing_sidechains_; }
bool ImportPoseOptions::read_fold_tree() const { return read_fold_tree_; }
bool ImportPoseOptions::rna() const { return rna_; }
bool ImportPoseOptions::skip_set_reasonable_fold_tree() const { return skip_set_reasonable_fold_tree_; }

std::string const & ImportPoseOptions::residue_type_set() const { return residue_type_set_; }

// mutators
void ImportPoseOptions::set_centroid( bool centroid ) { centroid_ = centroid; }
void ImportPoseOptions::set_fold_tree_io( bool fold_tree_io ) { fold_tree_io_ = fold_tree_io; }
void ImportPoseOptions::set_no_optH( bool no_optH ) { no_optH_ = no_optH; }
void ImportPoseOptions::set_pack_missing_sidechains( bool pack_missing_sidechains ) { pack_missing_sidechains_ = pack_missing_sidechains; }
void ImportPoseOptions::set_read_fold_tree( bool read_fold_tree ) { read_fold_tree_ = read_fold_tree; }
void ImportPoseOptions::set_rna( bool rna ) { rna_ = rna; }
void ImportPoseOptions::set_skip_set_reasonable_fold_tree_( bool skip_set_reasonable_fold_tree ) { skip_set_reasonable_fold_tree_ = skip_set_reasonable_fold_tree; }

void ImportPoseOptions::set_residue_type_set( std::string const & residue_type_set ) { residue_type_set_ = residue_type_set; }


void ImportPoseOptions::init_from_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	set_centroid( option[ in::file::centroid_input ]()
		|| option[ in::file::centroid ]()
		|| ( option[ in::file::fullatom ].user() && !option[ in::file::fullatom ]())
		|| ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set ]() == "centroid" ));
	
	// sanity check
	if ( centroid() &&
		( option[ in::file::fullatom ]()
			|| ( option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set ]() == "fa_standard" ))) {
		tr.Warning << "conflicting command line flags for centroid/full-atom input. Choosing fullatom!" << std::endl;
		set_centroid( false );
	}
	set_fold_tree_io( option[ inout::fold_tree_io ].user());
	set_no_optH( option[ packing::no_optH ]());
	set_pack_missing_sidechains( option[ packing::pack_missing_sidechains ].value());
	set_read_fold_tree( false ); // no option for this parameter - it can only be set to true if you call pose_from_pdd.
	set_rna(option[ in::file::residue_type_set ].user() && option[ in::file::residue_type_set]()  == "rna");
	set_skip_set_reasonable_fold_tree_( option[ run::skip_set_reasonable_fold_tree ].value());
	
	set_residue_type_set( option[ in::file::residue_type_set ]());
}

} // namespace import_pose
} // namespace core
