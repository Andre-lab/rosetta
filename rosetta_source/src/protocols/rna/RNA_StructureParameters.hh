// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 1.1.2.1 $
//  $Date: 2005/11/07 21:05:35 $
//  $Author: rhiju $
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @author Rhiju Das

#ifndef INCLUDED_protocols_rna_RNA_StructureParameters_HH
#define INCLUDED_protocols_rna_RNA_StructureParameters_HH

#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <protocols/rna/RNA_JumpLibrary.fwd.hh>
#include <protocols/toolbox/AllowInsert.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>

// C++ Headers
#include <string>
#include <map>
#include <iostream>
#include <list>

#include <utility/vector1.hh>



namespace protocols {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////
class RNA_Pairing {
public:
	core::Size pos1;
	core::Size pos2;
	char edge1; //W,H,S
	char edge2; //W,H,S
	char orientation; //A,P
};

typedef utility::vector1< RNA_Pairing > RNA_PairingList;

//////////////////////////////////////////////////////////////////////////////////////////////
class RNA_StructureParameters : public utility::pointer::ReferenceCount {
public:

	//constructor
	RNA_StructureParameters();
	~RNA_StructureParameters();
	void
	initialize(
		core::pose::Pose & pose,
		std::string const rna_params_file,
		std::string const jump_library_file,
		bool const ignore_secstruct );

	void
	setup_fold_tree_and_jumps_and_variants( core::pose::Pose & pose );

	bool
	random_jump_change( core::pose::Pose & pose ) const;

toolbox::AllowInsertOP
	allow_insert();

	void
	set_allow_insert(toolbox::AllowInsertOP allow_insert );

	bool
	check_base_pairs( core::pose::Pose & pose ) const;

	std::map< Size, Size >
	connections() const;

	std::list< Size >
	get_stem_residues(  core::pose::Pose const & pose ) const;

	void
	setup_base_pair_constraints( core::pose::Pose & pose );

	private:

	void
	initialize_secstruct( core::pose::Pose & pose );

	void
	override_secstruct( core::pose::Pose & pose );

	void
	append_virtual_anchor( core::pose::Pose & pose );

	void
	initialize_allow_insert( core::pose::Pose & pose  );

	void
	get_pairings_from_line(
												 std::istringstream & line_stream,
												 bool const in_stem );

	void
	read_chain_connection( std::istringstream & line_stream );

	void
	read_parameters_from_file( std::string const & pairing_file );


	void
	setup_jumps( core::pose::Pose & pose );

	void
	setup_chainbreak_variants( core::pose::Pose & pose );

	void
	setup_virtual_phosphate_variants( core::pose::Pose & pose );

	void
	setup_virtual_phosphate_variants_OLD( core::pose::Pose & pose ); //remove this by end of 2012

	void
	set_jump_library( RNA_JumpLibraryOP rna_jump_library );

	std::string const
	read_secstruct_from_file( std::string const & rna_secstruct_file );

	core::Size
	check_in_chain_connections( core::Size const & pos1, core::Size const & pos2 ) const;

	bool
	check_forward_backward(
		core::pose::Pose & pose,
		core::Size const jump_pos ) const;

	void
	add_new_RNA_jump(
		core::pose::Pose & pose,
		core::Size const & which_jump,
		bool & success ) const;

	void
	sample_alternative_chain_connection( core::pose::Pose & pose, core::Size const & which_jump ) const;

	void
	insert_base_pair_jumps( core::pose::Pose & pose, bool & success ) const;

	void
	fill_in_default_jump_atoms( core::kinematics::FoldTree & f, core::pose::Pose const & pose ) const;

private:
	RNA_JumpLibraryOP rna_jump_library_;
	RNA_PairingList rna_pairing_list_;

	utility::vector1 < utility::vector1 <core::Size > > obligate_pairing_sets_;
	utility::vector1 < utility::vector1 <core::Size > > possible_pairing_sets_;

	utility::vector1 < std::pair< utility::vector1 <core::Size >, utility::vector1 <core::Size > > > chain_connections_;

	//		int force_stems_;// deprecated.

	utility::vector1 <core::Size > cutpoints_open_;
	utility::vector1 <core::Size > cutpoints_closed_;
	utility::vector1 <core::Size > virtual_anchor_attachment_points_;

	bool secstruct_defined_;
	std::string rna_secstruct_;
	bool assume_non_stem_is_loop;

	bool add_virtual_anchor_;

	utility::vector1 < core::Size  > allow_insert_res_;
	toolbox::AllowInsertOP allow_insert_;

};

	typedef utility::pointer::owning_ptr< RNA_StructureParameters > RNA_StructureParametersOP;

} //rna
} //protocols

#endif
