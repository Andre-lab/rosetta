// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/Pose.hh
/// @brief  Pose class
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov

/// Be very careful with #includes in this file.  As almost every file in mini
/// #includes Pose.hh, any header file #included here will produce a huge
/// recompile when modified.  I pledge a personally delivered, hand crafted
/// ass whooping for any fool that #includes an unneccessary .hh file in Pose.hh.
/// Allowed .hh files are few: Residue, ResidueType, Conformation, Energies
/// AtomID and NamedAtomID.  Everything else should be forward included, or
/// not included at all.  The above pledge applies to any header file
/// unnecessarily #included in the .hh files #included here.
/// If Pose were to hold a ConformationOP and an EnergiesOP, and not inline
/// their access, then Residue.hh, ResidueType.hh, Conformation.hh and
/// Energies.hh could be replaced with .fwd.hh's instead.


#ifndef INCLUDED_core_pose_Pose_hh
#define INCLUDED_core_pose_Pose_hh


// type headers
#include <core/types.hh>

// Unit headers
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>

// AUTO-REMOVED #include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
// AUTO-REMOVED #include <core/id/NamedStubID.hh>
#include <core/id/TorsionID.fwd.hh>

#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>

#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/MetricValue.fwd.hh>

// //#include "cst_set.h"
// #include "jump_classes.h"
// #include "pose_param.h"
// #include "score_data.h"

// // ObjexxFCL Headers
// #include <ObjexxFCL/FArray3A.hh>
// #include <ObjexxFCL/FArray4D.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/vector1.fwd.hh>

// Numeric headers
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>

#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <utility/vector1.hh>


#ifdef WIN32
	#include <core/pose/signals/ConformationEvent.hh>
	#include <core/pose/signals/DestructionEvent.hh>
	#include <core/pose/signals/EnergyEvent.hh>
	#include <core/pose/signals/GeneralEvent.hh>
	#include <core/conformation/Residue.hh>
	#include <core/pose/datacache/CacheableObserver.hh>
	#include <core/id/AtomID.hh>
#endif


// C++ Headers
// #include <cassert>
// #include <cmath>
// #include <cstdlib>
// #include <iostream>
// #include <map>
// #include <string>


namespace core {
namespace pose {


	/// A molecular system including residues, kinematics, and energies

/**

  The Pose class represents a molecular system (protein-dna-ligand...)
as a container of Rosetta Residue objects together with
a Conformation object that defines how internal coordinate changes
propagate through the system and an Energies object that stores
information from the last energy evaluation.


The main responsibilities of the pose are:

@li  Kinematic:
  (a) to update the xyz coordinates in response to changes to internal
      degrees of freedom, and
  (b) to update internal coordinates when the user modifes the xyz
      (Cartesian) coords,

@li  Scoring:
  (a) to keep track of what parts of the structure have changed since
      the last score evaluation, and
  (b) to cache residue and residue-pair energies for efficient re-use

@li As a container:
  The pose provides a single object for passing
  a molecular system and for copying of entire molecules
  or stretches of molecules from one Pose object into another.


Output Methods:
Common Methods:
    Pose.assign
    Pose.atom_tree
    Pose.conformation
    Pose.dump_pdb
    Pose.energies
    Pose.fold_tree
    Pose.pdb_info
    Pose.residue
    Pose.sequence
    Pose.total_residue
**/

class Pose : public utility::pointer::ReferenceCount {

public:
	typedef id::AtomID AtomID;
	typedef id::NamedAtomID NamedAtomID;
	typedef id::TorsionID TorsionID;
	typedef id::DOF_ID DOF_ID;
	typedef conformation::Residue Residue;
	typedef conformation::Conformation Conformation;
	typedef conformation::ConformationOP ConformationOP;
	typedef pose::datacache::ObserverCache ObserverCache;
	typedef pose::datacache::ObserverCacheOP ObserverCacheOP;
	typedef pose::signals::ConformationEvent ConformationEvent;
	typedef pose::signals::DestructionEvent DestructionEvent;
	typedef pose::signals::EnergyEvent EnergyEvent;
	typedef pose::signals::GeneralEvent GeneralEvent;
	typedef scoring::constraints::ConstraintSetOP ConstraintSetOP;
	typedef scoring::constraints::ConstraintSetCOP ConstraintSetCOP;
	typedef basic::datacache::BasicDataCache BasicDataCache;
	typedef basic::datacache::BasicDataCacheOP BasicDataCacheOP;

public:

	/// @brief default constructor, builds an empty pose
    ///
    /// AtomTree      default   /bonding information
    /// Conformation  default   /change propagation
    /// Energies      default   /contains pose energy information
    /// FoldTree      default   /directs building of Pose
    /// Residue       default   /container of individual Residue objects
	Pose();

	/// @brief destructor -- > kill data on the heap
	~Pose();

	/// @brief copy constructor
	Pose( Pose const & src );

	/// @brief partial copy constructor
	Pose( Pose const & src, Size residue_begin, Size residue_end);

	/// @brief Construct pose from pdb file
	//Pose( std::string const & pdb_file );

	/// @brief Copies  <src>  into the pose
	///
	/// example(s):
	///     test_pose.assign(pose)
	/// See also:
	///     Pose
	Pose &
	operator=( Pose const & src );

	////////////////////////////////////////
	// tree builders / modifiers

	/// @brief Returns the pose Conformation (const-access)
    ///
    /// example(s):
    ///     pose.Conformation()
    /// See also:
    ///     Pose
    ///     Conformation
	Conformation const &
	conformation() const
	{
		return *conformation_;
	}

	/// @brief Returns the pose Conformation (non-const access)
	Conformation &
	conformation()
	{
		return *conformation_;
	}

	/// @brief Returns the pose FoldTree
	///
	/// example(s):
	///     pose.fold_tree()
	/// See also:
	///     Pose
	///     FoldTree
	kinematics::FoldTree const &
	fold_tree() const;

	/// @brief Sets the pose FoldTree to  <fold_tree_in>
	///
	/// example(s):
	///     pose.fold_tree( foldtree )
	/// See also:
	///     Pose
	///     pose.fold_tree
	///     FoldTree
	void
	fold_tree( kinematics::FoldTree const & fold_tree_in );

  /// @brief  Now that the conformation_ member data is an owning pointer,
	/// and we have derived classes of Conformation.
  void
  set_new_conformation( ConformationOP new_conformation );

	/// @brief  Now that the energies_ member data is an owning pointer,
  /// and we have derived classes of Energies.
  void
  set_new_energies_object( scoring::EnergiesOP energies );

	/// @brief Returns the pose AtomTree
	kinematics::AtomTree const &
	atom_tree() const;

	/// @brief Returns the chain number of residue  <seqpos>
	///
	/// example(s):
	///     pose.chain(3)
	/// See also:
	///     Pose
	///     Pose.annotated_sequence
	///	Pose.chain_sequence
	///     Pose.fold_tree
	///     Pose.sequence
	///     FoldTree
	int
	chain( Size const seqpos ) const;

	/// @brief Returns a vector of poses with one element per chain of the original pose
	utility::vector1< PoseOP >
	split_by_chain() const;

	/// @brief
	Pose
	split_by_chain(Size chain_id) const;

	/// @brief  Updates the pose chain IDs to match the pdb chain IDs.
	void
	update_pose_chains_from_pdb_chains();

	/// APL Removing illegal non-const accessors to residues which
	/// otherwise violate the data-integrity guarantees provided by
	/// class Conformation.
	/// accessors for iteration over residues
	/// THIS MUST BE REMOVED ASAP!
	//conformation::ResidueOPs::iterator res_begin();
	//conformation::ResidueOPs::iterator res_end  ();

	///////////////////////////////////////////
	/// @brief Returns the pose Energies (const-access)
	///
	/// example(s):
	///     pose.energies()
	/// See also:
	///     Pose
	///     Energies
	///     PDBInfo
	///     ScoreFunction
	///     create_score_function
	scoring::Energies const &
	energies() const {
		return *energies_;
	}

	/// @brief Returns the pose Energies (non-const access)
	scoring::Energies &
	energies() {
		return *energies_;
	}

	ConstraintSetCOP
	constraint_set() const;

	/// @brief adding a constraint is done by cloning the input constraint. A const copy is then returned
	scoring::constraints::ConstraintCOP
	add_constraint( scoring::constraints::ConstraintCOP cst );

	scoring::constraints::ConstraintCOPs
	add_constraints( scoring::constraints::ConstraintCOPs csts );

	/// @brief re object_comparison see comment for ConstraintSet::remove_constraint function
	bool
	remove_constraint(
		scoring::constraints::ConstraintCOP cst,
		bool object_comparison = false
	);

	/// @brief re object_comparison see comment for ConstraintSet::remove_constraint function
	bool
	remove_constraints(
		scoring::constraints::ConstraintCOPs csts,
		bool object_comparison = false
	);

	bool
	remove_constraints();

	void
	constraint_set( ConstraintSetOP );

	void transfer_constraint_set( const pose::Pose &pose );

	/// @brief Returns the pose PDBInfo (const)
	///
	/// example(s):
	///     pose.pdb_info()
	/// See also:
	///     Pose
	///     Energies
	///     PDBInfo
	///	ScoreFunction
	///	pose_from_pdb
	/// @return NULL if no PDBInfo instance exists, the pdb info instance otherwise
	PDBInfoCOP
	pdb_info() const;

	/// @brief Returns the pose PDBInfo
	/// @return NULL if no PDBInfo instance exists, the PDBInfo instance otherwise
	PDBInfoOP
	pdb_info();

	/// @brief Sets pose PDBInfo to <new_info>
	/// @param[in]  <new_info>  the new PDBInfo to copy, pass NULL
	/// if you want to zero the existence of PDBInfo inside this Pose
	/// @return the prior PDBInfo instance
	PDBInfoOP
	pdb_info( PDBInfoOP new_info );

	void
	metric( std::string const & calculator_name, std::string const & key, basic::MetricValueBase & val ) const;

	std::string
	print_metric( std::string const & calculator_name, std::string const & key ) const;

	////////////////////////////////////////
	// miscellaneous structure modification:

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// insert/append/delete residues
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	/// @brief Appends  <new_rsd>  (a residue) to pose by a new jump
	void
	append_residue_by_jump(
		conformation::Residue const & new_rsd,
		Size const jump_anchor_residue,
		std::string const& jump_anchor_atom = "",
		std::string const& jump_root_atom = "",
		bool const start_new_chain = false
	);

	/// @brief Appends  <new_rsd>  (a residue) to pose by a new bond
	///
	/// @details The default behavior is to append by a polymeric connection to the preceding residue
	/// If we want to connect via a non-polymer connection, we give the connection number, anchor residue
	/// and the connection number for the anchor residue. These connection numbers are wrt the connections_
	/// arrays in Residue and ResidueType
	///
	/// If build_ideal_bond is TRUE it will transform the coordinates of the new residue so that the bond
	/// geometry of the new bond is ideal according to the icoor_internal data in the residues.
	///
	/// Otherwise the incoming coordinates of new_rsd are preserved.
	void
	append_residue_by_bond(
		conformation::Residue const & new_rsd,
		bool const build_ideal_geometry = false,
		int const connection = 0,
		Size const anchor_residue = 0,
		int const anchor_connection = 0,
		bool const start_new_chain = false,
		bool const lookup_bond_length = false
	);

	///
	/// This code sorely belongs in Pose.cc
	/// @brief Adds  <new_rsd_in>  to pose at  <seqpos>
	void
	insert_residue_by_jump(
		Residue const & new_rsd_in,
		Size const seqpos, // desired seqpos of new_rsd
		Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
		std::string const& anchor_atomno = "",
		std::string const& root_atomno = ""
	);

	/// @brief Replaces the residue at  <seqpos>  with  <new_rsd_in>
	void
	replace_residue(
		Size const seqpos,
		Residue const & new_rsd_in,
		bool const orient_backbone
	);

	/// @brief Replaces the residue at  <seqpos>  with  <new_rsd>
	/// based on superposition on the specified input atom pairs
	/// NOTE: at the moment, only superposition on 3 atoms works
	/// This code sorely belongs in Pose.cc
	void
	replace_residue(
		int const seqpos,
		Residue const & new_rsd_in,
		utility::vector1< std::pair< std::string, std::string > > const & atom_pairs
	);

	/// @brief glues to seqpos and perhaps also seqpos+1
	void
	append_polymer_residue_after_seqpos(
		Residue const & new_rsd,
		Size const seqpos,
		bool const build_ideal_geometry
	);

	/// @brief glues to seqpos and perhaps also seqpos-1
	void
	prepend_polymer_residue_before_seqpos(
		Residue const & new_rsd,
		Size const seqpos,
		bool const build_ideal_geometry
	);

	void
	delete_polymer_residue( Size const seqpos );

	/// @brief Copy a stretch of coordinates/torsions from  <src>
	/// to pose
	void
	copy_segment(
		Size const size,
		Pose const & src,
		Size const begin,
		Size const src_begin
	);

	///////////////////////
	// miscellaneous access

	/// @brief BasicDataCache indexed by enum in core/pose/datacache/CacheableDataType.hh
	BasicDataCache const &
	data() const
	{
		return *data_cache_;
	}

	/// @brief BasicDataCache indexed by enum in core/pose/datacache/CacheableDataType.hh
	BasicDataCache &
	data()
	{
		return *data_cache_;
	}


	/// @brief ObserverCache indexed by enum in core/pose/datacache/CacheableObserverType.hh
	ObserverCache const &
	observer_cache() const
	{
		return *observer_cache_;
	}


	/// @brief ObserverCache indexed by enum in core/pose/datacache/CacheableObserverType.hh
	ObserverCache &
	observer_cache()
	{
		return *observer_cache_;
	}

	/// @brief Returns the number of residues in the pose conformation
	///
	/// example(s):
	///     pose.total_residue()
	/// See also:
	///     Pose
	///     Pose.n_residue
	///     Pose.sequence
	Size
	total_residue() const;

	/// @brief Returns the number of residues in the pose conformation
	/// example(s):
	///     pose.n_residue()
	/// See also:
	///     Pose
	///     Pose.sequence
	///     Pose.total_residue
	Size
	n_residue() const;

	/// @brief Returns true if there are no residues in the conformation
	///
	/// example(s):
	///     pose.empty()
	/// See also:
	///     Pose
	///     Pose.clear
	///	Pose.sequence
	///     Pose.total_residue
	bool
	empty() const;

	/// @brief Returns the number of jumps in the pose FoldTree
	///
	/// example(s):
	///     pose.num_jump()
	/// See also:
	///     Pose
	///     Pose.jump
	///     Pose.set_jump
	///     Pose.set_jump_now
	///     FoldTree
	///     Jump
	Size
	num_jump() const;

	/// @brief Returns the chemical::AA of the residue at  <seqpos>
	///
	/// example(s):
	///     pose.aa(17)
	/// See also:
	///     Pose
	///	Pose.Residue
	///     Pose.sequence
	///     Residue
	chemical::AA
	aa( Size const seqpos) const;

	/// @brief Returns the secondary structure of residue  <seqpos>
	/// this usually comes from fragments. The conformation object
	/// will not invoke DSSP to determine the secondary structure
	/// if e.g. it has not been made from fragments.
	/// 'H' = helical
	/// 'S' = strand or sheet
	/// 'E' = loop
	///
	/// example(s):
	///     pose.secstruct(3)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_secstruct
	///     Residue
	char
	secstruct( Size const seqpos ) const;

	/// @brief Returns a string representing pose secondary structure
	///
	/// example(s):
	///     pose.secstruct()
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_secstruct
	///     Residue
	std::string
	secstruct() const;

	/// @brief Assign the secondary structure of residue  <seqpos> to  <setting>
	///
	/// example(s):
	///     pose.set_secstruct(3,'H')
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.secstruct
	///     Residue
	void
	set_secstruct( Size const seqpos, char const setting );

	/// @brief Returns a string representing the 1-leter-coded
	/// sequence of the pose conformation
	///
	/// example(s):
	///     pose.sequence()
	/// See also:
	///     Pose
	///     Pose.annotated_sequence
	///     Pose.chain
	///     Pose.chain_sequence
	///     Pose.residue
	///     Pose.total_residue
	std::string
	sequence() const;

	/// @brief Returns the variant-tagged string representing the
	/// residue types that make up a conformation; e.g.
	/// M[MET:N-Terminus-Variant]CDH[HIS_D]LLR[ARG:C-Terminus-Variant]
	///
	/// example(s):
	///     pose.annotated_sequence()
	/// See also:
	///     Pose
	///     Pose.sequence
	///     Pose.total_residue
	///	Residue
	std::string
	annotated_sequence( bool show_all_variants = false ) const;

	/// @brief Returns the sequence for the chain  <chain_in>
	///
	/// example(s):
	///     pose.chain_sequence(1)
	/// See also:
	///     Pose
	///     Pose.chain
	///     Pose.residue
	///     Pose.sequence
	std::string
	chain_sequence( core::Size const chain_in ) const;

	/// @brief Returns the Residue at position  <seqpos>  (read access)
	/// Note: this method will trigger a refold if either the
	/// torsions or the coordinates are out-of-date
	///
	/// example(s):
	///     pose.residue(4)
	/// See also:
	///     Pose
	///     Pose.sequence
	///     Pose.total_residue
	///     Residue
	///     ResidueType
	Residue const &
	residue(
		Size const seqpos
	) const;

	/// @brief Returns the ResidueType at position  <seqpos>  (read access)
	/// Note: this method NOT will trigger a refold if either
	/// the torsions or the coordinates are out-of-date
	///
	/// example(s):
	///     pose.residue_type(5)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.sequence
	///     Pose.total_residue
	///     Residue
	///     ResidueType
	chemical::ResidueType const &
	residue_type(
		Size const seqpos
	) const;

	/// @brief Returns true if pose is ResidueType fullatom
	///
	/// @note convenience test for residue_type_set ( based on two
	/// middle residue -- to avoid hitting on ligands or pseudos )
	///
	/// example(s):
	///     pose.is_fullatom()
	/// See also:
	///     Pose
	///     Pose.is_centroid
	///     Residue
	///     ResidueType
	///@brief this is nt a good test --Doug
	bool is_fullatom() const;

	/// @brief Returns true if pose is ResidueType centroid
	///
	/// @note convenience test for residue_type_set ( based on two
	/// middle residue -- to avoid hitting on ligands or pseudos )
	///
	/// example(s):
	///     pose.is_centroid()
	/// See also:
	///     Pose
	///     Pose.is_fullatom
	///     Residue
	///     ResidueType
	///@brief this is nt a good test --Doug
	bool is_centroid() const;

	///////////////////////////////////////////////////////////////////////////
	// convenience access for protein residues



	/// @brief Returns the phi torsion angle of residue  <seqpos>
	/// @note assumes the residue is an amino acid
	///
	/// example(s):
	///     pose.phi(1)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_phi
	///     Residue
	Real
	phi( Size const seqpos ) const;

	/// @brief Sets the phi torsion angle of residue  <seqpos> to  <setting>
	/// @note  <setting>  must be in degrees, assumes residue is an amino acid
	///
	/// example(s):
	///     pose.set_phi(1)
	/// See also:
	///     Pose
	///     Pose.phi
	///     Pose.residue
	///     Residue
	void
	set_phi( Size const seqpos, Real const setting );

	/// @brief Returns the psi torsion angle of residue  <seqpos>
	/// Note: assumes the residue is an amino acid
	///
	/// example(s):
	///     pose.psi(2)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_psi
	///     Residue
	Real
	psi( Size const seqpos ) const;

	/// @brief Sets the psi torsion angle of residue  <seqpos>  to  <setting>
	/// @note  <setting>  must be in degrees, assumes residue is an amino acid
	///
	/// example(s):
	///     pose.set_psi(2)
	/// See also:
	///     Pose
	///     Pose.psi
	///     Pose.residue
	///     Residue
	void
	set_psi( Size const seqpos, Real const setting );

	/// @brief Returns the omega torsion angle of residue  <seqpos>
	/// @note assumes the residue is an amino acid
	///
	/// example(s):
	///     pose.omega(3)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_omega
	///     Residue
	Real omega( Size const seqpos ) const;

	/// @brief Sets the omega torsion angle of residue  <seqpos>  to  <setting>
	/// @note  <setting>  must be in degrees, assumes residue is an amino acid
	///
	/// example(s):
	///     pose.set_omega(3)
	/// See also:
	///     Pose
	///     Pose.omega
	///     Pose.residue
	///     Residue
	void
	set_omega( Size const seqpos, Real const setting );

	///////////////////////////////////////////////////////////////////////////
	// convenience access for nucleic acid residues



	/// @brief Returns the alpha torsion angle of residue  <seqpos>
	/// @note assumes the residue is an nucleic acid
	///
	/// example(s):
	///     pose.alpha(1)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_alpha
	///     Residue
	Real
	alpha( Size const pos ) const;

	/// @brief Sets the alpha torsion angle of residue  <seqpos>  to  <setting>
	/// Note:  <setting>  must be in degrees, assumes residue is an nucleic acid
	///
	/// example(s):
	///     pose.set_alpha(1)
	/// See also:
	///     Pose
	///     Pose.alpha
	///     Pose.residue
	///     Residue
	void
	set_alpha( Size const seqpos, Real const setting );

	/// @brief Returns the beta torsion angle of residue  <seqpos>
	/// @note assumes the residue is an nucleic acid
	///
	/// example(s):
	///     pose.beta(2)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_beta
	///     Residue
	Real
	beta( Size const seqpos ) const;

	/// @brief Sets the beta torsion angle of residue  <seqpos>  to  <setting>
	/// Note:  <setting>  must be in degrees, assumes residue is an nucleic acid
	///
	/// example(s):
	///     pose.set_beta(2)
	/// See also:
	///     Pose
	///     Pose.beta
	///     Pose.residue
	///     Residue
	void
	set_beta( Size const seqpos, Real const setting );

	/// @brief Returns the gamma torsion angle of residue  <seqpos>
	/// @note assumes the residue is an nucleic acid
	///
	/// example(s):
	///     pose.gamma(3)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_gamma
	///     Residue
	Real gamma( Size const seqpos ) const;

	/// @brief Sets the gamma torsion angle of residue  <seqpos>  to  <setting>
	/// Note:  <setting>  must be in degrees, assumes residue is an nucleic acid
	///
	/// example(s):
	///     pose.set_gamma(3)
	/// See also:
	///     Pose
	///     Pose.gamma
	///     Pose.residue
	///     Residue
	void
	set_gamma( Size const seqpos, Real const setting );

	/// @brief Returns the delta torsion angle of residue  <seqpos>
	/// @note assumes the residue is an nucleic acid
	///
	/// example(s):
	///     pose.delta(4)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_delta
	///     Residue
	Real
	delta( Size const pos ) const;

	/// @brief Sets the delta torsion angle of residue  <seqpos>  to  <setting>
	/// Note:  <setting>  must be in degrees, assumes residue is an nucleic acid
	///
	/// example(s):
	///     pose.set_delta(4)
	/// See also:
	///     Pose
	///     Pose.delta
	///     Pose.residue
	///     Residue
	void
	set_delta( Size const seqpos, Real const setting );

	/// @brief Returns the epsilon torsion angle of residue  <seqpos>
	/// @note assumes the residue is an nucleic acid
	///
	/// example(s):
	///     pose.epsilon(5)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_epsilon
	///     Residue
	Real
	epsilon( Size const seqpos ) const;

	/// @brief Sets the epsilon torsion angle of residue  <seqpos>  to  <setting>
	/// Note:  <setting>  must be in degrees, assumes residue is an nucleic acid
	///
	/// example(s):
	///     pose.set_epsilon(5)
	/// See also:
	///     Pose
	///     Pose.epsilon
	///     Pose.residue
	///     Residue
	void
	set_epsilon( Size const seqpos, Real const setting );

	/// @brief Returns the zeta torsion angle of residue  <seqpos>
	/// @note assumes the residue is an nucleic acid
	///
	/// example(s):
	///     pose.zeta(6)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_zeta
	///     Residue
	Real zeta( Size const seqpos ) const;

	/// @brief Sets the zeta torsion angle of residue  <seqpos>  to  <setting>
	/// Note:  <setting>  must be in degrees, assumes residue is an nucleic acid
	///
	/// example(s):
	///     pose.set_zeta(6)
	/// See also:
	///     Pose
	///     Pose.zeta
	///     Pose.residue
	///     Residue
	void
	set_zeta( Size const seqpos, Real const setting );

	/// @brief Returns the chi torsion angle of residue  <seqpos>
	/// @note assumes the residue is an nucleic acid
	///
	/// example(s):
	///     pose.chi(7)
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.set_chi
	///     Residue
	Real chi( Size const seqpos ) const;

	/// @brief Sets the chi torsion angle of residue  <seqpos>  to  <setting>
	/// Note:  <setting>  must be in degrees, assumes residue is an nucleic acid
	///
	/// example(s):
	///     pose.set_chi(7)
	/// See also:
	///     Pose
	///     Pose.chi
	///     Pose.residue
	///     Residue
	void
	set_chi( Size const seqpos, Real const setting );

	/////////////////////////////////////////////////////////////////////////////
	// jumps



	/// @brief Sets the pose FoldTree Jump  <jump_number>  to  <new_jump>
	///
	/// example(s):
	///     pose.set_jump(1,jump1)
	/// See also:
	///     Pose
	///     Pose.fold_tree
	///     Pose.jump
	///     Pose.set_jump_now
	///     FoldTree
	///     FoldTree.jump_edge
	///     Jump
	void
	set_jump(
		int const jump_number,
		const kinematics::Jump & new_jump
	);

	/// @brief Sets the pose FoldTree Jump  <jump_number>  to  <new_jump>
	///
	/// example(s):
	///     pose.set_jump_now(2,jump2)
	/// See also:
	///     Pose
	///     Pose.fold_tree
	///     Pose.jump
	///     Pose.set_jump
	///     FoldTree
	///     FoldTree.jump_edge
	///     Jump
	void
	set_jump_now(
		int const jump_number,
		const kinematics::Jump & new_jump
	);

	/// @brief Returns the pose FoldTree Jump  <jump_number>
	///
	/// example(s):
	///     pose.jump(1)
	/// See also:
	///     Pose
	///     Pose.fold_tree
	///     Pose.set_jump
	///     FoldTree
	///     FoldTree.jump_edge
	///     Jump
	kinematics::Jump const &
	jump( int const jump_number ) const;

	/// @brief Sets the pose FoldTree Jump  <id>  to  <new_jump>
	///
	/// example(s):
	///     pose.set_jump(1,jump1)
	/// See also:
	///     Pose
	///     Pose.fold_tree
	///     Pose.jump
	///     FoldTree
	///     FoldTree.jump_edge
	///     AtomID
	void
	set_jump(
		AtomID const & id,
		const kinematics::Jump & new_jump
	);

	/// @brief Returns the pose FoldTree Jump  <id>
	///
	/// example(s):
	///     pose.set_jump(R5N)
	/// See also:
	///     Pose
	///     Pose.fold_tree
	///     Pose.set_jump
	///     FoldTree
	///     FoldTree.jump_edge
	///     AtomID
	kinematics::Jump const &
	jump( AtomID const & id ) const;

	/////////////////////////////////////////////////////////////////////////////
	// generic torsion-angle access



	/// @brief Returns the  <chino>  chi torsion angle of residue  <seqpos>
	/// @note assumes the residue is an amino acid
	///
	/// example(s):
	///     pose.chi(1,7)
	/// See also:
	///     Pose
	///     Pose.set_chi
	///     Pose.residue
	///     Residue
	Real
	chi(
		int const chino,
		Size const seqpos
	) const;

	/// @brief Sets the  <chino>  chi torsion angle of residue  <seqpos>  to  <setting>
	/// @note  <setting>  must be in degrees, assumes residue is an amino acid
	///
	/// example(s):
	///     pose.set_chi(1,7,120)
	/// See also:
	///     Pose
	///     Pose.chi
	///     Pose.residue
	///     Residue
	void
	set_chi(
		int const chino,
		Size const seqpos,
		Real const setting
	);

	/// @brief Returns the Conformation torsion angle identified by  <id>
	///
	/// See also:
	///     Pose
	///     TorsionID
	Real
	torsion( TorsionID const & id ) const;

	/// @brief Sets the Conformation torsion angle identified
	/// by  <id>  to  <setting>
	///
	/// See also:
	///     Pose
	///     TorsionID
	void
	set_torsion( TorsionID const & id, Real const setting );

	///////////////////////////////////////////////////////////////////////////
	// access atomtree dof's



	/// @brief Returns the value of the AtomTree DOF  <id>
	///
	/// See also:
	///     Pose
	///     DOF_ID
	Real
	dof( DOF_ID const & id ) const;

	/// @brief Sets the value of the AtomTree DOF  <id>
	///
	/// See also:
	///     Pose
	///     DOF_ID
	void
	set_dof( DOF_ID const & id, Real const setting );

	/// @brief Returns true if pose has DOF  <id>
	///
	/// See also:
	///     Pose
	///     DOF_ID
	bool
	has_dof( DOF_ID const & id ) const;

	/// @brief Returns the location (xyz) of pose AtomID  <id>
	///
	/// example(s):
	///     atom = AtomID(1,1)
	///     pose.xyz(atom)
	/// See also:
	///     Pose
	///     Pose.residue
	///     AtomID
	///     Residue
	///     Residue.xyz
	PointPosition const &
	xyz( AtomID const & id ) const;

	/// @brief Returns the location (xyz) of pose NamedAtomID  <id>
	///
	/// Tutorial soon...
	/// See also:
	///     Pose
	///     Pose.residue
	///     NamedAtomID
	///     Residue
	///     Residue.xyz
	PointPosition const &
	xyz( NamedAtomID const & id ) const;

	/// @brief Sets the location (xyz) of pose AtomID  <id>  to
	/// the PointPosition  <point>
	///
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.xyz
	///     Residue
	void
	set_xyz( AtomID const & id, PointPosition const & point );

	/// @brief Sets the location (xyz) of pose NamedAtomID  <id>
	/// to the PointPosition  <point>
	///
	/// See also:
	///     Pose
	///     Pose.residue
	///     Pose.xyz
	///     Residue
	void
	set_xyz( NamedAtomID const & id, PointPosition const & point );

	/// @brief Sets the locations (xyz) of pose AtomIDs in  <ids>
	/// to mathcing PointPositions in  <points>
	void
	batch_set_xyz( utility::vector1< AtomID > const & ids, utility::vector1< PointPosition > const & points );

	/// @brief Gets the locations (xyz) of pose AtomIDs in  <ids>
	void
	batch_get_xyz( utility::vector1< AtomID > const & ids, utility::vector1< PointPosition > & points ) const;

	kinematics::Stub
	stub_from_id( id::NamedStubID const& id );

	/// @brief Sets pose coordinates such that the pose center is at the Euclidean origin
	void
	center();

	/// @brief Updates neighbor links in the pose Energies object
	void
	update_residue_neighbors();


	/// @brief Called by ScoreFunction at the beginning of scoring
	void
	scoring_begin(
		scoring::ScoreFunction const & info
	);

	/// @brief Called by ScoreFunction at the end of scoring
	void
	scoring_end( scoring::ScoreFunction const & scorefxn );

	/// @brief Called by PairEPotential to update the action coordinates for all residues
	void
	update_actcoords();

	/// @brief Updates the action coordinates for pose residue  <resid>
	void
	update_actcoord( Size resid );

	void
	update_orbital_coords(Size resid);

	/// @brief Apply a transform of the Rx + v form, where R is a
	/// rotation matrix and v is a translation vector.
	void
	apply_transform_Rx_plus_v(
		numeric::xyzMatrix< Real > const & R,
		Vector const & v
	);

	/// @brief Empty the pose contents
	///
	/// example(s):
	///     pose.clear()
	/// See also:
	///     Pose
	///     Pose.assign
	///     Pose.empty
	void
	clear();

	/// @brief Export pose data to the PDB file  <file_name>
	///
	/// example(s):
	///     pose.dump_pdb('new_01.pdb')
	/// See also:
	///     Pose
	///     pose_from_pdb
	bool
	dump_pdb(std::string const & file_name, std::string const & tag="1") const;

	void dump_pdb(std::ostream & out, std::string const & tag="1") const;

	/// @brief for writing a specified subset of residues in pdb format
	void
	dump_pdb(
		std::ostream & out,
		utility::vector1< core::Size > const & residue_indices,
		std::string const & tag="1"
	) const;


	/// @brief Export pose data to the PDB file  <file_name>,
	/// add some score output
	void
	dump_scored_pdb( std::string const & file_name, scoring::ScoreFunction const & scorefxn, std::string const & tag="1" );


public: // observer attach/detach


	/// @brief attach DestructionEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( DestructionEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection.
	/// @remarks DestructionEvent observers will only be notified upon destruction
	///  of the Pose
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_destruction_obs( MemFn fn, Ptr ptr ) {
		return destruction_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach DestructionEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( DestructionEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	/// @remarks DestructionEvent observers will only be notified upon destruction
	///  of the Pose
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_destruction_obs( MemFn fn, Ptr ptr ) const {
		return destruction_obs_hub_.disconnect( fn, ptr );
	}


	/// @brief attach GeneralEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( GeneralEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection.
	/// @remarks GeneralEvent observers will be notified whenever any signal
	///  derived from GeneralEvent occurs.
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_general_obs( MemFn fn, Ptr ptr ) {
		return general_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach GeneralEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( GeneralEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	/// @remarks GeneralEvent observers will be notified whenever any signal
	///  derived from GeneralEvent occurs.
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_general_obs( MemFn fn, Ptr ptr ) const {
		return general_obs_hub_.disconnect( fn, ptr );
	}


	/// @brief attach EnergyEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( EnergyEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection.
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_energy_obs( MemFn fn, Ptr ptr ) {
		return energy_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach EnergyEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( EnergyEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_energy_obs( MemFn fn, Ptr ptr ) const {
		return energy_obs_hub_.disconnect( fn, ptr );
	}


	/// @brief attach ConformationEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( ConformationEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection.
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_conformation_obs( MemFn fn, Ptr ptr ) {
		return conformation_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach ConformationEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( ConformationEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_conformation_obs( MemFn fn, Ptr ptr ) const {
		return conformation_obs_hub_.disconnect( fn, ptr );
	}


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////// internal methods ////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
private: //
	/// @brief initialize internal for a freshly build pose object. Use this function when you need
	///        to create new constructor.
	void init(void);


private: // observer notifications


	/// @brief notify DestructionEvent observers
	/// @remarks called only upon destruction of the Pose
	void
	notify_destruction_obs( DestructionEvent const & e );


	/// @brief notify GeneralEvent observers
	/// @remarks should only be called when there are no other suitable event types
	///  since specific event notifications will automatically fire a GeneralEvent signal
	void
	notify_general_obs( GeneralEvent const & e );


	/// @brief notify EnergyEvent observers
	/// @param e the event
	/// @param fire_general fire a GeneralEvent afterwards? default true
	void
	notify_energy_obs( EnergyEvent const & e, bool const fire_general = true );


	/// @brief notify ConformationEvent observers
	/// @param e the event
	/// @param fire_general fire a GeneralEvent afterwards? default true
	void
	notify_conformation_obs( ConformationEvent const & e, bool const fire_general = true );


private: // Pose as-an-observer methods


	/// @brief upon receiving a conformation::signals::XYZEvent
	void
	on_conf_xyz_change( core::conformation::signals::XYZEvent const & event );


private:



	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	/////////////////////////data ////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	/// the kinematics object, wraps a fold_tree and an atom_tree
	/**
		 responsible for updating Residue coords in response to changes
		 in the internal coordinates
		 also keeps track of what has changed since last update neighbors call
	**/
	ConformationOP conformation_;


	/// the cached energies object, stores total,rsd,and rsd-pair energies and neighbor info

	/// stores information from the last energy evaluation
	/// handles cacheing rsd and rsd-pair energies
	scoring::EnergiesOP energies_;

	// data which the pose can compute (akin to "decoystats" of Rosetta++)
	mutable metrics::PoseMetricContainerOP metrics_;

	/// @brief BasicDataCache indexed by enum in core/pose/datacache/CacheableDataType.hh
	/// @remarks contains data we can tuck inside the pose for convenience access
	/// @warning DataCache must always be initialized with the number of cacheable
	///  data types -- see the last enum entry.
	BasicDataCacheOP data_cache_;

	/// @brief ObserverCache indexed by enum in core/pose/datacache/CacheableObserverType.hh
	/// @warning ObserverCache must always be initialized with the number of cacheable
	///  observer types -- see the last enum entry.
	ObserverCacheOP observer_cache_;

	/// @brief pdb info
	PDBInfoOP pdb_info_;

	// constraint set
	ConstraintSetOP constraint_set_;

	/// @brief DestructionEvent observers
	/// @remarks notification only occurs when Pose object is destroyed
	mutable utility::signals::BufferedSignalHub< void, DestructionEvent > destruction_obs_hub_;

	/// @brief GeneralEvent observers
	/// @remarks GeneralEvent observers will be notified whenever any signal
	///  derived from GeneralEvent occurs.
	mutable utility::signals::BufferedSignalHub< void, GeneralEvent > general_obs_hub_;

	/// @brief EnergyEvent observers
	mutable utility::signals::BufferedSignalHub< void, EnergyEvent > energy_obs_hub_;

	/// @brief ConformationEvent observers
	/// @remarks fires when Conformation experiences a coordinate change (conformation::signals::XYZEvent)
	mutable utility::signals::BufferedSignalHub< void, ConformationEvent > conformation_obs_hub_;


}; // class Pose

/// @brief Test IO operator for debug and Python bindings
std::ostream & operator << ( std::ostream & os, Pose const & pose);

} // namespace pose
} // namespace core


#endif // INCLUDED_core_pose_Pose_HH
