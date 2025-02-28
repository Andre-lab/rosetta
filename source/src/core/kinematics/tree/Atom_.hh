// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/tree/Atom_.hh
/// @brief  Kinematics Atom abstract base class
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_tree_Atom__hh
#define INCLUDED_core_kinematics_tree_Atom__hh


// Package headers
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID.hh>

// Numeric headers


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace kinematics {
namespace tree {


/// @brief Kinematics Atom abstract base class
class Atom_ : public Atom
{
private: // Types

	typedef  Atom  Super;


public: // Types

	using Super::update_xyz_coords;
	using Super::update_internal_coords;
	using Super::insert_atom;


protected: // Creation

	/// @brief Default constructor
	inline
	Atom_() :
		raw_parent_( nullptr ),
		dof_refold_index_( 0 )
	{
		atoms_.reserve( 4 );
	}

	/// @brief Copy constructor
	/// @note  Copies value and type state but not context state
	inline
	Atom_( Atom_ const & atom ) :
		Super( atom ),
		raw_parent_( nullptr ),
		position_( atom.position_ ),
		dof_refold_index_( atom.dof_refold_index_ )
	{
		atoms_.reserve( 4 );
	}


public: // Creation

	/// @brief Destructor
	// does not clear pointers, use erase
	~Atom_() override
	{}


protected: // Assignment

	/// @brief Copy assignment
	/// @note  Copies value and type state but not context state
	inline
	Atom_ &
	operator =( Atom_ const & atom )
	{
		if ( this != &atom ) {
			position_ = atom.position_;
			dof_refold_index_ = atom.dof_refold_index_;
		}
		return *this;
	}


public: // Methods

	/// @brief copy atom with new memory allocation
	AtomOP
	clone( AtomAP parent_in, AtomPointer2D & atom_pointer ) const override;

	// assumes coords for our input stub are good
	/// @brief update xyz position of this atom and its offspring atoms
	void
	update_xyz_coords() override;

	/// @brief update internal coords of this atom and its offspring atoms (if recursive)
	void
	update_internal_coords(
		bool const recursive
	) override;

	/// @brief Update internal coordinates for this atom and possibly all children.
	/// @details If recursive is true, we update children, grandchildren, great-grandchildren,
	/// etc., but we don't use a recursive algorithm.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	update_internal_coords(
		Stub & stub,
		bool const recursive = true
	) override;

	/// @brief for DOFs which must be kept fixed due to topology of tree
	/** eg, phi of stub_atoms for jump_atoms */
	inline
	bool
	keep_dof_fixed(
		DOF_Type const  //type
	) const override
	{
		return false;
	}

	/// @brief dihedral angle between two bonded children to this atom
	Real
	dihedral_between_bonded_children(
		Atom const & child1,
		Atom const & child2
	) const override;

	/// @brief dump out AtomID for this atom, its parent and all its offspring
	void
	show() const override;

	/// @brief dump out AtomID for this atom, its parent and all its offspring up to n_level
	void
	show(int const & n_level) const override;

	///////////////////////////////////////////////////////////////////////////
	/// @brief update domain map
	void
	update_domain_map(
		int & current_color,
		int & biggest_color,
		DomainMap & domain_map,
		AtomID_Mask const & dof_moved,
		AtomID_Mask const & atom_moved
	) const override;


	///////////////////////////////////////////////////////////////////////////
	// manage atom_list

	/// @brief starting const iterator of the children atom list
	inline
	Atoms_ConstIterator
	atoms_begin() const override
	{
		return atoms_.begin();
	}

	/// @brief ending const iterator of the children atom list
	inline
	Atoms_ConstIterator
	atoms_end() const override
	{
		return atoms_.end();
	}

	/// @brief starting iterator of the children atom list
	inline
	Atoms_Iterator
	atoms_begin() override
	{
		return atoms_.begin();
	}

	/// @brief ending iterator of the children atom list
	inline
	Atoms_Iterator
	atoms_end() override
	{
		return atoms_.end();
	}

	/// @brief number of children atoms
	inline
	Size
	n_atom() const override
	{
		return atoms_.size();
	}

	/// @brief append an atom as this atom's child
	void
	append_atom( AtomOP ) override;

	/// @brief remove an atom from this atom's children
	void
	delete_atom( AtomOP ) override;

	/// @brief insert an atom as this atom's child
	void
	insert_atom( AtomOP ) override;


	/// @brief tries to insert at the position specified by the second argument
	void
	insert_atom( AtomOP, int const /*index*/ ) override;

	/// @brief replace the old atom by the new atom in the child atom list
	void
	replace_atom(
		AtomOP const old_atom,
		AtomOP const new_atom
	) override;

	/// @brief get non-jump atom by its index from the children atoms list
	AtomCOP
	get_nonjump_atom(
		Size const i
	) const override;

	/// @brief number of the child atoms
	Size
	n_children() const override;

	/// @brief number of the non-jump child atoms
	Size
	n_nonjump_children() const override
	{
		return atoms_end() - nonjump_atoms_begin();
	}

	/// @brief get a child atom by index (const method)
	AtomCOP
	child( Size const k ) const override;

	/// @brief get a child atom by index
	AtomOP
	child( Size const k ) override;

	/// @brief the atom-index of this child
	Size
	child_index( AtomCOP child ) const override;

	/// @brief the atom-index of this child
	Size
	raw_child_index( Atom const * child ) const override;

	/// @brief whether atom1 is downstream of this atom.
	bool
	downstream( AtomCOP atom1 ) const override;


public: // Properties

	/// @brief Atom identifier
	inline
	AtomID const &
	id() const override
	{
		return atom_id_;
	}


	/// @brief AtomID assignment
	inline
	void
	id( AtomID const & id_in ) override
	{
		atom_id_ = id_in;
	}


	/// @brief Atom identifier
	inline
	AtomID const &
	atom_id() const override
	{
		return atom_id_;
	}


	/// @brief Position
	inline
	Position const &
	position() const override
	{
		return position_;
	}


	/// @brief Position assignment
	inline
	void
	position( Position const & position_a ) override
	{
		debug_assert( position_a.is_finite() );
		position_ = position_a;
	}


	/// @brief Position
	inline
	Position const &
	xyz() const override
	{
		return position_;
	}


	/// @brief Position assignment
	inline
	void
	xyz( Position const & position_a ) override
	{
		position_ = position_a;
	}


	/// @brief x coordinate
	inline
	Length const &
	x() const override
	{
		return position_.x();
	}


	/// @brief y coordinate
	inline
	Length const &
	y() const override
	{
		return position_.y();
	}


	/// @brief z coordinate
	inline
	Length const &
	z() const override
	{
		return position_.z();
	}


	/// @brief Distance to an Atom
	inline
	Length
	distance( Atom const & atom ) const override
	{
		return position_.distance( atom.position() );
	}


	/// @brief Distance squared to an Atom
	inline
	Length
	distance_squared( Atom const & atom ) const override
	{
		return position_.distance_squared( atom.position() );
	}


	/// @brief Parent atom pointer
	inline
	AtomOP
	parent() override
	{
		return parent_.lock();
	}


	/// @brief Parent atom pointer
	inline
	AtomCOP
	parent() const override
	{
		return parent_.lock();
	}


	void
	parent( AtomAP parent_in ) override
	{
		parent_ = parent_in;
		AtomOP parent_op = parent_.lock();
		raw_parent_ = parent_op.get();
	}


	/// @brief stub centerd at this atom
	Stub
	get_stub() const override;

	/// @brief stub used to build this atom
	Stub
	get_input_stub() const override;

	/// @brief stub atom1 's id
	inline                                 // PHIL: These AtomID fxns could be faster by implementing analogs to the atom lookup calls at the cost of more near-duplicate code
	AtomID const &
	stub_atom1_id() const override
	{
		return stub_atom1()->id();
	}

	/// @brief stub atom2's id
	inline
	AtomID const &
	stub_atom2_id() const override
	{
		return raw_stub_atom2()->id();
	}

	/// @brief stub atom3's id
	inline
	AtomID const &
	stub_atom3_id() const override
	{
		return stub_atom3()->id();
	}

	/// @brief the center of the input stub for refolding this atom
	/** it is its parent*/
	inline
	AtomCOP
	input_stub_atom0() const override
	{
		return parent();
	}

	/// @brief the first atom to construct the input stub for refolding this atom
	/** it is its parent's stub_atom1, which normally the parent itself*/
	inline
	AtomCOP
	input_stub_atom1() const override
	{
		AtomCOP p = parent();
		debug_assert( p != nullptr );
		return p->stub_atom1();
	}

	/// @brief the second atom to construct the input stub for refolding this atom
	/** it is its parent's stub_atom2, which normally the parent's parent*/
	inline
	AtomCOP
	input_stub_atom2() const override
	{
		AtomCOP p = parent();
		debug_assert( p != nullptr );
		return p->stub_atom2();
	}

	/// @brief the third atom to construct the input stub for refolding this atom
	/** it is either its previous sibling or its parent's stub_atom3,*/
	inline
	AtomCOP
	input_stub_atom3() const override
	{
		AtomCOP parent_op = parent();
		AtomCOP sibling_op( previous_sibling() );
		if ( is_jump() || !sibling_op || sibling_op->is_jump() ||
				is_collinear( *(parent_op->stub_atom1()), *(parent_op->stub_atom2()), *sibling_op) ||
				( parent_op->is_jump() && sibling_op->id() == parent_op->stub_atom2_id() ) ) {
			return parent_op->stub_atom3();
		} else {
			return sibling_op;
		}
	}

	/// @brief input stub atom0's id
	inline
	AtomID const &
	input_stub_atom0_id() const override
	{
		return raw_input_stub_atom0()->id();
	}

	/// @brief input stub atom1's id
	inline
	AtomID const &
	input_stub_atom1_id() const override
	{
		return raw_input_stub_atom1()->id();
	}

	/// @brief input stub atom2's id
	inline
	AtomID const &
	input_stub_atom2_id() const override
	{
		return raw_input_stub_atom2()->id();
	}

	/// @brief input stub atom3's id
	inline
	AtomID const &
	input_stub_atom3_id() const override
	{
		return input_stub_atom3()->id();
	}


	/// @brief  routines for navigating the tree
	/// find the sibling atom before itself
	AtomCOP
	previous_sibling() const override;

	/// @brief find the child atom before this child in the list
	AtomCOP
	previous_child(
		AtomCOP child
	) const override;

	/// @brief find the child atom after this child in the list
	AtomOP
	next_child(
		AtomCOP child
	) override;

	/// @brief whether a Stub can be defined for this atom
	bool
	stub_defined() const override;


protected: // Methods

	/// @brief when subtrees have changed their coordinates
	void
	update_child_torsions(
		AtomOP const child
	) override;

	/// @brief constant iterator of the first non-jump (bonded) atom in the vector of children atoms.
	Atoms_ConstIterator
	nonjump_atoms_begin() const override;

	/// @brief iterator of the first non-jump (bonded) atom in the vector of children atoms.
	Atoms_Iterator
	nonjump_atoms_begin() override;

	/// @brief helper function to abort if something is wrong in atom tree
	void
	abort_bad_call() const;


	/// @brief Transform atom and children by linear transformation
	void
	transform_Ax_plus_b_recursive(
		Matrix const & A,
		Vector const & b,
		ResidueCoordinateChangeList & res_change_list
	) override;


	void
	get_path_from_root( utility::vector1< AtomCAP > & path ) const override;


	bool
	atom_is_on_path_from_root( AtomCOP atm ) const override;


	/// @brief Records this atom as having a changed DOF in the input list
	/// of Atoms with changed DOFs.  For use in output-sensitive refold subroutine.
	void
	note_dof_change(
		AtomDOFChangeSet & changset
	);

	/// @brief To ensure proper function of the output-senstive refold
	/// subroutine, derived classes must invoke this function during their
	/// update_xyz_coord subroutines.
	inline
	void
	note_xyz_uptodate()
	{
		dof_refold_index_ = 0;
	}


public:

	/// @brief base class implementation that traverses the subtree routed at this node
	/// in the depth-first traversal of the atoms requiring coordinate updates.
	void
	dfs(
		AtomDOFChangeSet & changeset,
		ResidueCoordinateChangeList & res_change_list,
		Size const start_atom_index
	) const override;


protected:

	/// @brief read access for derived classes
	inline
	Size
	dof_refold_index() const {
		return dof_refold_index_;
	}


protected: // Fields -- should be private...
	// private: // Fields

	/// @brief Atom ID
	AtomID atom_id_;

	/// @brief Associated conformation Atom
	//conformation::AtomAP conformation_atom_p_;

	/// @brief Parent atom pointer
	AtomAP parent_;

	/// @brief Workaround copy of the parent pointer for use in functions
	/// where locking and unlocking the parent_ pointer would be
	/// prohibitively slow
	Atom * raw_parent_;

	/// @brief xyz
	PointPosition position_;

	/// @brief Children atom pointers
	Atoms atoms_;

public:

	Atom const *
	raw_parent() const override;

	Atom const *
	raw_previous_sibling() const override;

	Atom const *
	raw_previous_child(
		Atom const * child
	) const override;

	Atom const *
	raw_input_stub_atom0() const override;

	Atom const *
	raw_input_stub_atom1() const override;

	Atom const *
	raw_input_stub_atom2() const override;

	Atom const *
	raw_input_stub_atom3() const override;

	Atom const *
	raw_get_nonjump_atom(
		Size const i
	) const override;


private:

	/// @brief Track my position in my owner's list of Atoms with modified DOFs.
	/// 0 when my dofs have not changed since the last update_coords.
	Size dof_refold_index_;

#ifdef    SERIALIZATION
	friend class cereal::access;
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // Atom_

} // namespace tree
} // namespace kinematics
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_kinematics_tree_Atom_ )
#endif // SERIALIZATION


#endif // INCLUDED_core_kinematics_Atom__HH
