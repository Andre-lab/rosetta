// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/AtomTree.cc
/// @brief  Atom tree class
/// @author Phil Bradley


// Unit headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/ResidueCoordinateChangeList.hh>

// Package headers
// AUTO-REMOVED #include <core/kinematics/DomainMap.hh>
// AUTO-REMOVED #include <core/kinematics/tree/BondedAtom.hh>
// AUTO-REMOVED #include <core/kinematics/tree/JumpAtom.hh>
// AUTO-REMOVED #include <core/kinematics/util.hh>

#include <basic/basic.hh> // periodic_range
#include <basic/prof.hh> // profiling
#include <basic/Tracer.hh> // profiling


// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray3D.hh>
//#include <ObjexxFCL/FArray4D.h>
//#include <ObjexxFCL/formatted.io.h>

// Numeric headers
// AUTO-REMOVED #include <numeric/all.fwd.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/assert.hh>
// AUTO-REMOVED #include <utility/io/orstream.hh>

// C++ headers
#include <cstdlib>
// AUTO-REMOVED #include <cstdio>

//Auto Headers
#include <core/kinematics/AtomWithDOFChange.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.hh>



namespace core {
namespace kinematics {

static basic::Tracer TR( "core.kinematics.AtomTree" );

/////////////////////////////////////////////////////////////////////////////
/// @details this will claim the tree as our own. new_root has information about its children,
/// and they have information about their children. From those atom positions, internal
/// coordinates can be updated and atom pointers will be added into the AtomTree map.
/// @note that we steal the atoms, ie it's incorporated into the AtomTree (by recording
/// their pointers),not cloned
AtomTree::AtomTree(
	AtomPointer2D const & new_atom_pointer,
	bool const from_xyz // = true
):
	root_( 0 ),
	atom_pointer_(), // default_setting_ = null pointer
	internal_coords_need_updating_( false ),
	xyz_coords_need_updating_( false ),
	topological_match_to_( 0 ),
	external_coordinate_residues_changed_( new ResidueCoordinateChangeList )
{
	replace_tree( new_atom_pointer, from_xyz );
	external_coordinate_residues_changed_->total_residue( new_atom_pointer.size() );
}

AtomTree::AtomTree():
	root_(0),
	atom_pointer_(),
	internal_coords_need_updating_( false ),
	xyz_coords_need_updating_( false ),
	topological_match_to_( 0 ),
	external_coordinate_residues_changed_( new ResidueCoordinateChangeList )
{}

/// @brief Destructor
AtomTree::~AtomTree()
{
	clear();
}



/////////////////////////////////////////////////////////////////////////////
/// @details copy ctor, uses operator=
AtomTree::AtomTree( AtomTree const & src ) :
	utility::pointer::ReferenceCount(),
	root_( 0 ), /// without this initialization, the destruction of this
	/// uninitialized pointer might have disasterous consequences
	atom_pointer_(), // default_setting_ = null pointer
	internal_coords_need_updating_( false ),
	xyz_coords_need_updating_( false ),
	topological_match_to_( 0 ),
	external_coordinate_residues_changed_( new ResidueCoordinateChangeList )
{
	*this = src;
}



/**
/////////////////////////////////////////////////////////////////////////////
///
///@li update xyz or internal coord before adding this atom.
///@li notify xyz or internal coord need to be updated after adding this atom.
///@li if atom id2 is not in the AtomID_Map and root atom is not set, atom id1
///    is set as the root of the tree (the parent of the root is 0).
///@li if atom id2 is in the AtomID_Map, atom id1 is added into the map either
///    as a bonded_atom or jump_atom and atom id2 is set as its parent.
///@li Throw an error if atom id2 is not in the map and root atom is already set.
///
void
AtomTree::add_atom(
	AtomID const & id1,
	AtomID const & id2,
	bool const add_bonded_atom,
	bool const from_xyz
)
{
	if ( from_xyz ) {
		update_xyz_coords();
	} else {
		update_internal_coords();
	}


	Atom* parent(0);
	if ( !id2.valid() || // root
		!atom_pointer_.has( id2 ) ||
		atom_pointer_[ id2 ] == 0 ) {
		//
		if ( root_ != 0 ) {
			utility_exit_with_message("add_atom: parent not in the tree");
		}
	} else {
		parent = atom_pointer_[ id2 ];
	}

	// create the new atom
	Atom* atom( add_bonded_atom ?
		static_cast< Atom* >( new BondedAtom() ) :
		static_cast< Atom* >( new JumpAtom() ) );
	atom->id( id1 );
	atom_pointer_.set( id1, atom );

	if ( parent ) {
		parent->append_atom( atom );
	} else {
		root_ = atom;
		atom->parent( 0 );
	}


	if ( from_xyz ) {
		internal_coords_need_updating_ = true;
	} else {
		xyz_coords_need_updating_ = true;
	}

}
**/


/////////////////////////////////////////////////////////////////////////////
void
AtomTree::find_root_from_atom_pointer()
{
	root_ = 0;
	for ( Size i=1; i<= atom_pointer_.size(); ++i ) {
		for ( Size j=1; j<= atom_pointer_[i].size(); ++j ) {
			assert( atom_pointer_[i][j] && atom_pointer_[i][j]->id() == AtomID( j,i ) );
			if ( atom_pointer_[i][j]->parent() == 0 ) {
				assert( !root_ );
				root_ = atom_pointer_[i][j]();
			}
		}
	}
}


/////////////////////////////////////////////////////////////////////////////
///
/// @details fill the AtomTree with a new tree of atoms by recording their pointers in
/// the map. Sync internal and xyz coords.
void
AtomTree::replace_tree(
	AtomPointer2D const & new_atom_pointer,
	bool const from_xyz // = true
)
{
	clear();

	atom_pointer_ = new_atom_pointer;

	find_root_from_atom_pointer();

	external_coordinate_residues_changed_->total_residue( atom_pointer_.size() );

	if ( from_xyz ) {
		internal_coords_need_updating_ = true;
		xyz_coords_need_updating_ = false;

		update_internal_coords();
	} else {
		xyz_coords_need_updating_ = true;
		internal_coords_need_updating_ = false;

		update_xyz_coords();
	}


	// anything that depends on the tree topology needs to be updated
	set_new_topology();
}

/**
void
AtomTree::setup_backrub_segment(
	utility::vector1< AtomID > const & mainchain,
	AtomID const & downstream_id, // mainchain child of last atom in mainchain vector
	utility::vector1< std::pair< Size, Size > > const & edges,
	Size const first_new_pseudo_residue
)
{

	// this assumes that the segment of mainchain has exactly one connection in and exactly one connection out, (to
	// the downstream_id atom) and all of the atoms to be replaced are bonded atoms
	//
	for ( Size i=1; i<= mainchain.size(); ++i ) {
		assert( !( atom_pointer_[ mainchain[i] ]->is_jump() ) );
	}
	AtomID const last_mainchain_id( mainchain[ mainchain.size() ] );
	assert( atom_pointer_[ downstream_id ]->parent()->id() == last_mainchain_id );

	// operation is done "by xyz" leaving the internal coords out of date
	update_xyz_coords();

	Atom * subtree_root
		( setup_backrub_atom_tree( mainchain, downstream_id, atom_pointer_, edges, first_new_pseudo_residue ) );
	assert( subtree_root->id() == mainchain[1] );

	Atom * old_root( atom_pointer_[ mainchain[1] ] );
	Atom * anchor( old_root->parent() );

	Atom * downstream_atom( atom_pointer_[ downstream_id ] );
	downstream_atom->parent()->delete_atom( downstream_atom ); // delete the old outgoing connection

	// update atom_pointer_
	subtree_root->update_atom_pointer( atom_pointer_, true );

	// insert the new tree
	anchor->replace_atom( old_root, subtree_root ); // incoming connection
	atom_pointer_[ last_mainchain_id ]->insert_atom( downstream_atom ); // outoing connection

	// erase the old tree
	// this call will erase and delete all of old_root's children
	old_root->erase();
	// free the last bit of old data
	delete old_root;

	internal_coords_need_updating_ = true;

}
**/


/// @details  This is a helper function to find a linear transform that when applied to the downstream stub
/// has the effect that RT( instub, transformed-downstream-stub) == target_rt
void
find_stub_transform(
	Stub const & stub1, // upstream stub
	Stub const & stub2, // downstream stub
	RT const & rt, // the target RT
	Stub::Matrix & A,
	Vector & b
)
{
	Stub::Matrix const & M1( stub1.M ), M2( stub2.M ), R( rt.get_rotation() );
	Vector const & v1( stub1.v ), v2( stub2.v ), t( rt.get_translation() );

	// look for a transformation of the form x |----> A*x + b
	//
	// this will change stub2 to stub2' with M2' = A * M2, v2' = A*v2 + b
	//
	// if we let (R,t) be the target RT, then we want
	//
	//  R = M1^T * M2' = M1^T * A * M2  ==> A = M1 * R * M2^T
	//
	//  t = M1^T * ( v2' - v1 ) ==> v2' = M1 * t + v1, which with b = v2' - A*v2 gives b = M1 * t + v1 - A * v2
	//

	A = M1 * R * M2.transposed();
	b = M1 * t + v1 - A * v2;
}





/////////////////////////////////////////////////////////////////////////////
/// assumes one incoming and at most one outgoing
/// Need fancier fxn for general case
///
void
AtomTree::delete_seqpos( Size const seqpos )
{
	Size const old_size( size() ), new_size( old_size - 1 );

	// find the anchor, root, and perhaps child atoms
	Size const natoms( atom_pointer_[seqpos].size() );
	Atom* anchor(0), *root(0), *child(0);
	for ( Size i=1; i<= natoms; ++i ) {
		Atom* atom( atom_pointer_[seqpos][i]() );
		if ( !atom ) continue;
		if ( Size(atom->parent()->id().rsd()) != seqpos ) {
			assert( !anchor );
			root = atom;
			anchor = atom->parent();
			// could break but debug 1 incoming connxn by continuing
		}
		for ( Atom::Atoms_ConstIterator iter=atom->atoms_begin(), iter_end = atom->atoms_end(); iter!= iter_end; ++iter ) {
			if ( Size((*iter)->id().rsd()) != seqpos ) {
				assert( !child );
				child = *iter;
				// could break but debug at most 1 outgoing connxn by continuing
			}
		}
	}

	if ( !anchor ) {
		utility_exit_with_message( "AtomTree::delete_seqpos can't handle deleting the root residue");
	}

	// rewire connxns
	if ( child ) {
		anchor->replace_atom( root, child );
	} else {
		anchor->delete_atom( root );
	}

// 	// now delete atoms
// 	for ( Size i=1; i<= natoms; ++i ) {
// 		delete atom_pointer_[seqpos][i];
// 	}

	atom_pointer_[ seqpos ].resize(0);

	// now renumber
	utility::vector1< int > old2new( old_size );
	for ( Size i=1; i<= old_size; ++i ) {
		if      ( i <  seqpos ) old2new[i] = i;
		else if ( i == seqpos ) old2new[i] = 0;
		else                    old2new[i] = i-1;
	}

	update_sequence_numbering( new_size, old2new );

	// anything that depends on the tree topology needs to be updated
	set_new_topology();
}


/// @note  This can also handle the case of inserting a residue, proved that we've already renumbered atomtree leaving an empty slot at seqpos.
/// @note  If the new residue contains the root of the tree, incoming.atom1 should be BOGUS_ATOM_ID
/// @note  Note that for all BondID's atom1 should be the parent and atom2 should be the child

/// Couple of special cases:
/// -- appending a residue
/// -- inserting a residue
/// -- replacing the root residue
/// -- inserting a new root residue

void
AtomTree::replace_residue_subtree(
	id::BondID const & incoming,
	utility::vector1< id::BondID > const & outgoing,
	AtomPointer1D const & new_atoms // this will be the new slice of our atom_pointer
)
{
	// general rule:
	//
	// before we set or get a type of coordinate, have to
	// do an update
	//
	// in this case, we are "setting" the xyz's for the new subtree
	// this will put the internal coords out of date, so we want to
	// ensure that the xyz's are up to date at the start, otherwise
	// we can end up in the disastrous situation of having both
	// internal and xyz coords out of data -- and thus no "reference"
	// state from which to update!
	//
	update_xyz_coords();

	Size const seqpos( incoming.atom2.rsd() ); // should be Size but AtomID::rsd() returns int
	AtomPointer1D const & old_atoms( atom_pointer_[ seqpos ] );

	// confirm that bond id's go from parent to child, atom1 to atom2
	for ( Size i=1; i<= outgoing.size(); ++i ) assert( outgoing[i].atom1.rsd() == seqpos );
	for ( Size i=1; i<= new_atoms.size(); ++i ) assert( new_atoms[i]->id() == AtomID( i, seqpos ) );

	//
	Atom * anchor_atom(0);
	Atom * old_root_atom(0);
	Atom * new_root_atom( new_atoms[ incoming.atom2.atomno() ]() );

	if ( incoming.atom1.valid() ) anchor_atom = atom_pointer( incoming.atom1 );

	// rewire the outgoing connections -- these are added to new atoms in the order that they come in the outgoing list
	// the children in these connections will all have parents in seqpos unless we are inserting into an empty slot,
	// ie unless old_atoms.empty()
	for ( Size i=1; i<= outgoing.size(); ++i ) {
		Atom * child( atom_pointer( outgoing[i].atom2 ) );
		Atom * old_parent( child->parent() );
		assert( child->id().rsd() != seqpos );
		if ( !old_parent ) {
			// we're becoming the new root residue
			assert( old_atoms.empty() && !incoming.atom1.valid() ); // implies anchor_atom == 0
			assert( !old_root_atom );
			old_root_atom = child;
		} else if ( old_parent->id().rsd() != seqpos ) {
			// we're inserting into a bond
			assert( old_atoms.empty() && old_parent->id() == incoming.atom1 );
			assert( !old_root_atom );
			old_root_atom = child;
		} else {
			// only necessary for debugging purposes
			old_parent->delete_atom( child );
		}
		Atom * new_parent( new_atoms[ outgoing[i].atom1.atomno() ]() );
		new_parent->insert_atom( child );
	}

	// potentially have to look for old_root_atom
	if ( old_root_atom ) {
		assert( old_atoms.empty() );
	} else {
		for ( Size i=1; i<= old_atoms.size(); ++i ) {
			Atom * old_atom( old_atoms[i]() );
			assert( old_atom ); // atom_pointer_ is ragged, always keep dimension equal to actual number of atoms
			if ( ! old_atom->parent() ) {
				// this was the root of the atomtree
				assert( !incoming.atom1.valid() );
				assert( !old_root_atom );
				old_root_atom = old_atom;
			} else if ( old_atom->parent()->id().rsd() != seqpos ) {
				// this is the root of the old tree
				assert( incoming.atom1 == old_atom->parent()->id() );
				assert( !old_root_atom );
				old_root_atom = old_atom;
			}
			// this is just debugging to confirm that outgoing vector actually contains all the outgoing connections
			for ( Size i=0; i< old_atom->n_children(); ++i ) assert( old_atom->child(i)->id().rsd() == seqpos );
		}
	}

	// rewire the incoming connection
	if ( anchor_atom ) {
		if ( old_root_atom ) {
			atom_pointer( incoming.atom1 )->replace_atom( old_root_atom, new_root_atom );
		} else {
			assert( outgoing.empty() && old_atoms.empty() );
			atom_pointer( incoming.atom1 )->insert_atom( new_root_atom );
		}

	} else {
		assert( root_ == old_root_atom );
		root_ = new_root_atom;
		new_root_atom->parent(0);
	}


	// now nuke the old atoms
	atom_pointer_[ seqpos ].clear();
	atom_pointer_[ seqpos ] = new_atoms;

	//
	// we've added the new atoms assuming that their xyz coords are
	// valid, but their internal coords (especially at the junctions
	// of the old and new) are likely to be messed up. This is also
	// true, eg, for any younger siblings of subtree_root
	//
	internal_coords_need_updating_ = true;

	// anything that depends on the tree topology needs to be updated
	set_new_topology();

}



/////////////////////////////////////////////////////////////////////////////
/// @brief retrieve a specific DOF by its ID.
Real
AtomTree::dof( DOF_ID const & id ) const
{
	update_internal_coords();
	return atom_pointer_[ id.atom_id() ]->dof( id.type() );
}

/////////////////////////////////////////////////////////////////////////////
/// @brief retrieve the xyz position of an atom in the tree by its AtomID
PointPosition const &
AtomTree::xyz( AtomID const & id ) const
{
	update_xyz_coords();
	// if ( !has(id) ) {
// 		std::cerr << "AtomTree::atom_pointer_ has not the atom " << id << std::endl;
// 	}
// 	runtime_assert( has( id ) );
	return atom_pointer_[ id ]->position();
}

/////////////////////////////////////////////////////////////////////////////
/// @brief retreive a kinematic Atom in the tree by its AtomID
tree::Atom const &
AtomTree::atom( AtomID const & id ) const
{
	// we don't know what kind of access the caller may perform:
	update_internal_coords();
	update_xyz_coords();

	return *(atom_pointer_[ id ]);
}

/////////////////////////////////////////////////////////////////////////////
/// @brief retreive a kinematic Atom in the tree by its AtomID -- no update!
tree::Atom const &
AtomTree::atom_dont_do_update( AtomID const & id ) const
{
	return *(atom_pointer_[ id ]);
}

/////////////////////////////////////////////////////////////////////////////
/// @brief retrieve a Jump in the tree by that JumpAtom's AtomID.
/// @note will abort if a BondedAtom's AtomID is passed in.
Jump const &
AtomTree::jump( AtomID const & id ) const
{
	update_internal_coords();
	return atom_pointer_[ id ]->jump();
}


/////////////////////////////////////////////////////////////////////////////
/// @brief get the DOF_ID of a torsion angle given those four atoms which define this torsion
///
/// @details an "offset" value is also calculated that torsion(id1,id2,id3,id4) = dof( dof_id ) + offset.
/// A BOGUS_DOF_ID will be returned if no proper DOF can be found for these four atoms.
/// offset is mainly for an atom with a previous sibling as the torsion(PHI) attached
/// to it is calculated as improper angle with respect to its sibling.
id::DOF_ID
AtomTree::torsion_angle_dof_id(
	AtomID const & atom1_in_id,
	AtomID const & atom2_in_id,
	AtomID const & atom3_in_id,
	AtomID const & atom4_in_id,
	Real & offset
) const
{
	using numeric::conversions::degrees;
	using numeric::constants::d::pi_2;
	using numeric::constants::d::pi;
	using numeric::dihedral;
	using numeric::dihedral_radians;

	bool const debug( false );

	// we use the internal dof's if necessary to calculate the offset
	update_internal_coords();

	if ( debug ) update_xyz_coords();

	assert( atom_pointer( atom1_in_id ) && atom_pointer( atom2_in_id ) &&
					atom_pointer( atom3_in_id ) && atom_pointer( atom4_in_id ) );

	// STUART -- (low priority) I'd like to be able to cache
	// the results of this calculation to allow faster access.


	Atom const
		* const atom1_in( atom_pointer( atom1_in_id ) ),
		* const atom2_in( atom_pointer( atom2_in_id ) ),
		* const atom3_in( atom_pointer( atom3_in_id ) ),
		* const atom4_in( atom_pointer( atom4_in_id ) );
	Atom const * atom1( atom1_in ), *atom2( atom2_in ),
		*atom3( atom3_in ), *atom4( atom4_in );
	// reorder the atoms if necessary
	// we want it to be the case that atom4 has
	// input_stub_atom1 == atom3 and
	// input_stub_atom2 == atom2
	//
	if ( !atom4->is_jump() &&
		!atom4->keep_dof_fixed( PHI ) && // not stub atom3 of a jump
		atom4->input_stub_atom1() == atom3 &&
		atom4->input_stub_atom2() == atom2 &&
		atom4->input_stub_atom3() == atom1 ) {
		// pass -- this is what we want, perfect match!!
	} else if ( !atom1->is_jump() &&
		!atom1->keep_dof_fixed( PHI ) && // not stub atom3 of a jump
		atom1->input_stub_atom1() == atom2 &&
		atom1->input_stub_atom2() == atom3 &&
		atom1->input_stub_atom3() == atom4 ) {
		// perfect case in reverse, want to reverse the order of the atoms
		atom1 = atom4_in;
		atom2 = atom3_in;
		atom3 = atom2_in;
		atom4 = atom1_in;
	} else if ( !atom4->is_jump() &&
		!atom4->keep_dof_fixed( PHI ) && // not stub atom3 of a jump
		atom4->input_stub_atom1() == atom3 &&
		atom4->input_stub_atom2() == atom2 ) {
		// pass -- this is what we want, not quite as perfect as 1st case though
	} else if ( !atom1->is_jump() &&
		!atom1->keep_dof_fixed( PHI ) && // not stub atom3 of a jump
		atom1->input_stub_atom1() == atom2 &&
		atom1->input_stub_atom2() == atom3 ) {
		// reverse the order of the atoms
		atom1 = atom4_in;
		atom2 = atom3_in;
		atom3 = atom2_in;
		atom4 = atom1_in;
	} else {
		// no good!
		return id::BOGUS_DOF_ID;
	}

	offset = 0.0; // initialize

	if ( atom4->input_stub_atom3() == atom1 ) {
		// perfect match
		return DOF_ID( atom4->id(), PHI );
	}

	assert( !atom4->is_jump() &&
		atom4->input_stub_atom0() == atom3 &&
		atom4->input_stub_atom1() == atom3 &&
		atom4->input_stub_atom2() == atom2 &&
		atom4->parent() == atom3 );

	// special case if atoms 1 and 4 are siblings -- not really well defined!
	if ( atom1->parent() == atom3 ) {
		Size const atom1_index( atom3->child_index( atom1 ) ),
			atom4_index( atom3->child_index( atom4 ) );
		if ( atom1_index < atom4_index ) {
			Real const current_value( atom3->dihedral_between_bonded_children( atom1, atom4 ) );
			offset = current_value - atom4->dof(PHI); // since torsion(id1...id4) = dof + offset;

			ASSERT_ONLY( Real const actual_current_value
				( dihedral_radians( atom4->xyz(), atom3->xyz(), atom2->xyz(), atom1->xyz() ) );)
			assert( std::abs( basic::subtract_radian_angles( actual_current_value, current_value ) ) < 1e-3 );
			assert( std::abs( basic::subtract_radian_angles( current_value, atom4->dof( PHI ) + offset ) ) < 1e-3 );

			return DOF_ID( atom4->id(), PHI );
		} else {
			Real const current_value( atom3->dihedral_between_bonded_children( atom1, atom4 ) );
			offset = current_value - atom1->dof( PHI ); // since torsion(id1...id4) = dof + offset;
			ASSERT_ONLY( Real const actual_current_value
				( dihedral_radians( atom4->xyz(), atom3->xyz(), atom2->xyz(), atom1->xyz() ) );)
			assert( std::abs( basic::subtract_radian_angles( actual_current_value, current_value ) ) < 1e-3 );
			assert( std::abs( basic::subtract_radian_angles( current_value, atom1->dof( PHI ) + offset ) ) < 1e-3 );
			return DOF_ID( atom1->id(), PHI );
		}
	}
	// atom4 is not the first sibling of atom3, get offset for that.
	if ( atom4 != atom3->get_nonjump_atom(0) ) {
		Atom const * new_atom4( atom3->get_nonjump_atom(0) );
		offset += atom3->dihedral_between_bonded_children( new_atom4, atom4 );

		if ( debug ) { // debugging
			ASSERT_ONLY( Real const actual_dihedral
				( dihedral_radians( atom4->xyz(), atom3->xyz(), atom2->xyz(),
					new_atom4->xyz() ) );)
			assert( std::abs( basic::subtract_radian_angles
					( actual_dihedral, offset ) ) < 1e-3 );
		}

		atom4 = new_atom4;
	}

	DOF_ID dof_id( atom4->id(), PHI );

	Atom const * dof_atom1( atom4->input_stub_atom3() );

	if ( dof_atom1 == atom1 ) {
		return dof_id;
	} else if ( atom1->parent() == atom2 && !atom1->is_jump() &&
		atom3->parent() == atom2 && !atom3->is_jump() ) {

		// handle offset between atom1 and dof_atom1
		// the only case we can do is if atom1->parent() == atom2
		//
		// this will happen, eg for chi1 if we are folding c2n:
		//
		// atom1 =  N    while    dof_atom1 = parent(atom2) =  C
		// atom2 = CA
		// atom3 = CB
		// atom4 = CG
		//

		// a little tricky: we don't want to access the positions,
		// since we can't be sure that they are up to date in a routine
		// that would be called inside dof setting and getting routines
		//
		//
		// need to use some spherical geometry
		//
		// we've got three unit vectors:
		//
		// u1 -- from atom2 to atom1
		// u2 -- from atom2 to dof_atom1
		// u3 -- from atom2 to atom3
		//
		// the angle between u2 and u1 = theta3 = pi - atom1(THETA)
		// the angle between u2 and u3 = theta1 = pi - atom3(THETA)
		// the torsion between u1 and u3 along u2 = phi2 = atom2->bonded_dih(1,3)
		//
		Real
			theta1( pi - atom3->dof( THETA ) ),
			theta3( pi - atom1->dof( THETA ) ),
			phi2( atom2->dihedral_between_bonded_children( atom1, atom3 ) ),
			sign_factor( 1.0 );
		phi2 = basic::periodic_range( phi2, pi_2 );
		if ( phi2 < 0 ) {
			sign_factor = -1;
			phi2 *= -1;
		}
		//If bond angle is varied, these angles can sometimes go out of range,
		// during exploration by minimizer:
		theta3 = basic::periodic_range( theta3, pi_2 );
		if ( theta3 < 0 ) {
			sign_factor *= -1;
			theta3 *= -1;
		}
		theta1 = basic::periodic_range( theta1, pi_2 );
		if ( theta1 < 0 ) {
			sign_factor *= -1;
			theta1 *= -1;
		}
		if ( !( 0 <= theta1 && theta1 <= pi &&
				0 <= theta3 && theta3 <= pi &&
				0 <=   phi2 &&   phi2 <= pi ) ) {
			std::cerr << "dof_atom1 " << dof_atom1->id() << std::endl;
			std::cerr << "atom1 " << atom1->id() << std::endl;
			std::cerr << "atom2 " << atom2->id() << std::endl;
			std::cerr << "atom3 " << atom3->id() << std::endl;
			std::cerr << "atom4 " << atom4->id() << std::endl;
			std::cerr << "THETA1 " << theta1 << std::endl;
			std::cerr << "THETA3 " << theta3 << std::endl;
			std::cerr << "PHI2 " << phi2 << std::endl;

			utility_exit_with_message
				( "AtomTree::torsion_angle_dof_id: angle range error" );

		}
		// want to solve for phi3 -- dihedral between u2 and u1 along axis
		// defined by u3
		//
		// spherical law of cosines says:
		//
		// cos(theta2) = cos(theta1) * cos(theta3) +
		//               sin(theta1) * sin(theta3) * cos(phi2)
		//
		// so we can solve for cos(theta2); then use formula again to solve for
		// cos(phi3).
		//
		Real const cos_theta2 = std::cos( theta1 ) * std::cos( theta3 ) +
			std::sin( theta1 ) * std::sin( theta3 ) * std::cos( phi2 );
		Real const theta2 = std::acos( cos_theta2 );
		assert( 0 <= theta2 && theta2 <= pi );

		// use formula again interchanging 2 and 3
		Real cos_phi3 = ( cos( theta3 ) - cos( theta1 ) * cos_theta2 ) /
			( sin( theta1 ) * sin( theta2 ) );

		// PHIL! really should handle degenerate cases more carefully
		// PHIL! also could reuse some of the trig calculations

		// dof-torsion: from dof_atom1 to atom4  ~  atom4 - dof_atom1
		// desired:     from atom1 to atom4  ~  atom4 - atom1
		//
		// atom4 - atom1 = atom4 - dof_atom1 + ( dof_atom1 - atom1 )
		//
		// so offset == dof_atom1 - atom1 ~ dihedral from atom1 to dof_atom1
		//

		//
		Real const dihedral_from_atom1_to_dof_atom1
			( sign_factor * std::acos( cos_phi3 ) );

		offset += dihedral_from_atom1_to_dof_atom1;

		if ( debug ) { // debugging
			using basic::periodic_range;
			using basic::subtract_radian_angles;
			Real const actual_dihedral_from_atom1_to_dof_atom1
				( dihedral_radians( atom1->xyz(), atom2->xyz(), atom3->xyz(),
					dof_atom1->xyz()) );

			Real const actual_dihedral_of_interest
				( dihedral_radians( atom1_in->xyz(), atom2_in->xyz(),
					atom3_in->xyz(), atom4_in->xyz()) );

			ASSERT_ONLY(Real const dev1
				( std::abs( subtract_radian_angles
					( actual_dihedral_from_atom1_to_dof_atom1,
						dihedral_from_atom1_to_dof_atom1 ) ) );)

			ASSERT_ONLY(Real const dev2
				( std::abs( subtract_radian_angles( actual_dihedral_of_interest,
						dof( dof_id ) + offset ) ) );)


			ASSERT_ONLY(Real const dev3
				( std::abs( subtract_radian_angles
					( dof( dof_id ),
						dihedral_radians( atom4->xyz(),
							atom4->input_stub_atom1()->xyz(),
							atom4->input_stub_atom2()->xyz(),
							atom4->input_stub_atom3()->xyz()))));)

			assert( dev1 < 1e-3 && dev2 < 1e-3 && dev3 < 1e-3 );

			if ( false ) {
				TR.Trace << "offset: " << offset << " sgn_fac: " << sign_factor <<
					' ' << dihedral_from_atom1_to_dof_atom1 << " =?= " <<
					actual_dihedral_from_atom1_to_dof_atom1 <<
					' ' << periodic_range( dof( dof_id ) + offset, pi_2 ) <<
					" =?= " << actual_dihedral_of_interest << std::endl;
			}
		}

		return dof_id;
	}

	// failure
	return id::BOGUS_DOF_ID;
}



/////////////////////////////////////////////////////////////////////////////
//
/// @details brief set a specific DOF in the tree
void
AtomTree::set_dof(
	DOF_ID const & id,
	Real const setting
)
{
	update_internal_coords();
	atom_pointer_[ id.atom_id() ]->set_dof( id.type(), setting, dof_changeset_ );
	xyz_coords_need_updating_ = true;
}

/////////////////////////////////////////////////////////////////////////////

void
AtomTree::set_xyz(
	AtomID const & id,
	PointPosition const & xyz
)
{
	update_xyz_coords();
	atom_pointer_[ id ]->position( xyz );
	internal_coords_need_updating_ = true;
}

/////////////////////////////////////////////////////////////////////////////

void
AtomTree::batch_set_xyz(
	utility::vector1<AtomID> const & ids,
	utility::vector1<PointPosition> const & xyzs
)
{
	runtime_assert( ids.size() == xyzs.size() );
	update_xyz_coords();
	for (core::Size i=1; i<=ids.size(); ++i)
		atom_pointer_[ ids[i] ]->position( xyzs[i] );
	internal_coords_need_updating_ = true;
}


/////////////////////////////////////////////////////////////////////////////
/// @note AtomID must point to a JumpAtom, otherwise it will die.
void
AtomTree::set_jump(
	AtomID const & id,
	Jump const & jump
)
{
	update_internal_coords();
	atom_pointer_[ id ]->jump( jump, dof_changeset_ );
	xyz_coords_need_updating_ = true;
}


/////////////////////////////////////////////////////////////////////////////
/// @note DEPRICATED.  AtomID must point to a JumpAtom, otherwise it will die.
/// Hopefully, the new output-sentitive refold will make this method unneccessary
void
AtomTree::set_jump_now(
	AtomID const & id,
	Jump const & jump
)
{
	/// set_jump_now no longer necessary.  Commented-out implementation below is buggy
	/// because the Conformation does not get its list of changed residues updated; so while
	/// the coordinates in the atom tree are updated, the coordinates in the Conformation are not.
	/// Instead of fixing this bug to handle this special case "set now", rely on the functional
	/// and more efficient output-sensitive "set_jump"
	set_jump( id, jump );

	// We use our parent atoms' positions to find our stub,
	// which we can't do if their positions aren't correct.
	//if( xyz_coords_need_updating_ ) {
	//	set_jump(id, jump);
	//	return;
	//}
	//update_internal_coords();
	//atom_pointer_[ id ]->jump( jump );
	//Stub stub( atom_pointer_[ id ]->get_input_stub() );
	//atom_pointer_[ id ]->update_xyz_coords( stub );
	// No change to xyz_coords_need_updating_ -- we've done our part,
	// but previous changes may have invalidated other coords.
	//xyz_coords_need_updating_ = true;
}


/////////////////////////////////////////////////////////////////////////////
/// @ details it is possible that no DOF can be obtained given these four atoms or the torsion
/// does not match exactly the DOF(PHI) angle. In the former case, a BOGUS_DOF_ID is
/// returned as indicator and in the latter case, an offset value is deducted from
/// input setting to set the real DOF value properly.
id::DOF_ID
AtomTree::set_torsion_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	AtomID const & atom4,
	Real const setting
)
{

	Real offset;
	DOF_ID const & id
		( torsion_angle_dof_id( atom1, atom2, atom3, atom4, offset ) );

	if ( !id.valid() ) {
		// couldnt find this angle
		return id;
	}

	// note: offset is defined (see torsion_angle_dof_id) so that
	//
	// torsion(atom1,atom2,atom3,atom4) = dof(id) + offset
	//

	set_dof( id, setting - offset );

	return id;
}

/////////////////////////////////////////////////////////////////////////////

id::DOF_ID
AtomTree::set_bond_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	Real const setting
)
{

	Real offset(0.0);
	DOF_ID const dof_id( bond_angle_dof_id( atom1, atom2, atom3, offset ) );

	if ( dof_id.valid() ) {
		assert( offset == 0.0 ); // not handling this case right now

		set_dof( dof_id, numeric::constants::d::pi - setting );
	}

	return dof_id;
}

/////////////////////////////////////////////////////////////////////////////

Real
AtomTree::bond_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3
) const
{

	Real offset(0.0);
	DOF_ID const dof_id( bond_angle_dof_id( atom1, atom2, atom3, offset ) );

	if ( dof_id.valid() ) {
		assert( offset == 0.0 ); // not handling this case right now
		return numeric::constants::d::pi - dof( dof_id );
	} else {
		TR << "unable to find DOF_ID for bond_angle: " << atom1 << ' ' << atom2 << ' ' << atom3 << std::endl;
	}

	return 0.0;
}

/////////////////////////////////////////////////////////////////////////////

id::DOF_ID
AtomTree::bond_angle_dof_id(
	AtomID const & atom1_in_id,
	AtomID const & atom2_in_id,
	AtomID const & atom3_in_id,
	Real & offset
) const
{
	offset = 0.0;

	// CASES
	// I.   the bond angle is the THETA dof of atom3 (atom3's parent is atom2 and atom2's parent is atom1)
	// II.  the bond angle is the THETA dof of atom1 (atom1's parent is atom2 and atom2's parent is atom3)
	// III. the bond angle is set by a torsion offset (atom1's parent is atom2 and atom3's parent is atom2)

	assert( atom_pointer( atom1_in_id ) && atom_pointer( atom2_in_id ) && atom_pointer( atom3_in_id ) );

	Atom const
		* const atom1_in( atom_pointer( atom1_in_id ) ),
		* const atom2_in( atom_pointer( atom2_in_id ) ),
		* const atom3_in( atom_pointer( atom3_in_id ) );

	Atom const * atom1( atom1_in ), *atom2( atom2_in ), *atom3( atom3_in );

	// reorder the atoms if necessary
	// not necessary at all in the current logic
	DOF_ID dof_id( id::BOGUS_DOF_ID );

	if ( !atom3->is_jump() &&
			 !atom3->keep_dof_fixed( THETA ) && // not stub atom2 of a jump
			 atom3->input_stub_atom1() == atom2 &&
			 atom3->input_stub_atom2() == atom1 ) {
		// case I
		dof_id = DOF_ID( atom3->id(), THETA );

	} else if ( !atom1->is_jump() &&
							!atom1->keep_dof_fixed( THETA ) && // not stub atom2 of a jump
							atom1->input_stub_atom1() == atom2 &&
							atom1->input_stub_atom2() == atom3 ) {
		// case II, perfect case in reverse, want to reverse the order of the atoms
		atom1 = atom3_in;
		atom2 = atom2_in;
		atom3 = atom1_in;
		dof_id = DOF_ID( atom3->id(), THETA );

	} else {
		// not handling other cases right now
	}

	return dof_id;
}

/////////////////////////////////////////////////////////////////////////////

id::DOF_ID
AtomTree::set_bond_length(
	AtomID const & atom1,
	AtomID const & atom2,
	Real const setting
)
{

	DOF_ID const dof_id( bond_length_dof_id( atom1, atom2 ) );

	if ( dof_id.valid() ) {
		set_dof( dof_id, setting );
	}

	return dof_id;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Real
AtomTree::bond_length(
	AtomID const & atom1,
	AtomID const & atom2
) const
{

	DOF_ID const dof_id( bond_length_dof_id( atom1, atom2 ) );

	if ( dof_id.valid() ) {
		return dof( dof_id );
	} else {
		TR << "unable to find DOF_ID for bond_length: " << atom1 << ' ' << atom2 << std::endl;
	}

	return 0.0;
}

/////////////////////////////////////////////////////////////////////////////

id::DOF_ID
AtomTree::bond_length_dof_id(
	AtomID const & atom1_id,
	AtomID const & atom2_id
) const
{

	assert( atom_pointer( atom1_id ) && atom_pointer( atom2_id ) );

	Atom const * atom1( atom_pointer( atom1_id ) ),* atom2( atom_pointer( atom2_id ) );

	if ( !atom2->is_jump() &&
			 atom2->input_stub_atom1() == atom1 ) {
		// case I
		return DOF_ID( atom2->id(), D );

	} else if ( !atom1->is_jump() &&
							atom1->input_stub_atom1() == atom2 ) {
		// case II
		return DOF_ID( atom1->id(), D );

	} else {
		// not handling other cases right now
	}

	return id::BOGUS_DOF_ID;
}

/////////////////////////////////////////////////////////////////////////////
///
/// @note this is NOT calculated straight from atom positions.Instead, it is
/// calculated from interal DOFs. If no DOF can be found from these four atoms,
/// 0.0 will be returned and a warning message is printed.
Real
AtomTree::torsion_angle(
	AtomID const & atom1,
	AtomID const & atom2,
	AtomID const & atom3,
	AtomID const & atom4
) const
{
 	update_internal_coords();

	// find the atomtree degree of freedom that corresponds to this torsion
	// angle
	Real offset;
	DOF_ID const & id
		( torsion_angle_dof_id( atom1, atom2, atom3, atom4, offset ) );

	if ( !id.valid() ) {
		// couldnt find this angle
		TR << "AtomTree::torsion_angle() cant find dof! " <<
			atom1 << ' ' << atom2 << ' ' << atom3 << ' ' << atom4 << std::endl;
		return 0.0;
	}

	// note: offset is defined (see torsion_angle_dof_id) so that
	//
	// torsion(atom1,atom2,atom3,atom4) = dof(id) + offset
	//

	return atom_pointer_[ id.atom_id() ]->dof( id.type() ) + offset;
}

/////////////////////////////////////////////////////////////////////////////

///
/// @details This is done by releasing memories of all atoms including the root atom
///, and clearing atom_pointer map
void
AtomTree::clear()
{
	atom_pointer_.clear();
	root_ = 0;
	external_coordinate_residues_changed_->total_residue(0);

	// anything that depends on the tree topology needs to be updated
	set_new_topology();

}

/////////////////////////////////////////////////////////////////////////////
/// @note this may trigger a coordinate update of src (access to root atom)
///

void
AtomTree::copy_coords(
	AtomTree const & src
)
{
	if ( !root_ ) {
		if ( src.root() ) utility_exit_with_message("AtomTree::copy_coords: I'm empty but src is not!");
		return;
	}
	internal_coords_need_updating_ = src.internal_coords_need_updating_;
	xyz_coords_need_updating_ = src.xyz_coords_need_updating_;

	root_->copy_coords( *(src.root()) );
	dof_changeset_ = src.dof_changeset_;
	(*external_coordinate_residues_changed_) = (*src.external_coordinate_residues_changed_);

}

////////////////////////////////////////////////////////////////////////////////////
/// @details domain map is residue-based and indicate which residues stay relatively rigid
/// to each other in one domain since last move and therefore their interaction
/// energy does not have to be recalculated. This function calls atom->update_domain_map
/// from the root atom of the tree.
void
AtomTree::update_domain_map(
	DomainMap & domain_map,
	AtomID_Mask const & dof_moved,
	AtomID_Mask const & xyz_moved
) const
{
	domain_map = -1;

	int current_color(1), biggest_color(1);
	root_->update_domain_map( current_color, biggest_color, domain_map,
		dof_moved, xyz_moved );
}


/////////////////////////////////////////////////////////////////////////////
// private
void
AtomTree::update_atom_ids_from_atom_pointer()
{
	for ( Size i=1, i_end = atom_pointer_.size(); i<= i_end; ++i ) {
		for ( Size j=1, j_end = atom_pointer_[i].size(); j<= j_end; ++j ) {
			assert( atom_pointer_[i][j] );
			if ( atom_pointer_[i][j] ) atom_pointer_[i][j]->id( AtomID(j,i) );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
void
AtomTree::update_sequence_numbering(
	Size const new_size,
	utility::vector1< int > const & old2new
)
{
	/// ResidueCoordinateChangeList is not setup to handle a remapping.  It must
	/// be empty, which means the Conformation its tracking moved data for must have
	/// already retrieved its moved data.
	assert( external_coordinate_residues_changed_->empty() );
	external_coordinate_residues_changed_->total_residue( new_size );

	atom_pointer_.update_sequence_numbering( new_size, old2new );

	update_atom_ids_from_atom_pointer();
}


/////////////////////////////////////////////////////////////////////////////
/// @details makes a complete copy of another AtomTree src by cloning the root_atom and
/// all its offspring atoms. Also update atom_pointer map.
/// @note clear the content of tree first to release memory before doing the assignment.
AtomTree &
AtomTree::operator=( AtomTree const & src )
{
	if (this == &src) {
		return *this;
	}

	if ( topological_match_to_ == &src || src.topological_match_to_ == this ) {
		/// Why include the second condition: src.topological_match_to_ == this ?
		/// As an optimization for the common case when atom trees are being copied
		/// back and forth into each other as happens in repeated calls to MC::boltzman.
		/// Moreover, "topological match to" is a commutative property.
		copy_coords( src );
	} else {
		// no memory leaks -- notify observers that the topology is changing.
		clear();

		// just to get the dimensions right:
		atom_pointer_ = src.atom_pointer_;

		// copy tree (recursive)
		if ( src.root_ ) {
			root_ = src.root_->clone( 0 /*parent=0*/, atom_pointer_ );
		}

		internal_coords_need_updating_ = src.internal_coords_need_updating_;
		xyz_coords_need_updating_ = src.xyz_coords_need_updating_;
		dof_changeset_ = src.dof_changeset_;

		(*external_coordinate_residues_changed_) = (*src.external_coordinate_residues_changed_);

	}

	if ( topological_match_to_ != & src ) {
		if ( topological_match_to_ != 0 ) {
			topological_match_to_->detatch_topological_observer( this );
		}
		src.attach_topological_observer( this );
	}
	return *this;

}

/////////////////////////////////////////////////////////////////////////////
/// @details update_internal_coords would be called in situations where we want
/// to ensure that the internal degrees are valid, eg when we are about
/// to set a new internal DOF.
/// @note private method, const to allow lazy updating of
/// @note see usage in AtomTree.cc
/// @note these guys are const because of the lazy updating scheme we are using:
/// they need to be called from within const accessing functions.

void
AtomTree::update_internal_coords() const
{
	// this would be bad:
	assert( ! ( xyz_coords_need_updating_ && internal_coords_need_updating_ ) );

	if ( internal_coords_need_updating_ ) {
		if ( !root_ ) utility_exit_with_message("phil how did we get here?-1");

		PROF_START( basic::ATOM_TREE_UPDATE_INTERNAL_COORDS ); // profiling
		root_->update_internal_coords( default_stub );
		PROF_STOP ( basic::ATOM_TREE_UPDATE_INTERNAL_COORDS );

		internal_coords_need_updating_ = false;
	}
}



/// @details  Set the transform between two stubs
/// Returns the atomid of the jump atom which moved.
/// @note  Requires that there be a jump in the atomtree between a pair of atoms in the stubs
id::AtomID
AtomTree::set_stub_transform(
	StubID const & stub_id1,
	StubID const & stub_id2,
	RT const & target_rt
)
{

	// look for the connection between these two stubs
	Atom * jump_atom(0);
	int dir(0);
	for ( Size i=1; i<= 3; ++i ) {
		Atom * atom1( atom_pointer( stub_id1.atom( i ) ) );
		for ( Size j=1; j<= 3; ++j ) {
			Atom * atom2( atom_pointer( stub_id2.atom( j ) ) );
			if ( atom1->is_jump() && atom1->parent() == atom2 ) {
				assert( !jump_atom );
				jump_atom = atom1;
				dir = -1;
			} else if ( atom2->is_jump() && atom2->parent() == atom1 ) {
				assert( !jump_atom );
				jump_atom = atom2;
				dir = 1;
			}
		}
	}
	if ( !jump_atom ) {
		utility_exit_with_message( "AtomTree::set_stub_transform: No jump between these atoms!" );
	}

	Stub const  instub( dir == 1 ? stub_from_id( stub_id1 ) : stub_from_id( stub_id2 ) );
	Stub const outstub( dir == 1 ? stub_from_id( stub_id2 ) : stub_from_id( stub_id1 ) );
	RT rt( target_rt );
	if ( dir == -1 ) rt.reverse();
	// now solve for a linear transform that will move outstub so that RT( instub, outstub ) == rt
	Stub::Matrix A;
	Vector b;
	find_stub_transform( instub, outstub, rt, A, b );
	jump_atom->transform_Ax_plus_b_recursive( A, b, *external_coordinate_residues_changed_ );
	internal_coords_need_updating_ = true; // this could be more efficient!

	// confirm that things worked
	assert( RT( stub_from_id( stub_id1 ), stub_from_id( stub_id2 ) ).distance_squared( target_rt ) < 1e-3 );
	return jump_atom->id();
}

/// @brief  get the transform between two stubs
RT
AtomTree::get_stub_transform(
	StubID const & stub_id1,
	StubID const & stub_id2
) const
{
	return RT( stub_from_id( stub_id1 ), stub_from_id( stub_id2 ) );
}
/////////////////////////////////////////////////////////////////////////////
// private
///\brief update xyz coordinates from internal cooridnates
void
AtomTree::update_xyz_coords() const
{
	// this would be bad:
	assert( ! ( xyz_coords_need_updating_ && internal_coords_need_updating_ ) );


	if ( xyz_coords_need_updating_ ) {
		if ( !root_ ) utility_exit_with_message("phil how did we get here-2?");

		PROF_START( basic::ATOM_TREE_UPDATE_XYZ_COORDS ); // profiling
		for ( Size ii = 1; ii <= dof_changeset_.size(); ++ii ) {
			if ( dof_changeset_[ ii ].reached_ ) continue;
			atom_pointer_[ dof_changeset_[ ii ].atomid_ ]->dfs( dof_changeset_, *external_coordinate_residues_changed_, ii );
		}

		for ( Size ii = 1; ii <= dof_changeset_.size(); ++ii ) {
			if ( dof_changeset_[ ii ].reached_ ) continue;
			//std::cout << "Refold from " << dof_changeset_[ ii ].atomid_.rsd() << std::endl; // << " " << dof_changeset_[ ii ].atomid_.atomno() << " " << dof_changeset_[ ii ].reached_ << std::endl;
			atom_pointer_[ dof_changeset_[ ii ].atomid_ ]->update_xyz_coords(); // it must find its own stub.
		}
		dof_changeset_.clear();
		PROF_STOP ( basic::ATOM_TREE_UPDATE_XYZ_COORDS );

		xyz_coords_need_updating_ = false;

		//std::cout << "REFOLD END" << std::endl;

	}
}

/// @brief The AtomTree provides to the Conformation object a list of residues
/// whose xyz coordinates have changed.  When the Conformation has finished reading off
/// residues that have changed from the AtomTree, and has copied the coordinates of
/// those residues into its conformation::Residue objects, it informs the AtomTree
/// to reset this list by a call to mark_changed_residues_registered
///
/// @details The list of which residues have had coordinate changes is unknown until
/// the DFS has completed.  The DFS must be triggered before iterators are given to the
/// Conformation object.
ResidueListIterator
AtomTree::residue_xyz_change_list_begin() const
{
	if ( xyz_coords_need_updating_ ) update_xyz_coords();
	return external_coordinate_residues_changed_->residues_moved_begin();
}

/// @details The list of which residues have had coordinate changes is unknown until
/// the DFS has completed.  The DFS must be triggered before iterators are given to the
/// Conformation object.
ResidueListIterator
AtomTree::residue_xyz_change_list_end() const
{
	if ( xyz_coords_need_updating_ ) update_xyz_coords();
	return external_coordinate_residues_changed_->residues_moved_end();
}


/// @brief The AtomTree provides a list of residues who's xyz coordinates have changed
/// to the Conformation object.  When the Conformation has finished reading off residues
/// that have changed from the AtomTree, and has copied the coordinates of those residues
/// into its conformation::Residue objects, it informs the AtomTree to reset this list
/// by a call to mark_changed_residues_registered
void
AtomTree::note_coordinate_change_registered() const
{
	external_coordinate_residues_changed_->clear();
}


/// @brief When an atom tree copies the topology of another atom tree, it must
/// register itself as a topological observer of that other tree.  When the other
/// tree changes its topology, then the tree being observed must notify its
/// observers that they are no longer a topological copy of this tree.  An atom
/// tree may only be the topological copy of a single other atom tree, though several
/// atom trees may be copies of a single atom tree.
void
AtomTree::attach_topological_observer( AtomTree const * observer ) const
{
	assert( observer->topological_match_to_ == 0 );
	bool resize_topo_observers_array_( false );
	for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
		if ( topological_observers_[ ii ] == 0 ) {
			resize_topo_observers_array_ = true;
			/// Make sure the observer is not already observing me
			assert( topological_observers_[ ii ] != observer );
		}
	}

	if ( resize_topo_observers_array_ ) {
		Size n_valid( 0 );
		for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
			if ( topological_observers_[ ii ] != 0 ) ++n_valid;
		}
		utility::vector1< AtomTree const * > valid_observers;
		valid_observers.reserve( n_valid + 1 );
		for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
			if ( topological_observers_[ ii ] != 0 ) valid_observers.push_back( topological_observers_[ ii ] );
		}
		topological_observers_.swap( valid_observers );
	}

	topological_observers_.push_back( observer );
	observer->topological_match_to_ = this;
}

/// @details The AtomTree being observed calls this notify_topological_change on all
/// of the AtomeTrees that are observing it.  The this object "detaches" itself from
/// the observee in this function.  In debug mode, this function ensures that this
/// object was actually observing the observee -- if it wasn't, then an internal error
/// has occurred.  The observee and the observer have gotten out of sync.
void
AtomTree::notify_topological_change( AtomTree const * ASSERT_ONLY( observee ) ) const
{
	assert( observee == topological_match_to_ );
	topological_match_to_ = 0;
}

/// @details When an AtomTree which was observing this AtomTree changes its
/// topology or is deleted, it must invoke this method on this AtomTree.
/// When this happens, this AtomTree marks the observer's position in
/// its list of observers as null.
void
AtomTree::detatch_topological_observer( AtomTree const * observer ) const
{
	assert( observer->topological_match_to_ == this );
	bool found( false );
	for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
		if ( topological_observers_[ ii ] == observer ) {
			found = true;
			topological_observers_[ ii ] = 0;
			observer->topological_match_to_ = 0;
			break;
		}
	}
	assert( found );
}

/// @details If this tree changes its topology, then its no longer a match in topology
/// to the tree it was copied from, nor are the trees that are a copy of it.  Detatch
/// this as an observer of the AtomTree it is observing (if any) and inform all the
/// trees observing this tree that the topology has changed before clearing the
/// topological_observers_ array.
void
AtomTree::set_new_topology()
{
	if ( topological_match_to_ != 0 ) {
		topological_match_to_->detatch_topological_observer( this );
	}
	for ( Size ii = 1; ii <= topological_observers_.size(); ++ii ) {
		if ( topological_observers_[ ii ] != 0 ) {
			topological_observers_[ ii ]->notify_topological_change( this );
		}
	}
	topological_observers_.clear();
}


////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// CARTESIAN-COORDINATE FRAGMENT INSERTION ROUTINES ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void
AtomTree::get_frag_atoms(
	StubID const & id,
	FragXYZ const & frag_xyz,
	Atom const * & frag_atom,
	Atom const * & nonfrag_atom // could be zero if atom1 is root of tree
) const
{
	bool const atom1_in_frag( frag_xyz.count( id.atom1 ) );
	bool const atom2_in_frag( frag_xyz.count( id.atom2 ) );

	Atom const * atom1( atom_pointer( id.atom1 ) );

	if ( atom1->is_jump() ) {
		if ( !atom1_in_frag && atom1->parent() && frag_xyz.count( atom1->parent()->atom_id() ) ) {
			frag_atom    = atom1->parent();
			nonfrag_atom = atom1;
			return;
		} else if ( atom1_in_frag && ( !atom1->parent() || !frag_xyz.count( atom1->parent()->atom_id() ) ) ) {
			frag_atom    = atom1;
			nonfrag_atom = atom1->parent();
			return;
		}
	}

	if ( atom1_in_frag && !atom2_in_frag ) {
		frag_atom    = atom_pointer( id.atom1 );
		nonfrag_atom = atom_pointer( id.atom2 );
		return;

	} else if ( atom2_in_frag && !atom1_in_frag ) {
		frag_atom    = atom_pointer( id.atom2 );
		nonfrag_atom = atom_pointer( id.atom1 );
		return;

	}
	utility_exit_with_message( "AtomTree::get_frag_atoms failed" );

}


/////////////////////////////////////////////////////////////////////////////
/// id is a frag atom
///
/// look for two more nearby frag atoms to define a pseudo-stub for getting coords for parents or children of this atom
///

id::StubID
AtomTree::get_frag_pseudo_stub_id(
	AtomID const & id,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	assert( frag_xyz.count( id ) );

	utility::vector1< AtomID > ids;

	Atom const * atom( atom_pointer( id ) );
	Atom const *parent( atom->parent() ), *child1( atom->get_nonjump_atom(0) ), *child2( atom->get_nonjump_atom(1) );

	if ( !atom->is_jump() && parent && frag_xyz.count( parent->atom_id() ) ) {
		ids.push_back( parent->atom_id() );
	}
	if ( child1 && frag_xyz.count( child1->atom_id() ) ) {
		ids.push_back( child1->atom_id() );
	}
	if ( child2 && frag_xyz.count( child2->atom_id() ) ) {
		ids.push_back( child2->atom_id() );
	}
	if ( child1 && frag_xyz.count( child1->atom_id() ) ) {
		Atom const * gchild1( child1->get_nonjump_atom(0) );
		if ( gchild1 && frag_xyz.count( gchild1->atom_id() ) ) {
			ids.push_back( gchild1->atom_id() );
		}
	}
	if ( child2 && frag_xyz.count( child2->atom_id() ) ) {
		Atom const * gchild2( child2->get_nonjump_atom(0) );
		if ( gchild2 && frag_xyz.count( gchild2->atom_id() ) ) {
			ids.push_back( gchild2->atom_id() );
		}
	}
	if ( !atom->is_jump() && parent && frag_xyz.count( parent->atom_id() ) ) {
		Atom const * gparent( parent->parent() );
		if ( !parent->is_jump() && gparent && frag_xyz.count( gparent->atom_id() ) ) {
			ids.push_back( gparent->atom_id() );
		}
		for ( Size i=0; i<parent->n_children(); ++i ) {
			Atom const * sibling( parent->child(i) );
			if ( sibling != atom && !sibling->is_jump() && frag_xyz.count( sibling->atom_id() ) ) {
				ids.push_back( sibling->atom_id() );
			}
		}
	}
	if ( ids.size() < 2 ) {
		fail = true;
		return StubID( id, id, id );
	}
	return StubID( id, ids[1], ids[2] );
}



/////////////////////////////////////////////////////////////////////////////
/// private helper for fragment insertion routines
///
Stub
AtomTree::get_frag_local_stub(
	StubID const & stubid,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	// first look for special case of stub of frag jump-child:
	AtomID const stub_atom1_id( stubid.atom1 ), stub_atom2_id( stubid.atom2 ), stub_atom3_id( stubid.atom3 );
	Atom const * stub_atom1( atom_pointer( stub_atom1_id ) );

	if ( !frag_xyz.count( stub_atom1_id ) && stub_atom1->is_jump() &&
			 stub_atom1->parent() && frag_xyz.count( stub_atom1->parent()->atom_id()) &&
			 stub_atom1->stub_atom1_id() == stub_atom1_id &&
			 stub_atom1->stub_atom2_id() == stub_atom2_id &&
			 stub_atom1->stub_atom3_id() == stub_atom3_id ) {

		/// special case: handled more easily this way
		StubID const instub_id( stub_atom1->input_stub_atom1_id(),
														stub_atom1->input_stub_atom2_id(),
														stub_atom1->input_stub_atom3_id() );
		assert( stub_atom1->input_stub_atom0_id() == stub_atom1->input_stub_atom1_id() );

		// current xyz transform:
		RT const current_rt( stub_from_id( instub_id ), stub_from_id( stubid ) );
		Stub const instub( get_frag_local_stub( instub_id, frag_xyz, fail ) ); // recursive call
		Stub outstub;
		TR.Trace << "get_frag_local_stub:: making jump: " <<
			stubid.atom1 << ' ' << stubid.atom2 << ' ' << stubid.atom3 << ' ' <<
			instub_id.atom1 << ' ' << instub_id.atom2 << ' ' << instub_id.atom3 << std::endl;

		current_rt.make_jump( instub, outstub );
		return outstub;
	}


	return Stub( get_frag_local_xyz( stub_atom1_id, frag_xyz, fail ),
							 get_frag_local_xyz( stub_atom2_id, frag_xyz, fail ),
							 get_frag_local_xyz( stub_atom3_id, frag_xyz, fail ) );
}

/////////////////////////////////////////////////////////////////////////////
/// private helper for fragment insertion routines
///
/// id is either in frag or a child or a gchild or a parent
///

Vector
AtomTree::get_frag_local_xyz(
	AtomID const & id,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	// easiest case:
	if ( frag_xyz.count( id ) ) return frag_xyz.find( id )->second;

	Atom const * atom( atom_pointer(id) );

	if ( ( atom->parent() && frag_xyz.count( atom->parent()->atom_id() ) ) ||
			 ( atom->parent() && atom->parent()->parent() && frag_xyz.count( atom->parent()->parent()->atom_id() ) ) ) {
		// child or grand child
		return get_frag_descendant_local_xyz( atom, frag_xyz, fail );
	}

	// now should be parent of frag
	Atom const * child( 0 );
	for ( Size i=0; i< atom->n_children(); ++i ) {
		if ( frag_xyz.count( atom->child(i)->atom_id() ) ) {
			assert( child == 0 );
			child = atom->child(i);
		}
	}
	assert( child ); // this is a known hole: could be asked for a grandparent of a frag, not just a parent. fix this phil
	return get_frag_parent_local_xyz( child, frag_xyz, fail );


}
/////////////////////////////////////////////////////////////////////////////
/// private helper for fragment insertion routines
///
Vector
AtomTree::get_frag_descendant_local_xyz(
	Atom const * atom,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	AtomID const id( atom->atom_id() );
	assert( !frag_xyz.count( id ) );
	bool const frag_child( atom->parent() && frag_xyz.count( atom->parent()->atom_id() ) );
	ASSERT_ONLY(bool const frag_gchild
		( atom->parent() && atom->parent()->parent() && frag_xyz.count( atom->parent()->parent()->atom_id() ) );)
	assert( frag_child || frag_gchild );

	AtomID id1( atom->input_stub_atom1_id() );
	AtomID id2( atom->input_stub_atom2_id() );
	AtomID id3( atom->input_stub_atom3_id() );
	assert( atom->input_stub_atom0_id() == id1 ); // cant handle this case yet

	if ( id == id1 || id == id2 || id == id3 ) {
		/// circular!!! potential for infinite loop!
		if ( frag_child ) {
			StubID tmp( get_frag_pseudo_stub_id( atom->parent()->atom_id(), frag_xyz, fail ) );
			id1 = tmp.atom1;
			id2 = tmp.atom2;
			id3 = tmp.atom3;
		} else {
			StubID tmp( get_frag_pseudo_stub_id( atom->parent()->parent()->atom_id(), frag_xyz, fail ) );
			id1 = tmp.atom1;
			id2 = tmp.atom2;
			id3 = tmp.atom3;
		}
		if ( fail ) return Vector(0.0);
	}

	Stub const current_stub( stub_from_id( StubID( id1, id2, id3) ) );
	Stub const   local_stub( get_frag_local_stub( StubID( id1, id2, id3 ), frag_xyz, fail ) ); // potentially recursive
	return local_stub.local2global( current_stub.global2local( xyz( id ) ) );
}

/////////////////////////////////////////////////////////////////////////////
/// private helper for fragment insertion routines
///
Vector
AtomTree::get_frag_parent_local_xyz(
	Atom const * child,
	FragXYZ const & frag_xyz,
	bool & fail
) const
{
	Atom const * parent( child->parent() );
	//assert( !parent->is_jump() ); // dont think we have to handle this case...
	assert( frag_xyz.count( child->atom_id() ) && ! frag_xyz.count( parent->atom_id() ) );

	// build a pseudo-stub for the parent
	StubID const pseudo_stubid( get_frag_pseudo_stub_id( child->atom_id(), frag_xyz, fail ) );
	if ( fail ) return Vector(0.0);
	Stub const current_stub( stub_from_id( pseudo_stubid ) );
	Stub const   local_stub( get_frag_local_stub( pseudo_stubid, frag_xyz, fail ) );
	assert( !fail ); // since pseudo stub atoms are all in the fragment
	return local_stub.local2global( current_stub.global2local( xyz( parent->atom_id() ) ) );
}



/////////////////////////////////////////////////////////////////////////////
// a private function
//
// the incoming stub should in fact be the true incoming connection
// and all the outgoing connections should be accounted for.
//
// This is all checked in assert statements.
//
void
AtomTree::insert_single_fragment(
	StubID const & instub_id,
	FragRT const & outstub_transforms,
	FragXYZ const & frag_xyz,
	utility::vector1< AtomID > & moving_atoms
)
{
	// get the incoming stub:
	Stub const instub( stub_from_id( instub_id ) );

	Atom const * instub_frag_atom(0), *instub_nonfrag_atom(0);
	get_frag_atoms( instub_id, frag_xyz, instub_frag_atom, instub_nonfrag_atom );

	assert( instub_frag_atom->parent() == instub_nonfrag_atom ); // sanity check

	utility::vector1< Atom const * > outstub_nonfrag_atoms; // just for debugging

	for ( FragRT::const_iterator it= outstub_transforms.begin(), ite= outstub_transforms.end(); it != ite; ++it ) {
		StubID const & outstub_id( it->first );
		Atom const * outstub_frag_atom(0), *outstub_nonfrag_atom(0);
		get_frag_atoms( outstub_id, frag_xyz, outstub_frag_atom, outstub_nonfrag_atom );
		outstub_nonfrag_atoms.push_back( outstub_nonfrag_atom ); // for debugging
		assert( outstub_nonfrag_atom->parent() == outstub_frag_atom );

		// now transform outstub_nonfrag_atom so that current transform becomes desired transform
		//
		Stub const & stub1( instub ); // to match names below
		Stub const   stub2( stub_from_id( outstub_id ) );
		RT const & rt( it->second );

		// now things are arranged so that when we transform outstub_nonfrag_atom and children, stub1 stays fixed and
		// stub2 moves, and we want RT( stub1, new_stub2 ) == target_rt

		Stub::Matrix const & M1( stub1.M ), M2( stub2.M ), R( rt.get_rotation() );
		Vector const & v1( stub1.v ), v2( stub2.v ), t( rt.get_translation() );

		// look for a transformation of the form x |----> A*x + b
		//
		// this will change stub2 to stub2' with M2' = A * M2, v2' = A*v2 + b
		//
		// if we let (R,t) be the target RT, then we want
		//
		//  R = M1^T * M2' = M1^T * A * M2  ==> A = M1 * R * M2^T
		//
		//  t = M1^T * ( v2' - v1 ) ==> v2' = M1 * t + v1, which with b = v2' - A*v2 gives b = M1 * t + v1 - A * v2
		//

		Stub::Matrix const A( M1 * R * M2.transposed() );
		Vector const b( M1 * t + v1 - A * v2 );

		atom_pointer_[ outstub_nonfrag_atom->atom_id() ]->
			transform_Ax_plus_b_recursive( A, b, *external_coordinate_residues_changed_ ); // get nonconst version
		internal_coords_need_updating_ = true;

		moving_atoms.push_back( outstub_nonfrag_atom->atom_id() );
	}


	for ( FragXYZ::const_iterator it=frag_xyz.begin(), ite= frag_xyz.end(); it != ite; ++it ) {
		AtomID const & id( it->first );
		Atom * atom( atom_pointer( id ) );

		// update xyz using instub
		atom->xyz( instub.local2global( it->second ) );

		{ // sanity check/debug
			// if parent not in frag, assert this is instub_frag_atom
			assert( ( atom->parent() && frag_xyz.count( atom->parent()->atom_id() ) ) || atom == instub_frag_atom );

			// if has a child not in frag, assert child is in outstub_nonfrag_atoms
			for ( Size i=0; i< atom->n_children(); ++i ) {
				assert( frag_xyz.count( atom->child(i)->atom_id() ) ||
								( std::find( outstub_nonfrag_atoms.begin(), outstub_nonfrag_atoms.end(), atom->child(i) ) !=
									outstub_nonfrag_atoms.end()));
			}
		}
		moving_atoms.push_back( id );

	}
	internal_coords_need_updating_ = true;


	// now debug transforms
	for ( FragRT::const_iterator it= outstub_transforms.begin(), ite= outstub_transforms.end(); it != ite; ++it ) {
		assert( RT( instub, stub_from_id( it->first )).distance_squared( it->second ) < 1e-3 );
	}


}


/////////////////////////////////////////////////////////////////////////////
void
AtomTree::set_jump_atom_stub_id(
	StubID const & id
)
{
	update_xyz_coords(); // since we will need to recalculate all the internal dofs when we're done

	Atom * atom1( atom_pointer( id.atom1 ) );
	Atom * atom2( atom_pointer( id.atom2 ) );
	Atom * atom3( atom_pointer( id.atom3 ) );
	if ( !atom1->is_jump() || atom2->is_jump() || atom3->is_jump() ||
			 atom2->parent() != atom1 || atom3->parent() != atom2 ) {
		utility_exit_with_message( "set_jump_atom_stub_id failed!" );
	}

	atom1->delete_atom( atom2 );
	atom1->insert_atom( atom2 ); // goes to head of line
	atom2->delete_atom( atom3 );
	atom2->insert_atom( atom3 ); // goes to head of line

	internal_coords_need_updating_ = true;

	assert( atom1->stub_atom1() == atom1 && atom1->stub_atom2() == atom2 && atom1->stub_atom3() == atom3 );

}


/////////////////////////////////////////////////////////////////////////////
// completely new plan
//
void
AtomTree::insert_fragment(
	StubID const & instub_id,
	FragRT const & outstub_transforms,
	FragXYZ const & frag_xyz,
	utility::vector1< AtomID > & moving_atoms
)
{
	// first compile a list of incoming/outgoing connections

	// look for more incoming and/or outgoing connections:

	utility::vector1< Atom const * > incoming_stub_frag_atoms, outgoing_stub_nonfrag_atoms;
	utility::vector1< Stub > incoming_local_stubs, outgoing_local_stubs;
	utility::vector1< StubID > incoming_stub_ids, outgoing_stub_ids;

	{
		Atom const * frag_atom, *nonfrag_atom;
		get_frag_atoms( instub_id, frag_xyz, frag_atom, nonfrag_atom );
		if ( frag_atom->parent() == nonfrag_atom ) {
			TR.Trace << "AtomTree::insert_fragment: instub is incoming" << std::endl;
			incoming_stub_frag_atoms.push_back( frag_atom );
			incoming_local_stubs.push_back( Stub() );
			incoming_stub_ids.push_back( instub_id );
		} else if ( nonfrag_atom && nonfrag_atom->parent() == frag_atom ) {
			TR.Trace << "AtomTree::insert_fragment: instub is outgoing" << std::endl;
			outgoing_stub_nonfrag_atoms.push_back( nonfrag_atom );
			outgoing_local_stubs.push_back( Stub() );
			outgoing_stub_ids.push_back( instub_id );
		}
	} // scope

	for ( FragRT::const_iterator it= outstub_transforms.begin(), ite= outstub_transforms.end(); it != ite; ++it ) {
		StubID const & outstub_id( it->first );
		Atom const * frag_atom, *nonfrag_atom;
		get_frag_atoms( outstub_id, frag_xyz, frag_atom, nonfrag_atom );
		if ( nonfrag_atom && nonfrag_atom->parent() == frag_atom ) {
			TR.Trace << "AtomTree::insert_fragment: outstub is outgoing" << std::endl;
			outgoing_stub_nonfrag_atoms.push_back( nonfrag_atom );
			outgoing_local_stubs.push_back( Stub( it->second ) );
			outgoing_stub_ids.push_back( outstub_id );
		} else if ( frag_atom->parent() == nonfrag_atom ) {
			TR.Trace << "AtomTree::insert_fragment: outstub is incoming" << std::endl;
			incoming_stub_frag_atoms.push_back( frag_atom );
			incoming_local_stubs.push_back( Stub( it->second ) );
			incoming_stub_ids.push_back( outstub_id );
		}
	}

	for ( FragXYZ::const_iterator it=frag_xyz.begin(), ite= frag_xyz.end(); it != ite; ++it ) {
		//AtomID const & id( it->first );
		Atom * atom( atom_pointer( it->first ) );
		if ( ( atom->parent() == 0 || !frag_xyz.count( atom->parent()->atom_id() ) ) &&
				 ( std::find( incoming_stub_frag_atoms.begin(), incoming_stub_frag_atoms.end(), atom ) ==
					 incoming_stub_frag_atoms.end() ) ) {
			// found a new incoming connection!
			TR.Trace << "AtomTree::insert_fragment: found new incoming connection: " << atom->atom_id() << std::endl;
			// now want to generate a StubID for this guy as well as a local stub
			if ( atom->is_jump() ) {
				// incoming jump-atom connection -- or root of tree??
				// get local coords for all the stub atoms
				bool fail( false );
				StubID const stubid( atom->stub_atom1_id(), atom->stub_atom2_id(), atom->stub_atom3_id() );
				Stub const stub( get_frag_local_stub( stubid, frag_xyz, fail ) );
				if ( !fail ) {
					incoming_stub_frag_atoms.push_back( atom );
					incoming_local_stubs.push_back( stub );
					incoming_stub_ids.push_back( stubid );
				}
				assert( !fail ); // could fail
			} else {
				// incoming bonded-atom connection
				bool fail( false );
				StubID const stubid( get_frag_pseudo_stub_id( atom->atom_id(), frag_xyz, fail ) );
				if ( !fail ) {
					incoming_stub_frag_atoms.push_back( atom );
					incoming_stub_ids.push_back( stubid );
					incoming_local_stubs.push_back( get_frag_local_stub( stubid, frag_xyz, fail ) );
					assert( !fail ); // shouldnt fail for a frag_pseudo_stub -- all atoms are already in frag
				}
				assert( !fail ); // could fail
			}
		}

		// look for outgoing connections:
		for ( Size i=0; i< atom->n_children(); ++i ) {
			Atom* child( atom->child(i) );
			if ( !frag_xyz.count( child->atom_id() ) &&
					 ( std::find( outgoing_stub_nonfrag_atoms.begin(), outgoing_stub_nonfrag_atoms.end(), child ) ==
						 outgoing_stub_nonfrag_atoms.end() ) ) {
				// new outgoing connection!
				TR.Trace << "AtomTree::insert_fragment: found new outgoing connection: " << atom->atom_id() << std::endl;
				bool fail( false );
				StubID const stubid( child->stub_atom1_id(), child->stub_atom2_id(), child->stub_atom3_id() );
				Stub const outstub( get_frag_local_stub( stubid, frag_xyz, fail ) );
				if ( !fail ) {
					outgoing_stub_nonfrag_atoms.push_back( child );
					outgoing_stub_ids.push_back( stubid );
					outgoing_local_stubs.push_back( outstub );
				}
				assert( !fail );
			}
		}
	}


	/// Now associate outgoing stubs and frag atoms with the incoming stubs, call insert_single_fragment once for each
	/// incoming stub.
	///
	for ( Size i=1; i<= incoming_stub_frag_atoms.size(); ++i ) {
		Atom const * instub_atom( incoming_stub_frag_atoms[i] );
		Stub const & local_instub( incoming_local_stubs[ i ] );
		FragXYZ new_frag_xyz;
		FragRT new_outstub_transforms;
		// get frag atoms that depend on this guy:
		for ( FragXYZ::const_iterator it=frag_xyz.begin(), ite= frag_xyz.end(); it != ite; ++it ) {
			AtomID const & id( it->first );
			Atom const * frag_atom( atom_pointer( id ) );
			if ( frag_atom->atom_is_on_path_from_root( instub_atom ) ) {
				new_frag_xyz[ id ] = local_instub.global2local( it->second );
			}
		}
		for ( Size j=1; j<= outgoing_stub_nonfrag_atoms.size(); ++j ) {
			Atom const * outstub_atom( outgoing_stub_nonfrag_atoms[j] );
			if ( outstub_atom->atom_is_on_path_from_root( instub_atom ) ) {
				new_outstub_transforms[ outgoing_stub_ids[j] ] = RT( local_instub, outgoing_local_stubs[j] );
			}
		}

		TR.Trace << "AtomTree::insert_fragment: inserting single fragment nout= " << new_outstub_transforms.size() <<
			" natoms= " << new_frag_xyz.size() << std::endl;

		insert_single_fragment( incoming_stub_ids[i], new_outstub_transforms, new_frag_xyz, moving_atoms );
	}


}


/// @details  Useful for guaranteeing that a stub remains within a single residue
void
AtomTree::promote_sameresidue_nonjump_child( AtomID const & parent_atom_id )
{
	update_xyz_coords(); // since we are about to invalidate the internal coords

	Atom * parent( atom_pointer( parent_atom_id ) );
	tree::Atom * sameresidue_child( 0 );
	for ( Size i=0; i< parent->n_nonjump_children(); ++i ) {
		tree::Atom * child( atom_pointer( parent->get_nonjump_atom( i )->id() ) ); // want nonconst, use atom_pointer
		if ( child->id().rsd() == parent->id().rsd() ) {
			sameresidue_child = child;
			break;
		}
	}
	if ( sameresidue_child ) {
		assert( !sameresidue_child->is_jump() );
		parent->delete_atom( sameresidue_child );
		parent->insert_atom( sameresidue_child );
	} else {
		TR.Warning << "promote_sameresidue_nonjump_child failed, parent has no non-jump, same-residue children!" <<
			std::endl;
	}

	internal_coords_need_updating_ = true;
	set_new_topology();
}




// 			while ( it != path1.end() &&
// 							std::find( stub1_ids.begin(), stub1_ids.end(), (*it)->parent()->atom_id() ) != stub1_ids.end() ) ++it;
// 			if ( it == path1.end() ) {
// 				it = path2.begin();
// 				while ( it != path2.end() &&
// 								std::find( stub2_ids.begin(), stub2_ids.end(), (*it)->parent()->atom_id() ) != stub2_ids.end() ) ++it;
// 				if ( it == path2.end() ) {
// 					utility_exit_with_message( "AtomTree::make_stub_transform: Unable to find bond to cut!" );
// 				}
// 			}
// 		}


} // namespace kinematics
} // namespace core
// 		// we use the internal dof's if necessary to calculate the offset
// 		assert( !internal_coords_need_updating_ );

// 		offset = 0.0;

// 		assert( atom_pointer_[ atom_id1 ] && atom_pointer_[ atom_id2 ] &&
// 						atom_pointer_[ atom_id3 ] && atom_pointer_[ atom_id4 ] );

// 		// STUART -- (low priority) I'd like to be able to cache
// 		// the results of this calculation to allow faster access.


// 		// reorder the atoms if necessary so that atom2 is the parent of atom3
// 		// and atom3 is the parent of atom4
// 		Atom *atom1( 0 ), *atom2( 0 ), *atom3( 0 ), *atom4( 0 );

// 		if ( atom_pointer_[ atom_id3 ]->parent()->atom_id() == atom_id2 ) {
// 			atom1 = atom_pointer_[ atom_id1 ];
// 			atom2 = atom_pointer_[ atom_id2 ];
// 			atom3 = atom_pointer_[ atom_id3 ];
// 			atom4 = atom_pointer_[ atom_id4 ];
// 		} else if ( atom_pointer_[ atom_id2 ]->parent()->atom_id() == atom_id3 ) {
// 			atom1 = atom_pointer_[ atom_id4 ];
// 			atom2 = atom_pointer_[ atom_id3 ];
// 			atom3 = atom_pointer_[ atom_id2 ];
// 			atom4 = atom_pointer_[ atom_id1 ];
// 		} else {
// 			// no good!
// 			return BOGUS_DOF_ID;
// 		}
// 		assert( atom3->parent() == atom2 );
// 		if ( atom4->parent() != atom3 ) {
// 			// also no good
// 			return BOGUS_DOF_ID;
// 		} else if ( atom4->is_jump() ) {
// 			// also no good
// 			return BOGUS_DOF_ID;
// 		}


// 		Atom const
// 			*dof_atom1( atom2->parent() ),
// 			*dof_atom4( atom3->get_nonjump_atom(0) ); // guaranteed not a jump

// 		// this is the answer:
// 		DOF_ID dof_id( dof_atom4->atom_id(), PHI );


// 		///////////////////////////////////////////////////////////////////////////
// 		// handle the easiest case: perfect match with an atomtree DOF
// 		//
// 		if ( dof_atom1 == atom1 &&
// 				 dof_atom4 == atom4 ) {
// 			assert( &(atom4->input_stub_atom0()) == atom3 &&
// 							&(atom4->input_stub_atom1()) == atom3 &&
// 							&(atom4->input_stub_atom2()) == atom2 &&
// 							&(atom4->input_stub_atom3()) == atom1 );

// 			offset = 0.0;
// 			return dof_id;
// 		}

// 		// now calculate the offset
// 		// by summing contributions from both ends
// 		offset = 0.0;

// 		// handle offset between atom4 and atom3's first child
// 		// note that we already know that atom4 is one of atom3's children
// 		// ( see above )
// 		//
// 		if ( dof_atom4 != atom4 ) {

// 			// recall: torsion(atom1-4) = dof(dof_id) + offset
// 			//
// 			offset += atom3->dihedral_between_bonded_children( dof_atom4, atom4 );
// 		}


// 		// handle offset between atom1 and true atom1
// 		// the only case we can do is if atom1->parent() == atom2
// 		//
// 		// this will happen, eg for chi1 if we are folding c2n:
// 		//
// 		// atom1 =  N    while    dof_atom1 = parent(atom2) =  C
// 		// atom2 = CA
// 		// atom3 = CB
// 		// atom4 = CG
// 		//
// 		if ( dof_atom1 != atom1 ) {
// 			if ( atom1->parent() != atom2 || atom1->is_jump() ) {
// 				// no good!
// 				return BOGUS_DOF_ID;
// 			}

// 			// a little tricky: we don't want to access the positions,
// 			// since we can't be sure that they are up to date in a routine
// 			// that would be called inside dof setting and getting routines
// 			//
// 			//
// 			// need to use some spherical geometry
// 			//
// 			// we've got three unit vectors:
// 			//
// 			// u1 -- from atom2 to atom1
// 			// u2 -- from dof_atom1 to atom2
// 			// u3 -- from atom2 to atom3
// 			//
// 			// the angle between u2 and u1 = theta1 = atom1(THETA)
// 			// the angle between u2 and u3 = theta3 = atom3(THETA)
// 			// the torsion between u1 and u3 along u2 = phi2 = atom2->bonded_dih(1,3)
// 			//
// 			Real
// 				theta1( atom1->dof( THETA ) ),
// 				theta3( atom3->dof( THETA ) ),
// 				phi2( atom2->dihedral_between_bonded_children( atom1, atom3 ) );

// 			// want to solve for phi3 -- dihedral between u2 and u1 along axis
// 			// defined by u3
// 			//
// 			// spherical law of cosines says:
// 			//
// 			// cos(theta2) = cos(theta1) * cos(theta3) +
// 			//               sin(theta1) * sin(theta3) * cos(phi2)
// 			//
// 			// so we can solve for cos(theta2); then use formula again to solve for
// 			// cos(phi3).
// 			//
// 			Real cos_theta2 = std::cos( theta1 ) * std::cos( theta3 ) +
// 				std::sin( theta1 ) * std::sin( theta3 ) * std::cos( phi2 );

// 			Real sin_theta2 = std::sqrt( 1 - cos_theta2 * cos_theta2 );

// 			// use formula again interchanging 2 and 3
// 			Real cos_phi3 = ( cos( theta3 ) - cos( theta1 ) * cos_theta2 ) /
// 				( sin( theta1 ) * sin_theta2 );

// 			// PHIL! really should handle degenerate cases more carefully
// 			// PHIL! also could reuse some of the trig calculations


// 			std::cout << "PHIL: calculate the sign factor in AtomTree::torsion_angle_dof_id!" << std::endl;

// 			Real sign_factor( 1.0 ); // this is not correct

// 			offset += sign_factor * std::acos( cos_phi3 );

// 		}

// 		// offset is passed out by reference
// 		return dof_id;
// 	}
// /////////////////////////////////////////////////////////////////////////////
// /// tmp hack

// class AtomBondChecker {
// public:
// 	typedef id::BondID BondID;

// public:

// 	AtomBondChecker( utility::vector1< BondID > const & bonds_in ):
// 		bonds_( bonds_in )
// 	{}


// 	bool
// 	operator()( Atom const * const & atom )
// 	{
// 		BondID id( atom->atom_id(), atom->parent()->atom_id() );
// 		if ( std::find( bonds_.begin(), bonds_.end(), id ) == bonds_.end() ) id.reverse();

// 		// return true if bond between atom and parent is in the list:
// 		return ( std::find( bonds_.begin(), bonds_.end(), id ) != bonds_.end() );
// 	}

// private:
// 	utility::vector1< BondID > bonds_;


// };

// 																 vector1< id::StubID > const & outstub_ids,
// 																 vector1< RT > const & outstub_transforms,
// 																 vector1< AtomID > const & frag_ids,
// 																 vector1< Vector > const & frag_xyz, // or use std::map< AtomID, Vector > ??
// 	AtomID id1( child.atom_id() ),id2,id3;
// 	if ( child->n_nonjump_children() == 0 ) {
// 		fail = true;
// 		return Vector(0.0);
// 	} else if ( child->n_nonjump_children() == 1 ) {
// 		Atom const * gchild( child->get_nonjump_atom( 0 ) ); // frag or child of frag
// 		id2 = gchild->atom_id();
// 		if ( gchild->n_nonjump_children() == 0 ) {
// 			fail = true;
// 			return Vector(0.0);
// 		} else {
// 			id3 = gchild->get_nonjump_atom(0)->atom_id(); // frag or child of frag or grandchild of frag
// 		}
// 	} else {
// 		id2 = child->get_nonjump_atom(0)->atom_id(); // frag or child of frag
// 		id3 = child->get_nonjump_atom(1)->atom_id(); // frag or child of frag
// 	}
// /////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////
// // private, called recursively
// //
// void
// AtomTree::insert_stub_transforms(
// 																 id::StubID const & instub_id,
// 																 utility::vector1< id::StubID > const & outstub_ids,
// 																 utility::vector1< RT > const & transforms,
// 																 id::AtomID_Mask & exclude,
// 																 utility::vector1< id::AtomID > & moved_atoms
// )
// {
// 	assert( outstub_ids.size() == transforms.size() );
// 	Size const nstub( outstub_ids.size() );

// 	// build vector of stub atom ids

// 	// trim stub atom ids to exclude overlapping atoms

// 	vector1< vector1< Atom const * > > paths( nstub );
// 	for ( Size i=1; i<= nstub; ++i ) {
// 		get_atom_path( out_atom_ids[i][1], in_atom_ids[1], paths[i] );
// 		sizes.push_back( paths[i].size () );
// 	}

// 	///
// 	Size const stub_index( argmin( sizes ) ); // see rotamer_trials.cc
// 	vector1< Atom const * > stub_path( paths[ stub_index ] );


// 	/// choose an atom to move ///////////////////////////////////////////////
// 	Atom * moving_atom( 0 );
// 	// first look for jump atom:
// 	for ( Size i=1; i<= stub_path.size(); ++i ) {
// 		AtomID const & id( stub_path[i]->atom_id() );
// 		if ( !excluded[ id ] && stub_path[i]->is_jump() ) {
// 			moving_atom = atom_pointer_[ id ];
// 			break;
// 		}
// 	}
// 	for ( Size i=1; i<= stub_path.size(); ++i ) {
// 		AtomID const & id( stub_path[i]->atom_id() );
// 		if ( !excluded[ id ] ) {
// 			excluded[ id ] = true;
// 			if ( !moving_atom ) moving_atom = atom_pointer_[ id ];
// 			break;
// 		}
// 	}
// 	if ( !moving_atom ) utility_exit_with_message( "AtomTree stub fragment insertion failed!!" );

// 	bool const moving_atom_on_incoming_root_path( my_path.back()->atom_is_on_root_path( moving_atom ) );
// 	assert( moving_atom_on_incoming_root_path || my_path.front()->atom_is_on_root_path( moving_atom ) );

// 	//// now figure out what transform we need to apply
// 	Stub stub1( xyz( stub1_id.atom1 ), xyz( stub1_id.atom2 ), xyz( stub1_id.atom3 ) );
// 	Stub stub2( xyz( stub2_id.atom1 ), xyz( stub2_id.atom2 ), xyz( stub2_id.atom3 ) );
// 	RT rt( target_rt );
// 	if ( moving_atom_on_path1 ) {
// 		// in this case when we transform the atom we will actually move the stub1 atoms, so reverse for consistency
// 		rt.reverse();
// 		Stub const tmp( stub1 );
// 		stub1 = stub2;
// 		stub2 = tmp;
// 	}

// 	// now things are arranged so that when we transform moving_atom and children, stub1 stays fixed and stub2 moves.
// 	// and we want RT( stub1, new_stub2 ) == target_rt

// 	Stub::Matrix const & M1( stub1.M ), M2( stub2.M ), R( rt.get_rotation() );
// 	Vector const & v1( stub1.v ), v2( stub2.v ), t( rt.get_translation() );

// 	// look for a transformation of the form x |----> A*x + b
// 	//
// 	// this will change stub2 to stub2' with M2' = A * M2, v2' = A*v2 + b
// 	//
// 	// if we let (R,t) be the target RT, then we want
// 	//
// 	//  R = M1^T * M2' = M1^T * A * M2  ==> A = M1 * R * M2^T
// 	//
// 	//  t = M1^T * ( v2' - v1 ) ==> v2' = M1 * t + v1, which with b = v2' - A*v2 gives b = M1 * t + v1 - A * v2
// 	//

// 	Stub::Matrix const A( M1 * R * M2.transposed() );
// 	Vector const b( M1 * t + v1 - A * v2 );

// 	update_xyz_coords();

// 	moving_atom->transform_Ax_plus_b_recursive( A, b );
// 	internal_coords_need_updating_ = true;


// 	{ // debug
// 		Stub const new_stub1( xyz( stub1_id.atom1 ), xyz( stub1_id.atom2 ), xyz( stub1_id.atom3 ) );
// 		Stub const new_stub2( xyz( stub2_id.atom1 ), xyz( stub2_id.atom2 ), xyz( stub2_id.atom3 ) );
// 		RT const new_rt( new_stub1, new_stub2 );
// 		std::cout << "debugRTD: " << new_rt.distance_squared( target_rt ) << std::endl;

// 		assert( new_rt.distance_squared( target_rt ) < 1e-3 );
// 	}

// 	//// tell the outside world which atom we moved (ie subtree rooted at this atom has moved)
// 	return moving_atom->atom_id();

// }


// /////////////////////////////////////////////////////////////////////////////

// id::AtomID
// AtomTree::make_stub_transform(
// 	id::StubID const & stub1_id, // triplet of atomids
// 	id::StubID const & stub2_id, // triplet of atomids
// 	RT const & target_rt,
// 	utility::vector1< id::BondID > const & preferred_bonds
// )
// {
// 	//using tree::Atom;

// 	//// get path between origin atoms of both stubs /////////////////////
// 	utility::vector1< Atom const * > path1, path2;

// 	Atom* const stub1_atom1( atom_pointer_[ stub1_id.atom1 ] );
// 	Atom* const stub2_atom1( atom_pointer_[ stub2_id.atom1 ] );


// 	stub1_atom1->get_path_from_root( path1 );
// 	stub2_atom1->get_path_from_root( path2 );

// 	assert( path1.front() == root_       && path2.front() == root_ &&
// 					path1.back () == stub1_atom1 && path2.back () == stub2_atom1 );

// 	//// Now remove all the common ancestors, and reverse the paths so they start at stub origin atoms
// 	{
// 		Size lcai(1); // last_common_ancestor_index
// 		Size const s1( path1.size() );
// 		Size const s2( path2.size() );
// 		while ( lcai < s1 && lcai < s2 ) {
// 			++lcai;
// 			if ( path1[ lcai ] != path2[ lcai ] ) {
// 				--lcai;
// 				break;
// 			}
// 		}
// 		path1.erase( path1.begin(), path1.begin() + lcai );
// 		path2.erase( path2.begin(), path2.begin() + lcai );
// 		assert( path1.empty() || path2.empty() ||
// 						( path1[ 1 ] != path2[ 1 ] && path1[ 1 ]->parent() == path2[ 1 ]->parent() ) );
// 		std::reverse( path1.begin(), path1.end() );
// 		std::reverse( path2.begin(), path2.end() );
// 		assert( path1.empty() || path1[ 1 ] == stub1_atom1 );
// 		assert( path2.empty() || path2[ 1 ] == stub2_atom1 );


// 		// trim backward to get past all stub atoms
// 		utility::vector1< AtomID > stub1_ids, stub2_ids;
// 		stub1_ids.push_back( stub1_id.atom1 );
// 		stub1_ids.push_back( stub1_id.atom2 );
// 		stub1_ids.push_back( stub1_id.atom3 );
// 		stub2_ids.push_back( stub2_id.atom1 );
// 		stub2_ids.push_back( stub2_id.atom2 );
// 		stub2_ids.push_back( stub2_id.atom3 );
//  		while ( !path1.empty() &&
//  						std::find( stub1_ids.begin(), stub1_ids.end(), path1.front()->parent()->atom_id() ) != stub1_ids.end() ) {
//  			path1.erase( path1.begin() );
//  		}
//  		while ( !path2.empty() &&
//  						std::find( stub2_ids.begin(), stub2_ids.end(), path2.front()->parent()->atom_id() ) != stub2_ids.end() ) {
//  			path2.erase( path2.begin() );
//  		}
// 	}


// 	// use this to find atoms whose bond to a parent is in preferred_bond_ids
// 	AtomBondChecker bond_checker( preferred_bonds );

// 	//// first look for a jump on either path that matches the preferred_bonds set

// 	// make paths with just the jump atoms
// 	utility::vector1< Atom const * > path1_jumps, path2_jumps;
// 	for ( Size i=1; i<= path1.size(); ++i ) if ( path1[i]->is_jump() ) path1_jumps.push_back( path1[i] );
// 	for ( Size i=1; i<= path2.size(); ++i ) if ( path2[i]->is_jump() ) path2_jumps.push_back( path2[i] );

// 	utility::vector1< Atom const * >::iterator it = find_if( path1_jumps.begin(), path1_jumps.end(), bond_checker );
// 	if ( it == path1_jumps.end() )             it = find_if( path2_jumps.begin(), path2_jumps.end(), bond_checker );
// 	if ( it == path2_jumps.end() )             it = find_if(       path1.begin(),       path1.end(), bond_checker );
// 	if ( it ==       path1.end() )             it = find_if(       path2.begin(),       path2.end(), bond_checker );
// 	if ( it ==       path2.end() ) {
// 		// no atoms matched the preferred bonds set
// 		// choose a jump if it exists
// 		// otherwise choose ...

// 		if ( !path1_jumps.empty() ) {
// 			it = path1_jumps.begin();
// 		} else if ( ! path2_jumps.empty() ) {
// 			it = path2_jumps.begin();
// 		} else if ( !path1.empty() ) {
// 			it = path1.begin();
// 		} else if ( !path2.empty() ) {
// 			it = path2.begin();
// 		} else {
// 			utility_exit_with_message( "AtomTree::make_stub_transform: Unable to find bond to cut!" );
// 		}
// 	}

// 	Atom* moving_atom( atom_pointer_[ (*it)->atom_id() ] ); // get nonconst version
// 	bool const moving_atom_on_path1( std::find( path1.begin(), path1.end(), *it ) != path1.end() );;

// 	//// now figure out what transform we need to apply
// 	Stub stub1( xyz( stub1_id.atom1 ), xyz( stub1_id.atom2 ), xyz( stub1_id.atom3 ) );
// 	Stub stub2( xyz( stub2_id.atom1 ), xyz( stub2_id.atom2 ), xyz( stub2_id.atom3 ) );
// 	RT rt( target_rt );
// 	if ( moving_atom_on_path1 ) {
// 		// in this case when we transform the atom we will actually move the stub1 atoms, so reverse for consistency
// 		rt.reverse();
// 		Stub const tmp( stub1 );
// 		stub1 = stub2;
// 		stub2 = tmp;
// 	}

// 	// now things are arranged so that when we transform moving_atom and children, stub1 stays fixed and stub2 moves.
// 	// and we want RT( stub1, new_stub2 ) == target_rt

// 	Stub::Matrix const & M1( stub1.M ), M2( stub2.M ), R( rt.get_rotation() );
// 	Vector const & v1( stub1.v ), v2( stub2.v ), t( rt.get_translation() );

// 	// look for a transformation of the form x |----> A*x + b
// 	//
// 	// this will change stub2 to stub2' with M2' = A * M2, v2' = A*v2 + b
// 	//
// 	// if we let (R,t) be the target RT, then we want
// 	//
// 	//  R = M1^T * M2' = M1^T * A * M2  ==> A = M1 * R * M2^T
// 	//
// 	//  t = M1^T * ( v2' - v1 ) ==> v2' = M1 * t + v1, which with b = v2' - A*v2 gives b = M1 * t + v1 - A * v2
// 	//

// 	Stub::Matrix const A( M1 * R * M2.transposed() );
// 	Vector const b( M1 * t + v1 - A * v2 );

// 	update_xyz_coords();

// 	moving_atom->transform_Ax_plus_b_recursive( A, b );
// 	internal_coords_need_updating_ = true;


// 	{ // debug
// 		Stub const new_stub1( xyz( stub1_id.atom1 ), xyz( stub1_id.atom2 ), xyz( stub1_id.atom3 ) );
// 		Stub const new_stub2( xyz( stub2_id.atom1 ), xyz( stub2_id.atom2 ), xyz( stub2_id.atom3 ) );
// 		RT const new_rt( new_stub1, new_stub2 );
// 		std::cout << "debugRTD: " << new_rt.distance_squared( target_rt ) << std::endl;

// 		assert( new_rt.distance_squared( target_rt ) < 1e-3 );
// 	}

// 	//// tell the outside world which atom we moved (ie subtree rooted at this atom has moved)
// 	return moving_atom->atom_id();

// }


