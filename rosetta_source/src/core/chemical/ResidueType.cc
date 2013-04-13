// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @file ResidueType.cc
///
/// @brief
/// A class for defining a type of residue
///
/// @details
/// This class contains the "chemical" information for residues. This does not contain the actual xyz coordinates of a
/// particular residue in a specific peptide.  (xyz coordinates are found in core/conformation/Residue.hh).  A residue
/// in Rosetta can be a ligand, DNA, amino acid, or basically anything.  A residue is read in through residue_io.cc and
/// read from parameter files, generally located in the database chemical/residue_types.  For ligands, or anything that
/// is not one of the natural 20 AAs, a parameter has to be provided to rosetta through the -extra_res_fa flag.
/// residue_io.cc sets private member data in ResidueType.  The primary data that are set are: atoms, mmatoms,
/// orbitals, and properties of the particular residue type.  These properties can be modified through patches, which
/// is controlled through PatchOperations.cc.  If the residue_type of a residue is modified, the indices of atoms and
/// mmatoms and everything associated with those indices must be redefined.  This reordering of indices is taken care
/// of with the function reorder_primary_data().
///
/// Setting of primary data and then reordering is important.  Primary data for the following are described:
///
/// Atoms: Setting of atoms includes indexing the atoms into vectors, saving their names into vectors/maps, saving the
/// associated mm_atom_type into a vector, saving bond connections into vectors, etc, etc.  Since everything is
/// allocated into vectors, it is easy to reorder those vectors.  On any given residue, the heavy atoms are put into
/// the vector first, (their indices are first,) and hydrogens are put in last.
///
/// Properties: Properties of a residue include things like DNA, PROTEIN, SC_ORBITALS, CHARGED, etc.  These properties
/// indicate the type of residue it is and what properties are associated with the residue.  They are set when read in.
/// Several lines of code must be modified to get them to work, all found here in ResidueType.cc.
///
/// Orbitals: Orbitals are indexed separately from atoms.  They function much the same way as atoms, except for some
/// key differences.  To find atoms bonded to orbitals, you must provide the atom index, not the orbital index.  (I
/// haven't figured out how to get the reverse to work because of the separate indices.)  Orbital xyz coordinates are
/// not updated when atom coordinates are.  This is to keep speed consistent with just having atoms.  To output the
/// orbitals, use the flag -output_orbitals.
///
/// @author
/// Phil Bradley
/// Steven Combs - these comments
////////////////////////////////////////////////////////////////////////

// Unit headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>

// Package Headers
#include <core/conformation/Residue.hh>

// Project Headers
#include <core/chemical/ResidueSupport.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/NumericTraits.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/PyAssert.hh>

// Options and Option key includes (needed for protonated versions of the residues - pH mode)
#include <basic/options/option.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/VariantType.hh>
#include <utility/vector1.hh>
#include <algorithm>

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end

namespace core {
namespace chemical {

using namespace ObjexxFCL;
using namespace ObjexxFCL::fmt;

static basic::Tracer tr("core.chemical.ResidueType");

/// must be a better place for this, probably already exists!
inline
std::string
strip_whitespace( std::string const & name )
{
	std::string trimmed_name( name );
	left_justify( trimmed_name ); trim( trimmed_name ); // simpler way to dothis?
	return trimmed_name;
}

ResidueType::ResidueType(
	AtomTypeSetCAP atom_types,
	ElementSetCAP elements,
	MMAtomTypeSetCAP mm_atom_types,
	orbitals::OrbitalTypeSetCAP orbital_types//, CSDAtomTypeSetCAP csd_atom_types kwk commenting out csd atom types until they are fully functional
) :
	utility::pointer::ReferenceCount(),
	atom_types_( atom_types ),
	elements_( elements ),
	mm_atom_types_( mm_atom_types ),
	orbital_types_( orbital_types),
	residue_type_set_( 0 ),
	graph_(),
	orbitals_(),
	natoms_(0),
	nheavyatoms_(0),
	n_hbond_acceptors_(0),
	n_hbond_donors_(0),
	n_orbitals_(0),
	n_backbone_heavyatoms_(0),
	first_sidechain_hydrogen_( 0 ),
	ndihe_( 0 ),
	rotamer_library_name_( "" ),
	use_ncaa_rotlib_( false ),
	ncaa_rotlib_n_rots_( 0 ),
	is_polymer_( false ),
	is_protein_( false ),
	is_charged_( false ),
	is_polar_( false ),
	has_sc_orbitals_(false),
	is_aromatic_( false ),
	is_DNA_( false ),
	is_RNA_( false ),
	is_NA_( false ),
	is_carbohydrate_( false ),
	is_ligand_( false ),
	is_surface_( false ),
	is_terminus_( false ),
	is_lower_terminus_( false ),
	is_upper_terminus_( false ),
	is_acetylated_nterminus_( false ),
	is_methylated_cterminus_( false ),
	is_coarse_( false ), //currently for coarse_RNA only
	is_adduct_( false ),
	aa_( aa_unk ),
	rotamer_aa_( aa_unk ),
	name_(),
	name3_(),
	name1_(),
	nbr_atom_(1),
	nbr_radius_( 0 ),
	force_nbr_atom_orient_(false),
	molecular_mass_(0),
	molar_mass_(0),
	n_actcoord_atoms_( 0 ),
	lower_connect_id_( 0 ),
	upper_connect_id_( 0 ),
	n_non_polymeric_residue_connections_( 0 ),
	n_polymeric_residue_connections_( 0 ),
	carbohydrate_info_(NULL),
	finalized_(false),
	nondefault_(false),
	base_restype_name_(""),
	serialized_(false)
{
}

ResidueType::~ResidueType()
{
	tr.Trace << "Residue dstor" << std::endl;
}

ResidueType::ResidueType(ResidueType const & residue_type):
	utility::pointer::ReferenceCount(),
	atom_types_( residue_type.atom_types_ ),
	elements_( residue_type.elements_ ),
	mm_atom_types_( residue_type.mm_atom_types_ ),
	orbital_types_( residue_type.orbital_types_ ),
	residue_type_set_( residue_type.residue_type_set_ ),
	graph_(residue_type.graph_),
	orbitals_(residue_type.orbitals_),
	natoms_(residue_type.natoms_),
	nheavyatoms_(residue_type.nheavyatoms_),
	n_hbond_acceptors_(residue_type.n_hbond_acceptors_),
	n_hbond_donors_(residue_type.n_hbond_donors_),
	n_orbitals_(residue_type.n_orbitals_),
	n_backbone_heavyatoms_(residue_type.n_backbone_heavyatoms_),
	first_sidechain_hydrogen_( residue_type.first_sidechain_hydrogen_ ),
	ndihe_( residue_type.ndihe_ ),
	orbital_bonded_neighbor_(residue_type.orbital_bonded_neighbor_),
	bonded_neighbor_(residue_type.bonded_neighbor_),
	bonded_neighbor_type_(residue_type.bonded_neighbor_type_),
	cut_bond_neighbor_(residue_type.cut_bond_neighbor_),
	atom_base_(residue_type.atom_base_),
	abase2_(residue_type.abase2_), // acceptors only
	attached_H_begin_(residue_type.attached_H_begin_),
	attached_H_end_(residue_type.attached_H_end_),
	parents_(residue_type.parents_),
	dihedral_atom_sets_(residue_type.dihedral_atom_sets_),
	dihedrals_for_atom_(residue_type.dihedrals_for_atom_),
	bondangle_atom_sets_(residue_type.bondangle_atom_sets_),
	bondangles_for_atom_(residue_type.bondangles_for_atom_),
	last_controlling_chi_(residue_type.last_controlling_chi_),
	atoms_last_controlled_by_chi_(residue_type.atoms_last_controlled_by_chi_),
	atoms_with_orb_index_(residue_type.atoms_with_orb_index_),
	Haro_index_(residue_type.Haro_index_),
	Hpol_index_(residue_type.Hpol_index_),
	accpt_pos_(residue_type.accpt_pos_),
	Hpos_polar_(residue_type.Hpos_polar_),
	Hpos_apolar_(residue_type.Hpos_apolar_),
	accpt_pos_sc_(residue_type.accpt_pos_sc_),
	Hpos_polar_sc_(residue_type.Hpos_polar_sc_),
	all_bb_atoms_(residue_type.all_bb_atoms_),
	all_sc_atoms_(residue_type.all_sc_atoms_),
	heavyatom_has_polar_hydrogens_(residue_type.heavyatom_has_polar_hydrogens_),
	heavyatom_is_an_acceptor_(residue_type.heavyatom_is_an_acceptor_),
	atom_is_polar_hydrogen_(residue_type.atom_is_polar_hydrogen_),
	mainchain_atoms_(residue_type.mainchain_atoms_),
	actcoord_atoms_(residue_type.actcoord_atoms_),
	chi_atoms_(residue_type.chi_atoms_),
	is_proton_chi_(residue_type.is_proton_chi_),
	proton_chis_(residue_type.proton_chis_),
	chi_2_proton_chi_(residue_type.chi_2_proton_chi_),
	proton_chi_samples_(residue_type.proton_chi_samples_),
	proton_chi_extra_samples_(residue_type.proton_chi_extra_samples_),
	path_distance_(residue_type.path_distance_),
	atom_graph_index_(), /// This must be regenerated below to hold the new new vertex_descriptors
	ordered_atoms_(), // This must be regenerated to hold the new vertex_descriptors
	orbitals_index_(residue_type.orbitals_index_),
	chi_rotamers_(residue_type.chi_rotamers_),
	rotamer_library_name_( residue_type.rotamer_library_name_ ),
	use_ncaa_rotlib_( residue_type.use_ncaa_rotlib_ ),
	ncaa_rotlib_path_( residue_type.ncaa_rotlib_path_),
	ncaa_rotlib_n_rots_( residue_type.ncaa_rotlib_n_rots_ ),
	ncaa_rotlib_n_bins_per_rot_(residue_type.ncaa_rotlib_n_bins_per_rot_),
	properties_(residue_type.properties_),
	is_polymer_( residue_type.is_polymer_ ),
	is_protein_( residue_type.is_protein_ ),
	is_charged_( residue_type.is_charged_ ),
	is_polar_( residue_type.is_polar_ ),
	has_sc_orbitals_(residue_type.has_sc_orbitals_),
	is_aromatic_( residue_type.is_aromatic_ ),
	is_DNA_( residue_type.is_DNA_ ),
	is_RNA_( residue_type.is_RNA_ ),
	is_NA_( residue_type.is_NA_ ),
	is_carbohydrate_( residue_type.is_carbohydrate_ ),
	is_ligand_( residue_type.is_ligand_ ),
	is_surface_( residue_type.is_surface_ ),
	is_terminus_( residue_type.is_terminus_ ),
	is_lower_terminus_( residue_type.is_lower_terminus_ ),
	is_upper_terminus_( residue_type.is_upper_terminus_ ),
	is_acetylated_nterminus_( residue_type.is_acetylated_nterminus_ ),
	is_methylated_cterminus_( residue_type.is_methylated_cterminus_ ),
	is_coarse_( residue_type.is_coarse_ ), //currently for coarse_RNA only
	is_adduct_( residue_type.is_adduct_ ),
	variant_types_( residue_type.variant_types_ ),
	aa_( residue_type.aa_ ),
	rotamer_aa_( residue_type.rotamer_aa_ ),
	name_( residue_type.name_),
	name3_( residue_type.name3_),
	name1_(residue_type.name1_),
	nbr_atom_(residue_type.nbr_atom_),
	nbr_radius_( residue_type.nbr_radius_ ),
	force_nbr_atom_orient_(residue_type.force_nbr_atom_orient_),
	molecular_mass_(residue_type.molecular_mass_),
	molar_mass_(residue_type.molar_mass_),
	n_actcoord_atoms_( residue_type.n_actcoord_atoms_ ),
	mol_data_(residue_type.mol_data_),
	residue_connections_(residue_type.residue_connections_),
	atom_2_residue_connection_map_(residue_type.atom_2_residue_connection_map_),
	atoms_within_one_bond_of_a_residue_connection_(residue_type.atoms_within_one_bond_of_a_residue_connection_),
	within1bonds_sets_for_atom_(residue_type.within1bonds_sets_for_atom_),
	atoms_within_two_bonds_of_a_residue_connection_(residue_type.atoms_within_two_bonds_of_a_residue_connection_),
	within2bonds_sets_for_atom_(residue_type.within2bonds_sets_for_atom_),
	lower_connect_id_( residue_type.lower_connect_id_ ),
	upper_connect_id_( residue_type.upper_connect_id_ ),
	n_non_polymeric_residue_connections_( residue_type.n_non_polymeric_residue_connections_ ),
	n_polymeric_residue_connections_( residue_type.n_polymeric_residue_connections_ ),
	branch_point_atoms_(residue_type.branch_point_atoms_),
	delete_atoms_(residue_type.delete_atoms_),
	force_bb_(residue_type.force_bb_),
	rna_residuetype_(residue_type.rna_residuetype_),
	carbohydrate_info_(residue_type.carbohydrate_info_),
	finalized_(residue_type.finalized_),
	defined_adducts_(residue_type.defined_adducts_),
	nondefault_(residue_type.nondefault_),
	base_restype_name_(residue_type.base_restype_name_),
	serialized_(residue_type.serialized_)
{
	// Setup the atom_graph_index_ and map old to new...
	std::map<VD, VD> old_to_new;
	for(
			VIterPair vp = boost::vertices(graph_), old_vp= boost::vertices(residue_type.graph_);
			vp.first != vp.second;
			++vp.first, ++old_vp.first
	){
		VIter v_iter= vp.first;
		VIter old_v_iter= old_vp.first;
		VD vd = *v_iter;
		VD old_vd = *old_v_iter;
		old_to_new[old_vd] = vd;//Assuming the boost::graph copy preserves ordering within the vertices list
		Atom & a = graph_[vd];
		assert( a == residue_type.graph_[old_vd]);
		NameVDInserted const name_vd_inserted = atom_graph_index_.insert(NameVDPair(a.name(), vd));
		assert(name_vd_inserted.second); // Don't add atoms with the same name
		atom_graph_index_.insert( NameVDPair( strip_whitespace( a.name() ), vd) );
		//assert(strip_name_vd_inserted.second); // If this is 4 chars, than it will be the same as before.
	}
	/// Setup the temporary ordered_atoms_ vector for refactor
	VDs::const_iterator begin = residue_type.ordered_atoms_.begin();
	VDs::const_iterator const end = residue_type.ordered_atoms_.end();
	for(; begin != end; ++begin){
		VD old_vd = *begin;
		VD vd = old_to_new[old_vd];
		ordered_atoms_.push_back(vd);
	}
}

///
ResidueTypeSet const &
ResidueType::residue_type_set() const
{
	if ( residue_type_set_ == 0 ) {
		utility_exit_with_message( "ResidueType::residue_type_set: pointer is not set!");
	}
	return *residue_type_set_;
}

//////////////////////////////////////////////////////////////////////////////

/// make a copy
ResidueTypeOP
ResidueType::clone() const
{
	ResidueTypeOP rsd_ptr( new ResidueType( *this ) );
	return rsd_ptr;
}

///
void
ResidueType::residue_type_set( ResidueTypeSetCAP set_in )
{
	residue_type_set_ = set_in;
}

Atom & ResidueType::atom(Size const atom_index){
	return graph_[ ordered_atoms_[atom_index] ];
}
Atom const & ResidueType::atom(Size const atom_index) const{
	return graph_[ ordered_atoms_[atom_index] ];
}
Atom & ResidueType::atom(std::string const & atom_name){
	return graph_[ atom_graph_index_[atom_name] ];
}
Atom const & ResidueType::atom(std::string const & atom_name) const{
	NameVDMap::const_iterator found = atom_graph_index_.find( atom_name );
	assert( found != atom_graph_index_.end());
	return graph_[  found->second ];
}
Orbital const & ResidueType::orbital(Size const orbital_index) const{
	return orbitals_[orbital_index];
}
Orbital const & ResidueType::orbital(std::string const & orbital_name) const{
	return orbitals_[ orbital_index(orbital_name) ];
}


/// @details set the atom which connects to the lower connection
void
ResidueType::set_lower_connect_atom( std::string const & atm_name )
{
	finalized_ = false;
	if ( atm_name == "NONE" ) {
		if ( lower_connect_id_ != 0 ) {
			tr.Debug << "ERASING LOWER_CONNECT: " << lower_connect_id_ << " lcid: " << upper_connect_id_ << std::endl;
			utility::vector1< ResidueConnection >::iterator to_erase( residue_connections_.begin() );
			to_erase += lower_connect_id_ - 1;
			residue_connections_.erase(  to_erase );
			update_residue_connection_mapping();
			assert( n_polymeric_residue_connections_ != 0 );
			--n_polymeric_residue_connections_;
			if ( lower_connect_id_ < upper_connect_id_ ) { --upper_connect_id_; }
			lower_connect_id_ = 0;

		}
	} else {
		if ( lower_connect_id_ == 0 ) {
			ResidueConnection rc( atom_index( atm_name ) );
			residue_connections_.push_back( rc );
			lower_connect_id_ = residue_connections_.size();
			++n_polymeric_residue_connections_;
		} else {
			residue_connections_[ lower_connect_id_ ].atomno( atom_index( atm_name ) );
		}
	}
	update_residue_connection_mapping();
}

/// set the atom which connects to the upper connection
void
ResidueType::set_upper_connect_atom( std::string const & atm_name )
{
	finalized_ = false;
	if ( atm_name == "NONE" ) {
		if ( upper_connect_id_ != 0 ) {
			tr.Debug << "ERASING UPPER_CONNECT: " << upper_connect_id_ << " lcid: " << lower_connect_id_  << std::endl;
			utility::vector1< ResidueConnection >::iterator to_erase( residue_connections_.begin() );
			to_erase += upper_connect_id_ - 1;
			residue_connections_.erase(  to_erase );
			assert( n_polymeric_residue_connections_ != 0 );
			--n_polymeric_residue_connections_;
			if ( upper_connect_id_ < lower_connect_id_ ) { --lower_connect_id_; }
			upper_connect_id_ = 0;
		}
	} else {
		if ( upper_connect_id_ == 0 ) {
			ResidueConnection rc( atom_index( atm_name ) );
			residue_connections_.push_back( rc );
			upper_connect_id_ = residue_connections_.size();
			++n_polymeric_residue_connections_;
		} else {
			residue_connections_[ upper_connect_id_ ].atomno( atom_index( atm_name ) );
		}
	}
	update_residue_connection_mapping();
}

ResidueConnection const &
ResidueType::upper_connect() const
{
	assert( is_polymer_ );
	assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ];
}

///
ResidueConnection const &
ResidueType::lower_connect() const
{
	assert( is_polymer_ );
	assert( lower_connect_id_ != 0 );
	return residue_connections_[ lower_connect_id_ ];
}

Size
ResidueType::lower_connect_atom() const {
	assert( is_polymer_ );
	assert( lower_connect_id_ != 0 );
	return residue_connections_[ lower_connect_id_ ].atomno();
}

/// @brief index number of the atom which connects to the upper connection
Size
ResidueType::upper_connect_atom() const
{
	assert( is_polymer_ );
	assert( upper_connect_id_ != 0 );
	return residue_connections_[ upper_connect_id_ ].atomno();
}

/// @brief number of ResidueConnections, counting polymeric residue connections
Size
ResidueType::n_residue_connections() const
{
	return residue_connections_.size();
}

Size
ResidueType::residue_connect_atom_index( Size const resconn_id ) const {
	return residue_connections_[ resconn_id ].atomno();
}


/// @brief get a ResidueConection
ResidueConnection const &
ResidueType::residue_connection( Size const i ) const
{
	return residue_connections_[i];
}

ResidueConnection &
ResidueType::residue_connection( Size const i )
{
	return residue_connections_[ i ];
}

///////////////////////////////////////////////////////////////////////////////
//

/// @brief Get the chemical atom_type for this atom by it index number in this residue
///
/// @details If we want the atom_type index (integer), we get this from
/// the conformation::Atom itself, as seen in the code below
AtomType const &
ResidueType::atom_type( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <= ordered_atoms_.size()), "ResidueType::atom_type( Size const atomno ): atomno is not in this ResidueType!");
	assert( (atomno > 0) && (atomno <= ordered_atoms_.size()) );
	return ( *atom_types_ )[ graph_[ ordered_atoms_[atomno] ].atom_type_index() ];
}

/// @brief Get the atom name by index
std::string const &
ResidueType::atom_name( Size const index ) const
{
	PyAssert((index > 0) && (index <= ordered_atoms_.size()), "ResidueType::atom_name( Size const index ): index is not in this ResidueType!");
	assert((index > 0) && (index <= ordered_atoms_.size()));
	return graph_[ ordered_atoms_[index] ].name();
}

/// @brief get index of an atom's base atom
Size
ResidueType::atom_base( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <= atom_base_.size()), "ResidueType::atom_base( Size const atomno ): atomno is not in this ResidueType!");
	return atom_base_[ atomno ];
}

/// @brief get index of an atom's second base atom
Size
ResidueType::abase2( Size const atomno ) const
{
	PyAssert((atomno > 0) && (atomno <= abase2_.size()), "ResidueType::abase2( Size const atomno ): atomno is not in this ResidueType!");
	return abase2_[ atomno ];
}

Size
ResidueType::number_bonded_heavyatoms( Size const atomno ) const
{
	return bonded_neighbor_[ atomno ].size () - number_bonded_hydrogens( atomno );
	/// Graph style
	Size count = 0;
	VD vd = ordered_atoms_[atomno];
	for(OutEdgeIterPair ep = boost::out_edges(vd, graph_); ep.first != ep.second; ++ep.first){
		OutEdgeIter e_iter= ep.first;
		ED ed = *e_iter;
		VD target = boost::target(ed, graph_);
		AtomType const& at = (*atom_types_)[graph_[target].atom_type_index()];
		if( at.is_heavyatom() ) ++count;
	}
	return count;
}

/// @brief indices of the bonded neighbors for an atom
AtomIndices const &
ResidueType::bonded_neighbor( Size const atomno ) const
{
	//return bonded_neighbor_[ atomno ];
	/// Graph style
	VD vd = ordered_atoms_[atomno];
	AtomIndices atoms;
	for(OutEdgeIterPair ep = boost::out_edges(vd, graph_); ep.first != ep.second; ++ep.first){
		OutEdgeIter e_iter= ep.first;
		ED ed = *e_iter;
		VD target = boost::target(ed, graph_);
		for(Size i=1; i <= ordered_atoms_.size(); ++i){
			if( ordered_atoms_[i] == target ) atoms.push_back( i );
		}
	}

	assert(bonded_neighbor_[atomno] == atoms);
	return bonded_neighbor_[ atomno ];
}

utility::vector1<BondName> const &
ResidueType::bonded_neighbor_types(Size const atomno) const
{
	//return bonded_neighbor_type_[atomno];
	/// Graph style
	utility::vector1<BondName> bond_names;
	VD vd = ordered_atoms_[atomno];
	for(OutEdgeIterPair ep = boost::out_edges(vd, graph_); ep.first != ep.second; ++ep.first){
		OutEdgeIter e_iter= ep.first;
		ResidueGraph::edge_descriptor ed = *e_iter;
		Bond b = graph_[ed];
		bond_names.push_back( b.bond_name() );
	}
	assert( bond_names == bonded_neighbor_type_[atomno]);
	return bonded_neighbor_type_[atomno];
}


/// @note this does not set xyz coordinates for the added atom
void
ResidueType::add_atom(
	std::string const & atom_name,
	std::string const & atom_type_name,
	std::string const & mm_atom_type_name,
	Real const charge
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	assert(atom_graph_index_.find(atom_name) == atom_graph_index_.end());
	assert(atom_graph_index_.find( strip_whitespace(atom_name)) == atom_graph_index_.end());

	// increment atom count
	++natoms_;
	// index lookup by name
	// store the atom types
	// the next calls will fail if the atom type name is unrecognized
	Size const type( atom_types_->atom_type_index( atom_type_name ) );
	Size const mm_type( mm_atom_types_->atom_type_index( mm_atom_type_name ) );

	Atom atom(
				atom_name,
				mm_atom_type_name,
				type, mm_type,
				charge,
				Vector(0.0)
				//,AtomIcoor(...) // Cannot set Icoor here because we have to first add the vertex to the ordered_atoms_ and atom_graph_index_
	);

	VD v = graph_.add_vertex( atom );
	atom_graph_index_[ atom_name ] = v;
	atom_graph_index_[ strip_whitespace( atom_name ) ] = v;
	ordered_atoms_.push_back(v);
	AtomAP graph_atom = &graph_[v];

	graph_atom->icoor( AtomICoor( 0.0, 0.0, 0.0, atom_name, atom_name, atom_name, *this ) );
	// store the name

	assert( ordered_atoms_.size() == natoms_ );

	if ( (*atom_types_)[type].is_acceptor() )
			++n_hbond_acceptors_;
	if ( (*atom_types_)[type].is_donor() )
			++n_hbond_donors_;

	// add the atomic weight of this atom
	if( elements_ )	{ // Be robust if elements_ isn't defined.
		std::string const & element_name= (*atom_types_)[type].element();
		int const element_index = elements_->element_index(element_name);
		molecular_mass_ += (*elements_)[element_index].mass();
		molar_mass_ += (*elements_)[element_index].weight();
	} else {
		tr.Warning << "WARNING Elements set undefined." << std::endl;
	}

	parents_.push_back(0);

	// allocate space for the new atom !!!!!!!!!!!!!!!!!!!!!!1
	// eg, in the atom/resconn map
	assert( atom_2_residue_connection_map_.size() == natoms_-1 &&
					bonded_neighbor_.size() == natoms_-1 &&
					atom_base_.size() == natoms_-1);
	atom_2_residue_connection_map_.resize( natoms_ );
	bonded_neighbor_.resize( natoms_ ); // added here to handle residues w/ no bonds
	bonded_neighbor_type_.resize( natoms_);
	cut_bond_neighbor_.resize( natoms_ ); // added here to handle residues w/ no bonds
	atom_base_.push_back( natoms_ ); // default value, should generally be changed

	orbital_bonded_neighbor_.resize(natoms_);

}

/// @brief flag an atom for deletion by adding its index to the delete_atom_ list
void
ResidueType::delete_atom( std::string const & name )
{
	assert( has( name ) );
	delete_atom( atom_index(name) );
}

/// @brief flag an atom for deletion by adding its index to the delete_atom_ list
void
ResidueType::delete_atom( Size const index )
{
	finalized_ = false;
	delete_atoms_.push_back( index );

	VD const vd = ordered_atoms_[index];
	graph_.clear_vertex(vd);
	graph_.remove_vertex(vd);
}

/// @brief set atom type
void
ResidueType::set_atom_type(
	std::string const & atom_name,
	std::string const & atom_type_name
)
{
	int const atom_type_index = atom_types_->atom_type_index( atom_type_name );
	Atom & a = graph_[ atom_graph_index_[atom_name] ];
	a.atom_type_index( atom_type_index );
}


/// @brief set mm atom type
void
ResidueType::set_mm_atom_type(
	std::string const & atom_name,
	std::string const & mm_atom_type_name
)
{
	int const mm_type_index = mm_atom_types_->atom_type_index( mm_atom_type_name );
	Atom & a = graph_[ atom_graph_index_[atom_name] ];
	a.mm_atom_type_index( mm_type_index );
}

/// @brief Get the MM atom_type for this atom by its index number in this residue
MMAtomType const &
ResidueType::mm_atom_type( Size const atomno ) const
{
	return ( *mm_atom_types_ )[graph_[ ordered_atoms_[atomno] ].mm_atom_type_index() ];
}

orbitals::OrbitalType const &
ResidueType::orbital_type(int const orbital_index)const
{
	return ( *orbital_types_ )[ orbitals_[ orbital_index ].orbital_type_index() ];
}

//core::Size
//ResidueType::orbital_type_index(Size const orb_index) const
//{
//	return orbitals_[orb_index].orbital_type_index();
//}

/// @note this does not set xyz coordiates for the added orbital but sets the index of the orbital and maps
/// it to the type of orbital.
void
ResidueType::add_orbital(
	std::string & orbital_name,
	std::string & orbital_type_name
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	// increment orbital count
	++n_orbitals_;

	// store the atom type
	// the next call will fail if the orbital type name is unrecognized
	Size type( orbital_types_->orbital_type_index( orbital_type_name ) );

	// store the name
	orbitals_.push_back(Orbital(orbital_name, type, Vector(0.0)));
	assert( orbitals_.size() == n_orbitals_ );

	orbitals_index_[ orbital_name ] = n_orbitals_;
	orbitals_index_[ strip_whitespace( orbital_name ) ] = n_orbitals_;
}

///////////////////////////////////////////////////////////////////////////////

/// @details add a bond between atom1 and atom2 and add a BondType object referencing the bond (default bond type of SingleBond)
/** update bonded_neighbor_ and resize it as necessary **/
void
ResidueType::add_bond(
	std::string const & atom_name1,
	std::string const & atom_name2
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
	   std::string message = "add_bond: atoms " + atom_name1 + " and " + atom_name2 + " dont exist!";
		utility_exit_with_message( message  );
	}

	/////// Standard Version /////////

	Size const i1( atom_index( atom_name1 ));
	Size const i2( atom_index( atom_name2 ) );

	if ( bonded_neighbor_.size() < Size( std::max( i1, i2) ) ) {
		// ensure enough space for nbrs
		utility_exit_with_message("ResidueType:: shouldnt get here -- resizing in add_atom");
		//bonded_neighbor_.resize( std::max( i1,i2 ) );
	} else {
		// check if bond already exists
		AtomIndices const & i1_nbrs( bonded_neighbor_[i1] );
		if ( std::find( i1_nbrs.begin(), i1_nbrs.end(), i2 ) != i1_nbrs.end() ) {
			utility_exit_with_message( "dont add residue bonds more than once!" );
		}
	}

	bonded_neighbor_[ i1 ].push_back( i2 );
	bonded_neighbor_[ i2 ].push_back( i1 );
	bonded_neighbor_type_[i1].push_back(SingleBond);
	bonded_neighbor_type_[i2].push_back(SingleBond);


	////// GRAPH VERSION //////////

	NameVDMap::const_iterator source = atom_graph_index_.find( atom_name1 );
	NameVDMap::const_iterator target = atom_graph_index_.find( atom_name2 );

	assert( source != atom_graph_index_.end());
	assert( target != atom_graph_index_.end());

	VD const vd_source= source->second;
	VD const vd_target = target->second;

	Atom a = graph_[vd_source];
	Atom b = graph_[vd_target];

	// check if bond already exists...
	if( boost::edge(vd_source, vd_target, graph_).second ){
		utility_exit_with_message( "dont add residue bonds more than once!" );
	}

	ResidueGraph::edge_descriptor e_added;
	bool added;
	boost::tie(e_added, added) = graph_.add_edge( vd_source, vd_target);
	assert(added);
}

/// @details add a bond between atom1 and atom2 and add a BondType object referencing the bond using the specified bondName
void ResidueType::add_bond(std::string const & atom_name1, std::string const & atom_name2, BondName bondLabel)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
		std::string message = "add_bond: atoms " + atom_name1 + " and " + atom_name2 + " dont exist!";
		utility_exit_with_message( message  );
	}

	/////// Standard Version /////////

	Size const i1( atom_index( atom_name1 ) );
	Size const i2( atom_index( atom_name2 ) );

	if ( bonded_neighbor_.size() < Size( std::max( i1, i2) ) ) {
		// ensure enough space for nbrs
		utility_exit_with_message("ResidueType:: shouldnt get here -- resizing in add_atom");
		//bonded_neighbor_.resize( std::max( i1,i2 ) );
	} else {
		// check if bond already exists
		AtomIndices const & i1_nbrs( bonded_neighbor_[i1] );
		if ( std::find( i1_nbrs.begin(), i1_nbrs.end(), i2 ) != i1_nbrs.end() ) {
			utility_exit_with_message( "dont add residue bonds more than once!" );
		}
	}

	bonded_neighbor_[ i1 ].push_back( i2 );
	bonded_neighbor_[ i2 ].push_back( i1 );
	bonded_neighbor_type_[i1].push_back(bondLabel);
	bonded_neighbor_type_[i2].push_back(bondLabel);
	//bondType_vector_.push_back(BondType(i1,i2,bondLabel));

	NameVDMap::const_iterator source = atom_graph_index_.find( atom_name1 );
	NameVDMap::const_iterator target = atom_graph_index_.find( atom_name2 );
	assert( source != atom_graph_index_.end());
	assert( target != atom_graph_index_.end());
	VD const vd_source = source->second;
	VD const vd_target = target->second;

	Atom a = graph_[vd_source];
	Atom b = graph_[vd_target];

	// check if bond already exists...
	if( boost::edge(vd_source, vd_target, graph_).second ){
		utility_exit_with_message( "dont add residue bonds more than once!" );
	}

	ResidueGraph::edge_descriptor e_added;
	bool added;
	boost::tie(e_added, added) = graph_.add_edge( vd_source, vd_target);
	assert(added);
}

///////////////////////////////////////////////////////////////////////////////
///@brief add an orbital bond between an atom and an orbital.
///@note NOTE!!!!! This is indexed based upon atoms, not orbitals. That means that in your params file
/// you must have the atom as the first and orbital as the second.
void
ResidueType::add_orbital_bond(
	std::string const & atom_name1,
	std::string const & orbital_name
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has_orbital( orbital_name ) ) {
	   std::string message = "add_bond: atoms " + atom_name1 + " and " + orbital_name + " dont exist!";
		utility_exit_with_message( message  );
	}

	Size const i2( orbitals_index_[ orbital_name ] );

	if(atom_graph_index_.find(atom_name1) == atom_graph_index_.end()){
		utility_exit_with_message("atom_name: " + atom_name1 +" not found. Improper params file!");

	}

	Size const i1( atom_index( atom_name1 ) );

	orbital_bonded_neighbor_[i1].push_back(i2);

}

orbitals::ICoorOrbitalData const &
ResidueType::orbital_icoor_data(Size const orbital_index) const{
	return orbitals_[orbital_index].icoor();
}


orbitals::ICoorOrbitalData const &
ResidueType::new_orbital_icoor_data(Size const orbital_index) const{
	return orbitals_[orbital_index].new_icoor();
}


/// @details add a cut_bond between atom1 and atom2, which disallows an atom-tree connection,
///            though the atoms are really bonded.
/** update cut_bond_ and resize it as necessary **/
void
ResidueType::add_cut_bond(
	std::string const & atom_name1,
	std::string const & atom_name2
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
	   std::string message = "add_cut_bond: atoms " + atom_name1 + " and " + atom_name2 + " dont exist!";
		utility_exit_with_message( message  );
	}

	Size const i1( atom_index( atom_name1 ) );
	Size const i2( atom_index( atom_name2 ) );

	if ( cut_bond_neighbor_.size() < Size( std::max( i1, i2) ) ) {
		// ensure enough space for nbrs
		utility_exit_with_message("ResidueType:: shouldnt get here -- resizing in add_atom");
		//bonded_neighbor_.resize( std::max( i1,i2 ) );
	} else {
		// check if bond already exists
		AtomIndices const & i1_nbrs( cut_bond_neighbor_[i1] );
		if ( std::find( i1_nbrs.begin(), i1_nbrs.end(), i2 ) != i1_nbrs.end() ) {
			utility_exit_with_message( "dont add residue bonds more than once!" );
		}
	}

	cut_bond_neighbor_[ i1 ].push_back( i2 );
	cut_bond_neighbor_[ i2 ].push_back( i1 );
}


///////////////////////////////////////////////////////////////////////////////

/// add a chi angle defined by four atoms
void
ResidueType::add_chi(
	Size const chino,
	std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4
	)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ||
			 !has( atom_name3 ) || !has( atom_name4 ) ) {
		utility_exit_with_message("ResidueType::add_chi: atoms dont exist!" );
	}

	AtomIndices atoms;
	atoms.push_back( atom_index( atom_name1 ) );
	atoms.push_back( atom_index( atom_name2 ) );
	atoms.push_back( atom_index( atom_name3 ) );
	atoms.push_back( atom_index( atom_name4 ) );
	if ( chi_atoms_.size() < chino ) {
		chi_atoms_.resize( chino );
		chi_rotamers_.resize( chino );
		chi_2_proton_chi_.resize( chino );
	}
	chi_atoms_[chino] = atoms;

	is_proton_chi_.push_back( false );
	chi_2_proton_chi_[ chino ] = 0;
} // add_chi


/// @details Describe proton behavior for residue type; where should rotamer samples be considered,
/// and if expanded rotamers are desired, what deviations from the original rotamer samples
/// should be included.
/// E.g. dihedral_samples of 60, -60, and 180 could have an extra_sample of
/// 20 which would produce rotamers at 40 60 & 80, -40 -60 & -80, and -160, 180 & 160.
/// Extra_samples at 10 and 20 would produce 15 different rotamer samples.
void
ResidueType::set_proton_chi(
	Size chino,
	utility::vector1< Real > dihedral_samples,
	utility::vector1< Real > extra_samples
)
{
	is_proton_chi_[ chino ] = true;
	proton_chis_.push_back( chino );
	proton_chi_samples_.push_back( dihedral_samples );
	proton_chi_extra_samples_.push_back( extra_samples );
	chi_2_proton_chi_[ chino ] = proton_chis_.size();
}

///////////////////////////////////////////////////////////////////////////////

/// @details add a rotamer bin for a given chi
/** a rotamer bin has the mean and standard deviation**/
void
ResidueType::add_chi_rotamer(
	Size const chino,
	Real const mean,
	Real const sdev
)
{
	if ( chi_rotamers_.size() < chino ) chi_rotamers_.resize( chino );
	chi_rotamers_[chino].push_back( std::make_pair( mean, sdev ) );
}


///////////////////////////////////////////////////////////////////////////////

/// @details sets atom_base_[atom1] = atom2
/** resize atom_base_ vector as necessary **/
void
ResidueType::set_atom_base(
	std::string const & atom_name1,
	std::string const & atom_name2
)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ) {
		utility_exit_with_message( "set_atom_base: atoms dont exist!" );
	}


	Size const i1( atom_index( atom_name1 ) );
	Size const i2( atom_index( atom_name2 ) );

	// debug connectivity
	AtomIndices const & i1_nbrs( bonded_neighbor_[i1] );
	AtomIndices const & i1_cut_nbrs( cut_bond_neighbor_[i1] );

	if ( ( std::find( i1_nbrs.begin(), i1_nbrs.end(), i2 ) == i1_nbrs.end() ) &&
			 !( i1 == 1 && i2 == 1 && natoms_ == 1 ) ) { // note we allow special exception for single-atom residue
		utility_exit_with_message( "set_atom_base: atoms must be bonded!" );
	}
	if ( atom_base_.size() < Size(i1) ) utility_exit_with_message("ResidueType:: shouldnt get here!");

	//Need to add a cut-bond check too!!!
	if ( std::find( i1_cut_nbrs.begin(), i1_cut_nbrs.end(), i2 ) != i1_cut_nbrs.end() )  {
		utility_exit_with_message( "set_atom_base: cut_bond disallows specification of atom base!" );
	}

	atom_base_[ i1 ] = i2;

}

///////////////////////////////////////////////////////////////////////////////

/// @details get all specified properties for this residue type
utility::vector1< std::string > const &
ResidueType::properties() const
{
	return properties_;
}

///////////////////////////////////////////////////////////////////////////////

/// @details add a property to this residue
/** update boolean property member data accordingly **/
void
ResidueType::add_property( std::string const & property )
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( property == "POLYMER" ) {
		is_polymer_ = true;
	} else if ( property == "PROTEIN" ) {
		is_protein_ = true;
		is_polymer_ = true;
	} else if ( property == "POLAR" ) {
		is_polar_ = true;
	} else if( property == "SC_ORBITALS"){
		has_sc_orbitals_ = true;
	} else if ( property == "CHARGED" ) {
		is_charged_ = true;
	} else if ( property == "AROMATIC" ) {
		is_aromatic_ = true;
	} else if ( property == "COARSE" ) {
		is_coarse_ = true; //currently only for RNA
	} else if ( property == "DNA" ) {
		is_DNA_ = true;
		is_NA_ = true;
		is_polymer_ = true;
	} else if ( property == "RNA" ) {
		is_RNA_ = true;
		is_NA_ = true;
		is_polymer_ = true;
	} else if ( property == "CARBOHYDRATE") {
		is_carbohydrate_ = true;
	} else if ( property == "LIGAND" ) {
		is_ligand_ = true;
	} else if ( property == "SURFACE" ) {
		is_surface_ = true;
	} else if ( property == "LOWER_TERMINUS" ) {
		is_terminus_ = true;
		is_lower_terminus_ = true;
	} else if ( property == "UPPER_TERMINUS" ) {
		is_terminus_ = true;
		is_upper_terminus_ = true;
	} else if ( property == "LOWERTERM_TRUNC" ) {
		is_terminus_ = true;
		is_lower_terminus_ = true;
	} else if ( property == "UPPERTERM_TRUNC" ) {
		is_terminus_ = true;
		is_upper_terminus_ = true;
	} else if ( property == "PHOSPHONATE" ) {
		is_polymer_ = true;
		is_phosphonate_ = true;
	} else if ( property == "PHOSPHONATE_UPPER" ) {
		is_terminus_ = true;
		is_upper_terminus_ = true;
		is_phosphonate_ = true;
		is_phosphonate_upper_ = true;
	} else if ( property == "TERMINUS" ) {
		is_terminus_ = true;
	} else if ( property == "ACETYLATED_NTERMINUS" ) {
		is_terminus_ = true;
		is_lower_terminus_ = true;
		is_acetylated_nterminus_ = true;
	} else if ( property == "METHYLATED_CTERMINUS" ) {
		is_terminus_ = true;
		is_upper_terminus_ = true;
		is_methylated_cterminus_ = true;
	} else if (property == "BRANCH_POINT") {
		;  // Null statement for now.... ~ Labonte
    } else if (
            (property == "ALDOSE") ||
            (property == "KETOSE") ||
            (property == "L_SUGAR") ||
            (property == "D_SUGAR") ||
            (property == "FURANOSE") ||
            (property == "PYRANOSE") ||
            (property == "SEPTANOSE") ||
            (property == "ALPHA_SUGAR") ||
            (property == "BETA_SUGAR") ||
            (property == "URONIC_ACID")) {
        ;  // Null statement -- these properties will be added to carbohydrate_info_ by update_derived_data().
	} else {
		tr.Warning << "WARNING:: unrecognized residue type property: " << property << std::endl;
	}

	properties_.push_back( property );
}

///////////////////////////////////////////////////////////////////////////////

/// @details delete a property to this residue
/** update boolean property member data accordingly **/
//    Added by Andy M. Chen in June 2009
//    This is needed for deleting properties, which occurs in certain PTM's

void
ResidueType::delete_property( std::string const & property )
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( property == "POLYMER" ) {
		is_polymer_ = false;
	} else if ( property == "PROTEIN" ) {
		is_protein_ = false;
	} else if ( property == "POLAR" ) {
		is_polar_ = false;
	}else if(property == "SC_ORBITALS"){
		has_sc_orbitals_ = false;
	}else if ( property == "CHARGED" ) {
		is_charged_ = false;
	} else if ( property == "AROMATIC" ) {
		is_aromatic_ = false;
	} else if ( property == "COARSE" ) {
		is_coarse_ = false;
	} else if ( property == "DNA" ) {
		is_DNA_ = false;
	} else if ( property == "RNA" ) {
		is_RNA_ = false;
	} else if ( property == "CARBOHYDRATE") {
		is_carbohydrate_ = false;
	} else if ( property == "LIGAND" ) {
		is_ligand_ = false;
	} else if ( property == "SURFACE" ) {
		is_surface_ = false;
	} else if ( property == "LOWER_TERMINUS" ) {
		// could add an is_lower_terminus_ bool if needed?
		is_lower_terminus_ = false;
	} else if ( property == "UPPER_TERMINUS" ) {
		is_upper_terminus_ = false;
	} else if ( property == "TERMINUS" ) {
		is_terminus_ = false;
	} else if ( property == "PHOSPHONATE" ) {
		is_phosphonate_ = false;
	} else if ( property == "PHOSPHONATE_UPPER" ) {
		is_phosphonate_upper_ = false;
	} else if ( property == "ACETYLATED_NTERMINUS" ) {
		is_acetylated_nterminus_ = false;
	} else if ( property == "METHYLATED_CTERMINUS" ) {
		is_methylated_cterminus_ = false;
	} else {
		tr.Warning << "WARNING:: unrecognized residuetype property: " << property << std::endl;
	}

	properties_.push_back( property );
}


///////////////////////////////////////////////////////////////////////////////

/// @details redefine a chi angle based on four atoms
//    Added by Andy M. Chen in June 2009
//    This function is almost an exact copy of the add_chi function except that vector resizing does NOT occur.
//    It is needed for certain PTM's that affects proton chis (e.g. phosphorylation and sulfation).

void
ResidueType::redefine_chi(
	Size const chino,
	std::string const & atom_name1,
	std::string const & atom_name2,
	std::string const & atom_name3,
	std::string const & atom_name4
	)
{
	// signal that we need to update the derived data
	finalized_ = false;

	if ( !has( atom_name1 ) || !has( atom_name2 ) ||
			 !has( atom_name3 ) || !has( atom_name4 ) ) {
		utility_exit_with_message("ResidueType::redefine_chi: atoms dont exist!" );
	}

	AtomIndices atoms;
	atoms.push_back( atom_index( atom_name1 ) );
	atoms.push_back( atom_index( atom_name2 ) );
	atoms.push_back( atom_index( atom_name3 ) );
	atoms.push_back( atom_index( atom_name4 ) );
	chi_atoms_[chino] = atoms;

	//Assumes that the redifined chi is NOT a proton chi.
	//  (This is adequate in most cases because PTMs tend to replace hydrogens
	//  with functional groups rather than the other way around.)
	is_proton_chi_[ chino ] = false;
	chi_2_proton_chi_[ chino ] = 0;


} // redefine_chi



/////////////////////////////////////////////////////////////////

/// @details add an atom to the list for calculating actcoord center
void
ResidueType::add_actcoord_atom( std::string const & atom )
{
	assert( is_protein() );
	finalized_ = false;
	tr.Trace << "adding act coord atom: " << name_ << ' ' << atom << std::endl;
	Size atomindex = atom_index( atom );
	actcoord_atoms_.push_back( atomindex );
	++n_actcoord_atoms_;
}

///////////////////////////////////////////////////////////////////////////////

///@details set up atom ordering map old2new, called by finalize()

/**
	 because some new heavy atoms are added by patching, some are flagged to be deleted
	 in delete_atoms_ and some are forced to be backbon atoms as in force_bb_

	 old2new[old_atom_index] = new_atom_index

	 sets natoms_, nheavyatoms_, and n_backbone_heavyatoms_

	 also fills attached_H_begin, attached_H_end

**/

void
ResidueType::setup_atom_ordering(
	AtomIndices & old2new
)
{
	/////////////////////////////////////////////////////////////////////////////
	// reorder!
	//
	// preserve order of heavyatoms in current list, modulo new backbone heavyatoms that have been added by patching
	// these guys need to be moved forward
	//
	// ** ensure that heavyatoms come before hydrogens
	// ** put all hydrogens attached to a heavyatom in a single row
	//
	// ** adding support for deleting atoms at this stage: calls to delete_atom( name )
	//    will append to delete_atoms_, must follow with call to finalize()


	// set atom counts

	// Graph Atoms are already deleted! Now all we need to do is sort them...

	Size const old_natoms( ordered_atoms_.size() );
	assert( natoms_ == old_natoms );
	assert(graph_.num_vertices() == old_natoms - delete_atoms_.size());
	natoms_ = old_natoms - delete_atoms_.size();

	// size and initialize the mapping from old to new indices
	old2new.clear();
	old2new.resize( old_natoms, 0 );

	// process delete_atoms_ to a friendlier form
	utility::vector1< bool > keep_me( old_natoms, true );
	for ( Size i=1; i<= old_natoms; ++i ) {
		if ( std::find( delete_atoms_.begin(), delete_atoms_.end(), i ) != delete_atoms_.end() ) {
			keep_me[i] = false;
		}
	}
	delete_atoms_.clear();

	// first fill a list of all the heavyatoms, insisting that backbone come first
	AtomIndices old_heavyatom_indices;
	for ( Size i=1; i<= old_natoms; ++i ) {
		if ( keep_me[ i ] && atom_type( i ).is_heavyatom() &&
				 ( ( i <= n_backbone_heavyatoms_ ) || // is backbone heavyatom by count from previous call to finalize()
					 ( std::find( force_bb_.begin(), force_bb_.end(), i ) != force_bb_.end() ) ) ) { // is forced-bb
			old_heavyatom_indices.push_back( i );
		}
	}

	// update the bb-heavy counter:
	n_backbone_heavyatoms_ = old_heavyatom_indices.size();
	force_bb_.clear();

	// now add sidechain heavyatoms
	for ( Size i=1; i<= old_natoms; ++i ) {
		if ( keep_me[ i ] && atom_type( i ).is_heavyatom() &&
				 ( std::find( old_heavyatom_indices.begin(), old_heavyatom_indices.end(), i ) == // not already done
					 old_heavyatom_indices.end() ) ) {
			old_heavyatom_indices.push_back( i );
		}
	}

	// all the heavyatoms
	nheavyatoms_ = old_heavyatom_indices.size();

	// setup old2new, also fill in attached_H_begin/end
	attached_H_begin_.clear();
	attached_H_end_.clear();
	attached_H_begin_.resize( nheavyatoms_, 0 );
	attached_H_end_.resize( nheavyatoms_, 0 /*do not change*/);

	Size new_H_index( nheavyatoms_ );
	for ( Size new_heavy_index = 1; new_heavy_index <= nheavyatoms_; ++new_heavy_index ) {
		Size const old_heavy_index( old_heavyatom_indices[ new_heavy_index ] );
		// mapping between indices
		old2new[ old_heavy_index ] = new_heavy_index;

		// now add the attached hydrogens
		attached_H_begin_[ new_heavy_index ] = new_H_index + 1;
		AtomIndices const & nbrs( bonded_neighbor_[ old_heavy_index ] );
		for ( Size j=1; j<= nbrs.size(); ++j ) {
			Size const old_H_index( nbrs[j] );
			if ( keep_me[ old_H_index ] && (*atom_types_)[ graph_[ordered_atoms_[ old_H_index ]].atom_type_index() ].is_hydrogen()) {
					++new_H_index;
					old2new[ old_H_index ] = new_H_index;
					attached_H_end_[ new_heavy_index ] = new_H_index;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////

/// reorder primary data in ResidueType given the old2new map, called by finalize()
/**
		update the rest private data in ResidueType object after old2new map is set up
		by calling setup_atom_ordering (some data have been updated there also)
**/
void
ResidueType::reorder_primary_data(
	AtomIndices const & old2new
)
{
	// now reorder using old2new  -- a bit of a nuisance
	// there must be a better way to handle this!!!
	Size const old_natoms( old2new.size() );
	assert( old_natoms == ordered_atoms_.size() );

	// copy the old per-atom data: note that attached_H_begin and attached_H_end have already been setup
	// and abase2_ is derived data setup down below
	VDs old_ordered_atoms( ordered_atoms_ );

	utility::vector1< utility::vector1 <core::Size> > old_orbital_bonded_neighbor(orbital_bonded_neighbor_);
	utility::vector1< AtomIndices > old_bonded_neighbor( bonded_neighbor_ );
	utility::vector1<utility::vector1<BondName> > old_bonded_neighbor_type(bonded_neighbor_type_);
	utility::vector1< AtomIndices > old_cut_bond_neighbor( cut_bond_neighbor_ );
	AtomIndices old_atom_base( atom_base_ );

	utility::vector1<Size> old_parents(parents_);

	if ( old_natoms != natoms_ ) { // because we deleted some atoms
		ordered_atoms_.resize( natoms_ );
		orbital_bonded_neighbor_.resize( natoms_ );
		bonded_neighbor_.resize( natoms_ );
		bonded_neighbor_type_.resize (natoms_);
		cut_bond_neighbor_.resize( natoms_ );
		atom_base_.resize( natoms_ );
		parents_.resize(natoms_);
		atom_2_residue_connection_map_.resize( natoms_ );
	}
	// fill in the new data
	for ( Size old_index=1; old_index<= old_natoms; ++old_index ) {
		Size const new_index( old2new[ old_index ] );
		if ( new_index == 0 ) continue; // deleted
		ordered_atoms_[ new_index] = old_ordered_atoms[ old_index];
		atom_base_[ new_index ] = old2new[ old_atom_base[ old_index ] ];
		assert( atom_base_[ new_index ] ); // this will fail if we deleted an atom which was the atom_base for another atom

		orbital_bonded_neighbor_[new_index] = old_orbital_bonded_neighbor[old_index];

		bonded_neighbor_[ new_index ].clear();
		for ( Size j=1, j_end = old_bonded_neighbor[ old_index ].size(); j<= j_end; ++j )
		{
			Size const old_nbr( old_bonded_neighbor[ old_index ][j] );
			if ( old2new[ old_nbr ] ) bonded_neighbor_[ new_index ].push_back( old2new[ old_nbr ] );
		}


		bonded_neighbor_type_[new_index].clear();
		for(Size j=1, j_end=old_bonded_neighbor_type[old_index].size(); j<=j_end; ++j)
		{
			Size const old_neighbor(old_bonded_neighbor[old_index][j]) ;
			if(old2new[old_neighbor] ) bonded_neighbor_type_[new_index].push_back(old_bonded_neighbor_type[old_index][j]);
		}

		cut_bond_neighbor_[ new_index ].clear();
		for ( Size j=1, j_end = old_cut_bond_neighbor[ old_index ].size(); j<= j_end; ++j )
		{
			Size const old_nbr( old_cut_bond_neighbor[ old_index ][j] );
			if ( old2new[ old_nbr ] ) cut_bond_neighbor_[ new_index ].push_back( old2new[ old_nbr ] );
		}

		parents_[ new_index ] = old_parents[ old_index ];

		VD vd = ordered_atoms_[ new_index ];
		Atom & a2 =  graph_[vd];
		for ( Size i=1; i<= 3; ++i ) {
			//ICoorAtomID & stub_atom( a1->icoor().stub_atom( i ) );
			ICoorAtomID & ordered_stub_atom( a2.icoor().stub_atom( i ) );
			//assert(stub_atom.type() == ordered_stub_atom.type());
			//assert(stub_atom.atomno()== ordered_stub_atom.atomno());
//			if ( stub_atom.type() == ICoorAtomID::INTERNAL ) {
//				stub_atom.atomno( old2new[ stub_atom.atomno() ] );
//				assert( stub_atom.atomno() ); // this will fail if we deleted a stub atom for some other atom
//			}

			if ( ordered_stub_atom.type() == ICoorAtomID::INTERNAL ) {
				ordered_stub_atom.atomno( old2new[ ordered_stub_atom.atomno() ] );
				assert( ordered_stub_atom.atomno() ); // this will fail if we deleted a stub atom for some other atom
			}
			//assert(*a1 == a2);
		}

		if(orbitals_.size() != 0){
				utility::vector1< core::Size > const orbital_indices(orbital_bonded_neighbor_[new_index]);
				for (
						utility::vector1< core::Size >::const_iterator
						orbital_index = orbital_indices.begin(),
						orbital_index_end = orbital_indices.end();
						orbital_index != orbital_index_end; ++orbital_index
				)
				{

					core::Size stub1( orbitals_[*orbital_index].new_icoor().get_stub1());
					core::Size stub2( orbitals_[*orbital_index].new_icoor().get_stub2());
					core::Size stub3( orbitals_[*orbital_index].new_icoor().get_stub3() );

					if ( stub1 == 0 || stub2 == 0 || stub3 == 0) {
						continue;
					}else{
						orbitals_[*orbital_index].new_icoor().replace_stub1( old2new[stub1]);
						orbitals_[*orbital_index].new_icoor().replace_stub2( old2new[stub2]);
						orbitals_[*orbital_index].new_icoor().replace_stub3( old2new[stub3]);
					}


				}

		}



	}

	atom_graph_index_.clear();
	for ( Size i=1; i<= natoms_; ++i ) {
		atom_graph_index_[ graph_[ ordered_atoms_[i]].name() ] = ordered_atoms_[i];
		atom_graph_index_[ strip_whitespace( graph_[ ordered_atoms_[i]].name()  ) ] = ordered_atoms_[i];
	}


	// chi_atoms_
	for ( Size i=1; i<= chi_atoms_.size(); ++i ) {
		for ( Size j=1; j<= 4; ++j ) {
			chi_atoms_[i][j] = old2new[ chi_atoms_[i][j] ];
		}
	}

	// mainchain_atoms_
	for ( Size i=1; i<= mainchain_atoms_.size(); ++i ) {
		mainchain_atoms_[i] = old2new[ mainchain_atoms_[i] ];
	}

	// actcoord_atoms_
	AtomIndices const old_actcoord_atoms = actcoord_atoms_;
	actcoord_atoms_.clear();
	//assert( n_actcoord_atoms_ == actcoord_atoms_.size() );
	for ( Size i=1; i<= old_actcoord_atoms.size(); ++i ) {
		Size new_index =  old2new[ old_actcoord_atoms[i]];
		if(new_index) actcoord_atoms_.push_back(new_index); // if atom hasn't been deleted
	}

	// nbr_atom_
	nbr_atom_ = nbr_atom_ ? old2new[ nbr_atom_ ] : nbr_atom_;


	/////////////////////////////////////////////////////////////////////////////
	// add additional reordering statements here for new data that you've added
	// to  ResidueType  that's sensitive to atom order

	///TODO fix this problem: deleting atoms can invalidate these residue_connections_
	for ( Size i=1; i<= residue_connections_.size(); ++i ) {
		residue_connections_[i].atomno( old2new[ residue_connections_[i].atomno() ] );
		AtomICoor new_icoor = residue_connections_[i].icoor();
		for ( Size j = 1; j <= 3; ++j ) {
			new_icoor.stub_atom( j ).atomno( old2new[ new_icoor.stub_atom( j ).atomno() ] );
		}
		residue_connections_[ i ].icoor( new_icoor );
		assert( residue_connections_[i].atomno() ); //this will fail if we deleted an atom involved in an inter-rsd connection
	}

	update_residue_connection_mapping();
}


///////////////////////////////////////////////////////////////////////////////

/// @details update derived data in ResidueType, called by finalize()
/**
		after primary data have been reordered, update derived data acoordingly,
		including\n, Hbond donor and acceptors, path_distance etc.
**/
void
ResidueType::update_derived_data()
{

	first_sidechain_hydrogen_ = natoms_ + 1;
	for ( Size i= n_backbone_heavyatoms_ + 1; i<= nheavyatoms_; ++i ) {
		if ( attached_H_begin_[i] <= attached_H_end_[i] ) {
			first_sidechain_hydrogen_ = attached_H_begin_[i];
			break;
		}
	}


	// compile atom-index lists of subsets of the atoms
	accpt_pos_.clear();
	accpt_pos_sc_.clear();
	Haro_index_.clear();
	Hpol_index_.clear();
	atoms_with_orb_index_.clear();
	Hpos_polar_.clear();
	Hpos_apolar_.clear();
	Hpos_polar_sc_.clear();
	all_bb_atoms_.clear();
	all_sc_atoms_.clear();



	for ( Size i=1; i<= natoms_; ++i ) {
		AtomType const & type( (*atom_types_)[ graph_[ ordered_atoms_[i]].atom_type_index() ] );
		//////////////////////////////////
		// info derived from the atom_type
		// hbond info
/*		if(type.is_lone_pair()){
			Orbitals_index_.push_back( i );
		}*/
		//get atoms with orbitals on it
		if(type.atom_has_orbital()){
			atoms_with_orb_index_.push_back(i);
		}

		//get aromatic hydrogens
		if(type.is_haro()){
			Haro_index_.push_back( i );

		}
		//get polar hydrogens
		if(type.is_polar_hydrogen()){
			Hpol_index_.push_back( i );
		}

		if ( type.is_acceptor() ) {
			accpt_pos_.push_back( i );
			if ( i > n_backbone_heavyatoms_ ) {
				accpt_pos_sc_.push_back( i );
			}

		}

		//if ( type.is_polar_hydrogen() &&   (std::abs( graph_[ordered_atoms_[ natoms_ ]].charge() ) > 1.0e-3) ) {
		if ( type.is_polar_hydrogen() &&   (!type.is_virtual() ) ) {
			Hpos_polar_.push_back( i );
			if ( i >= first_sidechain_hydrogen_ ) {
				Hpos_polar_sc_.push_back( i );
			}
		}

		if ( type.is_hydrogen() && ( !type.is_polar_hydrogen() ) ){
			Hpos_apolar_.push_back( i );
		}

		// Which atoms are backbone and which are sidechain; sometimes nice to just get
		// lists instead of iterating over the subranges.
		if ( type.is_hydrogen() ) {
			if ( i < first_sidechain_hydrogen_ ) {
				all_bb_atoms_.push_back( i );
			} else {
				all_sc_atoms_.push_back( i );
			}
		} else {
			if ( i <= n_backbone_heavyatoms_ ) {
				all_bb_atoms_.push_back( i );
			} else {
				all_sc_atoms_.push_back( i );
			}
		}

	}

	// donor heavy atoms, acceptor heavy atoms, donor hydrogen atoms setup.
	// Must be executed after Hpos_polar_ and accpt_pos_ have been updated.
	heavyatom_has_polar_hydrogens_.resize( natoms_ );
	heavyatom_is_an_acceptor_.resize( natoms_ );
	atom_is_polar_hydrogen_.resize( natoms_ );
	std::fill( heavyatom_has_polar_hydrogens_.begin(), heavyatom_has_polar_hydrogens_.end(), 0 );
	std::fill( heavyatom_is_an_acceptor_.begin(), heavyatom_is_an_acceptor_.end(), 0 );
	std::fill( atom_is_polar_hydrogen_.begin(), atom_is_polar_hydrogen_.end(), 0 );
	for ( Size ii = 1; ii <= Hpos_polar_.size(); ++ii ) {
		Size hind = Hpos_polar_[ ii ];
		atom_is_polar_hydrogen_[ ii ] = 1;
		Size base = atom_base_[ hind ];
		heavyatom_has_polar_hydrogens_[ base ] = 1;
	}
	for ( Size ii = 1; ii <= accpt_pos_.size(); ++ii ) {
		heavyatom_is_an_acceptor_[ accpt_pos_[ ii ] ] = 1;
	}


	// now setup abase2
	abase2_.clear();
	abase2_.resize( natoms_ );

	for ( Size ii=1, ii_end= accpt_pos_.size(); ii<= ii_end; ++ii ) {
		uint const i( accpt_pos_[ii] );
		uint const i_base( atom_base_[i] );
		assert( i_base != 0 );
		AtomIndices const & i_nbrs( bonded_neighbor_[i] );
		if ( i_nbrs.size() == 0 ) {
		        utility_exit_with_message( "failed to set abase2 for acceptor atom, it has no nbrs!" );
		} else if ( i_nbrs.size() == 1 ) {
			assert( i_nbrs[1] == i_base );
			abase2_[ i ] = atom_base_[ i_base ];
			//iwd  The first child of the root is root's atom_base.
			//iwd  But if that child has no children, it ends up as its own abase2.
			//iwd  So instead we use the second child of the parent,
			//iwd  which must exist if there are 3+ atoms in this tree.
			if( abase2_[ i ] == i ) {
				AtomIndices const & i_base_nbrs( bonded_neighbor_[i_base] );
				for(Size jj = 1, jj_end = i_base_nbrs.size(); jj <= jj_end; ++jj) {
					if(i_base_nbrs[ jj ] != i) {
						abase2_[ i ] = i_base_nbrs[ jj ];
						//TR << "Using " << atom_name(abase2_[ i ]) << " as abase2 for " << atom_name(i) << std::endl;
						break;
					}
				}
			}
			assert( abase2_[ i ] != i && abase2_[ i ] != i_base && abase2_[ i ] != 0 );
		} else if ( i_nbrs[1] == i_base ) {
			abase2_[ i ] = i_nbrs[2];
		} else {
			abase2_[ i ] = i_nbrs[1];
		}
	}


	// bond path distances
	FArray2D_int path_distances( get_residue_path_distances( *this ));
	path_distance_.resize( natoms_ );
	for ( Size ii = 1; ii <= natoms_; ++ii ) {
		path_distance_[ ii ].resize( natoms_ );
		for (Size jj = 1; jj <= natoms_; ++jj ) {
			path_distance_[ ii ][ jj ] = path_distances( ii, jj );
		}
	}

	// get dihedral angles
	dihedral_atom_sets_.clear();
	dihedrals_for_atom_.resize( natoms_ );
	for ( Size ii = 1; ii <= natoms_; ++ii ) dihedrals_for_atom_[ ii ].clear();

	// get for all pairs of atoms seperated by 1 bond
	for ( Size central_atom1 = 1; central_atom1 < natoms_; ++central_atom1 ) {
		for ( Size central_atom2 = central_atom1+1; central_atom2 <= natoms_; ++central_atom2 ) {
			if ( path_distance_[ central_atom1 ][ central_atom2 ] == 1 ) {

				// get all atoms seperated from central_atom1/2 by one bond that are not central_atom2/1
				utility::vector1< Size > ca1d1;
				utility::vector1< Size > ca2d1;

				// ca1
				for ( Size i = 1; i <= natoms_; ++i ) {
					if ( ( path_distance_[ central_atom1 ][ i ] == 1 ) && ( i != central_atom2 ) ) {
						ca1d1.push_back( i );
					}
				}
				// ca2
				for ( Size i = 1; i <= natoms_; ++i ) {
					if ( ( path_distance_[ central_atom2 ][ i ] == 1 ) && ( i != central_atom1 ) ) {
						ca2d1.push_back( i );
					}
				}

				// for each pair of dihedral angle start or end atoms create a dihedral angle using central atom
				for ( utility::vector1< Size >::iterator terminal_atom1 = ca1d1.begin();
						terminal_atom1 != ca1d1.end(); ++terminal_atom1 ) {
					for ( utility::vector1< Size >::iterator terminal_atom2 = ca2d1.begin();
							terminal_atom2 != ca2d1.end(); ++terminal_atom2 ) {
						dihedral_atom_set temp( *terminal_atom1, central_atom1, central_atom2, *terminal_atom2 );
						dihedral_atom_sets_.push_back( temp );
						Size const which_dihedral = dihedral_atom_sets_.size();
						dihedrals_for_atom_[ *terminal_atom1 ].push_back( which_dihedral );
						dihedrals_for_atom_[   central_atom1 ].push_back( which_dihedral );
						dihedrals_for_atom_[   central_atom2 ].push_back( which_dihedral );
						dihedrals_for_atom_[ *terminal_atom2 ].push_back( which_dihedral );
					}
				}

			}
		}
	}
	ndihe_ = dihedral_atom_sets_.size();



	// get bond angles
	bondangle_atom_sets_.clear();
	bondangles_for_atom_.resize( natoms_ );
	for ( Size ii = 1; ii <= natoms_; ++ii ) bondangles_for_atom_[ ii ].clear();

	// iterate over all atoms that could be a central atom
	for ( Size central_atom = 1; central_atom <= natoms_; ++central_atom ) {

		AtomIndices const & bonded_neighbors( bonded_neighbor_[central_atom] );
		Size const num_bonded_neighbors( bonded_neighbors.size() );

		// create all possible combinations of branching atoms
		for ( Size i = 1; i < num_bonded_neighbors; ++i ) {
			for ( Size j = i+1; j <= num_bonded_neighbors; ++j ) {
				bondangle_atom_set temp( bonded_neighbors[i], central_atom, bonded_neighbors[j] );
				bondangle_atom_sets_.push_back( temp );
				Size const which_angle = bondangle_atom_sets_.size();
				bondangles_for_atom_[ bonded_neighbors[i] ].push_back( which_angle );
				bondangles_for_atom_[        central_atom ].push_back( which_angle );
				bondangles_for_atom_[ bonded_neighbors[j] ].push_back( which_angle );
			}
		}
	}

	// Now for inter-residue connections, find the sets of atoms that are within one and within two bonds
	// of a residue connection point.  From these sets, all inter-residue bond angle and bond torsions may
	// be enumerated when evaluating residue pair energies.  Also compute the backwards mapping: a list for
	// each atom of the within-1-bond and within-2-bond sets that the atom is listed as being part of. These
	// lists are needed when evaluating atom derivatives wrt the bond dihedrals and angles.
	atoms_within_one_bond_of_a_residue_connection_.resize( residue_connections_.size() );
	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) atoms_within_one_bond_of_a_residue_connection_[ ii ].clear();

	within1bonds_sets_for_atom_.resize( natoms_ );
	for ( Size ii = 1; ii <= natoms_; ++ii ) within1bonds_sets_for_atom_[ ii ].clear();

	atoms_within_two_bonds_of_a_residue_connection_.resize( residue_connections_.size() );
	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) atoms_within_two_bonds_of_a_residue_connection_[ ii ].clear();

	within2bonds_sets_for_atom_.resize( natoms_ );
	for ( Size ii = 1; ii <= natoms_; ++ii ) within2bonds_sets_for_atom_[ ii ].clear();

	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) {
		Size const ii_resconn_atom = residue_connections_[ ii ].atomno();

		AtomIndices const & ii_bonded_neighbors( bonded_neighbor_[ ii_resconn_atom ] );
		Size const ii_num_bonded_neighbors( ii_bonded_neighbors.size() );

		for ( Size jj = 1; jj <= ii_num_bonded_neighbors; ++jj ) {
			Size const jj_atom = ii_bonded_neighbors[ jj ];

			// Record that ii_resconn_atom and jj_atom are within a single bond of residue connection ii.
			two_atom_set wi1( ii_resconn_atom, jj_atom );
			atoms_within_one_bond_of_a_residue_connection_[ ii ].push_back( wi1 );

			// For atoms ii_resconn_atom and jj_atom, mark residue connection ii as a
			// connection point the are within one bond of.
			Size const which_wi1 = atoms_within_one_bond_of_a_residue_connection_[ ii ].size();
			within1bonds_sets_for_atom_[ ii_resconn_atom ].push_back( std::make_pair( ii, which_wi1 ) );
			within1bonds_sets_for_atom_[ jj_atom ].push_back( std::make_pair( ii, which_wi1 ));

			AtomIndices const & jj_bonded_neighbors( bonded_neighbor_[ jj_atom ] );
			Size const jj_num_bonded_neighbors( jj_bonded_neighbors.size() );

			for ( Size kk = 1; kk <= jj_num_bonded_neighbors; ++kk ) {
				Size const kk_atom = jj_bonded_neighbors[ kk ];
				if ( kk_atom == ii_resconn_atom ) continue; // skip iiat->jjat->iiat

				three_atom_set wi2( ii_resconn_atom, jj_atom, kk_atom );
				atoms_within_two_bonds_of_a_residue_connection_[ ii ].push_back( wi2 );

				Size const which_wi2 = atoms_within_two_bonds_of_a_residue_connection_[ ii ].size();
				within2bonds_sets_for_atom_[ ii_resconn_atom ].push_back( std::make_pair( ii, which_wi2 ) );
				within2bonds_sets_for_atom_[ jj_atom ].push_back( std::make_pair( ii, which_wi2 ));
				within2bonds_sets_for_atom_[ kk_atom ].push_back( std::make_pair( ii, which_wi2 ));
			}
		}
	}


	if(is_RNA_){ //reinitialize and RNA derived data.
		//Reinitialize rna_residuetype_ object! This also make sure rna_residuetype_ didn't inherit anything from the previous update!
		//It appears that the rna_residuetype_ is shared across multiple ResidueType object, if the rna_residuetype_ is not reinitialized here!

		rna_residuetype_ = *( core::chemical::rna::RNA_ResidueTypeOP( new core::chemical::rna::RNA_ResidueType ) );

		//update_last_controlling_chi is treated seperately for RNA case. Parin Sripakdeevong, June 26, 2011
		rna_residuetype_.rna_update_last_controlling_chi(this, last_controlling_chi_, atoms_last_controlled_by_chi_);

		rna_residuetype_.update_derived_rna_data(this);
    } else if (is_carbohydrate_) {
        carbohydrate_info_ = new chemical::carbohydrates::CarbohydrateInfo(this);
        update_last_controlling_chi();
	} else {
		update_last_controlling_chi();
	}

}

///////////////////////////////////////////////////////////////////////////////

/// recalculate derived data, potentially reordering atom-indices
/**
	 This routine updates all the derived data.\n
	 Atom order will probably change after this call, so if you add a new
	 property that depends on atom-indices that will have to be updated below.
**/
/*
	 data that we have prior to calling this routine:

	 name                   type                 setting method
	 ----------------------------------------------------------
	 ordered_atoms_                 v1<AtomAP>            add_atom //from base class
	 atom_name_             v1<string>           add_atom
	 atom_graph_index_            map<string,VD>      add_atom
	 atomic_charge          v1<Real>             add_atom
	 bonded_neighbor_       v1<v1<int>>          add_bond
	 bonded_neighbor_type   v1<v1<BondName>>     add_bond
	 atom_base_             v1<int>              set_atom_base
	 chi_atoms_             v1<v1<int>>          add_chi
	 properties             bools                add_property
	 nbr_atom_              int                  nbr_atom( int )

	 This routine updates all the derived data.

	 Atoms_ order will probably change after this call, so if you add a new
	 property that depends on atom-indices that will have to be updated below.
*/
void
ResidueType::finalize()
{

	AtomIndices old2new;

	setup_atom_ordering( old2new );

	reorder_primary_data( old2new );

	update_derived_data();

	// debug -- these temporary arrays should have been cleared
	assert( force_bb_.empty() && delete_atoms_.empty() );

	// signal that derived data is up to date now
	finalized_ = true;

}

////////////////////////////////////////////////////////////////////


bool
ResidueType::variants_match( ResidueType const & other ) const
{

	if ( ! basic::options::option[ basic::options::OptionKeys::pH::pH_mode ].user() ) {
		for ( Size ii = 1; ii <= variant_types_.size(); ++ii ) {
			if ( ! other.has_variant_type( variant_types_[ ii ] ) ) {
				return false;
			}
		}
		return (variant_types_.size() == other.variant_types_.size());
	}

	//needed for protonated versions of the residues
	else {
		int this_variant_count_offset( 0 );
		for ( Size ii = 1; ii <= variant_types_.size(); ++ii ) {
			if ( variant_types_[ii] == PROTONATED || variant_types_[ii] == DEPROTONATED ) {
				this_variant_count_offset = 1;
				continue;
			}
			if ( ! other.has_variant_type( variant_types_[ ii ] ) ) {
				return false;
			}
		}

		int other_variant_count_offset( 0 );
		if( other.has_variant_type( PROTONATED ) || other.has_variant_type( DEPROTONATED ) ) {
			other_variant_count_offset = 1;
		}

		return ( ( variant_types_.size() - this_variant_count_offset ) ==
							  ( other.variant_types_.size() - other_variant_count_offset ) );
	}
}

bool
ResidueType::nonadduct_variants_match( ResidueType const & other ) const
{
	int this_variant_count_offset( 0 );
	for ( Size ii = 1; ii <= variant_types_.size(); ++ii ) {
		if ( variant_types_[ii] == ADDUCT ) {
			this_variant_count_offset = 1;
			continue;
		}
		if ( ! other.has_variant_type( variant_types_[ ii ] ) ) {
			return false;
		}
	}

		int other_variant_count_offset( 0 );
		if( other.has_variant_type( ADDUCT ) ) {
			other_variant_count_offset = 1;
		}

		return ( ( variant_types_.size() - this_variant_count_offset ) ==
							( other.variant_types_.size() - other_variant_count_offset ) );
}


bool
ResidueType::has_atom_name( std::string const & name ) const
{
	NameVDMap::const_iterator graph_iter
		( atom_graph_index_.find( name ) );
	if ( graph_iter == atom_graph_index_.end() ) {
		return false;
	}
	return true;
}

Size
ResidueType::atom_index( std::string const & name ) const
{
	//// NOTE: Currently we have to iterate twice because atom_graph_index stores vertex_descriptors not indices.
	//// A substantial change to the interface will fix this but everyone's code will need to switch to

	NameVDMap::const_iterator graph_iter
		( atom_graph_index_.find( name ) );
	if ( graph_iter == atom_graph_index_.end() ) {
#if defined BOINC
		// chu temporary graphic fix for boinc
		if ( name == "CA" && !is_protein() ) return 1;
#endif
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		utility_exit_with_message("unknown atom_name: " + name3() + "  " + name );
	}
	VD const & vd = graph_iter->second;

	Size ordered_index=0;
	Size i=1;
	for(; i <= ordered_atoms_.size(); ++i){
		if( std::find(delete_atoms_.begin(),delete_atoms_.end(), i) != delete_atoms_.end()) continue;//this atom is scheduled to be deleted.
		if( &graph_[ordered_atoms_[i]] == &graph_[vd]){
			ordered_index = i;
			break;
		}
	}

	if( ordered_index==0){
#if defined BOINC
		// chu temporary graphic fix for boinc
		if ( name == "CA" && !is_protein() ) return 1;
#endif
		tr.Error << "atom name : " << name << " not available in residue " << name3() << std::endl;
		show_all_atom_names( tr.Error );
		tr.Error << std::endl;
		utility_exit_with_message("unknown atom_name: " + name3() + "  " + name );
	}

	return ordered_index;
}

core::Size
ResidueType::orbital_index( std::string const & name ) const
{

	std::map< std::string, int >::const_iterator iter
		( orbitals_index_.find( name ) );
	if ( iter == orbitals_index_.end() ) {
		utility_exit_with_message("unknown orbital_name: " + name3() + "  " + name );
	}
	return iter->second;
}



void
ResidueType::set_backbone_heavyatom( std::string const & name )
{
	finalized_ = false;
	if ( n_backbone_heavyatoms_ ) {
		assert( force_bb_.empty() );
		for ( Size i=1; i<= n_backbone_heavyatoms_; ++i ) {
			force_bb_.push_back( i );
		}
		n_backbone_heavyatoms_ = 0;
	}
	force_bb_.push_back( atom_index( name ) );
}

/// @brief AtomICoord of an atom
AtomICoor const &
ResidueType::icoor( Size const atm ) const
{
	Atom const & a = graph_[ ordered_atoms_[atm] ]; // atm must become a ResidueGraph::vertex_descriptor
	return a.icoor();
}

Size
ResidueType::add_residue_connection( std::string const & atom_name )
{
	finalized_ = false;

	++n_non_polymeric_residue_connections_;
	residue_connections_.push_back( ResidueConnection( atom_index( atom_name ) ) );
	update_residue_connection_mapping();
	return residue_connections_.size();
}

///@details update actcoord
/** average geometrical center of the set of actcoord_atoms_ */
void
ResidueType::update_actcoord( conformation::Residue & rot ) const
{
	rot.actcoord().zero();
	if ( n_actcoord_atoms_ > 0 ) {
		for ( Size ii = 1; ii <= n_actcoord_atoms_; ++ii )
			{
				rot.actcoord() += rot.atoms()[ actcoord_atoms_[ ii ]].xyz();
			}
		rot.actcoord() /= n_actcoord_atoms_;
	}
}


/// @details set AtomICoor for an atom
///
/// will update the xyz coords as well if desired, useful inside a patching operation where new
/// atoms are being added.
///
void
ResidueType::set_icoor(
	Size const & index,
	std::string const & atm,
	Real const phi,
	Real const theta,
	Real const d,
	std::string const & stub_atom1,
	std::string const & stub_atom2,
	std::string const & stub_atom3,
	bool const update_xyz // = false
)
{
	ICoorAtomID id( atm, *this );
	AtomICoor const ic( index, phi, theta, d, stub_atom1, stub_atom2, stub_atom3, *this );

	Size const atomno( id.atomno() );
	switch ( id.type() ) {
	case ICoorAtomID::INTERNAL:
		if ( graph_.num_vertices() < atomno ) utility_exit_with_message("ResidueType:: shoudnt get here!");//icoor_.resize(atomno);
		if ( ordered_atoms_.size() < atomno ) utility_exit_with_message("ResidueType:: shoudnt get here!");//icoor_.resize(atomno);
		graph_[ordered_atoms_[ atomno ]].icoor( ic );

		// update atom_base?
		if ( ( stub_atom1 != atm ) && has( stub_atom1 ) &&
				 ( atom_base_.size() < atomno || atom_base_[ atomno ] == 0 || atom_base_[ atomno ] == atomno ) ) {
			set_atom_base( atm, stub_atom1 );
		}
		if ( update_xyz ) {
			set_ideal_xyz( atm, ic.build( *this ) );
			//std::cout << "building coords for atm " << name_ << ' ' << atm << ' ' <<
			//		ic.build(*this)(1) << ' ' <<
			//		ic.build(*this)(2) << ' ' <<
			//		ic.build(*this)(3) << std::endl;
		}
		break;
	case ICoorAtomID::CONNECT:
		residue_connections_[ atomno ].icoor( ic );
		break;
	case ICoorAtomID::POLYMER_LOWER:
		assert( lower_connect_id_ != 0 );
		residue_connections_[ lower_connect_id_ ].icoor( ic );
		break;
	case ICoorAtomID::POLYMER_UPPER:
		assert( upper_connect_id_ != 0 );
		residue_connections_[ upper_connect_id_ ].icoor( ic );
		break;
	default:
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
}

void
ResidueType::set_icoor(
	std::string const & atm,
	Real const phi,
	Real const theta,
	Real const d,
	std::string const & stub_atom1,
	std::string const & stub_atom2,
	std::string const & stub_atom3,
	bool const update_xyz // = false
)
{
	ICoorAtomID id( atm, *this );
	AtomICoor const ic( phi, theta, d, stub_atom1, stub_atom2, stub_atom3, *this );

	Size const atomno( id.atomno() );
	switch ( id.type() ) {
	case ICoorAtomID::INTERNAL:
		if ( ordered_atoms_.size() < atomno ) utility_exit_with_message("ResidueType:: shoudnt get here!");//icoor_.resize(atomno);
		graph_[ordered_atoms_[ atomno ]].icoor( ic );
		// update atom_base?
		if ( ( stub_atom1 != atm ) && has( stub_atom1 ) &&
				 ( atom_base_.size() < atomno || atom_base_[ atomno ] == 0 || atom_base_[ atomno ] == atomno ) ) {
			set_atom_base( atm, stub_atom1 );
		}
		if ( update_xyz ) {
			set_ideal_xyz( atm, ic.build( *this ) );
			//std::cout << "building coords for atm " << name_ << ' ' << atm << ' ' <<
			//		ic.build(*this)(1) << ' ' <<
			//		ic.build(*this)(2) << ' ' <<
			//		ic.build(*this)(3) << std::endl;
		}
		break;
	case ICoorAtomID::CONNECT:
		residue_connections_[ atomno ].icoor( ic );
		break;
	case ICoorAtomID::POLYMER_LOWER:
		assert( lower_connect_id_ != 0 );
		residue_connections_[ lower_connect_id_ ].icoor( ic );
		break;
	case ICoorAtomID::POLYMER_UPPER:
		assert( upper_connect_id_ != 0 );
		residue_connections_[ upper_connect_id_ ].icoor( ic );
		break;
	default:
		utility_exit_with_message( "unrecognized stub atom id type!" );
	}
}

//set the orbital icoor data.
void
ResidueType::set_orbital_icoor_id(
	std::string const & orbital,
	Real const phi,
	Real const theta,
	Real const d,
	std::string const & stub_atom1,
	std::string const & stub_atom2,
	std::string const & stub_atom3
)
{
	Size orb_indx(orbital_index(orbital));
	std::string stub1(stub_atom1);
	std::string stub2(stub_atom2);
	std::string stub3(stub_atom3);
	orbitals::ICoorOrbitalData icoor(phi, theta, d, stub1, stub2, stub3);
	orbitals_[ orb_indx ].icoor( icoor );

	core::Size s1(atom_index( stub_atom1 ));
	core::Size s2(atom_index( stub_atom2 ));
	core::Size s3(atom_index( stub_atom3 ));
	orbitals::ICoorOrbitalData new_icoor(phi, theta, d, s1, s2, s3);

	orbitals_[ orb_indx ].new_icoor( new_icoor );

}


void ResidueType::set_RotamerLibraryName( std::string const & filename )
{
	rotamer_library_name_ = filename;
}

/// @brief A residue parameter file can refer to a set of "pdb rotamers" that can be
/// superimposed onto a starting position for use in the packer.  These rotamers
/// are loaded into the pack::dunbrack::RotamerLibrary at the time of their first use.
std::string ResidueType::get_RotamerLibraryName() const
{
	return rotamer_library_name_;
}


void ResidueType::assign_neighbor_atom()
{
	//calculate the geometric center of all atoms in the residue

	Vector total(0.0,0.0,0.0);
	for(core::Size index = 1; index <= ordered_atoms_.size();++index)
	{
		total += graph_[ordered_atoms_[index]].ideal_xyz();
	}

	Vector center = total/ordered_atoms_.size();

	//locate the atom which is closest to the center
	core::Size min_index = 0;
	core::Real min_distance = 50000.0;

	for(core::Size index = 1; index <= ordered_atoms_.size();++index)
	{
		core::Real distance = center.distance(graph_[ordered_atoms_[index]].ideal_xyz());
		if( (distance < min_distance) && (!atom_is_hydrogen(index)) )
		{
			min_distance = distance;
			min_index = index;
		}
	}
	assert(min_index != 0);
	//set neighbor atom
	nbr_atom(graph_[ordered_atoms_[min_index]].name());
}

void ResidueType::assign_internal_coordinates()
{
	assign_internal_coordinates( atom_name(nbr_atom_) );
	AtomICoor nbr_icoor = icoor(nbr_atom_);
	set_icoor(atom_name(nbr_atom_),0.0,0.0,0.0,atom_name(nbr_icoor.stub_atom1().atomno()),
			atom_name(nbr_icoor.stub_atom2().atomno()),atom_name(nbr_icoor.stub_atom3().atomno()));

}

void ResidueType::assign_internal_coordinates(std::string const & current_atom)
{
	//%TODO: right now i'm ignoring M FRAG lines and M SPLT lines in molfiles
	core::Size current_atom_index = atom_index(current_atom);
	std::string parent_stub1;
	std::string parent_stub2;
	std::string parent_stub3;

	//the root atom has dummy stubs and icoords of 0
	if(current_atom_index == nbr_atom_)
	{
		core::Size first_child = bonded_neighbor_[current_atom_index][1];
		parent_stub1 = atom_name(current_atom_index);
		parent_stub2 = atom_name(first_child);

		if(bonded_neighbor_[first_child].size() >	0)
		{
			parent_stub3 = atom_name(bonded_neighbor_[first_child][1]);
		}else
		{
			parent_stub3 = atom_name(bonded_neighbor_[current_atom_index][2]);
		}
	}else
	{
		core::Size parent_index = parents_[current_atom_index];
		AtomICoor parent_icoor = icoor(parent_index);
		parent_stub1 = atom_name(current_atom_index);
		parent_stub2 = atom_name(parent_index);
		parent_stub3 = atom_name(parent_icoor.stub_atom2().atomno());
	}

	std::string previous_sibling = parent_stub3;
	AtomIndices children = bonded_neighbor_[current_atom_index];
	for(core::Size index = 1; index <children.size();++index)
	{
		core::Size child_index = children[index];

		if((parents_[child_index] != 0) && (child_index != nbr_atom_))
		{
			continue;
		}

		std::string child_stub1 = parent_stub1;
		std::string child_stub2 = parent_stub2;
		std::string child_stub3 = previous_sibling;
		parents_[child_index] = current_atom_index;
		if((current_atom_index == nbr_atom_) && (previous_sibling == parent_stub2))
		{
			child_stub3 = parent_stub3;
		}
		calculate_icoor(atom_name(child_index),child_stub1,child_stub2,child_stub3);
		//set_atom_base(atom_name(child_index),)
		assign_internal_coordinates(atom_name(child_index) );
		previous_sibling = atom_name(child_index);
	}
}

void ResidueType::calculate_icoor(std::string const & child,
		std::string const & stub_atom1,
		std::string const & stub_atom2,
		std::string const & stub_atom3)
{
	//std::cout <<child << " \""<<stub_atom1 << "\" \""<<stub_atom2<< "\" \""<<stub_atom3 << std::endl;
	// This is basically a direct port of calc_internal_coords()
	// found in /python/apps/public/molfile_to_params.py
	Vector const child_xyz = atom(atom_index(child)).ideal_xyz();
	Vector const stub1_xyz = atom(atom_index(stub_atom1)).ideal_xyz();
	Vector const stub2_xyz = atom(atom_index(stub_atom2)).ideal_xyz();
	Vector const stub3_xyz = atom(atom_index(stub_atom3)).ideal_xyz();

	core::Real distance = child_xyz.distance(stub1_xyz);
	core::Real theta = 0.0;
	core::Real phi = 0.0;
	if(distance <1e-2)
	{
		tr << "WARNING: extremely small distance=" << distance << " for " <<
				child << " ,using 0.0 for theta and phi."<<
				" If you were not expecting this warning, something is very wrong" <<std::endl;
	}else
	{
		theta = numeric::angle_radians<core::Real>(child_xyz,stub1_xyz,stub2_xyz);
		if( (theta < 1e-2) || (theta > numeric::NumericTraits<Real>::pi()-1e-2) )
		{
			phi = 0.0;
		}else
		{
			phi = numeric::dihedral_radians<core::Real>(child_xyz,stub1_xyz,stub2_xyz,stub3_xyz);
		}

	}
	//tr << child << " " << stub_atom1 << " "<< stub_atom2 << " " <<stub_atom3 << " " <<distance << " " << phi << " " << theta <<std::endl;
	set_icoor(child,phi,theta,distance,stub_atom1,stub_atom2,stub_atom3);
}

void
ResidueType::set_ideal_xyz(
	std::string const & atm,
	Vector const & xyz_in
)
{
	Size const index( atom_index(atm) );
	set_ideal_xyz(index,xyz_in);
}

void
ResidueType::set_ideal_xyz(
	Size index,
	Vector const & xyz_in
)
{
	if ( index > ordered_atoms_.size() ) ordered_atoms_.resize(index); ///TODO REMOVE THIS!!! This would create NULL vertex_descriptors
	Atom & a = graph_[ ordered_atoms_[index] ];
	a.ideal_xyz( xyz_in );
}

// Return the CarbohydrateInfo object containing sugar-specific properties for this residue.
core::chemical::carbohydrates::CarbohydrateInfoCOP
ResidueType::carbohydrate_info() const
{
    return carbohydrate_info_;
}

void
ResidueType::print_dihedrals() const
{
	tr.Debug << "START DIHEDRAL ANGLES ATOM NAMES" << std::endl;
	tr.Debug << "Number of dihe: " << ndihe_ << " " << dihedral_atom_sets_.size() << std::endl;
	for ( Size i = 1; i <= ndihe_; ++i )
		{

			AtomType at1 = atom_type( dihedral_atom_sets_[ i ].key1() );
			AtomType at2 = atom_type( dihedral_atom_sets_[ i ].key2() );
			AtomType at3 = atom_type( dihedral_atom_sets_[ i ].key3() );
			AtomType at4 = atom_type( dihedral_atom_sets_[ i ].key4() );
			MMAtomType at5 = mm_atom_type( dihedral_atom_sets_[ i ].key1() );
			MMAtomType at6 = mm_atom_type( dihedral_atom_sets_[ i ].key2() );
			MMAtomType at7 = mm_atom_type( dihedral_atom_sets_[ i ].key3() );
			MMAtomType at8 = mm_atom_type( dihedral_atom_sets_[ i ].key4() );

			tr.Debug << "PDB:" << "\t"
								<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key1() ]].name() << "\t"
								<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key2() ]].name() << "\t"
								<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key3() ]].name() << "\t"
								<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key4() ]].name() << "\t"
								<< "MM:" << "\t"
								<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key1() ]].mm_name() << "\t"
								<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key2() ]].mm_name() << "\t"
								<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key3() ]].mm_name() << "\t"
								<< graph_[ordered_atoms_[ dihedral_atom_sets_[ i ].key4() ]].mm_name() << "\t"
								<< "MM2:" << "\t"
								<< at5.name() << "\t"
								<< at6.name() << "\t"
								<< at7.name() << "\t"
								<< at8.name() << "\t"
								<< "ROS:" << "\t"
								<< at1.name() << "\t"
								<< at2.name() << "\t"
								<< at3.name() << "\t"
								<< at4.name() << "\t"
								<< std::endl;
		}
	tr.Debug << "END DIHEDRAL ANGLES ATOM NAMES" << std::endl;
}

void
ResidueType::print_bondangles() const
{
	tr.Debug << "START BOND ANGLES ATOM NAMES" << std::endl;
	tr.Debug << "Number of bond angles: " << bondangle_atom_sets_.size() << std::endl;
	for ( Size i = 1; i <= bondangle_atom_sets_.size(); ++i )
		{

			AtomType at1 = atom_type( bondangle_atom_sets_[ i ].key1() );
			AtomType at2 = atom_type( bondangle_atom_sets_[ i ].key2() );
			AtomType at3 = atom_type( bondangle_atom_sets_[ i ].key3() );
			MMAtomType at5 = mm_atom_type( bondangle_atom_sets_[ i ].key1() );
			MMAtomType at6 = mm_atom_type( bondangle_atom_sets_[ i ].key2() );
			MMAtomType at7 = mm_atom_type( bondangle_atom_sets_[ i ].key3() );

			tr.Debug << "PDB:" << "\t"
								<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key1() ]].name() << "\t"
								<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key2() ]].name() << "\t"
								<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key3() ]].name() << "\t"
								<< "MM:" << "\t"
								<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key1() ]].mm_name() << "\t"
								<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key2() ]].mm_name() << "\t"
								<< graph_[ordered_atoms_[ bondangle_atom_sets_[ i ].key3() ]].mm_name() << "\t"
								<< "MM2:" << "\t"
								<< at5.name() << "\t"
								<< at6.name() << "\t"
								<< at7.name() << "\t"
								<< "ROS:" << "\t"
								<< at1.name() << "\t"
								<< at2.name() << "\t"
								<< at3.name() << "\t"
								<< std::endl;
		}
	tr.Debug << "END BOND ANGLES ATOM NAMES" << std::endl;
}

void
ResidueType::print_pretty_path_distances() const
{
	tr.Debug << "START PATH DISTANCES" << std::endl;
	// print header line
	for ( Size i = 1; i <= natoms_; ++i )
		{
			tr.Debug << "\t" << graph_[ordered_atoms_[i]].name();
		}
	tr.Debug << std::endl;

	for ( Size j = 1; j <= natoms_; ++j )
		{
			tr.Debug << graph_[ordered_atoms_[j]].name() << "\t";
			for ( Size k = 1; k <= natoms_; ++k )
				{
					tr.Debug << path_distance_[j][k] << "\t";
				}
			tr.Debug << std::endl;
		}
	tr.Debug << "END PATH DISTANCES" << std::endl;
}

void
ResidueType::update_residue_connection_mapping()
{
	//std::fill( atom_2_residue_connection_map_.begin(), atom_2_residue_connection_map_.end(), 0 );
	for ( Size ii = 1; ii <= natoms_; ++ii ) { atom_2_residue_connection_map_[ ii ].clear(); }

	for ( Size ii = 1; ii <= residue_connections_.size(); ++ii ) {
		atom_2_residue_connection_map_[ residue_connections_[ ii ].atomno() ].push_back( ii );
		residue_connections_[ ii ].index( ii );
	}
}

void
ResidueType::update_last_controlling_chi() {
	last_controlling_chi_.resize( natoms_ );
	std::fill( last_controlling_chi_.begin(), last_controlling_chi_.end(), 0 );

	/// 1. First we have to mark all the atoms who are direct descendants of the 3rd
	/// atom in each chi; this prevents the note_chi_controls_atom recursion from overwriting
	/// the last-controlling chi for atoms descending from a particular chi.
	for ( Size ii = 1; ii <= nchi(); ++ii ) {
		Size const iiat3 = chi_atoms_[ ii ][ 3 ];
		// This may be unnecessary; I believe two atoms pair as each other's bases only at the mainchain.
		Size const iiat3base = atom_base_[ iiat3 ];
		AtomIndices const & ii_nbrs( bonded_neighbor_[ iiat3 ] );
		for ( Size jj = 1; jj <= ii_nbrs.size(); ++jj ) {
			Size const jj_atom = ii_nbrs[ jj ];
			if ( atom_base_[ jj_atom ] == iiat3 && iiat3base != jj_atom ) {
				last_controlling_chi_[ jj_atom ] = ii;
			}
		}
	}

	/// 2. Now, lets recurse through all the atoms that are not direct descendants
	/// of the 3rd atom in a chi.  E.g. chi2 in PHE controls several more atoms than
	/// just CD1 and CD2.
	for ( Size ii = nchi(); ii >= 1; --ii ) {
		/// Note children of atom 3 of chi_ii as being controlled by chi ii.
		Size const iiat3 = chi_atoms_[ ii ][ 3 ];
		// This may be unnecessary; I believe two atoms pair as each other's bases only at the mainchain.
		Size const iiat3base = atom_base_[ iiat3 ];
		AtomIndices const & ii_nbrs( bonded_neighbor_[ iiat3 ] );
		for ( Size jj = 1; jj <= ii_nbrs.size(); ++jj ) {
			Size const jj_atom = ii_nbrs[ jj ];
			if ( atom_base_[ jj_atom ] == iiat3 && iiat3base != jj_atom ) {
				note_chi_controls_atom( ii, jj_atom );
			}
		}
	}

	/// Now compute the atoms_last_controlled_by_chi_ arrays.

	/// get ready to allocate space in the atoms_last_controlled_by_chi_ arrays
	utility::vector1< Size > natoms_for_chi( nchi(), 0 );
	for ( Size ii = 1; ii <= natoms_; ++ii ) {
		if ( last_controlling_chi_[ ii ] != 0 ) {
			++natoms_for_chi[ last_controlling_chi_[ ii ] ];
		}
	}

	/// allocate space
	atoms_last_controlled_by_chi_.resize( nchi() );
	for ( Size ii = 1; ii <= nchi(); ++ii ) {
		atoms_last_controlled_by_chi_[ ii ].clear();
		atoms_last_controlled_by_chi_[ ii ].reserve( natoms_for_chi[ ii ] );
	}

	/// fill the arrays
	for ( Size ii = 1; ii <= natoms_; ++ii ) {
		if ( last_controlling_chi_[ ii ] != 0 ) {
			atoms_last_controlled_by_chi_[ last_controlling_chi_[ ii ]].push_back( ii );
		}
	}

}

/// @details O(N) recursive algorithm for determining the last chi for each atom.
/// Each atom is visited at most twice.
void
ResidueType::note_chi_controls_atom( Size chi, Size atomno )
{
	/// This should never be called on the "root" atom or it will enter an infinite loop
	assert( atom_base_[ atom_base_[ atomno ]] != atomno );

	/// End the recursion: this atom already has had it's last chi identified, and it's not
	/// the chi we're currently labeling atoms with.
	if ( last_controlling_chi_[ atomno ] != 0 && last_controlling_chi_[ atomno ] != chi ) return;

	last_controlling_chi_[ atomno ] = chi;

	AtomIndices const & nbrs( bonded_neighbor_[ atomno ] );
	for ( Size ii = 1; ii <= nbrs.size(); ++ii ) {
		/// descend into atoms who list atomno as their parent;
		/// atom_base_ defines a tree except at the root, where
		/// atom_base_[ atom_base_[ ii ]] == ii
		if ( atom_base_[ nbrs[ ii ]] == atomno ) {
			note_chi_controls_atom( chi, nbrs[ ii ] );
		}
	}

}

void
ResidueType::select_orient_atoms(
	Size & center,
	Size & nbr1,
	Size & nbr2
) const
{
	center = 0;
	nbr1 = 0;
	nbr2 = 0;

	// No backbone atoms, all backbone atoms, or orient mode explicitly set to nbr_atom
	if ( first_sidechain_atom() == 1 || first_sidechain_atom() > natoms() || force_nbr_atom_orient() ) {
		// If no backbone atoms (or all bb atoms), assume nbr_atom will be close to center-of-mass.
		center = nbr_atom();
		// If is hydrogen or too few neighbors, try trekking up the atom tree
		while( center > nheavyatoms() || bonded_neighbor(center).size() < 2 ) {
			center = atom_base(center);
		}
		AtomIndices const & nbrs( bonded_neighbor(center) );
		// First try to find two neighbors that are heavyatoms
		for( Size j=1; j<= nbrs.size(); ++j ) {
			Size const nbr( nbrs[j] );
			if( nbr <= nheavyatoms() ) {
				if ( nbr1 ) nbr2 = nbr;
				else nbr1 = nbr;
			}
		}
		// Failing that, just try for two neighbors!
		if( !( center && nbr1 && nbr2 ) ) {
			for( Size j=1; j<= nbrs.size(); ++j ) {
				Size const nbr( nbrs[j] );
				if ( nbr1 ) nbr2 = nbr;
				else nbr1 = nbr;
			}
		}
		if( !( center && nbr1 && nbr2 ) ) {
			// assert() isn't enough for these cases b/c they're typically ligands
			// and thus depend on user input -- need to be caught even in release mode.
			utility_exit_with_message("Cannot superimpose residues of type "+name());
		}
		//std::cout << "Superimposing on " << atom_name(center) << " " << atom_name(nbr1) << " " << atom_name(nbr2) << "\n";

	} else {
		// look for a backbone atom, one of whose neighbors is a sidechain atom
		// center will be this atom
		// nbr1 and nbr2 will be the backbone heavyatom nbrs of this atom
		// eg center = CA, nbr1 = N. nbr2 = C in the protein case
		for ( Size atom_index(1); atom_index <= natoms(); ++atom_index ) {
			if ( atom_is_backbone( atom_index ) ) {
				AtomIndices const & nbrs( bonded_neighbor( atom_index ) );
				center = 0; nbr1 = 0; nbr2 = 0;
				for ( Size nbr_index(1); nbr_index <= nbrs.size(); ++nbr_index ) {
					Size const nbr( nbrs[ nbr_index ] );
					if ( !atom_is_backbone( nbr ) && atom_base( nbr ) == atom_index ) {
						// nbr is a sidechain atom that branches from the atom at atom_index
						center = atom_index;
					} else if ( atom_is_backbone( nbr ) && nbr <= nheavyatoms() ) {
						// nbr is a backbone heavy atom neighbor of the atom at atom_index
						if ( nbr1 ) nbr2 = nbr;
						else nbr1 = nbr;
					}
				}
			} // atom_index is backbone
			if ( center && nbr1 && nbr2 ) break;
		} // atom_index
	}
}

void
ResidueType::report_adducts()
{
	if( defined_adducts_.size() == 0 ) return;

	for( Size ii = 1 ; ii <= defined_adducts_.size() ; ++ii) {
		Adduct & add( defined_adducts_[ii] );
		tr.Debug << "Residue: " << name3() << " Adduct: " << add.adduct_name() <<
			" Atom name: " << add.atom_name() << std::endl;
	}
}

void
ResidueType::debug_dump_icoor()
{

	tr.Debug << "ICoor for " << name3() << std::endl;
	for( Size ii = 1 ; ii <= natoms() ; ++ii) {
		tr.Debug << " Atom name: " << atom_name( ii ) << " ideal xyz " << atom(ii).ideal_xyz()[0] << "  " << atom(ii).ideal_xyz()[1] << "  " << atom(ii).ideal_xyz()[2] << std::endl;
	}

}


void
ResidueType::show_all_atom_names( std::ostream & out ) const {

	for( VIterPair vp = boost::vertices(graph_); vp.first != vp.second; ++vp.first){
		VIter v_iter= vp.first;
		VD vd = *v_iter;
		Atom a = graph_[vd];
		out << a.name() << std::endl;
	}

}

void
ResidueType::set_ncaa_rotlib_n_bin_per_rot( utility::vector1<Size> n_bins_per_rot )
{
	assert( ncaa_rotlib_n_rots_ == n_bins_per_rot.size() );
	ncaa_rotlib_n_bins_per_rot_.resize( ncaa_rotlib_n_rots_ );
	for( Size i = 1; i <= ncaa_rotlib_n_rots_; ++i ) {
		ncaa_rotlib_n_bins_per_rot_[i] = n_bins_per_rot[i];
	}
}

/// @brief  Check if atom is virtual.
bool
ResidueType::is_virtual( Size const & atomno ) const
{
	return ( atom_type( atomno ).is_virtual() );
}

/// @brief  Check if residue is 'VIRTUAL_RESIDUE'
bool
ResidueType::is_virtual_residue() const{
	return ( has_variant_type( "VIRTUAL_RESIDUE" ) );
}


///////////////////////////////////////////////////////////////////////////////

} // chemical
} // core
