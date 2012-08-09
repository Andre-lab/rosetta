// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin ChemicalManager
///
/// @brief
/// Chemical manager class
///
/// @detailed
/// The Chemical Manager is a singleton class, which means that it can only been initialized once (exist once in memory). Once initialized,
/// you can call it by simply access it via:
///
/// core::chemical::AtomTypeSetCAP atom_types =
/// core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
///
/// You can substitute AtomTypeSet, with whatever is seen below (residue_type_set, mm_atom_type_set, orbital_type_set).
/// In the below functions, the "tag_in" refers to fullatom, centroid, which basically tells what type of set to load in.
/// The chemical manager will call functions within the AtomTypeSet, MMAtomTypeSet, ResidueTypeSet, etc etc. The classes type set
/// reads in files from the database to create atom types, residue types, and mmatom types. The information from those files are stored
/// in the type class.
///
///
///
/// @authors
/// Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// Steven Combs - comments
///
///
/// @last_modified December 6 2010
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_core_chemical_ChemicalManager_hh
#define INCLUDED_core_chemical_ChemicalManager_hh

// Unit headers
#include <core/chemical/ChemicalManager.fwd.hh>

// Package headers
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/IdealBondLengthSet.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
//#include <core/chemical/CSDAtomTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

// C++
#include <map>

// Commented by inclean daemon #include <string>

namespace core {
namespace chemical {

/// @brief a class managing different sets of atom_type_set and residue_type_set
///
/// @details make it as a singleton class so that atom_type_set and residue_type_set are only
///input and initialized once. They can be later retrieved by querying this class.
class ChemicalManager
{
public:


public:
	static ChemicalManager * get_instance();

	/// @brief query atom_type_set by a name tag
	AtomTypeSetCAP
	atom_type_set( std::string const & tag );

	/// @brief query atom_type_set by a name tag
	ElementSetCAP
	element_set( std::string const & tag );

	/// @brief query ideal_bond_lengths
	IdealBondLengthSetCAP ideal_bond_length_set(std::string const & tag);

	/// @brief query mm_atom_type_set by a name tag
	MMAtomTypeSetCAP
	mm_atom_type_set( std::string const & tag );

	/// @brief query csd_atom_type_set by a name tag
//	CSDAtomTypeSetCAP
//	csd_atom_type_set( std::string const & tag );

	/// @brief query orbital_type_set by a name tag
	orbitals::OrbitalTypeSetCAP
	orbital_type_set(std::string const & tag);

	/// @brief query residue_type_set by a name tag
	ResidueTypeSetCAP
	residue_type_set( std::string tag );

	/// @brief query residue_type_set by a name tag
	ResidueTypeSet &
	nonconst_residue_type_set( std::string const & tag );



private:
	typedef std::map< std::string, AtomTypeSetOP > AtomTypeSets;
	typedef std::map< std::string, ElementSetOP > ElementSets;
	typedef std::map< std::string, IdealBondLengthSetOP> IdealBondLengthSets;
	typedef std::map< std::string, orbitals::OrbitalTypeSetOP > OrbitalTypeSets;
	typedef std::map< std::string, MMAtomTypeSetOP > MMAtomTypeSets;
//	typedef std::map< std::string, CSDAtomTypeSetOP > CSDAtomTypeSets;
	typedef std::map< std::string, ResidueTypeSetOP > ResidueTypeSets;

private:
	/// @brief private constructor
	ChemicalManager();
	/// @brief static data member holding pointer to the singleton class itself
	static ChemicalManager * instance_;

	/// @brief lookup map for querying atom_type_set by name tag
	AtomTypeSets atom_type_sets_;
	/// @brief lookup map for querying element_type_set by name tag
	ElementSets element_sets_;
	///@brief lookup map for querying orbital_type_set by name tag.
	OrbitalTypeSets orbital_type_sets_;
	/// @brief lookup map for querying mm_atom_type_set by name tag
	MMAtomTypeSets mm_atom_type_sets_;
	/// @brief lookup map for querying csd_atom_type_set by name tag
//	CSDAtomTypeSets csd_atom_type_sets_;
	/// @brief lookup map for querying residue_type_set by name tag
	ResidueTypeSets residue_type_sets_;

	IdealBondLengthSets ideal_bond_length_sets_;
};

} // namespace core
} // namespace chemical


#endif

// 	/// @brief  Duplicate a ResidueTypeSet, preparatory to modifying it in some way DOES NOT DUPLICATE ITS ATOMTYPESETS
// 	void
// 	copy_residue_type_set(
// 		std::string const & old_name,
// 		std::string const & new_name
// 	);

// 	/// @brief  Duplicate an AtomTypeSet, preparatory to modifying it in some way
// 	void
// 	copy_atom_type_set(
// 		std::string const & old_name,
// 		std::string const & new_name
// 	);
