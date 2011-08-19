// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/etrie/TrieCountPairAll.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/etable/etrie/TrieCountPairAll.hh>

// Package Headers
#include <core/scoring/trie/trie_vs_trie.hh>
#include <core/scoring/trie/trie_vs_path.hh>
#include <core/scoring/etable/etrie/CountPairData_1_1.hh>
#include <core/scoring/etable/etrie/CountPairData_1_2.hh>
#include <core/scoring/etable/etrie/CountPairData_1_3.hh>
#include <core/scoring/etable/etrie/CountPairDataGeneric.hh>

// Project Headers
#include <core/scoring/etable/EtableEnergy.hh>
//XRW_B_T1
//#include <core/scoring/etable/CoarseEtableEnergy.hh>
//XRW_E_T1
#include <core/scoring/hackelec/HackElecEnergy.hh>
#include <core/scoring/methods/MMLJEnergyInter.hh>

// STL Headers
#include <iostream>

//Auto Headers
#include <core/scoring/etable/etrie/EtableAtom.hh>
#include <core/scoring/trie/RotamerTrie.hh>


namespace core {
namespace scoring {
namespace etable {
namespace etrie {

using namespace trie;

TrieCountPairAll::~TrieCountPairAll() {}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::EtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::EtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}



//XRW_B_T1
/*
/////////////////////////// CoarseEtableEnergy //////////////////////////////////

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_1 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_2 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairData_1_3 > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< EtableAtom, CountPairDataGeneric > const & trie2,
	etable::CoarseEtableEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

*/
//XRW_E_T1

// HBONDS
void
TrieCountPairAll::resolve_trie_vs_trie(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & ,
	ObjexxFCL::FArray2D< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::TrieCountPairAll::resolve_trie_vs_trie reached with HBondEnergy" );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	RotamerTrie< hbonds::hbtrie::HBAtom, hbonds::hbtrie::HBCPData >  const & ,
	hbonds::HBondEnergy const & ,
	utility::vector1< core::PackerEnergy > & ,
	utility::vector1< core::PackerEnergy > & )
{
	utility_exit_with_message( "etable::etrie::TrieCountPairAll::resolve_trie_vs_path reached with HBondEnergy" );
}

/// Hack Elec E
void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void
TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_1 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_2 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairData_1_3 > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}


void
TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< hackelec::ElecAtom, CountPairDataGeneric > const & trie2,
	hackelec::HackElecEnergy const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

/////////////////////////////// MMLJEnergyInter ////////////////////////////

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_trie(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	ObjexxFCL::FArray2D< core::PackerEnergy > & pair_energy_table,
	ObjexxFCL::FArray2D< core::PackerEnergy > & temp_table
)
{
	trie_vs_trie( trie1, trie2, *this, sfxn, pair_energy_table, temp_table );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_1 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_2 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairData_1_3 > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void TrieCountPairAll::resolve_trie_vs_path(
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie1,
	trie::RotamerTrie< mm::mmtrie::MMEnergyTableAtom, CountPairDataGeneric > const & trie2,
	methods::MMLJEnergyInter const & sfxn,
	utility::vector1< core::PackerEnergy > & pair_energy_vector,
	utility::vector1< core::PackerEnergy > & temp_vector
)
{
	trie_vs_path( trie1, trie2, *this, sfxn, pair_energy_vector, temp_vector );
}

void
TrieCountPairAll::print()
{
	std::cout << "TrieCountPairAll" << std::endl;
}

} // namespace etrie
} // namespace etable
} // namespace scoring
} // namespace core

