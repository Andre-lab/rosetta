 // -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/etable/count_pair/CountPairIntraResC3.hh
/// @brief  Count pair for residue pairs connected with one bond, where the
/// crossover from excluding to counting atom pair interactions is at 3 bonds.
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_scoring_etable_count_pair_CountPairIntraResC3_hh
#define INCLUDED_core_scoring_etable_count_pair_CountPairIntraResC3_hh

#include <core/scoring/etable/count_pair/CountPairCrossover3.hh>

#include <core/types.hh>
#include <core/conformation/Atom.hh>

namespace core {
namespace scoring {
namespace etable {
namespace count_pair {

class CountPairIntraResC3 : public CountPairCrossover3
{
public:
	public:
	typedef CountPairCrossover3 parent;

public:
	CountPairIntraResC3(
		conformation::Residue const & res
	);

	virtual ~CountPairIntraResC3();

	///@brief function required by templated functions in atom_pair_energy_inline
	inline
	bool
	operator () (
		int const at1,
		int const at2,
		Real & weight,
		Size & path_dist
	) const
	{
		path_dist = path_dists_[ at1 ][ at2 ];
		return count_at_path_distance( path_dist, weight );
	}

	virtual
	bool
	count(
		int const at1,
		int const at2,
		Real &,
		Size & path_dist
	) const;

	/// Type Resolution Functions ///

	virtual
	void
	residue_atom_pair_energy(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		etable::EtableEnergy const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::EtableEnergy const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::EtableEnergy const &,
		EnergyMap &
	) const;

//XRW_B_T1
/*

	virtual
	void
	residue_atom_pair_energy(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		etable::CoarseEtableEnergy const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_backbone(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::CoarseEtableEnergy const &,
		EnergyMap &
	) const;


	virtual
	void
	residue_atom_pair_energy_sidechain_whole(
		conformation::Residue const &,
		conformation::Residue const &,
		etable::CoarseEtableEnergy const &,
		EnergyMap &
	) const;

*/
//XRW_E_T1

private:
	utility::vector1< utility::vector1< int > > const & path_dists_;

};

} // namespace count_pair
} // namespace etable
} // namespace scoring
} // namespace core

#endif
