// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/HolesEnergyRes.cc
/// @brief  Packing Score
/// @author Will Sheffler


//Unit headers
#include <core/scoring/packing/HolesEnergyRes.hh>

#include <core/scoring/packing/PoseBalls.hh>

//Package headers

//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/methods/EnergyMethodCreator.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/id/CacheableAtomID_MapVector.hh>
#include <basic/options/option.hh>

#include <core/scoring/packing/compute_holes_score_res.hh>

// AUTO-REMOVED #include <basic/prof.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

//utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

//C++ headers
#include <iostream>

//Auto Headers
#include <core/pose/util.hh>
//#include <core/import_pose/import_pose.hh>
#include <core/scoring/EnergyMap.hh>


namespace core {
namespace scoring {
namespace packing {


//////////////////////////////////////////////////////
//@brief THIS CLASS REQUIRES AN ENERGY METHOD CREATOR TO WORK PROPERLY
//////////////////////////////////////////////////////
HolesEnergyRes::HolesEnergyRes() :
	parent( 0 ) // THIS WILL FAIL
{
	//add_score_type( scoring::holes );
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
void
HolesEnergyRes::finalize_total_energy(
  pose::Pose & pose,
  ScoreFunction const &,
  EnergyMap & totals
) const
{
	totals[ holes ] = compute_holes_score_res(pose,params_);
	//
	//
	// pose::Pose nat;
	// core::import_pose::pose_from_pdb(nat,"hivp_high_1tw7.pdb.clean.align");
	//
	// Real rms = all_atom_rmsd(pose,nat);
	//
	// std::cout << "holes score: " << totals[holes]/(pose.total_residue()-2) << " rms " << rms << std::endl;
}

void
HolesEnergyRes::setup_for_derivatives(
	pose::Pose & pose,
	ScoreFunction const &
) const {
	std::cerr << "HolesEnergyRes::setup_for_derivatives" << std::endl;

	using namespace basic;
	using namespace datacache;
	using namespace id;
	using namespace numeric;
	if( !pose.data().has( core::pose::datacache::CacheableDataType::HOLES_POSE_INFO ) ) {
		pose.data().set( core::pose::datacache::CacheableDataType::HOLES_POSE_INFO, new CacheableAtomID_MapVector );
	}
	CacheableDataOP dat( pose.data().get_ptr( core::pose::datacache::CacheableDataType::HOLES_POSE_INFO ) );
	CacheableAtomID_MapVector *cachemap = (CacheableAtomID_MapVector*)dat();
	AtomID_Map<xyzVector<Real> > & derivs(cachemap->map());
	core::pose::initialize_atomid_map_heavy_only(derivs,pose);

	compute_holes_deriv_res( pose, params_, derivs );

}


void
HolesEnergyRes::eval_atom_derivative(
	id::AtomID const & aid,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	// std::cerr << "HolesEnergyRes::eval_atom_derivative " << aid << std::endl;
	using namespace basic;
	using namespace datacache;
	using namespace id;
	using namespace numeric;
	CacheableDataCOP dat( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::HOLES_POSE_INFO ) );
	CacheableAtomID_MapVector const *cachemap = (CacheableAtomID_MapVector const *)dat();
	AtomID_Map<xyzVector<Real> > const & derivs(cachemap->map());

	if( aid.rsd() > derivs.n_residue() || aid.atomno() > derivs.n_atom(aid.rsd()) ) {
		return;
	}
	// std::cerr << "eval_atom_derivative " << aid << " " << derivs[aid].x() << " " << derivs[aid].y() << " " << derivs[aid].z() << std::endl;
	// F2 += weights * derivs[aid];
	// F1 += weights * derivs[aid].cross(xyzVector<Real>(1,0,0));

  	numeric::xyzVector<core::Real> atom_x = pose.xyz(aid);
	numeric::xyzVector<core::Real> const f2( derivs[aid] );
	numeric::xyzVector<core::Real> const atom_y = atom_x - f2;   // a "fake" atom in the direcion of the gradient
	numeric::xyzVector<core::Real> const f1( atom_x.cross( atom_y ) );

	F1 += weights[ holes ] * f1;
	F2 += weights[ holes ] * f2;
}


core::Size
HolesEnergyRes::version() const
{
	return 1; //Initial version
}
} // packing
} // scoring
} // core
