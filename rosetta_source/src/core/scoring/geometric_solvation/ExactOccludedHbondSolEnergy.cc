// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @brief
/// @author jk

// Project Headers
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergyCreator.hh>

// AUTO-REMOVED #include <core/init.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/Tracer.hh>

#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/hbonds.hh> //ws
// Auto-header: duplicate removed #include <core/scoring/hbonds/types.hh> //ws
#include <core/scoring/hbonds/HBondTypeManager.hh> //ws

// AUTO-REMOVED #include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>

// C++ Headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/hbonds/FadeInterval.hh>
#include <core/scoring/hbonds/hbtrie/HBAtom.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/PyAssert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iosfwd>
#include <limits>
#include <list>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end




//Vector dummy_res_energy_vector_;

static basic::Tracer TR( "core.scoring.geometric_solvation.ExactOccludedHbondSolEnergy" );

namespace core {
namespace scoring {
namespace geometric_solvation {


/// @details This must return a fresh instance of the ExactOccludedHbondSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ExactOccludedHbondSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new ExactOccludedHbondSolEnergy(
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_skip_Hbonders ],
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_include_Hbond_contribution ],
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_pairwise ],
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_pairwise_by_res ],
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_split_between_res ],
		! basic::options::option[ basic::options::OptionKeys::score::exact_occ_self_res_no_occ ],
		basic::options::option[basic::options::OptionKeys::score::exact_occ_radius_scaling]
 );

}

ScoreTypes
ExactOccludedHbondSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( occ_sol_exact );
	return sts;
}


using namespace core;
using namespace core::scoring;
using namespace core::scoring::hbonds;

core::Real const geosol_kT = { 0.593 };
core::Real const max_possible_LK = { -5 };


// apply this weight to everything, so that scale will match LK
core::Real const LK_MATCHING_WEIGHT_EXACT = { 0.387829 };
core::Real const SKIP_HBONDER_CUT = { -0.1 };

GridInfo* GridInfo::instance_( 0 );
WaterWeightGridSet* WaterWeightGridSet::instance_( 0 );

//singleton class
GridInfo * GridInfo::get_instance()
{
	if ( instance_ == 0 )
	{
		 instance_ = new GridInfo();
	}
	return instance_;
}


// private constructor
GridInfo::GridInfo() {

	// Setup dimensions for water grids - use the same dimensions for all grids
	core::Real const water_grid_width = 10.;
	core::Real const water_grid_depth = 8.;
	// For speed use only 52 thousand points per grid - gives identical results...
	xnum_points_ = 41;
	ynum_points_ = 41;
	znum_points_ = 31;
	// Note: below gives 51 million points per grid
		//	xnum_points_ = 401;
		//	ynum_points_ = 401;
		//	znum_points_ = 321;

	xstep_ = water_grid_width / ( xnum_points_ - 1 );
	ystep_ = water_grid_width / ( ynum_points_ - 1 );
	zstep_ = water_grid_depth / ( znum_points_ - 1 );
	// Note: the point at the origin will NOT be considered in calculations - the grid starts AFTER the origin!!
	xorigin_ = -xstep_ * ( xnum_points_ + 1) / 2.;
	yorigin_ = -ystep_ * ( ynum_points_ + 1) / 2.;
	zorigin_ = 0.;

}


//singleton class
WaterWeightGridSet * WaterWeightGridSet::get_instance()
{
	if ( instance_ == 0 )
	{
		 instance_ = new WaterWeightGridSet();
	}
	return instance_;
}

// private constructor
WaterWeightGridSet::WaterWeightGridSet() :
  hb_database_(HBondDatabase::get_database()),
  hbondoptions_( new HBondOptions )
{
	using namespace hbonds;
	hbondoptions_->use_sp2_chi_penalty(false); // override command line settings

	// We need water grids for each donor and acceptor type
	// We could read them in, but let's just compute them from scratch instead
	// We'll store them as a map keyed on hbonds::HBEvalType, with each value a 3D vector of the weights

	TR << "computing and storing water weight grids for acceptor types." << std::endl;
	for(Size i = 1; i <= hbacc_MAX; i++){
		HBEvalType const hbe = HBEval_lookup(hbdon_H2O, HBAccChemType(i), seq_sep_other);
		sum_all_water_weights_[hbe] = fill_water_grid( all_water_weights_[hbe], hbe, *GridInfo::get_instance(), true /*water is donor*/);
	}

	TR << "computing and storing water weight grids for donor types." << std::endl;
	for(Size i = 1; i <= hbdon_MAX; i++){
		HBEvalType const hbe = HBEval_lookup(HBDonChemType(i), hbacc_H2O, seq_sep_other);
		sum_all_water_weights_[hbe] = fill_water_grid( all_water_weights_[hbe], hbe, *GridInfo::get_instance(), false /*water is acceptor*/);
	}
}


// Fill in the water grid
core::Real WaterWeightGridSet::fill_water_grid(
	std::vector < std::vector < std::vector <core::Real> > > & water_weights,
	hbonds::HBEvalType const & hbond_eval_type,
	GridInfo const & grid_info, bool const water_is_donor)
{

	static Vector const base_to_outer(0,0,1);
	core::Real const entropy_scaling = 1.0;
	core::Real const water_O_H_distance = 0.958;

	// Setup grid, initialize to zero
	water_weights.clear();
	water_weights.resize(grid_info.xnum_points());
	for (core::Size tx=0;tx<grid_info.xnum_points();tx++){
		water_weights[tx].resize(grid_info.ynum_points());
		for (core::Size ty=0;ty<grid_info.ynum_points();ty++){
			water_weights[tx][ty].resize(grid_info.znum_points(), 0.);
		}
	}

	// Fill in the water weight grid
	core::Real sum_grid_water_weight = 0.;
	core::Vector water_position(grid_info.xorigin(),grid_info.yorigin(),grid_info.zorigin());

	for (core::Size tx=0;tx<grid_info.xnum_points();tx++){
		water_position.x() += grid_info.xstep();
		water_position.y() = grid_info.yorigin();
		for (core::Size ty=0;ty<grid_info.ynum_points();ty++){
			water_position.y() += grid_info.ystep();
			water_position.z() = grid_info.zorigin();
			for (core::Size tz=0;tz<grid_info.znum_points();tz++){
				water_position.z() += grid_info.zstep();

				// Compute the current geometry
				core::Real AHdis, xD, xH;
				if ( water_is_donor ) {

					// water is the donor, give it perfect geometry
					xD = 0.9999;

					// compute the distance to the accepting water proton
					// subtract the water's OH distance to get the AHdis,
					// since the distance computed was from the acceptor to the water oxygen
					// note: water proton lies on the line between the acceptor and the water oxygen
					AHdis = water_position.length();
					AHdis -= water_O_H_distance; // water O-H distance

					// find cosine of the base-acceptor-water_proton angle (xH)
					// note: this is the same as the base-acceptor-water_oxygen angle
					// note: careful to normalize-by-copy, since water_position is reused around the loop
					xH = dot( base_to_outer, water_position.normalized() );

				} else {

					// water is the acceptor, give it perfect geometry
					xH = 1./3.;  // perfect geometry is cos( 180 - 109.5 degrees), which is 1/3

					// compute the distance to the accepting water
					AHdis = water_position.length();

					// find the cosine of the base-proton-water angle (xD)
					// note: careful to normalize-by-copy, since water_position is reused around the loop
					xD = dot( base_to_outer, water_position.normalized() );

				}

				if ( xH < MIN_xH ) continue;
				if ( xH > MAX_xH ) continue;
				if ( xD < MIN_xD ) continue;
				if ( xD > MAX_xD ) continue;
				if ( AHdis < MIN_R ) continue;
				if ( AHdis > MAX_R ) continue;

				// Get the Hbond energy
				core::Real curr_water_hbond;
				core::Real dummy_chi( 0.0 );
				assert( ! hbondoptions_->use_sp2_chi_penalty() ); // APL avoid the new sp2 chi term.
				hbond_compute_energy(*hb_database_, *hbondoptions_, hbond_eval_type,
					AHdis, xD, xH, dummy_chi, curr_water_hbond );

				// Save the Hbond energy
				curr_water_hbond *= entropy_scaling;
				if ( curr_water_hbond < 0 ) {
					core::Real curr_water_weight = exp( - curr_water_hbond / geosol_kT );
					water_weights[tx][ty][tz] = curr_water_weight;
					sum_grid_water_weight += curr_water_weight;
				}
			}
		}
	}

	return sum_grid_water_weight;
}


std::vector < std::vector < std::vector <core::Real> > > const & WaterWeightGridSet::get_water_weight_grid( hbonds::HBEvalType const & hbond_eval_type ) const {
	// Check that we have weights for this Hbond type
	all_water_weights_iterator curr_water_weights_iter = all_water_weights_.find( hbond_eval_type );
	if ( curr_water_weights_iter == all_water_weights_.end( ) ) {
		TR << "Could not look up map element" << std::endl;
		TR << "get_water_weight_grid hbond_eval_type " << hbond_eval_type << std::endl;
		assert(false);
		exit(1);
	}
	return curr_water_weights_iter->second;
}

core::Real WaterWeightGridSet::get_sum_water_weight_grid( hbonds::HBEvalType const & hbond_eval_type ) const {
	// Check that we have weights for this Hbond type
	sum_water_weights_iterator curr_sum_water_weights_iter = sum_all_water_weights_.find( hbond_eval_type );
	if ( curr_sum_water_weights_iter == sum_all_water_weights_.end( ) ) {
		TR << "Could not look up map element" << std::endl;
		TR << "get_sum_water_weight_grid hbond_eval_type " << hbond_eval_type << std::endl;
		assert(false);
		exit(1);
	}
	return curr_sum_water_weights_iter->second;
}


ExactOccludedHbondSolEnergy::~ExactOccludedHbondSolEnergy() {}

ExactOccludedHbondSolEnergy::ExactOccludedHbondSolEnergy(
	bool const exact_occ_skip_Hbonders,
	bool const exact_occ_include_Hbond_contribution,
	bool const exact_occ_pairwise,
	bool const exact_occ_pairwise_by_res,
	bool const exact_occ_split_between_res,
	bool const exact_occ_self_res_occ,
	core::Real const occ_radius_scaling,
	bool const verbose
) :
	parent( new ExactOccludedHbondSolEnergyCreator ),
	exact_occ_skip_Hbonders_( exact_occ_skip_Hbonders ),
	exact_occ_include_Hbond_contribution_( exact_occ_include_Hbond_contribution ),
	exact_occ_pairwise_( exact_occ_pairwise ),
	exact_occ_pairwise_by_res_( exact_occ_pairwise_by_res ),
	exact_occ_split_between_res_( exact_occ_split_between_res ),
	exact_occ_self_res_occ_( exact_occ_self_res_occ ),
	occ_radius_scaling_( occ_radius_scaling ),
	hb_database_(HBondDatabase::get_database()),
	hbondoptions_(new HBondOptions ),
	verbose_( verbose )
{
	hbondoptions_->use_sp2_chi_penalty( false ); // apl preserve old behavior

	if ( verbose_ ) TR <<"ExactOccludedHbondSolEnergy constructor" << std::endl;
	if ( exact_occ_split_between_res_ && ! exact_occ_pairwise_ ) {
		TR << "Error - cannot split occ energy between residues unless pairwise calculations are used!" << std::endl;
		exit(1);
	}

	// Keep a copy of the AtomTypeSet (to lookup atomic radii)
	atom_type_set_ptr_ = chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD );



	// Allocate memory for grid of occluded sites
	occluded_sites_.clear();
	occluded_sites_.resize(GridInfo::get_instance()->xnum_points());
	for (core::Size tx=0;tx<GridInfo::get_instance()->xnum_points();tx++){
		occluded_sites_[tx].resize(GridInfo::get_instance()->ynum_points());
		for (core::Size ty=0;ty<GridInfo::get_instance()->ynum_points();ty++){
			occluded_sites_[tx][ty].resize(GridInfo::get_instance()->znum_points());
		}
	}

}

ExactOccludedHbondSolEnergy::ExactOccludedHbondSolEnergy( ExactOccludedHbondSolEnergy const & src ):
	parent( src ),
	exact_occ_skip_Hbonders_( src.exact_occ_skip_Hbonders_ ),
	exact_occ_include_Hbond_contribution_( src.exact_occ_include_Hbond_contribution_ ),
	exact_occ_pairwise_( src.exact_occ_pairwise_ ),
	exact_occ_pairwise_by_res_( src.exact_occ_pairwise_by_res_ ),
	exact_occ_split_between_res_( src.exact_occ_split_between_res_ ),
	exact_occ_self_res_occ_( src.exact_occ_self_res_occ_ ),
	occ_radius_scaling_( src.occ_radius_scaling_ ),
  hb_database_(HBondDatabase::get_database()),
	verbose_( src.verbose_ ),
	atom_type_set_ptr_( src.atom_type_set_ptr_ )
{
	if ( verbose_ ) TR <<"ExactOccludedHbondSolEnergy constructor" << std::endl;
	if ( exact_occ_split_between_res_ && ! exact_occ_pairwise_ ) {
		TR << "Error - cannot split occ energy between residues unless pairwise calculations are used!" << std::endl;
		exit(1);
	}

	// Allocate memory for grid of occluded sites
	occluded_sites_.clear();
	occluded_sites_.resize(GridInfo::get_instance()->xnum_points());
	for (core::Size tx=0;tx<GridInfo::get_instance()->xnum_points();tx++){
		occluded_sites_[tx].resize(GridInfo::get_instance()->ynum_points());
		for (core::Size ty=0;ty<GridInfo::get_instance()->ynum_points();ty++){
			occluded_sites_[tx][ty].resize(GridInfo::get_instance()->znum_points());
		}
	}

}

methods::EnergyMethodOP
ExactOccludedHbondSolEnergy::clone() const
{
	return new ExactOccludedHbondSolEnergy( *this );
}

void
ExactOccludedHbondSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

void
ExactOccludedHbondSolEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	pose.update_residue_neighbors();
}

void
ExactOccludedHbondSolEnergy::setup_for_derivatives( pose::Pose & , ScoreFunction const & ) const
{
	TR << "Error - cannot compute derivatives for ExactOccludedHbondSolEnergy (occ_sol_exact)" << std::endl;
	assert(false);
	exit(1);
}

void
ExactOccludedHbondSolEnergy::setup_for_minimizing( pose::Pose & , ScoreFunction const & , kinematics::MinimizerMapBase const & ) const
{
	TR << "Error - cannot compute derivatives for ExactOccludedHbondSolEnergy (occ_sol_exact)" << std::endl;
	assert(false);
	exit(1);
}

Distance ExactOccludedHbondSolEnergy::atomic_interaction_cutoff() const { return 7.5; }

void ExactOccludedHbondSolEnergy::residue_energy(
	conformation::Residue const & polar_rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const {

	core::Size polar_resnum = (core::Size) polar_rsd.seqpos();
	core::Real residue_geosol(0.);

	// loop over donors in polar_rsd
	for ( chemical::AtomIndices::const_iterator
					hnum  = polar_rsd.Hpos_polar().begin(),
					hnume = polar_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const polar_atom( *hnum );
		Size const base_atom( polar_rsd.atom_base( polar_atom ) );
		HBEvalType const curr_hbond_eval_type = HBEval_lookup(get_hb_don_chem_type( base_atom, polar_rsd ), hbacc_H2O, seq_sep_other);

		// Figure out max LK energy
		std::string const base_atom_name = polar_rsd.atom_name( base_atom );
		core::Real max_possible_LK = (*atom_type_set_ptr_)[polar_rsd.atom_type_index(base_atom)].lj_radius();
		if ( ( base_atom_name == " N  " ) && polar_rsd.is_lower_terminus() ) max_possible_LK /= 3; // charged N-terminus
		if ( base_atom_name == " NZ " ) max_possible_LK /= 3; // Lys
		if ( base_atom_name == " ND2" ) max_possible_LK /= 2; // Asn
		if ( base_atom_name == " NE2" ) max_possible_LK /= 2; // Gln
		if ( base_atom_name == " NH1" ) max_possible_LK /= 2; // Arg
		if ( base_atom_name == " NH2" ) max_possible_LK /= 2; // Arg
		// Note: inner nitrogen of Arg (NE) is extra strong, since it's the same atom type as the other two but doesn't get
		// cut in half because there's only one proton...
		//			TR << "jk max LK for donor with base " << base_atom_name << " is  " << max_possible_LK << std::endl;

		// jk INSTEAD OF USING LK dG FREE, SET THEM ALL TO -5.0. THIS IS ALMOST TRUE ANYWAY, AND THE ONES THAT AREN'T SHOULD PROBABLY BE...
		max_possible_LK = -5.;

		// Compute Ebulk (using the LK energy)
		core::Real const Emax_weight = exp( max_possible_LK / geosol_kT );
		core::Real const sum_water_weights = WaterWeightGridSet::get_instance()->get_sum_water_weight_grid( curr_hbond_eval_type );
		core::Real const Ebulk_weight = ( sum_water_weights * Emax_weight ) / ( 1. - Emax_weight);

		// This grid constant is the denominator in computing solvation energies,
		// it depends on the grid dimensions, and sets the max possible solvation energy (in this case to match LK)
		core::Real const grid_constant = sum_water_weights + Ebulk_weight;

		core::Real polar_group_energy = 0.;
		if ( ! exact_occ_pairwise_ && ! exact_occ_pairwise_by_res_ ) {
			polar_group_energy = compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
				grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type ) );
		} else {
			// loop over all atoms of neighboring residues, INCLUDING SELF
			core::scoring::TenANeighborGraph const & graph = pose.energies().tenA_neighbor_graph();
			utility::vector1 <core::Size> neighborlist;
			neighborlist.push_back( polar_resnum);
			for ( core::graph::Graph::EdgeListConstIter
							neighbor_iter = graph.get_node( polar_resnum )->const_edge_list_begin(),
							neighbor_iter_end = graph.get_node( polar_resnum )->const_edge_list_end();
						neighbor_iter != neighbor_iter_end; ++neighbor_iter ) {
				neighborlist.push_back( (*neighbor_iter)->get_other_ind( polar_resnum ) );
			}
			for ( Size occ_inx = 1; occ_inx <= neighborlist.size(); ++occ_inx ) {
				core::Size const occ_resnum( neighborlist[occ_inx] );
				conformation::Residue const occ_rsd = pose.residue(occ_resnum);
				if ( exact_occ_pairwise_by_res_ ) {
					polar_group_energy += compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
						grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type ),
						true, occ_resnum );
				} else {
					for ( Size occ_atomno = 1; occ_atomno <= occ_rsd.natoms(); ++occ_atomno ) {
						polar_group_energy += compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
							grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type ),
							true, occ_resnum, true, occ_atomno );
					}
				}
			}
		}
		residue_geosol += polar_group_energy;

		//solvation energy for donor
		//		std::cout << "jk EXACT Donor " << base_atom_name << "  " << pose.residue(polar_resnum).aa() << " " << polar_resnum << "  " << polar_group_energy << std::endl;

	}

	// loop over acceptors in polar_rsd
	for ( chemical::AtomIndices::const_iterator
					anum  = polar_rsd.accpt_pos().begin(),
					anume = polar_rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const polar_atom( *anum );
		Size const base_atom ( polar_rsd.atom_base( polar_atom ) );
		hbonds::HBEvalType const curr_hbond_eval_type = HBEval_lookup( hbdon_H2O, get_hb_acc_chem_type( polar_atom, polar_rsd ), seq_sep_other);

		// Figure out max LK energy
		std::string const base_atom_name = polar_rsd.atom_name( base_atom );
		core::Real max_possible_LK = (*atom_type_set_ptr_)[ polar_rsd.atom_type_index( polar_atom ) ].lk_dgfree();
		//			TR << "jk max LK for acceptor " << polar_rsd.atom_name(polar_atom) << " is  " << max_possible_LK << std::endl;

		// jk INSTEAD OF USING LK dG FREE, SET THEM ALL TO -5.0. THIS IS ALMOST TRUE ANYWAY, AND THE ONES THAT AREN'T SHOULD PROBABLY BE...
		max_possible_LK = -5.;

		// Compute Ebulk (using the LK energy)
		core::Real const Emax_weight = exp( max_possible_LK / geosol_kT );
		core::Real const sum_water_weights = WaterWeightGridSet::get_instance()->get_sum_water_weight_grid( curr_hbond_eval_type );
		core::Real const Ebulk_weight = ( sum_water_weights * Emax_weight ) / ( 1. - Emax_weight);
		// This grid constant is the denominator in computing solvation energies,
		// it depends on the grid dimensions, and sets the max possible solvation energy (in this case to match LK)
		core::Real const grid_constant = sum_water_weights + Ebulk_weight;

		core::Real polar_group_energy = 0.;
		if ( ! exact_occ_pairwise_ && ! exact_occ_pairwise_by_res_ ) {
			polar_group_energy = compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
				grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type ) );
		} else {
			// loop over all atoms of neighboring residues, INCLUDING SELF
			core::scoring::TenANeighborGraph const & graph = pose.energies().tenA_neighbor_graph();
			utility::vector1 <core::Size> neighborlist;
			neighborlist.push_back( polar_resnum);
			for ( core::graph::Graph::EdgeListConstIter
							neighbor_iter = graph.get_node( polar_resnum )->const_edge_list_begin(),
							neighbor_iter_end = graph.get_node( polar_resnum )->const_edge_list_end();
						neighbor_iter != neighbor_iter_end; ++neighbor_iter ) {
				neighborlist.push_back( (*neighbor_iter)->get_other_ind( polar_resnum ) );
			}
			for ( Size occ_inx = 1; occ_inx <= neighborlist.size(); ++occ_inx ) {
				core::Size const occ_resnum( neighborlist[occ_inx] );
				conformation::Residue const occ_rsd = pose.residue(occ_resnum);
				if ( exact_occ_pairwise_by_res_ ) {
					polar_group_energy += compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
																															 grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type ),
																															 true, occ_resnum );
				} else {
					for ( Size occ_atomno = 1; occ_atomno <= occ_rsd.natoms(); ++occ_atomno ) {
						polar_group_energy += compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
																																 grid_constant,  WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type ),
																																 true, occ_resnum, true, occ_atomno );
					}
				}
			}
		}
		residue_geosol += polar_group_energy;

		//solvation energy for acceptor
		//		std::cout << "jk EXACT Acceptor " << base_atom_name << "  " << pose.residue(polar_resnum).aa() << " " << polar_resnum << "  " << polar_group_energy << std::endl;

	}

	if ( exact_occ_split_between_res_ ) {
		TR << "PAIRWISE OUTPUT FORMAT IS NOT YET SUPPORTED" << std::endl;
		// Here we need to add contributions from this residue occluding all neighboring polar groups then divide everything by two
		// note: this will make the code run twice as slow as it otherwise would, but that's okay because it's just for parameterization / analysis
		exit(1);
	}

	emap[ occ_sol_exact ] += LK_MATCHING_WEIGHT_EXACT * residue_geosol;

}


Real ExactOccludedHbondSolEnergy::compute_polar_group_sol_energy(
	pose::Pose const & pose,
	conformation::Residue const & polar_rsd,
  Size const polar_atom,
	bool const restrict_to_single_occluding_residue, // = false
	Size const single_occluding_resinx, // = 0
	bool const restrict_to_single_occluding_atom, // = false
	Size const single_occluding_atominx // = 0
) const {
	Size const base_atom( polar_rsd.atom_base( polar_atom ) );

	HBEvalType curr_hbond_eval_type;
	if ( polar_rsd.atom_type( base_atom ).is_donor()){
		curr_hbond_eval_type = HBEval_lookup(get_hb_don_chem_type( base_atom, polar_rsd ), hbacc_H2O, seq_sep_other);
	} else if( polar_rsd.atom_type( polar_atom).is_acceptor()){
		curr_hbond_eval_type = HBEval_lookup( hbdon_H2O, get_hb_acc_chem_type( polar_atom, polar_rsd ), seq_sep_other);
	} else {
		assert( false ); // Not a donor or an acceptor, don't know what to do.
	}

	Real const grid_constant = compute_grid_constant( curr_hbond_eval_type );
	return compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
																				grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type ),
																				restrict_to_single_occluding_residue,
																				single_occluding_resinx,
																				restrict_to_single_occluding_atom,
																				single_occluding_atominx);
}


core::Real ExactOccludedHbondSolEnergy::compute_polar_group_sol_energy(
	pose::Pose const & pose,
	conformation::Residue const & polar_rsd,
	core::Size const polar_atomno,
	GridInfo const & grid_info,
	core::Real const & grid_constant,
	std::vector < std::vector < std::vector <core::Real> > > const & water_weights,
	bool const restrict_to_single_occluding_residue, // = false
	core::Size const single_occluding_resinx, // = 0
	bool const restrict_to_single_occluding_atom, // = false
	core::Size const single_occluding_atominx // = 0
) const {

	// note: the "restrict_to_single_occluding_residue / atom" options are so that pairwise additivity can be enforced (for parameterization / analysis)
	core::Size const polar_resnum = (core::Size) polar_rsd.seqpos();

	// jumpout immediately if we're not including self but the single requested residue is self
	if ( restrict_to_single_occluding_residue && ! exact_occ_self_res_occ_ && ( polar_resnum == single_occluding_resinx ) ) {
		return 0;
	}

	// note: the "restrict_to_single_occluding_residue / atom" options are so that pairwise additivity can be enforced (for parameterization / analysis)

	Size const base_atom( polar_rsd.atom_base( polar_atomno ) );
	bool polar_group_is_acceptor =  polar_rsd.atom_type( polar_atomno).is_acceptor();
	bool polar_group_is_donor = polar_rsd.atom_type( base_atom ).is_donor();

	core::Real polar_group_hb_energy(0.);

	bool const hydrogens_can_occlude = false;
	core::Real const water_radius = 1.4;

	// Reset grid of occluded sites, by setting everything to false (for "not occluded")
	for (core::Size tx=0;tx<grid_info.xnum_points();tx++){
		for (core::Size ty=0;ty<grid_info.ynum_points();ty++){
			for (core::Size tz=0;tz<grid_info.znum_points();tz++){
				occluded_sites_[tx][ty][tz] = false;
			}
		}
	}

	// Find the transformation which puts the donor/acceptor of interest onto the existing grid
	// The plan is to apply this transformation to bring each occluding atom onto the existing grid (translation then matrix multiplication)
	core::Size const base_atomno( polar_rsd.atom_base( polar_atomno ) );
	core::Vector const & orig_polar_atom_xyz( polar_rsd.atom( polar_atomno ).xyz() );
	core::Vector const & orig_base_atom_xyz( polar_rsd.atom( base_atomno ).xyz() );
	core::Vector const translation_vector = -1.0 * orig_polar_atom_xyz;
	core::Vector const translated_base_atom = orig_base_atom_xyz + translation_vector;

	// We want to translate positions of occluding atoms _from_ a cartesian basis set into one using the reference frame of the polar group
	core::Vector const cartesian_x(1,0,0);
	core::Vector const cartesian_y(0,1,0);
	core::Vector const cartesian_z(0,0,1);

	// There's not a unique solution for this, so we'll arbitrarily pick a second basis vector, requiring that the dot product with desired_z is zero
	core::Vector desired_z = -1. * translated_base_atom.normalized();

	// Treat the case where polar_atom.z == base_atom.z ; this leads to desired_z.z == 0
	core::Real arbitrary_x, arbitrary_y, arbitrary_z;
	// jk note: we get numerical problems if z-coor is about zero (we divide by it below) - if so then use y to divide by
	if ( std::abs( desired_z.z() ) > 0.01 ) {
		arbitrary_x = 1;
		arbitrary_y = 1;
		arbitrary_z = -1. * ((arbitrary_x * desired_z.x()) + (arbitrary_y * desired_z.y())) / desired_z.z();
	} else {
		if ( std::abs( desired_z.y() ) > 0.01 ) {
			arbitrary_x = 1;
			arbitrary_z = 1;
			arbitrary_y = -1. * ((arbitrary_x * desired_z.x()) + (arbitrary_z * desired_z.z())) / desired_z.y();
		} else {
			// jk note: we get numerical problems if z-coor AND y-coor are about zero (it does happen!) - if so then use x to divide by
			arbitrary_y = 1;
			arbitrary_z = 1;
			arbitrary_x = -1. * ((arbitrary_y * desired_z.y()) + (arbitrary_z * desired_z.z())) / desired_z.x();
		}
	}
	core::Vector desired_x(arbitrary_x,arbitrary_y,arbitrary_z);
	desired_x.normalize();
	core::Vector desired_y = cross_product( desired_x, desired_z );

	// The transformation matrix to do this is i.i'  i.j', etc. where i,j,k are the unit vectors of the starting system
	// and i',j',k' are the unit vectors of the target system
	// for reference see http://kwon3d.com/theory/transform/transform.html
	numeric::xyzMatrix< Length > transformation_matrix;
	transformation_matrix.xx( desired_x.dot(cartesian_x) );
	transformation_matrix.xy( desired_x.dot(cartesian_y) );
	transformation_matrix.xz( desired_x.dot(cartesian_z) );
	transformation_matrix.yx( desired_y.dot(cartesian_x) );
	transformation_matrix.yy( desired_y.dot(cartesian_y) );
	transformation_matrix.yz( desired_y.dot(cartesian_z) );
	transformation_matrix.zx( desired_z.dot(cartesian_x) );
	transformation_matrix.zy( desired_z.dot(cartesian_y) );
	transformation_matrix.zz( desired_z.dot(cartesian_z) );

	// Double-check transformation matrix
	core::Vector new_base_atom_location = transformation_matrix * translated_base_atom;
	assert( std::abs(new_base_atom_location.normalized().x()) < 0.001 );
	assert( std::abs(new_base_atom_location.normalized().y()) < 0.001 );
	assert( std::abs(new_base_atom_location.normalized().z() + 1.) < 0.001 );

	utility::vector1 <core::Size> neighborlist;
	if ( restrict_to_single_occluding_residue ) {
		// consider only a single occluding residue (for pairwise additivity)
		assert ( single_occluding_resinx > 0 );
		neighborlist.push_back( single_occluding_resinx );
	} else {
		// loop over all atoms of neighboring residues, INCLUDING SELF
		core::scoring::TenANeighborGraph const & graph = pose.energies().tenA_neighbor_graph();
		neighborlist.push_back( polar_resnum);
		for ( core::graph::Graph::EdgeListConstIter
						neighbor_iter = graph.get_node( polar_resnum )->const_edge_list_begin(),
						neighbor_iter_end = graph.get_node( polar_resnum )->const_edge_list_end();
					neighbor_iter != neighbor_iter_end; ++neighbor_iter ) {
			neighborlist.push_back( (*neighbor_iter)->get_other_ind( polar_resnum ) );
		}
	}

	for ( Size occ_inx = 1; occ_inx <= neighborlist.size(); ++occ_inx ) {
		core::Size const occ_resnum( neighborlist[occ_inx] );
		if ( ! exact_occ_self_res_occ_ && ( polar_resnum == occ_resnum ) ) {
			continue;
		}
		conformation::Residue const occ_rsd = pose.residue(occ_resnum);

		// loop over all occluding atoms in this residue
		Size atom_startinx = 1;
		Size atom_lastinx = occ_rsd.natoms();
		if ( restrict_to_single_occluding_atom ) {
			// consider only a single occluding atom (for pairwise additivity)
			assert ( single_occluding_atominx > 0 );
			atom_startinx = single_occluding_atominx;
			atom_lastinx = single_occluding_atominx;
		}

		//		TR << "computing occlusion of polar residue " << polar_resnum << " by residue " << occ_resnum << std::endl;

		for ( Size occ_atomno = atom_startinx; occ_atomno <= atom_lastinx; ++occ_atomno ) {

			bool const occ_atom_is_hydrogen = occ_rsd.atom_is_hydrogen( occ_atomno );
			if ( occ_atom_is_hydrogen && ! hydrogens_can_occlude ) continue;

			// can be occluded by atoms directly bonded to this group, but not by self
			if ( polar_resnum == occ_resnum ) {
				if ( polar_atomno == occ_atomno ) continue;
				if ( base_atomno == occ_atomno ) continue;
			}

			core::Real occ_radius = (*atom_type_set_ptr_)[occ_rsd.atom_type_index( occ_atomno )].lj_radius();
			// catch proline NV here (and other virtual atoms, etc.)
			if ( occ_radius < 0.1 ) continue;
			occ_radius *= occ_radius_scaling_;
			core::Real const sq_dist_cut = ( occ_radius + water_radius ) * ( occ_radius + water_radius );

			// if we're not counting contributions from Hbonded atoms, and atom is Hbonded to our polar group, jump out here
			// note: this affects secondary structure a lot!
			// also note: this will not be smoothly differentiable, and will cause problems if ported directly to the fitted function.
			// for the fitted function, we need a better way to do this, eg. downweight the solvation penalty based on how good the Hbond is

			if ( exact_occ_include_Hbond_contribution_ || exact_occ_skip_Hbonders_ ) {
				bool occ_atom_is_Hbonded(false);

				// figure out if the occluding atom is a donor base
				for ( chemical::AtomIndices::const_iterator hnum = occ_rsd.Hpos_polar().begin(), hnume = occ_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
					Size const don_h_atom( *hnum );
					Size const base_atom ( occ_rsd.atom_base( don_h_atom ) );
					if ( occ_atomno == base_atom ) {
						// make sure the polar atom is an acceptor
						for ( chemical::AtomIndices::const_iterator anum = polar_rsd.accpt_pos().begin(), anume = polar_rsd.accpt_pos().end(); anum != anume; ++anum ) {
							Size const acc_atom( *anum );
						 	if ( polar_atomno == acc_atom ) {
								// If so, check if we have an Hbond to the polar group of interest
								HBDonChemType don_chem_type = get_hb_don_chem_type( occ_atomno, occ_rsd );
								HBAccChemType acc_chem_type = get_hb_acc_chem_type( polar_atomno, polar_rsd );
								// Note: this should really test if things are on different chains (comment from Hbond code)
								HBSeqSep seq_sep(get_seq_sep(don_chem_type, acc_chem_type, polar_resnum - occ_resnum));
								HBEvalType const hbe_type = HBEval_lookup(don_chem_type, acc_chem_type, seq_sep);
								Real hb_ener(0.);
								hb_energy_deriv( *hb_database_, *hbondoptions_, hbe_type,
									occ_rsd.atom( occ_atomno ).xyz(), occ_rsd.atom( don_h_atom ).xyz(),
									polar_rsd.atom( polar_atomno ).xyz(), polar_rsd.atom( polar_rsd.atom_base( polar_atomno ) ).xyz(),
									polar_rsd.atom( polar_rsd.abase2( polar_atomno ) ).xyz(), hb_ener);

								if ( hb_ener < 0. ) {
									switch ( get_hbond_weight_type(hbe_type) ) {
									case hbw_SC:
										hb_ener *= 1.1;
										break;
									case hbw_LR_BB:
									case hbw_SR_BB:
										hb_ener *= 1.17;
										break;
									default:
										break;
									}
									polar_group_hb_energy += hb_ener;
								}
								if ( hb_ener < SKIP_HBONDER_CUT ) {
									occ_atom_is_Hbonded = true;
								}
							}
						}
					}
				}
				if ( occ_atom_is_Hbonded && ! exact_occ_include_Hbond_contribution_ ) continue;

				// figure out if the occluding atom is an acceptor
				for ( chemical::AtomIndices::const_iterator anum = occ_rsd.accpt_pos().begin(), anume = occ_rsd.accpt_pos().end(); anum != anume; ++anum ) {
					Size const acc_atom( *anum );
					if ( occ_atomno == acc_atom ) {
						// make sure the polar atom is a donor
						for ( chemical::AtomIndices::const_iterator hnum = polar_rsd.Hpos_polar().begin(), hnume = polar_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
							Size const don_h_atom( *hnum );
							if ( polar_atomno == don_h_atom ) {
								Size const base_atom ( polar_rsd.atom_base( don_h_atom ) );
								// If so, check if we have an Hbond to the polar group of interest
								HBDonChemType don_chem_type = get_hb_don_chem_type( base_atom, polar_rsd );
								HBAccChemType acc_chem_type = get_hb_acc_chem_type( occ_atomno, occ_rsd );
								// Note: this should really test if things are on different chains (comment from Hbond code)
								HBSeqSep seq_sep(get_seq_sep(don_chem_type, acc_chem_type, occ_resnum - polar_resnum ));
								HBEvalType const hbe_type = HBEval_lookup(don_chem_type, acc_chem_type, seq_sep);
								Real hb_ener(0.);
								hb_energy_deriv( *hb_database_, *hbondoptions_, hbe_type, polar_rsd.atom( base_atom ).xyz(), polar_rsd.atom( polar_atomno ).xyz(),
									occ_rsd.atom( occ_atomno ).xyz(), occ_rsd.atom( occ_rsd.atom_base( occ_atomno ) ).xyz(),
									occ_rsd.atom( occ_rsd.abase2( occ_atomno ) ).xyz(), hb_ener);

								if ( hb_ener < 0. ) {
									switch ( get_hbond_weight_type(hbe_type) ) {
									case hbw_SC:
										hb_ener *= 1.1;
										break;
									case hbw_LR_BB:
									case hbw_SR_BB:
										hb_ener *= 1.17;
										break;
									default:
										break;
									}
									polar_group_hb_energy += hb_ener;
								}
								if ( hb_ener < SKIP_HBONDER_CUT ) {
									occ_atom_is_Hbonded = true;
								}
							}
						}
					}
				}
				if ( occ_atom_is_Hbonded && ! exact_occ_include_Hbond_contribution_ ) continue;

			} // done finding Hbonds for skip_Hbonders and/or including Hbond contribution

			// Apply the transformation to put this atom onto the current grid
			core::Vector const & orig_occ_atom_xyz( occ_rsd.atom( occ_atomno ).xyz() );
			core::Vector const translated_occ_atom_xyz = orig_occ_atom_xyz + translation_vector;
			core::Vector const transformed_occ_atom_xyz = transformation_matrix * ( orig_occ_atom_xyz + translation_vector );

			// Double-check transformations
			assert( std::abs(orig_polar_atom_xyz.distance( orig_occ_atom_xyz ) - transformed_occ_atom_xyz.magnitude()) < 0.001 );
			assert( std::abs(orig_base_atom_xyz.distance( orig_occ_atom_xyz ) - new_base_atom_location.distance( transformed_occ_atom_xyz )) < 0.001 );

			// Loop over all water positions, mark those which are occluded by this atom
			core::Vector water_position(grid_info.xorigin(),grid_info.yorigin(),grid_info.zorigin());
			for (core::Size wx=0;wx<grid_info.xnum_points();wx++){
				water_position.x() += grid_info.xstep();
				core::Real sq_xdist = ( water_position.x() - transformed_occ_atom_xyz.x() ) * ( water_position.x() - transformed_occ_atom_xyz.x() );
				if ( sq_xdist > sq_dist_cut ) continue;
				water_position.y() = grid_info.yorigin();
				for (core::Size wy=0;wy<grid_info.ynum_points();wy++){
					water_position.y() += grid_info.ystep();
					core::Real sq_ydist = ( water_position.y() - transformed_occ_atom_xyz.y() ) * ( water_position.y() - transformed_occ_atom_xyz.y() );
					if ( sq_ydist > sq_dist_cut ) continue;
					water_position.z() = grid_info.zorigin();
					for (core::Size wz=0;wz<grid_info.znum_points();wz++){
						water_position.z() += grid_info.zstep();
						core::Real sq_zdist = ( water_position.z() - transformed_occ_atom_xyz.z() ) * ( water_position.z() - transformed_occ_atom_xyz.z() );
						if ( sq_zdist > sq_dist_cut ) continue;
						core::Real sq_curr_dist = sq_xdist + sq_ydist + sq_zdist;
						if ( sq_curr_dist < sq_dist_cut ) {
							// this atom occludes this water site
							occluded_sites_[wx][wy][wz] = true;
						}
					}
				}
			}

		}
	}

	// Compute and store the solvation energy for the grid occluded by all nearby atoms
	// Compute the numerator (the sum of occluded weights)
	core::Real sum_occluded_weights(0.);
	for (core::Size tx=0;tx<grid_info.xnum_points();tx++){
		for (core::Size ty=0;ty<grid_info.ynum_points();ty++){
			for (core::Size tz=0;tz<grid_info.znum_points();tz++){
				if ( occluded_sites_[tx][ty][tz] ) {
					core::Real const curr_water_weight = water_weights[tx][ty][tz];
					sum_occluded_weights += curr_water_weight;
				}
			}
		}
	}

	// Compute the energetic cost of occluding this polar group
	core::Real geometric_solvation_energy = - geosol_kT * log( 1 - ( sum_occluded_weights / grid_constant ) );
	core::Real desired_hb_weight = 0.;

	if ( exact_occ_include_Hbond_contribution_ ) {

		if ( exact_occ_self_res_occ_ ) {

				if ( polar_group_is_acceptor ) {

					HBAccChemType acc_chem_type = get_hb_acc_chem_type( polar_atomno, polar_rsd );

					switch( acc_chem_type) {
					case hbacc_PBA:
						geometric_solvation_energy -= 0.775088;
						desired_hb_weight = geometric_solvation_energy *0.72908;
						break;
					case hbacc_CXA:
						geometric_solvation_energy -= 0.775088;
						desired_hb_weight = geometric_solvation_energy *0.94358;
						break;
					case hbacc_CXL:
						geometric_solvation_energy -=0.32359;
						desired_hb_weight = geometric_solvation_energy *0.83315;
						break;
					case hbacc_IMD:
						geometric_solvation_energy -= 1.07342;
						desired_hb_weight = geometric_solvation_energy *0.64448;
						break;
					case hbacc_IME:
						geometric_solvation_energy -= 1.67093;
						desired_hb_weight = geometric_solvation_energy *1.0758;
						break;
					case hbacc_AHX:
						geometric_solvation_energy -= 0.681323;
						desired_hb_weight = geometric_solvation_energy *0.97637;
						break;
					case hbacc_HXL:
						geometric_solvation_energy -= 0.752036;
						desired_hb_weight = geometric_solvation_energy *1.2094;
						break;
					default:
						std::cout << "Chemical acceptor type not found: " << acc_chem_type << std::endl;
						exit(1);
						break;
					}
				}
				if (polar_group_is_donor ) {

					HBDonChemType don_chem_type = get_hb_don_chem_type( base_atom, polar_rsd );

					switch (don_chem_type ) {
					case hbdon_PBA:
						geometric_solvation_energy -= 1.56141;
						desired_hb_weight = geometric_solvation_energy *0.59343;
						break;
					case hbdon_CXA:
						geometric_solvation_energy -= 0.4999125;
						desired_hb_weight = geometric_solvation_energy*.8906;
						break;
					case hbdon_IMD:
						geometric_solvation_energy -=0.957797;
						desired_hb_weight = geometric_solvation_energy *0.90386;
						break;
					case hbdon_IME:
						geometric_solvation_energy -= 0.591156;
						desired_hb_weight = geometric_solvation_energy *0.9971;
						break;
					case hbdon_IND:
						geometric_solvation_energy -=0.854121;
						desired_hb_weight = geometric_solvation_energy *0.97701;
						break;
					case hbdon_AMO:
						geometric_solvation_energy -= 0.33501;
						desired_hb_weight = geometric_solvation_energy*1.5641;
						break;
					case hbdon_GDE:
						geometric_solvation_energy -= 0.981806;
						desired_hb_weight = geometric_solvation_energy *0.8612;
						break;
					case hbdon_GDH:
						geometric_solvation_energy -=0.538248;
						desired_hb_weight = geometric_solvation_energy *1.142;
						break;
					case hbdon_AHX:
						geometric_solvation_energy -= 0.681433;
						desired_hb_weight = geometric_solvation_energy *1.1076;
						break;
					case hbdon_HXL:
						geometric_solvation_energy -= 0.785331;
						desired_hb_weight = geometric_solvation_energy *1.1066;
						break;
					default:
						std::cout << "Chemical donor type not found: " << don_chem_type << std::endl;
						exit(1);
						break;
					}
				}

		} else {

			// constant shift (brings non-Hbonded groups to zero), will match Hbonded groups at about zero
			if ( geometric_solvation_energy > 0 ) {
				if ( polar_group_is_acceptor ) {

					HBAccChemType acc_chem_type = get_hb_acc_chem_type( polar_atomno, polar_rsd );

					switch( acc_chem_type) {
					case hbacc_PBA:
						geometric_solvation_energy -= 0.624376;
						desired_hb_weight = geometric_solvation_energy *0.72832;
						break;
					case hbacc_CXA:
						geometric_solvation_energy -= 0.427747;
						desired_hb_weight = geometric_solvation_energy *0.89874;
						break;
					case hbacc_CXL:
						geometric_solvation_energy -= 0.286873;
						desired_hb_weight = geometric_solvation_energy *0.82413;
						break;
					case hbacc_IMD:
						geometric_solvation_energy -= 0.725837;
						desired_hb_weight = geometric_solvation_energy *0.85228;
						break;
					case hbacc_IME:
						geometric_solvation_energy -= 0.690253;
						desired_hb_weight = geometric_solvation_energy *0.91419;
						break;
					case hbacc_AHX:
						geometric_solvation_energy -= 0.64312;
						desired_hb_weight = geometric_solvation_energy *1.0239;
						break;
					case hbacc_HXL:
						geometric_solvation_energy -= 0.713553;
						desired_hb_weight = geometric_solvation_energy *1.0076;
						break;
					default:
						std::cout << "Chemical acceptor type not found: " << acc_chem_type << std::endl;
						exit(1);
						break;
					}
				}
				if (polar_group_is_donor ) {

					HBDonChemType don_chem_type = get_hb_don_chem_type( base_atom, polar_rsd );

					switch (don_chem_type ) {
					case hbdon_PBA:
						geometric_solvation_energy -= 0.794341;
						desired_hb_weight = geometric_solvation_energy *0.606;
						break;
					case hbdon_CXA:
						geometric_solvation_energy -= 0.294998;
						desired_hb_weight = geometric_solvation_energy*1.0568;
						break;
					case hbdon_IMD:
						geometric_solvation_energy -= 0.426234;
						desired_hb_weight = geometric_solvation_energy *0.82543;
						break;
					case hbdon_IME:
						geometric_solvation_energy -= 0.297841;
						desired_hb_weight = geometric_solvation_energy *1.1028;
						break;
					case hbdon_IND:
						geometric_solvation_energy -= 0.5463725;
						desired_hb_weight = geometric_solvation_energy *0.98083;
						break;
					case hbdon_AMO:
						geometric_solvation_energy -= 0.1029515;
						desired_hb_weight = geometric_solvation_energy*1.727;
						break;
					case hbdon_GDE:
						geometric_solvation_energy -= 0.2953195;
						desired_hb_weight = geometric_solvation_energy *0.83017;
						break;
					case hbdon_GDH:
						geometric_solvation_energy -= 0.22389;
						desired_hb_weight = geometric_solvation_energy *1.5544;
						break;
					case hbdon_AHX:
						geometric_solvation_energy -= 0.494319;
						desired_hb_weight = geometric_solvation_energy *0.98647;
						break;
					case hbdon_HXL:
						geometric_solvation_energy -= 0.451273;
						desired_hb_weight = geometric_solvation_energy *1.1331;
						break;
					default:
						std::cout << "Chemical donor type not found: " << don_chem_type << std::endl;
						exit(1);
						break;
					}
				}
			}

		}

	}

	geometric_solvation_energy = geometric_solvation_energy + (desired_hb_weight*polar_group_hb_energy );
	return geometric_solvation_energy;
}


core::Real ExactOccludedHbondSolEnergy::compute_grid_constant( HBEvalType const & hbond_eval_type ) const
{

	//		// Figure out max LK energy for DONOR
	//		std::string const base_atom_name = polar_rsd.atom_name( base_atom );
	//		core::Real max_possible_LK = etable_ptr_->lk_dgfree( polar_rsd.atom_type_index( base_atom ) );
	//		if ( ( base_atom_name == " N  " ) && polar_rsd.is_lower_terminus() ) max_possible_LK /= 3; // charged N-terminus
	//		if ( base_atom_name == " NZ " ) max_possible_LK /= 3; // Lys
	//		if ( base_atom_name == " ND2" ) max_possible_LK /= 2; // Asn
	//		if ( base_atom_name == " NE2" ) max_possible_LK /= 2; // Gln
	//		if ( base_atom_name == " NH1" ) max_possible_LK /= 2; // Arg
	//		if ( base_atom_name == " NH2" ) max_possible_LK /= 2; // Arg
	//		// Note: inner nitrogen of Arg (NE) is extra strong, since it's the same atom type as the other two but doesn't get
	//		// cut in half because there's only one proton...
	//		//			TR << "jk max LK for donor with base " << base_atom_name << " is  " << max_possible_LK << std::endl;
	//
	//    // --OR--
	//
	//		// Figure out max LK energy for ACCEPTOR
	//		std::string const base_atom_name = polar_rsd.atom_name( base_atom );
	//		core::Real max_possible_LK = etable_ptr_->lk_dgfree( polar_rsd.atom_type_index( polar_atom ) );
	//		//			TR << "jk max LK for acceptor " << polar_rsd.atom_name(polar_atom) << " is  " << max_possible_LK << std::endl;
	//
	//
	//		// jk INSTEAD OF USING LK dG FREE, SET THEM ALL TO -5.0. THIS IS ALMOST TRUE ANYWAY, AND THE ONES THAT AREN'T SHOULD PROBABLY BE...
	//		max_possible_LK = -5.;  // -> mjo moved magic number up to namespace level
	//
	//
	//    jk LK dG FREE is the energy difference between a completely exposed
	//    jk and a completely buried polar atom (in the LK model). They get
	//    jk these numbers from MD simulations in explicit water. The
	//    jk reference would be the original L+K paper (PMID 10223287), the
	//    jk numbers themselves are in Table 1 though in Rosetta we use
	//    jk slightly different numbers.
	//
	//	  jk I was originally using these numbers (since I too need to set
	//	  jk the max possible solvation energy for a given functional group),
	//	  jk but then realized that the LK numbers (or at least, their
	//	  jk Rosetta incarnation) are all basically 5 divided by the number
	//	  jk of protons / lone pairs (which I don't think L+K appreciated),
	//	  jk so I just started using this instead of LK dG FREE.
	//

	// Compute Ebulk (using the LK energy)
	core::Real const Emax_weight = exp( max_possible_LK / geosol_kT );
	core::Real const sum_water_weights = WaterWeightGridSet::get_instance()->get_sum_water_weight_grid( hbond_eval_type );
	core::Real const Ebulk_weight = ( sum_water_weights * Emax_weight ) / ( 1. - Emax_weight);
	// This grid constant is the denominator in computing solvation energies,
	// it depends on the grid dimensions, and sets the max possible solvation energy (in this case to match LK)
	core::Real const grid_constant = sum_water_weights + Ebulk_weight;

	return grid_constant;
}
core::Size
ExactOccludedHbondSolEnergy::version() const
{
	return 1; // Initial versioning
}



} // geometric_solvation
} // scoring
} // core
