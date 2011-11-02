// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/SecondaryStructurePotential.cc
/// @brief  Scoring manager class header
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Phil Bradley

// Unit header
#include <core/scoring/SecondaryStructurePotential.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <basic/database/open.hh>
#include <core/graph/DisjointSets.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/Energies.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

// numeric headers
#include <numeric/conversions.hh>
#include <numeric/trig.functions.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <basic/prof.hh>
// ObjexxFCL headers
//#include <ObjexxFCL/ObjexxFCL.hh>
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// C++ headers
#include <cmath>
#include <iostream>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
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
#include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/conformation/symmetry/SymDof.fwd.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
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
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/SS_Info.hh>
#include <core/scoring/SS_Killhairpins_Info.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/SecondaryStructureWeights.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/Methods.hh>
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
#include <utility/file/gzip_util.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/irstream.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/zipstream.hpp>
#include <utility/io/zipstream.ipp>
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
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/constants.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
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
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray4.fwd.hh>
#include <ObjexxFCL/FArray4.hh>
#include <ObjexxFCL/FArray4D.fwd.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
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
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <boost/pool/poolfwd.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end


namespace core {
namespace scoring {

/// @brief default constructor
SecondaryStructurePotential::SecondaryStructurePotential():
	iptsn_( 36 ),
	pts_( 2, 3, 36, 36 ),
	ds_( 6 ),
	idsn_( 6, idsn_initializer ),
	ids_( 6, ids_initializer ),
	ssdist_( 4, 2, ssdist_initializer ),
	hs_dp_( 10, hs_dp_initializer ),
	rsigma_dot_( 12, 18, 2, 2, rsigma_dot_initializer ),
	m_term_( 4, m_term_initializer )
{
	load_phi_theta_bins();
}


/// helper function
SS_Info const &
retrieve_const_ss_info_from_pose( pose::Pose const & pose )
{
	// ////using core::pose::datacache::CacheableDataType::SS_INFO;
	assert( pose.data().has( core::pose::datacache::CacheableDataType::SS_INFO ) );
	assert( dynamic_cast< SS_Info const *>( &( pose.data().get( core::pose::datacache::CacheableDataType::SS_INFO ))));
	return ( static_cast< SS_Info const &>(    pose.data().get( core::pose::datacache::CacheableDataType::SS_INFO )));
}

/// helper function
SS_Info &
retrieve_nonconst_ss_info_from_pose( pose::Pose & pose )
{
	// ////using core::pose::datacache::CacheableDataType::SS_INFO;

	if ( !pose.data().has( core::pose::datacache::CacheableDataType::SS_INFO ) ) {
		// create new one
		pose.data().set( core::pose::datacache::CacheableDataType::SS_INFO, new SS_Info() );
	}
	assert( pose.data().has( core::pose::datacache::CacheableDataType::SS_INFO ) );
	assert( dynamic_cast< SS_Info *>( &( pose.data().get( core::pose::datacache::CacheableDataType::SS_INFO ))));
	return ( static_cast< SS_Info &>(    pose.data().get( core::pose::datacache::CacheableDataType::SS_INFO )));
}

/// helper function
void
fill_bb_pos(
	pose::Pose const & pose,
	BB_Pos & bb_pos
)
{
	bb_pos.take_coordinates_from_pose( pose );
}


///
void
SecondaryStructurePotential::setup_for_scoring( pose::Pose & pose ) const
{
	basic::ProfileThis doit( basic::SECONDARY_STRUCTURE_ENERGY );

	SS_Info & ss_info( retrieve_nonconst_ss_info_from_pose( pose ) );

	// set the size (does nothing if already correct size)
	// note that this does NOT clear the results of the last calculation
	// that must be done by the individual routines using the cached data
	//
	ss_info.resize( pose.total_residue() );

	// fill the backbone position array for fast access to coords during calculation
	fill_bb_pos( pose, ss_info.bb_pos() );

	// identify strand and helix segments
	identify_ss( pose, ss_info.helices(), ss_info.strands() );

	//	std::cout << "identify_ss:\n" << ss_info.helices() << ss_info.strands();

}




/// @brief score secondary structure
void
SecondaryStructurePotential::score(
	pose::Pose const & pose,
	SecondaryStructureWeights const & weights,
	Real & hs_score,
	Real & ss_score,
	Real & rsigma_score,
	Real & sheet_score
) const
{
	basic::ProfileThis doit( basic::SECONDARY_STRUCTURE_ENERGY );
	ss_score = 0.0;
	hs_score = 0.0;
	rsigma_score = 0.0;

	hspair( pose, hs_score );

	sspair( pose, weights, ss_score, rsigma_score );

	sheet_score = sheets_from_dimers( pose ); // execute after scoring sspair
}


/// @brief score hspair
void
SecondaryStructurePotential::hspair(
	pose::Pose const & pose,
	Real & hs_score
) const
{
	basic::ProfileThis doit( basic::SECONDARY_STRUCTURE_HSPAIR_ENERGY );
	hs_score = 0.0;

	// retrieve the cached info on SSE's from the pose
	SS_Info const & ss_info( retrieve_const_ss_info_from_pose( pose ) );

	Helices const & helices( ss_info.helices() );
	Strands const & strands( ss_info.strands() );
	BB_Pos const & bb_pos( ss_info.bb_pos() );

	if ( strands.total_SS_dimer < 1 || helices.total_HH_dimer < 1 ) return; // helix-strand interactions

	// are we symmetric?
	bool symmetric=false;
	core::conformation::symmetry::SymmetryInfoOP symm_info = NULL;
	if (core::pose::symmetry::is_symmetric( pose ) ) {
		symmetric=true;
		SymmetricConformation const & SymmConf (
		    dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info()->clone();
	}

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( int ss1 = 1; ss1 <= helices.total_HH_dimer; ++ss1 ) {
		int const HH_resnum_ss1 = helices.HH_resnum(ss1);
		int const HH_helix_end_1ss1 = helices.HH_helix_end(1,ss1);
		int const HH_helix_end_2ss1 = helices.HH_helix_end(2,ss1);

		Vector pt1, pt2;
		helix_end( HH_resnum_ss1, bb_pos, pt1, pt2 );

		// loop over ALL the 3D neighbors of seqpos= HH_resnum_ss1 (ie not just standard loop over UPPER nbrs)
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node( HH_resnum_ss1 )->const_edge_list_begin(),
				irue = energy_graph.get_node( HH_resnum_ss1 )->const_edge_list_end();
				iru != irue; ++iru ) {

			int SS_resnum_ss2( (*iru)->get_second_node_ind() );

			//Edges always have first node < second node. Just in case we picked the wrong one:
			if (SS_resnum_ss2 == HH_resnum_ss1) SS_resnum_ss2 = (*iru)->get_first_node_ind();

			int const ss2( strands.SS_dimer( SS_resnum_ss2 ) );
			if ( ss2 == 0 ) continue;

			Vector const & pt3( bb_pos.N( SS_resnum_ss2   ) );
			Vector const & pt4( bb_pos.C( SS_resnum_ss2+1 ) );

			Vector vdist, cen1, cen2;
			Real dist;
			dist_pair(pt1,pt2,pt3,pt4,dist,cen1,cen2,vdist);

			if ( dist <= 12.0 ) {
				Real ph, th;
				spherical(pt2,pt4,ph,th,cen1,cen2,vdist);

				if ( ph < -180.0 ) {
					ph += 360.0f;
				} else if ( ph > 180.0 ) {
					ph -= 360.0f;
				}

				int iph = static_cast< int >(1+(ph+180.0f)/10);
				if ( iph > 36 ) {
					iph = 36;
				} else if ( iph < 1 ) {
					iph = 1;
				}

				int ith = static_cast< int >(1+(th/5));
				if ( ith > 36 ) {
					ith = 36;
				} else if ( ith < 1 ) {
					ith = 1;
				}

				int iseqsep;
				int itemp = std::min( get_foldtree_seqsep( pose, strands.SS_strand_end(2,ss2), HH_helix_end_1ss1 ) + 1,
															get_foldtree_seqsep( pose, HH_helix_end_2ss1, strands.SS_strand_end(1,ss2) ) + 1 );
//				int itemp = std::min( std::abs( strands.SS_strand_end(2,ss2) - HH_helix_end_1ss1 ) + 1,
//															 std::abs( HH_helix_end_2ss1 - strands.SS_strand_end(1,ss2) ) + 1 );
				if ( itemp >= 2 && itemp <= 10 ) {
					iseqsep = 2;
				} else {
					if ( itemp > 10 ) {
						iseqsep = 3;
					} else {
						iseqsep = 1;
					}
				}

				//				std::cout << "hs_intxn: " << HH_resnum_ss1 << ' ' << SS_resnum_ss2 << ' ' << pts_(1,iseqsep,iph,ith) <<
				//					"      " << cendist << " " <<
				//					std::endl;

				if (symmetric) {
					hs_score += pts_(1,iseqsep,iph,ith)*symm_info->score_multiply(HH_resnum_ss1,SS_resnum_ss2);
				} else {
					hs_score += pts_(1,iseqsep,iph,ith);
				}
			}
		} // loop over neighbors of HH_resnum_ss1
	} // ss1


	// modify by proper weighting
	hs_score *= 0.090;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
class DimerPairing : public utility::pointer::ReferenceCount {
public:
	DimerPairing(
		int const dimer1_in,
		int const dimer2_in,
		int const sign1_in,
		int const sign2_in,
		Real const score_in
	): dimer1_( dimer1_in ),
		 dimer2_( dimer2_in ),
		 sign1_( sign1_in ),
		 sign2_( sign2_in ),
		 score_( score_in ),
		 valid_( true )
	{}

	///
	int
	dimer1() const
	{
		return dimer1_;
	}

	///
	int
	dimer2() const
	{
		return dimer2_;
	}

	///
	int
	sign1() const
	{
		return sign1_;
	}

	///
	int
	sign2() const
	{
		return sign2_;
	}

	///
	Real
	score() const
	{
		return score_;
	}

	///
	bool
	valid() const
	{
		return valid_;
	}

	///
	void
	valid( bool const setting )
	{
		valid_ = setting;
	}

private:

	int dimer1_;
	int dimer2_;
	int sign1_;
	int sign2_;

	Real score_;

	bool valid_; // set to false if it conflicts with higher-scoring dimer

};

typedef utility::pointer::owning_ptr< DimerPairing > DimerPairingOP;

bool
dimer_pairing_pointer_sorter( DimerPairingOP const & a, DimerPairingOP const & b )
{
	return ( a->score() < b->score() );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief score sspair

void
SecondaryStructurePotential::sspair(
	pose::Pose const & pose,
	SecondaryStructureWeights const & wts,
	Real & ss_score,
	Real & rsigma_score
) const
{
	basic::ProfileThis doit( basic::SECONDARY_STRUCTURE_SSPAIR_ENERGY );
	 ////using core::pose::datacache::CacheableDataType::core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO;
	ss_score = 0.0;
	rsigma_score = 0.0;

	//
	int const lowstrand( wts.get_ss_lowstrand() );
	int const cutoff( wts.get_ss_cutoff() );

	SS_Info const & ss_info( retrieve_const_ss_info_from_pose( pose ) );

	Strands const & strands( ss_info.strands() );
	BB_Pos const & bb_pos( ss_info.bb_pos() );

	if ( strands.total_SS_dimer < 1 ) return;

	// We need to initialize some things for symmetry scoring. Default we are assymetric
	bool symmetric=false;
	core::conformation::symmetry::SymmetryInfoOP symm_info = NULL;
	if (core::pose::symmetry::is_symmetric( pose ) ) {
		symmetric=true;
		SymmetricConformation const & SymmConf (
		    dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info()->clone();
	}

	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	// new plan: keep a list of dimer pairs with good interaction energy
	// sort by energy after compiling

	typedef utility::vector1< DimerPairingOP > DimerPairings;
	DimerPairings dimer_pairs;

	 // dimer to dimer score (<6.5A)
	static FArray1D_real const SS_penalty( 11, ss_penalty_initializer );

//car local
	Real const ssdist_12 = ssdist_(1,2);

	//car initialize
	// these are used in the sheet score (?)
	for ( int ss1 = 1; ss1 <= strands.total_SS_dimer; ++ss1 ) {
		strands.dimer_neighbor(1,ss1) = 0;
		strands.dimer_neighbor(2,ss1) = 0;
	}



//car ss1 is the first dimer in the possible pair
	for ( int ss1 = 1; ss1 < strands.total_SS_dimer; ++ss1 ) {

		int const SS_resnum_ss1 = strands.SS_resnum(ss1);
		int const SS_strand_end_1ss1 = strands.SS_strand_end(1,ss1);
		int const SS_strand_end_2ss1 = strands.SS_strand_end(2,ss1);
		int & dimer_neighbor_1ss1( strands.dimer_neighbor(1,ss1) );
		int & dimer_neighbor_2ss1( strands.dimer_neighbor(2,ss1) );
		int const SS_strand_ss1 = strands.SS_strand(ss1);

//  SS_resnum(ss1) IS THE SEQUENCE POSITION OF STRAND ss1
		Vector const & pt1( bb_pos.N( SS_resnum_ss1   ) );
		Vector const & pt2( bb_pos.C( SS_resnum_ss1+1 ) );

//car ss2 is second dimer, always C term to first
		// loop over the UPPER neighbors of seqpos= SS_resnum_ss1
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node( SS_resnum_ss1 )->const_upper_edge_list_begin(),
				irue = energy_graph.get_node( SS_resnum_ss1 )->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			int const SS_resnum_ss2( (*iru)->get_second_node_ind() );
			int const ss2( strands.SS_dimer( SS_resnum_ss2 ) );

			if ( ss2 == 0 ) continue; // skip if this position isnt 1st residue of a dimer

//  IF TWO DIMERS HAVE A GOOD PHI/THETA AND DOT PRODUCT, ss_orient=1
			int ss_orient( 0 );

//car calculate the sequence separation between the two dimers
//car note that this is not the strand separation that phitheta is
//car  conditioned on;  the dimer separation is used to decide which
//car pair to count in ss_score based on the value of cutoff

			int const dimer_seqsep = get_foldtree_seqsep( pose, SS_resnum_ss2, SS_resnum_ss1 );
//			int const dimer_seqsep = SS_resnum_ss2 - SS_resnum_ss1;

			int const SS_strand_end_1ss2 = strands.SS_strand_end(1,ss2);

			if ( SS_strand_end_1ss1 != SS_strand_end_1ss2 ) { // not in same strand
				Vector const pt3( bb_pos.N( SS_resnum_ss2   ) );
				Vector const pt4( bb_pos.C( SS_resnum_ss2+1 ) );

//car find vector and distance between midpoints
				Vector cen1, cen2, vdist;
				Real dist;
				dist_pair(pt1,pt2,pt3,pt4,dist,cen1,cen2,vdist);

				Real strand_dist_cutoff = wts.get_strand_dist_cutoff(); // Default is 6.5 Angstroms

				bool const stretch_strand_dist_cutoff = wts.get_stretch_strand_dist_cutoff();
				if ( stretch_strand_dist_cutoff ) {
					Real const seq_sep_scale  = wts.get_seq_sep_scale();
					Real const max_strand_dist_cutoff  = wts.get_max_strand_dist_cutoff();
					strand_dist_cutoff = 6.5 +
						(max_strand_dist_cutoff - 6.5) *(dimer_seqsep / seq_sep_scale);
					if (strand_dist_cutoff > max_strand_dist_cutoff) strand_dist_cutoff = max_strand_dist_cutoff;
				}

				//				std::cout << "dist " << ss1 << ' ' << ss2 << ' ' << dist << ' ' << strand_dist_cutoff << std::endl;

				if ( dist < strand_dist_cutoff ) {
//car find phi and theta
					Real phi_ss, theta;
					spherical(pt2,pt4,phi_ss,theta,cen1,cen2,vdist);

//car find sequence separation between STRANDS  (not between dimers)
					int const strand_seqsep = get_foldtree_seqsep( pose, SS_strand_end_1ss2, SS_strand_end_2ss1 ) + 1;
//					int const strand_seqsep = SS_strand_end_1ss2 - SS_strand_end_2ss1 + 1;

//car bin all these values...
					int istrand_seqsep;
					if ( strand_seqsep >= 2 && strand_seqsep <= 10 ) {
						istrand_seqsep = 2;
					} else {
						if ( strand_seqsep > 10 ) {
							istrand_seqsep = 3;
						} else {
							istrand_seqsep = 1;
						}
					}

					if ( phi_ss > 180.0 ) {
						phi_ss -= 360.0f;
					} else if ( phi_ss < -180.0 ) {
						phi_ss += 360.0f;
					}

					int iphi = static_cast< int >(1+(phi_ss+180.0)/10);
					if ( iphi > 36 ) {
						iphi = 36;
					} else if ( iphi < 1 ) {
						iphi = 1;
					}

					int itheta = static_cast< int >(1+(theta/5));
					if ( itheta > 36 ) {
						itheta = 36;
					} else if ( itheta < 1 ) {
						itheta = 1;
					}

//car find dp_all and sign1,sign2  (convert vdist to unit vec first)
					Vector const unit_vdist( vdist.normalized() );

					Real dpall;
					int sign1, sign2;
					pair_dp( SS_resnum_ss1, SS_resnum_ss2, bb_pos, dpall, unit_vdist, sign1, sign2 );

//car evaluate the first two scoring terms:
					Real phithetascore = pts_(2,istrand_seqsep,iphi,itheta);
					Real distscore = 0.0;

//rhiju Extra rewards/bonuses for parallel and antiparallel terms. Only affects
//         long-range pairings.
					if (theta<90 & istrand_seqsep==3) phithetascore *= wts.get_parallel_weight();
					if (theta>90 & istrand_seqsep==3) phithetascore *= wts.get_antiparallel_weight();

//car save the total score for the pair and the signs
//car dimer_pair_score is the total for the dimer pair
//     LOCAL STRANDS DON'T GET SO MUCH OF A SCORE BONUS
//car add terms to pair only if pair separation > cutoff
//Objexx: Assumes these arrays have same dimensions as ss_orient
					Real dimer_pair_score( 0.0 );

					if ( dimer_seqsep >= cutoff ) {
						if ( lowstrand > 0.5 ) distscore = ssdist_12;
						dimer_pair_score += phithetascore + distscore;
					}

					//					std::cout << "pts: " << ss1 << ' ' << ss2 << ' ' << phithetascore << ' ' << cutoff << std::endl;

					if ( phithetascore < 0.0 ) {
//car bin dpall to get idot
//km changed bins 5 and 6
						int idot;
						if ( dpall > 0.33 && dpall <= 0.66 ) {
							idot = 2;
						} else if ( dpall > 0.66 && dpall <= 1.00 ) {
							idot = 3;
						} else if ( dpall > 1.00 && dpall <= 1.33 ) {
							idot = 4;
						} else if ( dpall > 1.33 && dpall <= 1.60 ) {
							idot = 5;
						} else if ( dpall > 1.60 && dpall <= 1.80 ) {
							idot = 6;
						} else if ( dpall > 1.80 && dpall <= 2.00 ) {
							idot = 5;
						} else {
							idot = 1;
						}

//car lookup the dotscore
						Real dotscore = ds_( idot );

//car again, if the distance between pairs is too small, don't add this
//car term to the total
//     LOCAL STRANDS DON'T GET SO MUCH OF A SCORE BONUS
						if ( dimer_seqsep >= cutoff ) dimer_pair_score += dotscore; // dimer_pair_score(ss1,ss2)

//car note that ss_orient is 1 if phithetascore<0 and dotscore<0 and dist<6.5
//car but these terms not in dimer_pair_score if dimer_seqsep < cutoff
						if ( dotscore < 0.0 ) ss_orient = 1; // ss_orient(ss1,ss2)

//js quick fix for rsigma dependence on sign1,sign2 by evaluating here.
//js
						Real sig;
						sigma( pt2, cen1, vdist, sig );
						if ( sig > 179.0 ) {
							sig = 179.0;
						} else if ( sig < 0.0 ) {
							sig = 0.0;
						}
						int isig = static_cast< int >( sig / 10 ) + 1;
						if ( dist > 6.4 ) {
							dist = 6.4;
						} else if ( dist < 3.5 ) {
							dist = 3.5;
						}
						int idist = static_cast< int >( ( dist - 3.5 ) / 0.25 ) + 1;

							//tempscore_rsigma = rsigma_dot_(idist,isig,sign1,sign2);

							// The definition of dimer signs (sign1,sign2)
							// appears inverted (1 should be 2, vice versa).
						Real tempscore_rsigma = rsigma_dot_(idist, isig, 3 - sign1, 3 - sign2);

							// Modify sigma potential to no longer give an rsigma bonus
							//  to strands in wrong register.
						if (sign1 == 1 & sign2 == 1 & sig < 110. & sig > 70.) tempscore_rsigma = 0.0;
						if (sign1 == 1 & sign2 == 2 & (sig < 75. | sig > 95.)) tempscore_rsigma = 0.0;;
						if (sign1 == 2 & sign2 == 1 & (sig < 90. | sig > 110.)) tempscore_rsigma = 0.0;
						if (sign1 == 2 & sign2 == 2 & sig < 120. & sig > 80.) tempscore_rsigma = 0.0;

						if ( tempscore_rsigma > 0.0 ) {
//car add in all unfavorable rsigma scores (these pairs may not contribute
//car to ss_score if they have a favorable dimer_pair_score)
//car note there are no positive scores in structure_db.cc as of 7/03
							if (symmetric) { // multiply with score factors for the edge
								rsigma_score += tempscore_rsigma*symm_info->score_multiply(SS_resnum_ss1, SS_resnum_ss2);
							} else {
								rsigma_score += tempscore_rsigma;
							//							std::cout << "rsigma: " << ss1 << ' ' << ss2 << ' ' << tempscore_rsigma << std::endl;
							}
						} else {
//car and favorable ones if phitheta and dot score favorable & dist<6.5
//car note that rsigma is not subject to the favorable-dimer-interactions-must
//car be-consistent-pairwise-strands rule that is applied to the ss_score
//car below
							if ( ss_orient == 1 ) {
								if (symmetric) { // multiply with score factors for the edge
									rsigma_score += tempscore_rsigma*symm_info->score_multiply(SS_resnum_ss1, SS_resnum_ss2);
								} else {
									rsigma_score += tempscore_rsigma;
								}
							}
								 // ss_orient(ss1,ss2)
						}

						// NOTE: there was a barcode section here in the original classic rosetta code that was taken out

						//////////////////////////////
						// barcode sspair constraints:
						//
						if (pose.data().has( core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO ) ) {
							// these are scored regardless of dimer_seqsep
							// (eg to penalize hairpins even in score5)
							// but only if the interaction is favorable
							// ie dimer_score < 0.0
							//
							// the score is equal to
							// -1 * weight * dimer_score
							//
							// the -1 is so that a positive weight in the barcode file
							// will penalize occurrence of the feature
							//
							float tmp_distscore( lowstrand > 0.5 ? ssdist_12 : 0.0 );
							// reconstruct dimer_score w/o the dimer_seqsep cutoff:
							float const dimer_score
								( phithetascore + tmp_distscore + dotscore );
							if ( dimer_score < 0.0 ) {
								// score is added directly to sspair score

								// offset the residue position for the second dimer
								// if they are anti-parallel. This has the effect that
								// if the dimers are aligned, then the residues passed
								// in will in fact be aligned as well.
								int const dimer2_offset( theta > 90.0 ? 1 : 0 );

								core::Size const resnum_ss1( SS_resnum_ss1 );
								core::Size const resnum_ss2_offset( SS_resnum_ss2+dimer2_offset);

								//static_cast< SS_KILLHAIRPINS_Info const &>( pose.data().get( core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO )));

								float const score_delta( hairpin_killing_score( pose, resnum_ss1,
																																resnum_ss2_offset,
																																theta, dimer_score ));
								if (symmetric) { // multiply with score factors for the edge
									ss_score += dimer_pair_score*symm_info->score_multiply(SS_resnum_ss1, SS_resnum_ss2); // [ l ]; // dimer_pair_score(ss1,ss2)
									ss_score += score_delta*symm_info->score_multiply(SS_resnum_ss1, SS_resnum_ss2);
								} else {
									ss_score += score_delta;
								}
								//if ( false && std::abs( score_delta ) > 0.1 )
								if ( false && std::abs( score_delta ) > 0.1 )
									std::cout << "Apply score: " << score_delta << ' ' <<
										SS_resnum_ss1 << ' ' << SS_resnum_ss2 << std::endl;
							}
						}

						//BARCODE END

						// rhiju Allow penalty for locally paired strands to be set by user.
						int sequence_separation = std::abs( strands.SS_resnum(ss2) - strands.SS_resnum(ss1) );
						if ( sequence_separation <= 11)	{
							dimer_pair_score += wts.get_localstrandpair_penalty() *
								SS_penalty(sequence_separation);
						}

//car now add  pairs that interact favorably to a list for later trimming...
//car pairs that interact unfavorably count toward the ss_score
//car and then are zeroed out.

//							std::cout << "rsigma: " << ss1 << ' ' << ss2 << ' ' << tempscore_rsigma << std::endl;

						if ( dimer_pair_score < 0.0 ) {
//car here we should also build a list of what dimers ss1 and ss2 are close
//car to so that we can speed up the check for disallowed pairs at the end
//car of this function
							dimer_pairs.push_back( new DimerPairing( ss1, ss2, sign1, sign2, dimer_pair_score ) );

							if (symmetric) {  // multiply with score factors for the edge
									// We need to have all possible dimer pairs present for the sheet detection to work. The energy graph does not contain all residues, only the ones in the scoring subunit and edges to the scoring subunit. So we need to add dimer pairongs for no-scoring subunits to make the list complete.

									// Map a residue pair to all other symmetric pairs and store in a vector of pairs
									std::vector < std::pair < Size, Size > > symm_ss_pairs (symm_info->map_symmetric_res_pairs( SS_resnum_ss1, SS_resnum_ss2 ) );
									std::vector< std::pair < Size, Size > >::const_iterator it_begin = symm_ss_pairs.begin();
									std::vector< std::pair < Size, Size > >::const_iterator it_end = symm_ss_pairs.end();
									std::vector< std::pair < Size, Size > >::const_iterator it;
									for ( it = it_begin; it != it_end ; ++it ) {
										Size clone_SS_resnum_ss1 ( it->first );
										Size clone_SS_resnum_ss2 ( it->second );
										int symm_ss1_clone ( strands.SS_dimer ( clone_SS_resnum_ss1 ) );
										int symm_ss2_clone ( strands.SS_dimer ( clone_SS_resnum_ss2 ) );
										if ( symm_info->bb_follows( clone_SS_resnum_ss1 ) == 0 || symm_info->bb_follows( clone_SS_resnum_ss2 ) == 0 ) continue; // Do not add we are in a scoring subunit. These are already added...
										if ( symm_ss1_clone == 0 || symm_ss2_clone == 0 ) continue; // if there are no strands at clone postions...
										if ( clone_SS_resnum_ss1 > clone_SS_resnum_ss2 ) continue; // strands have to be in sequence order

										dimer_pairs.push_back( new DimerPairing( symm_ss1_clone, symm_ss2_clone, sign1, sign2, dimer_pair_score ) );
								}
							}

						} else {
							if (symmetric) { // multiply with the score_mulptipy for the edge
								if ( symm_info->bb_follows(SS_resnum_ss1) !=0 && symm_info->bb_follows(SS_resnum_ss2) != 0 ) continue; // we have already added this score
								ss_score += dimer_pair_score*symm_info->score_multiply(SS_resnum_ss1, SS_resnum_ss2); // [ l ]; // dimer_pair_score(ss1,ss2)
						} else {
							ss_score += dimer_pair_score; // [ l ]; // dimer_pair_score(ss1,ss2)
						}
						}
					}      // good phi/th ( phithetascore<0.0 )

//js collect data on neighbors of dimers for making sheets later
//js somewhat odd method of selecting the two neighbors for determining what
//js makes a sheet.  Selects the first and last dimer that follow the current
//js dimer in sequence.  Not completely sure it is even symmetric.
//js There has been discussion of replacing this with something like
//js a) a check that these dimers are also well oriented according to SS_strand_score
//js b) or possibly using the good two that are selected with later.


					if ( dist <= 5.5 ) { // && ss_orient(i,j) == 1} ) {
						if ( dimer_neighbor_1ss1 == 0 ) {
							dimer_neighbor_1ss1 = ss2;
						} else if ( strands.SS_strand(dimer_neighbor_1ss1) != strands.SS_strand(ss2) ) {
							dimer_neighbor_2ss1 = ss2;
						}
						if ( ( dimer_neighbor_1ss1 != ss2 ) &&
								 ( dimer_neighbor_2ss1 != ss2 ) ) {
							if ( strands.dimer_neighbor(1,ss2) == 0 ) {
								strands.dimer_neighbor(1,ss2) = ss1;
							} else if ( strands.SS_strand(strands.dimer_neighbor(1,ss2)) != SS_strand_ss1 ) {
								strands.dimer_neighbor(2,ss2) = ss1;
							}
						}
					}      // within sheet distance cutoff (5.5 Angstrom)

					if (symmetric) {  // multiply with score factors for the edge
									// We need to have all possible dimer pairs present for the sheet detection to work. The energy graph does not contain all residues, only the ones in the scoring subunit and edges to the scoring subunit. So we need to add dimer neighbors for no-scoring subunits to make the list complete.


						std::vector < std::pair < Size, Size > > symm_ss_pairs (symm_info->map_symmetric_res_pairs( SS_resnum_ss1, SS_resnum_ss2 ) );
						std::vector< std::pair < Size, Size > >::const_iterator it_begin = symm_ss_pairs.begin();
						std::vector< std::pair < Size, Size > >::const_iterator it_end = symm_ss_pairs.end();
						std::vector< std::pair < Size, Size > >::const_iterator it;
						for ( it = it_begin; it != it_end ; ++it ) {
							Size clone_SS_resnum_ss1 ( it->first );
							Size clone_SS_resnum_ss2 ( it->second );
							int symm_ss1_clone ( strands.SS_dimer ( clone_SS_resnum_ss1 ) );
							int symm_ss2_clone ( strands.SS_dimer ( clone_SS_resnum_ss2 ) );

							if ( symm_info->bb_follows( clone_SS_resnum_ss1 ) == 0 || symm_info->bb_follows( clone_SS_resnum_ss2 ) == 0 ) continue;
							if ( symm_ss1_clone == 0 || symm_ss2_clone == 0 ) continue;
							if ( clone_SS_resnum_ss1 > clone_SS_resnum_ss2 ) continue;

							int & dimer_neighbor_1ss1_clone( strands.dimer_neighbor(1,symm_ss1_clone) );
							int & dimer_neighbor_2ss1_clone( strands.dimer_neighbor(2,symm_ss1_clone) );

							if ( dist <= 5.5 ) { // && ss_orient(i,j) == 1} ) {
								if ( dimer_neighbor_1ss1_clone == 0 ) {
									dimer_neighbor_1ss1_clone = symm_ss2_clone;
								} else if ( strands.SS_strand(dimer_neighbor_1ss1_clone) != strands.SS_strand(symm_ss2_clone) ) {
									dimer_neighbor_2ss1_clone = symm_ss2_clone;
								}
								if ( ( dimer_neighbor_1ss1_clone != symm_ss2_clone ) &&
										 ( dimer_neighbor_2ss1_clone != symm_ss2_clone ) ) {
									if ( strands.dimer_neighbor(1,symm_ss2_clone) == 0 ) {
										strands.dimer_neighbor(1,symm_ss2_clone) = symm_ss1_clone;
									} else if ( strands.SS_strand(strands.dimer_neighbor(1,symm_ss2_clone)) != strands.SS_strand( symm_ss1_clone ) ) {
										strands.dimer_neighbor(2,symm_ss2_clone) = symm_ss1_clone;
									}
								}
							}      // within sheet distance cutoff (5.5 Angstrom)
						}
					}

				}         // within distance cutoff (6.5 Angstrom)
			}            // not the same strand
		} // nbrs of SS_resnum_ss1
	} // ss1


//car okay, we've now scored all dimer pairs and we have to figure out
//car if there are inconsistent pairs on our list (ie each dimer can
//car interact with at most two strands, one on each side; note that
//car a dimer can interact with multiple dimers on a single side, as long
//car as they're in the same strand)
//car
//car What follows below is very inefficient, but here's the idea:
//car  loop through the list of favorable pairs to find the best score;
//car  this pair is 'allowed' and added to the total.  This pair now
//car  defines which strand is to one side of dimer1 and which strand is
//car  to one side of dimer2 (the side is defined by signdimer).  So
//car  now we loop through all possible dimer pairs and if they are
//car  not in
//car  then mark them as disallowed by setting the score for that pair
//car   (ie dimer_pair_score) to zero.
//car  Now repeat this for the next-best scoring pair that remains in the
//car  the favorable pairs list.  When there are no more allowed pairs,
//car  escape.
//car this will be sped up enormously by
//car   1. sort the allowed pairs list so we just have to loop through it 1x
//car   2. save a list (above) of every dimer that a particular dimer is
//car      near so that marking disallowed pairs is fast

//  TOTAL SS_SCORE (DONT LET A DIMER INTERACT WITH MORE THAN 2 OTHERS)
//    UNLESS THE SCORE IS UNFAVORABLE, THEN IT STILL GETS THE PENALTY
//car    (note the penalties have already been added in)



	std::sort( dimer_pairs.begin(), dimer_pairs.end(), dimer_pairing_pointer_sorter );

	for ( DimerPairings::iterator it= dimer_pairs.begin(), ite= dimer_pairs.end(); it != ite; ++it ) {
		DimerPairing const & pairing( **it );
		if ( !pairing.valid() ) continue;
		if (symmetric) { // multiply with score factors for the edge
			int dim1 ( strands.SS_resnum( pairing.dimer1()) );
			int dim2 ( strands.SS_resnum( pairing.dimer2()) );
			if ( symm_info->bb_follows(dim1) !=0 && symm_info->bb_follows(dim2) != 0 ) continue;
			ss_score += pairing.score()*symm_info->score_multiply( dim1, dim2 );
		} else {
			ss_score += pairing.score();
		}

		int const dimer1( pairing.dimer1() );
		int const dimer2( pairing.dimer2() );

		int const sign1( pairing.sign1() );
		int const sign2( pairing.sign2() );

		assert( dimer1 < dimer2 );

//		std::cout << "ss_dimer: " << strands.SS_resnum(dimer1) << ' ' << strands.SS_resnum(dimer2) << ' ' <<
//			pairing.score() << std::endl;


//DB  BONUS FOR NONLOCAL PAIRS!!
//car this was added to correct for an excess of local strand pairs
//car and really should probably go above where each dimer is calculated...
//car  ask DB about this...

		if ( std::abs( strands.SS_resnum(dimer2) - strands.SS_resnum(dimer1) ) > 11 &&
				 std::abs( strands.SS_strand_end(1,dimer2) - strands.SS_strand_end(2,dimer1) ) > 8 ) {
			if (symmetric) { // multiply with score factors for the edge
				int dim1 ( strands.SS_resnum( pairing.dimer1()) );
				int dim2 ( strands.SS_resnum( pairing.dimer2()) );
				if ( symm_info->bb_follows(dim1) !=0 && symm_info->bb_follows(dim2) != 0 ) continue; // we have already added these
				ss_score += -0.2*symm_info->score_multiply(strands.SS_resnum(dimer1), strands.SS_resnum(dimer2) );
			} else {
				ss_score -= 0.2;
			}
		}

		int const SS_strand_dimer1 = strands.SS_strand( dimer1 );
		int const SS_strand_dimer2 = strands.SS_strand( dimer2 );

// ARE THERE OTHER DIMERS INTERACTING WITH THE BEST PAIR?

		DimerPairings::iterator it2( it );
		++it2;

		for ( ; it2 != ite; ++it2 ) {

			DimerPairing & other( **it2 );
			if ( !other.valid() ) continue;

			int const other_strand1( strands.SS_strand( other.dimer1() ) );
			int const other_strand2( strands.SS_strand( other.dimer2() ) );

//car if dimer1 and ss2 interact favorably
//car and ss2 and dimer2 are in different strands
//car and ss2 is on the same side of dimer1 as dimer2 is...
//car then mark this pair as dissallowed

			if ( ( other.dimer1() == dimer1 && other_strand2 != SS_strand_dimer2 && other.sign1() == sign1 ) ||
					 ( other.dimer2() == dimer1 && other_strand1 != SS_strand_dimer2 && other.sign2() == sign1 ) ||
					 ( other.dimer1() == dimer2 && other_strand2 != SS_strand_dimer1 && other.sign1() == sign2 ) ||
					 ( other.dimer2() == dimer2 && other_strand1 != SS_strand_dimer1 && other.sign2() == sign2 ) ) {
				other.valid( false );
			}
		} // it2
	} // it

	// modify by proper weighting
	ss_score *= 0.498 * 0.75;
	rsigma_score *= 0.1;
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
float
SecondaryStructurePotential::hairpin_killing_score(
	pose::Pose const & pose,
	Size const & pos1,
	Size const & pos2,
	Real const & theta,
	float const & ss_score
) const
{
	 ////using core::pose::datacache::CacheableDataType::core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO;
	basic::ProfileThis doit( basic::SECONDARY_STRUCTURE_ENERGY );
	float total_score(0.0);
	// will this be too slow? called inside sspair score
	// for each dimer pair... well -- actually only for dimer_pairs
	// with favorable dimer scores... probably not too many then

	//ITERATE OVER HAIRPINS
	//for ( constraint_const_iterator it = cst_list.begin(),
	//				it_end = cst_list.end(); it != it_end; ++it ) {
	if (pose.data().has( core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO ) ) {

		runtime_assert( dynamic_cast< SS_Killhairpins_Info const *>( &( pose.data().get( core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO ))));
		SS_Killhairpins_Info kill_hairpin_info( static_cast< SS_Killhairpins_Info const &>( pose.data().get( core::pose::datacache::CacheableDataType::SS_KILLHAIRPINS_INFO )));

		if ( theta < 0.0 || theta > 180.0 ) {
			//CHANGE TO RUNTIME ASSERT
			std::cout << "WARNING:: sspair_constraint_score: theta should lie " <<
				" between 0 and 180 degrees! theta = "<< theta << std::endl;
			return 0.0;
		}

		//ALL HAIRPINS ARE ANTIPARALLEL
		//int const orientation( int_value );

		if ( (kill_hairpin_info.check_hairpin(pos1, pos2)) && kill_hairpin_info.kill_antiparallel() &&
				 ( theta > 90.0 ) ) { //theta > 90.0 means the pairing is antiparallel
			float const penalty(100.0); //PENALTY VALUE, MAYBE PUT IN DATA CACHE???
			// penalty constraint or unmatched bonus constraint
			total_score += -1 * penalty * ss_score;
		}

		if ( (kill_hairpin_info.check_hairpin(pos1, pos2)) && kill_hairpin_info.kill_parallel() &&
				 ( theta <= 90.0 ) ) { //theta <= 90.0 means the pairing is parallel
			float const penalty(100.0); //PENALTY VALUE, MAYBE PUT IN DATA CACHE???
			// penalty constraint or unmatched bonus constraint
			total_score += -1 * penalty * ss_score;
		}

	}
	return total_score;
}

//////////////////////////////////////
//////////////////////////////////////


/// @details
/// apl This function reads the strand pairing information for all dimer
/// apl neighbors to determine how many strands are in each sheet.  It
/// apl uses a fast conected-component detection data structure to do so.
/// apl Once the number of strands for each sheet is known, a per-sheet "poker
/// apl hand score" is computed and these scores summed.
/// apl Old comments follow.
///
///js This function takes a set of dimer neighbors, and determines how
///js many sheets there, and how many strands are in each sheet
///js This information is then used to calculate the "poker hand" score,
///js which reflects to probability of that distribution of strands and
///js sheets.
///js In current version, it seems to simply penalize sheets with fewer
///js strands compared to those with more strands.
///
///js This function looks at a list of dimers, which contains up to
///js two neighbors for each dimer.  In priniciple these neighbors would
///js be hydrogen bond partners in neighboring strands.  This function
///js will take this list, however it is made.
///
///js Currently, dimer neighbors are defined somewhat arbitrarily.
///js If I understand the code correctly, the first and last dimers in
///js sequence that follow the current dimer, and that are within 6.5
///js angstroms, are the neighbors.  There is no orientation dependence
///js on what counts as a strand neighbor.
///
///js A sheet is then loosely defined by all the strands that are connected
///js by one of these neighbors.  This amounts to "single-linkage clustering."
///js A sheet is determined by the set of strands that can be linked by
///js the dimer neighbors.  Note that a one neighbor is enough to say that
///js two strands are in a sheet.
///
///js The final score is put into sheet_score, and is determined by the
///js number of sheets of each size.
///
///js   Basic strategy:  go through all dimers, finding neigboring strands
///js   Because each dimer can only have two neighbors we can first search
///js   all the way down one side, and then the other.  This will give a set
///js   of dimers that are connected.  The strands in which these dimers reside
///js   are then considered connected into sheets.  By going through all
///js   dimers with the proper bookkeeping, we can determine how many sheets
///js   there are, and how many strands in each.

Real
SecondaryStructurePotential::sheets_from_dimers(
	pose::Pose const & pose
) const
{
	basic::ProfileThis doit( basic::SECONDARY_STRUCTURE_SHEETS_FROM_DIMERS_ENERGY );
	SS_Info const & ss_info( retrieve_const_ss_info_from_pose( pose ) );
	Strands const & strands( ss_info.strands() );
	//std::cout << "strands.total_strands " <<  strands.total_strands << std::endl;

// We need to set up some stuff for symmetry
	bool symmetric=false;
	core::conformation::symmetry::SymmetryInfoOP symm_info = NULL;
	if (core::pose::symmetry::is_symmetric( pose ) ) {
		symmetric=true;
		SymmetricConformation const & SymmConf (
		    dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
		symm_info = SymmConf.Symmetry_Info()->clone();
	}


	graph::DisjointSets sheet_sets( strands.total_strands );
	for ( int ii = 1; ii <= strands.total_SS_dimer; ++ii ) {
		//std::cout << "strands.SS_strand(" << ii << ") = " << strands.SS_strand( ii ) << " with neighbs: " << strands.dimer_neighbor( 1, ii ) << " " << strands.dimer_neighbor( 2, ii ) << std::endl;
		int const ii_sheet = strands.SS_strand( ii );
		for ( int direction = 1; direction <= 2; ++direction ) {
			int const ii_neighbor = strands.dimer_neighbor( direction, ii );
			if ( ii_neighbor != 0 ) {
				int const ii_neighbor_sheet = strands.SS_strand( ii_neighbor );
				sheet_sets.ds_union( ii_sheet, ii_neighbor_sheet );
			}
		}
	}


	/// This code duplicates r++ behavior -- but it's not at all clear this is what should be done
	/*utility::vector1< bool > visited( strands.total_SS_dimer, false );
	graph::DisjointSets sheet_sets( strands.total_strands );

	for ( int ii = 1; ii <= strands.total_SS_dimer; ++ii ) {
		if ( visited[ ii ] ) continue;
		visited[ ii ] = true;
		for ( int direction = 1; direction <= 2; ++direction ) {
			int searching = ii;
			while ( true ) {
				int neighbor = strands.dimer_neighbor( direction, ii );
				if ( neighbor != 0 ) {
					if ( visited[ neighbor ] ) break;
					int const searching_sheet = strands.SS_strand( searching );
					int const neighbor_sheet = strands.SS_strand( neighbor );
					sheet_sets.ds_union( searching_sheet, neighbor_sheet );
					visited[ neighbor ] = true;
					searching = neighbor;
				} else {
					break;
				}
			}
		}
	}*/

	if ( symmetric ) {
		// The basic idea is to go through all sheets and weigh the sheet_score by the multiply_score factors.
		// If a sheet is across a subunit interface the score gets weighted by the score factors for that interface.
		// Otherwise the weight is equal to the weight of the scoring subunit

		// We need to be able to know in which subunit a particular strand resides. So we make a map of residue
		// numbers of the end of a strand to the strand number.
		utility::vector1< Size > sheet_sizes_sym;
		utility::vector1< Size > strand_to_SS_resnum (strands.total_strands, 1 );
		for ( int i = 2; i <= strands.total_SS_dimer; ++i ) {
			if ( strands.SS_strand_end(1,i) != strands.SS_strand_end(1,i-1) ) {
				strand_to_SS_resnum[ strands.SS_strand( i ) ] = strands.SS_strand_end(1,i);
			}
		}

		// Create a reduced set of strands that have scoring weights > 0.
		utility::vector1< Size > weight_sheets;
		graph::DisjointSets reduced_sheet_sets;
		std::map< Size, utility::vector1< Size > > const set_and_nodes ( sheet_sets.sets() );
		std::map< Size, utility::vector1< Size > >::const_iterator it_start = set_and_nodes.begin();
		std::map< Size, utility::vector1< Size > >::const_iterator it_end = set_and_nodes.end();
		// Loop over all sheets
		for ( std::map< Size, utility::vector1< Size > >::const_iterator it = it_start; it != it_end; ++it ) {
			utility::vector1< Size > const node_list ( sheet_sets.nodes_in_set( it->first ) );
			utility::vector1< Size >::const_iterator itn_start = node_list.begin();
			 utility::vector1< Size >::const_iterator itn_end = node_list.end();
			Size strand_res_native = 0;
			Size weight = 0;
			// Loop over all strands in sheet
			for ( utility::vector1< Size >::const_iterator itn = itn_start; itn != itn_end; ++itn ) {
				Size strand_res (strand_to_SS_resnum[ *itn ] );
				// determine the weight. Two cases: in scoring subunit or across an inteface to the scoring
				// subunit
				if ( symm_info->bb_is_independent( strand_res ) )  {
				 strand_res_native = strand_res;
				 weight =  symm_info->score_multiply(strand_res_native,strand_res_native);
				} else {
					if ( strand_res !=0 ) continue; //if we don't have a strand in the scoring subunit
					Size new_weight (symm_info->score_multiply(strand_res_native, strand_res) );
					if (  new_weight > weight ) weight = new_weight;
				}
			}
			// Only add sheets with weights > 0
			if ( weight > 0 ) {
				sheet_sizes_sym.push_back( it->second.size() );
				weight_sheets.push_back( weight );
			}
		}

		Real sheet_score( 0.0 );
		Size const size_four( 4 );
		for ( Size ii = 1; ii <= sheet_sizes_sym.size(); ++ii ) {
			sheet_score += weight_sheets[ ii ]*m_term_( std::min( sheet_sizes_sym[ ii ], size_four ) );
		}
		return sheet_score * 2.019; // pre-weighting from r++::structure.cc:733
	}

	utility::vector1< Size > sheet_sizes = sheet_sets.disjoint_set_sizes();

	Real sheet_score( 0.0 );
	Size const size_four( 4 );
	for ( Size ii = 1; ii <= sheet_sizes.size(); ++ii ) {
		sheet_score += m_term_( std::min( sheet_sizes[ ii ], size_four ) );
	}

	/// APL -- to be ported in the future
	//if ( get_handedness_score_flag() && files_paths::use_filter(files_paths::sheet_type) ){
	//	int result = 0;
	//	sheet_filter::SheetFilter sf(position_, secstruct_, total_residue_);
	//	sf.compute_result(result); // Trigger evaluation if Ingo's sheet filter, including handedness checks.
	//	sheet_score += sf.get_handedness_score();
	//}

	//std::cout << " SHEET-SCORE: " << sheet_score << std::endl;
	return sheet_score * 2.019; // pre-weighting from r++::structure.cc:733
}





/// @details identifies secondary structure at the start of scoring and loads the Strands/Helices
/// data
void
SecondaryStructurePotential::identify_ss(
	pose::Pose const & pose,
	Helices & helices,
	Strands & strands
) const
{

	conformation::Conformation const & conf( pose.conformation() ); // get secstruct info from this guy
	int const total_residue( pose.total_residue() );

	FArray1D_int dimer_type( total_residue ); // what type of dimer each position is
	// 0-none, 1-SS 2-HH

	//car find SSdimers and HHdimers
	//car find N-terminal ends of SS elements
	//car save a map of the dimers to aid finding C-termini later
	int lastH = 1;
	int lastE = 1;
	int lastL = 1;
	helices.total_HH_dimer = 0;
	strands.total_SS_dimer = 0;
	for ( int i = 1; i <= total_residue; ++i ) strands.SS_dimer(i) = 0; // initialize

	//std::cout << "secstruct: ";
	for ( int i = 1; i <= total_residue; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( "REPLONLY" ) ){
				continue;
		}
		char const c = conf.secstruct(i);
		//std::cout << conf.secstruct(i);
		if ( c == 'L' ) {
			lastL = i;
		} else if ( c == 'E' ) {
			lastE = i;
		} else if ( c == 'H' ) {
			lastH = i;
		}

		dimer_type(i) = 0;
		if ( i <= total_residue-1 ) {
			if ( conf.secstruct(i) == 'E' && conf.secstruct(i+1) == 'E' ) {
				++strands.total_SS_dimer;
				strands.SS_resnum(strands.total_SS_dimer) = i;
				strands.SS_dimer ( i ) = strands.total_SS_dimer; // reverse mapping -- add to ctor,operator=,etc
				strands.SS_strand_end(1, strands.total_SS_dimer) = std::max(lastH,lastL);
				dimer_type(i) = 1;
			}
		}
		if ( i >= 2 && i <= total_residue-2 ) {
			if ( conf.secstruct(i-1) == 'H' && conf.secstruct(i) == 'H' && conf.secstruct(i+1) == 'H' &&
			 conf.secstruct(i+2) == 'H' ) {
				++helices.total_HH_dimer;
				helices.HH_resnum(helices.total_HH_dimer) = i;
				//helices.HH_dimer ( i ) = helices.total_HH_dimer; // reverse mapping, is this needed?
				helices.HH_helix_end(1, helices.total_HH_dimer) = std::max(lastE,lastL);
				dimer_type(i) = 2;
			}
		}
	}
	//std::cout << std::endl;

	lastH = total_residue;
	lastE = total_residue;
	lastL = total_residue;
	int iHH = helices.total_HH_dimer;
	int iSS = strands.total_SS_dimer;

	for ( int i = total_residue; i >= 1; --i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( "REPLONLY" ) ){
				continue;
		}
		char const c = conf.secstruct(i);
		if ( c == 'L' ) {
			lastL = i;
		} else if ( c == 'E' ) {
			lastE = i;
		} else if ( c == 'H' ) {
			lastH = i;
		}

		if ( dimer_type(i) == 1 ) {
			strands.SS_strand_end(2, iSS) = std::min(lastH, lastL);
			--iSS;
		} else if ( dimer_type(i) == 2 ) {
			helices.HH_helix_end(2, iHH) = std::min(lastE, lastL);
			--iHH;
		}
	}

	//car assign strand numbers to ss_dimers
	if ( strands.total_SS_dimer > 0 ) {
		strands.SS_strand(1) = 1; // first ss-dimer belongs to strand 1
		strands.total_strands = 1;

		for ( int i = 2; i <= strands.total_SS_dimer; ++i ) {
			if ( strands.SS_strand_end(1,i) == strands.SS_strand_end(1,i-1) ) {
				strands.SS_strand(i) = strands.SS_strand(i-1); // same strand
			} else {
				++strands.total_strands;
				strands.SS_strand(i) = strands.total_strands;
			}
		}
	} else {
		strands.total_strands = 0;
	}
	///std::cout << " strands.total_strands " <<  strands.total_strands;// << std::endl;

}


/// @brief function reads in two points in sequence and returns two points in space,
/// @brief the endpoints of the axis through an alpha-helix
void
SecondaryStructurePotential::helix_end(
	int const & pos1,
	BB_Pos const & bb_pos,
	Vector & p1,
	Vector & p2
) const
{
	int const s1 = pos1-1;
	int const s2 = pos1;
	int const s3 = pos1+1;
	int const s4 = pos1+2;

	static Real const eleven_inv = 1.0 / 11.0;

	Vector const Epos_sum( (                  bb_pos.CA( s1 ) + bb_pos.C( s1 ) ) +
												 ( bb_pos.N( s2 ) + bb_pos.CA( s2 ) + bb_pos.C( s2 ) ) +
												 ( bb_pos.N( s3 ) + bb_pos.CA( s3 ) + bb_pos.C( s3 ) ) +
												 ( bb_pos.N( s4 ) + bb_pos.CA( s4 )                  ) );

	p1 = ( Epos_sum + bb_pos.N( s1 ) ) * eleven_inv;
	p2 = ( Epos_sum + bb_pos.C( s4 ) ) * eleven_inv;

// 	for ( int i = 1; i <= 3; ++i ) {
// 		Real const Epos_sum =
// 		 Eposition(i,2,s1) + Eposition(i,4,s1) +
// 		 Eposition(i,1,s2) + Eposition(i,2,s2) + Eposition(i,4,s2) +
// 		 Eposition(i,1,s3) + Eposition(i,2,s3) + Eposition(i,4,s3) +
// 		 Eposition(i,1,s4) + Eposition(i,2,s4);
// 		p1(i) = ( Epos_sum + Eposition(i,1,s1) ) * eleven_inv;
// 		p2(i) = ( Epos_sum + Eposition(i,4,s4) ) * eleven_inv;
// 	}
}


/// @brief calculate sum of dot product of the co vectors of strand dimers ss1 and ss2
/// @brief with the vector connecting the midpoints of the dimer vectors (vdist)
/// @brief also determine return the sign of the dot products for each dimer
/// @brief to determine which direction the CO groups point
void
SecondaryStructurePotential::pair_dp(
	int const & ss1,
	int const & ss2,
	BB_Pos const & bb_pos,
	Real & dp,
	Vector const & vdist,
	int & sign1,
	int & sign2
) const
{

//car parameters
//	static Real const dist_co = { 1.231015f }; // length of C=O bond
	static Real const dist_co_inv = { 1.0f / 1.231015f };

//car local
	Vector temp;

	dp = 0.0;
	Real dp1 = 0.0;
	Real dp2 = 0.0;
	Real sdp1 = 0.0;
	Real sdp2 = 0.0;

	for ( int i = ss1; i <= ss1+1; ++i ) {
		if ( i == ss1+1 ) {
			temp = dist_co_inv * ( bb_pos.C(i) - bb_pos.O(i) );
			// 			for ( int j = 1, l = l0; j <= 3; ++j, ++l ) {
			// 				temp(j) = -( Eposition[ l+3 ] - Eposition[ l ] ) * dist_co_inv; // 5=O, 4=C
			// // 				temp(j) = -( Eposition[ l+3 ] - Eposition[ l ] ) * dist_co_inv; // 5=O, 4=C
			// 				 //       -( Eposition(j,5,i) - Eposition(j,4,i) )
			// 			}
		} else {
			temp = dist_co_inv * ( bb_pos.O(i) - bb_pos.C(i) );
			// 			for ( int j = 1, l = l0; j <= 3; ++j, ++l ) {
			// 				temp(j) = ( Eposition[ l+3 ] - Eposition[ l ] ) * dist_co_inv; // 5=O, 4=C
			// 				 //       ( Eposition(j,5,i) - Eposition(j,4,i) )
			// 			}
		}
		//if ( vdist(1) != 10.0 ) {   // why is this checked?  (car)
		Real const tempdot = temp.dot( vdist );
		dp1 += std::abs(tempdot);
		sdp1 += tempdot;
			//}
	}
	dp1 *= 0.5;

	for ( int i = ss2; i <= ss2+1; ++i ) {
		if ( i == ss2+1 ) {
			temp = dist_co_inv * ( bb_pos.C(i) - bb_pos.O(i) );
			// 			for ( int j = 1, l = l0; j <= 3; ++j, ++l ) {
			// 				temp(j) = -( Eposition[ l+3 ] - Eposition[ l ] ) * dist_co_inv; // 5=O, 4=C
			// 				 //       -( Eposition(j,5,i) - Eposition(j,4,i) )
			// 			}
		} else {
			temp = dist_co_inv * ( bb_pos.O(i) - bb_pos.C(i) );
			// 			for ( int j = 1, l = l0; j <= 3; ++j, ++l ) {
			// 				temp(j) = ( Eposition[ l+3 ] - Eposition[ l ] ) * dist_co_inv; // 5=O, 4=C
			// 				 //       ( Eposition(j,5,i) - Eposition(j,4,i) )
			// 			}
		}
		//if ( vdist(1) != 10.0 ) {
		Real const tempdot = temp.dot( vdist );
		dp2 += std::abs(tempdot);
		sdp2 += tempdot;
		//}
	}
	dp2 *= 0.5;

	// 	if ( vdist(1) == 10.0 ) {
	// 		dp = 10.0;
	// 	} else {
	dp = dp1 + dp2;
//js These signs tell whether the first c=o bond vector of a dimer points
//js at the other dimer.  sign1 = 1 means that the first c=o bond of dimer1
//js points at dimer2.  sign2 = 1 means that dimer2 points at dimer1.  When
//js sign1 or sign2 equals 2, that dimer points away from the other dimer.
	sign1 = ( sdp1 > 0.0 ? 2 : 1 );
	sign2 = ( sdp2 < 0.0 ? 2 : 1 );
}


/// @brief identifies the sequence separation along the fold tree
/// lin the old way to calculate the sequence separation takes an asumption of no-break chain
/// lin when there is chain break between pos1 and pos2, we add a gap to make a correct calculation in ss energy
int
SecondaryStructurePotential::get_foldtree_seqsep(
	pose::Pose const & pose,
	int pos1,
	int pos2,
	int gap_size
) const
{
	if( pose.fold_tree().is_simple_tree() ) return std::abs( pos1 - pos2 );

	int begin ( std::min(pos1,pos2) );
	int end   ( std::max(pos1,pos2) );

	bool is_break ( false );

	for ( int i = begin; i < end; ++i ) {
		//if( pose.fold_tree().is_cutpoint(i) ) { is_break=true; break; }
		if( pose.residue_type(i).is_terminus() ) { is_break=true; break; }
	}

	if( is_break ) {
		return end-begin+gap_size;
	} else {
		return end-begin;
	}

}


/// @brief find the vector connecting the midpoints of two line segments
/// @brief defined by a1,a2 and a3,a4
/// @param[out] dist length of v21
/// @param[out] v21  vector connecting midpoints
void
SecondaryStructurePotential::dist_pair(
	Vector const & a1,
	Vector const & a2,
	Vector const & a3,
	Vector const & a4,
	Real & dist, // length of v21
	Vector & cen1,
	Vector & cen2,
	Vector & v21 // vector connecting midpoints
)
{
	//car find midpoint coordinates
	cen1 = Real( 0.5 )*( a1 + a2 );
	cen2 = Real( 0.5 )*( a3 + a4 );

	//car find distance between midpoint coordinates
	v21 = cen2 - cen1;
	dist = v21.length();
	// 	cen1(1) = ( a1(1) + a2(1) )*0.5f;
	// 	cen1(2) = ( a1(2) + a2(2) )*0.5f;
	// 	cen1(3) = ( a1(3) + a2(3) )*0.5f;

	// 	cen2(1) = ( a3(1) + a4(1) )*0.5f;
	// 	cen2(2) = ( a3(2) + a4(2) )*0.5f;
	// 	cen2(3) = ( a3(3) + a4(3) )*0.5f;

// 	v21(1) = cen2(1) - cen1(1);
// 	v21(2) = cen2(2) - cen1(2);
// 	v21(3) = cen2(3) - cen1(3);
// 	dist = std::sqrt(
// 	 v21(1) * v21(1) +
// 	 v21(2) * v21(2) +
// 	 v21(3) * v21(3) );
}


/// @details define a coordinate system with Z axis along cen1->a2 (v1),
/// xz plane defined by cen1->a2 and v21. and origin at cen1.
/// Find the spherical coordinates phi,and theta of point a4 if
/// vector cen2->a4 was moved such that cen2=cen1
///
/// @param[out] v21  vector connecting midpoints
void
SecondaryStructurePotential::spherical(
	Vector const & a2,
	Vector const & a4,
	Real & phi,
	Real & theta,
	Vector const & cen1,
	Vector const & cen2,
	Vector const & v21
)
{
	using numeric::conversions::degrees;
	using numeric::sin_cos_range;

	//car v1 vector from center to end of dimer1 vector
	//car v2                                   2

	Vector v1( a2 - cen1 );
	// 	v1(1) = a2(1) - cen1(1);
	// 	v1(2) = a2(2) - cen1(2);
	// 	v1(3) = a2(3) - cen1(3);
	Vector v2( a4 - cen2 );
	// 	v2(1) = a4(1) - cen2(1);
	// 	v2(2) = a4(2) - cen2(2);
	// 	v2(3) = a4(3) - cen2(3);

//car find unit vector along v1 = uz
	Vector const uz( v1.normalized_or_zero() );

//car find unit vector perpendicular v21 v1 plane
	Vector const uy( uz.cross( v21 ).normalized_or_zero() );

//car third unit vector to define coordinate system
	Vector const ux( uy.cross( uz ).normalized_or_zero() );

//car find projection of v2 onto each of these azes
	Real const v2x = v2.dot( ux ); //v2(1)*ux(1) + v2(2)*ux(2) + v2(3)*ux(3); // v2x=v2.ux
	Real const v2y = v2.dot( uy ); //v2(1)*uy(1) + v2(2)*uy(2) + v2(3)*uy(3); // v2y=v2.uy
	Real const v2z = v2.dot( uz ); //v2(1)*uz(1) + v2(2)*uz(2) + v2(3)*uz(3); // v2z=v2.uz

//car and length of v2
	Real r1 = v2.length();

//car unit vector along v21
	Vector const u21( v21.normalized_or_zero() );

//car projection of u21 on uz
// pb -- this doesnt seem to be used
	//Real const u21z = u21.dot( uz ); //u21(1)*uz(1) + u21(2)*uz(2) + u21(3)*uz(3); // u21z=u21.uz


//car the same thing in function calls (is this expensive?)
//$$$      unitvec(v1,uz);
//$$$      cros(uz,v21,uy);
//$$$      unitvec(uy,uy);
//$$$      cros(uy,uz,ux);
//$$$      unitvec(ux,ux);
//$$$
//$$$      v2x = dotprod(v2,ux);
//$$$      v2y = dotprod(v2,uy);
//$$$      v2z = dotprod(v2,uz);
//$$$      r1 = std::sqrt( ( v2(1) * v2(1) )+( v2(2) * v2(2) )+( v2(3) * v2(3) ) );
//$$$
//$$$      unitvec(v21,u21);
//$$$      u21z = dotprod(u21,uz);


	phi = 200.0; // why set to 200?  if v2y = 0, v2 lies in xz plane and phi
	// is zero; if v2x=0, v2 lies in yz plane and phi is 90 or -90
	if ( v2y != 0.0 && v2x != 0.0 ) {
		phi = degrees( std::atan2( v2y, v2x ) );
	}

	Real const v2z_r1 = v2z / r1;
	if ( r1 != 0.0 && std::abs( v2z_r1 ) <= 1.0 ) {
		theta = degrees( numeric::arccos( v2z_r1 ) ); //std::acos( sin_cos_range( v2z_r1 ) ) );
	} else if ( v2z_r1 > 1.0 ) {
		theta = 0.0;
	} else {
		theta = 180.0;
	}
}


/// @brief find angle sigma between vectors cen1->a2 and v21
/// @param[out] sig sigma
void
SecondaryStructurePotential::sigma(
	Vector const & a2,
	Vector const & cen1,
	Vector const & v21,
	Real & sig
)
{
	//car this could be done in spherical

	// 	subvec(a2,cen1,v1);
	// 	unitvec(v1,uz);
	// 	unitvec(v21,u21);
	Real const u21z = ( a2 - cen1 ).normalized().dot( v21.normalized() ); //dotprod(u21,uz);

	sig = numeric::conversions::degrees( numeric::arccos( u21z ) ); //std::acos( sin_cos_range( u21z ) );

	// 	sig = 200.0; // why 200? should be 0 or 180
	// 	if ( std::abs(u21z) <= 1.0 ) { //Objexx:SGM This logic is hackery that should be cleaned up
	// 		if ( std::abs(u21z) < 1.0 ) {
	// 		}
	// 		to_degrees( sig );
	// 	}
}


/// @brief load phi/theta bins for use in secondary structure scoring

void
SecondaryStructurePotential::load_phi_theta_bins(
	std::string const & hs_filename,
	std::string const & ss_filename
)
{
	using ObjexxFCL::fmt::skip;

	// local
	int isep,iph,itheta,idot;
	Real totn,tot;

	typedef ObjexxFCL::FArray3D< Real > FArray3D_real;
	FArray3D_real pts_HS( 36, 36, 3 );
	FArray3D_real pts_SS( 36, 36, 3 );
	FArray1D_real pts_correct( 3 );
//------------------------------------------------------------------------------

	for ( itheta = 1; itheta <= 36; ++itheta ) {
//        radians = pi/180.0;
//        angle1 = (itheta+17)*5*radians;
//        angle2 = (itheta+18)*5*radians;
//        std::cout << SS( angle1 ) << SS( angle2 ) << std::endl;
//        iptsn_(itheta) = 10000*std::abs(std::sin(angle1)-std::sin(angle2));
//        std::cout << SS( itheta ) << SS( iptsn_(itheta) ) << std::endl;
//  FOR PRESMOOTHED/BINNED COUNTS FROM INGO USE NO THETA NORM
		iptsn_(itheta) = 100;
	}

	// FIXME: need equivalent to open_data_file() function here
	//utility::io::izstream & HS_stream( open_data_file( hs_filename ) );
	utility::io::izstream HS_stream;
	basic::database::open( HS_stream, hs_filename );
	for ( isep = 1; isep <= 3; ++isep ) {
		for ( itheta = 1; itheta <= 36; ++itheta ) {
			for ( iph = 1; iph <= 36; ++iph ) {
				HS_stream >> pts_HS(itheta,iph,isep) >> skip;
			}
		}
	}
	HS_stream.close();
	HS_stream.clear();

	for ( isep = 1; isep <= 3; ++isep ) {
		tot = 0.0;
		totn = 0.0;
		for ( iph = 1; iph <= 36; ++iph ) {
			for ( itheta = 1; itheta <= 36; ++itheta ) {
				//  SMALL COUNTS CORRECTION
				pts_HS(itheta,iph,isep) += iptsn_(itheta)*0.000001f;
				tot += pts_HS(itheta,iph,isep);
				totn += iptsn_(itheta);
			}
		}
		for ( iph = 1; iph <= 36; ++iph ) {
			for ( itheta = 1; itheta <= 36; ++itheta ) {
				pts_(1,isep,iph,itheta) = -std::log(pts_HS(itheta,iph,isep)/tot) +
				 std::log(iptsn_(itheta)/totn);
			}
		}
	}

	// FIXME: need equivalent to open_data_file() function here
	utility::io::izstream SS_stream;
	basic::database::open( SS_stream, ss_filename );
	for ( isep = 1; isep <= 3; ++isep ) {
		for ( itheta = 1; itheta <= 36; ++itheta ) {
			for ( iph = 1; iph <= 36; ++iph ) {
				SS_stream >> pts_SS(itheta,iph,isep) >> skip;
			}
		}
		if ( isep == 1 ) SS_stream.seek_beg();
	}
	SS_stream.close();
	SS_stream.clear();

	for ( isep = 1; isep <= 3; ++isep ) {
		tot = 0.0;
		totn = 0.0;
		for ( iph = 1; iph <= 36; ++iph ) {
			for ( itheta = 1; itheta <= 36; ++itheta ) {
				//  SMALL COUNTS CORRECTION
				pts_SS(itheta,iph,isep) += iptsn_(itheta)*0.000001f;
				tot += pts_SS(itheta,iph,isep);
				totn += iptsn_(itheta);
			}
		}
		for ( iph = 1; iph <= 36; ++iph ) {
			for ( itheta = 1; itheta <= 36; ++itheta ) {
				pts_(2,isep,iph,itheta) = -std::log(pts_SS(itheta,iph,isep)/tot) +
				 std::log(iptsn_(itheta)/totn);
			}
		}
	}

	tot = 0.0;
	totn = 0.0;
	for ( idot = 1; idot <= 6; ++idot ) {
		tot += ids_(idot);
		totn += idsn_(idot);
	}
	for ( idot = 1; idot <= 6; ++idot ) {
		if ( ids_(idot) != 0 ) {
			ds_(idot) = -std::log(ids_(idot)/tot) + std::log(idsn_(idot)/totn);
		} else {
			ds_(idot) = 0.0;
		}
	}

	for ( isep = 1; isep <= 3; ++isep ) {
		pts_correct(isep) = 0.0;
		for ( iph = 1; iph <= 36; ++iph ) {
			for ( itheta = 1; itheta <= 36; ++itheta ) {
				if ( pts_(2,isep,iph,itheta) > pts_correct(isep) )
				 pts_correct(isep) = pts_(2,isep,iph,itheta);
			}
		}
	}
}


void
SecondaryStructurePotential::idsn_initializer(
	FArray1D_int & idsn
)
{
	// triangle-2 random numbers
	//     data idsn/56,167,278,278,167,56/
	// sort of triangle-4 random numbers
	int i = 0;
	idsn( ++i ) = 5596;
	idsn( ++i ) = 16581;
	idsn( ++i ) = 27823;
	idsn( ++i ) = 27823;
	idsn( ++i ) = 16581;
	idsn( ++i ) = 5596;
}


void
SecondaryStructurePotential::ids_initializer(
	FArray1D_int & ids
)
{
	int i = 0;
	ids( ++i ) = 1;
	ids( ++i ) = 48;
	ids( ++i ) = 368;
	ids( ++i ) = 2378;
	ids( ++i ) = 7141;
	ids( ++i ) = 8904;
}


void
SecondaryStructurePotential::ssdist_initializer(
	FArray2D_real & ssdist
)
{
	//  DATA 0-12A, GOES BACK TO PROTEINS 1999 PAPER
	//     data ssdist/2.652527,0.7284873,0.0176830,-0.2566608,
	//    #           -1.471609,0.0104174,0.0679096, 0.4667910/

	//  DATA 0-12A+, CREATED JANUARY 30, 1999
	FArray1A_real ssdist1d( ssdist, ssdist.size() ); // 1D view
	int i = 0;
	ssdist1d( ++i ) = 2.3962;
	ssdist1d( ++i ) = 0.56921;
	ssdist1d( ++i ) = -0.20262;
	ssdist1d( ++i ) = -0.55172;
	ssdist1d( ++i ) = -1.6408;
	ssdist1d( ++i ) = -0.63196;
	ssdist1d( ++i ) = -0.57115;
	ssdist1d( ++i ) = -0.26221;
}


void
SecondaryStructurePotential::hs_dp_initializer(
	FArray1D_real & hs_dp
)
{
	int i = 0;
	hs_dp( ++i ) = 0.416;
	hs_dp( ++i ) = -0.412;
	hs_dp( ++i ) = -0.542;
	hs_dp( ++i ) = -0.489;
	hs_dp( ++i ) = -0.351;
	hs_dp( ++i ) = -0.104;
	hs_dp( ++i ) = 0.211;
	hs_dp( ++i ) = 0.494;
	hs_dp( ++i ) = 0.942;
	hs_dp( ++i ) = 1.897;
}


void
SecondaryStructurePotential::rsigma_dot_initializer(
	FArray4D_real & rsigma_dot
)
{
	// section 12
	//js --------------------------
	//js new rsigma stats that take into account whether the first
	//js c=o bond vector points towards away from the other dimer
	rsigma_dot(  1,  1,  1,  2 ) = -1.038100;
	rsigma_dot(  1,  2,  1,  2 ) = -1.038100;
	rsigma_dot(  1,  3,  1,  2 ) = -1.038100;
	rsigma_dot(  1,  4,  1,  2 ) = -1.038100;
	rsigma_dot(  1,  5,  1,  2 ) = -1.038100;
	rsigma_dot(  1,  6,  1,  2 ) = -1.038100;
	rsigma_dot(  1,  7,  1,  2 ) = -1.038100;
	rsigma_dot(  1,  8,  1,  2 ) = -2.984000;
	rsigma_dot(  1,  9,  1,  2 ) = -3.746110;
	rsigma_dot(  1, 10,  1,  2 ) = -3.746110;
	rsigma_dot(  1, 11,  1,  2 ) = -2.647500;
	rsigma_dot(  1, 12,  1,  2 ) = -1.038100;
	rsigma_dot(  1, 13,  1,  2 ) = -1.038100;
	rsigma_dot(  1, 14,  1,  2 ) = -1.038100;
	rsigma_dot(  1, 15,  1,  2 ) = -1.038100;
	rsigma_dot(  1, 16,  1,  2 ) = -1.038100;
	rsigma_dot(  1, 17,  1,  2 ) = -1.038100;
	rsigma_dot(  1, 18,  1,  2 ) = -1.038100;
	rsigma_dot(  2,  1,  1,  2 ) = -0.973500;
	rsigma_dot(  2,  2,  1,  2 ) = -0.973500;
	rsigma_dot(  2,  3,  1,  2 ) = -0.973500;
	rsigma_dot(  2,  4,  1,  2 ) = -0.973500;
	rsigma_dot(  2,  5,  1,  2 ) = -0.973500;
	rsigma_dot(  2,  6,  1,  2 ) = -0.973500;
	rsigma_dot(  2,  7,  1,  2 ) = -2.072100;
	rsigma_dot(  2,  8,  1,  2 ) = -0.973500;
	rsigma_dot(  2,  9,  1,  2 ) = -4.943820;
	rsigma_dot(  2, 10,  1,  2 ) = -5.178220;
	rsigma_dot(  2, 11,  1,  2 ) = -4.823670;
	rsigma_dot(  2, 12,  1,  2 ) = -2.919400;
	rsigma_dot(  2, 13,  1,  2 ) = -2.583000;
	rsigma_dot(  2, 14,  1,  2 ) = -0.973500;
	rsigma_dot(  2, 15,  1,  2 ) = -0.973500;
	rsigma_dot(  2, 16,  1,  2 ) = -0.973500;
	rsigma_dot(  2, 17,  1,  2 ) = -0.973500;
	rsigma_dot(  2, 18,  1,  2 ) = -0.973500;
	rsigma_dot(  3,  1,  1,  2 ) = -0.912900;
	rsigma_dot(  3,  2,  1,  2 ) = -0.912900;
	rsigma_dot(  3,  3,  1,  2 ) = -0.912900;
	rsigma_dot(  3,  4,  1,  2 ) = -0.912900;
	rsigma_dot(  3,  5,  1,  2 ) = -0.912900;
	rsigma_dot(  3,  6,  1,  2 ) = -0.912900;
	rsigma_dot(  3,  7,  1,  2 ) = -2.011500;
	rsigma_dot(  3,  8,  1,  2 ) = -4.048400;
	rsigma_dot(  3,  9,  1,  2 ) = -5.981810;
	rsigma_dot(  3, 10,  1,  2 ) = -7.006470;
	rsigma_dot(  3, 11,  1,  2 ) = -6.696730;
	rsigma_dot(  3, 12,  1,  2 ) = -4.131780;
	rsigma_dot(  3, 13,  1,  2 ) = -0.912900;
	rsigma_dot(  3, 14,  1,  2 ) = -0.912900;
	rsigma_dot(  3, 15,  1,  2 ) = -0.912900;
	rsigma_dot(  3, 16,  1,  2 ) = -0.912900;
	rsigma_dot(  3, 17,  1,  2 ) = -0.912900;
	rsigma_dot(  3, 18,  1,  2 ) = -0.912900;
	rsigma_dot(  4,  1,  1,  2 ) = -0.855700;
	rsigma_dot(  4,  2,  1,  2 ) = -0.855700;
	rsigma_dot(  4,  3,  1,  2 ) = -0.855700;
	rsigma_dot(  4,  4,  1,  2 ) = -0.855700;
	rsigma_dot(  4,  5,  1,  2 ) = -0.855700;
	rsigma_dot(  4,  6,  1,  2 ) = -2.465200;
	rsigma_dot(  4,  7,  1,  2 ) = -3.420700;
	rsigma_dot(  4,  8,  1,  2 ) = -3.900270;
	rsigma_dot(  4,  9,  1,  2 ) = -6.315330;
	rsigma_dot(  4, 10,  1,  2 ) = -8.116270;
	rsigma_dot(  4, 11,  1,  2 ) = -8.329380;
	rsigma_dot(  4, 12,  1,  2 ) = -5.298390;
	rsigma_dot(  4, 13,  1,  2 ) = -1.954400;
	rsigma_dot(  4, 14,  1,  2 ) = -0.855700;
	rsigma_dot(  4, 15,  1,  2 ) = -0.855700;
	rsigma_dot(  4, 16,  1,  2 ) = -0.855700;
	rsigma_dot(  4, 17,  1,  2 ) = -0.855700;
	rsigma_dot(  4, 18,  1,  2 ) = -0.855700;
	rsigma_dot(  5,  1,  1,  2 ) = -0.801700;
	rsigma_dot(  5,  2,  1,  2 ) = -0.801700;
	rsigma_dot(  5,  3,  1,  2 ) = -0.801700;
	rsigma_dot(  5,  4,  1,  2 ) = -0.801700;
	rsigma_dot(  5,  5,  1,  2 ) = -0.801700;
	rsigma_dot(  5,  6,  1,  2 ) = -1.900300;
	rsigma_dot(  5,  7,  1,  2 ) = -3.634900;
	rsigma_dot(  5,  8,  1,  2 ) = -5.119160;
	rsigma_dot(  5,  9,  1,  2 ) = -6.381410;
	rsigma_dot(  5, 10,  1,  2 ) = -9.110860;
	rsigma_dot(  5, 11,  1,  2 ) = -9.615860;
	rsigma_dot(  5, 12,  1,  2 ) = -6.566870;
	rsigma_dot(  5, 13,  1,  2 ) = -2.411100;
	rsigma_dot(  5, 14,  1,  2 ) = -1.900300;
	rsigma_dot(  5, 15,  1,  2 ) = -0.801700;
	rsigma_dot(  5, 16,  1,  2 ) = -0.801700;
	rsigma_dot(  5, 17,  1,  2 ) = -0.801700;
	rsigma_dot(  5, 18,  1,  2 ) = -0.801700;
	rsigma_dot(  6,  1,  1,  2 ) = -0.750400;
	rsigma_dot(  6,  2,  1,  2 ) = -0.750400;
	rsigma_dot(  6,  3,  1,  2 ) = -0.750400;
	rsigma_dot(  6,  4,  1,  2 ) = -0.750400;
	rsigma_dot(  6,  5,  1,  2 ) = -0.750400;
	rsigma_dot(  6,  6,  1,  2 ) = -2.359800;
	rsigma_dot(  6,  7,  1,  2 ) = -4.511580;
	rsigma_dot(  6,  8,  1,  2 ) = -5.325090;
	rsigma_dot(  6,  9,  1,  2 ) = -6.768980;
	rsigma_dot(  6, 10,  1,  2 ) = -8.613260;
	rsigma_dot(  6, 11,  1,  2 ) = -9.405770;
	rsigma_dot(  6, 12,  1,  2 ) = -7.032650;
	rsigma_dot(  6, 13,  1,  2 ) = -3.885880;
	rsigma_dot(  6, 14,  1,  2 ) = -3.315300;
	rsigma_dot(  6, 15,  1,  2 ) = -0.750400;
	rsigma_dot(  6, 16,  1,  2 ) = -0.750400;
	rsigma_dot(  6, 17,  1,  2 ) = -0.750400;
	rsigma_dot(  6, 18,  1,  2 ) = -0.750400;
	rsigma_dot(  7,  1,  1,  2 ) = -0.701600;
	rsigma_dot(  7,  2,  1,  2 ) = -0.701600;
	rsigma_dot(  7,  3,  1,  2 ) = -0.701600;
	rsigma_dot(  7,  4,  1,  2 ) = -0.701600;
	rsigma_dot(  7,  5,  1,  2 ) = -0.701600;
	rsigma_dot(  7,  6,  1,  2 ) = -2.311000;
	rsigma_dot(  7,  7,  1,  2 ) = -4.365150;
	rsigma_dot(  7,  8,  1,  2 ) = -5.446520;
	rsigma_dot(  7,  9,  1,  2 ) = -6.795160;
	rsigma_dot(  7, 10,  1,  2 ) = -7.145720;
	rsigma_dot(  7, 11,  1,  2 ) = -7.527050;
	rsigma_dot(  7, 12,  1,  2 ) = -6.885740;
	rsigma_dot(  7, 13,  1,  2 ) = -4.508250;
	rsigma_dot(  7, 14,  1,  2 ) = -1.800200;
	rsigma_dot(  7, 15,  1,  2 ) = -2.647500;
	rsigma_dot(  7, 16,  1,  2 ) = -0.701600;
	rsigma_dot(  7, 17,  1,  2 ) = -0.701600;
	rsigma_dot(  7, 18,  1,  2 ) = -0.701600;
	rsigma_dot(  8,  1,  1,  2 ) = -0.655100;
	rsigma_dot(  8,  2,  1,  2 ) = -0.655100;
	rsigma_dot(  8,  3,  1,  2 ) = -0.655100;
	rsigma_dot(  8,  4,  1,  2 ) = -0.655100;
	rsigma_dot(  8,  5,  1,  2 ) = -0.655100;
	rsigma_dot(  8,  6,  1,  2 ) = -3.699600;
	rsigma_dot(  8,  7,  1,  2 ) = -4.089060;
	rsigma_dot(  8,  8,  1,  2 ) = -5.417250;
	rsigma_dot(  8,  9,  1,  2 ) = -6.293430;
	rsigma_dot(  8, 10,  1,  2 ) = -6.663890;
	rsigma_dot(  8, 11,  1,  2 ) = -6.788470;
	rsigma_dot(  8, 12,  1,  2 ) = -6.368800;
	rsigma_dot(  8, 13,  1,  2 ) = -4.732610;
	rsigma_dot(  8, 14,  1,  2 ) = -1.753700;
	rsigma_dot(  8, 15,  1,  2 ) = -0.655100;
	rsigma_dot(  8, 16,  1,  2 ) = -0.655100;
	rsigma_dot(  8, 17,  1,  2 ) = -0.655100;
	rsigma_dot(  8, 18,  1,  2 ) = -0.655100;
	rsigma_dot(  9,  1,  1,  2 ) = -0.610600;
	rsigma_dot(  9,  2,  1,  2 ) = -0.610600;
	rsigma_dot(  9,  3,  1,  2 ) = -0.610600;
	rsigma_dot(  9,  4,  1,  2 ) = -0.610600;
	rsigma_dot(  9,  5,  1,  2 ) = -1.709200;
	rsigma_dot(  9,  6,  1,  2 ) = -3.906460;
	rsigma_dot(  9,  7,  1,  2 ) = -3.906460;
	rsigma_dot(  9,  8,  1,  2 ) = -5.355550;
	rsigma_dot(  9,  9,  1,  2 ) = -6.070210;
	rsigma_dot(  9, 10,  1,  2 ) = -5.971910;
	rsigma_dot(  9, 11,  1,  2 ) = -6.070210;
	rsigma_dot(  9, 12,  1,  2 ) = -5.852370;
	rsigma_dot(  9, 13,  1,  2 ) = -4.901080;
	rsigma_dot(  9, 14,  1,  2 ) = -2.807800;
	rsigma_dot(  9, 15,  1,  2 ) = -2.220100;
	rsigma_dot(  9, 16,  1,  2 ) = -0.610600;
	rsigma_dot(  9, 17,  1,  2 ) = -1.709200;
	rsigma_dot(  9, 18,  1,  2 ) = -0.610600;
	rsigma_dot( 10,  1,  1,  2 ) = -0.568100;
	rsigma_dot( 10,  2,  1,  2 ) = -0.568100;
	rsigma_dot( 10,  3,  1,  2 ) = -0.568100;
	rsigma_dot( 10,  4,  1,  2 ) = -2.177500;
	rsigma_dot( 10,  5,  1,  2 ) = -2.765300;
	rsigma_dot( 10,  6,  1,  2 ) = -4.002050;
	rsigma_dot( 10,  7,  1,  2 ) = -4.802170;
	rsigma_dot( 10,  8,  1,  2 ) = -5.598500;
	rsigma_dot( 10,  9,  1,  2 ) = -5.809810;
	rsigma_dot( 10, 10,  1,  2 ) = -5.809810;
	rsigma_dot( 10, 11,  1,  2 ) = -5.732850;
	rsigma_dot( 10, 12,  1,  2 ) = -5.572010;
	rsigma_dot( 10, 13,  1,  2 ) = -4.858520;
	rsigma_dot( 10, 14,  1,  2 ) = -3.935360;
	rsigma_dot( 10, 15,  1,  2 ) = -1.666700;
	rsigma_dot( 10, 16,  1,  2 ) = -1.666700;
	rsigma_dot( 10, 17,  1,  2 ) = -0.568100;
	rsigma_dot( 10, 18,  1,  2 ) = -1.666700;
	rsigma_dot( 11,  1,  1,  2 ) = -0.527200;
	rsigma_dot( 11,  2,  1,  2 ) = -0.527200;
	rsigma_dot( 11,  3,  1,  2 ) = -0.527200;
	rsigma_dot( 11,  4,  1,  2 ) = -1.625900;
	rsigma_dot( 11,  5,  1,  2 ) = -3.746110;
	rsigma_dot( 11,  6,  1,  2 ) = -4.190800;
	rsigma_dot( 11,  7,  1,  2 ) = -4.761350;
	rsigma_dot( 11,  8,  1,  2 ) = -5.897880;
	rsigma_dot( 11,  9,  1,  2 ) = -5.850250;
	rsigma_dot( 11, 10,  1,  2 ) = -5.789930;
	rsigma_dot( 11, 11,  1,  2 ) = -5.371430;
	rsigma_dot( 11, 12,  1,  2 ) = -5.608640;
	rsigma_dot( 11, 13,  1,  2 ) = -5.306360;
	rsigma_dot( 11, 14,  1,  2 ) = -4.288440;
	rsigma_dot( 11, 15,  1,  2 ) = -3.092200;
	rsigma_dot( 11, 16,  1,  2 ) = -0.527200;
	rsigma_dot( 11, 17,  1,  2 ) = -0.527200;
	rsigma_dot( 11, 18,  1,  2 ) = -0.527200;
	rsigma_dot( 12,  1,  1,  2 ) = -0.488000;
	rsigma_dot( 12,  2,  1,  2 ) = -0.488000;
	rsigma_dot( 12,  3,  1,  2 ) = -0.488000;
	rsigma_dot( 12,  4,  1,  2 ) = -3.623500;
	rsigma_dot( 12,  5,  1,  2 ) = -4.722120;
	rsigma_dot( 12,  6,  1,  2 ) = -4.692710;
	rsigma_dot( 12,  7,  1,  2 ) = -5.197550;
	rsigma_dot( 12,  8,  1,  2 ) = -5.544260;
	rsigma_dot( 12,  9,  1,  2 ) = -5.740290;
	rsigma_dot( 12, 10,  1,  2 ) = -5.478450;
	rsigma_dot( 12, 11,  1,  2 ) = -5.347830;
	rsigma_dot( 12, 12,  1,  2 ) = -5.363220;
	rsigma_dot( 12, 13,  1,  2 ) = -5.103140;
	rsigma_dot( 12, 14,  1,  2 ) = -4.722120;
	rsigma_dot( 12, 15,  1,  2 ) = -3.984530;
	rsigma_dot( 12, 16,  1,  2 ) = -3.196100;
	rsigma_dot( 12, 17,  1,  2 ) = -0.488000;
	rsigma_dot( 12, 18,  1,  2 ) = -0.488000;

	// section 11
	rsigma_dot(  1,  1,  1,  1 ) = -0.552300;
	rsigma_dot(  1,  2,  1,  1 ) = -0.552300;
	rsigma_dot(  1,  3,  1,  1 ) = -0.552300;
	rsigma_dot(  1,  4,  1,  1 ) = -0.552300;
	rsigma_dot(  1,  5,  1,  1 ) = -0.552300;
	rsigma_dot(  1,  6,  1,  1 ) = -0.552300;
	rsigma_dot(  1,  7,  1,  1 ) = -0.552300;
	rsigma_dot(  1,  8,  1,  1 ) = -0.552300;
	rsigma_dot(  1,  9,  1,  1 ) = -0.552300;
	rsigma_dot(  1, 10,  1,  1 ) = -1.651000;
	rsigma_dot(  1, 11,  1,  1 ) = -0.552300;
	rsigma_dot(  1, 12,  1,  1 ) = -0.552300;
	rsigma_dot(  1, 13,  1,  1 ) = -0.552300;
	rsigma_dot(  1, 14,  1,  1 ) = -0.552300;
	rsigma_dot(  1, 15,  1,  1 ) = -0.552300;
	rsigma_dot(  1, 16,  1,  1 ) = -0.552300;
	rsigma_dot(  1, 17,  1,  1 ) = -0.552300;
	rsigma_dot(  1, 18,  1,  1 ) = -0.552300;
	rsigma_dot(  2,  1,  1,  1 ) = -0.487800;
	rsigma_dot(  2,  2,  1,  1 ) = -0.487800;
	rsigma_dot(  2,  3,  1,  1 ) = -0.487800;
	rsigma_dot(  2,  4,  1,  1 ) = -0.487800;
	rsigma_dot(  2,  5,  1,  1 ) = -0.487800;
	rsigma_dot(  2,  6,  1,  1 ) = -0.487800;
	rsigma_dot(  2,  7,  1,  1 ) = -0.487800;
	rsigma_dot(  2,  8,  1,  1 ) = -0.487800;
	rsigma_dot(  2,  9,  1,  1 ) = -0.487800;
	rsigma_dot(  2, 10,  1,  1 ) = -1.586400;
	rsigma_dot(  2, 11,  1,  1 ) = -2.097200;
	rsigma_dot(  2, 12,  1,  1 ) = -1.586400;
	rsigma_dot(  2, 13,  1,  1 ) = -0.487800;
	rsigma_dot(  2, 14,  1,  1 ) = -0.487800;
	rsigma_dot(  2, 15,  1,  1 ) = -0.487800;
	rsigma_dot(  2, 16,  1,  1 ) = -0.487800;
	rsigma_dot(  2, 17,  1,  1 ) = -0.487800;
	rsigma_dot(  2, 18,  1,  1 ) = -0.487800;
	rsigma_dot(  3,  1,  1,  1 ) = -0.427200;
	rsigma_dot(  3,  2,  1,  1 ) = -0.427200;
	rsigma_dot(  3,  3,  1,  1 ) = -0.427200;
	rsigma_dot(  3,  4,  1,  1 ) = -0.427200;
	rsigma_dot(  3,  5,  1,  1 ) = -0.427200;
	rsigma_dot(  3,  6,  1,  1 ) = -1.525800;
	rsigma_dot(  3,  7,  1,  1 ) = -2.036600;
	rsigma_dot(  3,  8,  1,  1 ) = -2.036600;
	rsigma_dot(  3,  9,  1,  1 ) = -0.427200;
	rsigma_dot(  3, 10,  1,  1 ) = -2.373100;
	rsigma_dot(  3, 11,  1,  1 ) = -1.525800;
	rsigma_dot(  3, 12,  1,  1 ) = -2.992100;
	rsigma_dot(  3, 13,  1,  1 ) = -2.373100;
	rsigma_dot(  3, 14,  1,  1 ) = -0.427200;
	rsigma_dot(  3, 15,  1,  1 ) = -0.427200;
	rsigma_dot(  3, 16,  1,  1 ) = -0.427200;
	rsigma_dot(  3, 17,  1,  1 ) = -0.427200;
	rsigma_dot(  3, 18,  1,  1 ) = -0.427200;
	rsigma_dot(  4,  1,  1,  1 ) = -0.370000;
	rsigma_dot(  4,  2,  1,  1 ) = -0.370000;
	rsigma_dot(  4,  3,  1,  1 ) = -0.370000;
	rsigma_dot(  4,  4,  1,  1 ) = -0.370000;
	rsigma_dot(  4,  5,  1,  1 ) = -0.370000;
	rsigma_dot(  4,  6,  1,  1 ) = -0.370000;
	rsigma_dot(  4,  7,  1,  1 ) = -2.315900;
	rsigma_dot(  4,  8,  1,  1 ) = -2.567200;
	rsigma_dot(  4,  9,  1,  1 ) = -2.767900;
	rsigma_dot(  4, 10,  1,  1 ) = -3.588900;
	rsigma_dot(  4, 11,  1,  1 ) = -3.414500;
	rsigma_dot(  4, 12,  1,  1 ) = -4.033580;
	rsigma_dot(  4, 13,  1,  1 ) = -3.804010;
	rsigma_dot(  4, 14,  1,  1 ) = -0.370000;
	rsigma_dot(  4, 15,  1,  1 ) = -0.370000;
	rsigma_dot(  4, 16,  1,  1 ) = -0.370000;
	rsigma_dot(  4, 17,  1,  1 ) = -0.370000;
	rsigma_dot(  4, 18,  1,  1 ) = -0.370000;
	rsigma_dot(  5,  1,  1,  1 ) = -0.316000;
	rsigma_dot(  5,  2,  1,  1 ) = -0.316000;
	rsigma_dot(  5,  3,  1,  1 ) = -0.316000;
	rsigma_dot(  5,  4,  1,  1 ) = -0.316000;
	rsigma_dot(  5,  5,  1,  1 ) = -0.316000;
	rsigma_dot(  5,  6,  1,  1 ) = -1.925400;
	rsigma_dot(  5,  7,  1,  1 ) = -3.683200;
	rsigma_dot(  5,  8,  1,  1 ) = -4.122610;
	rsigma_dot(  5,  9,  1,  1 ) = -3.812460;
	rsigma_dot(  5, 10,  1,  1 ) = -3.683200;
	rsigma_dot(  5, 11,  1,  1 ) = -3.360500;
	rsigma_dot(  5, 12,  1,  1 ) = -4.969910;
	rsigma_dot(  5, 13,  1,  1 ) = -4.869830;
	rsigma_dot(  5, 14,  1,  1 ) = -3.149200;
	rsigma_dot(  5, 15,  1,  1 ) = -0.316000;
	rsigma_dot(  5, 16,  1,  1 ) = -0.316000;
	rsigma_dot(  5, 17,  1,  1 ) = -0.316000;
	rsigma_dot(  5, 18,  1,  1 ) = -0.316000;
	rsigma_dot(  6,  1,  1,  1 ) = -0.264700;
	rsigma_dot(  6,  2,  1,  1 ) = -0.264700;
	rsigma_dot(  6,  3,  1,  1 ) = -0.264700;
	rsigma_dot(  6,  4,  1,  1 ) = -0.264700;
	rsigma_dot(  6,  5,  1,  1 ) = -0.264700;
	rsigma_dot(  6,  6,  1,  1 ) = -2.829600;
	rsigma_dot(  6,  7,  1,  1 ) = -5.680760;
	rsigma_dot(  6,  8,  1,  1 ) = -5.076840;
	rsigma_dot(  6,  9,  1,  1 ) = -4.271990;
	rsigma_dot(  6, 10,  1,  1 ) = -4.234950;
	rsigma_dot(  6, 11,  1,  1 ) = -3.761170;
	rsigma_dot(  6, 12,  1,  1 ) = -5.567960;
	rsigma_dot(  6, 13,  1,  1 ) = -6.108200;
	rsigma_dot(  6, 14,  1,  1 ) = -4.956010;
	rsigma_dot(  6, 15,  1,  1 ) = -1.874100;
	rsigma_dot(  6, 16,  1,  1 ) = -0.264700;
	rsigma_dot(  6, 17,  1,  1 ) = -0.264700;
	rsigma_dot(  6, 18,  1,  1 ) = -0.264700;
	rsigma_dot(  7,  1,  1,  1 ) = -0.215900;
	rsigma_dot(  7,  2,  1,  1 ) = -0.215900;
	rsigma_dot(  7,  3,  1,  1 ) = -0.215900;
	rsigma_dot(  7,  4,  1,  1 ) = -0.215900;
	rsigma_dot(  7,  5,  1,  1 ) = -0.215900;
	rsigma_dot(  7,  6,  1,  1 ) = -5.011660;
	rsigma_dot(  7,  7,  1,  1 ) = -7.313420;
	rsigma_dot(  7,  8,  1,  1 ) = -6.580620;
	rsigma_dot(  7,  9,  1,  1 ) = -4.726730;
	rsigma_dot(  7, 10,  1,  1 ) = -4.748470;
	rsigma_dot(  7, 11,  1,  1 ) = -2.923900;
	rsigma_dot(  7, 12,  1,  1 ) = -5.700670;
	rsigma_dot(  7, 13,  1,  1 ) = -6.799280;
	rsigma_dot(  7, 14,  1,  1 ) = -6.059410;
	rsigma_dot(  7, 15,  1,  1 ) = -2.161800;
	rsigma_dot(  7, 16,  1,  1 ) = -0.215900;
	rsigma_dot(  7, 17,  1,  1 ) = -0.215900;
	rsigma_dot(  7, 18,  1,  1 ) = -0.215900;
	rsigma_dot(  8,  1,  1,  1 ) = -0.169300;
	rsigma_dot(  8,  2,  1,  1 ) = -0.169300;
	rsigma_dot(  8,  3,  1,  1 ) = -0.169300;
	rsigma_dot(  8,  4,  1,  1 ) = -0.169300;
	rsigma_dot(  8,  5,  1,  1 ) = -0.169300;
	rsigma_dot(  8,  6,  1,  1 ) = -6.132930;
	rsigma_dot(  8,  7,  1,  1 ) = -8.735140;
	rsigma_dot(  8,  8,  1,  1 ) = -7.789560;
	rsigma_dot(  8,  9,  1,  1 ) = -5.199790;
	rsigma_dot(  8, 10,  1,  1 ) = -4.981530;
	rsigma_dot(  8, 11,  1,  1 ) = -3.882920;
	rsigma_dot(  8, 12,  1,  1 ) = -5.059700;
	rsigma_dot(  8, 13,  1,  1 ) = -7.329420;
	rsigma_dot(  8, 14,  1,  1 ) = -7.217740;
	rsigma_dot(  8, 15,  1,  1 ) = -4.723230;
	rsigma_dot(  8, 16,  1,  1 ) = -0.169300;
	rsigma_dot(  8, 17,  1,  1 ) = -0.169300;
	rsigma_dot(  8, 18,  1,  1 ) = -0.169300;
	rsigma_dot(  9,  1,  1,  1 ) = -0.124900;
	rsigma_dot(  9,  2,  1,  1 ) = -0.124900;
	rsigma_dot(  9,  3,  1,  1 ) = -0.124900;
	rsigma_dot(  9,  4,  1,  1 ) = -1.223500;
	rsigma_dot(  9,  5,  1,  1 ) = -3.069300;
	rsigma_dot(  9,  6,  1,  1 ) = -6.756900;
	rsigma_dot(  9,  7,  1,  1 ) = -8.997660;
	rsigma_dot(  9,  8,  1,  1 ) = -7.292710;
	rsigma_dot(  9,  9,  1,  1 ) = -5.704630;
	rsigma_dot(  9, 10,  1,  1 ) = -4.442390;
	rsigma_dot(  9, 11,  1,  1 ) = -4.202430;
	rsigma_dot(  9, 12,  1,  1 ) = -4.567550;
	rsigma_dot(  9, 13,  1,  1 ) = -7.550850;
	rsigma_dot(  9, 14,  1,  1 ) = -8.291970;
	rsigma_dot(  9, 15,  1,  1 ) = -5.617960;
	rsigma_dot(  9, 16,  1,  1 ) = -1.223500;
	rsigma_dot(  9, 17,  1,  1 ) = -1.223500;
	rsigma_dot(  9, 18,  1,  1 ) = -0.124900;
	rsigma_dot( 10,  1,  1,  1 ) = -0.082300;
	rsigma_dot( 10,  2,  1,  1 ) = -0.082300;
	rsigma_dot( 10,  3,  1,  1 ) = -0.082300;
	rsigma_dot( 10,  4,  1,  1 ) = -1.180900;
	rsigma_dot( 10,  5,  1,  1 ) = -3.932490;
	rsigma_dot( 10,  6,  1,  1 ) = -6.858840;
	rsigma_dot( 10,  7,  1,  1 ) = -8.186740;
	rsigma_dot( 10,  8,  1,  1 ) = -6.228670;
	rsigma_dot( 10,  9,  1,  1 ) = -5.720690;
	rsigma_dot( 10, 10,  1,  1 ) = -4.548250;
	rsigma_dot( 10, 11,  1,  1 ) = -4.451790;
	rsigma_dot( 10, 12,  1,  1 ) = -4.548250;
	rsigma_dot( 10, 13,  1,  1 ) = -7.242410;
	rsigma_dot( 10, 14,  1,  1 ) = -8.619920;
	rsigma_dot( 10, 15,  1,  1 ) = -6.349540;
	rsigma_dot( 10, 16,  1,  1 ) = -3.026800;
	rsigma_dot( 10, 17,  1,  1 ) = -1.180900;
	rsigma_dot( 10, 18,  1,  1 ) = -1.180900;
	rsigma_dot( 11,  1,  1,  1 ) = -0.041500;
	rsigma_dot( 11,  2,  1,  1 ) = -0.041500;
	rsigma_dot( 11,  3,  1,  1 ) = -0.041500;
	rsigma_dot( 11,  4,  1,  1 ) = -1.987400;
	rsigma_dot( 11,  5,  1,  1 ) = -4.275620;
	rsigma_dot( 11,  6,  1,  1 ) = -6.582550;
	rsigma_dot( 11,  7,  1,  1 ) = -6.915710;
	rsigma_dot( 11,  8,  1,  1 ) = -6.166200;
	rsigma_dot( 11,  9,  1,  1 ) = -5.483930;
	rsigma_dot( 11, 10,  1,  1 ) = -4.595390;
	rsigma_dot( 11, 11,  1,  1 ) = -4.616230;
	rsigma_dot( 11, 12,  1,  1 ) = -4.931860;
	rsigma_dot( 11, 13,  1,  1 ) = -6.246070;
	rsigma_dot( 11, 14,  1,  1 ) = -8.296560;
	rsigma_dot( 11, 15,  1,  1 ) = -6.734840;
	rsigma_dot( 11, 16,  1,  1 ) = -3.337400;
	rsigma_dot( 11, 17,  1,  1 ) = -1.140100;
	rsigma_dot( 11, 18,  1,  1 ) = -0.041500;
	rsigma_dot( 12,  1,  1,  1 ) = -0.002300;
	rsigma_dot( 12,  2,  1,  1 ) = -0.002300;
	rsigma_dot( 12,  3,  1,  1 ) = -0.002300;
	rsigma_dot( 12,  4,  1,  1 ) = -2.199500;
	rsigma_dot( 12,  5,  1,  1 ) = -4.729680;
	rsigma_dot( 12,  6,  1,  1 ) = -6.131350;
	rsigma_dot( 12,  7,  1,  1 ) = -5.729140;
	rsigma_dot( 12,  8,  1,  1 ) = -5.604410;
	rsigma_dot( 12,  9,  1,  1 ) = -5.211780;
	rsigma_dot( 12, 10,  1,  1 ) = -4.892640;
	rsigma_dot( 12, 11,  1,  1 ) = -4.992730;
	rsigma_dot( 12, 12,  1,  1 ) = -5.275290;
	rsigma_dot( 12, 13,  1,  1 ) = -5.945090;
	rsigma_dot( 12, 14,  1,  1 ) = -7.716970;
	rsigma_dot( 12, 15,  1,  1 ) = -6.693140;
	rsigma_dot( 12, 16,  1,  1 ) = -3.972590;
	rsigma_dot( 12, 17,  1,  1 ) = -1.948200;
	rsigma_dot( 12, 18,  1,  1 ) = -0.002300;

	// section 22
	rsigma_dot(  1,  1,  2,  2 ) = -0.550000;
	rsigma_dot(  1,  2,  2,  2 ) = -0.550000;
	rsigma_dot(  1,  3,  2,  2 ) = -0.550000;
	rsigma_dot(  1,  4,  2,  2 ) = -0.550000;
	rsigma_dot(  1,  5,  2,  2 ) = -0.550000;
	rsigma_dot(  1,  6,  2,  2 ) = -0.550000;
	rsigma_dot(  1,  7,  2,  2 ) = -0.550000;
	rsigma_dot(  1,  8,  2,  2 ) = -2.159400;
	rsigma_dot(  1,  9,  2,  2 ) = -0.550000;
	rsigma_dot(  1, 10,  2,  2 ) = -2.159400;
	rsigma_dot(  1, 11,  2,  2 ) = -0.550000;
	rsigma_dot(  1, 12,  2,  2 ) = -0.550000;
	rsigma_dot(  1, 13,  2,  2 ) = -0.550000;
	rsigma_dot(  1, 14,  2,  2 ) = -0.550000;
	rsigma_dot(  1, 15,  2,  2 ) = -0.550000;
	rsigma_dot(  1, 16,  2,  2 ) = -0.550000;
	rsigma_dot(  1, 17,  2,  2 ) = -0.550000;
	rsigma_dot(  1, 18,  2,  2 ) = -0.550000;
	rsigma_dot(  2,  1,  2,  2 ) = -0.485500;
	rsigma_dot(  2,  2,  2,  2 ) = -0.485500;
	rsigma_dot(  2,  3,  2,  2 ) = -0.485500;
	rsigma_dot(  2,  4,  2,  2 ) = -0.485500;
	rsigma_dot(  2,  5,  2,  2 ) = -0.485500;
	rsigma_dot(  2,  6,  2,  2 ) = -1.584100;
	rsigma_dot(  2,  7,  2,  2 ) = -3.050400;
	rsigma_dot(  2,  8,  2,  2 ) = -2.094900;
	rsigma_dot(  2,  9,  2,  2 ) = -2.431400;
	rsigma_dot(  2, 10,  2,  2 ) = -1.584100;
	rsigma_dot(  2, 11,  2,  2 ) = -0.485500;
	rsigma_dot(  2, 12,  2,  2 ) = -0.485500;
	rsigma_dot(  2, 13,  2,  2 ) = -0.485500;
	rsigma_dot(  2, 14,  2,  2 ) = -0.485500;
	rsigma_dot(  2, 15,  2,  2 ) = -0.485500;
	rsigma_dot(  2, 16,  2,  2 ) = -0.485500;
	rsigma_dot(  2, 17,  2,  2 ) = -0.485500;
	rsigma_dot(  2, 18,  2,  2 ) = -0.485500;
	rsigma_dot(  3,  1,  2,  2 ) = -0.424800;
	rsigma_dot(  3,  2,  2,  2 ) = -0.424800;
	rsigma_dot(  3,  3,  2,  2 ) = -0.424800;
	rsigma_dot(  3,  4,  2,  2 ) = -0.424800;
	rsigma_dot(  3,  5,  2,  2 ) = -0.424800;
	rsigma_dot(  3,  6,  2,  2 ) = -3.258000;
	rsigma_dot(  3,  7,  2,  2 ) = -4.231500;
	rsigma_dot(  3,  8,  2,  2 ) = -3.980180;
	rsigma_dot(  3,  9,  2,  2 ) = -2.370700;
	rsigma_dot(  3, 10,  2,  2 ) = -1.523400;
	rsigma_dot(  3, 11,  2,  2 ) = -2.822700;
	rsigma_dot(  3, 12,  2,  2 ) = -0.424800;
	rsigma_dot(  3, 13,  2,  2 ) = -0.424800;
	rsigma_dot(  3, 14,  2,  2 ) = -0.424800;
	rsigma_dot(  3, 15,  2,  2 ) = -0.424800;
	rsigma_dot(  3, 16,  2,  2 ) = -0.424800;
	rsigma_dot(  3, 17,  2,  2 ) = -0.424800;
	rsigma_dot(  3, 18,  2,  2 ) = -0.424800;
	rsigma_dot(  4,  1,  2,  2 ) = -0.367700;
	rsigma_dot(  4,  2,  2,  2 ) = -0.367700;
	rsigma_dot(  4,  3,  2,  2 ) = -0.367700;
	rsigma_dot(  4,  4,  2,  2 ) = -0.367700;
	rsigma_dot(  4,  5,  2,  2 ) = -3.312100;
	rsigma_dot(  4,  6,  2,  2 ) = -3.412200;
	rsigma_dot(  4,  7,  2,  2 ) = -5.509340;
	rsigma_dot(  4,  8,  2,  2 ) = -4.128880;
	rsigma_dot(  4,  9,  2,  2 ) = -3.075700;
	rsigma_dot(  4, 10,  2,  2 ) = -2.765600;
	rsigma_dot(  4, 11,  2,  2 ) = -2.765600;
	rsigma_dot(  4, 12,  2,  2 ) = -2.765600;
	rsigma_dot(  4, 13,  2,  2 ) = -0.367700;
	rsigma_dot(  4, 14,  2,  2 ) = -0.367700;
	rsigma_dot(  4, 15,  2,  2 ) = -0.367700;
	rsigma_dot(  4, 16,  2,  2 ) = -0.367700;
	rsigma_dot(  4, 17,  2,  2 ) = -0.367700;
	rsigma_dot(  4, 18,  2,  2 ) = -0.367700;
	rsigma_dot(  5,  1,  2,  2 ) = -0.313600;
	rsigma_dot(  5,  2,  2,  2 ) = -0.313600;
	rsigma_dot(  5,  3,  2,  2 ) = -0.313600;
	rsigma_dot(  5,  4,  2,  2 ) = -1.412200;
	rsigma_dot(  5,  5,  2,  2 ) = -3.021700;
	rsigma_dot(  5,  6,  2,  2 ) = -4.779520;
	rsigma_dot(  5,  7,  2,  2 ) = -6.825360;
	rsigma_dot(  5,  8,  2,  2 ) = -4.657420;
	rsigma_dot(  5,  9,  2,  2 ) = -3.358100;
	rsigma_dot(  5, 10,  2,  2 ) = -3.609400;
	rsigma_dot(  5, 11,  2,  2 ) = -4.245440;
	rsigma_dot(  5, 12,  2,  2 ) = -3.021700;
	rsigma_dot(  5, 13,  2,  2 ) = -1.412200;
	rsigma_dot(  5, 14,  2,  2 ) = -0.313600;
	rsigma_dot(  5, 15,  2,  2 ) = -0.313600;
	rsigma_dot(  5, 16,  2,  2 ) = -0.313600;
	rsigma_dot(  5, 17,  2,  2 ) = -0.313600;
	rsigma_dot(  5, 18,  2,  2 ) = -0.313600;
	rsigma_dot(  6,  1,  2,  2 ) = -0.262300;
	rsigma_dot(  6,  2,  2,  2 ) = -0.262300;
	rsigma_dot(  6,  3,  2,  2 ) = -0.262300;
	rsigma_dot(  6,  4,  2,  2 ) = -1.871800;
	rsigma_dot(  6,  5,  2,  2 ) = -2.208200;
	rsigma_dot(  6,  6,  2,  2 ) = -6.862190;
	rsigma_dot(  6,  7,  2,  2 ) = -7.463490;
	rsigma_dot(  6,  8,  2,  2 ) = -4.837030;
	rsigma_dot(  6,  9,  2,  2 ) = -4.436700;
	rsigma_dot(  6, 10,  2,  2 ) = -4.916280;
	rsigma_dot(  6, 11,  2,  2 ) = -5.687270;
	rsigma_dot(  6, 12,  2,  2 ) = -4.405450;
	rsigma_dot(  6, 13,  2,  2 ) = -2.208200;
	rsigma_dot(  6, 14,  2,  2 ) = -0.262300;
	rsigma_dot(  6, 15,  2,  2 ) = -1.360900;
	rsigma_dot(  6, 16,  2,  2 ) = -0.262300;
	rsigma_dot(  6, 17,  2,  2 ) = -0.262300;
	rsigma_dot(  6, 18,  2,  2 ) = -0.262300;
	rsigma_dot(  7,  1,  2,  2 ) = -0.213500;
	rsigma_dot(  7,  2,  2,  2 ) = -0.213500;
	rsigma_dot(  7,  3,  2,  2 ) = -0.213500;
	rsigma_dot(  7,  4,  2,  2 ) = -0.213500;
	rsigma_dot(  7,  5,  2,  2 ) = -4.256580;
	rsigma_dot(  7,  6,  2,  2 ) = -8.696750;
	rsigma_dot(  7,  7,  2,  2 ) = -7.275720;
	rsigma_dot(  7,  8,  2,  2 ) = -4.746130;
	rsigma_dot(  7,  9,  2,  2 ) = -4.145350;
	rsigma_dot(  7, 10,  2,  2 ) = -4.923060;
	rsigma_dot(  7, 11,  2,  2 ) = -6.725270;
	rsigma_dot(  7, 12,  2,  2 ) = -5.620700;
	rsigma_dot(  7, 13,  2,  2 ) = -2.611400;
	rsigma_dot(  7, 14,  2,  2 ) = -0.213500;
	rsigma_dot(  7, 15,  2,  2 ) = -0.213500;
	rsigma_dot(  7, 16,  2,  2 ) = -0.213500;
	rsigma_dot(  7, 17,  2,  2 ) = -0.213500;
	rsigma_dot(  7, 18,  2,  2 ) = -0.213500;
	rsigma_dot(  8,  1,  2,  2 ) = -0.167000;
	rsigma_dot(  8,  2,  2,  2 ) = -0.167000;
	rsigma_dot(  8,  3,  2,  2 ) = -0.167000;
	rsigma_dot(  8,  4,  2,  2 ) = -1.265600;
	rsigma_dot(  8,  5,  2,  2 ) = -6.383610;
	rsigma_dot(  8,  6,  2,  2 ) = -9.098560;
	rsigma_dot(  8,  7,  2,  2 ) = -6.418910;
	rsigma_dot(  8,  8,  2,  2 ) = -4.561460;
	rsigma_dot(  8,  9,  2,  2 ) = -4.174340;
	rsigma_dot(  8, 10,  2,  2 ) = -4.979190;
	rsigma_dot(  8, 11,  2,  2 ) = -7.462060;
	rsigma_dot(  8, 12,  2,  2 ) = -7.337890;
	rsigma_dot(  8, 13,  2,  2 ) = -3.722400;
	rsigma_dot(  8, 14,  2,  2 ) = -0.167000;
	rsigma_dot(  8, 15,  2,  2 ) = -0.167000;
	rsigma_dot(  8, 16,  2,  2 ) = -0.167000;
	rsigma_dot(  8, 17,  2,  2 ) = -0.167000;
	rsigma_dot(  8, 18,  2,  2 ) = -0.167000;
	rsigma_dot(  9,  1,  2,  2 ) = -0.122600;
	rsigma_dot(  9,  2,  2,  2 ) = -0.122600;
	rsigma_dot(  9,  3,  2,  2 ) = -0.122600;
	rsigma_dot(  9,  4,  2,  2 ) = -1.221200;
	rsigma_dot(  9,  5,  2,  2 ) = -7.827370;
	rsigma_dot(  9,  6,  2,  2 ) = -8.067400;
	rsigma_dot(  9,  7,  2,  2 ) = -5.809530;
	rsigma_dot(  9,  8,  2,  2 ) = -4.966740;
	rsigma_dot(  9,  9,  2,  2 ) = -4.327250;
	rsigma_dot(  9, 10,  2,  2 ) = -4.884730;
	rsigma_dot(  9, 11,  2,  2 ) = -7.072410;
	rsigma_dot(  9, 12,  2,  2 ) = -8.543900;
	rsigma_dot(  9, 13,  2,  2 ) = -5.252450;
	rsigma_dot(  9, 14,  2,  2 ) = -1.732000;
	rsigma_dot(  9, 15,  2,  2 ) = -0.122600;
	rsigma_dot(  9, 16,  2,  2 ) = -0.122600;
	rsigma_dot(  9, 17,  2,  2 ) = -0.122600;
	rsigma_dot(  9, 18,  2,  2 ) = -0.122600;
	rsigma_dot( 10,  1,  2,  2 ) = -0.080000;
	rsigma_dot( 10,  2,  2,  2 ) = -0.080000;
	rsigma_dot( 10,  3,  2,  2 ) = -1.689400;
	rsigma_dot( 10,  4,  2,  2 ) = -2.913200;
	rsigma_dot( 10,  5,  2,  2 ) = -7.579420;
	rsigma_dot( 10,  6,  2,  2 ) = -6.907620;
	rsigma_dot( 10,  7,  2,  2 ) = -5.422330;
	rsigma_dot( 10,  8,  2,  2 ) = -5.161400;
	rsigma_dot( 10,  9,  2,  2 ) = -4.590850;
	rsigma_dot( 10, 10,  2,  2 ) = -4.970340;
	rsigma_dot( 10, 11,  2,  2 ) = -6.331900;
	rsigma_dot( 10, 12,  2,  2 ) = -8.504850;
	rsigma_dot( 10, 13,  2,  2 ) = -6.881280;
	rsigma_dot( 10, 14,  2,  2 ) = -2.644900;
	rsigma_dot( 10, 15,  2,  2 ) = -0.080000;
	rsigma_dot( 10, 16,  2,  2 ) = -1.178600;
	rsigma_dot( 10, 17,  2,  2 ) = -0.080000;
	rsigma_dot( 10, 18,  2,  2 ) = -0.080000;
	rsigma_dot( 11,  1,  2,  2 ) = -0.039200;
	rsigma_dot( 11,  2,  2,  2 ) = -1.137800;
	rsigma_dot( 11,  3,  2,  2 ) = -2.437100;
	rsigma_dot( 11,  4,  2,  2 ) = -4.356660;
	rsigma_dot( 11,  5,  2,  2 ) = -6.159470;
	rsigma_dot( 11,  6,  2,  2 ) = -6.109910;
	rsigma_dot( 11,  7,  2,  2 ) = -5.082600;
	rsigma_dot( 11,  8,  2,  2 ) = -5.203960;
	rsigma_dot( 11,  9,  2,  2 ) = -5.029610;
	rsigma_dot( 11, 10,  2,  2 ) = -5.180840;
	rsigma_dot( 11, 11,  2,  2 ) = -5.633880;
	rsigma_dot( 11, 12,  2,  2 ) = -7.850740;
	rsigma_dot( 11, 13,  2,  2 ) = -8.236990;
	rsigma_dot( 11, 14,  2,  2 ) = -3.258000;
	rsigma_dot( 11, 15,  2,  2 ) = -2.604100;
	rsigma_dot( 11, 16,  2,  2 ) = -1.648600;
	rsigma_dot( 11, 17,  2,  2 ) = -0.039200;
	rsigma_dot( 11, 18,  2,  2 ) = -0.039200;
	rsigma_dot( 12,  1,  2,  2 ) = 0.000000;
	rsigma_dot( 12,  2,  2,  2 ) = -2.197200;
	rsigma_dot( 12,  3,  2,  2 ) = -2.708000;
	rsigma_dot( 12,  4,  2,  2 ) = -5.093700;
	rsigma_dot( 12,  5,  2,  2 ) = -5.802070;
	rsigma_dot( 12,  6,  2,  2 ) = -5.587200;
	rsigma_dot( 12,  7,  2,  2 ) = -5.198450;
	rsigma_dot( 12,  8,  2,  2 ) = -5.459540;
	rsigma_dot( 12,  9,  2,  2 ) = -5.068860;
	rsigma_dot( 12, 10,  2,  2 ) = -5.272950;
	rsigma_dot( 12, 11,  2,  2 ) = -5.370590;
	rsigma_dot( 12, 12,  2,  2 ) = -6.597100;
	rsigma_dot( 12, 13,  2,  2 ) = -7.631870;
	rsigma_dot( 12, 14,  2,  2 ) = -5.043380;
	rsigma_dot( 12, 15,  2,  2 ) = -3.218800;
	rsigma_dot( 12, 16,  2,  2 ) = -2.833200;
	rsigma_dot( 12, 17,  2,  2 ) = -2.197200;
	rsigma_dot( 12, 18,  2,  2 ) = 0.000000;

	// section 21
	rsigma_dot(  1,  1,  2,  1 ) = -0.974100;
	rsigma_dot(  1,  2,  2,  1 ) = -0.974100;
	rsigma_dot(  1,  3,  2,  1 ) = -0.974100;
	rsigma_dot(  1,  4,  2,  1 ) = -0.974100;
	rsigma_dot(  1,  5,  2,  1 ) = -0.974100;
	rsigma_dot(  1,  6,  2,  1 ) = -0.974100;
	rsigma_dot(  1,  7,  2,  1 ) = -0.974100;
	rsigma_dot(  1,  8,  2,  1 ) = -2.920000;
	rsigma_dot(  1,  9,  2,  1 ) = -3.171400;
	rsigma_dot(  1, 10,  2,  1 ) = -0.974100;
	rsigma_dot(  1, 11,  2,  1 ) = -2.583600;
	rsigma_dot(  1, 12,  2,  1 ) = -0.974100;
	rsigma_dot(  1, 13,  2,  1 ) = -0.974100;
	rsigma_dot(  1, 14,  2,  1 ) = -0.974100;
	rsigma_dot(  1, 15,  2,  1 ) = -0.974100;
	rsigma_dot(  1, 16,  2,  1 ) = -0.974100;
	rsigma_dot(  1, 17,  2,  1 ) = -0.974100;
	rsigma_dot(  1, 18,  2,  1 ) = -0.974100;
	rsigma_dot(  2,  1,  2,  1 ) = -0.909600;
	rsigma_dot(  2,  2,  2,  1 ) = -0.909600;
	rsigma_dot(  2,  3,  2,  1 ) = -0.909600;
	rsigma_dot(  2,  4,  2,  1 ) = -0.909600;
	rsigma_dot(  2,  5,  2,  1 ) = -0.909600;
	rsigma_dot(  2,  6,  2,  1 ) = -0.909600;
	rsigma_dot(  2,  7,  2,  1 ) = -0.909600;
	rsigma_dot(  2,  8,  2,  1 ) = -4.343580;
	rsigma_dot(  2,  9,  2,  1 ) = -6.015540;
	rsigma_dot(  2, 10,  2,  1 ) = -4.205430;
	rsigma_dot(  2, 11,  2,  1 ) = -2.008200;
	rsigma_dot(  2, 12,  2,  1 ) = -2.519000;
	rsigma_dot(  2, 13,  2,  1 ) = -0.909600;
	rsigma_dot(  2, 14,  2,  1 ) = -0.909600;
	rsigma_dot(  2, 15,  2,  1 ) = -0.909600;
	rsigma_dot(  2, 16,  2,  1 ) = -0.909600;
	rsigma_dot(  2, 17,  2,  1 ) = -0.909600;
	rsigma_dot(  2, 18,  2,  1 ) = -0.909600;
	rsigma_dot(  3,  1,  2,  1 ) = -0.849000;
	rsigma_dot(  3,  2,  2,  1 ) = -0.849000;
	rsigma_dot(  3,  3,  2,  1 ) = -0.849000;
	rsigma_dot(  3,  4,  2,  1 ) = -0.849000;
	rsigma_dot(  3,  5,  2,  1 ) = -0.849000;
	rsigma_dot(  3,  6,  2,  1 ) = -0.849000;
	rsigma_dot(  3,  7,  2,  1 ) = -3.046200;
	rsigma_dot(  3,  8,  2,  1 ) = -5.990630;
	rsigma_dot(  3,  9,  2,  1 ) = -7.627760;
	rsigma_dot(  3, 10,  2,  1 ) = -6.013760;
	rsigma_dot(  3, 11,  2,  1 ) = -2.794900;
	rsigma_dot(  3, 12,  2,  1 ) = -2.458400;
	rsigma_dot(  3, 13,  2,  1 ) = -0.849000;
	rsigma_dot(  3, 14,  2,  1 ) = -0.849000;
	rsigma_dot(  3, 15,  2,  1 ) = -0.849000;
	rsigma_dot(  3, 16,  2,  1 ) = -0.849000;
	rsigma_dot(  3, 17,  2,  1 ) = -0.849000;
	rsigma_dot(  3, 18,  2,  1 ) = -0.849000;
	rsigma_dot(  4,  1,  2,  1 ) = -0.791800;
	rsigma_dot(  4,  2,  2,  1 ) = -0.791800;
	rsigma_dot(  4,  3,  2,  1 ) = -0.791800;
	rsigma_dot(  4,  4,  2,  1 ) = -0.791800;
	rsigma_dot(  4,  5,  2,  1 ) = -0.791800;
	rsigma_dot(  4,  6,  2,  1 ) = -2.401300;
	rsigma_dot(  4,  7,  2,  1 ) = -3.736250;
	rsigma_dot(  4,  8,  2,  1 ) = -6.512120;
	rsigma_dot(  4,  9,  2,  1 ) = -8.856760;
	rsigma_dot(  4, 10,  2,  1 ) = -7.235940;
	rsigma_dot(  4, 11,  2,  1 ) = -4.010690;
	rsigma_dot(  4, 12,  2,  1 ) = -2.401300;
	rsigma_dot(  4, 13,  2,  1 ) = -1.890400;
	rsigma_dot(  4, 14,  2,  1 ) = -0.791800;
	rsigma_dot(  4, 15,  2,  1 ) = -0.791800;
	rsigma_dot(  4, 16,  2,  1 ) = -0.791800;
	rsigma_dot(  4, 17,  2,  1 ) = -0.791800;
	rsigma_dot(  4, 18,  2,  1 ) = -0.791800;
	rsigma_dot(  5,  1,  2,  1 ) = -0.737700;
	rsigma_dot(  5,  2,  2,  1 ) = -0.737700;
	rsigma_dot(  5,  3,  2,  1 ) = -0.737700;
	rsigma_dot(  5,  4,  2,  1 ) = -0.737700;
	rsigma_dot(  5,  5,  2,  1 ) = -0.737700;
	rsigma_dot(  5,  6,  2,  1 ) = -2.683700;
	rsigma_dot(  5,  7,  2,  1 ) = -4.171730;
	rsigma_dot(  5,  8,  2,  1 ) = -7.436010;
	rsigma_dot(  5,  9,  2,  1 ) = -9.773380;
	rsigma_dot(  5, 10,  2,  1 ) = -8.510500;
	rsigma_dot(  5, 11,  2,  1 ) = -4.669570;
	rsigma_dot(  5, 12,  2,  1 ) = -3.445800;
	rsigma_dot(  5, 13,  2,  1 ) = -1.836400;
	rsigma_dot(  5, 14,  2,  1 ) = -0.737700;
	rsigma_dot(  5, 15,  2,  1 ) = -0.737700;
	rsigma_dot(  5, 16,  2,  1 ) = -0.737700;
	rsigma_dot(  5, 17,  2,  1 ) = -0.737700;
	rsigma_dot(  5, 18,  2,  1 ) = -0.737700;
	rsigma_dot(  6,  1,  2,  1 ) = -0.686500;
	rsigma_dot(  6,  2,  2,  1 ) = -0.686500;
	rsigma_dot(  6,  3,  2,  1 ) = -0.686500;
	rsigma_dot(  6,  4,  2,  1 ) = -0.686500;
	rsigma_dot(  6,  5,  2,  1 ) = -0.686500;
	rsigma_dot(  6,  6,  2,  1 ) = -0.686500;
	rsigma_dot(  6,  7,  2,  1 ) = -4.297370;
	rsigma_dot(  6,  8,  2,  1 ) = -7.670240;
	rsigma_dot(  6,  9,  2,  1 ) = -9.274480;
	rsigma_dot(  6, 10,  2,  1 ) = -8.780530;
	rsigma_dot(  6, 11,  2,  1 ) = -5.395980;
	rsigma_dot(  6, 12,  2,  1 ) = -2.632400;
	rsigma_dot(  6, 13,  2,  1 ) = -3.084300;
	rsigma_dot(  6, 14,  2,  1 ) = -0.686500;
	rsigma_dot(  6, 15,  2,  1 ) = -0.686500;
	rsigma_dot(  6, 16,  2,  1 ) = -0.686500;
	rsigma_dot(  6, 17,  2,  1 ) = -0.686500;
	rsigma_dot(  6, 18,  2,  1 ) = -0.686500;
	rsigma_dot(  7,  1,  2,  1 ) = -0.637700;
	rsigma_dot(  7,  2,  2,  1 ) = -0.637700;
	rsigma_dot(  7,  3,  2,  1 ) = -0.637700;
	rsigma_dot(  7,  4,  2,  1 ) = -0.637700;
	rsigma_dot(  7,  5,  2,  1 ) = -1.736300;
	rsigma_dot(  7,  6,  2,  1 ) = -3.682200;
	rsigma_dot(  7,  7,  2,  1 ) = -4.645000;
	rsigma_dot(  7,  8,  2,  1 ) = -6.967380;
	rsigma_dot(  7,  9,  2,  1 ) = -7.766160;
	rsigma_dot(  7, 10,  2,  1 ) = -7.608390;
	rsigma_dot(  7, 11,  2,  1 ) = -5.416790;
	rsigma_dot(  7, 12,  2,  1 ) = -3.773160;
	rsigma_dot(  7, 13,  2,  1 ) = -3.345700;
	rsigma_dot(  7, 14,  2,  1 ) = -2.834900;
	rsigma_dot(  7, 15,  2,  1 ) = -2.247100;
	rsigma_dot(  7, 16,  2,  1 ) = -0.637700;
	rsigma_dot(  7, 17,  2,  1 ) = -0.637700;
	rsigma_dot(  7, 18,  2,  1 ) = -0.637700;
	rsigma_dot(  8,  1,  2,  1 ) = -0.591100;
	rsigma_dot(  8,  2,  2,  1 ) = -0.591100;
	rsigma_dot(  8,  3,  2,  1 ) = -0.591100;
	rsigma_dot(  8,  4,  2,  1 ) = -0.591100;
	rsigma_dot(  8,  5,  2,  1 ) = -2.788400;
	rsigma_dot(  8,  6,  2,  1 ) = -4.087650;
	rsigma_dot(  8,  7,  2,  1 ) = -4.702020;
	rsigma_dot(  8,  8,  2,  1 ) = -6.084200;
	rsigma_dot(  8,  9,  2,  1 ) = -6.877140;
	rsigma_dot(  8, 10,  2,  1 ) = -6.827510;
	rsigma_dot(  8, 11,  2,  1 ) = -5.980210;
	rsigma_dot(  8, 12,  2,  1 ) = -4.087650;
	rsigma_dot(  8, 13,  2,  1 ) = -3.635700;
	rsigma_dot(  8, 14,  2,  1 ) = -4.202060;
	rsigma_dot(  8, 15,  2,  1 ) = -2.537100;
	rsigma_dot(  8, 16,  2,  1 ) = -0.591100;
	rsigma_dot(  8, 17,  2,  1 ) = -0.591100;
	rsigma_dot(  8, 18,  2,  1 ) = -0.591100;
	rsigma_dot(  9,  1,  2,  1 ) = -0.546700;
	rsigma_dot(  9,  2,  2,  1 ) = -0.546700;
	rsigma_dot(  9,  3,  2,  1 ) = -0.546700;
	rsigma_dot(  9,  4,  2,  1 ) = -0.546700;
	rsigma_dot(  9,  5,  2,  1 ) = -3.254700;
	rsigma_dot(  9,  6,  2,  1 ) = -4.102040;
	rsigma_dot(  9,  7,  2,  1 ) = -4.916140;
	rsigma_dot(  9,  8,  2,  1 ) = -5.788440;
	rsigma_dot(  9,  9,  2,  1 ) = -5.777800;
	rsigma_dot(  9, 10,  2,  1 ) = -6.126420;
	rsigma_dot(  9, 11,  2,  1 ) = -5.997730;
	rsigma_dot(  9, 12,  2,  1 ) = -4.721080;
	rsigma_dot(  9, 13,  2,  1 ) = -4.102040;
	rsigma_dot(  9, 14,  2,  1 ) = -4.157610;
	rsigma_dot(  9, 15,  2,  1 ) = -3.111600;
	rsigma_dot(  9, 16,  2,  1 ) = -1.645300;
	rsigma_dot(  9, 17,  2,  1 ) = -0.546700;
	rsigma_dot(  9, 18,  2,  1 ) = -0.546700;
	rsigma_dot( 10,  1,  2,  1 ) = -0.504100;
	rsigma_dot( 10,  2,  2,  1 ) = -0.504100;
	rsigma_dot( 10,  3,  2,  1 ) = -0.504100;
	rsigma_dot( 10,  4,  2,  1 ) = -0.504100;
	rsigma_dot( 10,  5,  2,  1 ) = -3.548700;
	rsigma_dot( 10,  6,  2,  1 ) = -4.873580;
	rsigma_dot( 10,  7,  2,  1 ) = -5.176960;
	rsigma_dot( 10,  8,  2,  1 ) = -5.756400;
	rsigma_dot( 10,  9,  2,  1 ) = -5.573030;
	rsigma_dot( 10, 10,  2,  1 ) = -5.702630;
	rsigma_dot( 10, 11,  2,  1 ) = -5.668920;
	rsigma_dot( 10, 12,  2,  1 ) = -5.078840;
	rsigma_dot( 10, 13,  2,  1 ) = -4.678520;
	rsigma_dot( 10, 14,  2,  1 ) = -4.310790;
	rsigma_dot( 10, 15,  2,  1 ) = -4.217700;
	rsigma_dot( 10, 16,  2,  1 ) = -2.113600;
	rsigma_dot( 10, 17,  2,  1 ) = -0.504100;
	rsigma_dot( 10, 18,  2,  1 ) = -0.504100;
	rsigma_dot( 11,  1,  2,  1 ) = -0.463300;
	rsigma_dot( 11,  2,  2,  1 ) = -0.463300;
	rsigma_dot( 11,  3,  2,  1 ) = -0.463300;
	rsigma_dot( 11,  4,  2,  1 ) = -2.409200;
	rsigma_dot( 11,  5,  2,  1 ) = -4.313460;
	rsigma_dot( 11,  6,  2,  1 ) = -5.225480;
	rsigma_dot( 11,  7,  2,  1 ) = -5.493750;
	rsigma_dot( 11,  8,  2,  1 ) = -5.353660;
	rsigma_dot( 11,  9,  2,  1 ) = -5.154660;
	rsigma_dot( 11, 10,  2,  1 ) = -5.397780;
	rsigma_dot( 11, 11,  2,  1 ) = -5.726000;
	rsigma_dot( 11, 12,  2,  1 ) = -5.480590;
	rsigma_dot( 11, 13,  2,  1 ) = -4.018660;
	rsigma_dot( 11, 14,  2,  1 ) = -4.313460;
	rsigma_dot( 11, 15,  2,  1 ) = -4.574180;
	rsigma_dot( 11, 16,  2,  1 ) = -2.861200;
	rsigma_dot( 11, 17,  2,  1 ) = -0.463300;
	rsigma_dot( 11, 18,  2,  1 ) = -1.561900;
	rsigma_dot( 12,  1,  2,  1 ) = -0.424100;
	rsigma_dot( 12,  2,  2,  1 ) = -0.424100;
	rsigma_dot( 12,  3,  2,  1 ) = -2.822000;
	rsigma_dot( 12,  4,  2,  1 ) = -3.132100;
	rsigma_dot( 12,  5,  2,  1 ) = -5.115440;
	rsigma_dot( 12,  6,  2,  1 ) = -5.697090;
	rsigma_dot( 12,  7,  2,  1 ) = -5.747100;
	rsigma_dot( 12,  8,  2,  1 ) = -5.386930;
	rsigma_dot( 12,  9,  2,  1 ) = -5.019210;
	rsigma_dot( 12, 10,  2,  1 ) = -5.186260;
	rsigma_dot( 12, 11,  2,  1 ) = -5.775950;
	rsigma_dot( 12, 12,  2,  1 ) = -5.600240;
	rsigma_dot( 12, 13,  2,  1 ) = -4.394380;
	rsigma_dot( 12, 14,  2,  1 ) = -3.257300;
	rsigma_dot( 12, 15,  2,  1 ) = -3.920600;
	rsigma_dot( 12, 16,  2,  1 ) = -3.643000;
	rsigma_dot( 12, 17,  2,  1 ) = -1.522700;
	rsigma_dot( 12, 18,  2,  1 ) = -0.424100;
}

void
SecondaryStructurePotential::m_term_initializer(
	FArray1D_real & m_term
)
{
	m_term(1) = 1.87;
	m_term(2) = .61;
	m_term(3) = .74;
	m_term(4) = .17; // score for 4 or more
}

/// @brief Penalty for pairing strand dimers that are close in sequence.
/// @details Inferred from the log ratio of pairing probabilities of strands
/// @details in the PDB vs. in Rosetta decoys. Calculated as a function of
/// @details strand separation.
void
SecondaryStructurePotential::ss_penalty_initializer(
	FArray1D_real & SS_penalty
)
{
	// For strand separations less than 5, statistics become small, so
	// set penalty to be a constant.
	SS_penalty( 1) = 1.13386;
	SS_penalty( 2) = 1.13386;
	SS_penalty( 3) = 1.13386;
	SS_penalty( 4) = 1.13386;
	SS_penalty( 5) = 1.13386;
	SS_penalty( 6) = 0.70241;
	SS_penalty( 7) = 0.57908;
	SS_penalty( 8) = 0.44451;
	SS_penalty( 9) = 0.31653;
	SS_penalty(10) = 0.22074;
	SS_penalty(11) = 0.14869;
}


} // ns scoring
} // ns core



/// @brief This function takes a set of dimer neighbors, and determines how
/// @brief many sheets there, and how many strands are in each sheet
/// @brief This information is then used to calculate the "poker hand" score,
/// @brief which reflects to probability of that distribution of strands and
/// @brief sheets.
//
//js This function takes a set of dimer neighbors, and determines how
//js many sheets there, and how many strands are in each sheet
//js This information is then used to calculate the "poker hand" score,
//js which reflects to probability of that distribution of strands and
//js sheets.
//js In current version, it seems to simply penalize sheets with fewer
//js strands compared to those with more strands.
//
//js This function looks at a list of dimers, which contains up to
//js two neighbors for each dimer.  In priniciple these neighbors would
//js be hydrogen bond partners in neighboring strands.  This function
//js will take this list, however it is made.
//
//js Currently, dimer neighbors are defined somewhat arbitrarily.
//js If I understand the code correctly, the first and last dimers in
//js sequence that follow the current dimer, and that are within 6.5
//js angstroms, are the neighbors.  There is no orientation dependence
//js on what counts as a strand neighbor.
//
//js A sheet is then loosely defined by all the strands that are connected
//js by one of these neighbors.  This amounts to "single-linkage clustering."
//js A sheet is determined by the set of strands that can be linked by
//js the dimer neighbors.  Note that a one neighbor is enough to say that
//js two strands are in a sheet.
//
//js The final score is put into sheet_score, and is determined by the
//js number of sheets of each size.
//
//js   Basic strategy:  go through all dimers, finding neigboring strands
//js   Because each dimer can only have two neighbors we can first search
//js   all the way down one side, and then the other.  This will give a set
//js   of dimers that are connected.  The strands in which these dimers reside
//js   are then considered connected into sheets.  By going through all
//js   dimers with the proper bookkeeping, we can determine how many sheets
//js   there are, and how many strands in e


/**

NOT PORTING IN FIRST PASS (PB)

void
SecondaryStructurePotential::sheets_from_dimers(
	Real & sheet_score
)
{
	int const & total_residue = *total_residue_; // yab: misc removal

	static FArray1D_bool searched( MAX_RES() );
	static FArray1D_int strand_sheet( MAX_RES() );
	 // 40 is the maximum number of strands
	static FArray1D_int num_of_strands( MAX_RES() );
	 // 11 is the maxmumber number of sheets
	static FArray2D_int strand_sheet_list( MAX_RES(), MAX_RES() );
	static FArray1D_real const m_term( 4, m_term_initializer );

	for ( int current_dimer = 1; current_dimer <= strands.total_SS_dimer;
	 ++current_dimer ) {
//js      Set all dimers as unchecked.
		searched(current_dimer) = false;
//js      set all sheet locations as null
//js         dimer_sheet(current_dimer) = 0
	}
	for ( int current_strand = 1; current_strand <= strands.total_strands;
	 ++current_strand ) {
//js      Set the sheets of all strands as null
		strand_sheet(current_strand) = 0;
	}
	for ( int current_sheet = 1; current_sheet <= total_residue;
	 ++current_sheet ) {
		num_of_strands(current_sheet) = 0;
	}

//js Find the neighbors of each dimer.  Some will be found during
//js the search, some will be initial nodes.  That is why we keep
//js track of whether it has been searched.
	int num_of_sheets = 0;
	int current_sheet = 0;
	for ( int current_dimer = 1; current_dimer <= strands.total_SS_dimer;
	 ++current_dimer ) {
		if ( !searched(current_dimer) ) {
//js we need to check this one
			searched(current_dimer) = true;
			int current_strand = strands.SS_strand(current_dimer);
			// place node strand in sheet
			if ( strand_sheet(current_strand) == 0 ) { // it is not in a sheet, so:
				// make new sheet
				++num_of_sheets;
				current_sheet = num_of_sheets;

				// place strand in current sheet
				strand_sheet(current_strand) = current_sheet;
				++num_of_strands(current_sheet);
				strand_sheet_list(num_of_strands(current_sheet),current_sheet) =
				 current_strand;
			} else {
				current_sheet = strand_sheet(current_strand);
			}
			for ( int direction = 1; direction <= 2; ++direction ) {
			 // the two directions of searching
				int neighbor = strands.dimer_neighbor(direction,current_dimer);
				while ( neighbor != 0 ) {
//               if ( neighbor != 0 ) {
//js                  if ( !searched(neighbor) ) {
					searched(neighbor) = true;

					current_strand = strands.SS_strand(neighbor);
					if ( strand_sheet(current_strand) == 0 ) {
// js  if neighbor strand is not in a sheet already, put it in the working sheet
						strand_sheet(current_strand) = current_sheet;
						++num_of_strands(current_sheet);

						strand_sheet_list(num_of_strands(current_sheet),current_sheet) =
						 current_strand;
					} else if ( strand_sheet(current_strand) != current_sheet ) {
// js if neighbor strand is  already in a different sheet, merge the sheets.
// js the sheet of the new strand must have a lower sheet number, so give the
// js strands of the current sheet to its sheet
						int const new_sheet = strand_sheet(current_strand);
						int & num_of_strandsnew_sheet( num_of_strands(new_sheet) );
						for ( int merge = 1,
						 lss = strand_sheet_list.index(merge,current_sheet),
						 mergee = num_of_strands(current_sheet);
						 merge <= mergee; ++merge, ++lss ) {
							++num_of_strandsnew_sheet;
							int const strand_sheet_list_mc = strand_sheet_list[ lss ];
							 // strand_sheet_list(merge,current_sheet)
							strand_sheet_list(num_of_strandsnew_sheet,new_sheet) =
							 strand_sheet_list_mc;
							strand_sheet(strand_sheet_list_mc) = new_sheet;
						}
						--num_of_sheets;
// rhiju After merging one sheet with another, need to erase traces
// rhiju of sheet that got eaten up, and reorder other sheets.
						num_of_strands(current_sheet) = 0;
						for (int shiftsheet = current_sheet; shiftsheet <= num_of_sheets;
								 ++shiftsheet){
							num_of_strands(shiftsheet) = num_of_strands(shiftsheet+1);
							for (int i = 1; i <= num_of_strands(shiftsheet); ++i){
								int strandtoshift =	strand_sheet_list(i,shiftsheet+1);
								strand_sheet_list(i,shiftsheet) = strandtoshift;
								strand_sheet(strandtoshift) = shiftsheet;
							}
						}
						current_sheet = new_sheet;
					}
//js                  }
					neighbor = strands.dimer_neighbor(direction,neighbor);
				}
			}
		}
	}

//js calculate score, based on the number of sheets of each size
	Real sheet_score_sum = 0.0;
	for ( int current_sheet = 1; current_sheet <= num_of_sheets; ++current_sheet ) {
		sheet_score_sum += m_term( std::min( num_of_strands(current_sheet), 4 ) );
	}
	sheet_score = sheet_score_sum;
	//	sheet_score *+ get_sheet_wt();

	// FIXME: keep Ingo's sheet filter?  the sheet filter class is in classic rosetta
	//        sheet_filter.h/sheet_filter.cc and should be self contained, so it's
	//        directly liftable into mini without much modification
	// FIXME: need equivalent for files::paths below
	if ( get_handedness_score_flag() && files_paths::use_filter(files_paths::sheet_type) ){
		int result = 0;
		sheet_filter::SheetFilter sf(position_, secstruct_, total_residue_);
		sf.compute_result(result); // Trigger evaluation if Ingo's sheet filter, including handedness checks.
		sheet_score += sf.get_handedness_score();
	}

	// modify by proper weighting
	sheet_score *= 2.019 * get_sheet_weight();
}

**/

