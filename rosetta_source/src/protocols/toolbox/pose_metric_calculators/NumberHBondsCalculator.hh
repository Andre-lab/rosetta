// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/protocols/toolbox/PoseMetricCalculators/NumberHBondsCalculator.hh
/// @brief
/// @author Florian Richter


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_NumberHBondsCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_NumberHBondsCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/id/AtomID_Map.hh>
//#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>

// AUTO-REMOVED #include <utility/vector1.hh>

#include <set>

#include <utility/vector1.hh>


namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {

class NumberHBondsCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:


  NumberHBondsCalculator() :
    hb_database( core::scoring::hbonds::HBondDatabase::get_database( "standard_params" ) )
  {};

  NumberHBondsCalculator( std::set< core::Size > special_region ) :
    hb_database( core::scoring::hbonds::HBondDatabase::get_database( "standard_params" ) ),
    special_region_(special_region)
  {};


  core::pose::metrics::PoseMetricCalculatorOP clone() const { return new NumberHBondsCalculator(); };

  static core::Real
  sum_Hbond_terms( core::scoring::EnergyMap const & emap );

protected:

  virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
  virtual std::string print( std::string const & key ) const;
  virtual void recompute( core::pose::Pose const & this_pose );

  void
  determine_res_to_recompute(
    core::pose::Pose const & pose,
    utility::vector1< bool > & res_to_recompute
  );

  void
  compute_Hbonds_for_residue(
    core::pose::Pose const & pose,
    core::Size i,
    utility::vector1< bool > const & res_to_recompute,
    core::scoring::hbonds::HBondSet & hb_set
  );

private:

  core::scoring::hbonds::HBondDatabaseCOP hb_database;
  core::Size all_Hbonds_;
  core::Size special_region_Hbonds_;

  core::id::AtomID_Map< core::Size > atom_Hbonds_;
  utility::vector1< core::Size > residue_Hbonds_;

  //holds the calculated energies to prevent unnecessary recalculation
  utility::vector1< core::Real > ref_residue_total_energies_;

  std::set< core::Size > special_region_;
};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#endif
