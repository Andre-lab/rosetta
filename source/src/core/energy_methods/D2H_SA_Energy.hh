// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/D2H_SA_Energy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Björn Wallner


#ifndef INCLUDED_core_energy_methods_D2H_SA_Energy_hh
#define INCLUDED_core_energy_methods_D2H_SA_Energy_hh


// Package headers
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

// Utility headers


namespace core {
namespace energy_methods {



class D2H_SA_Energy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy  parent;
public:

	/// @brief call sasa.cc to calculate the surface area
	D2H_SA_Energy();

	/// clone
	core::scoring::methods::EnergyMethodOP
	clone() const override;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	void
	finalize_total_energy(
		pose::Pose & pose,
		core::scoring::ScoreFunction const &,
		core::scoring::EnergyMap & totals
	) const override;
	core::Size version() const override;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const override {}

	//void
	//D2H_SA_Energy::set_rsa_range(Size s,Size e) {};

private:
	//utility::vector1< core::Real > all_data_;
	//utility::vector1< core::Real > all_group_;
	utility::vector1< core::Real > data_;
	utility::vector1< core::Real > position_;
	utility::vector1< core::Size > group_;


	// utility::vector1< core::Real > data2_;
	//utility::vector1< core::Real > position2_;
	//utility::vector1< core::Size > group2_;
	// ObjexxFCL::FArray2D< core::Real > time_data_;
	utility::vector1< core::Real >  time_data_mean_;
	utility::vector1< core::Real >  time_data_std_;
	Real mean_;
	Real sd_;
	Size chain_; // chain to use for surface calculation, usually the one in the middle...
	Size rsa_index_start_; //avoid searching for chain each scoring...
	Size rsa_index_end_;
	bool HDX_data_defined_;
	Real reweight_;
};


}
}

#endif // INCLUDED_core_energy_methods_D2H_SA_Energy_HH
