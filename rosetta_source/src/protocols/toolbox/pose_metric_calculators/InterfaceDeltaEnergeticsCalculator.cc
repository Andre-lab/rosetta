// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metrics/SasaCalculator.cc
/// @brief  SasaCalculator class
/// @author John Karanicolas

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/InterfaceDeltaEnergeticsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/methods/Methods.hh> //for long range energies
#include <core/scoring/LREnergyContainer.hh> //long range energies

// Utility headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <cassert>

using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;
using namespace utility;

namespace protocols{
namespace toolbox {
namespace pose_metric_calculators {


// preferred constructor - use an existing InterfaceNeighborDefinitionCalculator
InterfaceDeltaEnergeticsCalculator::InterfaceDeltaEnergeticsCalculator( std::string const & NameOfInterfaceNeighborDefinitionCalculator ) :
	EnergyDependentCalculator(),
	name_of_InterfaceNeighborDefinitionCalculator_(NameOfInterfaceNeighborDefinitionCalculator)
{
	if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( name_of_InterfaceNeighborDefinitionCalculator_ ) ) {
		basic::Error() << "Tried to tie InterfaceDeltaEnergeticsCalculator to InterfaceNeighborDefinitionCalculator " <<
			name_of_InterfaceNeighborDefinitionCalculator_ << " but this calculator does not exist." << std::endl;
		utility_exit();
	}
}


// less preferred constructor - creates a new InterfaceNeighborDefinitionCalculator
InterfaceDeltaEnergeticsCalculator::InterfaceDeltaEnergeticsCalculator( Size const chain1_number, Size const chain2_number ) :
	EnergyDependentCalculator()
{
	name_of_InterfaceNeighborDefinitionCalculator_ = "IEC_private_IDC_" + to_string(chain1_number) + "_" + to_string(chain2_number);
	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( name_of_InterfaceNeighborDefinitionCalculator_ ) ) {
		basic::Error() << "InterfaceDeltaEnergeticsCalculator cannot create a new InterfaceNeighborDefinitionCalculator named " <<
			name_of_InterfaceNeighborDefinitionCalculator_ << std::endl;
		utility_exit();
	}
	core::pose::metrics::PoseMetricCalculatorOP int_calculator =
		new InterfaceNeighborDefinitionCalculator(chain1_number,chain2_number);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( name_of_InterfaceNeighborDefinitionCalculator_, int_calculator );
}


// less preferred constructor - creates a new InterfaceNeighborDefinitionCalculator
InterfaceDeltaEnergeticsCalculator::InterfaceDeltaEnergeticsCalculator( char const chain1_letter, char const chain2_letter ) :
	EnergyDependentCalculator()
{
	name_of_InterfaceNeighborDefinitionCalculator_ = "IEC_private_IDC_" + to_string(chain1_letter) + "_" + to_string(chain2_letter);
	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( name_of_InterfaceNeighborDefinitionCalculator_ ) ) {
		basic::Error() << "InterfaceDeltaEnergeticsCalculator cannot create a new InterfaceNeighborDefinitionCalculator named " <<
			name_of_InterfaceNeighborDefinitionCalculator_ << std::endl;
		utility_exit();
	}
	core::pose::metrics::PoseMetricCalculatorOP int_calculator =
		new InterfaceNeighborDefinitionCalculator(chain1_letter,chain2_letter);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( name_of_InterfaceNeighborDefinitionCalculator_, int_calculator );
}


void InterfaceDeltaEnergeticsCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const {

	if ( key == "weighted_total" ) {
		basic::check_cast( valptr, &weighted_total_, "weighted_total expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( weighted_total_ );

	} else {
		scoring::ScoreType requested_scoretype( scoring::ScoreTypeManager::score_type_from_name( key ) );
		Real scoreval = delta_energies_unweighted_[ requested_scoretype ] * weights_[ requested_scoretype ];
		basic::check_cast( valptr, &scoreval, "Interface "+key+" expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( scoreval );
	}

}


std::string InterfaceDeltaEnergeticsCalculator::print( std::string const & key ) const {

	if ( key == "weighted_total" ) {
		return utility::to_string( weighted_total_ );
	} else {
		scoring::ScoreType requested_scoretype( scoring::ScoreTypeManager::score_type_from_name( key ) );
		Real scoreval = delta_energies_unweighted_[ requested_scoretype ] * weights_[ requested_scoretype ];
		return utility::to_string( scoreval );
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}


void InterfaceDeltaEnergeticsCalculator::recompute( Pose const & this_pose ) {

	// JK MAKE SURE THAT THE GRAPH STATE HERE IS "GOOD"

	// Get the first and last resnum of each chain, using name_of_InterfaceNeighborDefinitionCalculator_
	basic::MetricValue<Size> mv_size;
	this_pose.metric(name_of_InterfaceNeighborDefinitionCalculator_,"first_chain_first_resnum",mv_size);
	Size ch1_begin_num = mv_size.value();
	this_pose.metric(name_of_InterfaceNeighborDefinitionCalculator_,"first_chain_last_resnum",mv_size);
	Size ch1_end_num = mv_size.value();
	this_pose.metric(name_of_InterfaceNeighborDefinitionCalculator_,"second_chain_first_resnum",mv_size);
	Size ch2_begin_num = mv_size.value();
	this_pose.metric(name_of_InterfaceNeighborDefinitionCalculator_,"second_chain_last_resnum",mv_size);
	Size ch2_end_num = mv_size.value();

	// Clear the energy-holders, get the (unweighted) energies from the pose
	delta_energies_unweighted_.clear();
	scoring::EnergyGraph const & energy_graph( this_pose.energies().energy_graph() );

	// Loop over interactions across the interface
	for ( Size i = ch1_begin_num; i <= ch1_end_num; ++i ) {
		for ( graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
						irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
			const scoring::EnergyEdge * edge( static_cast< const scoring::EnergyEdge *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			if ( ( j >= ch2_begin_num ) && ( j <= ch2_end_num ) ) {
				delta_energies_unweighted_ += edge->fill_energy_map();
			}
		}
	}

	// Graph is asymmetric, so switch i/j and redo
	for ( Size i = ch2_begin_num; i <= ch2_end_num; ++i ) {
		for ( graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
						irue = energy_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {
			const scoring::EnergyEdge * edge( static_cast< const scoring::EnergyEdge *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			if ( ( j >= ch1_begin_num ) && ( j <= ch1_end_num ) ) {
				delta_energies_unweighted_ += edge->fill_energy_map();
			}
		}
	}

	//let's not forget the long range energies
	for( Size lr = 1; lr <= scoring::methods::n_long_range_types; lr++){
		scoring::methods::LongRangeEnergyType lr_type = scoring::methods::LongRangeEnergyType( lr );
		scoring::LREnergyContainerCOP lrec = this_pose.energies().long_range_container( lr_type );
		//runtime_assert( lrec );
		if( !lrec ) continue;
		if( lrec->empty() ) continue;

		for( Size i = ch1_begin_num; i <= ch1_end_num; ++i ){
			for ( scoring::ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( i );
					*rni != *( lrec->const_upper_neighbor_iterator_end( i ) );
					++(*rni) ) {
				Size j = rni->upper_neighbor_id();
				if( ( j >= ch2_begin_num ) && ( j <= ch2_end_num ) ) {
					scoring::EnergyMap emap;
					rni->retrieve_energy( emap );
					delta_energies_unweighted_ += emap;
				}
			}
		} //loop over chain 1 residues

		for( Size i = ch2_begin_num; i <= ch2_end_num; ++i ){
			for ( scoring::ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( i );
					*rni != *( lrec->const_upper_neighbor_iterator_end( i ) );
					++(*rni) ) {
				Size j = rni->upper_neighbor_id();
				if( ( j >= ch1_begin_num ) && ( j <= ch1_end_num ) ) {
					scoring::EnergyMap emap;
					rni->retrieve_energy( emap );
					delta_energies_unweighted_ += emap;
				}
			}
		} //loop over chain 2 residues

	} //loop over long range energy types

	// Save the most recently used weights, as well as the total (weighted) delta score
	weights_ = this_pose.energies().weights();
	weighted_total_ = delta_energies_unweighted_.dot(weights_);

	//debug stuff
	//std::cerr << "intef Delta E calc has following non-zero values: ";
	//for( Size ii = 1; ii <= scoring::n_score_types; ii++ ){
	//	scoring::ScoreType scotype = scoring::ScoreType( ii );
	//	if( weights_[scotype] != 0 ){
	//		std::cerr << scoring::name_from_score_type( scotype ) << " " << delta_energies_unweighted_[ scotype ] * weights_[ scotype ] << "; ";
	//	}
	//}
	//std::cerr << " you happy?" << std::endl;
	//debug stuff over
}


} // PoseMetricCalculators
} // toolbox
} // protocols
