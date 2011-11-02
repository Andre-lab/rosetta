// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/flxbb/SelectResiduesByLayer.cc
/// @brief select residues depending on layer: core, boundary, and surface
/// The layer of residues are defined by accessible sufrace area of CB + mainchain
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/flxbb/SelectResiduesByLayer.hh>

// Project Headers
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/sasa.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <fstream>

using basic::T;
using basic::Error;
using basic::Warning;


static basic::Tracer TR("protocols.flxbb.SelectResiduesByLayer");

namespace protocols {
namespace flxbb {

/// @brief default constructor
SelectResiduesByLayer::SelectResiduesByLayer():
	pick_core_( false ),
	pick_boundary_( false ),
	pick_surface_( false ),
	pore_radius_( 2.0 ),
	make_rasmol_format_file_( false )
{
	initialize( 20.0, 40.0 );
}


/// @brief value constructor
SelectResiduesByLayer::SelectResiduesByLayer( bool const pick_core, bool const pick_boundary, bool const pick_surface ):
	pick_core_( pick_core ),
	pick_boundary_( pick_boundary ),
	pick_surface_( pick_surface ),
	pore_radius_( 2.0 ),
	make_rasmol_format_file_( false )
{
	initialize( 20.0, 40.0 );
}

/// @brief value constructor
SelectResiduesByLayer::SelectResiduesByLayer( String const pick ):
	pick_core_( false ),
	pick_boundary_( false ),
	pick_surface_( false ),
	pore_radius_( 2.0 ),
	make_rasmol_format_file_( false )
{
	initialize( 20.0, 40.0 );
	utility::vector1< String > layers( utility::string_split( pick, '_' ) );
	for( utility::vector1< String >::const_iterator iter = layers.begin(); iter != layers.end() ; ++iter) {
		String layer(*iter);
		if ( layer == "core" ) {
			pick_core_ = true;
		} else if ( layer == "surface" ) {
			pick_surface_ = true;
		} else if ( layer == "boundary" ) {
			pick_boundary_ = true;
		} else {
			TR << "Error!, wrong specification of layer_mode " << layer << std::endl;
			runtime_assert( false );
		}
	} // utility::vector1
}

/// @brief destructor
SelectResiduesByLayer::~SelectResiduesByLayer() {}

/// @brief copy constructor
SelectResiduesByLayer::SelectResiduesByLayer( SelectResiduesByLayer const & rval ):
	pick_core_( rval.pick_core_ ),
	pick_boundary_( rval.pick_boundary_ ),
	pick_surface_( rval.pick_surface_ ),
	pore_radius_( rval.pore_radius_ ),
	burial_( rval.burial_ ),
	surface_( rval.surface_ ),
	excluded_aatypes_for_selection_( rval.excluded_aatypes_for_selection_ ),
	selected_core_residues_( rval.selected_core_residues_ ),
	selected_boundary_residues_( rval.selected_boundary_residues_ ),
	selected_surface_residues_( rval.selected_surface_residues_ ),
	make_rasmol_format_file_( rval.make_rasmol_format_file_ ),
	rsd_sasa_( rval.rsd_sasa_ ),
	rsd_layer_( rval.rsd_layer_ )
{}

/// @brief accessible surface for evaluating residues are in surface or not
void
SelectResiduesByLayer::sasa_surface( Real const r, String const ss )
{
	if ( ss == "" ) {
		surface_[ 'H' ] = r;
		surface_[ 'L' ] = r;
		surface_[ 'E' ] = r;
	} else {
		runtime_assert( ss.length() == 1 );
		runtime_assert( ss.at( 0 ) == 'L' || ss.at( 0 ) == 'E' || ss.at( 0 ) == 'H' );
		surface_[ ss.at( 0 ) ] = r;
	}
}

/// @brief accessible surface for evaluating residues are in core or not
void
SelectResiduesByLayer::sasa_core( Real const r, String const ss )
{
	if ( ss == "" ) {
		burial_[ 'H' ] = r;
		burial_[ 'L' ] = r;
		burial_[ 'E' ] = r;
	} else {
		runtime_assert( ss.length() == 1 );
		runtime_assert( ss.at( 0 ) == 'L' || ss.at( 0 ) == 'E' || ss.at( 0 ) == 'H' );
		burial_[ ss.at( 0 ) ] = r;
	}
}

/// @brief
void
SelectResiduesByLayer::initialize( Real burial, Real surface )
{
	surface_.insert( std::map< char, Real >::value_type( 'H', surface ) );
	surface_.insert( std::map< char, Real >::value_type( 'L', surface ) );
	surface_.insert( std::map< char, Real >::value_type( 'E', surface ) );
	burial_.insert( std::map< char, Real >::value_type( 'H', burial ) );
	burial_.insert( std::map< char, Real >::value_type( 'L', burial ) );
	burial_.insert( std::map< char, Real >::value_type( 'E', burial ) );
}


/// @brief amino acid types excluded for selection
void
SelectResiduesByLayer::exclude_aatypes_for_selection( utility::vector1< AA > const & aas )
{
	excluded_aatypes_for_selection_ = aas;
}


/// @brief amino acid types restricted for selection
void
SelectResiduesByLayer::restrict_aatypes_for_selection( utility::vector1< AA > const & aas )
{
	restricted_aatypes_for_selection_ = aas;
}


/// @brief accessbile surface are of each residue
core::Real
SelectResiduesByLayer::rsd_sasa( Size const i ) const
{
	return rsd_sasa_[ i ];
}

/// @brief return defined layer for each residue
std::string
SelectResiduesByLayer::layer( Size const i ) const
{
	return rsd_layer_[ i ];
}

/// @brief selected residues on boundary
utility::vector1< core::Size > const &
SelectResiduesByLayer::selected_boundary_residues() const
{
	return selected_boundary_residues_;
}


/// @brief selected residues in core
utility::vector1< core::Size > const &
SelectResiduesByLayer::selected_core_residues() const
{
	return selected_core_residues_;
}


/// @brief selected residues on surface
utility::vector1< core::Size > const &
SelectResiduesByLayer::selected_surface_residues() const
{
	return selected_surface_residues_;
}


/// @brief return accessible surface area for each residue
utility::vector1< core::Real > const
SelectResiduesByLayer::calc_rsd_sasa( Pose const & pose ) const {

	// define atom_map for main-chain and CB
	core::id::AtomID_Map< bool > atom_map;
	core::pose::initialize_atomid_map( atom_map, pose, false );
	for ( Size ir = 1; ir <= pose.total_residue(); ++ir ) {
		for ( Size j = 1; j<=5; ++j ) {
			core::id::AtomID atom( j, ir );
			atom_map.set( atom, true );
		}
	}

	// calc sasa
	core::id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, pore_radius_, false, atom_map );

	return rsd_sasa;
} // calc_residue_sasa


/// @brief
utility::vector1< core::Size > const
SelectResiduesByLayer::compute( Pose const & pose, String secstruct )
{
	// define secstruct if it is empty
	if( secstruct == "" ) {
		secstruct = pose.secstruct();
	}
	runtime_assert( pose.total_residue() == secstruct.length() );

	// clear
	rsd_sasa_.clear();
	rsd_layer_.clear();
	selected_core_residues_.clear();
	selected_surface_residues_.clear();
	selected_boundary_residues_.clear();

	std::ofstream output;
	if( make_rasmol_format_file_ ) {
		output.open( "srb.ras" ,std::ios::out );
	}

	TR << " pore_radius : " <<  pore_radius_ << std::endl;
	TR << " core ( E, L, H ): " << burial_[ 'E' ] << ' ' << burial_[ 'L' ] << ' ' << burial_[ 'H' ] << std::endl;
	TR << " surface (E, L, H ): " << surface_[ 'E' ] << ' ' << surface_[ 'L' ] << ' ' << surface_[ 'H' ] << std::endl;

	rsd_sasa_ = calc_rsd_sasa( pose );
	rsd_layer_.resize( pose.total_residue() );

	utility::vector1< Size > selected_residues;
	for( Size iaa=1; iaa<=pose.total_residue(); iaa++ ) {

		char ss = secstruct.at( iaa-1 );
		runtime_assert( ss == 'L' || ss =='E' || ss=='H' );

		rsd_layer_[ iaa ] = "";

		bool flag( false );
		for( Size i=1; i<=excluded_aatypes_for_selection_.size(); i++ ) {
			if( pose.aa( iaa ) == excluded_aatypes_for_selection_[ i ] ) {
				flag = true;
				break;
			}
		}
		if( flag ) continue;

		if( restricted_aatypes_for_selection_.size() > 0 ) flag = true;
		for( Size i=1; i<=restricted_aatypes_for_selection_.size(); i++ ) {
			if( pose.aa( iaa ) == restricted_aatypes_for_selection_[ i ] ) {
				flag = false;
				break;
			}
		}
		if( flag ) continue;

		if ( rsd_sasa_[ iaa ] <= burial_[ ss ] && pick_core_ ) {

			selected_core_residues_.push_back( iaa );
			selected_residues.push_back( iaa );

			if( make_rasmol_format_file_ ) {
				output << "select " << iaa << std::endl;
				output << "color blue " << std::endl;
			}

			rsd_layer_[ iaa ] = "core";

		} else if ( rsd_sasa_[ iaa ] >= surface_[ ss ] && pick_surface_ ) {

			selected_surface_residues_.push_back( iaa );
			selected_residues.push_back( iaa );

			if( make_rasmol_format_file_ ) {
				output << "select " << iaa << std::endl;
				output << "color red " << std::endl;
			}

			rsd_layer_[ iaa ] = "surface";

		} else if ( rsd_sasa_[ iaa ] < surface_[ ss ] && rsd_sasa_[ iaa ] > burial_[ ss ] && pick_boundary_ ) {

			selected_boundary_residues_.push_back( iaa );
			selected_residues.push_back( iaa );

			if( make_rasmol_format_file_ ) {
				output << "select " << iaa << std::endl;
				output << "color green " << std::endl;
			}

			rsd_layer_[ iaa ] = "boundary";

		}

	}

	return selected_residues;

} // compute



} // flxbb
} // protocols
