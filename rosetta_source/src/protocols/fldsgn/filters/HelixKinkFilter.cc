// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/HelixKinkFilter.cc
/// @brief
/// @detailed filter structures out, which have kink helices
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/HelixKinkFilter.hh>
#include <protocols/fldsgn/filters/HelixKinkFilterCreator.hh>

// Package Headers
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/util.hh>

// Project Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

// Parser headers
#include <utility/tag/Tag.hh>

// AUTO-REMOVED #include <protocols/fldsgn/topology/HSSTriplet.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#ifdef WIN32
	#include <protocols/fldsgn/topology/HSSTriplet.hh>
#endif


//// C++ headers
static basic::Tracer TR("protocols.fldsgn.filters.HelixKinkFilter");

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
HelixKinkFilter::HelixKinkFilter():
	Filter( "HelixKink" ),
	bend_angle_( 20.0 ),
	secstruct_( "" )
{}


// @brief copy constructor
HelixKinkFilter::HelixKinkFilter( HelixKinkFilter const & rval ):
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	bend_angle_( rval.bend_angle_ ),
	secstruct_( rval.secstruct_ )
{}


/// @brief
bool
HelixKinkFilter::apply( Pose const & pose ) const
{
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;
	using protocols::fldsgn::topology::Helix;
	using protocols::fldsgn::topology::Helices;
	using protocols::fldsgn::topology::check_kink_helix;
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;

	// set SS_Info
	String secstruct( pose.secstruct() );
	if( secstruct_ != "" ) {
		secstruct = secstruct_;
	}
	SS_Info2_OP  ss_info = new SS_Info2( pose, secstruct );
	Helices const & helices( ss_info->helices() );

	if( helices.size() < 1 ) {
		TR << "There is no helix definition in pose. " << std::endl;
		return true;
	}

	if ( ! pose.energies().data().has( HBOND_SET ) ) {
		TR << " Pose does not have HBOND_SET. Checking hbonds will be skipped. " << std::endl;
	}

	// check kink
	for( Size ii=1; ii<=helices.size(); ++ii ) {

		TR << "Helix " << ii << ", res " << helices[ ii ]->begin() << "-" << helices[ ii ]->end() << ", ";
		// check helix bend
		if ( helices[ ii ]->bend() > bend_angle_ ) {
			TR << "is bended angle=" << helices[ ii ]->bend() << std::endl;
			return false;
		}

		// check broken hydrogen within helix
		/// @brief check kink of helix, return number of loosen hydrogen
		if ( pose.energies().data().has( HBOND_SET ) ) {
			Size broken_hbonds( check_kink_helix( pose, helices[ ii ]->begin()-1, helices[ ii ]->end()-5 ) );
			if( broken_hbonds > 0 ) {
				TR << "is kinked, hbonds are broken. " << std::endl;
				return false;
			}
		}

		TR << "is OK." << std::endl;
	}

	// check broken hydrogen bond
	TR << " Filter success ! " << std::endl;

	return true;

}


/// @brief parse xml
void
HelixKinkFilter::parse_my_tag(
	TagPtr const tag,
	DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	bend_angle_  = tag->getOption<Real>( "bend",  20 );
	// secondary strucuture info
	String const blueprint = tag->getOption<String>( "blueprint", "" );
	if( blueprint != "" ) {
		protocols::jd2::parser::BluePrint blue( blueprint );
		secstruct_ = blue.secstruct();
	}
}

protocols::filters::FilterOP
HelixKinkFilterCreator::create_filter() const { return new HelixKinkFilter; }

std::string
HelixKinkFilterCreator::keyname() const { return "HelixKink"; }

} // filters
} // fldsgn
} // protocols
