// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

// Unit headers
#include <protocols/moves/SetTorsion.hh>
#include <protocols/moves/SetTorsionCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/util.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/ResidueFactory.hh>
//parsing
#include <utility/tag/Tag.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


// Utility Headers

// Unit Headers

// C++ headers

namespace protocols {
namespace moves {

using namespace core;
using namespace core::chemical;
using namespace std;

using core::pose::Pose;
using core::conformation::Residue;

static basic::Tracer TR( "protocols.moves.SetTorsion" );

std::string
SetTorsionCreator::keyname() const
{
	return SetTorsionCreator::mover_name();
}

protocols::moves::MoverOP
SetTorsionCreator::create_mover() const {
	return new SetTorsion;
}

std::string
SetTorsionCreator::mover_name()
{
	return "SetTorsion";
}

SetTorsion::~SetTorsion() {}

///@brief default ctor
SetTorsion::SetTorsion() :
	parent(),
	angle_( 0 ),
	resnum_( 0 ),
	torsion_name_( "" )
{}

void SetTorsion::apply( Pose & pose ) {
	runtime_assert( resnum() > 0 );
	runtime_assert( resnum() <= pose.total_residue() );
	runtime_assert( torsion_name() == "phi" || torsion_name() == "psi" );

	if( torsion_name() == "phi" )
		pose.set_phi( resnum(), angle() );
	else if( torsion_name() == "psi" )
		pose.set_psi( resnum(), angle() );
}

std::string
SetTorsion::get_name() const {
	return SetTorsionCreator::mover_name();
}

void SetTorsion::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		Pose const & pose)
{
	angle( tag->getOption< core::Real >( "angle" ) );
	resnum( protocols::rosetta_scripts::parse_resnum( tag->getOption< std::string >( "resnum" ), pose ) );
	torsion_name( tag->getOption< std::string >( "torsion_name" ) );
	TR<<"Set torsion "<<torsion_name_<<" at residue "<< resnum_ << " to "<<angle_<<std::endl;
}


} // moves
} // protocols
