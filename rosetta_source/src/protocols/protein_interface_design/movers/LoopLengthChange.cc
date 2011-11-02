// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/LoopLengthChange.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChangeCreator.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace std;
using namespace core::scoring;

static basic::Tracer TR( "protocols.protein_interface_design.movers.LoopLengthChange" );

std::string
LoopLengthChangeCreator::keyname() const
{
	return LoopLengthChangeCreator::mover_name();
}

protocols::moves::MoverOP
LoopLengthChangeCreator::create_mover() const {
	return new LoopLengthChange;
}

std::string
LoopLengthChangeCreator::mover_name()
{
	return "LoopLengthChange";
}

LoopLengthChange::LoopLengthChange() :
	Mover( LoopLengthChangeCreator::mover_name() ),
	loop_start_( 0 ), loop_end_( 0 ), delta_( 0 )
{
}


LoopLengthChange::~LoopLengthChange() {}

void
LoopLengthChange::apply( core::pose::Pose & pose )
{
	TR<<"Changing loop "<<loop_start()<<"-"<<loop_end()<<" by "<<delta()<<std::endl;
	runtime_assert( delta() != 0 );
  runtime_assert( loop_end() > loop_start() );
  runtime_assert( loop_end() + delta() >= loop_start() );
  if( delta() < 0 ){
    for( int del(-1); del>=delta(); --del )
      pose.delete_polymer_residue( loop_start() + 1 );
  }
  else if( delta() > 0 ){
    using namespace core::chemical;
    using namespace core::conformation;

    ResidueTypeSet const & residue_set( pose.residue( 1 ).residue_type_set() ); // residuetypeset is noncopyable
    ResidueCOP new_res = ResidueFactory::create_residue( residue_set.name_map( name_from_aa( aa_from_oneletter_code( 'A' ) ) ) );
    for( core::Size leng(1); leng<=(core::Size) delta(); ++leng )
      pose.conformation().safely_append_polymer_residue_after_seqpos( *new_res, loop_start() + 1, true/*build_ideal
_geometry*/ );
  }
  pose.update_residue_neighbors();
}

std::string
LoopLengthChange::get_name() const {
	return LoopLengthChangeCreator::mover_name();
}

void
LoopLengthChange::parse_my_tag( TagPtr const tag, protocols::moves::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{
	using namespace protocols::rosetta_scripts;

	loop_start( parse_resnum( tag->getOption< std::string >( "loop_start" ), pose ) );
	loop_end( parse_resnum( tag->getOption< std::string >( "loop_end" ), pose ) );
	delta( tag->getOption< int >( "delta" ) );

  runtime_assert( delta() != 0 );
  runtime_assert( loop_end() > loop_start() );
  runtime_assert( loop_end() + delta() >= loop_start() );

	TR<<"LoopLengthChange with loop "<<loop_start()<<"-"<<loop_end()<<" and delta "<<delta()<<std::endl;
}

protocols::moves::MoverOP
LoopLengthChange::clone() const {
    return( protocols::moves::MoverOP( new LoopLengthChange( *this ) ));
}

void
LoopLengthChange::loop_start( core::Size const loop_start ){
	loop_start_ = loop_start;
}

core::Size
LoopLengthChange::loop_start() const{
  return( loop_start_ );
}

void
LoopLengthChange::loop_end( core::Size const l ){
  loop_end_ = l;
}

core::Size
LoopLengthChange::loop_end() const{
  return( loop_end_ );
}

void
LoopLengthChange::delta( int const d ){
  delta_ = d;
}

int
LoopLengthChange::delta() const{
  return( delta_ );
}

} //movers
} //protein_interface_design
} //protocols
