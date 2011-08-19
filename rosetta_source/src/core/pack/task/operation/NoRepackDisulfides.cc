// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/NoRepackDisulfides.cc
/// @brief  prevent disulfides from being repacked; assumes disulfide info in
///         Pose is up-to-date
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/NoRepackDisulfidesCreator.hh>

// package headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// C++ headers
// AUTO-REMOVED #include <set>


namespace core {
namespace pack {
namespace task {
namespace operation {


// static
static basic::Tracer TR( "core.pack.task.operation.NoRepackDisulfides" );


/// @brief default constructor
NoRepackDisulfides::NoRepackDisulfides() :
	Super()
{}


/// @brief copy constructor
NoRepackDisulfides::NoRepackDisulfides( NoRepackDisulfides const & rval ) :
	Super( rval )
{}


/// @brief default destructor
NoRepackDisulfides::~NoRepackDisulfides() {}

TaskOperationOP NoRepackDisulfidesCreator::create_task_operation() const
{
	return new NoRepackDisulfides;
}

/// @brief clone this object
NoRepackDisulfides::TaskOperationOP NoRepackDisulfides::clone() const {
	return new NoRepackDisulfides( *this );
}

/// @brief apply operations to PackerTask
void NoRepackDisulfides::apply( Pose const & pose, PackerTask & task ) const {
	using core::Size;
	using core::chemical::aa_cys;
	using core::chemical::DISULFIDE;
	using core::conformation::Residue;

	for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
		Residue const & res = pose.residue( i );

		// residue appears to be a disulfide...?
		if ( res.aa() == aa_cys && res.has_variant_type( DISULFIDE ) ) {

			bool no_repack = true;

			// try and double-check if possible
			if ( res.type().has_atom_name( "SG" ) ) {
				// check its partner to see if it's really the case
				Size const sg_index = res.atom_index( "SG" );
				Size const sg_conn_index = res.type().residue_connection_id_for_atom( sg_index );
				Residue const & partner_res = pose.residue( res.residue_connection_partner( sg_conn_index ) );

				if ( partner_res.aa() != aa_cys || !partner_res.has_variant_type( DISULFIDE ) ) {
					no_repack = false;
				}
			}

			// set repack status
			if ( no_repack ) {
				TR.Debug << "found disulfide residue " << i << ", preventing repack at this position" << std::endl;
				task.nonconst_residue_task( i ).prevent_repacking();
			} else {
				TR.Warning << "WARNING: residue " << i << " marked as disulfide but has no partner, allowing repack at this position" << std::endl;
			}

		}

	} // foreach residue

}


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core

