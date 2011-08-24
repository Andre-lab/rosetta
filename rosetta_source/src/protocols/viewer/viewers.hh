// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_protocols_viewer_viewers_hh
#define INCLUDED_protocols_viewer_viewers_hh

// Unit headers

// Package headers
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

#include <protocols/viewer/GraphicsState.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/viewer/triangle.hh>

// Project headers
#include <utility/vector1.hh>

// C++ Headers
#include <string>

//Auto Headers
#include <iostream>


namespace protocols {
namespace viewer {

typedef void* (*VoidFunc)(void*);

static std::string empty_string("");

#ifndef GL_GRAPHICS ///////////////////////////////////////////////////////

inline
void
add_conformation_viewer(
	core::conformation::Conformation &,
	std::string const = empty_string,
	int const = 900,
	int const = 900,
	bool = false
)
{
}

inline
void
add_monte_carlo_viewer(
	moves::MonteCarlo &,
	std::string const = empty_string,
	int const = 900,
	int const = 900,
	bool = false
)
{
}

inline
int
viewer_main( VoidFunc worker_main ){ worker_main(NULL); return 0; }

void set_bg_color( core::Vector new_bg_color );


#else // GL_GRAPHICS ////////////////////////////////////////////////////

void
add_conformation_viewer(
	core::conformation::Conformation & conformation,
	std::string const name_in = empty_string,
	int const length = 900,
	int const width  = 900,
	bool debug_pause=false
);

void
add_monte_carlo_viewer(
	moves::MonteCarlo & mc,
	std::string const name_in = empty_string,
	int const length = 900,
	int const width  = 900,
	bool debug_pause=false
);


int
viewer_main( VoidFunc worker_main );

#endif

#if defined GL_GRAPHICS || BOINC_GRAPHICS

void
display_residues(
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	core::id::AtomID const & anchor_id
);

void
display_residues_wireframe(
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	core::id::AtomID const & anchor_id
);

void
display_residues_wireframe(
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	core::Vector const & center
);

void set_bg_color( core::Vector new_bg_color );


void draw_pose(const core::pose::Pose & pose,
                GraphicsState & gs, bool centered);

void draw_conformation( utility::vector1< core::conformation::ResidueCOP > const & residues,
                        utility::vector1< char > const & ss,
                        GraphicsState & gs, core::Vector const & center);

void
render_density(
			GraphicsState &gs,
			utility::vector1< triangle > &triangles );

void
display_density(
			GraphicsState &gs,
			utility::vector1< triangle > &triangles ) ;

void
draw_conformation_and_density(
			utility::vector1< core::conformation::ResidueCOP > const & residues,
			utility::vector1< char > const & ss,
			utility::vector1< triangle > &triangles,
			GraphicsState & gs,
			core::Vector const & center);


#endif ////////////////////////////////////////////////////////////////

void
add_monte_carlo_silent_viewer(
	moves::MonteCarlo & mc,
	std::string const name_in,
	bool fullatom
);

/// @brief Allows for graceful exit of graphics viewers.
void
clear_conformation_viewers();


} // viewer
} // protocols


#endif
