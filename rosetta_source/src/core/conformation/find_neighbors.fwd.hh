// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/find_neighbors.fwd.hh
/// @brief  forward headers for find_neighbors.hh
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
///
/// @remarks Thanks to Will Sheffler for his ideas on refining this and extending it to atom neighbors
/// @remarks Adapting libRosetta code for generalized neighbor detection

#ifndef INCLUDED_core_conformation_find_neighbors_fwd_hh
#define INCLUDED_core_conformation_find_neighbors_fwd_hh

// Package Headers
// AUTO-REMOVED #include <core/conformation/PointGraph.fwd.hh>
#include <core/types.hh>



// Numeric headers
// AUTO-REMOVED #include <numeric/numeric.functions.hh>
// AUTO-REMOVED #include <numeric/xyzTriple.hh>
// AUTO-REMOVED #include <numeric/xyzVector.hh>

// ObjexxFCL headers
//#include <ObjexxFCL/KeyFArray1D.hh>
//#include <ObjexxFCL/KeyFArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray3D.hh>

// Utility headers
//#include <utility/pointer/access_ptr.hh>

// boost headers
// AUTO-REMOVED #include <boost/unordered_map.hpp>

// C++ headers
#include <cassert>
// AUTO-REMOVED #include <cmath>
// AUTO-REMOVED #include <cstdlib>
#include <limits>
// AUTO-REMOVED #include <map>
#include <vector>

#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/xyzTriple.fwd.hh>


namespace core {
namespace conformation {

enum Strategy {
	NAIVE,
	AUTOMATIC,
	OCTREE,
	THREEDGRID
};

// move the following typedef to top of file instead of
// within find_neighbors()
typedef numeric::xyzTriple< core::Size > CubeKey;
typedef numeric::xyzVector< core::Real > PointPosition;

template <class Vertex, class Edge>
void
find_neighbors_naive(
	utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff
);

template <class Vertex, class Edge>
void
find_neighbors_octree(
	utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	Strategy strategy
);

template <class Vertex, class Edge>
void
find_neighbors_3dgrid(
	utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff
);

template <class Vertex, class Edge>
void
find_neighbors_naive_restricted(
	utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const & residue_selection
);

template <class Vertex, class Edge>
void
find_neighbors_octree_restricted(
	utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const & residue_selection,
	Strategy strategy
);


template <class Vertex, class Edge>
void
find_neighbors_3dgrid_restricted(
	utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
	core::Real neighbor_cutoff,
	utility::vector1< bool > const &  residue_selection
);

template <class Vertex, class Edge>
core::Size
get_nearest_neighbor(
		utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
		core::Size node_id,
		core::Real neighbor_cutoff,
		Strategy strategy = AUTOMATIC
);

/*
// Commented out to make clang compile - duplication of default arguments and forward declaration of template functions confuses the compiler (and me!)
// Brian Weitzner and Sergey Lyskov 3/5/2011
 
template <class Vertex, class Edge>
void
find_neighbors(
    utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
    core::Real neighbor_cutoff,
    Strategy strategy = AUTOMATIC
);
 
template <class Vertex, class Edge>
void
find_neighbors_restricted(
    utility::pointer::owning_ptr<graph::UpperEdgeGraph<Vertex, Edge> > point_graph,
    core::Real neighbor_cutoff,
    utility::vector1< bool > const & residue_selection,
    Strategy strategy = AUTOMATIC
);
*/


}
}

#endif
