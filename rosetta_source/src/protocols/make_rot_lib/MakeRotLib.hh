// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2008 University of Washington
// (C) 199x-2008 University of California Santa Cruz
// (C) 199x-2008 University of California San Francisco
// (C) 199x-2008 Johns Hopkins University
// (C) 199x-2008 University of North Carolina, Chapel Hill
// (C) 199x-2008 Vanderbilt University

#ifndef INCLUDED_protocols_make_rot_lib_makerotlib_HH
#define INCLUDED_protocols_make_rot_lib_makerotlib_HH

// core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// unit headers
#include <protocols/make_rot_lib/RotData.hh>

// c++ headers
#include <string>

namespace protocols {
namespace MakeRotLib {

// typedefs
typedef utility::vector1<RotData> RotVec;

// hack hack hack
void
asp_corrections( RotVec & );

void
glu_corrections( RotVec & );

void
phe_tyr_corrections( RotVec & );

void
peptoid_trans_hack( RotVec & );

void
peptoid_cis_hack( RotVec & );

// minimize side chain dihedral angles of each rotamer
void
min_rotamers (RotVec &, core::scoring::ScoreFunctionOP, std::string);

// fill in centroid with cmd-line opts
void
init_rotamers_centroids(RotVec &, RotVec &, core::Size &, std::string, std::string &, bool, core::Real, core::Real );

// itterates over distances to centroids, assigns to closest centroid
// assigns rots to clusters and determins if clusters changed
bool
calc_rotamer_clusters (RotVec &);

// find new centroids
bool
calc_centroids (RotVec &, RotVec &);

//calculates distance b/w 2 RotData objects
core::Real
calc_dist (RotData &, RotData &);

//calculates avg distance between centroid and all points in its cluster
core::Real
avg_cluster_cen_dist (RotVec &, core::Size &);

//finds all distances between all rotamers and all centroids
void
calc_all_dist (RotVec &, RotVec &);

// pull out best rots from rotamers and add them to final_rotamers
void
get_final_rots(RotVec &, RotVec &, core::Size &);

// calc probabilites for final rots & normilize
void
get_final_rot_probs( RotVec &);

void
calc_std_dev (RotVec &, core::scoring::ScoreFunctionOP, std::string );

// prints out rotamers
void
pretty_print_rd( RotData & );

// print out rotamers acording to the Dunbrack format
void
dunbrack_print( RotVec &, RotVec &, std::string );

} // namespace MakeRotLib
} // namespace protocols

#endif // INCLUDED_protocols_makerotlib_makerotlib_HH
