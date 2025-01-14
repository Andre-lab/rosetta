// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/VoxelGrid.cc
/// @brief      Neighbour class
/// @details    ... TO FILL
/// @author     Mads Jeppesen

#include <core/scoring/shape/Neighbour.hh>

int Neighbour::x() {
    return _x;
}
int Neighbour::y(){
    return _y;
}
int Neighbour::z() {
    return _z;
}

core::Real Neighbour::equal_grid_distance() {
    return equal_grid_distance_;
}

core::Real Neighbour::grid_resolution_distance(const core::Real & grid_resolution) {
    return grid_resolution * equal_grid_distance_;
}
