// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/VoxelGrid.hh
/// @brief      Neighbour class
/// @details    ... TO FILL
/// @author     Mads Jeppesen

#ifndef INCLUDED_core_scoring_shape_Neighbour_hh
#define INCLUDED_core_scoring_shape_Neighbour_hh

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/pointer/owning_ptr.hh>

// core
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/chemical/AtomType.hh>

// protocols
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// stl
#include <math.h>       /* ceil  and fabs */
#include <algorithm>   /* sort */
#include <fstream>     /* ofstream */
#include <chrono>
#include <ctime>
#include <vector>

/// @brief Neighbour object that is used to store 1 neighbour position as well as the distance to the neighbour.
class Neighbour {

public:

    /// @brief constructor.
    Neighbour(int x, int y, int z, core::Real distance) {
        _x = x,
        _y = y,
        _z = z;
        equal_grid_distance_ = distance;
    };

    /// @brief destructor
    ~Neighbour() = default;

    /// @brief x position.
    int x();

    /// @brief y position.
    int y();

    /// @brief z position.
    int z();

    core::Real
    equal_grid_distance();

    /// @brief distance to neighbour.
    core::Real
    grid_resolution_distance(const core::Real & grid_resolution);

private:
    int _x,_y,_z;
    core::Real equal_grid_distance_;
};


#endif // INCLUDED_core_scoring_shape_Neighbour_hh
