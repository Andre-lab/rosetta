// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/VoxelGrid.cc
/// @brief      VoxelGrid class
/// @details    ... TO FILL
/// @author     Mads Jeppesen
/// @author     Ingemar Andre

// mads
#include <core/scoring/shape/VoxelGrid.hh>
#include <core/scoring/shape/zernike/ZernikeDescriptor.hh>
#include <core/scoring/shape/zernike/Zernike2Ddescriptor.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
//#include <core/scoring/shape/Neighbour.hh>
#include <core/types.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
//#include <math.h>       /* pow */
#include <stack>
#include <cmath>
//#include <external/pybind11/include/pybind11/pybind11.h>
//#include <external/pybind11/include/pybind11/numpy.h>
#include <core/kinematics/Jump.hh>

#include <basic/Tracer.hh>
#include <vector>
#include <algorithm>
#include <iterator>
/// ingemar
//#include <core/types.hh>
//#include <queue>
//#include <basic/options/keys/edensity.OptionKeys.gen.hh>
//#include <core/scoring/electron_density/ElectronDensity.hh>
#include <basic/options/keys/zernike_descriptor.OptionKeys.gen.hh>
#include <iomanip>
#include <ctime>
#include <chrono>

namespace core {
namespace scoring {
namespace shape {

static basic::Tracer TR( "core.scoring.shape.VoxelGrid" );

// todo: so this is called when we initialize the farray - but can we initialize w. below?
Voxel::Voxel() = default;

Voxel::Voxel(int x, int y, int z):
        filled_(false),
        visited_by_fs_(false),
        surface_by_fs_(false),
        visited_by_edt_(false),
        x_(x),
        y_(y),
        z_(z),
        edt_root_x_(x),
        edt_root_y_(y),
        edt_root_z_(z),
        edt_dist_(std::numeric_limits<int>::max())
    {}

void
Voxel::fill() {
    filled_ = true;
}

void
Voxel::unfill() {
    filled_ = false;
}

bool
Voxel::is_filled() {
    return filled_;
}

void
Voxel::visit_by_fs() {
    visited_by_fs_ = true;
}

bool
Voxel::is_visited_by_fs() {
    return visited_by_fs_;
}

void
Voxel::unvisit_by_fs() {
    visited_by_fs_ = false;
}

void
Voxel::surface_by_fs() {
    surface_by_fs_ = true;
}

void
Voxel::unsurface_by_fs() {
    surface_by_fs_ = false;
}

bool
Voxel::is_surface_by_fs() {
    return surface_by_fs_;
}

void
Voxel::visit_by_edt() {
    visited_by_edt_ = true;
}

bool
Voxel::is_visited_by_edt() {
    return visited_by_edt_;
}

void
Voxel::unvisit_by_edt() {
    visited_by_edt_ = false;
}

void
Voxel::surface_by_edt() {
    surface_by_edt_ = true;
}

void
Voxel::unsurface_by_edt() {
    surface_by_edt_ = false;
}

bool
Voxel::is_surface_by_edt() {
    return surface_by_edt_;
}

int Voxel::x() const{
    return x_;
}

int Voxel::y() const {
    return y_;
}

int Voxel::z() const {
    return z_;
}

core::Real
Voxel::edt_dist() const {
    return edt_dist_;
}

int
Voxel::edt_root_x() const {
    return edt_root_x_;
}

int
Voxel::edt_root_y() const {
    return edt_root_y_;
}

int
Voxel::edt_root_z() const {
    return edt_root_z_;
}

void
Voxel::edt_root_from_voxel(Voxel const & center) {
    edt_root_x_ = center.edt_root_x();
    edt_root_y_ = center.edt_root_y();
    edt_root_z_ = center.edt_root_z();
}

core::Real
Voxel::edt_root_distance(Voxel const & neighbour) const {
    core::Real xd = edt_root_x_ - neighbour.x();
    core::Real yd = edt_root_y_ - neighbour.y();
    core::Real zd = edt_root_z_ - neighbour.z();
    return sqrt(xd*xd + yd*yd + zd*zd);
}

void
Voxel::set_edt_dist(core::Real dst) {
    edt_dist_ = dst;
}


void
Voxel::clear() {
    unfill();
    unvisit_by_fs();
    unsurface_by_fs();
    unvisit_by_edt();
    unsurface_by_edt();
    edt_root_x_ = x_;
    edt_root_y_ = y_;
    edt_root_z_ = z_;
    edt_dist_ = std::numeric_limits<int>::max();
}

bool
Voxel::operator==(const Voxel v) const {
	return (v.x() == x_ && v.y() == y_ && v.z() == z_ );
}

VoxelGrid::VoxelGrid(
        int size,
        std::string const & surface_type,
        core::Real probe,
        core::Real shell,
        bool recalculate,
        bool scale_on,
        core::Real scaling,
        bool never_leave_neighbour):
        grid_dimension_(size),
        actual_grid_dimension_(size + 2),
        geometric_center_at_grid_center_((size % 2) / 2),
        scaling_(scaling),
        scale_on_(scale_on),
        probe_radius_(probe),
        shell_thickness_(shell),
		never_leave_neighbour_(never_leave_neighbour),
        recalculate_(recalculate),
        padding_(1.1){
    amount_of_voxel_neighbours_ = std::floor(size / 2);
    set_surface_type(surface_type);
    initialize_grid();
    }


void
VoxelGrid::initialize_grid() {
    generate_neighbours();
    grid_.dimension(actual_grid_dimension_, actual_grid_dimension_, actual_grid_dimension_); // + 2 because we need to have space for the visited voxels
    for (int x = 1; x <= actual_grid_dimension_; x++) {
        for (int y = 1; y <= actual_grid_dimension_; y++) {
            for (int z = 1; z <= actual_grid_dimension_; z++)
                grid_(x,y,z) = Voxel(x,y,z);
        }
    }
}


void
VoxelGrid::set_scaling(core::Real scaling) {
    assert(scaling_ <= 1);
    scaling_ = scaling;
    scale_on_ = true;
}

void
VoxelGrid::set_grid_size(int size) {
    grid_dimension_ = size;
    actual_grid_dimension_ = size + 2;
    geometric_center_at_grid_center_ = (size % 2) / 2;
    initialize_grid();
}

void
VoxelGrid::set_surface_type(std::string const & surface_type) {
    if (surface_type == "MS")
        surface_type_ = 1;
    else if (surface_type == "SAS")
        surface_type_ = 2;
    else if (surface_type == "VDW")
        surface_type_ = 3;
    else
        utility_exit_with_message("Did not understand surface_type. 'MS', 'SAS' and 'VDW' is understood");
}

core::scoring::shape::Voxel &
VoxelGrid::operator()(int x, int y, int z) {
    return grid_(x, y, z);
}

void
VoxelGrid::voxelize(std::vector<std::vector<core::Real>> atom_coords, std::vector<core::Real> ljs) {

	core::Real padding = padding_;

	// gogo label if grid is not big enough
	// todo: probably not considered good to use goto so can be depreciated in the future - can make a wrapper function that returns a bool if fail!
	beginning:

	// clean the grid (cleans the grid_ and surface_ objects)
	clean_grid();

    // calculate the geometric center of the pose
	int n_atoms = atom_coords.size(); // len of vector
	core::Real total_x = 0.0;
	core::Real total_y = 0.0;
	core::Real total_z = 0.0;
    for (const auto& coord : atom_coords) {
        assert(coord.size() == 3 && "Inner vector does not have exactly 3 elements");
        total_x += coord[0];
        total_y += coord[1];
        total_z += coord[2];
    }
    std::tuple<core::Real, core::Real, core::Real> geometric_center(total_x / n_atoms, total_y / n_atoms,
                                                                    total_z / n_atoms);

    // 1. calculate vectors from all atom_coordinates (x, y, z) to the geometric center
    // 2. find max/min distance in x, y and z direction
    std::vector<std::tuple<core::Real, core::Real, core::Real>> atom_coords_wrt_geometric_center;
    core::Real max = 0;
    core::Real x_in_angstrom;
    core::Real y_in_angstrom;
    core::Real z_in_angstrom;
    for (const auto& atom_coord : atom_coords) {

        core::Real x = atom_coord[0];
        core::Real y = atom_coord[1];
        core::Real z = atom_coord[2];

        // alternative is valarray
        x_in_angstrom = x - std::get<0>(geometric_center);
        y_in_angstrom = y - std::get<1>(geometric_center);
        z_in_angstrom = z - std::get<2>(geometric_center);
        atom_coords_wrt_geometric_center.emplace_back(x_in_angstrom, y_in_angstrom, z_in_angstrom);

        max = std::max(std::abs(x_in_angstrom), max);
        max = std::max(std::abs(y_in_angstrom), max);
        max = std::max(std::abs(z_in_angstrom), max);
    }

    // find max grid_size_in_angstrom - this is the size of the grid if it was grid_size_in_voxels and not if it was actual_grid_size_in_voxels.
    core::Real half_grid_size_in_angstrom = max;

    // set the maximum atom radius depending on the pose residue types.
    core::Real max_atom_radius = *std::max_element(ljs.begin(), ljs.end());

    // increase the scaling to accomodate extra voxelization from taking the surface_type into account.
    half_grid_size_in_angstrom += max_atom_radius;
    if (surface_type_ != 3)
        half_grid_size_in_angstrom += probe_radius_;

    // if the user want to scale the voxelization, then this is set here.
    if (scale_on_) {
        half_grid_size_in_angstrom /= scaling_;
    }

    // setting the grid_resolution (see @details for more info)
    grid_res_ = (half_grid_size_in_angstrom / ((grid_dimension_ / 2) - 1 - geometric_center_at_grid_center_)) * padding;

    // iterate through
    std::vector<core::Real>::iterator lj_it;
    std::vector<std::tuple<core::Real, core::Real, core::Real>>::iterator atom_coord_wrt;
    std::vector<int>::iterator number_of_identical_neighbours;
    std::vector<Neighbour>::iterator neighbour;
    core::pose::Pose::iterator res;
    core::Real max_atom_acting_distance;
    int middle_grid_voxel = actual_grid_dimension_ / 2;
    int x_in_voxel, y_in_voxel, z_in_voxel, i;
    for (atom_coord_wrt = atom_coords_wrt_geometric_center.begin(), lj_it = ljs.begin(); lj_it != ljs.end(); atom_coord_wrt++, lj_it++) {
        // 1 put atom in the middle of the grid in euclid space (middle_grid_coord).
        // 2 put the atom at the place it was with with respect to the geometric center (which would correspond to putting the geometric center at middle_grid_coord).
        // 3 scale with the grid resolution
        // ceil has the same effect as floor in a 1 indexed space (which the grid is!)
        x_in_voxel = std::ceil( middle_grid_voxel + std::get<0>(*atom_coord_wrt) / grid_res_);
        y_in_voxel = std::ceil( middle_grid_voxel + std::get<1>(*atom_coord_wrt) / grid_res_);
        z_in_voxel = std::ceil( middle_grid_voxel + std::get<2>(*atom_coord_wrt) / grid_res_);

        grid_(x_in_voxel, y_in_voxel, z_in_voxel).fill();

        if (surface_type_ == 1 || surface_type_ == 2)
            max_atom_acting_distance = (*lj_it + probe_radius_) / grid_res_;
        else
            max_atom_acting_distance = *lj_it / grid_res_;

        // fill the surrounding neighbours
        for (number_of_identical_neighbours = identical_neighbours_.begin(), neighbour = neighbours_.begin();
             number_of_identical_neighbours != identical_neighbours_.end(); ++number_of_identical_neighbours) {
            core::Real distance_to_neighbour = neighbour->equal_grid_distance();
            if (distance_to_neighbour <= max_atom_acting_distance) {
                for (i = 0; i != *number_of_identical_neighbours; ++neighbour, ++i) {
                    int x_in_voxel_plus_neighbour = x_in_voxel + neighbour->x();
                    int y_in_voxel_plus_neighbour = y_in_voxel + neighbour->y();
                    int z_in_voxel_plus_neighbour = z_in_voxel + neighbour->z();
                    if(is_xyz_at_grid_edge(x_in_voxel_plus_neighbour, y_in_voxel_plus_neighbour, z_in_voxel_plus_neighbour)) {
                        std::string error = "(" + std::to_string(x_in_voxel_plus_neighbour) + "," + std::to_string(y_in_voxel_plus_neighbour) + ","
                                + std::to_string(z_in_voxel_plus_neighbour) + ") is filled at boundary and shouldn't be!";
                        if (!recalculate_){
                            utility_exit_with_message(error);
                        }
                        TR << error << " - Recalculating!" << std::endl;
                        padding +=  0.1;
                        goto beginning;
                    }
                    if(is_xyz_out_of_grid(x_in_voxel_plus_neighbour, y_in_voxel_plus_neighbour, z_in_voxel_plus_neighbour)) {
                        std::string error = "(" + std::to_string(x_in_voxel_plus_neighbour) + "," + std::to_string(y_in_voxel_plus_neighbour) + ","
                                + std::to_string(z_in_voxel_plus_neighbour) + ") is filled out of the grid and shouldn't be!.";
                        if (!recalculate_){
                            utility_exit_with_message(error);
                        }
                        TR << error << " - Recalculating!" << std::endl;
                        goto beginning;
                        padding +=  0.1;
                    }
                    grid_(x_in_voxel_plus_neighbour, y_in_voxel_plus_neighbour, z_in_voxel_plus_neighbour).fill();
                }
            } else
                break;
        }
    }
}


void
VoxelGrid::voxelize(core::pose::Pose const & pose) {
    // find x, y, z and lj of every atom of the pose
    std::vector<std::vector<core::Real>> atom_coords;
    std::vector<core::Real> ljs;
    core::Real lj;
    for (core::pose::Pose::const_iterator res = pose.begin(); res != pose.end(); ++res) {
        for (core::Size atom_index = 1; atom_index <= res->natoms(); ++atom_index) {
            // skip if atom is Hydrogen or if residue is virtual.
            if (!res->atom_is_hydrogen(atom_index && res->aa() != core::chemical::aa_vrt)) {
                core::Real x = res->xyz(atom_index).x();
                core::Real y = res->xyz(atom_index).y();
                core::Real z = res->xyz(atom_index).z();
                // only way I know of to test for a centroid atom
                if (res->atom_type(atom_index).lj_radius() == 0.0)
                    lj = res->nbr_radius();
                else
                    lj = res->atom_type(atom_index).lj_radius();
                atom_coords.emplace_back(x);
                atom_coords.emplace_back(y);
                atom_coords.emplace_back(z);
                ljs.emplace_back(lj);
            }
        }
    }
	voxelize(atom_coords, ljs);
}
//
// void VoxelGrid::voxel_grid_from_file( std::string filename)
// {
// 	 using namespace basic::options;
//    using namespace basic::options::OptionKeys;
//
// 	std::string type = option[zernike_descriptor::zernike_transform_type].value();
//
//    clean_grid();
//    grid_res_ = 2.0;
//     std::ifstream infile(filename.c_str());
//     std::string line;
//     while ( getline(infile, line) ) {
//       utility::vector1< std::string > tokens ( utility::split( line ) );
//       int x = std::stoi(tokens[1]);
//       int y = std::stoi(tokens[2]);
//       int z = std::stoi(tokens[3]);
// 			// put 2D slice at z = 0
// 			if ( type == "2D"  || type == "2D_3D" ) {
// //					TR << "Grid: " << x << " " << y << std::endl;
// 					for (int i = -5; i <= 5; i++ ) {
// //						z = int(grid_dimension_/2) + i;
// 						z = int(actual_grid_dimension_/2) + i;
// 						grid_(x,y,z).fill();
// 					}
//     	} else {
// 				grid_(x,y,z).fill();
// 			}
// 		}
// //   } else std::cout << "Unable to open file" << filename << std::endl;
//
//   return;
// }

void VoxelGrid::set_voxel_grid( std::vector< std::vector<int> > slice_for_ZT )
{
	 using namespace basic::options;
   using namespace basic::options::OptionKeys;

	std::string type = option[zernike_descriptor::zernike_transform_type].value();

   clean_grid();
   grid_res_ = 2.0;
   for (core::Size x=0; x< slice_for_ZT[0].size(); x++ ) {
    for (core::Size y=0; y< slice_for_ZT[0].size(); y++ ) {
			if ( slice_for_ZT[x][y] == 1 ) {
				for (int i = -5; i <= 5; i++ ) {
					int z = int(actual_grid_dimension_/2) + i;
					grid_(x,y,z).fill();
				}
			}
		}
	}
	return;
}

std::vector< std::vector<int> >
VoxelGrid::slice_outline_grid_from_file( std::string filename)
{
	 using namespace basic::options;
   using namespace basic::options::OptionKeys;

	std::string type = option[zernike_descriptor::zernike_transform_type].value();
	int grid_size = option[zernike_descriptor::grid_size].value();

	  std::vector< std::vector<int> > slice_for_ZT;
  for (int i = 0; i< grid_size; i++ ) {
    slice_for_ZT.push_back( std::vector<int>(grid_size) );
  }
    std::ifstream infile(filename.c_str());
    std::string line;
    while ( getline(infile, line) ) {
      utility::vector1< std::string > tokens ( utility::split( line ) );
      int x = std::stoi(tokens[1]);
      int y = std::stoi(tokens[2]);
//      int z = std::stoi(tokens[3]);
			slice_for_ZT[x][y] = 1;
		}
//   } else std::cout << "Unable to open file" << filename << std::endl;

  return slice_for_ZT;
}


bool
VoxelGrid::is_xyz_at_grid_edge(int x, int y, int z) {
    return x == 1 || x == actual_grid_dimension_ ||
           y == 1 || y == actual_grid_dimension_ ||
           z == 1 || z == actual_grid_dimension_;
}

bool
VoxelGrid::is_xyz_out_of_grid(int x, int y, int z) {
    return x < 1 || x > actual_grid_dimension_ ||
           y < 1 || y > actual_grid_dimension_ ||
           z < 1 || z > actual_grid_dimension_;
}

bool
VoxelGrid::is_utilizing_the_whole_grid() {
    for (int i = 2; i < actual_grid_dimension_ ; i++) {
        for (int k = 2; k < actual_grid_dimension_; k++) {
            for (int j = 2; j < actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_filled() && (i == 2 || i == actual_grid_dimension_ - 1 ||
                                                k == 2 || k == actual_grid_dimension_ - 1 ||
                                                j == 2 || j == actual_grid_dimension_ - 1)) {
                    return true;
                }
            }
        }
    }
    return false;
}
void
VoxelGrid::generate_neighbours() {

    core::Real x_shift;
    core::Real y_shift;
    core::Real z_shift;
    for (int x = -amount_of_voxel_neighbours_; x != amount_of_voxel_neighbours_ + 1; x++) {
        if (x < 0) {x_shift = 0.5;} else if (x==0) {x_shift = 0;} else {x_shift = -0.5;}
        for (int y = -amount_of_voxel_neighbours_; y != amount_of_voxel_neighbours_ + 1; y++) {
            if (y < 0) {y_shift = 0.5;} else if (y==0) {y_shift = 0;} else {y_shift = -0.5;}
            for (int z = -amount_of_voxel_neighbours_; z != amount_of_voxel_neighbours_ + 1; z++) {
                if (z < 0) {z_shift = 0.5;} else if (z==0) {z_shift = 0;} else {z_shift = -0.5;}
                if (x == 0 && y==0 && z==0) continue;
                core::Real equal_grid_distance = std::sqrt(((x+x_shift) * (x+x_shift) + (y+y_shift) * (y+y_shift) + (z+z_shift) * (z+z_shift)));
                neighbours_.emplace_back(Neighbour(x, y, z, equal_grid_distance));
            }
        }
    }

    // sort neighbours into ascending order.
    std::sort(neighbours_.begin(), neighbours_.end(), [](Neighbour a, Neighbour b) -> bool {
        return a.equal_grid_distance() < b.equal_grid_distance();
    });

    // find the number of neighbours that have similar values.
    int count = 1;
    // maybe assert that the vector is not empty which could be the case and then you get a bad acces error.
    for (auto it = neighbours_.begin() + 1; it != neighbours_.end(); ++it) {
        // there must be a better comparison method than this? // you can use max(of a or b) on the right as well..
        // could aso use the permutations of x,y,z. 1,0,1 is the same as 1,1,0 etc.
        if (std::fabs((it - 1)->equal_grid_distance() - it->equal_grid_distance())  < 0.0001 * std::fabs((it - 1)->equal_grid_distance())) {
            count += 1;
        }
        else {
            identical_neighbours_.push_back((count));
            count = 1;
        }
    }
}

bool
VoxelGrid::sanity_checks() {
    bool checks = true;
    if (!grid_not_filled_at_boundary())
        checks = false;
    return checks;
}

bool
VoxelGrid::grid_not_filled_at_boundary() {
    // for fixed x
    bool check = true;
    for (int i : {1, actual_grid_dimension_} ) {
        for (int k = 1; k <= actual_grid_dimension_; k++) {
            for (int j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_filled()) {
                    TR << "(" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ") is filled." << std::endl;
                    check = false;
                }
            }
        }
    }
    // for fixed y
    for (int k : {1, actual_grid_dimension_} ) {
        for (int i = 1; i <= actual_grid_dimension_; i++) {
            for (int j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_filled()) {
                    TR << "(" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ") is filled."
                       << std::endl;
                    check = false;
                }
            }
        }
    }
    // for fixed z
    for (int j : {1, actual_grid_dimension_} ) {
        for (int i = 1; i <= actual_grid_dimension_; i++) {
            for (int k = 1; k <= actual_grid_dimension_; k++) {
                if (grid_(i, k, j).is_filled()) {
                    TR << "(" + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ") is filled."
                       << std::endl;
                    check = false;
                }
            }
        }
    }
    if (check)
        TR << "grid_not_filled_at_boundary: OK" << std::endl;
    else
        TR << "grid_not_filled_at_boundary: NOT OK" << std::endl;
    return check;
}

void VoxelGrid::mark_surface(int x, int y, int z) {
    surface_voxels_frontier_.push(& grid_(x, y, z));

    // only import for debugging:
    grid_(x, y, z).surface_by_fs();
}

void VoxelGrid::find_surface_with_stack(){
    std::stack<std::vector<int>> stack;

    // start from the corner of the grid
    stack.push({1, 1, 1});

    // search
    while (!stack.empty()) {
        std::vector<int> index = stack.top();
        int x = index[0];
        int y = index[1];
        int z = index[2];
        stack.pop();

        // mark voxel as visited
        grid_(x, y, z).visit_by_fs();

        for(int i = 0; i < 26; ++i){
            // retrieve a neighbour
            int x_to_check = x + neighbours_[i].x();
            int y_to_check = y + neighbours_[i].y();
            int z_to_check = z + neighbours_[i].z();

            // if the neighbour is out of the grid; continue
            if (is_xyz_out_of_grid(x_to_check, y_to_check, z_to_check))
                continue;

            // if the nieghbour is filled, mark it it as surface if it is within the adjacent voxels
            if (grid_(x_to_check, y_to_check, z_to_check).is_filled() && i < 6)
                    mark_surface(x_to_check, y_to_check, z_to_check);

            // add adjacent neighbour voxels to the stack if the following conditions DO NOT apply:
            //  1. Voxel is already visited
            //  2. Voxel is filled.
            if ( !( grid_(x_to_check, y_to_check, z_to_check).is_visited_by_fs() ||
                    grid_(x_to_check, y_to_check, z_to_check).is_filled() )) {
                stack.push({x_to_check, y_to_check, z_to_check});
            }
        }
    }
}

void
VoxelGrid::find_surface_with_recursion(){
    std::tuple<int, int, int> edge;
    try {
        edge = find_edge();
    } catch(const char * err) {
        std::cerr << err << std::endl;
    }
    int x = std::get<0>(edge);
    int y = std::get<1>(edge);
    int z = std::get<2>(edge);
    fillsurface_recursion(x, y, z);
}

void
VoxelGrid::fillsurface_recursion(int x, int y, int z){

	// check if out of the grid:
	if (is_xyz_out_of_grid(x,y,z))
		return;

	// return if voxel is visited or filled
	if (grid_(x, y, z).is_visited_by_fs() || grid_(x, y, z).is_filled())
		return;

	// mark voxel as visited
	grid_(x, y, z).visit_by_fs();

	if (detect_neighbours(x, y, z)) {
		fillsurface_recursion(x, y, z + 1);
		fillsurface_recursion(x, y, z - 1);
		fillsurface_recursion(x + 1, y, z);
		fillsurface_recursion(x - 1, y, z);
		fillsurface_recursion(x, y + 1, z);
		fillsurface_recursion(x, y - 1, z);
	}

}

void
VoxelGrid::find_surface() {
    if (never_leave_neighbour_) {
        find_surface_with_recursion();
    } else {
        find_surface_with_stack();
    }
}

bool
VoxelGrid::search_in_z_dimension(int xi, int yi, int & zi){
    for (zi = 1; zi <= actual_grid_dimension_; ++zi) {
        if (grid_(xi, yi, zi).is_filled()) {
            return true;
        }
    }
    return false;
}

std::tuple<int, int, int>
VoxelGrid::find_edge() {
    int mid = actual_grid_dimension_ / 2;
    int zi = 0;

    // 1st quadrant
    for (int xi = mid; xi < actual_grid_dimension_; xi++) {
        for (int yi = mid; yi < actual_grid_dimension_; yi++) {
            if (search_in_z_dimension(xi, yi, zi))
                return std::make_tuple(xi, yi, zi - 1);
        }
    }

    // 2nd quadrant
    for (int xi = mid - 1; xi > 1; xi--) {
        for (int yi = mid; yi < actual_grid_dimension_; yi++) {
            if (search_in_z_dimension(xi, yi, zi))
                return std::make_tuple(xi, yi, zi - 1);
        }
    }

    // 3rd quadrant
    for (int xi = mid; xi < actual_grid_dimension_; xi++) {
        for (int yi = mid - 1; yi > 1; yi--) {
            if (search_in_z_dimension(xi, yi, zi))
                return std::make_tuple(xi, yi, zi - 1);
        }
    }

    // 4rd quadrant
    for (int xi = mid - 1; xi > 1; xi--) {
        for (int yi = mid - 1; yi > 1; yi--) {
            if (search_in_z_dimension(xi, yi, zi))
                return std::make_tuple(xi, yi, zi - 1);
        }
    }

     utility_exit_with_message("Did not find an edge!");

    //todo: stil?
    // compiler wants me to return a value outside of the for loops
//    return std::make_tuple(0,0,0);
}

// we go through all 26
bool
VoxelGrid::detect_neighbours(int x, int y, int z) {

    bool found_filled_voxel = false;
    for(int i = 0; i < 26; ++i){
        // retrieve a neighbour
        int x_to_check = x + neighbours_[i].x();
        int y_to_check = y + neighbours_[i].y();
        int z_to_check = z + neighbours_[i].z();
        // check if neighbour is out of the grid:
        if (is_xyz_out_of_grid(x_to_check, y_to_check, z_to_check))
            continue;
        // check if neighbour is filled
        if (grid_(x_to_check, y_to_check, z_to_check).is_filled()) {
            found_filled_voxel = true;
            // if neighbour is filled and within the 6 closets voxels then mark it as surface
            if (i < 6) {
                mark_surface(x_to_check, y_to_check, z_to_check);
            }
        }
    }
    return found_filled_voxel;
}

int
VoxelGrid::actual_grid_size() {
    return actual_grid_dimension_;
}

/// todo: in the future you can change it to py::array_t<double> to output a numpy array
/// todo: also the indexing in the for loop should be to + 1 !!!!!
/// todo: make it visible with in pyrosetta - somehow this function is not????
//std::map<std::string, pybind11::array_t<core::Real> >
//pybind11::array_t<int>
//VoxelGrid::get_grid_map() {
//    // collect data in vectors first!
////    std::map<std::string, pybind11::array_t<core::Real> > map_to_data;
//    std::vector<std::tuple<int,int,int>> filldata;
////    std::vector<std::tuple<core::Real,core::Real, core::Real>> surfacedata;
////    std::vector<std::tuple<core::Real,core::Real, core::Real>> visitdata;
//    int i, j, k;
//    for (i = 1; i <=  actual_grid_size_in_voxels_; i++) {
//        for (k = 1; k <= actual_grid_size_in_voxels_; k++) {
//            for (j = 1; j <= actual_grid_size_in_voxels_; j++) {
//                if (grid_(i, k, j).is_filled())
//                    filldata.emplace_back(i,j,k);
////                else if (grid_(i, k, j).visited_by_fs())
////                    visitdata.emplace_back(i,j,k);
////                else if (grid_(i, k, j).is_surface_by_fs())
////                    surfacedata.emplace_back(i,j,k);
//            }
//        }
//    }
//
//    int fillsize = 3*filldata.size();
//    int *fill = new int[fillsize];
//    for (std::tuple<int,int,int> xyz : filldata)
//        fill[i] = std::get<0>(xyz);
//    for (std::tuple<int,int,int> xyz : filldata)
//        fill[i] = std::get<1>(xyz);
//    for (std::tuple<int,int,int> xyz : filldata)
//        fill[i] = std::get<2>(xyz);
//
////    double** surfacedata_array =  new double**[surfacedata.size()][3];
////    double filldata_array[3][filldata.size()];
////    double filldata_array[3][filldata.size()];
//
////    map_to_data["filldata"] = filldata_array;
////    map_to_data["visitdata"] = visitdata_array;
////    map_to_data["surfacedata"] = surfacedata_array;
//
//    return pybind11::array_t<int>(
//                    {3, fillsize}, // shape
////                    {1000*1000*8, 1000*8, 4}, // C-style contiguous strides for double
//                    fill); // the data pointer
//}


/// Does not show up in pyrosetta - I think it is a problem with the std:tuple
//std::map<std::string,  std::vector<std::tuple<int, int, int>>>
//VoxelGrid::get_grid_map() {
//    std::map<std::string, std::vector<std::tuple<int, int, int>>> map_to_data;
//    std::vector<std::tuple<int, int, int>> filldata;
//    std::vector<std::tuple<int, int, int>> surface_by_fs_data;
//    std::vector<std::tuple<int, int, int>> visit_by_fs_data;
//    std::vector<std::tuple<int, int, int>> visit_by_edt_data;
//    std::vector<std::tuple<int, int, int>> surface_by_edt_data;
//    int i, j, k;
//    for (i = 1; i <= actual_grid_size_in_voxels_; i++) {
//        for (k = 1; k <= actual_grid_size_in_voxels_; k++) {
//            for (j = 1; j <= actual_grid_size_in_voxels_; j++) {
//                if (grid_(i, k, j).is_filled())
//                    filldata.emplace_back(std::tuple<int, int, int>(i, j, k));
//                else if (grid_(i, k, j).visited_by_fs())
//                    visit_by_fs_data.emplace_back(std::tuple<int, int, int>(i, j, k));
//                else if (grid_(i, k, j).is_surface_by_fs())
//                    surface_by_fs_data.emplace_back(std::tuple<int, int, int>(i, j, k));
//                else if (grid_(i, k, j).is_visited_by_edt())
//                    visit_by_edt_data.emplace_back(std::tuple<int, int, int>(i, j, k));
//                else if (grid_(i, k, j).is_surface_by_edt())
//                    surface_by_edt_data.emplace_back(std::tuple<int, int, int>(i, j, k));
//            }
//        }
//    }
//    map_to_data["fill"] = filldata;
//    map_to_data["visit_by_fs"] = visit_by_fs_data;
//    map_to_data["surface_by_fs"] = surface_by_fs_data;
//    map_to_data["surface_by_edt"] = surface_by_edt_data;
//    map_to_data["visit_by_edt"] = visit_by_edt_data;
//    return map_to_data;
//}

//std::map<std::string,  std::vector<std::vector<int>>>
//VoxelGrid::get_grid_mapping() {
//    std::map<std::string, std::vector<std::vector<int>>> map_to_data;
//
//    std::vector<std::vector<int>> filldata;
//    std::vector<int> filldata_x;
//    std::vector<int> filldata_y;
//    std::vector<int> filldata_z;
//
//    std::vector<std::vector<int>> surface_by_fs_data;
//    std::vector<int> surface_by_fs_x;
//    std::vector<int> surface_by_fs_y;
//    std::vector<int> surface_by_fs_z;
//
//    std::vector<std::vector<int>> visit_by_fs_data;
//    std::vector<int> visit_by_fs_x;
//    std::vector<int> visit_by_fs_y;
//    std::vector<int> visit_by_fs_z;
//
//    std::vector<std::vector<int>> visit_by_edt_data;
//    std::vector<int> visit_by_edt_x;
//    std::vector<int> visit_by_edt_y;
//    std::vector<int> visit_by_edt_z;
//
//    std::vector<std::vector<int>> surface_by_edt_data;
//    std::vector<int> surface_by_edt_x;
//    std::vector<int> surface_by_edt_y;
//    std::vector<int> surface_by_edt_z;
//    int i, j, k;
//    for (i = 1; i <= actual_grid_size_in_voxels_; i++) {
//        for (k = 1; k <= actual_grid_size_in_voxels_; k++) {
//            for (j = 1; j <= actual_grid_size_in_voxels_; j++) {
//                if (grid_(i, k, j).is_filled()) {
//                    filldata_x.emplace_back(i);
//                    filldata_y.emplace_back(j);
//                    filldata_z.emplace_back(k);
//                }
//                else if (grid_(i, k, j).visited_by_fs()) {
//                    visit_by_fs_x.emplace_back(i);
//                    visit_by_fs_y.emplace_back(j);
//                    visit_by_fs_z.emplace_back(k);
//                }
//                else if (grid_(i, k, j).is_surface_by_fs()) {
//                    surface_by_fs_x.emplace_back(i);
//                    surface_by_fs_y.emplace_back(j);
//                    surface_by_fs_z.emplace_back(k);
//                }
//                else if (grid_(i, k, j).is_visited_by_edt()) {
//                    visit_by_edt_x.emplace_back(i);
//                    visit_by_edt_y.emplace_back(j);
//                    visit_by_edt_z.emplace_back(k);
//                }
//                else if (grid_(i, k, j).is_surface_by_edt()) {
//                    surface_by_edt_x.emplace_back(i);
//                    surface_by_edt_y.emplace_back(j);
//                    surface_by_edt_z.emplace_back(k);
//                }
//            }
//        }
//    }
//
//    filldata.emplace_back(filldata_x);
//    filldata.emplace_back(filldata_y);
//    filldata.emplace_back(filldata_z);
//
//    visit_by_fs_data.emplace_back(visit_by_fs_x);
//    visit_by_fs_data.emplace_back(visit_by_fs_y);
//    visit_by_fs_data.emplace_back(visit_by_fs_z);
//
//    surface_by_fs_data.emplace_back(surface_by_fs_x);
//    surface_by_fs_data.emplace_back(surface_by_fs_y);
//    surface_by_fs_data.emplace_back(surface_by_fs_z);
//
//    visit_by_edt_data.emplace_back(visit_by_edt_x);
//    visit_by_edt_data.emplace_back(visit_by_edt_y);
//    visit_by_edt_data.emplace_back(visit_by_edt_z);
//
//    surface_by_edt_data.emplace_back(surface_by_edt_x);
//    surface_by_edt_data.emplace_back(surface_by_edt_y);
//    surface_by_edt_data.emplace_back(surface_by_edt_z);
//
//    map_to_data["fill"] = filldata;
//    map_to_data["visit_by_fs"] = visit_by_fs_data;
//    map_to_data["surface_by_fs"] = surface_by_fs_data;
//    map_to_data["visit_by_edt"] = visit_by_edt_data;
//    map_to_data["surface_by_edt"] = surface_by_edt_data;
//
//    return map_to_data;
//}


std::vector<std::vector<int>>
VoxelGrid::get_filled_voxels() {
    std::vector<std::vector<int>> filldata;
    std::vector<int> filldata_x;
    std::vector<int> filldata_y;
    std::vector<int> filldata_z;
    int i, j, k;
    for (i = 1; i <= actual_grid_dimension_; i++) {
        for (k = 1; k <= actual_grid_dimension_; k++) {
            for (j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_filled()) {
                    filldata_x.emplace_back(i);
                    filldata_y.emplace_back(j);
                    filldata_z.emplace_back(k);
                }
            }
        }
    }
    filldata.emplace_back(filldata_x);
    filldata.emplace_back(filldata_y);
    filldata.emplace_back(filldata_z);
    return filldata;
}

std::vector<std::vector<int>>
VoxelGrid::get_visit_by_fs_voxels() {
    std::vector<std::vector<int>> visit_by_fs;
    std::vector<int> visit_by_fs_x;
    std::vector<int> visit_by_fs_y;
    std::vector<int> visit_by_fs_z;
    int i, j, k;
    for (i = 1; i <= actual_grid_dimension_; i++) {
        for (k = 1; k <= actual_grid_dimension_; k++) {
            for (j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_visited_by_fs()) {
                    visit_by_fs_x.emplace_back(i);
                    visit_by_fs_y.emplace_back(j);
                    visit_by_fs_z.emplace_back(k);
                }
            }
        }
    }
    visit_by_fs.emplace_back(visit_by_fs_x);
    visit_by_fs.emplace_back(visit_by_fs_y);
    visit_by_fs.emplace_back(visit_by_fs_z);
    return visit_by_fs;
}

std::vector<std::vector<int>>
VoxelGrid::get_surface_by_fs_voxels() {
    std::vector<std::vector<int>> surface_by_fs;
    std::vector<int> surface_by_fs_x;
    std::vector<int> surface_by_fs_y;
    std::vector<int> surface_by_fs_z;
    int i, j, k;
    for (i = 1; i <= actual_grid_dimension_; i++) {
        for (k = 1; k <= actual_grid_dimension_; k++) {
            for (j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_surface_by_fs()) {
                    surface_by_fs_x.emplace_back(i);
                    surface_by_fs_y.emplace_back(j);
                    surface_by_fs_z.emplace_back(k);
                }
            }
        }
    }
    surface_by_fs.emplace_back(surface_by_fs_x);
    surface_by_fs.emplace_back(surface_by_fs_y);
    surface_by_fs.emplace_back(surface_by_fs_z);
    return surface_by_fs;
}

std::vector<std::vector<int>>
VoxelGrid::get_visit_by_edt_voxels() {
    std::vector<std::vector<int>> visit_by_edt;
    std::vector<int> visit_by_edt_x;
    std::vector<int> visit_by_edt_y;
    std::vector<int> visit_by_edt_z;
    int i, j, k;
    for (i = 1; i <= actual_grid_dimension_; i++) {
        for (k = 1; k <= actual_grid_dimension_; k++) {
            for (j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_visited_by_edt()) {
                    visit_by_edt_x.emplace_back(i);
                    visit_by_edt_y.emplace_back(j);
                    visit_by_edt_z.emplace_back(k);
                }
            }
        }
    }
    visit_by_edt.emplace_back(visit_by_edt_x);
    visit_by_edt.emplace_back(visit_by_edt_y);
    visit_by_edt.emplace_back(visit_by_edt_z);
    return visit_by_edt;
}

std::vector<std::vector<int>>
VoxelGrid::get_surface_by_edt_voxels() {
    std::vector<std::vector<int>> surface_by_edt;
    std::vector<int> surface_by_edt_x;
    std::vector<int> surface_by_edt_y;
    std::vector<int> surface_by_edt_z;
    int i, j, k;
    for (i = 1; i <= actual_grid_dimension_; i++) {
        for (k = 1; k <= actual_grid_dimension_; k++) {
            for (j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_surface_by_edt()) {
                    surface_by_edt_x.emplace_back(i);
                    surface_by_edt_y.emplace_back(j);
                    surface_by_edt_z.emplace_back(k);
                }
            }
        }
    }
    surface_by_edt.emplace_back(surface_by_edt_x);
    surface_by_edt.emplace_back(surface_by_edt_y);
    surface_by_edt.emplace_back(surface_by_edt_z);
    return surface_by_edt;
}

void
VoxelGrid::write_json_data(std::ofstream & file, std::string voxeltype, std::vector<std::vector<int>> & voxelvector) {
    file << "\t\"" << voxeltype << "\": [";
    for (size_t i = 0; i < voxelvector.size(); ++i) {
        const auto& xyz = voxelvector[i];
        file << "[";
        for (size_t j = 0; j < xyz.size(); ++j) {
            file << xyz[j];
            if (j < xyz.size() - 1) {
                file << ", ";
            }
        }
        file << "]";
        if (i < voxelvector.size() - 1) {
            file << ", ";
        }
    }
    file << "]";
}

void
VoxelGrid::output_grid_json(std::string name) {
    // setup stream
    std::ofstream file;
    file.open(name);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << name << std::endl;
        return;
    }

    std::vector<std::vector<int>> filldata;
    std::vector<std::vector<int>> visited;
    std::vector<std::vector<int>> surface;
    std::vector<std::vector<int>> surface_edt;
    std::vector<std::vector<int>> visited_edt;

    int i, j, k;
    for (i = 1; i <= actual_grid_dimension_; i++) {
        for (k = 1; k <= actual_grid_dimension_; k++) {
            for (j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_filled())
                    filldata.emplace_back(std::vector<int>{i, j, k});
                if (grid_(i, k, j).is_visited_by_fs())
                    visited.emplace_back(std::vector<int>{i, j, k});
                if (grid_(i, k, j).is_surface_by_fs())
                    surface.emplace_back(std::vector<int>{i, j, k});
                if (grid_(i, k, j).is_surface_by_edt())
                    surface_edt.emplace_back(std::vector<int>{i, j, k});
                if (grid_(i, k, j).is_visited_by_edt())
                    visited_edt.emplace_back(std::vector<int>{i, j, k});
            }
        }
    }
    file << "{\n";
    write_json_data(file, "filled", filldata);
    file << ",\n";
    write_json_data(file, "visited", visited);
    file << ",\n";
    write_json_data(file, "surface", surface);
    file << ",\n";
    write_json_data(file, "visited_edt", visited_edt);
    file << ",\n";
    write_json_data(file, "surface_edt", surface_edt);
    file << "\n}\n";

    file.close();
}

// todo: 1 index these badboys instead
void
VoxelGrid::output_grid_csv(std::string name) {
    // write visit and fill data
    std::ofstream filldata;
    std::ofstream visit_by_fs_data;
    std::ofstream surfacedata;
    std::ofstream surface_by_edt_data;
    std::ofstream visit_by_edt_data;
    filldata.open(name + "filled.voxels");
    filldata << "i,k,j,boxlength" << std::endl;
    visit_by_fs_data.open(name + "visit_by_fs.voxels");
    visit_by_fs_data << "i,k,j,boxlength" << std::endl;
    surfacedata.open(name + "surface.voxels");
    surfacedata << "i,k,j,boxlength" << std::endl;
    surface_by_edt_data.open(name + "surface_by_edt.voxels");
    surface_by_edt_data << "i,k,j,boxlength" << std::endl;
    visit_by_edt_data.open(name + "visit_by_edt.voxels");
    visit_by_edt_data << "i,k,j,boxlength" << std::endl;
    int i, j, k;
    for (i = 1; i <= actual_grid_dimension_; i++) {
        for (k = 1; k <= actual_grid_dimension_; k++) {
            for (j = 1; j <= actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_filled())
                    filldata << i << "," << k << "," << j << "," << actual_grid_dimension_ << "\n";
                if (grid_(i, k, j).is_visited_by_fs())
                    visit_by_fs_data << i << "," << k << "," << j << "," << actual_grid_dimension_ << "\n";
                if (grid_(i, k, j).is_surface_by_fs())
                    surfacedata << i << "," << k << "," << j << "," << actual_grid_dimension_ << "\n";
                if (grid_(i, k, j).is_surface_by_edt())
                    surface_by_edt_data << i << "," << k << "," << j << "," << actual_grid_dimension_ << "\n";
                if (grid_(i, k, j).is_visited_by_edt())
                    visit_by_edt_data << i << "," << k << "," << j << "," << actual_grid_dimension_ << "\n";
            }
        }
    }
    // write surface data
    filldata.close();
    visit_by_fs_data.close();
    surfacedata.close();
    visit_by_edt_data.close();
    surface_by_edt_data.close();
}

/// @brief: Writes the grid to disk according to the binvoxformat (run length encoding compression).
void VoxelGrid::output_grid_binvox(std::string name) {
    std::ofstream binvoxdata =  std::ofstream(name, std::ios::out | std::ios::binary);
    // write header
    const char *header1 = "#binvox 1\n";
    binvoxdata.write(header1, strlen(header1)); // char is always 1 byte. second argument tells how many bytes the file pointer ahead. strlen gives 1*length of array
    std::string temp("dim");
    for (int i = 0; i != 3; ++i) {
        temp += " ";
        temp += std::to_string(actual_grid_dimension_);
    }
    temp += "\n";
    const char *header2 = temp.c_str();
    binvoxdata.write(header2, strlen(header2));
    const char *header3 = "translate 0 0 0\n";
    binvoxdata.write(header3, strlen(header3));
    const char *header4 = "scale 1\n";
    binvoxdata.write(header4, strlen(header4));
    const char *header5 = "data\n";
    binvoxdata.write(header5, strlen(header5));

    // This is the voxel method according to http://www.patrickmin.com/binvox/binvox.html
    int i,j,k;
    std::vector<int> flat(actual_grid_dimension_ * actual_grid_dimension_ * actual_grid_dimension_);
    for (i = 0; i < actual_grid_dimension_; i++) {
        for (k = 0; k < actual_grid_dimension_; k++) {
            for (j = 0; j < actual_grid_dimension_; j++) {
                if (grid_(i, k, j).is_filled())
                    flat[i * actual_grid_dimension_ * actual_grid_dimension_ + j * actual_grid_dimension_ + k] = 1;
            }
        }
    }

    // write data
    int state = flat[0];
    int count = 0;
    for (int index = 0; index != actual_grid_dimension_ * actual_grid_dimension_ * actual_grid_dimension_; ++index) {
        if (flat[index] == state) {
            count += 1;
            if (count == 255) {
                binvoxdata.write(reinterpret_cast<char *>(&state), 1);
                binvoxdata.write(reinterpret_cast<char *>(&count), 1);
                count = 0;
            }
        } else {
            binvoxdata.write(reinterpret_cast<char *>(&state), 1);
            binvoxdata.write(reinterpret_cast<char *>(&count), 1);
            state = flat[index];
            count = 1;
        }
    }
    if (count > 0) {
        binvoxdata.write(reinterpret_cast<char *>(&state), 1);
        binvoxdata.write(reinterpret_cast<char *>(&count), 1);}

    binvoxdata.close();
}

void VoxelGrid::clean_grid() {
    // before we begin we need to make sure that we have no filled/visited voxels in our grid
    // we do not need to clean the surface member variable since the nsurface variable keeps track on the # of surface voxels
    // but we do need to set it to 0.
    // todo: are we cleaning the last voxels?
    int x,y,z;
    for (x = 1; x <= actual_grid_dimension_; ++x) {
        for (y = 1; y <= actual_grid_dimension_; ++y) {
            for (z = 1; z <= actual_grid_dimension_; ++z)
                grid_(x, y, z).clear();
        }
    }
//   // Nuke arrays
    // Doesn't need to be nuked necessarily because the code goes through the queue until it is empty (see edt())
    surface_voxels_frontier_ = {};
    shell_ = {};
}

int
VoxelGrid::grid_size() {
    return grid_dimension_;
}

std::string
VoxelGrid::surface_type(){
    if (surface_type_ == 1)
        return "MS";
    else if (surface_type_ == 2)
        return "SAS";
    else // surface_type == 3
        return "VDW";
}

core::Real
VoxelGrid::grid_resolution() {
    return grid_res_;
}

core::Real
VoxelGrid::probe_radius() {
    return probe_radius_;
}

core::Real
VoxelGrid::shell_thickness() {
    return shell_thickness_;
}

void
VoxelGrid::set_surface_vectors( numeric::xyzVector<Real>& surface_vector1,
				numeric::xyzVector<Real>& surface_vector2 ){
	slice_vector1_ = surface_vector1;
	slice_vector2_ = surface_vector2;
}

numeric::xyzVector<Real>
VoxelGrid::slice_surface_vector1() {
    return slice_vector1_;
}

numeric::xyzVector<Real>
VoxelGrid::slice_surface_vector2() {
    return slice_vector2_;
}

core::Real
VoxelGrid::slice_average_x() {
	return slice_average_x_;
}

core::Real
VoxelGrid::slice_average_y() {
	return slice_average_y_;
}

void
VoxelGrid::set_probe_radius( core::Real radius ) {
    probe_radius_ = radius;
}

void
VoxelGrid::set_shell_thickness( core::Real thickness ) {
    shell_thickness_ = thickness;
}

// todo you should also consider the euclidian  space (aka grid_res)
// todo: incorporate the changes that ingemar made for the linearization.
void
VoxelGrid::edt() {

    set_edt_search_depth();

    Voxel *v_center;
    Voxel *v_neighbour;
    int lim  = std::numeric_limits<int>::max();
    int x, y, z;
    while (!surface_voxels_frontier_.empty()) {

        // Get one voxel from the edt frontier (starting with the surface voxels)
        // v_center is the center voxel, we consider the neighbours of it below
        v_center = surface_voxels_frontier_.front();
        surface_voxels_frontier_.pop();


        // First we look at v_center

        // mark as visit - THIS IS ONLY A DEBUG FEATURE AND DOES NOT DO ANYTHING TO THE CODE
        v_center->visit_by_edt();

        // if it is part of the shell then mark it as so and add it to shell voxels.
        if (!v_center->is_surface_by_edt() && surface_by_edt_criteria(0)) {
            v_center->surface_by_edt();
            shell_.push(v_center);
        }

        // Now we consider the neighbours of v_center (v_neighbour)

        // iterate over all  26 neighbours
        for (int i = 0; i < 26; ++i) {
            // retrieve the neighbour x,y,z position in the grid
            // todo: put this neighbour retrieving in the voxelgrid class?
            x = v_center->x() + neighbours_[i].x();
            y = v_center->y() + neighbours_[i].y();
            z = v_center->z() + neighbours_[i].z();

            // We consider the neighbour ONLY if it is:
            // - filled: because then that is part of the protein - the "inside".
            // - NOT visited_by_fs: because this marks what is not the protein - the "outside".
            // - NOT surface_by_fs: because these are already in the surface_voxels queue and are dealt with as v_center above.
            // - NOT surface_by_edt: because these are already put in the surface_voxels_by_edt queue
            //   and therefore does not need to be reconsidered.
            if ( grid_(x,y,z).is_filled() &&
                !grid_(x,y,z).is_visited_by_fs() &&
                !grid_(x,y,z).is_surface_by_fs() &&
                !grid_(x,y,z).is_surface_by_edt()) {
                // retrieve a neighbour from the grid
                v_neighbour = & grid_(x, y, z);
                // mark as visit - THIS IS ONLY A DEBUG FEATURE AND DOES NOT DO ANYTHING TO THE CODE
                v_neighbour->visit_by_edt();
                // calculate the distance in voxel space from the root coordinates of v_center to v_neighbour.
                core::Real dist_in_voxel_space = v_center->edt_root_distance(*v_neighbour);
                // dont do more than nescesarry => stop when dist_in_voxel_space is more than max_edt_search
                // max_edt_search depends on the surface type and is set at the top of the edt() call.
                if ( dist_in_voxel_space > max_edt_search_)
                    continue;
                // get the current recorded distance (found in another iteration) or defaults to lim (big number)
                core::Real stored_dist_in_voxel_space = v_neighbour->edt_dist();
                // The following does:
                // true if the new distance found is closer than what was previously found
                // - Then reset the edt_dist (stored_dist_in_voxel_space) and reset the root xyz position for that voxel
                if (dist_in_voxel_space < stored_dist_in_voxel_space) {
                    v_neighbour->set_edt_dist(dist_in_voxel_space);
                    v_neighbour->edt_root_from_voxel(*v_center);
                    // if the voxel is never visited before add it to the surface_voxel queue
                    // true if newer visited before because stored_dist_in_angstrom is initialized as very big (stored_dist_in_angstrom == lim).
                    if (stored_dist_in_voxel_space == lim)
                        surface_voxels_frontier_.push(v_neighbour);
                }

                // if the v_neighbour should be part of the shell, then put in the shell.
                if (surface_by_edt_criteria(dist_in_voxel_space)) {
                    v_neighbour->surface_by_edt();
                    shell_.push(v_neighbour);
                }
            }
        }
    }
}

void
VoxelGrid::set_edt_search_depth() {
    if (surface_type_ == 1)
        max_edt_search_ = shell_thickness_ + std::ceil(probe_radius_ / grid_res_);
    else
        max_edt_search_ = shell_thickness_;
}

bool
VoxelGrid::surface_by_edt_criteria(core::Real dist_in_voxel_space) {
    if (surface_type_ == 1) {
        core::Real dist_in_angstrom = dist_in_voxel_space * grid_res_;
        return dist_in_angstrom > probe_radius_;
    }
    // surface_type 2 or 3
    return true;
}

void
VoxelGrid::zernike_transform(int order) {
    double* zernike_array = new double[grid_dimension_ * grid_dimension_ * grid_dimension_] {0};
    while (!shell_.empty()) {
        Voxel *v = shell_.front();
        shell_.pop();
        // linear index from ingemar andre commit 98f7f35 on madsjeppesen/shape.
        // int linear_index( ( (z-1)*gridsize  + (y-1) )*gridsize + x-1);
        // -3 because the actual grid size is dim+2 (see the VoxelGrid description).
        int linear_index = ((v->z() - 2) * grid_dimension_ + (v->y() - 2)) * grid_dimension_ + v->x() - 2;
        zernike_array[linear_index] = 1;
    }
    zernike_descriptor_ = zernike::ZernikeDescriptor<double, double>(zernike_array, grid_dimension_, order);
    delete [] zernike_array;
}

void
VoxelGrid::zernike_transform_3D_slice(int order, std::vector< std::vector< std::vector<int> > > slice) {

	using namespace basic::options;
   	using namespace basic::options::OptionKeys;

        bool debug = option[zernike_descriptor::debug].value();

	if (debug) {
   		std::ofstream savefile32dd("pose.slice.grid.2d.3d.dat");
	}

	double* zernike_array = new double[grid_dimension_ * grid_dimension_ * grid_dimension_] {0};
		std::vector< Voxel > slice_voxel;
		int slice_dim( slice.size());
		for (int i = 0; i< slice_dim; i++ ) {
			for (int j = 0; j< slice_dim; j++ ) {
				std::vector<int> s(slice[i][j]);
				slice_voxel.push_back( Voxel(s[0], s[1], s[2]) );
			}
		}

	int num_surface_voxels = 0;
	int max_dim = 0;
    while (!shell_.empty()) {
        Voxel *v = shell_.front();
        shell_.pop();
        // linear index from ingemar andre commit 98f7f35 on madsjeppesen/shape.
        // int linear_index( ( (z-1)*gridsize  + (y-1) )*gridsize + x-1);
        // -3 because the actual grid size is dim+2 (see the VoxelGrid description).
				std::vector<Voxel>::const_iterator b = slice_voxel.begin();
				std::vector<Voxel>::const_iterator e = slice_voxel.end();
				if (std::find(b, e, *v) != e) {
					if (debug) {
    						std::ofstream savefile32dd("pose.slice.grid.2d.3d.dat", std::ofstream::app);
						savefile32dd << v->x() << " " << v->y() << std::endl;
					}
					num_surface_voxels += 1;
					if (v->x() > max_dim) max_dim = v->x();
					if (v->y() > max_dim) max_dim = v->y();

        	int linear_index = ((v->z() - 2) * grid_dimension_ + (v->y() - 2)) * grid_dimension_ + v->x() - 2;
        	zernike_array[linear_index] = 1;
				}
    }
		num_surface_2d_pixels_ = num_surface_voxels;
		max_surface_dimension_ = max_dim;

    zernike_descriptor_ = zernike::ZernikeDescriptor<double, double>(zernike_array, grid_dimension_, order);
    delete [] zernike_array;
}


/*
void
VoxelGrid::sample_random_surface_vectors()
{
	numeric::xyzVector<Real> surface_vector1(numeric::random::rg().uniform(),numeric::random::rg().uniform(),numeric::random::rg().uniform());
	numeric::xyzVector<Real> surface_vector2(numeric::random::rg().uniform(),numeric::random::rg().uniform(),numeric::random::rg().uniform());
	surface_vector1.normalize();
	surface_vector2.normalize();
	numeric::xyzVector<Real> perp (cross(surface_vector1, surface_vector2 ) );
	surface_vector2 = cross(  perp, surface_vector1 );
	surface_vector2.normalize();
	slice_vector1_ = surface_vector1;
	slice_vector2_ = surface_vector2;
	return;
}
*/

// Copy paste programming. This function lives in protocols, so I cannot call it in core
numeric::xyzMatrix_real
VoxelGrid::random_reorientation_matrix(const core::Real phi_range, const core::Real psi_range)
{
	// a genuine rotation matrix which will randomly reorient the coord sys.
	// from Euler theorem
	const core::Real phi( phi_range * numeric::random::rg().uniform() ); // degrees
	const core::Real psi( psi_range * numeric::random::rg().uniform() ); // degrees
	const core::Real theta(
		numeric::conversions::degrees( std::acos(numeric::sin_cos_range( 1.0 - 2.0 * numeric::random::rg().uniform() ) ) )
	); // degrees

	//TR << "random_reorientation_matrix phi: " << phi << " psi: " << psi << " theta: " << theta << std::endl;
	return
		numeric::z_rotation_matrix_degrees(  psi   ) *
		numeric::y_rotation_matrix_degrees(  theta ) *
		numeric::z_rotation_matrix_degrees(  phi   );
}


void
VoxelGrid::sample_random_surface_vectors()
{
	numeric::xyzVector<Real> surface_vector1(1,0,0);
	numeric::xyzVector<Real> surface_vector2(0,1,0);
	numeric::xyzMatrix< Real > rot = random_reorientation_matrix();
	surface_vector1 = rot*surface_vector1;
	surface_vector2 = rot*surface_vector2;
	surface_vector1.normalize();
	surface_vector2.normalize();
	slice_vector1_ = surface_vector1;
	slice_vector2_ = surface_vector2;
	return;
}


core::pose::PoseOP
VoxelGrid::rotate_pose( core::pose::Pose const & pose) {
	core::pose::PoseOP rotated_pose = rotate_pose(pose, slice_vector1_, slice_vector2_);
	return rotated_pose;
}

core::pose::PoseOP
VoxelGrid::rotate_pose( core::pose::Pose const & pose, numeric::xyzVector<Real> surface_vector1, numeric::xyzVector<Real> surface_vector2 ) {
    numeric::xyzVector<Real> rot_x = surface_vector1.normalize();
    numeric::xyzVector<Real> rot_y = surface_vector2.normalize();
    numeric::xyzVector<Real> rot_z = cross(rot_x, rot_y ).normalize();
    numeric::xyzVector<Real> rot_x_new = cross(rot_y, rot_z ).normalize();
    //numeric::xyzVector<Real> rot_x_new = cross(rot_z, rot_y ).normalize();
    rot_x = rot_x_new;
    numeric::xyzMatrix<core::Real> rot_mat = numeric::xyzMatrix<core::Real>::rows(rot_x, rot_y, rot_z);
    //TR << "DET Voxelgrid " << rot_mat.det() << std::endl;
    core::pose::PoseOP pose_rotated(pose.clone());
    core::Size nres (pose_rotated->total_residue() );
    addVirtualResAsRoot(*pose_rotated);
    core::conformation::ResidueOP VRT_base = core::conformation::ResidueFactory::create_residue(*core::pose::get_restype_for_pose( *pose_rotated, "VRT" ));
    pose_rotated->append_residue_by_jump(*VRT_base, nres+1);
		kinematics::FoldTree newF( pose_rotated->fold_tree() );
    newF.reorder( nres+2 );
    pose_rotated->fold_tree( newF );
    numeric::xyzVector<core::Real> com(0,0,0); //core::pose::get_center_of_mass(*pose);
    core::kinematics::Jump j = pose_rotated->jump(2);
    j.set_translation(com);
    pose_rotated->set_jump(2, j);
    j.set_rotation(rot_mat);
    pose_rotated->set_jump(2, j);
    return pose_rotated;
}

std::vector< std::vector< std::vector<int> > >
VoxelGrid::xy_plane(core::Size N) {
	int max_len = N;
	//numeric::xyzVector<Real> center(N/2-1,N/2-1,N/2-1);
	numeric::xyzVector<Real> center(N/2,N/2,N/2);
	std::vector< std::vector< std::vector<int> > > grid;
	for (int i=0; i< max_len; i++) {
		grid.push_back(std::vector< std::vector<int> >(max_len));
	}
	numeric::xyzVector<Real> xvec(1,0,0);
	numeric::xyzVector<Real> yvec(0,1,0);

	for ( int dir1=-max_len/2; dir1 < max_len/2; dir1++) {
		for ( int dir2=-max_len/2; dir2 < max_len/2; dir2++) {
    	numeric::xyzVector<Real> v = center + xvec*dir1 + yvec*dir2;
			int x = round(v[0]);
			int y = round(v[1]);
			int z = round(v[2]);
        if ( x >= 0 && y >= 0 && z >= 0 &&
					 x < max_len && y < max_len && z < max_len ) {
						int grid_x = dir1 + max_len/2;
						int grid_y = dir2 + max_len/2;
						std::vector<int> grid_pos;
						grid_pos.push_back(x);
						grid_pos.push_back(y);
						grid_pos.push_back(z);
            grid[grid_x][grid_y] =grid_pos;
				}
		}
	}
	return grid;

}

std::vector< std::vector< std::vector<int> > >
VoxelGrid::sample_slice_plane( core::Size N) {
	numeric::xyzVector<Real> surface_vector1(numeric::random::rg().uniform(),numeric::random::rg().uniform(),numeric::random::rg().uniform());
	numeric::xyzVector<Real> surface_vector2(numeric::random::rg().uniform(),numeric::random::rg().uniform(),numeric::random::rg().uniform());
	surface_vector1.normalize();
	surface_vector2.normalize();
	numeric::xyzVector<Real> perp (cross(surface_vector1, surface_vector2 ) );
	surface_vector2 = cross(  perp, surface_vector1 );
  surface_vector2.normalize();
  slice_vector1_ = surface_vector1;
  slice_vector2_ = surface_vector2;
	//numeric::xyzVector<Real> center(N/2-1,N/2-1,N/2-1);
	numeric::xyzVector<Real> center(N/2,N/2,N/2);
	int max_len = int(N*sqrt(2));
	std::vector< std::vector< std::vector<int> > > grid;
	for (int i=0; i< max_len; i++) {
		grid.push_back(std::vector< std::vector<int> >(max_len));
	}
	for ( int dir1=-max_len/2; dir1 < max_len/2; dir1++) {
		for ( int dir2=-max_len/2; dir2 < max_len/2; dir2++) {
    		numeric::xyzVector<Real> v = center + surface_vector1*dir1 + surface_vector2*dir2;
				int x = round(v[0]);
				int y = round(v[1]);
				int z = round(v[2]);
        if ( x >= 0 && y >= 0 && z >= 0 &&
					 x < max_len && y < max_len && z < max_len ) {
						int grid_x = dir1 + max_len/2;
						int grid_y = dir2 + max_len/2;
						std::vector<int> grid_pos;
						grid_pos.push_back(x);
						grid_pos.push_back(y);
						grid_pos.push_back(z);
            grid[grid_x][grid_y] =grid_pos;
				}
		}
	}
	return grid;
}

std::vector< std::vector< std::vector<int> > >
VoxelGrid::slice_plane_from_surface_vectors( numeric::xyzVector<Real> surface_vector1, numeric::xyzVector<Real> surface_vector2, core::Size N) {
	surface_vector1.normalize();
	surface_vector2.normalize();
	numeric::xyzVector<Real> perp (cross(surface_vector1, surface_vector2 ) );
	surface_vector2 = cross(  perp, surface_vector1 );
  surface_vector2.normalize();

  slice_vector1_ = surface_vector1;
  slice_vector2_ = surface_vector2;
	//numeric::xyzVector<Real> center(N/2-1,N/2-1,N/2-1);
	numeric::xyzVector<Real> center(N/2,N/2,N/2);
	int max_len = int(N*sqrt(2));
	std::vector< std::vector< std::vector<int> > > grid;
	for (int i=0; i< max_len; i++) {
		grid.push_back(std::vector< std::vector<int> >(max_len));
	}
	for ( int dir1=-max_len/2; dir1 < max_len/2; dir1++) {
		for ( int dir2=-max_len/2; dir2 < max_len/2; dir2++) {
    		numeric::xyzVector<Real> v = center + surface_vector1*dir1 + surface_vector2*dir2;
				int x = round(v[0]);
				int y = round(v[1]);
				int z = round(v[2]);
        if ( x >= 0 && y >= 0 && z >= 0 &&
					 x < max_len && y < max_len && z < max_len ) {
						int grid_x = dir1 + max_len/2;
						int grid_y = dir2 + max_len/2;
						std::vector<int> grid_pos;
						grid_pos.push_back(x);
						grid_pos.push_back(y);
						grid_pos.push_back(z);
            grid[grid_x][grid_y] =grid_pos;
				}
		}
	}
	return grid;
}

std::vector< std::vector<int> >
VoxelGrid::shrink_2D_grid( std::vector< std::vector<int> > slice_for_ZT )
{
	int N( slice_for_ZT.size());
	int max_x = 0;
	int max_y = 0;
	int min_x = 10000000;
	int min_y = 10000000;
	float x_sum = 0;
	float y_sum = 0;
	Size num_voxels = 0;
	for (int i = 0; i< N; i++ ) {
		for (int j = 0; j< N; j++ ) {
			if  (slice_for_ZT[i][j] == 1) {
				x_sum += i;
				y_sum += j;
				num_voxels++;
				if (i > max_x)  max_x = i;
				if (j > max_y) max_y = j;
				if (i < min_x) min_x = i;
				if (j < min_y) min_y = j;
			}
		}
	}
	int average_x = round(x_sum/num_voxels);
	int average_y = round(y_sum/num_voxels);

	int max_dim = round((max_x - average_x)*2);
	if (max_dim < round(average_x - min_x)*2)
		max_dim = round(average_x - min_x)*2;
	if (max_dim < round(max_y - average_y)*2)
		max_dim = round(max_y - average_y)*2;
	if (max_dim < round(average_y - min_y)*2)
		max_dim = round(average_y - min_y)*2;

//	TR << "Average pos: " << average_x << " " << average_y << std::endl;
//	TR << "MIN, MAX: " << max_x << " " << max_y << " " << min_x << " " << min_y  << std::endl;

	int pad=4;
	int N_new = int(max_dim)+pad;
	std::vector< std::vector<int> > slice_for_ZT_shrink;
	for (int i = 0; i< N_new; i++ ) {
		slice_for_ZT_shrink.push_back( std::vector<int>(N_new) );
	}
	for (int i = 0; i< N; i++ ) {
		for (int j = 0; j< N; j++ ) {
			if  (slice_for_ZT[i][j] == 1) {
				int x=int(i - average_x + N_new/2);
				int y=int(j - average_y + N_new/2);
				slice_for_ZT_shrink[x][y]=1;
			}
		}
	}
	TR << "Size of shrinked grid: " << N_new << std::endl;
	return slice_for_ZT_shrink;
}

void
VoxelGrid::setup_2D_descriptor(int actual_grid_size, int order)
{
	zernike_2D_descriptor_ = zernike::Zernike2Ddescriptor(actual_grid_size, order);
}

void
VoxelGrid::zernike2D_transform(int order) {

	 using namespace basic::options;
   using namespace basic::options::OptionKeys;

	std::string type = option[zernike_descriptor::zernike_transform_type].value();

        bool debug = option[zernike_descriptor::debug].value();

	std::vector< std::vector< std::vector<int> > > slice = xy_plane( actual_grid_dimension_ );
	int slice_dim( slice.size());

	if ( type == "2D_3D" ) {
		zernike_transform_3D_slice(order, slice);
		return;

	}

	// Clear arrays
	edt_surface_vector_.clear();
	std::vector< std::vector<int> > slice_coords;
	std::map< std::vector<int>, std::vector<int> > map;
	for (int i = 0; i< slice_dim; i++ ) {
		for (int j = 0; j< slice_dim; j++ ) {
			std::vector<int> s(slice[i][j]);
			if (s.size() > 0 ) {
				slice_coords.push_back(s);
				std::vector<int> xy;
				xy.push_back(i);
				xy.push_back(j);
				map.insert( std::pair< std::vector<int>, std::vector<int> > (s,xy) );
			}
		}
	}
	std::vector< std::vector<int> > slice_for_ZT;
	for (int i = 0; i< slice_dim; i++ ) {
		slice_for_ZT.push_back( std::vector<int>(slice_dim) );
	}
	int i, j, k;
        for (i = 1; i <= actual_grid_dimension_; i++) {
  	   for (k = 1; k <= actual_grid_dimension_; k++) {
    	      for (j = 1; j <= actual_grid_dimension_; j++) {
      	          if (grid_(i, k, j).is_surface_by_edt() ) {
		     std::vector<int> v = {i,k,j};
                     edt_surface_vector_.push_back(v);
//					std::cout << v[0] << "," << v[1] << "," << v[2]  << std::endl;
				}
			}
		}
	}
    std::sort(edt_surface_vector_.begin(), edt_surface_vector_.end());
    std::sort(slice_coords.begin(), slice_coords.end());
    std::vector< std::vector<int> > v_intersection;

    std::set_intersection(edt_surface_vector_.begin(), edt_surface_vector_.end(),
                          slice_coords.begin(), slice_coords.end(),
                          std::back_inserter(v_intersection));


	float sum_x(0);
	float sum_y(0);
	int num_pixels(0);

  std::vector< std::vector<int> > intersection_plane;
	for(std::vector<int> v : v_intersection) {
			std::vector<int> xy = map[v];
			intersection_plane.push_back(xy);
				sum_x += xy[0];
				sum_y += xy[1];
				num_pixels +=1;
	}

	slice_average_x_ = sum_x/num_pixels;
	slice_average_y_ = sum_y/num_pixels;
	//ADD for alignment
    //std::ofstream savefile("pose.slice.grid.dat");
    int max_dim = 0;
    for(std::vector<int> v : intersection_plane) {
				int x = round( v[0] - slice_average_x_ + slice_dim/2 );
				int y = round( v[1] - slice_average_y_ + slice_dim/2 );
				slice_for_ZT[x][y] = 1;
				if (x > max_dim) max_dim = x;
				if (y > max_dim) max_dim = y;
				if (debug) {
                                	std::ofstream savefile("pose.slice.grid.dat", std::ofstream::app);
					savefile << x << " " << y << std::endl;
				}
    }
		num_surface_2d_pixels_ = v_intersection.size();
		max_surface_dimension_ = max_dim;


		TR.Debug << "Calling setup for zernike_2D_descriptor_ VoxelGrid::zernike2D_transform" << std::endl;
		zernike_2D_descriptor_ = zernike::Zernike2Ddescriptor(slice_dim, order);

		zernike_2D_descriptor_.Transform(slice_for_ZT);
}

void
VoxelGrid::zernike2D_transform_from_stored_slice(int order, numeric::xyzVector<Real> surface_vector1, numeric::xyzVector<Real> surface_vector2 ) {

	using namespace basic::options;
   	using namespace basic::options::OptionKeys;

        bool debug = option[zernike_descriptor::debug].value();

	std::vector< std::vector< std::vector<int> > > slice = slice_plane_from_surface_vectors( surface_vector1, surface_vector2, actual_grid_dimension_ );
	std::vector< std::vector<int> > slice_coords;
	std::map< std::vector<int>, std::vector<int> > map;
	int slice_dim( slice.size());
	for (int i = 0; i< slice_dim; i++ ) {
		for (int j = 0; j< slice_dim; j++ ) {
			std::vector<int> s(slice[i][j]);
			if (s.size() > 0 ) {
				slice_coords.push_back(s);
				std::vector<int> xy;
				xy.push_back(i);
				xy.push_back(j);
				map.insert( std::pair< std::vector<int>, std::vector<int> > (s,xy) );
			}
		}
	}
	std::vector< std::vector<int> > slice_for_ZT;
	for (int i = 0; i< slice_dim; i++ ) {
		slice_for_ZT.push_back( std::vector<int>(slice_dim) );
	}

    if ( edt_surface_vector_.size() == 0 ) {
			int i, j, k;
    	for (i = 1; i <= actual_grid_dimension_; i++) {
        	for (k = 1; k <= actual_grid_dimension_; k++) {
            	for (j = 1; j <= actual_grid_dimension_; j++) {
                	if (grid_(i, k, j).is_surface_by_edt() ) {
											std::vector<int> v = {i,k,j};
                    	edt_surface_vector_.push_back(v);
									}
								}
						}
				}
		}
    std::sort(edt_surface_vector_.begin(), edt_surface_vector_.end());
    std::sort(slice_coords.begin(), slice_coords.end());
    std::vector< std::vector<int> > v_intersection;

    std::set_intersection(edt_surface_vector_.begin(), edt_surface_vector_.end(),
                          slice_coords.begin(), slice_coords.end(),
                          std::back_inserter(v_intersection));

	float sum_x(0);
	float sum_y(0);
	int num_pixels(0);

  std::vector< std::vector<int> > intersection_plane;
	for(std::vector<int> v : v_intersection) {
			std::vector<int> xy = map[v];
			intersection_plane.push_back(xy);
				sum_x += xy[0];
				sum_y += xy[1];
				num_pixels +=1;
	}

	slice_average_x_ = sum_x/num_pixels;
	slice_average_y_ = sum_y/num_pixels;

//    std::ofstream savefile("pose.slice.grid.dat");
    int max_dim = 0;
    for(std::vector<int> v : intersection_plane ) {
				int x = round( v[0] - slice_average_x_ + slice_dim/2 );
				int y = round( v[1] - slice_average_y_ + slice_dim/2 );
				slice_for_ZT[x][y] = 1;
				if (x > max_dim) max_dim = x;
				if (y > max_dim) max_dim = y;
				if (debug) {
                                	std::ofstream savefile("pose.slice.grid.dat", std::ofstream::app);
					savefile << x << " " << y << std::endl;
				}

//				savefile << x << " " << y << std::endl;
		}

	num_surface_2d_pixels_ = v_intersection.size();
	max_surface_dimension_ = max_dim;

    ///transform
//    auto finish = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = finish - start;
//    std::cout << "Time per transform: " << elapsed.count()  << std::endl;
		if ( !zernike_2D_descriptor_.initialized() ) {
				TR << "Calling setup for zernike_2D_descriptor_ VoxelGrid::zernike2D_transform_from_stored_slice " << std::endl;
//    	zernike_2D_descriptor_ = zernike::Zernike2Ddescriptor(slice_for_ZT, slice_dim, order);
//    	zernike_2D_descriptor_ = zernike::Zernike2Ddescriptor(N, order);
		}

		TR << "Transforming with grid length of: " << slice_dim << std::endl;
    zernike_2D_descriptor_ = zernike::Zernike2Ddescriptor( slice_dim, order);
		zernike_2D_descriptor_.Transform(slice_for_ZT);

}


void
VoxelGrid::zernike2D_transform_from_slice(int order) {

	using namespace basic::options;
   	using namespace basic::options::OptionKeys;

	std::string type = option[zernike_descriptor::zernike_transform_type].value();

	std::vector< std::vector< std::vector<int> > > slice = xy_plane( actual_grid_dimension_ );

	if ( type == "2D_3D" ) {
		TR << "Slice transformed with 3D zernike transform: " << std::endl;
		zernike_transform_3D_slice(order, slice);
		return;
	}

	std::vector< std::vector<int> > slice_coords;
	int center_grid ( int(actual_grid_dimension_/2) );
	for (int i = 1; i <= actual_grid_dimension_; i++ ) {
		for (int j = 1; j <= actual_grid_dimension_; j++ ) {
			if ( grid_(i,j,center_grid).is_filled() ) {
				std::vector<int> s{i,j,center_grid};
				slice_coords.push_back(s);
			}
		}
	}

	int slice_dim( actual_grid_dimension_ );
	std::vector< std::vector<int> > slice_for_ZT;
	for (int i = 0; i< slice_dim; i++ ) {
		slice_for_ZT.push_back( std::vector<int>(slice_dim) );
	}

    if ( edt_surface_vector_.size() == 0 ) {
			int i, j, k;
    	for (i = 1; i <= actual_grid_dimension_; i++) {
        	for (k = 1; k <= actual_grid_dimension_; k++) {
            	for (j = 1; j <= actual_grid_dimension_; j++) {
                	if (grid_(i, k, j).is_surface_by_edt() ) {
											std::vector<int> v = {i,k,j};
                    	edt_surface_vector_.push_back(v);
									}
								}
						}
				}
		}
    std::sort(edt_surface_vector_.begin(), edt_surface_vector_.end());
    std::sort(slice_coords.begin(), slice_coords.end());
    std::vector< std::vector<int> > v_intersection;

    std::set_intersection(edt_surface_vector_.begin(), edt_surface_vector_.end(),
                          slice_coords.begin(), slice_coords.end(),
                          std::back_inserter(v_intersection));

	float sum_x(0);
	float sum_y(0);
  int num_pixels(0);

	for(std::vector<int> v : v_intersection) {
				sum_x += v[0];
				sum_y += v[1];
				num_pixels +=1;
	}

	slice_average_x_ = sum_x/num_pixels;
	slice_average_y_ = sum_y/num_pixels;

    //std::ofstream savefile3("grid.slice.dat");
    int max_dim = 0;
    for(std::vector<int> v : v_intersection) {
				int x = round( v[0] - slice_average_x_ + slice_dim/2 );
				int y = round( v[1] - slice_average_y_ + slice_dim/2 );
				slice_for_ZT[x][y] = 1;
				if (x > max_dim) max_dim = x;
				if (y > max_dim) max_dim = y;

		//			savefile3 << x << " " << y << std::endl;
		}

		num_surface_2d_pixels_ = v_intersection.size();
		max_surface_dimension_ = max_dim;

    ///transform
//    auto finish = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = finish - start;
//    std::cout << "Time per transform: " << elapsed.count()  << std::endl;

		TR << "Transforming with grid length of: " << slice_dim << std::endl;
    zernike_2D_descriptor_ = zernike::Zernike2Ddescriptor( slice_dim, order);
		zernike_2D_descriptor_.Transform(slice_for_ZT);

}


void
VoxelGrid::zernike2D_transform_from_outline_slice(std::vector< std::vector<int> > slice_for_ZT, int order) {

	int slice_dim( slice_for_ZT.size());
	zernike_2D_descriptor_ = zernike::Zernike2Ddescriptor( slice_dim, order);
	zernike_2D_descriptor_.Transform(slice_for_ZT);
}


zernike::ZernikeDescriptor<double, double> &
VoxelGrid::get_zernike_descriptor()  {
    return zernike_descriptor_;
}

zernike::Zernike2Ddescriptor &
VoxelGrid::get_zernike_2D_descriptor()  {
    return zernike_2D_descriptor_;
}


void
VoxelGrid::save_invariants( std::string filename ) {
    zernike_descriptor_.SaveInvariants( filename.c_str() );
}


} // namespace shape
} // namespace scoring
} // namespace core
