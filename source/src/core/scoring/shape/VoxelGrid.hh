// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/VoxelGrid.hh
/// @brief      The VoxelGrid class voxelizes a pose and from it produces zernike descriptors (or invariants)
/// @details    The VoxelGrid class main functions and purposes are:
///             1: voxelize(pose) - Voxelizes a pose in a grid.
///             2: find_surface() - Finds the surface voxels of the voxelized pose.
///             3: edt() - Does a Euclidian Distance Transform (EDT) from the surface voxels and marks shell voxels.
///             4: zernike_transform(order) - from the shell voxels calculates zernike descriptors with order.
///             5: invariants() - returns the invariants calculated.
///
///             The main objects of the class are:
///             grid_: A FArray3D that contains a (DIM+2)*(DIM+2)*(DIM+2)  object that acts like the grid.
///             The reason for having +2 in the grid DIM is that the find_surface() function needs at least 1 voxel
///             on each side to move around the voxelized pose.
///             surface_voxels_frontier_: A Queue that initially contains the surface voxels found by
///             find_surface() but is later on populated and trimmed by edt().
///             shell_: A Queue that marks the shell voxels used by zernike_transform() to caluclate zernike descriptors
///             zernike_descriptor_: A core::scoring::shape::zernike::ZernikeDescriptor<double, double>  object that
///             holds the zernike_descriptors calculated from the grid. Can return invariants through invariants().
///
///             Each function and there algorithm is described in detail below.
///
/// @author     Mads Jeppesen
/// @author     Ingemar Andre

#ifndef INCLUDED_core_scoring_shape_VoxelGrid_hh
#define INCLUDED_core_scoring_shape_VoxelGrid_hh

// Utility headers
#include <ObjexxFCL/FArray3D.hh>

// core
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/shape/zernike/ZernikeDescriptor.hh>
#include <core/scoring/shape/zernike/Zernike2Ddescriptor.hh>


// C++ headers
#include <vector>
#include <string>
#include <fstream>
#include <cstring>

// STL
#include <queue>

// other
#include <core/scoring/shape/Neighbour.hh>
//#include <external/pybind11/include/pybind11/pybind11.h>
//#include <external/pybind11/include/pybind11/numpy.h>


namespace core {
namespace scoring {
namespace shape {

/// @brief: Voxel object stored in a VoxelGrid object. Can be filled and visited.
class Voxel {
    // todo: maybe its just better for the VoxelGrid class to be friend with this badboy and just have it acces these directly?
public:

    Voxel();

    ///@brief: constructor with parameters.
    Voxel(int x, int y, int z);

    ///@brief: default destructor.
    ~Voxel() = default;

    ///@brief: mark the voxel as filled.
    void fill();

    ///@brief: mark the voxel as unfilled.
    void unfill();

    ///@brief: check if the voxel has been filled.
    bool is_filled();

    ///@brief: mark the voxel as visited.
    void visit_by_fs();

    ///@brief: check if the voxel has been visited.
    bool is_visited_by_fs();

    ///@brief: mark the voxel as unvisited.
    void unvisit_by_fs();

    ///@brief: mark the voxel as visited.
    void visit_by_edt();

    ///@brief: check if the voxel has been visited.
    bool visited_by_edt();

    ///@brief: mark the voxel as unvisited.
    void unvisit_by_edt();

    ///@brief: mark the voxel as surface.
    void surface_by_fs();

    ///@brief: mark the voxel as surface
    void unsurface_by_fs();

    ///@brief: check if the voxel has been labelled as surface.
    bool is_surface_by_fs();

    ///@brief: mark the voxel as surface.
    void surface_by_edt();

    ///@brief: mark the voxel as surface
    void unsurface_by_edt();

    ///@brief: check if the voxel has been labelled as surface.
    bool is_surface_by_edt();

    ///@brief: returns the x coordinate
    int x() const;

    ///@brief: returns the y coordinate
    int y() const;

    ///@brief: returns the z coordinate
    int z() const;

    ///@brief: returns the edt distance
    core::Real
    edt_dist() const;

    ///@brief: sets the edt distance
    void
    set_edt_dist(core::Real dst);

    ///@brief: returns the edt root x coordinate
    int
    edt_root_x() const;

    ///@brief: returns the edt root y coordinate
    int
    edt_root_y() const;

    ///@brief: returns the edt root z coordinate
    int
    edt_root_z() const;

    ///@brief: sets the edt_root_x, edt_root_y and edt_root_z and from the voxel edt_root_x, edt_root_y and edt_root_z
    void
    edt_root_from_voxel(Voxel const & center);

    ///@brief: returns the distance between the edt_root_x, edt_root_y and edt_root_z and the voxel x, y, and z
    core::Real
    edt_root_distance(Voxel const & neighbour) const;

    ///@brief: clears all data (all bool member variables are set to false)
    void clear();

    ///@brief: comparison operator.
    bool operator == (const core::scoring::shape::Voxel v) const;

private:
    bool filled_;
    bool visited_by_fs_;
    bool surface_by_fs_;
    bool visited_by_edt_;
    bool surface_by_edt_;
    int x_;
    int y_;
    int z_;
    int edt_root_x_;
    int edt_root_y_;
    int edt_root_z_;
    core::Real edt_dist_;
};

///@brief: 3D voxel grid object. Can take voxelize a pose and find which voxels are surface voxels.
class VoxelGrid {
public:


    ///@brief: constructor with parameters
    /// @param surface_type the surface type to use (NOTE: Currently only MS,SAS and VDW surface types are supported. Would like to use electrondensity etc. in the future).
    /// @param probe the probe radii to use.
    VoxelGrid( int size = 64,
               std::string const & surface_type = "MS",
               core::Real probe = 1.4,
               core::Real shell = 2,
               bool recalculate = true,
               bool scale_on = false,
               core::Real scaling = 1);

    ///@brief: default destructor.
    ~VoxelGrid() = default;

    ///@brief: index operator used to find the voxel at position x, y, z.
    core::scoring::shape::Voxel & operator () (int x, int y, int z);

    ///@brief: voxelizes a pose.
    ///
    ///@details: Algorithm:
    /// 1. The pose grid is cleaned for any labelling (filled and visited).
    /// 2. The geometric center of the pose is found and the maximum distance to any atom in the pose to the geometric center is recorded.
    /// 3. This distance is used together with the probe size and maximum lennard jones radii,
    /// to set the grid size in angstrom, and is used to set the grid resolution (aka the size of a voxel in angstrom).
    /// the following formula sets the grid resolution:
//    / OLD:
//    / = (maximum distance to an atom from geometrical center + maximu LJ radii + probe) / (dim/2 - 0.5 - 1 - 0.5 (if dim % 2 != 0))
//    / - dim/2 becuase we consider the geometrical center to atom with the maximum distance (hence half the grid dim).
//    / - dim and not dim+2 becuase we want the pose to only consider the 2 - (dim-1) gridvoxels (see the description at the top of the source file)
//    / - -2 to give some extra cushioning because:
//    /         The distance to the nearest voxel is only half a voxel, and to the next one 1.5 voxel (see the Neighbour class).
//    /         and because 1 voxel should be filled right at either 1 or dim+2 according to the equation above - so we cushion with 1 voxel.
//    / - -0.5 (if dim % 2 != 0) becuase if the bool is true, then the geometrical center will be put at the midpoint of a voxel,
    /// NEW
    /// pretty much same as before but added 1 voxel extra (problem with 6qfj) + removed the geometric_center_at_grid_center // todo. REMOVE THIS VARIABLE IN THE FUTURE
    /// = (maximum distance to an atom from geometrical center + maximu LJ radii + probe) / (dim/2 - 0.5 - 1 - 0.5 (if dim % 2 != 0))
    /// and on that note have +0.5 advantage.
    /// 4. The pose is voxelized by putting the geometrical center at the grid mid point and then iterating through each atom
    /// and filling the voxels that are (atom distance to geomtrical center / grid resolution) away.
    /// 5. A solid is made by filling in neighbouring atom acording to the surface type (probe/no probe) and LJ radii.
    /// All the voxelized voxels at this point are marked "filled=true"
    ///
    ///@param pose the pose to voxelize
    ///
    /// example(s):
    ///     voxelgrid.voxelize(pose)
    void voxelize(core::pose::Pose const & pose);
    /// @brief function to initialize a grid stored from file
    void voxel_grid_from_file( std::string filename);
    void set_voxel_grid( std::vector< std::vector<int> > slice_for_ZT );
    /// @brief
    std::vector< std::vector<int> >
    slice_outline_grid_from_file( std::string filename);

    ///@brief: finds the surface of a voxelized pose.
    ///@Details: Algorithm:
    /// 1. Finds a surface voxel from the voxels marked filled in voxelize() by searching in one
    /// direction starting from the grid edge.
    /// 2. Recursively walks around in 6 directions (not 26 directions, which is much slower) and marks any surface
    /// voxels it hits while also never leaving to far away from the surface.
    /// This approach is a faster modification of a full flood fill algorithm.
    ///
    /// example(s):
    ///     voxelgrid.find_surface()
    void find_surface();


    ///@brief: finds the shell voxels through an edt algorithm.
    ///@Details: Algorithm:
    /// 1. Starts from the surface voxels found in find_surface() that are considered the root of the EDT algorithm.
    /// 2. Propagating inwards ("propagiting the frontier") taking 1 voxel at the time, it sets the x, y and z
    /// coordinate of the voxel in the grid to the closest root voxels and calculates the distances of all neighbours
    /// with respect to the root x,y, and z.
    /// 3. It adds voxels to the shell if allowed by the surface type (MS, SAS, VDW etc.)
    /// 4. It stops the algorithm after a max edt seach has been conducted with again depends on the surface type.
    ///
    /// example(s):
    ///     voxelgrid.edt()
    void
    edt();

    ///@brief: Does a zernike transform
    ///@Details: Algorithm:
    /// 1. top secret :O ....
    ///
    /// example(s):
    ///     voxelgrid.zernike_transform()
    void
    zernike_transform(int order);

	void
	zernike_transform_3D_slice(int order, std::vector< std::vector< std::vector<int> > > slice);

   ///@brief: Does a zernike transform
    ///@Details: Algorithm:
    /// 1. top secret :O ....
    ///
    /// example(s):
    ///     voxelgrid.zernike_transform()

	std::vector< std::vector<int> >
	shrink_2D_grid( std::vector< std::vector<int> > slice_for_ZT );

	void
	setup_2D_descriptor(int actual_grid_size, int order);

    void
    zernike2D_transform(int order);

	void
	zernike2D_transform_from_slice(int order);

 void
  zernike2D_transform_from_stored_slice(
      int order,
      numeric::xyzVector<Real> surface_vector1,
      numeric::xyzVector<Real> surface_vector2 );

	void
	zernike2D_transform_from_outline_slice(std::vector< std::vector<int> > slice_for_ZT, int order);
    ///@brief: returns the filled voxels
    std::vector<std::vector<int>>
    get_filled_voxels();

    ///@brief: returns the visit_by_fs voxels
    std::vector<std::vector<int>>
    get_visit_by_fs_voxels();

    ///@brief: returns the surface_by_fs voxels
    std::vector<std::vector<int>>
    get_surface_by_fs_voxels();

    ///@brief: returns the visit_by_edt voxels
    std::vector<std::vector<int>>
    get_visit_by_edt_voxels();

    ///@brief: returns the surface_by_edt voxels
    std::vector<std::vector<int>>
    get_surface_by_edt_voxels();


    ///@brief: outputs all voxel information data (the same as the 5 functions above) into a csv format.
    void output_grid_csv(std::string name);

    //TODO: SHOULD BE EXTENDED AND TESTED MORE!
    ///@brief: outputs voxel information data of the filled voxels.
    void output_grid_binvox(std::string name);

//    std::map<std::string, std::vector<std::tuple<core::Real,core::Real, core::Real>>>
//    get_grid_map();

    int
    actual_grid_size();

    int
    grid_size();

    std::string
    surface_type();

    core::Real
    grid_resolution();

    core::Real
    probe_radius();

    core::Real
    shell_thickness();

    void
    set_surface_vectors( numeric::xyzVector<Real>& surface_vector1,
													numeric::xyzVector<Real>& surface_vector2 );

    numeric::xyzVector<Real>
    slice_surface_vector1();

    numeric::xyzVector<Real>
    slice_surface_vector2();

		core::Real
		slice_average_x();

		core::Real
		slice_average_y();

    void
    set_scaling (core::Real scaling);

    void
    set_grid_size(int);

    void
    set_surface_type(std::string const & surface_type);


    void
    set_probe_radius( core::Real radius );

    void
    set_shell_thickness( core::Real thickness );

    bool
    sanity_checks();

    bool
    is_utilizing_the_whole_grid();

//    pybind11::array_t<int>
//    get_grid_map();

    // not visible in pyrosetta
//    std::map<std::string,  std::vector<std::tuple<int, int, int>>>
//    get_grid_map();
//
//    std::map<std::string,  std::vector<std::vector<int>>>
//    get_grid_mapping();




    void
    save_invariants( std::string filename );

    zernike::ZernikeDescriptor<double, double> &
    get_zernike_descriptor();

    zernike::Zernike2Ddescriptor &
    get_zernike_2D_descriptor();

    numeric::xyzMatrix_real
    random_reorientation_matrix(const core::Real phi_range= 360.0, const core::Real psi_range= 360.0);

    void sample_random_surface_vectors();

		core::pose::PoseOP
		rotate_pose( core::pose::Pose const & pose);

    core::pose::PoseOP
    rotate_pose( core::pose::Pose const & pose, numeric::xyzVector<Real> surface_vector1,
                                          numeric::xyzVector<Real> surface_vector2 );

		std::vector< std::vector< std::vector<int> > >
		xy_plane(core::Size N);

    std::vector< std::vector< std::vector<int> > >
    sample_slice_plane( core::Size N);

    std::vector< std::vector< std::vector<int> > >
    slice_plane_from_surface_vectors(
        numeric::xyzVector<Real> surface_vector1,
        numeric::xyzVector<Real> surface_vector2,
       core::Size N);

    ///@brief: return invariants
    std::vector<double> &
    invariants()  {
        return zernike_descriptor_.GetInvariants();
    };

    ///@brief: return invariants 2D
    std::vector<double> &
    invariants2D()  {
        return zernike_2D_descriptor_.GetInvariants();
    };

   ///@brief: How many surface 2D pixels do we have?
   core::Size
   num_surface_pixels() {
       return num_surface_2d_pixels_;
   };

   ///@brief: How many surface 2D pixels do we have?
   core::Size
   max_surface_dim() {
       return max_surface_dimension_;
   };

private:

    void
    initialize_grid();

    void
    generate_neighbours();

    ///@brief: unfills and unvisists all voxels in the grid.
    void
    clean_grid();

    ///@brief: detects if neighbours (all 26 of them) adjacent to the voxel at postion x, y, z are filled.
    bool
    detect_neighbours(int x, int y, int z);

    ///@brief: mark neighbours (6 of the closest onces) adjacent to the voxel at postion x, y, z as filled.
    void
    mark_surface(int x, int y, int z);

    ///@brief: fills voxels that are on the surface. Only marks the 6 closest ones to the closest unfilled voxel.
    void
    fillsurface(int x, int y, int z);

    ///@brief: finds 1 surface voxel. Is used to initialize the start of the fillsurface method.
    std::tuple<int, int, int>
    find_edge();

    /// @brief checks if any filled voxels are in 1 or dim+2 voxels. Returns true if not.
    bool
    grid_not_filled_at_boundary();

    bool
    is_xyz_at_grid_edge(int x, int y, int z);

    bool
    is_xyz_out_of_grid(int x, int y, int z);

    bool
    surface_by_edt_criteria(core::Real dist);

    void
    set_edt_search_depth();

    bool
    search_in_z_dimension(int xi, int yi, int & zi);

    ObjexxFCL::FArray3D<Voxel> grid_;
    std::queue<Voxel*> surface_voxels_frontier_;
    std::queue<Voxel*> shell_;
    int grid_dimension_;
    int actual_grid_dimension_;
    int geometric_center_at_grid_center_;
    std::vector<Neighbour> neighbours_;
    std::vector<int> identical_neighbours_;
    core::Real scaling_;
    bool scale_on_;
    int amount_of_voxel_neighbours_;
    core::Real grid_res_;
    core::Real probe_radius_;
    core::Real shell_thickness_;
    bool recalculate_;
    int surface_type_; // MS = 1, SAS = 2, VDW = 3
    core::Real max_edt_search_; // a constant that will limit the EDT search in the edt() function.
    zernike::ZernikeDescriptor<double, double>  zernike_descriptor_; //todo: check why this is a double, double? thought it was 1d??
    zernike::Zernike2Ddescriptor zernike_2D_descriptor_;
    numeric::xyzVector<Real> slice_vector1_, slice_vector2_;
    core::Real slice_average_x_, slice_average_y_;
    std::vector< std::vector<int> > edt_surface_vector_;
    core::Real padding_;
    core::Size num_surface_2d_pixels_;
    core::Size max_surface_dimension_;
//    double* zernike_array_; // new double[gridsize*gridsize*gridsize];
};

} // namespace shape
} // namespace scoring
} // namespace core

#endif //INCLUDED_core_scoring_shape_VoxelGrid_hh
