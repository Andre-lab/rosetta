// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/shape/ZernikeDescriptorEnergy.hh
/// @brief
/// @details
/// @author Ingmear Andre
/// @author Mads Jeppesen


#ifndef INCLUDED_core_scoring_shape_ZernikeDescriptorEnergy_hh
#define INCLUDED_core_scoring_shape_ZernikeDescriptorEnergy_hh

#include <complex>
#include <core/scoring/shape/ZernikeDescriptorData.hh>
#include <core/scoring/shape/VoxelGrid.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <basic/datacache/CacheableData.hh>



namespace core {
namespace scoring {
namespace shape {

static basic::Tracer TR( "core.scoring.shape.ZernikeDescriptorEnergy"  );

class ZernikeDescriptorCalculator {
public:

   ZernikeDescriptorCalculator(
        int order=20,
        int size=64,
        std::string const & surface_type = "MS",
        core::Real probe = 1.4,
        core::Real shell = 2,
        int num_slices = 10000);

    std::vector<double>
    invariants_from_pose(core::pose::Pose const & pose);

  std::vector<std::complex<double>>
  moments_from_pose(core::pose::Pose const & pose);

    std::vector<double>
    invariants_slice_from_pose_and_surface_vectors(
       core::pose::Pose const & pose,
       numeric::xyzVector<Real> surface_vector1,
       numeric::xyzVector<Real> surface_vector2 );

	std::vector<double>
    invariants_from_grid_file(std::string filename);


    std::vector<double>
    invariants_from_outline_grid_file(std::string filename);

    std::vector< std::vector<int> >
    read_slice_grid(std::string filename);

    std::vector<double>
    invariants_slice_from_pose(core::pose::Pose const & pose);

    std::vector<double>
    invariants_slice_from_pose_oriented(core::pose::Pose const & pose);

    std::vector<double>
    invariants_from_file(std::string const & filename);

   std::vector<double>
   invariants_from_grid( std::vector< std::vector<int> > slice_for_ZT );

    //  not found in pyrosetta
//    std::vector<double> &&
//    invariants_from_pose(core::pose::Pose const & pose);

    // todo: compile with pybind11 library
//    pybind11::array_t<double>
//    invariants_from_pose_numpy(core::pose::Pose const & pose);

        //  not found in pyrosetta
//    std::vector<double> &&
//    invariants_from_file(std::string const & filename);

    void
    invariants_from_pose_to_file(core::pose::Pose const & pose, std::string const & filename);

	void
    invariants_2D_3D_from_pose_to_file(core::pose::Pose const & pose, std::string const & filename);

   void
    invariants_from_oriented_pose_to_file(core::pose::Pose const & pose, std::string const & filename);

    void
    invariants_2D_from_pose_and_slice(
        core::pose::Pose const & pose,
        numeric::xyzVector<Real> surface_vector1,
        numeric::xyzVector<Real> surface_vector2,
        std::string const & filename );

    void
    invariants_2D_3D_from_pose_and_slice(
        core::pose::Pose const & pose,
        numeric::xyzVector<Real> surface_vector1,
        numeric::xyzVector<Real> surface_vector2,
        std::string const & filename );

    void
    invariants_2D_from_pose_to_file(core::pose::Pose const & pose, std::string const & filename);

    void
    invariants_from_grid_file_to_file(std::string input_filename, std::string const & filename);

    void
    invariants_from_grid_outline_file_to_file(std::string input_filename, std::string const & filename);

   std::vector<core::Real> get_3d_descriptors();
   std::vector<core::Real> get_2d_descriptors();

    std::vector<std::complex<double>> get_3d_moments();
    std::vector<std::complex<double>> get_2d_moments();

    VoxelGrid &
    grid();

private:
    int order_;
    core::scoring::shape::VoxelGrid voxelgrid_;
    int num_slices_;
};

class ZernikeDescriptorEnergy : public methods::WholeStructureEnergy  {
public:

    typedef methods::WholeStructureEnergy  parent;

public:

    ZernikeDescriptorEnergy();

    //clone
    methods::EnergyMethodOP
    clone() const override;

    /////////////////////////////////////////////////////////////////////////////
    // scoring
    /////////////////////////////////////////////////////////////////////////////

    // todo: has changed this from  const
    void
    finalize_total_energy(
            pose::Pose & pose,
            ScoreFunction const &,
            EnergyMap & totals
    ) const override;

    void
    setup_for_derivatives(
            pose::Pose &,
            ScoreFunction const &
    )
    const override;

    void
    setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const override;

    void
    eval_atom_derivative(
            id::AtomID const &,
            pose::Pose const &,
            kinematics::DomainMap const &,
            ScoreFunction const &,
            EnergyMap const &,
            Vector &,// F1,
            Vector & // F2
    ) const override;

    void
    setup_zernike_descriptor_reference( core::pose::Pose & pose ) const;

    core::Real
    compute(core::pose::Pose const & pose, bool report) const;

    core::Size version() const override;

    void
    indicate_required_context_graphs(
            utility::vector1< bool > & /*context_graphs_required*/
    ) const override {};

    std::vector<double> const &
    zernike_reference_data_pose(core::pose::Pose const & pose) const;

private:
    void
    report_on_current_invariants(std::vector<double> const & zd_current) const {
        for ( auto i : zd_current)
            TR << i << std::endl;
    };

    void
    report_on_current_L2_norm(core::Real L2_norm) const {
        TR << "Zernike: L2_norm = " << sqrt(L2_norm) << std::endl;
    };

private:
    mutable ZernikeDescriptorCalculator zd_calc_;
    mutable ZernikeDescriptorData zd_reference_;
};

} //shape
} //scoring
} //core

#endif // INCLUDED_core_scoring_shape_ZernikeDescriptorEnergy_HH
