//
// Created by Mads Jeppesen on 4/25/24.
//

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <string>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <sstream>

// zernike
#include <core/scoring/shape/ZernikeDescriptorEnergy.hh>
#include <core/scoring/shape/VoxelGrid.hh>
#include <core/scoring/shape/zernike/ZernikeMoments.hh>
#include <core/scoring/shape/zernike/ZernikeDescriptor.hh>

// stolen from ingemars app (zernike align)

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/string_util.hh>

// core
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/util.hh>
#include <core/scoring/shape/VoxelGrid.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <iostream>


// protocols
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/scoring/shape/ZernikeDescriptorEnergy.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/zernike_descriptor.OptionKeys.gen.hh>

static basic::Tracer TR( "apps.pilot.zernike_reconstruction" );
OPT_KEY( String, reconstruction)
OPT_KEY( String, voxelgrid)
OPT_KEY( String, descriptors)
OPT_KEY( Boolean, never_leave_neighbour)

int
main( int argc, char ** argv ) {

    using namespace core;
    using namespace scoring;
    using namespace shape;
    using namespace utility;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;


    NEW_OPT( reconstruction, "Name of the resulting reconstruction file (in csv format).", "");
    NEW_OPT( voxelgrid, "Name the resulting internal voxelgrid layout (in json format).", "");
    NEW_OPT( descriptors, "Name the resulting Zernike descriptors (in csv format)", "");

    // todo: this should be a key under zernike_descriptor like the others
    NEW_OPT( never_leave_neighbour, "Use the never leave neighbour algorithm", false);

    // zernike outname
    devel::init( argc, argv );

    // read in pose
    vector1<std::string> files = option[in::file::s]();
    std::string file = files[1]; // i dont care about any other files specified here
    core::pose::PoseOP pose = core::import_pose::pose_from_file(file);

    // read in options for the zernike calculation
    int order = option[zernike_descriptor::order].value();
    int gridsize = option[zernike_descriptor::grid_size].value();
    std::string surface_type = option[zernike_descriptor::surface_type].value();
    core::Real probe = option[zernike_descriptor::probe_radius].value();
    core::Real shell = option[zernike_descriptor::shell_thickness].value();

    // read in options for output information
    std::string reconstruction_file_name = option[reconstruction];
    std::string voxelgrid_file_name = option[voxelgrid];
    std::string descriptors_file_name = option[descriptors];
    bool never_leave_neighbour_opt = option[never_leave_neighbour];

    // define calculator
    core::scoring::shape::ZernikeDescriptorCalculator zc =
            core::scoring::shape::ZernikeDescriptorCalculator(order, gridsize, surface_type, probe, shell, 1000, never_leave_neighbour_opt);

    if (!descriptors_file_name.empty()) {
        zc.invariants_from_pose_to_file(*pose, descriptors_file_name);
        std::cout << "Output the descriptors to: " << descriptors_file_name << std::endl;
    }

    if (!reconstruction_file_name.empty()) {
        // now we have moments in there
        //core::scoring::shape::zernike::ZernikeDescriptor<double, double>::ComplexT3D grid = zc.shape_reconstruction<double,double>(*pose);
        core::scoring::shape::zernike::ZernikeDescriptor<double, double>::ComplexT3D grid = zc.shape_reconstruction(
                *pose);
        core::scoring::shape::zernike::ZernikeMoments<double, double>().SaveGrid(grid, reconstruction_file_name);
        std::cout << "Output the reconstruction to: " << reconstruction_file_name << std::endl;
    }

    if (!voxelgrid_file_name.empty()) {
        // calculate the invariants (which sets the grid)
        zc.invariants_from_pose(*pose);
        // output the voxelgrid information
        shape::VoxelGrid vx = zc.get_voxelgrid();
        vx.output_grid_json(voxelgrid_file_name);
        std::cout << "Output the voxelgrid to: " << voxelgrid_file_name << std::endl;
    }


    return 0;
}
