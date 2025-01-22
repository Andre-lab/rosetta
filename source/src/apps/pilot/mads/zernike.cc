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

static basic::Tracer TR("apps.pilot.zernike_reconstruction");
OPT_KEY(String, xyzfile)
OPT_KEY(String, reconstruction_out)
OPT_KEY(String, voxelgrid_out)
OPT_KEY(String, descriptor_out)
OPT_KEY(Boolean, never_leave_neighbour)

void readPointCloud(const std::string& file,
                    std::vector<std::vector<double>>& atom_coords,
                    std::vector<double>& lj) {
    std::ifstream infile(file);
    if (!infile.is_open()) {
        throw std::runtime_error("Error opening file: " + file);
    }

    std::string line;
    bool is_header = true;

    while (std::getline(infile, line)) {
        if (is_header) {
            // Skip the header line
            is_header = false;
            continue;
        }

        std::istringstream ss(line);
        std::string token;
        std::vector<double> coords(3);
        double lj_value;

        // Parse x, y, z, lj
        for (int i = 0; i < 3; ++i) {
            if (!std::getline(ss, token, ',')) {
                throw std::runtime_error("Error parsing line: " + line);
            }
            coords[i] = std::stod(token); // Convert to double
        }

        if (!std::getline(ss, token, ',')) {
            throw std::runtime_error("Error parsing lj value: " + line);
        }
        lj_value = std::stod(token); // Convert to double

        // Add to the vectors
        atom_coords.push_back(coords);
        lj.push_back(lj_value);
    }
}



int
main(int argc, char** argv)
{
    using namespace core;
    using namespace scoring;
    using namespace shape;
    using namespace utility;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    NEW_OPT(xyzfile, "Name of the file", "");
    NEW_OPT(reconstruction_out, "Name of the resulting reconstruction_out file (in csv format).", "");
    NEW_OPT(voxelgrid_out, "Name the resulting internal voxelgrid_out layout (in json format).", "");
    NEW_OPT(descriptor_out, "Name the resulting Zernike descriptor_out (in csv format)", "");
    NEW_OPT(never_leave_neighbour, "Use the never leave neighbour algorithm", false);

    // zernike outname
    devel::init(argc, argv);

    // read in xyz coordinate file
    std::string file = option[xyzfile];
    std::vector<std::vector<double>> atom_coords;
    std::vector<double> ljs;
    try {
        readPointCloud(file, atom_coords, ljs);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // read in options for the zernike calculation
    int order = option[zernike_descriptor::order].value();
    int gridsize = option[zernike_descriptor::grid_size].value();
    std::string surface_type = option[zernike_descriptor::surface_type].value();
    core::Real probe = option[zernike_descriptor::probe_radius].value();
    core::Real shell = option[zernike_descriptor::shell_thickness].value();

    // read in options for output information
    std::string reconstruction_file_name = option[reconstruction_out];
    std::string voxelgrid_file_name = option[voxelgrid_out];
    std::string descriptors_file_name = option[descriptor_out];
    bool never_leave_neighbour_opt = option[never_leave_neighbour];

    // define calculator
    core::scoring::shape::ZernikeDescriptorCalculator zc =
        core::scoring::shape::ZernikeDescriptorCalculator(order, gridsize, surface_type, probe, shell, 1000,
                                                          never_leave_neighbour_opt);

    // calculate and output the zernike descriptors:
    if (!descriptors_file_name.empty()) {
        zc.invariants_from_vector_to_file(atom_coords, ljs, descriptors_file_name);
        std::cout << "Output the descriptor to: " << descriptors_file_name << std::endl;
    }

    // output the voxelgrid
    if (!voxelgrid_file_name.empty()) {
        // calculate the invariants (which sets the grid)
        zc.invariants_from_vector(atom_coords, ljs);
        // output the voxelgrid_out information
        shape::VoxelGrid vx = zc.get_grid();
        vx.output_grid_json(voxelgrid_file_name);
        std::cout << "Output the voxelgrid file to: " << voxelgrid_file_name << std::endl;
    }

    // output the reconstruction file
    if (!reconstruction_file_name.empty()) {
        // now we have moments in there
        //core::scoring::shape::zernike::ZernikeDescriptor<double, double>::ComplexT3D grid = zc.shape_reconstruction<double,double>(*pose);
        core::scoring::shape::zernike::ZernikeDescriptor<double, double>::ComplexT3D grid = zc.reconstruct_shape_from_vector(atom_coords, ljs);
        core::scoring::shape::zernike::ZernikeMoments<double, double>().SaveGrid(grid, reconstruction_file_name);
        std::cout << "Output the reconstruction file to: " << reconstruction_file_name << std::endl;
    }

    return 0;
}