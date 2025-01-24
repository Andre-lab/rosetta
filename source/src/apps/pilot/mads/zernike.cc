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

OPT_KEY(String, xyzlist)

OPT_KEY(String, reconstruction_out)

OPT_KEY(String, voxelgrid_out)

OPT_KEY(String, descriptor_out)

OPT_KEY(Boolean, never_leave_neighbour)

void readPointCloud(const std::string& file,
                    std::vector<std::vector<double>>& atom_coords,
                    std::vector<double>& lj)
{
    std::ifstream infile(file);
    if (!infile.is_open())
    {
        throw std::runtime_error("Error opening file: " + file);
    }

    std::string line;
    bool is_header = true;

    while (std::getline(infile, line))
    {
        if (is_header)
        {
            // Skip the header line
            is_header = false;
            continue;
        }

        std::istringstream ss(line);
        std::string token;
        std::vector<double> coords(3);
        double lj_value;

        // Parse x, y, z, lj
        for (int i = 0; i < 3; ++i)
        {
            if (!std::getline(ss, token, ','))
            {
                throw std::runtime_error("Error parsing line: " + line);
            }
            coords[i] = std::stod(token); // Convert to double
        }

        if (!std::getline(ss, token, ','))
        {
            throw std::runtime_error("Error parsing lj value: " + line);
        }
        lj_value = std::stod(token); // Convert to double

        // Add to the vectors
        atom_coords.push_back(coords);
        lj.push_back(lj_value);
    }
}

std::string modifyFileName(const std::string& original_name, size_t index)
{
    // Find the last dot in the file name to identify the suffix
    size_t last_dot = original_name.find_last_of('.');

    // Remove the suffix if it exists
    std::string base_name;
    std::string extension;

    if (last_dot != std::string::npos)
    {
        base_name = original_name.substr(0, last_dot); // Get the base name without suffix
        extension = original_name.substr(last_dot); // Get the extension (including dot)
    }
    else
    {
        base_name = original_name; // If no dot, use the whole name
        extension = ""; // No extension
    }

    // Create the new file name
    return base_name + "_" + std::to_string(index) + extension;
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
    NEW_OPT(xyzlist, "Name of the file for a list that contains xyzfiles to use", "");
    NEW_OPT(reconstruction_out, "Name of the resulting reconstruction_out file (in csv format).", "");
    NEW_OPT(voxelgrid_out, "Name the resulting internal voxelgrid_out layout (in json format).", "");
    NEW_OPT(descriptor_out, "Name the resulting Zernike descriptor_out (in csv format)", "");
    NEW_OPT(never_leave_neighbour, "Use the never leave neighbour algorithm", false);

    // zernike outname
    devel::init(argc, argv);

    // read in xyz coordinate file
    std::string file = option[xyzfile];
    std::string list = option[xyzlist];
    assert((!file.empty() ^ !list.empty()) && "Exactly one of xyzfile or xyzlist must be set");

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

    // Print all variables to the screen
    std::cout << "Running Zernike calculations with the following options: " << file << std::endl;
    std::cout << "xyzfile: " << file << std::endl;
    std::cout << "xyzlist: " << list << std::endl;
    std::cout << "order: " << order << std::endl;
    std::cout << "gridsize: " << gridsize << std::endl;
    std::cout << "surface_type: " << surface_type << std::endl;
    std::cout << "probe: " << probe << std::endl;
    std::cout << "shell: " << shell << std::endl;
    std::cout << "reconstruction_out: " << reconstruction_file_name << std::endl;
    std::cout << "voxelgrid_out: " << voxelgrid_file_name << std::endl;
    std::cout << "descriptors_out: " << descriptors_file_name << std::endl;
    std::cout << "never_leave_neighbour_opt: " << (never_leave_neighbour_opt ? "true" : "false") << std::endl;

    // define calculator
    core::scoring::shape::ZernikeDescriptorCalculator zc =
        core::scoring::shape::ZernikeDescriptorCalculator(order, gridsize, surface_type, probe, shell, 1000,
                                                          never_leave_neighbour_opt);

    // generate output names
    // Declare files as a vector of vector of strings
    std::vector<std::vector<std::string>> file_names;
    std::vector<std::vector<std::vector<double>>> all_atom_coords;
    std::vector<std::vector<double>> all_ljs;

    // Generate output names
    if (!file.empty())
    {
        // read in files
        std::vector<std::vector<double>> atom_coords;
        std::vector<double> ljs;
        try
        {
            readPointCloud(file, atom_coords, ljs);
        }
        catch (const std::exception& e)
        {
            std::cerr << e.what() << std::endl;
            return 1;
        }
        all_atom_coords.push_back(atom_coords);
        all_ljs.push_back(ljs);
        // read in files
        file_names.push_back({descriptors_file_name, voxelgrid_file_name, reconstruction_file_name});
    }
    else
    {
        // Open the file containing the list of file paths
        std::ifstream infile(list);
        if (!infile.is_open())
        {
            std::cerr << "Failed to open file: " << list << std::endl;
            return 1;
        }

        int count = 0;
        std::string line;
        while (std::getline(infile, line))
        {
            if (line.empty())
            {
                continue; // Skip empty lines
            }

            std::vector<std::vector<double>> atom_coords;
            std::vector<double> ljs;

            try
            {
                // Read data from the file
                readPointCloud(line, atom_coords, ljs);
                // Append the results to the aggregated data
                // all_atom_coords.insert(all_atom_coords.end(), atom_coords.begin(), atom_coords.end());
                // all_ljs.insert(all_ljs.end(), ljs.begin(), ljs.end());
                all_atom_coords.push_back(atom_coords);
                all_ljs.push_back(ljs);
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error processing file " << line << ": " << e.what() << std::endl;
            }
            // add new file names
            std::string new_descriptors_file_name;
            std::string new_voxelgrid_file_name;
            std::string new_reconstruction_file_name;
            if (!descriptors_file_name.empty()) {
                new_descriptors_file_name = modifyFileName(descriptors_file_name, count);
            } else {
                new_descriptors_file_name = descriptors_file_name;
            }
            if (!voxelgrid_file_name.empty()) {
                new_voxelgrid_file_name = modifyFileName(voxelgrid_file_name, count);
            } else            {
                new_voxelgrid_file_name = voxelgrid_file_name;
            }
            if (!reconstruction_file_name.empty()) {
                new_reconstruction_file_name = modifyFileName(reconstruction_file_name, count);
            } else {
                new_reconstruction_file_name = reconstruction_file_name;
            }

            count += 1;
            file_names.push_back({new_descriptors_file_name, new_voxelgrid_file_name, new_reconstruction_file_name});
        }
        infile.close();
    }

    //  calculate things:
    // Ensure all sizes match
    if (file_names.size() != all_atom_coords.size() || file_names.size() != all_ljs.size())
    {
        std::cerr << "Error: Size mismatch between files, all_atom_coords, and all_ljs!" << std::endl;
        return 1;
    }

    // Loop over the outer vectors
    for (size_t i = 0; i < file_names.size(); ++i)
    {
        // get all objects
        std::vector<std::vector<double>> atom_coords = all_atom_coords[i];
        std::vector<double> ljs = all_ljs[i];
        std::string descriptors_file_name = file_names[i][0];
        std::string voxelgrid_file_name = file_names[i][1];
        std::string reconstruction_file_name = file_names[i][2];

        // calculate and output the zernike descriptors:
        if (!descriptors_file_name.empty())
        {
            zc.invariants_from_vector_to_file(atom_coords, ljs, descriptors_file_name);
            std::cout << "Output the descriptor to: " << descriptors_file_name << std::endl;
        }

        // output the voxelgrid
        if (!voxelgrid_file_name.empty())
        {
            // calculate the invariants (which sets the grid)
            zc.invariants_from_vector(atom_coords, ljs);
            // output the voxelgrid_out information
            shape::VoxelGrid vx = zc.get_grid();
            vx.output_grid_json(voxelgrid_file_name);
            std::cout << "Output the voxelgrid file to: " << voxelgrid_file_name << std::endl;
        }

        // output the reconstruction file
        if (!reconstruction_file_name.empty())
        {
            // now we have moments in there
            //core::scoring::shape::zernike::ZernikeDescriptor<double, double>::ComplexT3D grid = zc.shape_reconstruction<double,double>(*pose);
            core::scoring::shape::zernike::ZernikeDescriptor<double, double>::ComplexT3D grid = zc.
                reconstruct_shape_from_vector(atom_coords, ljs);
            core::scoring::shape::zernike::ZernikeMoments<double, double>().SaveGrid(grid, reconstruction_file_name);
            std::cout << "Output the reconstruction file to: " << reconstruction_file_name << std::endl;
        }
    }

    return 0;
}
