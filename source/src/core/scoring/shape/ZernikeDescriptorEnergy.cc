// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/shape/ZernikeDescriptorEnergy.cc
/// @brief
/// @author Ingemar Andre
/// @author Mads Jeppesen


//Unit headers
#include <core/scoring/shape/ZernikeDescriptorData.hh>
#include <core/scoring/shape/ZernikeDescriptorEnergy.hh>
#include <core/scoring/shape/ZernikeDescriptorEnergyCreator.hh>
#include <core/scoring/shape/VoxelGrid.hh>

//Package headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <core/id/AtomID.hh>
#include <basic/options/keys/zernike_descriptor.OptionKeys.gen.hh>

/// Utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>

//Auto Headers
/// #include <core/pose/util.tmpl.hh>

#include <utility/io/izstream.hh>
#include <utility/file/gzip_util.hh>

//STL
#include <iomanip>
#include <ctime>
#include <chrono>

#include <core/import_pose/import_pose.hh>

//#include <external/pybind11/include/pybind11/pybind11.h>
//#include <external/pybind11/include/pybind11/numpy.h>
//#include <external/pybind11/include/pybind11/cast.h>

// term specific headers


using namespace core::scoring;
using namespace core::scoring::methods;

namespace core {
namespace scoring {
namespace shape {

ZernikeDescriptorCalculator::ZernikeDescriptorCalculator(
    int order,
    int size,
    std::string const & surface_type,
    core::Real probe,
    core::Real shell,
    int num_slices,
    bool never_leave_neighbour):
    order_(order),
    voxelgrid_(size, surface_type, probe, shell, never_leave_neighbour),
    num_slices_(num_slices) {}

std::vector<double>
ZernikeDescriptorCalculator::invariants_from_pose(core::pose::Pose const & pose){
    voxelgrid_.voxelize(pose);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
    voxelgrid_.zernike_transform(order_);
    return voxelgrid_.invariants();
}

std::vector<std::complex<double>>
ZernikeDescriptorCalculator::moments_from_pose(core::pose::Pose const & pose){
    voxelgrid_.voxelize(pose);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
    voxelgrid_.zernike_transform(order_);
    return voxelgrid_.get_zernike_descriptor().GetMoments();
}

std::vector<double>
ZernikeDescriptorCalculator::invariants_from_vector(std::vector<std::vector<core::Real>> atom_coords,
std::vector<core::Real> ljs) {
    voxelgrid_.voxelize(atom_coords, ljs);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
    voxelgrid_.zernike_transform(order_);
    return voxelgrid_.invariants();
}

std::vector<std::complex<double>>
    ZernikeDescriptorCalculator::moments_from_vector(std::vector<std::vector<core::Real>> atom_coords,
    std::vector<core::Real> ljs) {
    voxelgrid_.voxelize(atom_coords, ljs);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
    voxelgrid_.zernike_transform(order_);
    return voxelgrid_.get_zernike_descriptor().GetMoments();
}

std::vector<std::vector<std::vector<std::complex<double>>>>
ZernikeDescriptorCalculator::reconstruct_shape_from_pose(core::pose::Pose const & pose){
    voxelgrid_.voxelize(pose);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
    voxelgrid_.zernike_transform(order_);
    return voxelgrid_.get_zernike_descriptor().ReconstructionWrapper();
}

std::vector<std::vector<std::vector<std::complex<double>>>>
ZernikeDescriptorCalculator::reconstruct_shape_from_vector(std::vector<std::vector<core::Real>> atom_coords,
    std::vector<core::Real> ljs){
    voxelgrid_.voxelize(atom_coords, ljs);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
    voxelgrid_.zernike_transform(order_);
    return voxelgrid_.get_zernike_descriptor().ReconstructionWrapper();
}


std::vector<double>
ZernikeDescriptorCalculator::invariants_slice_from_pose_and_surface_vectors(core::pose::Pose const & pose, numeric::xyzVector<Real> surface_vector1, numeric::xyzVector<Real> surface_vector2){
		using namespace basic::options;
   		using namespace basic::options::OptionKeys;
	
        	bool debug = option[zernike_descriptor::debug].value();

		TR.Debug << "Sampling slize..." << std::endl;
		TR.Debug << "Updating voxelgrid..." << std::endl;
		voxelgrid_.set_surface_vectors( surface_vector1, surface_vector2 );
		core::pose::PoseOP rotated_pose = voxelgrid_.rotate_pose(pose);
    // trigger refolding of chain. UGLY!!!
		numeric::xyzVector< Real > coord = rotated_pose->residue(1).atom(rotated_pose->residue(1).atom_index("CA")).xyz();
		if (debug) {
			rotated_pose->dump_pdb("initial_rot_pose.pdb");
		}
		voxelgrid_.voxelize(*rotated_pose);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
    voxelgrid_.output_grid_csv("voxelgrid.slice.");
		voxelgrid_.zernike2D_transform(order_);
    return voxelgrid_.invariants2D();
}

std::vector<double>
ZernikeDescriptorCalculator::invariants_from_outline_grid_file(std::string filename){
    TR.Debug << "Reading outline grid..." << std::endl;
    std::vector< std::vector<int> > slice_for_ZT = voxelgrid_.slice_outline_grid_from_file(filename);
    voxelgrid_.zernike2D_transform_from_outline_slice(slice_for_ZT, order_);
    return voxelgrid_.invariants2D();
}

std::vector<double>
ZernikeDescriptorCalculator::invariants_from_grid( std::vector< std::vector<int> > slice_for_ZT ){
    voxelgrid_.set_voxel_grid(slice_for_ZT);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
		//voxelgrid_.output_grid_csv("voxelgrid.sample_grid");
    voxelgrid_.zernike2D_transform_from_slice(order_);
    return voxelgrid_.invariants2D();
}

std::vector<double>
ZernikeDescriptorCalculator::invariants_slice_from_pose(core::pose::Pose const & pose) {
   using namespace basic::options;
   using namespace basic::options::OptionKeys;
   std::string type = option[zernike_descriptor::zernike_transform_type].value();
   bool debug = option[zernike_descriptor::debug].value();

    TR.Debug << "Sampling slize..." << std::endl;
    TR.Debug << "Rotating pose..." << std::endl;
    voxelgrid_.sample_random_surface_vectors();
    core::pose::PoseOP rotated_pose = voxelgrid_.rotate_pose(pose);
    TR.Debug << "Updating voxelgrid..." << std::endl;
    // trigger refolding of chain. UGLY!!!
    numeric::xyzVector<Real> coord = rotated_pose->residue(1).atom(rotated_pose->residue(1).atom_index("CA")).xyz();
    voxelgrid_.voxelize(*rotated_pose);
    voxelgrid_.find_surface();
    voxelgrid_.edt();
    voxelgrid_.zernike2D_transform(order_);
    // For alignment
    if (debug) {
    	rotated_pose->dump_pdb("rotated.pdb");
    }
		if ( type == "2D_3D" ) {
			return voxelgrid_.invariants();
		} else {
		return voxelgrid_.invariants2D();
	}

}

std::vector<double>
ZernikeDescriptorCalculator::invariants_slice_from_pose_oriented(core::pose::Pose const & pose){
 	using namespace basic::options;
	using namespace basic::options::OptionKeys;
 	std::string type = option[zernike_descriptor::zernike_transform_type].value();

		TR.Debug << "Updating voxelgrid..." << std::endl;
		voxelgrid_.voxelize(pose);
		voxelgrid_.find_surface();
		voxelgrid_.edt();
		voxelgrid_.setup_2D_descriptor(voxelgrid_.actual_grid_size(), order_);
		//voxelgrid_.output_grid_csv("voxelgrid.pose");
		voxelgrid_.zernike2D_transform(order_);
		if ( type == "2D_3D" ) {
			return voxelgrid_.invariants();
		} else {
		        return voxelgrid_.invariants2D();
	}
}

// std::vector<double>
// ZernikeDescriptorCalculator::invariants_from_grid_file(std::string filename){
//
// 	 using namespace basic::options;
//    using namespace basic::options::OptionKeys;
//
//     std::cout << "Reading voxelgrid..." << std::endl;
//     voxelgrid_.voxel_grid_from_file(filename);
//     voxelgrid_.find_surface();
//     voxelgrid_.edt();
// 		//voxelgrid_.output_grid_csv("voxelgrid");
// 		std::string type = option[zernike_descriptor::zernike_transform_type].value();
// 		if ( type == "2D" ) {
//       std::cout << "2D zernike transfrom of grid..." << std::endl;
// 			voxelgrid_.zernike2D_transform_from_slice(order_);
// 			return voxelgrid_.invariants2D();
// 		} else {
//       std::cout << "2D 3D zernike transform of grid..." << std::endl;
// 			voxelgrid_.zernike2D_transform_from_slice(order_);
// 			return voxelgrid_.invariants();
// 	}
// }


std::vector< std::vector<int> >
ZernikeDescriptorCalculator::read_slice_grid(std::string filename){
    std::vector< std::vector<int> > slice_for_ZT = voxelgrid_.slice_outline_grid_from_file(filename);
    return slice_for_ZT;
}

std::vector<double>
ZernikeDescriptorCalculator::invariants_from_file(std::string const & filename){
    std::vector<double> zd;
    std::ifstream infile(filename.c_str());
    std::string line;
    int counter(0);
    while ( getline(infile, line) ) {
        counter++;
        std::istringstream line_stream(line);
        if ( counter == 1 ) {
            TR << "Reading header: " << std::endl;
            TR << line << std::endl;
            // do nothing with the header
            TR << "Reading descriptors: " << std::endl;
            continue;
        }
        core::Real val;
        line_stream >> val;
        TR << val << std::endl;
        zd.push_back(val);
    }
    // todo: change TR when class has its own file.
    TR << "Read " << zd.size() << " zernike descriptors" << std::endl;
    return zd;
}

void
ZernikeDescriptorCalculator::invariants_from_vector_to_file(std::vector<std::vector<core::Real>> atom_coords,
    std::vector<core::Real> ljs, std::string const & filename){
    std::ofstream zd_file;
    zd_file.open(filename);

    // header
    // time_t now = time(0);
    // tm *ltm = localtime(&now);
    zd_file << "Zernike invariants | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << " | Surface type: " << voxelgrid_.surface_type() <<
    " | Probe radius: " <<  std::setprecision(3) << voxelgrid_.probe_radius() << " | Shell thickness: " << std::setprecision(3) << voxelgrid_.shell_thickness() << std::endl;

    // descriptors
    std::vector<double> invariants = this->invariants_from_vector(atom_coords, ljs);
    for (double descriptor : invariants) {
        zd_file << descriptor << std::endl;
    }
    zd_file.close();
}

void
ZernikeDescriptorCalculator::invariants_from_pose_to_file(core::pose::Pose const & pose, std::string const & filename){
    std::ofstream zd_file;
    zd_file.open(filename);

    // header
    time_t now = time(0);
    tm *ltm = localtime(&now);
    zd_file << "Invariants for " << pose.pdb_info()->name() << " | date: " << ltm->tm_mday << "-" << 1 + ltm->tm_mon << "-" << 1900 + ltm->tm_year;
    zd_file << " | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << " | Surface type: " << voxelgrid_.surface_type() <<
    " | Probe radius: " <<  std::setprecision(3) << voxelgrid_.probe_radius() << " | Shell thickness: " << std::setprecision(3) << voxelgrid_.shell_thickness() << std::endl;

    // descriptors
    std::vector<double> invariants = this->invariants_from_pose(pose);
    for (double descriptor : invariants) {
        zd_file << descriptor << std::endl;
    }
    zd_file.close();
}

void
ZernikeDescriptorCalculator::invariants_from_oriented_pose_to_file(core::pose::Pose const & pose, std::string const & filename){
    std::ofstream zd_file;
    zd_file.open(filename);

    // header
    time_t now = time(0);
    tm *ltm = localtime(&now);
    zd_file << "Invariants for " << pose.pdb_info()->name() << " | date: " << ltm->tm_mday << "-" << 1 + ltm->tm_mon << "-" << 1900 + ltm->tm_year;
    zd_file << " | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << " | Surface type: " << voxelgrid_.surface_type() <<
    " | Probe radius: " <<  std::setprecision(3) << voxelgrid_.probe_radius() << " | Shell thickness: " << std::setprecision(3) << voxelgrid_.shell_thickness() << std::endl;

    // descriptors
    std::vector<double> invariants = this->invariants_slice_from_pose_oriented(pose);
    for (double descriptor : invariants) {
        zd_file << descriptor << std::endl;
    }
    zd_file.close();
}

void
ZernikeDescriptorCalculator::invariants_2D_from_pose_to_file(core::pose::Pose const & pose, std::string const & filename){
    std::ofstream zd_file;
    zd_file.open(filename);
		time_t now = time(0);
		tm *ltm = localtime(&now);
    voxelgrid_.setup_2D_descriptor(voxelgrid_.actual_grid_size(), order_);
    for (int slice=1; slice <= num_slices_; slice++ ) {
/// OBS!!!!
//     		voxelgrid_ = VoxelGrid(voxelgrid_.grid_size(), voxelgrid_.surface_type(),voxelgrid_.probe_radius(), voxelgrid_.shell_thickness());
        // descriptors
        std::vector<double> invariants = this->invariants_slice_from_pose(pose);

        zd_file << "Invariants for " << pose.pdb_info()->name() << " | date: " << ltm->tm_mday << "-" << 1 + ltm->tm_mon << "-" << 1900 + ltm->tm_year;
        zd_file << " | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << " | Surface type: " << voxelgrid_.surface_type() <<
        " | Probe radius: " <<  std::setprecision(3) << voxelgrid_.probe_radius() << " | Shell thickness: " << std::setprecision(3) << voxelgrid_.shell_thickness() <<
        " | Slice surface vector 1: " << voxelgrid_.slice_surface_vector1()(1) <<  "," << voxelgrid_.slice_surface_vector1()(2) << "," << voxelgrid_.slice_surface_vector1()(3) <<
        " | Slice surface vector 2: " << voxelgrid_.slice_surface_vector2()(1) <<  "," << voxelgrid_.slice_surface_vector2()(2) << "," << voxelgrid_.slice_surface_vector2()(3) <<
        " | Num surface pixels: " << voxelgrid_.num_surface_pixels() << 
	" | Max surface dimension: " << voxelgrid_.max_surface_dim() << std::endl;
               for (double descriptor : invariants) {
            zd_file << descriptor << std::endl;
        }
    }
   zd_file.close();
	utility::file::gzip( filename, true );
}

void
ZernikeDescriptorCalculator::invariants_2D_3D_from_pose_to_file(core::pose::Pose const & pose, std::string const & filename){
	  using namespace basic::options;
    using namespace basic::options::OptionKeys;
    std::string type = option[zernike_descriptor::zernike_transform_type].value();

    std::ofstream zd_file;
    zd_file.open(filename);
		time_t now = time(0);
		tm *ltm = localtime(&now);
    voxelgrid_.setup_2D_descriptor(voxelgrid_.actual_grid_size(), order_);
    for (int slice=1; slice <= num_slices_; slice++ ) {
/// OBS!!!!
//     		voxelgrid_ = VoxelGrid(voxelgrid_.grid_size(), voxelgrid_.surface_type(),voxelgrid_.probe_radius(), voxelgrid_.shell_thickness());
        // descriptors
        //std::vector<double> invariants = this->invariants_3D_slice_from_pose(pose);
        std::vector<double> invariants = this->invariants_slice_from_pose(pose);

        zd_file << "Invariants for " << pose.pdb_info()->name() << " | date: " << ltm->tm_mday << "-" << 1 + ltm->tm_mon << "-" << 1900 + ltm->tm_year;
        zd_file << " | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << " | Surface type: " << voxelgrid_.surface_type() <<
        " | Probe radius: " <<  std::setprecision(3) << voxelgrid_.probe_radius() << " | Shell thickness: " << std::setprecision(3) << voxelgrid_.shell_thickness() <<
        " | Slice surface vector 1: " << voxelgrid_.slice_surface_vector1()(1) <<  "," << voxelgrid_.slice_surface_vector1()(2) << "," << voxelgrid_.slice_surface_vector1()(3) <<
        " | Slice surface vector 2: " << voxelgrid_.slice_surface_vector2()(1) <<  "," << voxelgrid_.slice_surface_vector2()(2) << "," << voxelgrid_.slice_surface_vector2()(3) <<
        " | Num surface pixels: " << voxelgrid_.num_surface_pixels() << " | Max surface dimension: " << voxelgrid_.max_surface_dim() << " | " << type  << std::endl;
               for (double descriptor : invariants) {
            zd_file << descriptor << std::endl;
        }
    }
   zd_file.close();
	utility::file::gzip( filename, true );
}


void
ZernikeDescriptorCalculator::invariants_2D_from_pose_and_slice(core::pose::Pose const & pose, numeric::xyzVector<Real> surface_vector1, numeric::xyzVector<Real> surface_vector2, std::string const & filename){
    std::ofstream zd_file;
    zd_file.open(filename);
        // header
        time_t now = time(0);
        tm *ltm = localtime(&now);
  //     auto start = std::chrono::high_resolution_clock::now();
       // descriptors
        std::vector<double> invariants = this->invariants_slice_from_pose_and_surface_vectors(pose, surface_vector1, surface_vector2);
 std::cout << "filename: " << filename << std::endl;
        zd_file << "Invariants for " << pose.pdb_info()->name() << " | date: " << ltm->tm_mday << "-" << 1 + ltm->tm_mon << "-" << 1900 + ltm->tm_year;
        zd_file << " | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << " | Surface type: " << voxelgrid_.surface_type() <<
        " | Probe radius: " <<  std::setprecision(3) << voxelgrid_.probe_radius() << " | Shell thickness: " << std::setprecision(3) << voxelgrid_.shell_thickness() <<
        " | Slice surface vector 1: " << voxelgrid_.slice_surface_vector1()(1) <<  "," << voxelgrid_.slice_surface_vector1()(2) << "," << voxelgrid_.slice_surface_vector1()(3) <<
        " | Slice surface vector 2: " << voxelgrid_.slice_surface_vector2()(1) <<  "," << voxelgrid_.slice_surface_vector2()(2) << "," << voxelgrid_.slice_surface_vector2()(3)  <<
        " | Num surface pixels: " << voxelgrid_.num_surface_pixels() << 
	" | Max surface dimension: " << voxelgrid_.max_surface_dim() << std::endl;
   //     auto finish = std::chrono::high_resolution_clock::now();
   //     std::chrono::duration<double> elapsed = finish - start;
   //     TR.Debug << "Time for ZD2: " << elapsed.count() << std::endl;
        for (double descriptor : invariants) {
            zd_file << descriptor << std::endl;
        }
   zd_file.close();
}

void
ZernikeDescriptorCalculator::invariants_2D_3D_from_pose_and_slice(core::pose::Pose const & pose, numeric::xyzVector<Real> surface_vector1, numeric::xyzVector<Real> surface_vector2, std::string const & filename){
    std::ofstream zd_file;
    zd_file.open(filename);
        // header
        time_t now = time(0);
        tm *ltm = localtime(&now);
  //     auto start = std::chrono::high_resolution_clock::now();
       // descriptors
        std::vector<double> invariants = this->invariants_slice_from_pose_and_surface_vectors(pose, surface_vector1, surface_vector2);
 std::cout << "filename: " << filename << std::endl;
        zd_file << "Invariants for " << pose.pdb_info()->name() << " | date: " << ltm->tm_mday << "-" << 1 + ltm->tm_mon << "-" << 1900 + ltm->tm_year;
        zd_file << " | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << " | Surface type: " << voxelgrid_.surface_type() <<
        " | Probe radius: " <<  std::setprecision(3) << voxelgrid_.probe_radius() << " | Shell thickness: " << std::setprecision(3) << voxelgrid_.shell_thickness() <<
        " | Slice surface vector 1: " << voxelgrid_.slice_surface_vector1()(1) <<  "," << voxelgrid_.slice_surface_vector1()(2) << "," << voxelgrid_.slice_surface_vector1()(3) <<
        " | Slice surface vector 2: " << voxelgrid_.slice_surface_vector2()(1) <<  "," << voxelgrid_.slice_surface_vector2()(2) << "," << voxelgrid_.slice_surface_vector2()(3)  <<
        " | Num surface pixels: " << voxelgrid_.num_surface_pixels() << 
	" | Max surface dimension: " << voxelgrid_.max_surface_dim() << std::endl;
   //     auto finish = std::chrono::high_resolution_clock::now();
   //     std::chrono::duration<double> elapsed = finish - start;
   //     TR.Debug << "Time for ZD2: " << elapsed.count() << std::endl;
        for (double descriptor : invariants) {
            zd_file << descriptor << std::endl;
        }
   zd_file.close();
}


// void
// ZernikeDescriptorCalculator::invariants_from_grid_file_to_file(std::string input_filename, std::string const & output_filename){
//     using namespace basic::options;
//     using namespace basic::options::OptionKeys;
//
//     std::ofstream zd_file;
//     zd_file.open(output_filename);
//     // header
//     time_t now = time(0);
//     tm *ltm = localtime(&now);
// //    auto start = std::chrono::high_resolution_clock::now();
//     // descriptors
//     std::vector<double> invariants = this->invariants_from_grid_file(input_filename);
//     std::string type = option[zernike_descriptor::zernike_transform_type].value();
//
//     zd_file << "Invariants for " << input_filename << " | date: " << ltm->tm_mday << "-" << 1 + ltm->tm_mon << "-" << 1900 + ltm->tm_year;
//     zd_file << " | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << " | Surface type: " << voxelgrid_.surface_type() <<
//     " | Probe radius: " <<  std::setprecision(3) << voxelgrid_.probe_radius() << " | Shell thickness: " << std::setprecision(3) << voxelgrid_.shell_thickness()
//     << " | " << type << std::endl;
// //    auto finish = std::chrono::high_resolution_clock::now();
// //    std::chrono::duration<double> elapsed = finish - start;
// //    std::cout << "Time for ZD: " << elapsed.count()  << std::endl;
//     for (double descriptor : invariants) {
//         zd_file << descriptor << std::endl;
//     }
//    zd_file.close();
// }

void
ZernikeDescriptorCalculator::invariants_from_grid_outline_file_to_file(std::string input_filename, std::string const & output_filename){
    std::ofstream zd_file;
    zd_file.open(output_filename);
    // header
    time_t now = time(0);
    tm *ltm = localtime(&now);
//    auto start = std::chrono::high_resolution_clock::now();
    // descriptors
    std::vector<double> invariants = this->invariants_from_outline_grid_file(input_filename);

    zd_file << "Invariants for " << input_filename << " | date: " << ltm->tm_mday << "-" << 1 + ltm->tm_mon << "-" << 1900 + ltm->tm_year;
    zd_file << " | Parameters: Order:" << order_ << " | Grid dimension: " << voxelgrid_.grid_size() << std::endl;
//    auto finish = std::chrono::high_resolution_clock::now();
//    std::chrono::duration<double> elapsed = finish - start;
//    std::cout << "Time for ZD: " << elapsed.count()  << std::endl;
    for (double descriptor : invariants) {
        zd_file << descriptor << std::endl;
    }
   zd_file.close();
}

std::vector<core::Real>
ZernikeDescriptorCalculator::get_3d_descriptors() {
  return voxelgrid_.get_zernike_descriptor().GetInvariants();
}

std::vector<core::Real>
ZernikeDescriptorCalculator::get_2d_descriptors() {
  return voxelgrid_.get_zernike_2D_descriptor().GetInvariants();
}

std::vector<std::complex<double>>
ZernikeDescriptorCalculator::get_3d_moments() {
	return voxelgrid_.get_zernike_descriptor().GetMoments();
}


std::vector<std::complex<double>>
ZernikeDescriptorCalculator::get_2d_moments() {
	return voxelgrid_.get_zernike_2D_descriptor().GetMoments();
}

VoxelGrid &
ZernikeDescriptorCalculator::get_grid() {
    return voxelgrid_;
}

/// constructor with parameters - also the default constructor
ZernikeDescriptorEnergy::ZernikeDescriptorEnergy() : parent(utility::pointer::make_shared<ZernikeDescriptorEnergyCreator>()) {
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    // todo make these options the can be parsed in through the energymethodsoptions
    int order = option[zernike_descriptor::order].value();
    int size = option[zernike_descriptor::grid_size].value();
    std::string surface_type = option[zernike_descriptor::surface_type].value();
    core::Real probe = option[zernike_descriptor::probe_radius].value();
    core::Real shell = option[zernike_descriptor::shell_thickness].value();
    std::string type = option[zernike_descriptor::zernike_transform_type].value();
    int num_slices = option[zernike_descriptor::num_2d_slices].value();

    zd_calc_ = ZernikeDescriptorCalculator(order, size, surface_type, probe, shell, num_slices);

    TR << "Parameters" << std::endl;
    TR << "- Order: " << std::setw(5+15) << order << std::endl;
    TR << "- Grid dimension: " << std::setw(5+6) << size << std::endl;
    TR << "- Grid surface type: " << std::setw(5+3) << surface_type << std::endl;
    TR << "- Grid probe radius: " << std::setw(5+3) << probe  << std::endl;
    TR << "- Grid shell thickness: " << std::setw(5) << shell << std::endl;
    TR << "- Transform type: " << std::setw(5) << type << std::endl;
    TR << "- Number of 2D slices: " << std::setw(5) << num_slices << std::endl;

    if (option[zernike_descriptor::zernike_descriptor_file].user()) {
        std::string file_name = option[zernike_descriptor::zernike_descriptor_file].value();
        zd_reference_ = ZernikeDescriptorData(zd_calc_.invariants_from_file(file_name));
    }
    //  todo: to get this to work  in the future - move code from core.3 to core.6
//    else if (option[zernike_descriptor::zernike_descriptor_pose].user()) {
//        std::string const & pose_name = option[zernike_descriptor::zernike_descriptor_pose].name();
//        core::pose::PoseOP pose_reference = core::import_pose::pose_from_file( pose_name);
//        zd_reference_ = ZernikeDescriptorData(zd_calc_.invariants_from_pose(*pose_reference));
//    }
    else
        utility_exit_with_message("supply either pose or file as a reference");

}

/// @details This must return a fresh instance of the ZernikeDescriptorEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ZernikeDescriptorEnergyCreator::create_energy_method(
        methods::EnergyMethodOptions const &) const {
    return utility::pointer::make_shared<ZernikeDescriptorEnergy>();
}

ScoreTypes
ZernikeDescriptorEnergyCreator::score_types_for_method() const {
    ScoreTypes sts;
    sts.push_back(zernike_descriptor);
    return sts;
}

void
ZernikeDescriptorEnergy::finalize_total_energy(
        pose::Pose &pose,
        ScoreFunction const &,
        EnergyMap &totals) const {
    core::Real result = compute(pose, false); // before: const_cast< const pose::Pose &>( pose )
    totals[zernike_descriptor] = result;
}

EnergyMethodOP
ZernikeDescriptorEnergy::clone() const {
    return utility::pointer::make_shared< ZernikeDescriptorEnergy >( *this );
    }

void
ZernikeDescriptorEnergy::setup_for_derivatives(pose::Pose &, ScoreFunction const &) const {}


void
ZernikeDescriptorEnergy::setup_for_scoring(pose::Pose & pose, ScoreFunction const &) const{
    setup_zernike_descriptor_reference( pose );
}

void
ZernikeDescriptorEnergy::eval_atom_derivative(
    id::AtomID const & aid,
    pose::Pose const & pose,
    kinematics::DomainMap const &,
    ScoreFunction const &,
    EnergyMap const &,
    Vector &,
    Vector & ) const {
    conformation::Residue const & rsd = pose.residue( aid.rsd() );
    if ( ! rsd.is_protein() )
        return;
    Vector f1(0.0), f2(0.0);
}

void
ZernikeDescriptorEnergy::setup_zernike_descriptor_reference( core::pose::Pose & pose ) const {
    pose.data().set( core::pose::datacache::CacheableDataType::ZERNIKE_DESCRIPTOR_DATA, ZernikeDescriptorDataOP(new ZernikeDescriptorData(zd_reference_)));
}

core::Real
ZernikeDescriptorEnergy::compute( core::pose::Pose const & pose, bool report ) const {

    /// TODO copies ? we want references!
    std::vector<double> const zd_current = zd_calc_.invariants_from_pose(pose);
    std::vector<double> const zd_reference = zernike_reference_data_pose(pose);

    if (report)
        report_on_current_invariants(zd_current);

    if (zd_current.size() != zd_reference.size() )
        utility_exit_with_message("Zernike Descriptors have different sizes. "
        "Did you parse the zernike_file commandline flag? "
        "If not there's no reference shape to score against.");

    // Calculate the Zernike score as a L2 norm between the computed descriptors and the reference descriptor.
    core::Real L2_norm = 0;
    for (core::Size i=0; i != zd_current.size(); ++i )
        L2_norm += (zd_current[i] - zd_reference[i] ) * (zd_current[i] - zd_reference[i] );
    core::Real sqrt_L2_norm = sqrt(L2_norm);

    if (report)
        report_on_current_L2_norm(L2_norm);

    // Debug
    if ( sqrt_L2_norm > 100000 ) {
        for (core::Size i=0; i<zd_current.size(); ++i ) {
            std::cout << i << " : " << zd_current[i] << " vs " << zd_reference[i] << std::endl;
        }
        utility_exit_with_message("Zernike Descriptors weird!!!");
    }

    // todo: can delete this later on.
    TR << sqrt_L2_norm  << std::endl;

    return sqrt_L2_norm;
}

core::Size
ZernikeDescriptorEnergy::version() const
{
return 1; // Initial versioning
}

std::vector<double> const &
ZernikeDescriptorEnergy::zernike_reference_data_pose(core::pose::Pose const & pose) const {
    // no need to check  if it is in pose - because an error will be thrown internally
    // if not set in pose then set it now.

    return utility::pointer::static_pointer_cast<const core::scoring::shape::ZernikeDescriptorData>(
            pose.data().get_const_ptr(core::pose::datacache::CacheableDataType::ZERNIKE_DESCRIPTOR_DATA))->invariants();
}

} // shape
} // scoring
} // core
