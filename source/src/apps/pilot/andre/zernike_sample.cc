
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


static basic::Tracer TR( "apps.pilot.zernike_align" );
    OPT_KEY( String, grid_file_pdb)
//    OPT_KEY( String, surface_vector1_string)
    OPT_KEY( String, surface_vectors)
//    OPT_KEY( String, surface_vector2_string)
    OPT_KEY( String, grid_file_sampled_shape)

float
zeal_score( std::vector<std::complex<double>> moments1, std::vector<std::complex<double>> moments2 ) {

     core::Size moment_size1( moments1.size() );
     core::Size moment_size2( moments1.size() );
     if ( moment_size1 != moment_size2 ) {
         TR << "Moments have different size!!! Exit" << std::endl;
          std::exit(0);
    }
     std::complex<double> prod(0);
     std::complex<double> prodA(0);
     std::complex<double> prodB(0);
     for ( core::Size index=0; index < moment_size1; index++ ) {
	prod += moments1[index]*std::conj(moments2[index]);
	prodA += moments1[index]*std::conj(moments1[index]);
	prodB += moments2[index]*std::conj(moments2[index]);
//        TR << "Moments1: " << moments1[index] << " " << moments2[index] << std::endl;
     }
     float real_prod( prod.real() );
     float real_root_prodA ( std::sqrt( prodA.real() ) );
     float real_root_prodB ( std::sqrt( prodB.real() ) );
     float zeal_score = real_prod/(real_root_prodA*real_root_prodB);
    return zeal_score;
}

std::vector< std::vector<int> >
get_nm(int n) {
    std::vector<std::vector<int> > nm;
    for (int ni = 0; ni <=n; ni++) {
        std::vector<int> mv;
//        for (int mi = -ni; mi <= ni; mi++) {
        for (int mi = 0; mi <= ni; mi++) {
            int absm = abs(mi);
            if ((ni - absm) % 2 == 0) {
                mv.push_back(mi);
            }
        }
        nm.push_back(mv);
    }
    return nm;
}

numeric::xyzVector<core::Real>
get_center_of_pose(core::pose::Pose & pose) {

    core::Real total_x = 0;
    core::Real total_y = 0;
    core::Real total_z = 0;
    int n_atoms = 0;
    for (core::pose::Pose::const_iterator res = pose.begin(); res != pose.end(); ++res) {
        for (core::Size atom_index = 1; atom_index <= res->natoms(); ++atom_index) {
            // skip if atom is Hydrogen or if residue is virtual.
            if (!res->atom_is_hydrogen(atom_index && res->aa() != core::chemical::aa_vrt)) {
                core::Real x = res->xyz(atom_index).x();
                core::Real y = res->xyz(atom_index).y();
                core::Real z = res->xyz(atom_index).z();
		if (atom_index == 1 ) TR <<  "res " << *res << " " << x << " " << y << " " << z << std::endl;
                total_x += x;
                total_y += y;
                total_z += z;
                n_atoms += 1;
            }
        }
    }

    // calculate the geometric center of the pose
    numeric::xyzVector<core::Real> geometric_center(total_x / n_atoms, total_y / n_atoms,
                                                                    total_z / n_atoms);
   return geometric_center;
}

core::pose::PoseOP
rotate_pose( core::pose::Pose const & pose, numeric::xyzVector<core::Real> surface_vector1, numeric::xyzVector<core::Real> surface_vector2 ) {
    numeric::xyzVector<core::Real> rot_x = surface_vector1.normalize();
    numeric::xyzVector<core::Real> rot_y = surface_vector2.normalize();
    numeric::xyzVector<core::Real> rot_z = cross(surface_vector1, surface_vector2 ).normalize();
    numeric::xyzMatrix<core::Real> rot_mat = numeric::xyzMatrix<core::Real>::rows(rot_x, rot_y, rot_z);
    core::pose::PoseOP pose_rotated(pose.clone());
    core::Size nres (pose_rotated->total_residue() );
    addVirtualResAsRoot(*pose_rotated);
    core::conformation::ResidueOP VRT_base = core::conformation::ResidueFactory::create_residue(*core::pose::get_restype_for_pose( *pose_rotated, "VRT" ));
    pose_rotated->append_residue_by_jump(*VRT_base, nres+1);
    core::kinematics::FoldTree newF( pose_rotated->fold_tree() );
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


int
main( int argc, char ** argv ) {

    using namespace core;
    using namespace utility;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

//    NEW_OPT( grid_file_pdb, "grid file for pose to align", "");
//    NEW_OPT( surface_vector1_string, "surface vector1 for slice", "");
    NEW_OPT( surface_vectors, "surface vectors file", "");
//    NEW_OPT( surface_vector2_string, "surface vector1 for slice", "");
    NEW_OPT( grid_file_sampled_shape, "grid file for sampled shape file to align to", "");

    devel::init( argc, argv );

    // options for zernike calculations
    int order = option[zernike_descriptor::order].value();
    int size = option[zernike_descriptor::grid_size].value();
    std::string surface_type = option[zernike_descriptor::surface_type].value();
    core::Real probe = option[zernike_descriptor::probe_radius].value();
    core::Real shell = option[zernike_descriptor::shell_thickness].value();
    std::string zernike_transform_type = option[zernike_descriptor::zernike_transform_type].value();
//    bool outline_transform = option[zernike_descriptor::transform_outline].value();
//    int num_slices = option[zernike_descriptor::num_2d_slices].value();

    vector1<std::string> files;
    files = option[in::file::s]();
    std::string file = files[1];
    TR << file << std::endl;
    core::pose::PoseOP pose = core::import_pose::pose_from_file(file);
    std::string output = pose->pdb_info()->name();
    std::string pdb_prefix(utility::string_split(utility::string_split(output, '/').back(), '.').front());
    TR << "prefix: " << pdb_prefix << std::endl;

   // sets the pose to centroid if not the fullatom flag is set
   if (!basic::options::option[in::file::fullatom].user()) {
        protocols::simple_moves::SwitchResidueTypeSetMover to_centroid(core::chemical::CENTROID);
        to_centroid.apply(*pose);
   }

    if ( zernike_transform_type == "2D" ) {
         TR << "2D zernike transform" << std::endl;

        core::scoring::shape::ZernikeDescriptorCalculator zc_calc = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell,1);
	//zc_calc.invariants_slice_from_pose_to_file(*pose );
	zc_calc.invariants_2D_from_pose_to_file( *pose, basic::options::option[out::file::o].value() + "/sample." + pdb_prefix + ".inv");
    }

    if ( zernike_transform_type == "2D_3D" ) {
         TR << "2D_3D zernike transform" << std::endl;

        core::scoring::shape::ZernikeDescriptorCalculator zc_calc = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell,1);
	//zc_calc.invariants_slice_from_pose_to_file(*pose );
	zc_calc.invariants_2D_3D_from_pose_to_file( *pose, basic::options::option[out::file::o].value() + "/sample." + pdb_prefix + ".inv");
    }

  return 0;
}

