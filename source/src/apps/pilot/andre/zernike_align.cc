
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
     core::Size moment_size2( moments2.size() );
     if ( moment_size1 != moment_size2 ) {
         TR << "Moments have different size!!! Exit " << moment_size1 << " vs " << moment_size2 << std::endl;
          std::exit(0);
    }
     std::complex<double> prod(0);
     std::complex<double> prodA(0);
     std::complex<double> prodB(0);
     for ( core::Size index=0; index < moment_size1; index++ ) {
	prod += moments1[index]*std::conj(moments2[index]);
	prodA += moments1[index]*std::conj(moments1[index]);
	prodB += moments2[index]*std::conj(moments2[index]);
      //  TR << "Moments1: " << moments1[index] << " " << moments2[index] << std::endl;
     }
     float real_prod( prod.real() );
     float real_root_prodA ( std::sqrt( prodA.real() ) );
     float real_root_prodB ( std::sqrt( prodB.real() ) );
     float zeal_score = real_prod/(real_root_prodA*real_root_prodB);
    return zeal_score;
}

float
zdiff_score( std::vector<core::Real> descriptor1, std::vector<core::Real> descriptor2 ) {

     core::Size descriptor_size1( descriptor1.size() );
     core::Size descriptor_size2( descriptor2.size() );
     if ( descriptor_size1 != descriptor_size2 ) {
         TR << "Descriptors have different size!!! Exit" << std::endl;
          std::exit(0);
    }
     double diff_sum(0);
     for ( core::Size index=0; index < descriptor_size1; index++ ) {
	diff_sum += (descriptor1[index] - descriptor2[index])*(descriptor1[index] - descriptor2[index]);
     }
     float zdiff ( std::sqrt( diff_sum ) );
    return zdiff;
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
    numeric::xyzVector<core::Real> rot_x_new = cross(rot_y, rot_z ).normalize();
    rot_x = rot_x_new;
    numeric::xyzMatrix<core::Real> rot_mat = numeric::xyzMatrix<core::Real>::rows(rot_x, rot_y, rot_z);
    TR << "Determinant rotation is: " << rot_mat.det() << std::endl;
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

numeric::xyzMatrix<core::Real>
get_rotation_matrix_z_rad( float angle_rad ) {

    numeric::xyzVector<core::Real> x_vec( cos(angle_rad), -sin(angle_rad), 0);
    numeric::xyzVector<core::Real> y_vec( sin(angle_rad), cos(angle_rad), 0);
    numeric::xyzVector<core::Real> z_vec( 0, 0, 1);

    numeric::xyzMatrix<core::Real> rot_mat_align = numeric::xyzMatrix<core::Real>::rows(x_vec, y_vec, z_vec);

    return rot_mat_align;
}

numeric::xyzMatrix<core::Real>
get_rotation_matrix_z( float angle) {
    float angle_rad = float(angle/180.0*3.1415);
    return get_rotation_matrix_z_rad(angle_rad);
}

numeric::xyzMatrix<core::Real>
get_rotation_axis_angle_rad( numeric::xyzVector<core::Real> axis, float angle_rad ) {

    core::Real xx_ = axis(1)*axis(1)*(1-cos(angle_rad)) + cos(angle_rad);
    core::Real xy_ = axis(1)*axis(2)*(1-cos(angle_rad)) - axis(3)*sin(angle_rad);
    core::Real xz_ = axis(1)*axis(3)*(1-cos(angle_rad)) + axis(2)*sin(angle_rad);

    core::Real yx_ = axis(1)*axis(2)*(1-cos(angle_rad)) + axis(3)*sin(angle_rad);
    core::Real yy_ = axis(2)*axis(2)*(1-cos(angle_rad)) + cos(angle_rad);
    core::Real yz_ = axis(2)*axis(3)*(1-cos(angle_rad)) - axis(1)*sin(angle_rad);

    core::Real zx_ = axis(1)*axis(3)*(1-cos(angle_rad)) - axis(2)*sin(angle_rad);
    core::Real zy_ = axis(2)*axis(3)*(1-cos(angle_rad)) + axis(1)*sin(angle_rad);
    core::Real zz_ = axis(3)*axis(3)*(1-cos(angle_rad)) + cos(angle_rad);
	    
    numeric::xyzVector<core::Real> x_vec( xx_,xy_,xz_);
    numeric::xyzVector<core::Real> y_vec( yx_, yy_, yz_);
    numeric::xyzVector<core::Real> z_vec( zx_,zy_, zz_);

    numeric::xyzMatrix<core::Real> rot_mat_align = numeric::xyzMatrix<core::Real>::rows(x_vec, y_vec, z_vec);

    return rot_mat_align;
}

numeric::xyzMatrix<core::Real>
get_rotation_axis_angle( numeric::xyzVector<core::Real> axis, float angle) {
    float angle_rad = float(angle/180.0*3.1415);
    return get_rotation_matrix_z_rad(angle_rad);
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
    bool debug = option[zernike_descriptor::debug].value();

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

    std::string surface_filename = option[surface_vectors];
    std::ifstream infile(surface_filename.c_str());
    std::string line;

    numeric::xyzVector<Real> slice_vector1(0,0,0);
    numeric::xyzVector<Real> slice_vector2(0,0,0);
    while ( getline(infile, line) ) {
       TR << line << std::endl;
    	utility::vector1< std::string > identifier = utility::string_split(line, ':');
    	utility::vector1< std::string > items = utility::string_split(identifier[2], ',');
    	Real v1(utility::string2float(items[1]));
    	Real v2(utility::string2float(items[2]));
    	Real v3(utility::string2float(items[3]));
        if ( identifier[1] == "surface_vector1" ) {
    		numeric::xyzVector<Real> vec1(v1,v2,v3);
       		slice_vector1 = vec1;
        }
       if ( identifier[1] == "surface_vector2" ) {
                numeric::xyzVector<Real> vec2(v1,v2,v3);
                slice_vector2 = vec2;
        }
    }

     std::string grid_file = option[grid_file_sampled_shape];
     core::scoring::shape::ZernikeDescriptorCalculator zc_calc_grid = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
     std::vector< std::vector<int> > slice_for_ZT(zc_calc_grid.read_slice_grid(grid_file));
     zc_calc_grid.invariants_from_grid(slice_for_ZT);
     core::scoring::shape::ZernikeDescriptorCalculator zc_grid_output = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
     zc_grid_output.invariants_from_grid_file_to_file(grid_file, basic::options::option[out::file::o].value() + "/" + grid_file + ".grid.inv");

    std::vector<std::complex<double>> moments_pose;
    std::vector<std::complex<double>>   moments_grid;
    std::vector<core::Real> invariants_pose; 
    std::vector<core::Real> invariants_grid;

    if ( zernike_transform_type == "2D" || zernike_transform_type == "2D_3D" ) {
        core::scoring::shape::ZernikeDescriptorCalculator zc_calc_pose = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
    	if ( zernike_transform_type == "2D" ) {
         	TR << "2D zernike transform" << std::endl;
         	zc_calc_pose.invariants_2D_from_pose_and_slice(*pose, slice_vector1, slice_vector2, basic::options::option[out::file::o].value() + "/" + pdb_prefix + ".align.inv" );
     		moments_pose = zc_calc_pose.get_2d_moments();
           	moments_grid = zc_calc_grid.get_2d_moments();
      	        invariants_pose = zc_calc_pose.get_2d_descriptors();
                invariants_grid = zc_calc_grid.get_2d_descriptors();
	}
	else if ( zernike_transform_type == "2D_3D" ) {
         TR << "2D_3D zernike transform" << std::endl;
         	zc_calc_pose.invariants_2D_3D_from_pose_and_slice(*pose, slice_vector1, slice_vector2, basic::options::option[out::file::o].value() + "/" + pdb_prefix + ".align.inv" );
		moments_pose = zc_calc_pose.get_3d_moments();
           	moments_grid = zc_calc_grid.get_3d_moments();
                invariants_pose = zc_calc_pose.get_3d_descriptors();
                invariants_grid = zc_calc_grid.get_3d_descriptors();
	}

//     std::vector<std::complex<double>>   moments2 = zc_calc2.get_2d_moments();
//     float zeal_score2 = zeal_score(moments1, moments2 );
//     TR << "Zeal score: " << zeal_score2 << std::endl;

    
//     std::vector<std::complex<double>>   moments_grid = zc_calc_grid.get_2d_moments();
     float zeal_score_pose_grid = zeal_score(moments_pose, moments_grid );
     TR << "Zeal score: " << zeal_score_pose_grid << std::endl;
      float zdiff = zdiff_score( invariants_pose, invariants_grid );
     TR << "Zeal L2 norm: " << zdiff << std::endl;

////////////////////////
/*       core::Real deg_rot_90  = float(90)/180.0*3.1415;
       std::vector< std::vector<int> > rotate_grid_90;
       core::Size grid_size(slice_for_ZT[0].size());
       for (core::Size i = 0; i< grid_size; i++ ) {
           rotate_grid_90.push_back( std::vector<int>(grid_size) );
       }

        for (core::Size x=0; x< grid_size; x++ ) {
	    for (core::Size y=0; y< grid_size; y++ ) {
                if (slice_for_ZT[x][y] == 1) {
                    int translated_coords_x = x - grid_size/2;
                    int translated_coords_y = y - grid_size/2;
                    int rotated_x = int(translated_coords_x*cos(deg_rot_90) - translated_coords_y*sin(deg_rot_90));
                    int rotated_y = int(translated_coords_x*sin(deg_rot_90) + translated_coords_y*cos(deg_rot_90));
                    rotated_x += grid_size/2;
                    rotated_y += grid_size/2;
                   rotate_grid_90[rotated_x][rotated_y] = slice_for_ZT[x][y];
               }
           }
       }
       core::scoring::shape::ZernikeDescriptorCalculator zc_calc_90 = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
       zc_calc_90.invariants_from_grid(rotate_grid_90);
       std::vector<std::complex<double>>   moments_90 = zc_calc_90.get_2d_moments();

*/
///////////////////////
    numeric::xyzVector<Real> slice_x = slice_vector1.normalize();
    numeric::xyzVector<Real> slice_y = slice_vector2.normalize();
    numeric::xyzVector<Real> slice_z = cross(slice_x, slice_y ).normalize();
 //   numeric::xyzVector<Real> slice_x_new = cross(slice_z, slice_y ).normalize();
//    slice_x = slice_x_new;
    TR << "slice_x: " << slice_x(1) << " " << slice_x(2) << " " << slice_x(3) << std::endl;
    numeric::xyzMatrix<core::Real> rot_mat_to_plane = numeric::xyzMatrix<core::Real>::rows(slice_x, slice_y, slice_z);
    TR << "MAT: " << rot_mat_to_plane(1,1) << " " << rot_mat_to_plane(1,2) << " " << rot_mat_to_plane(1,3) << std::endl
	      << "MAT: " << rot_mat_to_plane(2,1) << " " << rot_mat_to_plane(2,2) << " " << rot_mat_to_plane(2,3) << std::endl
	      << "MAT: " << rot_mat_to_plane(3,1) << " " << rot_mat_to_plane(3,2) << " " << rot_mat_to_plane(3,3) << std::endl;

   TR << "Determinant is: " << rot_mat_to_plane.det() << std::endl;
   numeric::xyzVector<Real> rot_slice_vector1 = rot_mat_to_plane*slice_x;
   numeric::xyzVector<Real> rot_slice_vector2 = rot_mat_to_plane*slice_y;
   numeric::xyzVector<Real> rot_slice_vector3 = rot_mat_to_plane*slice_z;

TR << "After rot x: " << rot_slice_vector1(1) << " " << rot_slice_vector1(2) << " " << rot_slice_vector1(3) << std::endl;
TR << "After rot y: " << rot_slice_vector2(1) << " " << rot_slice_vector2(2) << " " << rot_slice_vector2(3) << std::endl;
TR << "After rot z: " << rot_slice_vector3(1) << " " << rot_slice_vector3(2) << " " << rot_slice_vector3(3) << std::endl;

// These are the correct facet unit vectors
//    numeric::xyzVector<Real> facet_x(0.1546028147050133, 0.9401951216827611, -0.30353105747060816);
//    numeric::xyzVector<Real> facet_y(-0.7957753207953651, 0.30058090157287454, 0.5257306919161859);
//    numeric::xyzVector<Real> facet_z(0.585525070768384, 0.16026307988297306, 0.7946547279970735);

// Here I have modified the facet vectors to compensate for a bug within Rosetta in how the virtual residue code
// It is likely a problem with the reference frame constructed for jumps. The corrections are
// facet_x: y-> -y
// facet y: x-> -x, z->-z 
    numeric::xyzVector<Real> facet_x(0.1546028147050133, -0.9401951216827611, -0.30353105747060816);
    numeric::xyzVector<Real> facet_y(0.7957753207953651, 0.30058090157287454, -0.5257306919161859);
    numeric::xyzVector<Real> facet_z(0.585525070768384, 0.16026307988297306, 0.7946547279970735);

    numeric::xyzMatrix<core::Real> rot_mat_to_facet = numeric::xyzMatrix<core::Real>::rows(facet_x, facet_y, facet_z);
   TR << "Facet determinant is: " << rot_mat_to_facet.det() << std::endl;

 numeric::xyzVector<Real> cart_x(1,0,0);
 numeric::xyzVector<Real> cart_y(0,1,0);
 numeric::xyzVector<Real> cart_z(0,0,1);

   numeric::xyzVector<Real> rot_facet1 = rot_mat_to_facet*cart_x;
   numeric::xyzVector<Real> rot_facet2 = rot_mat_to_facet*cart_y;
   numeric::xyzVector<Real> rot_facet3 = rot_mat_to_facet*cart_z;

TR << "After facet rot x: " << rot_facet1(1) << " " << rot_facet1(2) << " " << rot_facet1(3) << std::endl;
TR << "After facet rot y: " << rot_facet2(1) << " " << rot_facet2(2) << " " << rot_facet2(3) << std::endl;
TR << "After facet rot z: " << rot_facet3(1) << " " << rot_facet3(2) << " " << rot_facet3(3) << std::endl;
TR << "facet_x: " << facet_x(1) << " " << facet_x(2) << " " << facet_x(3) << std::endl;
TR << "facet_y: " << facet_y(1) << " " << facet_y(2) << " " << facet_y(3) << std::endl;
TR << "facet_z: " << facet_z(1) << " " << facet_z(2) << " " << facet_z(3) << std::endl;

core::pose::PoseOP rotated_pose = rotate_pose(*pose, slice_vector1, slice_vector2);

TR << rotated_pose->fold_tree() << std::endl;
if (debug) {
    rotated_pose->dump_pdb("poseVRT.pdb");
}

TR << "CENTER: " << zc_calc_pose.grid().slice_average_x() << " " << zc_calc_pose.grid().slice_average_y() << std::endl;
core::Real offset_x = zc_calc_pose.grid().slice_average_x() - zc_calc_pose.grid().actual_grid_size()/2;
core::Real offset_y = zc_calc_pose.grid().slice_average_y() - zc_calc_pose.grid().actual_grid_size()/2;
offset_x *= zc_calc_pose.grid().grid_resolution();
offset_y *= zc_calc_pose.grid().grid_resolution();

TR << "OFFSET Angstrom: " << offset_x << " " << offset_y << std::endl;
numeric::xyzVector<Real> offset_vec(-offset_x,-offset_y,0);
//numeric::xyzVector<Real> offset_vec(0,0,0);

core::conformation::ResidueOP VRT_coord = core::conformation::ResidueFactory::create_residue(*core::pose::get_restype_for_pose( *rotated_pose, "VRT" ));
rotated_pose->replace_residue( rotated_pose->total_residue() -1 , *VRT_coord, false );
numeric::xyzVector< Real > coord = rotated_pose->residue(1).atom(rotated_pose->residue(1).atom_index("CA")).xyz();

if (debug) {
   rotated_pose->dump_pdb("recentered_and_rotated.pdb");
}
    int optimal_angle(0);
    float max_corr(0);
    float max_corr_analytical(0);
    for (core::Size step = 0; step < 360; step=step+2 ) {
	core::Real rot_angle  = float(step)/180.0*3.1415;
	numeric::xyzMatrix<core::Real> rot_mat = get_rotation_matrix_z_rad( rot_angle );
	core::pose::PoseOP pose_test(rotated_pose->clone());
	core::kinematics::Jump j2 = pose_test->jump(2);
	j2.reverse();
	j2.set_rotation(j2.get_rotation()*rot_mat);
	pose_test->set_jump(2, j2);
	//UGLY make sure its get refolded!!!
	coord = pose_test->residue(1).atom(pose_test->residue(1).atom_index("CA")).xyz();
	core::scoring::shape::ZernikeDescriptorCalculator zc_calc_test = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
	zc_calc_test.invariants_from_oriented_pose_to_file(*pose_test, basic::options::option[out::file::o].value() + "/" + pdb_prefix + ".test.inv" );
	 std::vector<std::complex<double>> moments_test;
	if ( zernike_transform_type == "2D" )
		moments_test = zc_calc_test.get_2d_moments();
	else if ( zernike_transform_type == "2D_3D" )
		moments_test = zc_calc_test.get_3d_moments();

	float zeal_score_pose_test = zeal_score(moments_grid, moments_test );

/*       std::vector< std::vector<int> > rotate_grid;
       core::Size grid_size(slice_for_ZT[0].size());
       for (core::Size i = 0; i< grid_size; i++ ) {
           rotate_grid.push_back( std::vector<int>(grid_size) );
       }

        for (core::Size x=0; x< grid_size; x++ ) {
	    for (core::Size y=0; y< grid_size; y++ ) {
                if (slice_for_ZT[x][y] == 1) {
                    int translated_coords_x = x - grid_size/2;
                    int translated_coords_y = y - grid_size/2;
                    int rotated_x = int(translated_coords_x*cos(rot_angle) - translated_coords_y*sin(rot_angle));
                    int rotated_y = int(translated_coords_x*sin(rot_angle) + translated_coords_y*cos(rot_angle));
                    rotated_x += grid_size/2;
                    rotated_y += grid_size/2;
                   // TR << x << " " << y << " " << rotated_x << " " << rotated_y << std::endl;
                   rotate_grid[rotated_x][rotated_y] = slice_for_ZT[x][y];
               }
           }
       }
       core::scoring::shape::ZernikeDescriptorCalculator zc_calc4 = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
 */      
	
	if ( zernike_transform_type == "2D" ) {
		std::vector< std::vector<int> > nm = get_nm(order);
 		std::vector<std::complex<double>> rotated_moments( moments_grid );
       		Size counter=0;
       		for (int ni=0; ni< int(nm.size());ni++ ) {
	       		std::cout << "ni: " << ni << std::endl;
          		for (int j=0; j< int(nm[ni].size());j++ ) {
             		int m = nm[ni][j];
             		std::complex<double> imval(0, rot_angle * double(m) );
             		rotated_moments[counter] = moments_grid[counter]*exp(imval);
//             TR << "ROT: " << rotated_moments[counter] << " " << moments1[counter] << " " << exp(imval) << std::endl;
             		counter++;
           		}
        	}
       		float zeal_score_test = zeal_score(moments_pose, rotated_moments );
    	}
//       zc_calc4.invariants_from_grid(rotate_grid);
//       std::vector<std::complex<double>>   moments4 = zc_calc4.get_2d_moments();
//       float zeal_score4 = zeal_score(moments1, moments4 );
//       if (zeal_score4 > max_corr ) {
//       if (zeal_score_test > max_corr ) {
//           max_corr = zeal_score_test;
//           optimal_angle = step;
//       }
       if (zeal_score_pose_test > max_corr ) {
           max_corr = zeal_score_pose_test;
	   //max_corr_analytical = zeal_score_test;
           optimal_angle = step;
       }

       //TR << "Angle: " << step << " Zeal score: " << zeal_score_test << " pose: " << zeal_score_pose_test << std::endl;
       TR << "Angle: " << step << " pose: " << zeal_score_pose_test << std::endl;
    }
//TR << "Optimal angle: " << " " << optimal_angle << " zeal score: " << max_corr << "Analytical: " << max_corr_analytical << std::endl;
    TR << "Optimal angle: " << " " << optimal_angle << " zeal score: " << max_corr << std::endl;

if (debug) {
    float optimal_angle_rad_inv = float(optimal_angle)/180.0*3.1415;
     std::ofstream myfile ("opt_rot_grid.dat");
       std::vector< std::vector<int> > rotate_grid_opt;
       core::Size grid_size(slice_for_ZT[0].size());
       for (core::Size i = 0; i< grid_size; i++ ) {
           rotate_grid_opt.push_back( std::vector<int>(grid_size) );
       }

        for (core::Size x=0; x< grid_size; x++ ) {
	    for (core::Size y=0; y< grid_size; y++ ) {
                if (slice_for_ZT[x][y] == 1) {
                    int translated_coords_x = x - grid_size/2;
                    int translated_coords_y = y - grid_size/2;
                    int rotated_x = int(translated_coords_x*cos(optimal_angle_rad_inv) - translated_coords_y*sin(optimal_angle_rad_inv));
                    int rotated_y = int(translated_coords_x*sin(optimal_angle_rad_inv) + translated_coords_y*cos(optimal_angle_rad_inv));
		    rotated_x += grid_size/2;
                    rotated_y += grid_size/2;
                    rotate_grid_opt[rotated_x][rotated_y] = slice_for_ZT[x][y];
                    myfile << rotated_x << " " << rotated_y << std::endl;
               }
           }
       }
       core::scoring::shape::ZernikeDescriptorCalculator zc_calc_opt = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
       zc_calc_opt.invariants_from_grid(rotate_grid_opt);
       
      if ( zernike_transform_type == "2D" ) {
          std::vector<std::complex<double>>   moments_opt = zc_calc_opt.get_2d_moments();
          float zeal_score_opt_again = zeal_score(moments_pose, moments_opt );
	  TR << "Optimal angle: " << " " << optimal_angle << " zeal score again: " << zeal_score_opt_again  << std::endl;
	}

      if ( zernike_transform_type == "2D_3D" ) {
          std::vector<std::complex<double>>   moments_opt = zc_calc_opt.get_3d_moments();
          float zeal_score_opt_again = zeal_score(moments_pose, moments_opt );
	  TR << "Optimal angle: " << " " << optimal_angle << " zeal score again: " << zeal_score_opt_again  << std::endl;
	}
}

    float optimal_angle_rad = float(optimal_angle)/180.0*3.1415;
    numeric::xyzMatrix<core::Real> rot_mat_align = get_rotation_matrix_z_rad( optimal_angle_rad );

    rot_mat_to_facet.transpose();

    //numeric::xyzVector<Real> center_of_pose = core::pose::get_center_of_mass(*pose);
    //numeric::xyzVector<Real> center_of_pose = get_center_of_pose(*pose);
    //TR << "COM: " << center_of_pose(1) << " " << center_of_pose(2) << " " << center_of_pose(3) << std::endl;

core::kinematics::Jump j2 = rotated_pose->jump(2);
j2.reverse();
j2.set_rotation(j2.get_rotation()*rot_mat_align);
rotated_pose->set_jump(2, j2);
//UGLY make sure its get refolded!!!
coord = rotated_pose->residue(1).atom(rotated_pose->residue(1).atom_index("CA")).xyz();

if (debug) {
   rotated_pose->dump_pdb("pose_align.pdb");
   core::scoring::shape::ZernikeDescriptorCalculator zc_calc_final = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
   zc_calc_final.invariants_from_oriented_pose_to_file(*rotated_pose, basic::options::option[out::file::o].value() + "/" + pdb_prefix + ".align.final.inv" );

   std::vector<std::complex<double>>   moments_opt_final;
   if ( zernike_transform_type == "2D" ) {
       moments_opt_final = zc_calc_final.get_2d_moments();
   } else if ( zernike_transform_type == "2D_3D" ) {
       moments_opt_final = zc_calc_final.get_3d_moments();
   }

   float zeal_score_opt_final = zeal_score(moments_grid, moments_opt_final );
   TR << "Optimal angle: " << " " << optimal_angle << " zeal score final: " << zeal_score_opt_final  << std::endl;
}

rotated_pose->replace_residue( rotated_pose->total_residue() -1 , *VRT_coord, false );

j2 = rotated_pose->jump(2);
j2.reverse();
j2.set_rotation(rot_mat_to_facet);
rotated_pose->set_jump(2, j2);
//UGLY make sure its get refolded!!!
coord = rotated_pose->residue(1).atom(rotated_pose->residue(1).atom_index("CA")).xyz();

if (debug) {
   numeric::xyzMatrix<core::Real> R = j2.get_rotation();
    TR << "Rf: " << R(1,1) << " " << R(1,2) << " " << R(1,3) << std::endl
       << "Rf: " << R(2,1) << " " << R(2,2) << " " << R(2,3) << std::endl
       << "Rf: " << R(3,1) << " " << R(3,2) << " " << R(3,3) << std::endl;

   numeric::xyzVector<Real> r1 = R*cart_x;
   numeric::xyzVector<Real> r2 = R*cart_y;
   numeric::xyzVector<Real> r3 = R*cart_z;
   TR << "facet_x final: " << r1(1) << " " << r1(2) << " " << r1(3) << std::endl;
   TR << "facet_y final: " << r2(1) << " " << r2(2) << " " << r2(3) << std::endl;
   TR << "facet_z final: " << r3(1) << " " << r3(2) << " " << r3(3) << std::endl;
}

rotated_pose->dump_pdb("to_facet.pdb");

}
  return 0;
}

