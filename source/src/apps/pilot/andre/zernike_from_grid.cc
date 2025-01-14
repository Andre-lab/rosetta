
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
#include <core/scoring/shape/VoxelGrid.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <iostream>


// protocols
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/scoring/shape/ZernikeDescriptorEnergy.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/zernike_descriptor.OptionKeys.gen.hh>


static basic::Tracer TR( "apps.pilot.zernike" );

int
main( int argc, char ** argv ) {

    using namespace core;
    using namespace utility;
    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    devel::init( argc, argv );

    // reading in structure
    vector1<std::string> files;
    vector1<file::FileName> list;
    if ( basic::options::option[in::file::s].user() ) {
        TR << "using -s option" << std::endl;
        files=option[in::file::s]();
    } else if ( option[in::file::l].user() ) {
        TR << "using -l option " << std::endl;
        list = basic::options::option[ in::file::l ]();
        for ( file::FileName const & listfile : list ) {
            utility::io::izstream grids(listfile);
            std::string fname;
            while ( grids >> fname ) {
                files.push_back(fname);
            }
        }
    }


    // options for zernike calculations
    int order = option[zernike_descriptor::order].value();
    int size = option[zernike_descriptor::grid_size].value();
    std::string surface_type = option[zernike_descriptor::surface_type].value();
    core::Real probe = option[zernike_descriptor::probe_radius].value();
    core::Real shell = option[zernike_descriptor::shell_thickness].value();
    std::string zernike_transform_type = option[zernike_descriptor::zernike_transform_type].value();
    bool outline_transform = option[zernike_descriptor::transform_outline].value();

    // processing each file in turn outputting the zernike descriptors
    for ( auto const & file : files ) {
        TR << "Processing file: " << file << std::endl;
        core::scoring::shape::ZernikeDescriptorCalculator zc_calc = core::scoring::shape::ZernikeDescriptorCalculator(order,size, surface_type, probe,shell);
        TR << "printing inv" << std::endl;
        if ( outline_transform) {
        	zc_calc.invariants_from_grid_outline_file_to_file(file, basic::options::option[out::file::o].value() + "/" + file + ".inv");
        } else {
		zc_calc.invariants_from_grid_file_to_file(file, basic::options::option[out::file::o].value() + "/" + file + ".inv");
       	}
    }

    return 0;
}
