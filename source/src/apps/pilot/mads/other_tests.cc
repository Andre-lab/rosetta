#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/task_operations/SelectBySASAOperation.hh>
#include <string>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/TaskFactory.hh>

// zernike
#include <core/scoring/shape/ZernikeDescriptorEnergy.hh>

int
main( int argc, char ** argv ) {

    std::string pdb_file("/home/shared/databases/capsids/T1/prepared/unrelaxed/native_pdb/1stm.pdb");
    devel::init(argc, argv);


    // read stuff in
    utility::vector1<std::string> filenames = basic::options::option[basic::options::OptionKeys::in::file::s].value();

    // create symmetrical pose
    std::string sym_file("/home/shared/databases/capsids/T1/prepared/symmdef/1stm.symm");
    core::pose::PoseOP pose = core::import_pose::pose_from_file(pdb_file);
    protocols::symmetry::SetupForSymmetryMover setup_symmetry = protocols::symmetry::SetupForSymmetryMover(
            sym_file);
    setup_symmetry.apply(*pose);

    // put it into pymol
    protocols::moves::PyMOLMover pmm = protocols::moves::PyMOLMover(); // std::string const & address, unsigned int port, unsigned int max_packet_siz
    pmm.keep_history(true);
    pmm.apply(*pose);

    // create sasa selection
    protocols::task_operations::SelectBySASAOperation sasaop = protocols::task_operations::SelectBySASAOperation(
            "sc", "delta");
    core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task(*pose);
    sasaop.set_layers(true, true, true);
    sasaop.set_layer_asa(-30, -5);
    sasaop.set_verbose(true);
    sasaop.sym_dof_names("JUMP5fold1,JUMP5fold111");
    sasaop.apply(*pose, *dummy_task);
    //

    return 0;
}