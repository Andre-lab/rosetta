#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/task_operations/SelectByDeltaSASAOperation.hh>
#include <string>
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>

// zernike
#include <core/scoring/shape/ZernikeDescriptorEnergy.hh>

int
main( int argc, char ** argv ) {
    // init Rosetta
    devel::init(argc, argv);

    // initialize the SelectByDeltaSASAOperation class and set options
    protocols::task_operations::SelectByDeltaSASAOperation sasaop = protocols::task_operations::SelectByDeltaSASAOperation();
    sasaop.set_layers(true, true, true);
//    sasaop.set_layer_asa(-30, -5);
    sasaop.set_verbose(true);
    sasaop.set_output(true);

    // create symmetrical pose
    sasaop.set_output(true);
    std::string pdb_file("input/1stm.pdb");
    std::string sym_file("input/1stm.symm");
    core::pose::PoseOP pose_sym = core::import_pose::pose_from_file(pdb_file);
    protocols::symmetry::SetupForSymmetryMover setup_symmetry = protocols::symmetry::SetupForSymmetryMover(
            sym_file);
    std::string name_before = pose_sym->pdb_info()->name();
    setup_symmetry.apply(*pose_sym);
    pose_sym->pdb_info()->name(name_before) ;

    // TEST SYMMETRICAL CASE
    sasaop.sym_dof_names("JUMP5fold1,JUMP5fold111");
    core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task(*pose_sym);
    sasaop.apply(*pose_sym, *dummy_task);

    // create non-symmetrical_pose
    core::pose::PoseOP pose_non_sym = core::import_pose::pose_from_file("input/grafted_pose.pdb");

    // TEST NON SYMMETRICAL CASE
    sasaop.sym_dof_names("");
    sasaop.set_jump_nums("2");
    dummy_task = core::pack::task::TaskFactory::create_packer_task(*pose_non_sym);
    sasaop.apply(*pose_non_sym, *dummy_task);

    return 0;
    }


