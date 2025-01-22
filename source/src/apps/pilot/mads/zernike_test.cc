//
// Created by Mads Jeppesen on 3/7/22.
//

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
#include <sstream>


// zernike
#include <core/scoring/shape/ZernikeDescriptorEnergy.hh>

#include <string>
#include <vector>
#include <utility/file/FileName.hh>


bool
comparison_test(std::vector<double> a, std::vector<double> b, int precision = 5) {
  auto a_it = a.begin();
  auto b_it = b.begin();
  assert(a.size() == b.size());
  double eps = std::pow(10, -precision);
  while (a_it != a.end()) {
      if ( std::abs(*a_it - *b_it) > eps)
          return false;
      else {
          a_it++;
          b_it++;
      }

  }
  return true;
}

template <typename T>
std::stringstream
vector_to_stream(std::vector<T> vec) {
    std::stringstream ss;
    std::copy(vec.begin(), vec.end(),std::ostream_iterator<T>(ss,","));
    return ss;
}

template <typename T>
std::stringstream
set_to_stream(std::set<T> set) {
    std::stringstream ss;
    std::copy(set.begin(), set.end(),std::ostream_iterator<T>(ss,","));
    return ss;
}

int
main( int argc, char ** argv ) {
    ////// PARAMETER TO INDUCE FAILURE
    bool induce_fail = false;
    ////// PARAMETER TO INDUCE FAILURE

    devel::init(argc, argv);

    int order = 20;
    int size = 64;
    std::string surface_type = "MS";
    core::Real probe = 1.4;
    core::Real shell = 2;

    // define calcutor
    core::scoring::shape::ZernikeDescriptorCalculator zc =
            core::scoring::shape::ZernikeDescriptorCalculator(order, size, surface_type, probe, shell);

    // read in reference data and the pdb list into respective vector = 20;
    std::string line, value;

    std::stringstream floats;
    floats << probe << std::setprecision(2) << "_" << shell << std::setprecision(2);
    std::string identifier = std::to_string(order) + "_" + std::to_string(size) + "_" + surface_type + "_" + floats.str();
    std::string inv_file_str = "reference/expected_results_" + identifier + "_inv.csv" ;
    std::ifstream inv_file(inv_file_str);
    std::string real_file_str = "reference/expected_results_" + identifier + "_imag.csv";
    std::ifstream imag_file(real_file_str);
    std::string imag_file_str = "reference/expected_results_" + identifier + "_real.csv";
    std::ifstream real_file(imag_file_str);
    std::string pdb_file_str = "reference/pdblist.txt";
    std::ifstream pdblist(pdb_file_str);
    // assert these files exist!
    try {
        if (inv_file.fail())
            utility_exit_with_message(inv_file_str + " does not exist!");
        if (real_file.fail())
            utility_exit_with_message(real_file_str + " does not exist!");
        if (imag_file.fail())
            utility_exit_with_message(imag_file_str + " does not exist!");
        if (imag_file.fail())
            utility_exit_with_message(pdb_file_str + " does not exist!");
    }
    catch ( utility::excn::Exception & e) {
        std::cout << "Make sure you are in the shape/tests repository!" << std::endl;
        return 1;
    }
    std::vector<std::vector<double>> imag, real, inv;
    std::vector<std::string> pdbs;
    if (inv_file.is_open()) {
        while ( getline (inv_file, line) ) {
            std::stringstream str(line);
            inv.emplace_back(std::vector<double>());
            while (std::getline(str, value, *","))
                inv.back().push_back(std::stod(value));
        }
        inv_file.close();
    }
    if (real_file.is_open()) {
        while ( getline (real_file, line) ) {
            std::stringstream str(line);
            real.emplace_back(std::vector<double>());
            while (std::getline(str, value, *","))
                real.back().push_back(std::stod(value));
        }
        real_file.close();
    }
    if (imag_file.is_open()) {
        while (getline(imag_file, line)) {
            std::stringstream str(line);
            imag.emplace_back(std::vector<double>());
            while (std::getline(str, value, *","))
                imag.back().push_back(std::stod(value));
        }
        imag_file.close();
    }
    if (pdblist.is_open()) {
        getline(pdblist, line); // skip the header
        while ( getline(pdblist, line) )
            pdbs.push_back(line);
        pdblist.close();
    }

    // calculate moments and invs from the reference pdbs
    std::vector<std::vector<double>> imag_new, real_new, inv_new;
    for (std::string pdb : pdbs) {
        core::pose::PoseOP pose = core::import_pose::pose_from_file(pdb);
        inv_new.push_back(zc.invariants_from_pose(*pose));
        std::vector<std::complex<double>> moments = zc.moments_from_pose(*pose);
        std::vector<double> real_moments;
        std::vector<double> imag_moments;
        for (std::complex<double> moment : moments) {
            real_moments.push_back(moment.real());
            imag_moments.push_back(moment.imag());
        }
        ////// CODE IMPLEMENTED TO INDUCE FAILURE
        if (induce_fail) {
            std::cout << "inducing fail" << std::endl;
            real_moments[2] = 300;
            imag_moments[2] = 5500;
            inv[0][2] = 3005;
        }
        ///// CODE IMPLEMENTED TO INDUCE FAILURE
        real_new.push_back(real_moments);
        imag_new.push_back(imag_moments);
    }

    // check if the references are the same
    auto inv_it = inv.begin();
    auto real_it = real.begin();
    auto imag_it = imag.begin();
    auto inv_new_it = inv_new.begin();
    auto real_new_it = real_new.begin();
    auto imag_new_it = imag_new.begin();
    auto pdb_it = pdbs.begin();
    int failures = 0;
    std::set<std::string> failed_pdbs;
    while (pdb_it != pdbs.end()) {
        std::cout << "==================" << std::endl;
        std::cout << *pdb_it << std::endl;
        if (comparison_test(*inv_it, *inv_new_it))
            std::cout << "inv: True" << std::endl;
        else {
            std::cout << "inv: False" << std::endl;
            std::cout << "Expected: " << vector_to_stream(*inv_it).str() << std::endl;
            std::cout << "Got:      " << vector_to_stream(*inv_new_it).str() << std::endl;
            failures += 1;
            failed_pdbs.emplace(*pdb_it);
        }
        if (comparison_test(*real_it, *real_new_it))
            std::cout << "real: True" << std::endl;
        else {;
            std::cout << "real: False" << std::endl;
            std::cout << "Expected: " << vector_to_stream(*real_it).str() << std::endl;
            std::cout << "Got:      " << vector_to_stream(*real_new_it).str() << std::endl;
            failures += 1;
            failed_pdbs.emplace(*pdb_it);
        }
        if (comparison_test(*imag_it, *imag_new_it))
            std::cout << "imag: True" << std::endl;
        else {
            std::cout << "imag: False" << std::endl;
            std::cout << "Expected: " << vector_to_stream(*imag_it).str() << std::endl;
            std::cout << "Got:      " << vector_to_stream(*imag_new_it).str() << std::endl;
            failures += 1;
            failed_pdbs.emplace(*pdb_it);
        }
        pdb_it++; // increment the pdbs
        inv_it++; real_it++; imag_it++;  // increment the reference values
        inv_new_it++; real_new_it++; imag_new_it++; // increment the calculated values
//        break; /// DELETE!
    }
    if (failures == 0) {
        std::cout << "ALL TESTS PASSED" << std::endl;
    }
    else {
        std::cout << std::to_string(failures) << " TESTS FAILED!" << std::endl;
        std::cout << "THE FOLLOWING PDBS FAILED:" << std::endl;
        std::cout << set_to_stream(failed_pdbs).str() << std::endl;
    }

    // now compare them

    return 0;
}