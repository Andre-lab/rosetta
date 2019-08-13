// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#include <numeric/constants.hh>

#include <core/pose/Pose.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>

#include <utility/graph/Graph.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <devel/init.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/util.hh>//option.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.fwd.hh>
#include <protocols/comparative_modeling/LoopRelaxMover.hh>
#include <core/kinematics/FoldTree.hh>
#include <numeric/random/random.hh>
#include <protocols/loops/Loops.hh>
#include <core/chemical/VariantType.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pose/disulfide_util.hh>
#include <core/conformation/util.hh>
#include <protocols/loops/LoopRelaxMover.fwd.hh>
#include <devel/init.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/disulfides/FullatomDisulfidePotential.hh>
#include <protocols/jobdist/standard_mains.hh>


//#include <protocols/jobdist/standard_mains.hh>
//#include <protocols/moves/Mover.hh>
#include <core/io/raw_data/ScoreFileData.hh>

#include <ObjexxFCL/format.hh>
#define INF 99999999
// C++ headers
//#include <cstdlib>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/JobDistributors.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>

#define CYSPOS 2
using basic::Warning;
using basic::Error;
using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
namespace ss_loop_builder {
basic::options::IntegerOptionKey start_cys("ss_loop_builder:start");
basic::options::IntegerOptionKey   end_cys("ss_loop_builder:end");
basic::options::IntegerOptionKey   ncycles("ss_loop_builder:ncycles");
}
std::string int2string(int n){
	std::stringstream Num;
	std::string str;
	Num << n;
	str = Num.str();
	return str;
}

void modelSSLoop(Size startCys, Size endCys, pose::Pose& workpose, core::scoring::ScoreFunctionOP scorefxn){
	protocols::comparative_modeling::LoopRelaxMover loop_relax_mover;
	utility::vector1< std::pair<Size, Size> > dss;
	std::pair<Size,Size> bond(startCys,endCys);
	dss.push_back(bond);
	workpose.conformation().fix_disulfides(dss);
	//loop_relax_mover.scorefxns(scorefxn,scorefxn);
	loop_relax_mover.remodel("perturb_kic");
	loop_relax_mover.refine("refine_kic");
	loop_relax_mover.fa_scorefxn(scorefxn);
	loop_relax_mover.copy_sidechains(true);
	protocols::loops::Loops loop;
	std::cout<<"MODELING: \n";
	for ( Size i=startCys+1; i<endCys; ++i ) {
		std::cout<<workpose.residue(i).name();
	}
	std::cout<<std::endl;
	std::cout<<"TOTAL PEPTIDE: "<<workpose.sequence()<<std::endl;

	loop.add_loop(startCys+1,endCys-1,(startCys+endCys)/2);
	loop_relax_mover.loops(loop);
	loop_relax_mover.apply(workpose);
}
void printEnergies(pose::Pose& workpose) {
	for ( Size i=1; i<=workpose.size(); ++i ) {
		std::cout<<"=========Residue "<< i<<": "<< workpose.residue(i).name()<<"============"<<std::endl;
		workpose.energies().residue_total_energies(i).print();
	}
}

void rotateUntilCys(pose::Pose& workpose, Size untilCys){

	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,workpose.size());
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,1);

	// while ( untilCys  %  workpose.size() != 1){

	Size counter = 1;
	pose::Pose tempPose(workpose,untilCys,untilCys);
	core::pose::remove_lower_terminus_type_from_pose_residue(tempPose,tempPose.size());
	core::pose::remove_upper_terminus_type_from_pose_residue(tempPose,1);
	for ( Size j=untilCys+1; j%workpose.size()!=untilCys; ++j ) {
		//handle the case where j&n_residue == 0.. we actually want it to point to the last residue
		Size resi = (j%workpose.size() == 0)?workpose.size():j%workpose.size();
		std::cout<<resi<<std::endl;
		tempPose.append_polymer_residue_after_seqpos(workpose.residue(resi),counter++,false);
		std::cout<<"Counter: "<<counter<<std::endl;
	}
	//  std::cout<<"BEFORE====================="<<std::endl;
	//
	//  core::pose::remove_lower_terminus_type_from_pose_residue(tempPose,1);
	//  core::pose::remove_upper_terminus_type_from_pose_residue(tempPose,1);
	//  std::cout<<"AFTER====================="<<std::endl;
	//
	//  for( Size j=1; j<workpose.size(); ++j){
	//   tempPose.conformation().safely_append_polymer_residue_after_seqpos(workpose.residue(j),j,false);
	//   //tempPose.append_polymer_residue_after_seqpos(workpose.residue(j),j,false);
	//  }
	//
	//  std::cout<<tempPose.sequence()<<std::endl;
	//  core::pose::remove_lower_terminus_type_from_pose_residue(tempPose,workpose.size());
	//  core::pose::remove_upper_terminus_type_from_pose_residue(tempPose,1);
	workpose = tempPose;
	core::pose::remove_lower_terminus_type_from_pose_residue(workpose,workpose.size());
	core::pose::remove_upper_terminus_type_from_pose_residue(workpose,1);
	//  untilCys++;
	workpose.dump_pdb("tempPdb.pdb");
	workpose.conformation().declare_chemical_bond(1,"N",workpose.size(),"C");
	//  workpose.set_omega(w,-180.0);
	// }

}
void updateSSAtoms(pose::Pose& workpose, Size upNum, utility::vector1_int& vec, Size pepsize) {
	Size toAdd = pepsize+1-upNum;
	utility::vector1<std::pair<Size,Size> > dss;
	for ( Size i=1; i<=vec.size(); ++i ) {
		vec[i] = (vec[i] + toAdd) % pepsize;
	}
	std::sort(vec.begin(),vec.end());
	for ( Size i=1; i<=vec.size(); i+=2 ) {
		std::pair<Size,Size> p(vec[i],vec[i+1]);
		dss.push_back(p);
	}

	workpose.conformation().fix_disulfides(dss);
	// for (Size i=1; i<=vec.size(); i+=2 ){
	////  core::conformation::change_cys_state(vec[i],"CYD",workpose.conformation());
	////  core::conformation::change_cys_state(vec[i+1],"CYD",workpose.conformation());
	//  std::cout<<"BEFORE FORMIMG: "<<workpose.residue(vec[i]).name()<<" "<<workpose.residue(vec[i+1]).name()<<std::endl;
	//  if (!core::conformation::is_disulfide_bond(workpose.conformation(),vec[i],vec[i+1])){
	//   //core::conformation::form_disulfide(workpose.conformation(),vec[i],vec[i+1]);
	//  }
	//  std::cout<<"AFTER FORMIMG: "<<workpose.residue(vec[i]).name()<<std::endl;
	// }
}

core::scoring::constraints::ConstraintSetOP createConstraintSet(pose::Pose& workpose, utility::vector1 < std::pair< Size, Size > >& dss){
	core::scoring::constraints::ConstraintSetOP cst_set( new core::scoring::constraints::ConstraintSet() );
	core::scoring::disulfides::SG_Dist_FuncOP spring = new core::scoring::disulfides::SG_Dist_Func();
	conformation::Conformation const & conformation( workpose.conformation() );
	for ( Size i=1; i<=dss.size(); ++i ) {
		core::pose::Pose::AtomID atomsS1 (workpose.residue(dss[i].first).atom_index("SG"),dss[i].first);
		core::pose::Pose::AtomID atomsS2 (workpose.residue(dss[i].second).atom_index("SG"),dss[i].second);
		cst_set->add_constraint(new core::scoring::constraints::AtomPairConstraint(atomsS1,atomsS2,spring));
	}
	return cst_set;
}

void printScoreTermsPose(pose::Pose& workpose,std::string filename){
	std::cout<<"============SCORE:================="<<std::endl;
	std::string scorefile_name = "score";
	std::map < std::string, core::Real > score_map;
	score_map = protocols::jobdist::get_score_map(workpose);
	core::io::raw_data::ScoreFileData sfd(scorefile_name);
	sfd.write_pose(workpose, score_map,filename);
}
core::scoring::constraints::ConstraintSetOP createMinimizationConstraint(pose::Pose& workpose){
	core::scoring::constraints::ConstraintSetOP cst_set( new core::scoring::constraints::ConstraintSet() );
	core::pose::Pose::AtomID atom1(workpose.residue(1).atom_index("N"),1);
	core::pose::Pose::AtomID atom2(workpose.residue(workpose.size()).atom_index("C"),workpose.size());
	Real dist = workpose.residue(1).atom("N").xyz().distance(workpose.residue(workpose.size()).atom("C").xyz());
	core::scoring::constraints::HarmonicFuncOP spring = new core::scoring::constraints::HarmonicFunc(dist,0.00001);

	cst_set->add_constraint(new core::scoring::constraints::AtomPairConstraint(atom1,atom2,spring));
	return cst_set;
}

void deleteEdge( core::scoring::EnergyGraph&  g){
	utility::graph::EdgeListIterator it = g.edge_list_begin();
	while ( it != g.edge_list_end() ) {
		if ( (*it)->get_first_node_ind() == 1 && (*it)->get_second_node_ind() == 14 ) {
			g.delete_edge(*it);
			return;
		}
		++it;
	}
}

core::scoring::constraints::ConstraintSetOP createDihedralConstraint(pose::Pose& workpose) {
	core::scoring::constraints::ConstraintSetOP cst_set( new core::scoring::constraints::ConstraintSet() );
	core::scoring::constraints::HarmonicFuncOP spring = new core::scoring::constraints::HarmonicFunc(180,1);
	id::AtomID atom1(workpose.residue(workpose.size()).atom_index("CA"),workpose.size());
	id::AtomID atom2(workpose.residue(workpose.size()).atom_index("C"),workpose.size());
	id::AtomID atom3(workpose.residue(1).atom_index("N"),1);
	id::AtomID atom4(workpose.residue(1).atom_index("CA"),1);

	cst_set->add_constraint(new core::scoring::constraints::DihedralConstraint(atom1,atom2,atom3,atom4,spring));
	return cst_set;
}

void performMinimization(pose::Pose& workpose, core::scoring::ScoreFunctionOP scorefxn){
	core::kinematics::MoveMapOP tempMap = new core::kinematics::MoveMap();

	tempMap->set_bb_true_range(1,workpose.size());
	protocols::minimization_packing::MinMover minimizer(tempMap, scorefxn, "lbfgs_armijo_atol", 0.0001, true /*nb_list*/ );
	std::cout<<"Score before minimization: "<<scorefxn->score(workpose)<<std::endl;
	workpose.dump_pdb("before_minimization.pdb");
	minimizer.apply(workpose);
	workpose.dump_pdb("after_minimization.pdb");

	std::cout<<"Score after minimization: "<<scorefxn->score(workpose)<<std::endl;
}

void packRotamers( pose::Pose& workpose, core::scoring::ScoreFunctionOP scorefxn){
	pack::task::PackerTaskOP task;
	task = pack::task::TaskFactory::create_packer_task( workpose );
	task->initialize_from_command_line().restrict_to_repacking();//flags+no-design
	pack::task::operation::NoRepackDisulfides noRepackDisulf;
	pack::task::operation::PreventRepacking prevent;
	prevent.include_residue(1);
	prevent.include_residue(workpose.size());
	prevent.apply(workpose,*task);
	noRepackDisulf.apply(workpose, *task);
	protocols::minimization_packing::PackRotamersMoverOP packer = new protocols::minimization_packing::PackRotamersMover(scorefxn, task);
	packer->apply(workpose);
	std::cout<<"Score after side chain repacking: "<<scorefxn->score(workpose)<<std::endl;
}

void writeScores(pose::Pose& workpose, pose::Pose& nativePose, const std::string filename, Real score) {
	Real rmsd = protocols::simple_filters::native_CA_rmsd(nativePose,workpose);
	std::ofstream writer("scorefile.out",std::ios::app);
	writer<<filename<<"\t"<<score<<"\t"<<rmsd<<std::endl;
	writer.close();
}

int main (int argc, char** argv) {
	try {
		devel::init(argc,argv);
		basic::Tracer TR( "protocols.moves.CycPep" );
		Size nstruct = option[out::nstruct];
		pose::Pose nativePose;
		core::import_pose::pose_from_file(nativePose, basic::options::start_file(), core::import_pose::PDB_file);
		for ( Size s=0; s<nstruct; ++s ) {
			pose::Pose workpose = nativePose;
			core::pose::remove_lower_terminus_type_from_pose_residue(workpose,1);
			core::pose::remove_upper_terminus_type_from_pose_residue(workpose,workpose.size());

			option[ out::file::fullatom ].value(true);
			core::scoring::ScoreFunctionOP scorefxn(
				core::scoring::get_score_function());
			utility::vector1 < std::pair< Size, Size > > dss;

			utility::vector1_int ssAtoms;
			workpose.dump_pdb("before_minimization_first.pdb");
			core::conformation::disulfide_bonds(workpose,dss);
			workpose.conformation().declare_chemical_bond(1,"N",workpose.size(),"C");
			//
			performMinimization(workpose,scorefxn);
			packRotamers(workpose,scorefxn);
			workpose.dump_pdb("minimized_repacked.pdb");
			TR<<"Number of d-Sulfides: "<<dss.size()<<std::endl;
			for ( Size i=1; i<=dss.size(); ++i ) {
				TR<<"Cys bond at "<<dss[i].first<<" "<<dss[i].second<<std::endl;
				ssAtoms.push_back(dss[i].first);
				ssAtoms.push_back(dss[i].second);
			}
			std::sort(ssAtoms.begin(),ssAtoms.end());

			//  for (Size i=1; i<=ssAtoms.size(); ++i){
			//   workpose.dump_pdb("before_rotating"+int2string(i)+".pdb");
			//   rotateUntilCys(workpose,ssAtoms[i]);
			//   workpose.dump_pdb("after_rotating"+int2string(i)+".pdb");
			//   updateSSAtoms(workpose,ssAtoms[i],ssAtoms,workpose.size());
			//   performMinimization(workpose,scorefxn);
			//   workpose.dump_pdb("workpose_"+int2string(i)+".pdb");
			//  }

			for ( Size i=1; i<=ssAtoms.size(); ++i ) {
				rotateUntilCys(workpose,ssAtoms[i]);
				updateSSAtoms(workpose,ssAtoms[i],ssAtoms,workpose.size());
				workpose.dump_pdb("after_rotating"+int2string(i)+".pdb");
				core::scoring::constraints::ConstraintSetOP cst_set = createDihedralConstraint(workpose);
				workpose.add_constraints(cst_set->get_all_constraints());
				scorefxn->set_weight(core::scoring::dihedral_constraint,1);
				performMinimization(workpose,scorefxn);

				TR<<"Score after minimization (2): " <<scorefxn->score(workpose)<<std::endl;
				dss.clear();
				core::conformation::disulfide_bonds(workpose,dss);
				modelSSLoop(ssAtoms[1],ssAtoms[2],workpose,scorefxn);
				performMinimization(workpose,scorefxn);
				packRotamers(workpose,scorefxn);
				scorefxn->set_weight(core::scoring::dihedral_constraint,0);
				TR<<"Score after loop modeling+minimization: "<<scorefxn->score(workpose)<<std::endl;
				//   printEnergies(workpose);
				//   std::cout<<"---------------------Trying to connect "<<workpose.residue(ssAtoms[1]).name()<<", "<<workpose.residue(ssAtoms[2]).name();
				//   std::cout<<"At places "<<ssAtoms[1]<<" "<<ssAtoms[2]<<std::endl;
				//   for(Size j=1; j<=dss.size(); ++j){
				//    core::conformation::form_disulfide(workpose.conformation(),dss[j].first,dss[j].second);
				//    core::pose::rebuild_disulfide(workpose,dss);
				//
				//   }
				workpose.dump_pdb("workpose_loopm_"+int2string(i)+".pdb");
				workpose.remove_constraints(cst_set->get_all_constraints());
			}
			writeScores(workpose,nativePose,"Job_"+int2string(s)+".pdb",scorefxn->score(workpose));
			workpose.dump_pdb("Job_"+int2string(s)+".pdb");
		}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}


//int main(int argc, char** argv) {
// using namespace core;
// using namespace basic::options;
// using namespace basic::options::OptionKeys;
//
// devel::init(argc,argv);
// core: :scoring::ScoreFunctionOP scorefxn_ = core::scoring::get_score_function();
// protocols::loops::Loops myloop;
// myloop.add_loop(2,6,4);
// protocols::loops::LoopRelaxMover loop_relax_mover;
// loop_relax_mover.remodel("perturb_kic");
// loop_relax_mover.refine("refine_kic");
// loop_relax_mover.loops(myloop);
// loop_relax_mover.copy_sidechains(false);
// basic::Tracer TR("protocols.moves.CycPep");
// pose::Pose mypose;
// core::import_pose::pose_from_file(mypose, options::start_file(), core::import_pose::PDB_file);
// Size pepsize = mypose.size();
// pose::Pose temp = mypose;
//// temp.dump_pdb("temp_0.pdb");
// core::pose::remove_lower_terminus_type_from_pose_residue(temp,1);
// core::pose::remove_upper_terminus_type_from_pose_residue(temp,pepsize);
//// temp.dump_pdb("temp_1.pdb");
////
//// for (Size i=1; i<=pepsize-1; ++i){
////  temp.append_polymer_residue_after_seqpos(mypose.residue(i),i,false);
////  temp.dump_pdb("temp_"+int2string(i)+".pdb");
//// }
////
// TR << "INITIAL Score: " << scorefxn_->score(temp)<<std::endl;
// return 0;
// temp.dump_pdb("start.pdb");
// for (Size i=1; i<=pepsize; ++i){
////  core::pose::remove_lower_terminus_type_from_pose_residue(temp,pepsize);
////  core::pose::remove_upper_terminus_type_from_pose_residue(temp,1);
//  pack::task::PackerTaskOP task;
//  task = pack::task::TaskFactory::create_packer_task( temp );
//  pack::task::operation::NoRepackDisulfides noRepackDisulf;
//  noRepackDisulf.apply(temp, *task);
//  loop_relax_mover.apply(temp);
//  //temp.dump_pdb("test"+int2string(i)+".pdb");
//  //mypose.clear();
//  //mypose.append_residue_by_bond(temp.residue(pepsize));
////  temp.dump_pdb("test_temp_before.pdb");
//
////  if (i==1)
////   temp.dump_pdb("test_temp_after.pdb");
//  pose::Pose tmp2(temp, pepsize, pepsize);
//  core::pose::remove_lower_terminus_type_from_pose_residue(tmp2,1);
//  core::pose::remove_upper_terminus_type_from_pose_residue(tmp2,1);
//  for( Size j=1; j<=pepsize-1; ++j){
//   tmp2.append_polymer_residue_after_seqpos(temp.residue(j),j,false);
//  }
//  //TR << tmp2.sequence() << std::endl;
//
//  core::pose::remove_lower_terminus_type_from_pose_residue(tmp2,pepsize);
//  core::pose::remove_upper_terminus_type_from_pose_residue(tmp2,1);
//  tmp2.dump_pdb("tst2"+int2string(i)+".pdb");
////  for (Size i=1; i<=pepsize; ++i){
////   TR << "res: " << i << " lower: " << tmp2.residue(i).has_variant_type(core::chemical::LOWER_TERMINUS) <<std::endl;
////   TR << "res: " << i << " upper: " << tmp2.residue(i).has_variant_type(core::chemical::UPPER_TERMINUS) <<std::endl;
////  }
//  temp = tmp2;
//  break;
// }
//
//// for (Size i=1; i<=pepsize; ++i){
////  //loop_relax_mover.apply(mypose);
////
////  mypose.replace_residue(1,temp.residue(pepsize),true);
////  std::cerr<<"TTTTTTT\n";
////  for( Size j=1; j<=pepsize-1; ++j){
////   mypose.replace_residue(j+1,temp.residue(j),true);
////   mypose.conformation().update_polymeric_connection(j);
////   mypose.conformation().update_polymeric_connection(j+1);
////  }
////  temp = mypose;
//// }
//
// /**
//  * This is the actual definition of a new edge + cut required!!!
//  */
//// fd.delete_edge(e);
//// fd.add_edge(14,1,-1);
//// fd.new_jump(5,7,6);
//// mypose.fold_tree(fd);
//// for (Size i =1; i<2; ++i){
////  std::cout<<i<<" "<<mypose.size()-i<<" "<<(mypose.size()-i)/2<<std::endl;
////  myloop.add_loop(3,13,8);
////
////  loop_relax_mover.loops(myloop);
////  loop_relax_mover.apply(mypose);
////
//// }
// io::pdb::dump_pdb(mypose,"test.pdb");
//
//}
