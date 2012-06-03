// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c)University of Washington UW TechTransfer, email:license@u.washington.edu.

/// @file protocols/antibody2/ModelCDRH3.cc
/// @brief models CDR H3 loop using loop modeling
/// @detailed
///// @author Jianqing Xu ( xubest@gmail.com )
//


#include <protocols/antibody2/ModelCDRH3.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/loops/LoopMoverFactory.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <utility/exit.hh>
#include <protocols/antibody2/AntibodyUtil.hh>
#include <protocols/antibody2/H3PerturbCCD.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/H3CterInsert.hh>
#include <protocols/antibody2/RefineCDRH1Centroid.hh>
#include <protocols/moves/PyMolMover.hh>


static basic::Tracer TR("protocols.antibody2.ModelCDRH3");

using namespace core;

namespace protocols {
namespace antibody2 {

ModelCDRH3::ModelCDRH3() : Mover(){}

    
ModelCDRH3::~ModelCDRH3() {}
    
ModelCDRH3::ModelCDRH3( AntibodyInfoOP antibody_info) : Mover(){
    user_defined_ = false;
    ab_info_ = antibody_info;
    
    init();
}
    

    
ModelCDRH3::ModelCDRH3( AntibodyInfoOP antibody_info,                 
                        core::scoring::ScoreFunctionCOP lowres_scorefxn) : Mover()
{
	user_defined_ = true;
    ab_info_ = antibody_info;
    lowres_scorefxn_  = new core::scoring::ScoreFunction(*lowres_scorefxn);
    
	init();
}



    
void ModelCDRH3::init( )
{
	Mover::type( "ModelCDRH3" );

	set_default();

    //TODO:
    //JQX: need to deal with this
    if( is_camelid_ && !ab_info_->is_extended() && !ab_info_->is_kinked() ){
        c_ter_stem_ = 0;
    }
    
    h3_cter_insert_mover_ = new H3CterInsert(ab_info_, is_camelid_);
    h3_perturb_ccd_build_ = new H3PerturbCCD(ab_info_, lowres_scorefxn_);   
}

    
void ModelCDRH3::set_default()
{
    benchmark_          = false;
    is_camelid_         = false;
    do_cter_insert_     = true;
    loops_flag_         = true;
    dle_flag_           = true;
    use_pymol_diy_      = false;

    remodel_            = "legacy_perturb_ccd";
        
    c_ter_stem_ = 3;
    max_cycle_ = 20;
        
    cen_cst_ = 10.0;

        
    if(!user_defined_)
    {
        lowres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
            lowres_scorefxn_->set_weight( scoring::chainbreak, 10./3. );
            lowres_scorefxn_->set_weight( scoring::atom_pair_constraint, cen_cst_ );
    }
} 

    

void ModelCDRH3::set_lowres_score_func(scoring::ScoreFunctionCOP lowres_scorefxn ){
    lowres_scorefxn_ = new core::scoring::ScoreFunction(*lowres_scorefxn);
}
    
    
void ModelCDRH3::set_task_factory(pack::task::TaskFactoryCOP tf){
    tf_ = new pack::task::TaskFactory(*tf);
}    

    

void ModelCDRH3::turn_off_H3_filter(){
    h3_perturb_ccd_build_->turn_off_H3_filter();
}    

void ModelCDRH3::turn_on_and_pass_the_pymol(moves::PyMolMoverOP pymol){        
    use_pymol_diy_ = true;
    pymol_ = pymol;
    h3_cter_insert_mover_->turn_on_and_pass_the_pymol(pymol);
    h3_perturb_ccd_build_->turn_on_and_pass_the_pymol(pymol);
}
    

void ModelCDRH3::apply( pose::Pose & pose_in )
{


    TR << "Applying CDR H3 modeler" << std::endl;

    using namespace core::pose;
    using namespace core::scoring;
    using namespace protocols::moves;


    pose::Pose start_pose = pose_in;



    Size framework_loop_begin( ab_info_->get_CDR_loop("h3")->start() );
    Size framework_loop_end  ( ab_info_->get_CDR_loop("h3")->stop()  );
    Size cutpoint = ab_info_->get_CDR_loop("h3")->cut() ; // keep the cutpoint unchanged
    Size framework_loop_size = (framework_loop_end - framework_loop_begin) + 1;

    loops::Loop cdr_h3( framework_loop_begin, framework_loop_end, cutpoint, 0, true );
    loops::Loop trimmed_cdr_h3(framework_loop_begin, framework_loop_end - c_ter_stem_, cutpoint, 0, true );
    loops::Loop input_loop;
    
    if( framework_loop_size <= 6 ) {
        do_cter_insert_ = false;
        TR<<"loop_size <= 6, AUTOMATICALLY TURNING OFF THE C_TERMINAL INSERT"<<std::endl;
    }
    
    if (do_cter_insert_){
        //JQX: the h3 loop removing the cterminal 3 residues
        input_loop = trimmed_cdr_h3;
    }
    else{
        //JQX: the original h3 loop
        input_loop = cdr_h3;
    }
        

    simple_one_loop_fold_tree( pose_in, cdr_h3 );

    // switching to centroid mode
    simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
    simple_moves::SwitchResidueTypeSetMover to_full_atom( chemical::FA_STANDARD );

    
    
    // Building centroid mode loop
    to_centroid.apply( pose_in );

    if(remodel_=="legacy_perturb_ccd"){
    // some initialization before you do h3 loop modeling
//    my_LoopMover xxx ;
//    xxx.set_extended_torsions( pose_in, cdr_h3 );
       set_extended_torsions( pose_in, cdr_h3 );
       //pose_in.dump_pdb("extended_idealized_centroid.pdb");
        //JQX:  this function is in loops_main.cc file
        //      firstly, idealize the loop (indealize bonds as well)
        //      phi(-150),  all the residue, except the first one
        //      psi(150),   all the residue, except the last one
        //      omega(180), all the residue, except the first & last one
        //JQX:  in R2: the function is called "insert_init_frag", which is 
        //      in the file "jumping_util.cc". All the phi, psi, omega are 
        //      assigned to all the residues. "L" secondary structure is 
        //      also assinged. The bonds are idealized using 
        //      framework_pose.insert_ideal_bonds(begin-1, end)
        
        
        h3_perturb_ccd_build_->pass_the_loop(input_loop);

    }
    else{
        // use H3 to define a loops object
        loops::LoopsOP pass_loops = new loops::Loops(); 
        pass_loops->add_loop( input_loop   );
        pass_loops->set_extended(true); // the IndepdentLoopMover will extend the loop
        
        // create a LoopMover type based on the string remode_
        remodel_mover_ = static_cast< loops::loop_mover::IndependentLoopMover * >
        ( loops::LoopMoverFactory::get_instance()->create_loop_mover(remodel_, pass_loops).get() ) ;
        if ( !remodel_mover_ ) { utility_exit_with_message( "Error: no remodel mover defined!" );}
        
        // deal with the fragment files if the remodel_ type is not KIC
        if ( remodel_ != "perturb_kic"  ){
            utility::vector1< core::fragment::FragSetOP > frag_libs;
            loops::read_loop_fragments( frag_libs );
            runtime_assert( frag_libs.size() > 0 );
            for ( Size i = 1; i <= frag_libs.size(); ++i ) { 
                remodel_mover_->add_fragments( frag_libs[i]) ; 
            }
        }
        
        // if you have native structure to compare, do this
        if (get_native_pose()) remodel_mover_->set_native_pose(get_native_pose()) ;
        
        // scoring function, their default scoring function is the same as specified here, but put it anyway
        remodel_mover_->set_scorefxn( lowres_scorefxn_ );

    }
    
    
    /*  JQX: the following code is probably not ncessary*/
    
    Size unaligned_cdr_loop_begin(0), unaligned_cdr_loop_end(0);
    std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
    core::import_pose::pose_from_pdb( hfr_pose_, path+"hfr.pdb" );
    std::string cdr_name = "h3";
    AntibodyInfoOP hfr_info =  new AntibodyInfo ( hfr_pose_, cdr_name );
    unaligned_cdr_loop_begin = hfr_info->current_start;
    unaligned_cdr_loop_end   = hfr_info->current_end;
    
    if(framework_loop_size > 4){  //JQX: add this if statement to match R2_antibody
        pose_in.set_psi  (framework_loop_begin - 1, hfr_pose_.psi( unaligned_cdr_loop_begin - 1 )   );
        pose_in.set_omega(framework_loop_begin - 1, hfr_pose_.omega( unaligned_cdr_loop_begin - 1 ) );
    }

        //pose_in.dump_pdb("after_copying_nter.pdb");

    
    antibody2::AntibodyInfoOP starting_antibody;
    starting_antibody = new AntibodyInfo(*ab_info_);
    bool closed_cutpoints( false );
    
    
    Size cycle ( 1 );
    while( !closed_cutpoints && cycle < max_cycle_) {
        ab_info_ = starting_antibody;
        if (do_cter_insert_){
            h3_cter_insert_mover_->apply(pose_in);
        }
        
        //pose_in.dump_pdb("after_c_insert.pdb");

        if(remodel_=="legacy_perturb_ccd"){
            h3_perturb_ccd_build_->apply(pose_in);
        }
        else{
            remodel_mover_->apply(pose_in);
        }

        
        closed_cutpoints = cutpoints_separation( pose_in, ab_info_ );
        ++cycle;
    } // while( ( cut_separation > 1.9 )
    
    TR <<  "Finished Modeling Centroid CDR H3 loop" << std::endl;
    

    
    
    //#############################  //JQX: this should not be here
    if( is_camelid_ ){
        RefineCDRH1Centroid refine_cdr_centroid( ab_info_->get_CDR_loop("h1") );
        refine_cdr_centroid.apply(pose_in);
    }
    //#############################
    
    
    
    
    to_full_atom.apply( pose_in );

    utility::vector1<bool> allow_chi_copy( pose_in.total_residue(), true );
    for( Size ii=ab_info_->get_CDR_loop("h3")->start(); ii<=ab_info_->get_CDR_loop("h3")->stop(); ii++ ){
        allow_chi_copy[ii] = false;
    }
    //recover sidechains from starting structures except H3
    protocols::simple_moves::ReturnSidechainMover recover_sidechains( start_pose, allow_chi_copy );
    recover_sidechains.apply( pose_in );
    
    TR << "Finished applying CDR H3 modeler" << std::endl;
    

    return;
} // ModelCDRH3::apply()

    
    
    
std::string ModelCDRH3::get_name() const {
	return "ModelCDRH3";
}


//basic::Tracer & my_LoopMover::tr() const{
//    return TR;
//}






} // namespace antibody2
} // namespace protocols



