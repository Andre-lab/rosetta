// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/methods/util.cc
/// @brief  miscellaneous utility functions for forge
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Possu Huang (possu@u.washington.edu)

// unit headers
#include <protocols/forge/methods/util.hh>

// package headers
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>

// project headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/graph/DisjointSets.hh>
#include <core/id/types.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/loops/Loops.hh>

#include <core/pose/PDBPoseMap.hh>

// numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <set>

#include <utility/vector1.hh>


// REQUIRED FOR WINDOWS
#ifdef _WIN32
#include <cstdio>
#include <ctype.h>
#endif

using basic::T;

namespace protocols {
namespace forge {
namespace methods {


// static
static basic::Tracer TR( "protocols.forge.methods.util" );
static numeric::random::RandomGenerator RG( 2211990 ); // magic number, don't change


/// @brief perform union( root, i ) for all 'i' within the closed interval
///  [left, right]
/// @param[in] root position to union with; can be any number, does not have
///  to be a true root of a set
/// @param[in] left start of the interval
/// @param[in] right end of the interval
/// @param[in,out] uf
void
union_interval(
	core::Size const root,
	core::Size const left,
	core::Size const right,
	core::graph::DisjointSets & uf
)
{
	using core::Size;

	assert( left <= right );
	assert( root <= uf.n_nodes() );
	assert( left <= uf.n_nodes() );
	assert( right <= uf.n_nodes() );

	for ( Size i = left; i <= right; ++i ) {
		uf.ds_union( root, i );
	}
}

/// @brief moving left to right, find the first true cutpoint within specified
///  extent
/// @return the cutpoint position, otherwise 0 if not found
/// @details 0, pose.n_residue(), and any cutpoint with lower/upper terminus
///  are not counted as cutpoints
core::Size
find_cutpoint(
	core::pose::Pose const & pose,
	core::Size left,
	core::Size right
)
{
	using core::Size;
	using core::kinematics::FoldTree;

	FoldTree const & ft = pose.fold_tree();

	for ( Size i = left; i <= right; ++i ) {
		if ( ft.is_cutpoint( i ) && i < pose.n_residue() &&
		     !pose.residue( i ).is_lower_terminus() &&
		     !pose.residue( i ).is_upper_terminus()
		) {
			return i;
		}
	}

	return 0;
}


/// @brief moving left to right, count the number of true cutpoints within the
///  specified extent
/// @details 0, pose.n_residue(), and any cutpoint with lower/upper terminus
///  are not counted as cutpoints
core::Size
count_cutpoints(
	core::pose::Pose const & pose,
	core::Size left,
	core::Size right
)
{
	using core::Size;
	using core::kinematics::FoldTree;

	Size n = 0;
	FoldTree const & ft = pose.fold_tree();

	for ( Size i = left; i <= right; ++i ) {
		if ( ft.is_cutpoint( i ) && i < pose.n_residue() &&
		     !pose.residue( i ).is_lower_terminus() &&
		     !pose.residue( i ).is_upper_terminus()
		) {
			++n;
		}
	}

	return n;
}


/// @brief set omega to 180 for a range of residues [left, right]
void
trans_omega(
	core::Size const left,
	core::Size const right,
	core::pose::Pose & pose
)
{
	using core::Size;

	for ( Size i = left; i <= right; ++i ) {
		pose.set_omega( i, 180.0 );
	}
}


/// @brief create Loop object w/ random cutpoint from an Interval
protocols::loops::Loop
interval_to_loop( protocols::forge::build::Interval const & interval ) {
	using core::Size;
	using protocols::loops::Loop;

	// pick a cutpoint fully inside the loop so that there is always
	// at least one moveable residue on each side
	Size const cut = RG.random_range( interval.left, interval.right - 1 );

	return Loop( interval.left, interval.right, cut );
}


/// @brief create fold tree from loops
/// @remarks This is a generic replacement function for the one in protocols::loops
///  and will be moved there in the near future.
core::kinematics::FoldTree
fold_tree_from_loops(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops
)
{
	using core::Size;
	using core::kinematics::FoldTree;
	using core::kinematics::MoveMap;
	using protocols::loops::Loop;
	using protocols::loops::Loops;

	using protocols::forge::methods::fold_tree_from_pose;

	// setup movemap, mark only loops moveable, everything else fixed
	MoveMap mm;
	loops.switch_movemap( mm, core::id::BB , true );
	loops.switch_movemap( mm, core::id::CHI , true );

	// Generate initial fold tree from pose.  If existing root is a virtual
	// residue use it, otherwise use the existing root if non-moveable,
	// otherwise use a random non-moveable residue.
	Size root;
	if ( pose.residue( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		root = pose.fold_tree().root();
	} else if ( !mm.get_bb( pose.fold_tree().root() ) ) {
		root = pose.fold_tree().root();
	} else { // need to pick a random non-moveable residue

		utility::vector1< Size > fixed_bb;
		for ( Size i = 1, ie = pose.n_residue(); i <= ie; ++i ) {
			if ( !mm.get_bb( i ) ) {
				fixed_bb.push_back( i );
			}
		}

		if ( fixed_bb.size() == 0 ) {
			root = 1;
		} else {
			root = fixed_bb[ RG.random_range( 1, fixed_bb.size() ) ];
		}
	}

	FoldTree ft = fold_tree_from_pose( pose, root, mm );

	// track chain termini to distinguish internal vs terminal loops
	std::set< Size > lower_termini;
	std::set< Size > upper_termini;
	for ( Size i = 1, ie = pose.conformation().num_chains(); i <= ie; ++i ) {
		lower_termini.insert( pose.conformation().chain_begin( i ) );
		upper_termini.insert( pose.conformation().chain_end( i ) );
	}

	// post-modify tree with new loop jump/cuts for internal loop modeling
	for ( Loops::const_iterator i = loops.begin(), ie = loops.end(); i != ie; ++i ) {
		Loop const & loop = *i;

		if ( lower_termini.find( loop.start() ) == lower_termini.end() && upper_termini.find( loop.stop() ) == upper_termini.end() ) {
			ft.new_jump( loop.start() - 1, loop.stop() + 1, loop.cut() );
		}
	}

	return ft;
}


/// @brief set a single loop fold tree
/// @remarks This is a generic replacement function for the one in protocols::loops
///  and will be moved there in the near future.
void
set_single_loop_fold_tree(
	core::pose::Pose & pose,
	protocols::loops::Loop const & loop
)
{
	using core::kinematics::FoldTree;
	using protocols::loops::Loops;

	Loops loops;
	loops.add_loop( loop );

	FoldTree ft = fold_tree_from_loops( pose, loops );
	pose.fold_tree( ft );
}

	// duplicated code from resfile reader... unfortunately there's no easy way around.
	// question: if no chain is supplied should it be accepted?
	// yes just pass ' ' for the chain
	// how if a symbol is a chain or not?
	// all commands begin with something in the command map, if it's not a command treat it as a chain

utility::vector1< bool >
parse_resfile_string_with_no_lockdown( core::pose::Pose const & pose, core::pack::task::PackerTask & the_task, std::string const & resfile_string )// throw(ResfileReaderException)
{
	using namespace std;
	using namespace core::pack::task;
	using namespace core;
	istringstream resfile(resfile_string);

	bool have_read_start_token = false;

	utility::vector1< bool > non_default_lines( the_task.total_residue(), false );
	utility::vector1< std::string > default_tokens;
	utility::vector1< Size > origin_lines_of_default_tokens;

	core::uint lineno = 0;
	while ( resfile ) {
		map< string, ResfileCommandOP > command_map = create_command_map();
		utility::vector1< string > tokens( tokenize_line( resfile ));
		++lineno;

		// for debug
		//std::cout << "line->";
		//for( Size i=1; i <= tokens.size(); i++){
		//	std::cout << tokens[ i ] << ", ";
		//}
		//std::cout << std::endl;

		Size ntokens( tokens.size() );
		if ( ntokens == 0 ) continue;
		if ( comment_begin( tokens, 1 ) ) continue; // ignore the rest of this line

		if ( have_read_start_token ) {
			Size which_token = 1;

			// expected format: <residue identifier> <chain identifier> <commands*>
			//( the res/chain combo is used to get the pose's resid)

			// PDB numbering can be negative
			std::string const PDBnum_token = get_token( which_token, tokens );
			int PDBnum;
			char icode = ' ';
#ifdef _WIN32
			if ( isalpha( *PDBnum_token.rbegin() ) ) { // REQUIRED FOR WINDOWS
#else
			if ( std::isalpha( *PDBnum_token.rbegin() ) ) {
#endif
				PDBnum = atoi( PDBnum_token.substr( 0, PDBnum_token.length() - 1 ).c_str() );
				icode = *PDBnum_token.rbegin();
			} else { // no insertion code
				PDBnum = atoi( PDBnum_token.c_str() );
			}
			++which_token;
			char chain;
			chain = get_token( which_token, tokens )[ 0 ];
			if (chain == '_') chain = ' ';
			++which_token;

			Size resid(0);
			if(pose.pdb_info()){
				resid = pose.pdb_info()->pdb2pose().find( chain, PDBnum, icode );
			} else {
				if((1 <= PDBnum) && (PDBnum <= pose.total_residue())) {
					resid = PDBnum;
				}
			}
			if (resid == 0){
				std::stringstream err_msg;
				err_msg  << "On line " << lineno << ", the pose does not have residue (" << chain << ", " << PDBnum << ").";
				onError( err_msg.str());
			}
			non_default_lines[ resid ] = true;

			while ( which_token <= ntokens ) {
				if ( comment_begin( tokens, which_token ) ) break; // ignore the rest of this line
				if ( command_map.find( get_token( which_token, tokens ) ) == command_map.end() ) {
					std::stringstream err_msg;
					err_msg  << "On line " << lineno << " command '" << get_token( which_token, tokens) <<"' is not recognized.";
					onError(err_msg.str());
					which_token++;
					continue;
				}

				ResfileCommandOP command = command_map[ get_token( which_token, tokens ) ];

				try{
					command->initialize_from_tokens( tokens, which_token, resid );
					command->residue_action( the_task, resid );
				} catch ( ResfileReaderException() ){
					// there was a problem with this command.  If we're doing error recovery skip to next command.
					while( which_token <= ntokens && command_map.find( get_token(which_token, tokens ) ) == command_map.end() )
						which_token++;
					continue;
				}
			}

		} else { // the start token has not been read
			// read in default behaviors, store them, process them later
			if ( get_token( 1, tokens) == "START" ) {
				have_read_start_token = true;
			} else {
				for ( Size ii = 1; ii <= ntokens; ++ii ) {
					if ( comment_begin( tokens, ii ) ) break; // ignore the rest of this line
					default_tokens.push_back( get_token( ii, tokens ) );
					origin_lines_of_default_tokens.push_back( lineno );
				}
			}

		}
	}

	if ( ! have_read_start_token ) {
		T("core.pack.task.ResfileReader") << "RESFILE WARNING: reached the end of resfile without finding a 'start' token." << std::endl;
		T("core.pack.task.ResfileReader") << "RESFILE WARNING: No residue-specific behavior specified in resfile" << std::endl;
	}

	return non_default_lines;

	// now process default behaviors
/*
	for ( Size ii = 1; ii <= non_default_lines.size(); ++ii ) {
		if ( ! non_default_lines[ ii ] ) {
			Size which_token = 1, ntokens = default_tokens.size();

			while( which_token <= ntokens ){
				if ( command_map.find( get_token( which_token, default_tokens ) ) == command_map.end() ) {
					std::stringstream err_msg;
					err_msg  << "The default  command '" << get_token( which_token, default_tokens) <<"' is not recognized.";
					onError(err_msg.str());
					which_token++;
					continue;
				}

				ResfileCommandOP command = command_map[ get_token( which_token, default_tokens ) ];

				try{
					command->residue_action( default_tokens, which_token, the_task, ii );
				} catch ( ResfileReaderException() ){
					// there was a problem with this command.  If we're doing error recovery skip to next command.
					while( which_token <= ntokens && command_map.find( get_token(which_token, default_tokens ) ) == command_map.end() )
						which_token++;
					continue;
				}
			}
		}
	}
*/
}

core::pack::task::TaskFactoryOP
remodel_generic_taskfactory(){
  using core::pack::task::operation::IncludeCurrent;
  using core::pack::task::operation::InitializeFromCommandline;
  using core::pack::task::operation::NoRepackDisulfides;
  using protocols::toolbox::task_operations::LimitAromaChi2Operation;

  core::pack::task::TaskFactoryOP TF = new core::pack::task::TaskFactory();

  TF->push_back( new InitializeFromCommandline() ); // also inits -ex options
  TF->push_back( new IncludeCurrent() ); // enforce keeping of input sidechains
  TF->push_back( new NoRepackDisulfides() );
  if (!basic::options::option[basic::options::OptionKeys::remodel::design::allow_rare_aro_chi].user()){
  	TF->push_back( new LimitAromaChi2Operation() );
	}

	return TF;
}

void
fill_non_loop_cst_set(
  core::pose::Pose & pose,
	  protocols::loops::Loops loops)
{
	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring::constraints;

	core::Real const coord_sdev( 2.0 );
	core::Size const my_anchor(1);

	ConstraintSetOP cst_set = pose.constraint_set()->clone();

	std::set<core::Size> loopRange;

 	for ( protocols::loops::Loops::const_iterator it = loops.begin(), ite = loops.end(); it != ite; ++it ) {
     protocols::loops::Loop const & loop = *it;
			for (core::Size i = loop.start(); i<= loop.stop(); i++){

				loopRange.insert(i);

			}
	}
	core::Size const nres( pose.total_residue());
	for (core::Size i =1 ; i<=nres ; ++i){
		if (loopRange.find(i) != loopRange.end()){ //value exist(check!)
			continue;
		}
		else {
			Residue const & i_rsd( pose.residue(i));

			for (core::Size ii = 1; ii <= i_rsd.nheavyatoms(); ++ii) {
				cst_set->add_constraint( new CoordinateConstraint(AtomID(ii,i), AtomID(1, my_anchor), i_rsd.xyz(ii), new HarmonicFunc(0.0, coord_sdev)));
			}
		}
	}

	pose.constraint_set( cst_set );
}
	void fixH(core::pose::Pose & pose) {
		using namespace core;
		using core::id::AtomID;
		for(Size i = 1; i <= pose.n_residue(); ++i) {
			numeric::xyzVector<Real> n  = pose.residue(i).xyz("N");
			numeric::xyzVector<Real> ca = pose.residue(i).xyz("CA");
			Size in = i-1;
			if(in == 0) in = pose.n_residue();
			numeric::xyzVector<Real> c  = pose.residue(in).xyz("C");
			numeric::xyzVector<Real> h  = n + (n-(ca+c)/2.0).normalized()*1.01;
			if (pose.residue(i).name3() != "PRO"){
				pose.set_xyz(AtomID(pose.residue(i).atom_index("H"),i), h );
			}
		}
	}


void cyclize_pose(core::pose::Pose & pose) {
	using namespace core;
	using namespace core::pose;
	using namespace core::scoring::constraints;
	using namespace chemical;
	using core::id::AtomID;
	//pose.conformation().show_residue_connections();
  Size N = pose.n_residue();
  for(Size i = 1; i <= N; ++i) {
    if(pose.residue(i).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(pose,i);
    if(pose.residue(i).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(pose,i);
    if(pose.residue(i).has_variant_type(CUTPOINT_UPPER)) core::pose::remove_variant_type_from_pose_residue(pose,CUTPOINT_UPPER,i);
    if(pose.residue(i).has_variant_type(CUTPOINT_LOWER)) core::pose::remove_variant_type_from_pose_residue(pose,CUTPOINT_LOWER,i);
  }
  if(!pose.residue(1).has_variant_type(CUTPOINT_UPPER)) {
		core::pose::add_variant_type_to_pose_residue(pose,CUTPOINT_UPPER,1);
		/*
		//make connection
		conformation::ResidueOP cloneRes = new conformation::Residue(*pose.residue(1).clone());

		cloneRes->residue_connection_partner(1, N, 2);
		cloneRes->residue_connection_partner(2, 2, 1);

		pose.replace_residue(1, *cloneRes, false);
		*/

	}
  if(!pose.residue(N).has_variant_type(CUTPOINT_LOWER)) {
		core::pose::add_variant_type_to_pose_residue(pose,CUTPOINT_LOWER,N);
		/*
		//make connection
		conformation::ResidueOP cloneRes = new conformation::Residue(*pose.residue(N).clone());

		cloneRes->residue_connection_partner(1, N-1, 2);
		cloneRes->residue_connection_partner(2, 1, 1);

		pose.replace_residue(N, *cloneRes, false);
		*/
	}
  pose.conformation().declare_chemical_bond( 1, "N", N, "C" );
  fixH(pose);

	pose.conformation().update_polymeric_connection(1);

//	conformation::ResidueOP cloneRes2 = new conformation::Residue(*pose.residue(2).clone());
//	cloneRes2->residue_connection_partner(1, 1, 2);
//	cloneRes2->residue_connection_partner(2, 3, 1);
//std::cout << "resconn1CR2: " << cloneRes2->connected_residue_at_resconn( 1 ) << " ";
//        std::cout << cloneRes2->residue_connection_conn_id(1);
//        std::cout << " resconn2: " << cloneRes2->connected_residue_at_resconn( 2 ) << " ";
//        std::cout << cloneRes2->residue_connection_conn_id(2);
//        std::cout << " seqpos: " << cloneRes2->seqpos() << std::endl;
//
//	pose.replace_residue(2, *cloneRes2, false);
//	pose.conformation().update_polymeric_connection(2);
//std::cout << "resconn1PR2: " << pose.residue(2).connected_residue_at_resconn( 1 ) << " ";
//         std::cout << pose.residue(2).residue_connection_conn_id(1);
//         std::cout << " resconn2: " << pose.residue(2).connected_residue_at_resconn( 2 ) << " ";
//         std::cout << pose.residue(2).residue_connection_conn_id(2);
//         std::cout << " seqpos: " << pose.residue(2).seqpos() << std::endl;

	/*
	//make connection
	conformation::ResidueOP cloneRes1 = new conformation::Residue(*pose.residue(1).clone());
//	conformation::ResidueOP cloneResN = new conformation::Residue(*pose.residue(N).clone());

	cloneRes1->residue_connection_partner(1, N, 2);
	cloneRes1->residue_connection_partner(2, 2, 1);
//	cloneRes1->residue_connection_partner(2, N, 2);
	cloneRes2->residue_connection_partner(1, 1, 2);
	cloneRes2->residue_connection_partner(2, 3, 1);
//	cloneResN->residue_connection_partner(1, N-1, 2);
//	cloneResN->residue_connection_partner(2, 1, 1);

	pose.replace_residue(1, *cloneRes1, false);
	pose.replace_residue(2, *cloneRes2, false);
//	pose.replace_residue(N, *cloneResN, false);
*/

//	pose.conformation().show_residue_connections();

  using namespace core::scoring::constraints;
  AtomID a1( pose.residue(1).atom_index(   "N"), 1 ), a2( pose.residue(pose.n_residue()).atom_index("OVL1"), pose.n_residue() );
  AtomID b1( pose.residue(1).atom_index(  "CA"), 1 ), b2( pose.residue(pose.n_residue()).atom_index("OVL2"), pose.n_residue() );
  AtomID c1( pose.residue(1).atom_index("OVU1"), 1 ), c2( pose.residue(pose.n_residue()).atom_index(   "C"), pose.n_residue() );
//  pose.remove_constraints();
  pose.add_constraint(new AtomPairConstraint(a1,a2,new HarmonicFunc(0.0,0.1)));
  pose.add_constraint(new AtomPairConstraint(b1,b2,new HarmonicFunc(0.0,0.1)));
  pose.add_constraint(new AtomPairConstraint(c1,c2,new HarmonicFunc(0.0,0.1)));
}




//void lockdown_task_with_default_behaviors(){
//}



} // namespace methods
} // namespace forge
} // namespace protocols
