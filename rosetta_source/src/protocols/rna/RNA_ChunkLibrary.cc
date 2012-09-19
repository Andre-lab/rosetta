 // -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
 // vi: set ts=2 noet:
 //  CVS information:
 //  $Revision: 1.1.2.1 $
 //  $Date: 2005/11/07 21:05:35 $
 //  $Author: rhiju $
 // (c) Copyright Rosetta Commons Member Institutions.
 // (c) This file is part of the Rosetta software suite and is made available under license.
 // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
 // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
 // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Rosetta Headers
#include <protocols/rna/RNA_ChunkLibrary.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/MiniPose.fwd.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/rms_util.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <numeric/random/random.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
// AUTO-REMOVED #include <core/kinematics/AtomTree.hh>
// AUTO-REMOVED #include <core/kinematics/tree/Atom.hh>
#include <core/conformation/Residue.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Numeric headers
// AUTO-REMOVED #include <numeric/constants.hh>

// C++ headers
// AUTO-REMOVED #include <fstream>
#include <iostream>

#include <utility/vector1.hh>



static numeric::random::RandomGenerator RG(2380934);  // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.rna.rna_chunk_library" ) ;

namespace protocols{
namespace rna{

	using namespace core;
	using namespace ObjexxFCL;

	using core::Size;
	using core::Real;

	using core::pose::ResMap;

	///////////////////////////////////////////////////////////////////////
	ChunkSet::ChunkSet( utility::vector1< core::pose::MiniPoseOP > const & mini_pose_list,
											ResMap const & res_map ) {

		mini_pose_list_ = mini_pose_list;

		res_map_ = res_map;

		// not much information in mini_pose --> assume that all atoms are OK for copying.
		core::pose::MiniPose const & mini_pose = *(mini_pose_list[ 1 ]);
		for ( Size i = 1; i <= mini_pose.total_residue(); i++ ){
			for ( Size j = 1; j <= mini_pose.coords()[i].size(); j++ ){
				atom_id_mask_[ core::id::AtomID( j, i ) ] = true;
			}
		}

	}

	///////////////////////////////////////////////////////////////////////
	ChunkSet::ChunkSet( utility::vector1< core::pose::PoseOP > const & pose_list,
											ResMap const & res_map ) {
		for ( Size n = 1; n <= pose_list.size(); n++ ) {
			mini_pose_list_.push_back( core::pose::MiniPoseOP( new core::pose::MiniPose( *(pose_list[n]) ) ) );
		}


		res_map_ = res_map;

		core::pose::Pose const & pose = *( pose_list[1] );
		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			core::conformation::Residue rsd = pose.residue( i );
			for ( Size j = 1; j <= rsd.natoms(); j++ ){
				atom_id_mask_[ core::id::AtomID( j, i ) ] = !rsd.is_virtual( j );
			}
		}


	}

	///////////////////////////////////////////////////////////////////////
	ChunkSet::~ChunkSet()	{}


	///////////////////////////////////////////////////////////////////////
	void
	ChunkSet::insert_chunk_into_pose( core::pose::Pose & pose, Size const & chunk_pose_index,toolbox::AllowInsertOP const & allow_insert ) const{

		using namespace core::pose;
		using namespace core::id;

		core::pose::MiniPose const & scratch_pose ( *(mini_pose_list_[ chunk_pose_index ]) );

		//		TR << "SCRATCH_POSE " << scratch_pose.sequence() << ' ' << scratch_pose.fold_tree() << std::endl;

		std::map< AtomID, AtomID > atom_id_map = get_atom_id_map( pose, allow_insert );

		copy_dofs( pose, scratch_pose, atom_id_map  );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	std::map< id::AtomID, id::AtomID >
	ChunkSet::get_atom_id_map(  core::pose::Pose & pose, toolbox::AllowInsertOP const & allow_insert ) const{

		std::map< id::AtomID, id::AtomID > atom_id_map;

		allow_insert->calculate_atom_id_map( pose, res_map_, mini_pose_list_[1]->fold_tree(), atom_id_map );

		// This should prevent copying dofs for virtual phosphates, if they are tagged as such in the input silent files.
		filter_atom_id_map_with_mask( atom_id_map );

		return atom_id_map;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	ChunkSet::filter_atom_id_map_with_mask( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ) const{

		using namespace core::id;

		std::map< AtomID, AtomID > atom_id_map_new;

		for ( std::map< AtomID, AtomID >::const_iterator
						it=atom_id_map.begin(), it_end = atom_id_map.end(); it != it_end; ++it ) {

			AtomID const & insert_atom_id = it->first;
			AtomID const & source_atom_id = it->second;

			std::map< AtomID, bool >::const_iterator it_mask = atom_id_mask_.find( source_atom_id );
			if ( it_mask == atom_id_mask_.end() ) utility_exit_with_message( "Some problem with atom_id_mask in defining atom_id_map " );
			if ( !it_mask->second ) continue; // this source_atom_id is not allowed by mask, probably came from a virtual phosphate.

			atom_id_map_new[ insert_atom_id ] = source_atom_id;
		}

		atom_id_map = atom_id_map_new;

	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	core::pose::MiniPoseOP const
	ChunkSet::mini_pose( Size const idx ) const {
		return mini_pose_list_[ idx ];
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////
	RNA_ChunkLibrary::RNA_ChunkLibrary(){
		// currently nothing.
		chunk_coverage_ = 0.0;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	// constructor -- needs a list of silent files. Each silent file
	//  has solutions for a particular piece of the desired pose.
	//  THIS SHOULD BE DEPRECATED SOON -- no longer in use?
	//  Better to explicitly specify '-chunk_res'.
	RNA_ChunkLibrary::RNA_ChunkLibrary(
								utility::vector1 < std::string > const & silent_files,
								core::pose::Pose const & pose,
								std::map< Size, Size > const & connections_in_big_pose /* to figure out mapping to big pose*/ )
	{

		std::string const & sequence_of_big_pose( pose.sequence() );
		coarse_rna_ = pose.residue( 1 ).is_coarse();

		// allow_insert keeps track of where chunks are placed -- only allow
		// fragment insertions *outside* these regions.
		allow_insert_ = new toolbox::AllowInsert( pose );
		covered_by_chunk_.dimension( sequence_of_big_pose.size(), false );

		utility::vector1< Size > input_res;
		Size chunk_res_count( 0 );

		for ( Size n = 1; n <= silent_files.size(); n++ ) {

			utility::vector1< pose::PoseOP > pose_list;
			process_input_file( silent_files[n], pose_list );

			core::pose::Pose const & scratch_pose( *(pose_list[1]) );

			utility::vector1< ResMap > res_maps;

			// There may be more than one part of the pose to which this sequence maps.
			figure_out_possible_res_maps( res_maps, scratch_pose, sequence_of_big_pose, connections_in_big_pose );

			for (Size k = 1; k <= res_maps.size(); k++ )  {
				check_res_map( res_maps[ k ], *(pose_list[1]), sequence_of_big_pose );

				ChunkSetOP chunk_set( new ChunkSet( pose_list, res_maps[ k ] ) );
				chunk_sets_.push_back( chunk_set );

				zero_out_allow_insert( res_maps[ k ], pose, scratch_pose, n );
			}

			for ( ResMap::const_iterator
							it=res_maps[1].begin(), it_end = res_maps[1].end(); it != it_end; ++it ) {
				input_res.push_back( it->first );
			}

		}

		figure_out_chunk_coverage();

		//std::cout << "INPUT_RES: ";
		//		for ( Size n = 1; n <= input_res.size(); n++ ) std::cout << ' ' << input_res[ n ];
		//		std::cout << std::endl;

	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	// deprecate soon?
	RNA_ChunkLibrary::RNA_ChunkLibrary(
								utility::vector1 < std::string > const & silent_files,
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & input_res )
	{
		utility::vector1< std::string > pdb_files_BLANK;
		initialize_rna_chunk_library( pdb_files_BLANK, silent_files, pose, input_res );
	}


	//////////////////////////////////////////////////////////////////////////////////////////////
	// constructor -- needs a list of silent files. Each silent file
	//  has solutions for a particular piece of the desired pose.
	RNA_ChunkLibrary::RNA_ChunkLibrary(
								utility::vector1 < std::string > const & pdb_files,
								utility::vector1 < std::string > const & silent_files,
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & input_res )
	{
		initialize_rna_chunk_library( pdb_files, silent_files, pose, input_res );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::initialize_rna_chunk_library(
								utility::vector1 < std::string > const & pdb_files,
								utility::vector1 < std::string > const & silent_files,
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & input_res )
	{
		std::string const & sequence_of_big_pose( pose.sequence() );
		coarse_rna_ = pose.residue( 1 ).is_coarse();

		// allow_insert keeps track of where chunks are placed -- only allow
		// fragment insertions *outside* these regions.
		allow_insert_ = new toolbox::AllowInsert( pose );
		covered_by_chunk_.dimension( sequence_of_big_pose.size(), false );

		utility::vector1< std::string > all_input_files;
		utility::vector1< bool > is_pdb_file;
		for ( Size n = 1; n <= pdb_files.size(); n++ ){
			all_input_files.push_back( pdb_files[n] );
			is_pdb_file.push_back( true );
		}
		for ( Size n = 1; n <= silent_files.size(); n++ ){
			all_input_files.push_back( silent_files[n] );
			is_pdb_file.push_back( false );
		}

		Size count( 0 );
		for ( Size n = 1; n <= all_input_files.size(); n++ ) {

			utility::vector1< pose::PoseOP > pose_list;
			process_input_file( all_input_files[n], pose_list,  is_pdb_file[n] );

			core::pose::Pose const & scratch_pose( *(pose_list[1]) );

			// There may be more than one part of the pose to which this sequence maps.
			ResMap res_map;

			for ( Size i = 1; i <= scratch_pose.sequence().size(); i++ ) {
				count++;
				if ( sequence_of_big_pose[ input_res[ count ] -1 ] != scratch_pose.sequence()[ i - 1 ] ){
					std::cout << "Problem with input_file: " << all_input_files[n] << std::endl;
					std::cout << "mismatch in sequence   in  big pose: " << sequence_of_big_pose[ input_res[ count ] -1 ] << input_res[count] <<
						"  in input pose: " << scratch_pose.sequence()[ i - 1 ]  << i << std::endl;
					utility_exit_with_message( "mismatch in input_res sequence" );
				}
				res_map[ input_res[count ] ] = i;
			}

			ChunkSetOP chunk_set( new ChunkSet( pose_list, res_map ) );
			chunk_sets_.push_back( chunk_set );

			zero_out_allow_insert( res_map, pose, scratch_pose, n );

			//check_fold_tree_OK( res_map, pose, scratch_pose );

		}
		if ( count != input_res.size() ){
			utility_exit_with_message( "Number of input res does not match total res in input silent files!" );
		}

		figure_out_chunk_coverage();

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::add_chunk_set(
												std::string const & silent_file,
												ResMap const & res_map,
												pose::Pose const & big_pose )
	{

		utility::vector1< pose::PoseOP > pose_list;

		process_input_file( silent_file, pose_list );
		check_res_map( res_map, *(pose_list[1]), big_pose.sequence() );

		ChunkSetOP chunk_set( new ChunkSet( pose_list, res_map ) );
		chunk_sets_.push_back( chunk_set );

	}

	////////////////////////////////////////////////////////////////////////////
	void RNA_ChunkLibrary::insert_chunk_into_pose(
					 pose::Pose & pose,
					 Size const & chunk_list_index,
					 Size const & chunk_pose_index ) const
	{
		chunk_sets_[ chunk_list_index ]->insert_chunk_into_pose( pose, chunk_pose_index, allow_insert_ );
	}


	//////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChunkLibrary::random_chunk_insertion( core::pose::Pose & pose ) const{

		Size const chunk_set_index = static_cast <int> ( RG.uniform() * num_chunk_sets() ) + 1;

		ChunkSet const & chunk_set( *chunk_sets_[ chunk_set_index ] );

		if ( chunk_set.num_chunks() < 2 )  return false;

		Size const chunk_index = static_cast <int> ( RG.uniform() * chunk_set.num_chunks() ) + 1;

		chunk_set.insert_chunk_into_pose( pose, chunk_index, allow_insert_ );

		//		TR << "INSERTED CHUNK " << chunk_index << " FROM SET " << chunk_set_index << std::endl;

		return true;
	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::zero_out_allow_insert( ResMap const & res_map,
																					 core::pose::Pose const & pose,
																					 core::pose::Pose const & scratch_pose,
																					 core::Size const domain_num )
	{
		using namespace core::id;
		using namespace core::conformation;

		// connected doesn't do anything anymore...
		FArray1D< bool > connected( pose.total_residue(), false );
		covered_by_chunk_ = false;

		for ( ResMap::const_iterator
						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

			Size const i = it->first; //Index in big pose.
			Size const i_scratch = it->second; //Index in scratch pose (chunk).

			covered_by_chunk_( i ) = true;

			Residue const & rsd_i = pose.residue(i);
			for ( Size j = 1; j <= rsd_i.natoms(); j++ ){

				std::string const & atomname = rsd_i.atom_name( j );
				Residue const & scratch_rsd = scratch_pose.residue(i_scratch);

				if ( scratch_rsd.has( atomname ) ) {
					Size const & scratch_index = scratch_pose.residue( i_scratch ).atom_index( atomname );
					if ( !scratch_rsd.is_virtual( scratch_index ) ) {
						allow_insert_->set_domain( AtomID(j,i), domain_num);
					}
				}
			}

			//We don't trust phosphates at the beginning of chains!
			// MAKE THIS AN OPTION?
			//			if ( i_scratch == 1 || scratch_pose.fold_tree().is_cutpoint( i_scratch - 1 ) ) allow_insert_->set_phosphate( i, pose, true );

		}

	}


	//////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChunkLibrary::check_fold_tree_OK( pose::Pose const & pose ){

		for (Size k = 1; k <= chunk_sets_.size(); k++ )  {
			ChunkSet & chunk_set = *(chunk_sets_[k]);
			bool const OK = chunk_set.check_fold_tree_OK( pose );
			if (!OK){
				std::cout << "Problem with pose fold tree -- not enough jumps to handle the number of chains in chunk set " << k << std::endl;
				utility_exit_with_message( "FoldTree in pose does not have the right number of jumps to match chunk_res" );
			}
		}

		// A return value is required. Please fix this.
#ifdef WIN32
		return true;
#endif
	}


	//////////////////////////////////////////////////////////////////////////////
	bool
	ChunkSet::check_fold_tree_OK( pose::Pose const & pose ){

		// Check where the chunk is mapped to in the big pose.
		// There should be at least the same number of jumps in the big pose
		//  as there are chains in the scratch_pose.
		utility::vector1< bool > is_chunk_res( pose.total_residue(), false );
		for ( ResMap::const_iterator
						it=res_map_.begin(), it_end = res_map_.end(); it != it_end; ++it ) {
			Size const i = it->first; //Index in big pose.
			is_chunk_res[ i ] = true;
		}

		Size const num_jumps_scratch = mini_pose_list_[1]->fold_tree().num_jump(); // number of chains - 1

		Size num_jumps_in_big_pose_in_scratch_region( 0 );
		for ( Size n = 1; n <= pose.num_jump(); n++ ) {
			if (! is_chunk_res[ pose.fold_tree().upstream_jump_residue( n ) ] ) continue;
			if (! is_chunk_res[ pose.fold_tree().downstream_jump_residue( n ) ] ) continue;
			num_jumps_in_big_pose_in_scratch_region++;
		}

		if ( num_jumps_scratch > num_jumps_in_big_pose_in_scratch_region ){
			std::cout << "Number of jumps in chunk pose               : " << num_jumps_scratch << std::endl;
			std::cout << "Number of jumps in full pose in chunk region: " << num_jumps_in_big_pose_in_scratch_region  << "  out of total jumps " << pose.num_jump() << std::endl;
			return false;
		}

		if ( num_jumps_scratch < num_jumps_in_big_pose_in_scratch_region ){
			std::cout << "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!" << std::endl;
			std::cout << "Number of jumps in chunk pose               : " << num_jumps_scratch << std::endl;
			std::cout << "Does not match:" << std::endl;
			std::cout << "Number of jumps in full pose in chunk region: " << num_jumps_in_big_pose_in_scratch_region  << "  out of total jumps " << pose.num_jump() << std::endl;
			std::cout << "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!" << std::endl;
			// Just a warning
			//return false;
		}

		return true;

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::figure_out_chunk_coverage()
	{

		Size const tot_res( allow_insert_->nres() );
		Size num_chunk_res( 0 );
		Size num_other_res( 0 );

		for ( Size n = 1; n <= tot_res; n++ ) {
			// Allow insert keeps track of where the chunk *aren't*, and
			// where other moves (fragments, jumps) can be carried out.
			if ( covered_by_chunk_(n) ){
				num_chunk_res++;
			} else {
				num_other_res++;
			}
		}
		chunk_coverage_ = Real( 3 * num_chunk_res ) / ( 3 * num_chunk_res +  tot_res );

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::get_component_sequences(
							 core::pose::Pose const & pose,
							 utility::vector1< std::string > & sequences,
							 utility::vector1< core::Size > & chain_id,
							 utility::vector1< core::Size > & sequence_start ) const{

		chain_id.clear();
		sequences.clear();

		std::string sequence = "";
		Size count( 1 );
		sequence_start.push_back( 1 );

		for ( Size i = 1; i <= pose.total_residue(); i++ ) {

			sequence += pose.residue(i).name1();
			chain_id.push_back( count );

			if ( pose.fold_tree().is_cutpoint( i ) ) {
				sequences.push_back( sequence );
				sequence = "";
				count++;
				if ( i < pose.total_residue() ) sequence_start.push_back( i+1 );
			}

		}


		for (Size n = 1; n <= sequences.size(); n++ ) {
			TR << "SEQUENCE " << n << " " << sequences[ n ] << std::endl;
		}

	}

	//////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::figure_out_possible_res_maps(
					utility::vector1< ResMap > & res_maps,
					pose::Pose const & scratch_pose,
					std::string const & sequence_of_big_pose,
					std::map< Size, Size > const & connections_in_big_pose ) const
	{

		// Note -- using zero-indexed vectors -- easier to do modulo, etc.
		utility::vector1< std::string > scratch_sequences;
		utility::vector1< core::Size >  chain_id;
		utility::vector1< core::Size >  scratch_sequence_start;

		get_component_sequences( scratch_pose, scratch_sequences, chain_id, scratch_sequence_start );

		// Go through each sequence and look for matches
		utility::vector1< utility::vector1< Size > > matches_to_each_scratch_sequence;
		get_sequence_matches( matches_to_each_scratch_sequence, scratch_sequences, sequence_of_big_pose );

		// now actually find some good res maps.
		find_res_maps( chain_id, scratch_sequence_start, scratch_sequences, matches_to_each_scratch_sequence, scratch_pose, connections_in_big_pose, res_maps );
	}

	////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::find_res_maps(
								utility::vector1< Size > const & chain_id,
								utility::vector1< Size > const & scratch_sequence_start,
								utility::vector1< std::string > const & scratch_sequences,
								utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
								core::pose::Pose const & scratch_pose,
								std::map< Size, Size > const & connections_in_big_pose,
								utility::vector1< ResMap > & res_maps ) const
	{
		res_maps.clear();

		// Loop over matches for chain 1 -- if there are any chunks that match, they should be taggable by
		// their match to chain 1.

		Size num_chain( 1 );
		for (Size k = 1; k <= matches_to_each_scratch_sequence[ num_chain ].size(); k++ ) {
			ResMap res_map;
			fill_res_map( res_map, matches_to_each_scratch_sequence[ num_chain ][ k ],
										scratch_sequence_start[ num_chain ] /*should be 1*/,
										scratch_sequences[ num_chain ].size()  );

			// Dig deep into each connection from this chain.
			check_connections( num_chain, res_map,
												 chain_id, scratch_sequence_start, scratch_sequences, matches_to_each_scratch_sequence, scratch_pose, connections_in_big_pose, res_maps );
		}

		TR << "Number of matches found:  " << res_maps.size() << std::endl;
		if ( res_maps.size() == 0 )  utility_exit_with_message(  "Could not match silent file with sequence "+scratch_pose.sequence() );

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::get_sequence_matches( 	utility::vector1< utility::vector1< Size > > & matches_to_each_scratch_sequence,
																						utility::vector1< std::string > const &	scratch_sequences,
																						std::string const & sequence_of_big_pose ) const
	{
		//		Size tot_matches = 1;
		for ( Size n = 1; n <= scratch_sequences.size(); n++ ) {
			utility::vector1< Size > matches;
			std::string const scratch_sequence( scratch_sequences[n] );
			Size const scratch_sequence_length = scratch_sequence.size();
			for ( Size i = 0; i <= sequence_of_big_pose.size() - scratch_sequence_length; i++ ) {
				bool does_it_match( true );
				for (Size offset = 0; offset < scratch_sequence_length; offset++ ) {
					if ( sequence_of_big_pose[ i + offset ] != scratch_sequence[ offset ] ) {
						does_it_match = false;
						break;
					}
				}
				if (does_it_match) {
					matches.push_back( i+1 ); // convert numbering to start with 1
					TR << "Found match to scratch_sequence " << n <<
						//						"  starting at scratch pose position " << scratch_sequence_start[ n ] << " " <<
						"   at big pose position: " << i+1 << std::endl;
				}
			}

			matches_to_each_scratch_sequence.push_back( matches );

			if ( matches.size() < 1 ) 	 utility_exit_with_message(  "Could not find match to sequence" );

		}
}
	////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::check_connections( Size const & num_chain, ResMap & res_map,
																			 utility::vector1< Size > const & chain_id,
																			 utility::vector1< Size > const & scratch_sequence_start,
																			 utility::vector1< std::string > const & scratch_sequences,
																			 utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
																			 core::pose::Pose const & scratch_pose,
																			 std::map< Size, Size > const & connections_in_big_pose,
																			 utility::vector1< ResMap > & res_maps ) const
	{
		// We might already be done!
		if ( res_map.size() == scratch_pose.total_residue() ) { //everything assigned!
			res_maps.push_back( res_map );
			return;
		}

		// Look for connections coming off the current chain
		Size const n_jump( scratch_pose.num_jump() ) ;

		for ( Size n = 1; n <= n_jump; n++ ) {

			//			TR << "CHECK OUT JUMP  " << n << " for chain " << num_chain << std::endl;
			Size const res1 =  scratch_pose.fold_tree().upstream_jump_residue( n );
			Size const res2 =  scratch_pose.fold_tree().downstream_jump_residue( n );

			if ( chain_id[ res1 ] == num_chain ) {
				test_matches( res1, res2, res_map,
											chain_id, scratch_sequence_start,
											scratch_sequences, matches_to_each_scratch_sequence,
											scratch_pose, connections_in_big_pose, res_maps );
			} else if ( chain_id[ res2 ] == num_chain ) {
				test_matches( res2, res1, res_map,
											chain_id, scratch_sequence_start,
											scratch_sequences, matches_to_each_scratch_sequence,
											scratch_pose, connections_in_big_pose, res_maps );
			}

		}

	}

	//////////////////////////////////////////////
	void
	RNA_ChunkLibrary::test_matches( Size const & res1, Size const & res2, ResMap & res_map,
																	utility::vector1< Size > const & chain_id,
																	utility::vector1< Size > const & scratch_sequence_start,
																	utility::vector1< std::string > const & scratch_sequences,
																	utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
																	core::pose::Pose const & scratch_pose,
																	std::map< Size, Size > const & connections_in_big_pose,
																	utility::vector1< ResMap > & res_maps ) const {

		//		TR << " TESTING : " << res1 << " " << res2 << "  from chains: " << chain_id[ res1 ] << " " << chain_id[ res2 ] <<
		//			". Length of resmap " << res_map.size() << " out of " << scratch_pose.total_residue() <<  std::endl;

		Size const next_chain = chain_id[ res2 ];
		// res2 may already be "taken care of" -- no need to check it out.
		for ( ResMap::const_iterator
						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {
			if ( it->second == res2 ) {
				//				TR << " ALREADY TESTED " << std::endl;
				return;
			}
		}

		// Cycle through potential matches -- are any good?
		for (Size k = 1; k <= matches_to_each_scratch_sequence[ next_chain ].size(); k++ ) {

			ResMap res_map_test( res_map );

			bool const res_map_ok = fill_res_map( res_map_test, matches_to_each_scratch_sequence[ next_chain ][ k ],
																						scratch_sequence_start[ next_chain ],scratch_sequences[next_chain].size()  );

			if (!res_map_ok) continue;

			// Are the res1 and res2 connected in the big pose?
			// Need to figure out where they map to in big pose.
			Size res1_map( 0 ), res2_map( 0 );
			for ( ResMap::const_iterator
							it=res_map_test.begin(), it_end = res_map_test.end(); it != it_end; ++it ) {
				if ( it->second == res1 ) res1_map = it->first;
				if ( it->second == res2 ) res2_map = it->first;
			}
			if (res1_map == 0 || res2_map == 0 ){
				TR << res1 << " " << res1_map << "        " << res2 << " " << res2_map << std::endl;
				utility_exit_with_message( "SHOULD NOT BE HERE! " );
			}


			bool connection_ok( false );
			for ( ResMap::const_iterator
							it=connections_in_big_pose.begin(), it_end = connections_in_big_pose.end(); it != it_end; ++it ) {

				Size const res1_in_big_pose = it->first;
				Size const res2_in_big_pose = it->second;

				if ( res1_in_big_pose == res1_map &&
						 res2_in_big_pose == res2_map ) {
					connection_ok = true; break;
				}

				if ( res2_in_big_pose == res1_map &&
						 res1_in_big_pose == res2_map ) {
					connection_ok = true; break;
				}
			}

			if (!connection_ok) {
				//				TR << "DENIED!!  " << std::endl;
				continue;
			} else {
				//				TR << "OK!!  " << std::endl;
				res_map = res_map_test;
				check_connections( next_chain, res_map,
													 chain_id, scratch_sequence_start,
													 scratch_sequences, matches_to_each_scratch_sequence,
													 scratch_pose, connections_in_big_pose, res_maps );
				break;
			}
		}


	}


	///////////////////////////////////////////////////////////
	bool
	RNA_ChunkLibrary::fill_res_map( ResMap & res_map, Size const & match_pos, Size const & scratch_start_pos, Size const & scratch_sequence_length ) const
	{
		bool one_to_one( true );
		for (Size offset = 0; offset < scratch_sequence_length; offset++ ) {
			Size const big_pose_pos = match_pos + offset;
			Size const scratch_pos = scratch_start_pos + offset;
			if ( res_map.find( big_pose_pos) == res_map.end() ) {
				res_map[ big_pose_pos ] = scratch_pos;
				TR << "MAPPING " << match_pos+offset << " --> " << scratch_start_pos+offset << std::endl;
			} else {
				// this is a problem -- expect res_maps to be one-to-one.
				one_to_one = false;
				break;
			}
		}
		return one_to_one;
	}



	///////////////////
	// DELETE FOLLOWING AFTER STUFF WORKS!
	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::check_res_map_recursively( ResMap const & res_map_old,
																							 utility::vector1< std::string > const & scratch_sequences,
																							 utility::vector1< utility::vector1< Size > > const & matches_to_each_scratch_sequence,
																							 pose::Pose const & scratch_pose,
																							 std::map< Size, Size > const & connections_in_big_pose,
																							 utility::vector1< core::Size > const & chain_id,
																							 Size const & num_sequence, Size const & num_match, utility::vector1< ResMap > & res_maps ) const{
		ResMap res_map( res_map_old );

		Size const match_pos = matches_to_each_scratch_sequence[ num_sequence ][ num_match ];
		//		TR << "HEY MATCH_POS " << match_pos << std::endl;
		Size const scratch_sequence_length = scratch_sequences[num_sequence].size();
		Size i( res_map.size() );
		for (Size offset = 0; offset < scratch_sequence_length; offset++ ) {
			res_map[ match_pos + offset + 1 ] = i + 1; // Add back in 1 to match pose numbering.
			i++;
		}

		bool const jump_match = check_jump_match( scratch_pose, connections_in_big_pose, res_map, chain_id );

		//		for ( ResMap::const_iterator
		//						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {
		//			TR << it->first << " mapped to " << it->second << std::endl;
		//		}

		//		TR << " CHECKING: " << num_sequence << " -- " << num_match << " : " << res_map.size() << " --> " << jump_match << std::endl;

		if ( !jump_match ) return;

		if ( num_sequence == matches_to_each_scratch_sequence.size() ) { //Success!
			res_maps.push_back( res_map );
		}

		//		TR << "ADVANCING to sequence " << num_sequence << std::endl;
		Size const num_sequence_next = num_sequence + 1;

		if ( num_sequence_next <= matches_to_each_scratch_sequence.size() ) {

			for ( Size k = 1; k <= matches_to_each_scratch_sequence[ num_sequence_next ].size(); k++ ) {
				//				TR << "ABOUT TO TRY  " << num_sequence_next << " " << k << std::endl;
				check_res_map_recursively( res_map, scratch_sequences, matches_to_each_scratch_sequence, scratch_pose, connections_in_big_pose, chain_id, num_sequence_next, k, res_maps );
			}

		}

		return;

	}



	///////////////////////////////////////////////////////////////////////////
	// Any jump inside the scratch pose (which came from a user-inputted silent file with its own fold tree)
	// must correspond to a pairing defined in the new pose of interest. Could make this a little less restrictive,
	// just asking for any connection between the chain segments that are putatively matching.
	bool
	RNA_ChunkLibrary::check_jump_match(
									 pose::Pose const & scratch_pose,
									 std::map< Size, Size > const & connections_in_big_pose,
									 ResMap const & res_map,
									 utility::vector1< Size > const & chain_id ) const
	{

		Size const n_jump( scratch_pose.num_jump() ) ;

		for ( Size n = 1; n <= n_jump; n++ ) {

			Size const res1 =  scratch_pose.fold_tree().upstream_jump_residue( n );
			Size const res2 =  scratch_pose.fold_tree().downstream_jump_residue( n );

			// First check that res1 and res2 are in the resmap -- otherwise no point in looking for a matching
			// jump in the big pose.
			bool found_res1( false ), found_res2( false );
			for ( ResMap::const_iterator it = res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {
				if ( res1 == it->second ) found_res1 = true;
				if ( res2 == it->second ) found_res2 = true;
			}
			if ( !found_res1 || !found_res2 ) {
				//consider it OK if the res_map isn't complete yet.
				continue;
			}

			bool connection_ok( false );

			for ( ResMap::const_iterator
							it=connections_in_big_pose.begin(), it_end = connections_in_big_pose.end(); it != it_end; ++it ) {

				Size const res1_in_big_pose = it->first;
				Size const res2_in_big_pose = it->second;

				//				TR << "Checking connection: " << res1_in_big_pose << " - " << res2_in_big_pose << std::endl;

				ResMap::const_iterator res1_map_id = res_map.find( res1_in_big_pose );
				ResMap::const_iterator res2_map_id = res_map.find( res2_in_big_pose );

				if ( res1_map_id != res_map.end() &&
						 res2_map_id != res_map.end() ) {

					Size const res1_in_scratch_pose = res1_map_id->second;
					Size const res2_in_scratch_pose = res2_map_id->second;
					//					TR << " Checking connection in scratch_pose: " << res1_in_scratch_pose << " - " << res2_in_scratch_pose << std::endl;

					if ( ( chain_id[ res1 ] == chain_id[ res1_in_scratch_pose ] &&
								 chain_id[ res2 ] == chain_id[ res2_in_scratch_pose ] ) ||
							 ( chain_id[ res1 ] == chain_id[ res2_in_scratch_pose ] &&
								 chain_id[ res2 ] == chain_id[ res1_in_scratch_pose ] ) ) {
						connection_ok = true;
						break;
					}

				}
			}

			if (!connection_ok) {
				//				TR << "FAIL! Could not find a match to " << res1 << " -- " << res2 << std::endl;
				return false;
			}

		}

		return true;
	}


	//////////////////////////////////////////////////////////////////////////////
	bool
	RNA_ChunkLibrary::check_res_map( ResMap const & res_map, pose::Pose const & scratch_pose, std::string const & sequence ) const{

		// CHECK SEQUENCE HERE!!! EXIT IF NO MATCH!!!!!!!

		for ( ResMap::const_iterator
						it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

			// For now, just do bonded atoms...update later to do jumps too.
			Size const i = it->first; //Index in big pose.
			Size const i_scratch_pose = it->second; // Index in the little "chunk" or "scratch" pose

			if ( sequence[ i-1 ] != scratch_pose.residue( i_scratch_pose ).name1() ){
				utility_exit_with_message(  "Mismatched sequence!!" );
				return false;
			}

		}
		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::process_input_file( std::string const & input_file,
																				utility::vector1< pose::PoseOP > & pose_list,
																				bool is_pdb /*= false*/ ) const
	{
		using namespace core::io::silent;
		using namespace protocols::rna;

		core::chemical::ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );

		if ( is_pdb ){

			pose::PoseOP pose_op( new pose::Pose );
			core::import_pose::pose_from_pdb( *pose_op, *rsd_set, input_file );
			ensure_phosphate_nomenclature_matches_mini( *pose_op );
			figure_out_reasonable_rna_fold_tree( *pose_op );
			pose_list.push_back( pose_op );

		} else { //its a silent file.

			SilentFileData silent_file_data;
			silent_file_data.read_file( input_file );
			for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
							end = silent_file_data.end(); iter != end; ++iter ) {
				pose::PoseOP pose_op( new pose::Pose );
				iter->fill_pose( *pose_op );
				pose_list.push_back( pose_op );
			}

		}

		// further cleanup.
		for (Size n = 1; n <= pose_list.size(); n++ ){

			pose::PoseOP pose_op = pose_list[ n ];

			remove_cutpoints_closed( *pose_op );

			if ( coarse_rna_ && !pose_op->residue(1).is_coarse() ){
				pose::Pose coarse_pose;
				make_coarse_pose( *pose_op, coarse_pose );
				*pose_op = coarse_pose;
			}

			virtualize_5prime_phosphates( *pose_op );
		}

		// std::cout << "DONE: " << input_file << std::endl;

		if ( pose_list.size() < 1)  {
			utility_exit_with_message(  "No structure found in input file  " + input_file );
		}

	}

	////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::initialize_random_chunks( pose::Pose & pose, bool const dump_pdb /* = false */) const{
		for ( Size n = 1; n <= num_chunk_sets(); n++ ) {

			ChunkSet const & chunk_set( *chunk_sets_[ n ] );

			Size chunk_index = static_cast<int>( RG.uniform() * chunk_set.num_chunks() ) + 1;

			// JUST FOR TESTING
			if ( dump_pdb ) chunk_index = 1;

			//TR << "NUM_CHUNKS " << chunk_index << " " << chunk_set.num_chunks() << std::endl;
			chunk_set.insert_chunk_into_pose( pose, chunk_index, allow_insert_ );


			// useful for tracking homology modeling: perhaps we can align to first chunk as well -- 3D alignment of Rosetta poses are
			// arbitrarily set to origin (except in special cases with virtual residues...)
			if ( n==1  /*&&  pose.residue( pose.total_residue() ).name3() != "VRT"*/ ) align_to_chunk( pose, chunk_set, chunk_index  );

			if ( dump_pdb ) pose.dump_pdb( "start_"+string_of(n)+".pdb" );

		}

		//exit( 0 );

	}

	///////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::superimpose_to_first_chunk( pose::Pose & pose ) const{
		runtime_assert( chunk_sets_.size() > 0 );
		ChunkSet const & chunk_set( *chunk_sets_[ 1 ] );
		align_to_chunk( pose, chunk_set,  1  );
	}

	///////////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::align_to_chunk( pose::Pose & pose, ChunkSet const & chunk_set, Size const chunk_index ) const{

		using namespace core::id;

		std::map< AtomID, AtomID > atom_id_map = chunk_set.get_atom_id_map( pose, allow_insert_ );

		id::AtomID_Map< id::AtomID >  alignment_atom_id_map; // weird alternative format needed for superimpose_pose
		core::pose::initialize_atomid_map( alignment_atom_id_map, pose, id::BOGUS_ATOM_ID );
		for ( std::map< AtomID, AtomID >::const_iterator
						it=atom_id_map.begin(), it_end = atom_id_map.end(); it != it_end; ++it ) {
			alignment_atom_id_map.set( it->first, it->second );
		}

		core::scoring::superimpose_pose( pose, *(chunk_set.mini_pose( chunk_index )), alignment_atom_id_map );
	}



	////////////////////////////////////////////////////////////////
	void
	RNA_ChunkLibrary::set_allow_insert(toolbox::AllowInsertOP allow_insert ){
		allow_insert_ = allow_insert;
	}

}
}

