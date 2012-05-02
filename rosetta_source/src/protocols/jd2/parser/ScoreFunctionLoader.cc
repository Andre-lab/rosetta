// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/parser/DataLoader.cc
/// @brief  Implementation of the XML parser's DataLoader base class (ctor & dstor)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <protocols/jd2/parser/ScoreFunctionLoader.hh>
#include <protocols/jd2/parser/StandardLoaderCreators.hh>

// Project Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <protocols/moves/DataMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#define foreach BOOST_FOREACH

namespace protocols {
namespace jd2 {
namespace parser {

static basic::Tracer TR( "protocols.jd2.parser.ScoreFunctionLoader" );

ScoreFunctionLoader::ScoreFunctionLoader() {}
ScoreFunctionLoader::~ScoreFunctionLoader() {}

void ScoreFunctionLoader::load_data(
	core::pose::Pose const &,
	utility::tag::TagPtr const tag,
	moves::DataMap & data
) const
{
	using namespace utility::tag;
	typedef utility::vector0< TagPtr > TagPtrs;

	TagPtrs const scorefxn_tags( tag->getTags() );

	foreach(TagPtr scorefxn_tag, scorefxn_tags){
		using namespace core::scoring;
		using namespace core::scoring::symmetry;

		ScoreFunctionOP in_scorefxn;
		std::string const scorefxn_name( scorefxn_tag->getName() );
		std::string const scorefxn_weights( scorefxn_tag->getOption<std::string>( "weights", "standard" ) );
		if(  scorefxn_tag->hasOption( "weights" ) && scorefxn_tag->hasOption( "patch" ) ) {
			std::string const scorefxn_patch( scorefxn_tag->getOption<std::string>( "patch" ) );
			in_scorefxn = ScoreFunctionFactory::create_score_function( scorefxn_weights, scorefxn_patch);
			TR << "defined score function \"" << scorefxn_name << "\" with weights \""
				<< scorefxn_weights << "\" and patch \"" << scorefxn_patch << "\"\n";
		} else if ( scorefxn_tag->hasOption( "weights" ) ) {
			in_scorefxn = ScoreFunctionFactory::create_score_function( scorefxn_weights );
			TR << "defined score function \"" << scorefxn_name << "\" with weights \""
				<< scorefxn_weights << "\"\n";
		} else {
			in_scorefxn = new ScoreFunction;
			in_scorefxn->reset();
			TR << "***WARNING***: No weights/patch defined. Defining " << scorefxn_name << " with all-zero weights.\n";
		}
		foreach(TagPtr mod_tag, scorefxn_tag->getTags()){
			if( mod_tag->getName() == "Reweight" ) {
				std::string const scoretype_name( mod_tag->getOption<std::string>( "scoretype" ) );
				core::Real const weight( mod_tag->getOption<core::Real>( "weight" ) );
				TR<<"setting "<<scorefxn_name<<" weight " << scoretype_name << " to " << weight<<'\n';
				core::scoring::ScoreType const type = score_type_from_name( scoretype_name );
				in_scorefxn->set_weight( type, weight );
			}

			// Set energy method options:
			if( mod_tag->getName() == "Set" ){
				core::scoring::methods::EnergyMethodOptions emoptions( in_scorefxn->energy_method_options() );
				core::scoring::hbonds::HBondOptionsOP hboptions( emoptions.hbond_options() );

				if( mod_tag->hasOption( "softrep_etable" )) {
					if ( mod_tag->getOption<bool>( "softrep_etable" )) {
						emoptions.etable_type( core::scoring::FA_STANDARD_SOFT );

					}
				}

				if( mod_tag->hasOption( "hack_elec_min_dis" )) {
					emoptions.hackelec_min_dis( mod_tag->getOption<core::Real>( "hack_elec_min_dis" ) );
				}
				if( mod_tag->hasOption( "hack_elec_max_dis" )) {
					emoptions.hackelec_max_dis( mod_tag->getOption<core::Real>( "hack_elec_max_dis" ) );
				}
				if( mod_tag->hasOption( "hack_elec_dielectric" )) {
					emoptions.hackelec_die( mod_tag->getOption<core::Real>( "hack_elec_dielectric" ) );
				}
				if( mod_tag->hasOption( "hack_elec_no_dis_dep_die" )) {
					emoptions.hackelec_no_dis_dep_die( mod_tag->getOption<bool>( "hack_elec_no_dis_dep_die" ) );
				}
				if( mod_tag->hasOption( "exclude_protein_protein_hack_elec" )) {
					emoptions.exclude_protein_protein_hack_elec( mod_tag->getOption<bool>( "exclude_protein_protein_hack_elec" ) );
				}
				if( mod_tag->hasOption( "exclude_DNA_DNA" )) {
					emoptions.exclude_DNA_DNA( mod_tag->getOption<bool>( "exclude_DNA_DNA" ) );
				}
				if( mod_tag->hasOption( "exclude_DNA_DNA_hbond" )) {
					hboptions->exclude_DNA_DNA( mod_tag->getOption<bool>( "exclude_DNA_DNA_hbond" ) );
				}
				if( mod_tag->hasOption( "use_hb_env_dep_DNA" )) {
					hboptions->use_hb_env_dep_DNA( mod_tag->getOption<bool>( "use_hb_env_dep_DNA" ) );
				}
				if( mod_tag->hasOption( "use_hb_env_dep" )) {
					hboptions->use_hb_env_dep( mod_tag->getOption<bool>( "use_hb_env_dep" ) );
				}
				if( mod_tag->hasOption( "smooth_hb_env_dep" )) {
					hboptions->smooth_hb_env_dep( mod_tag->getOption<bool>( "smooth_hb_env_dep" ) );
				}
				if( mod_tag->hasOption( "decompose_bb_hb_into_pair_energies" )) {
					hboptions->decompose_bb_hb_into_pair_energies( mod_tag->getOption<bool>( "decompose_bb_hb_into_pair_energies" ) );
				}

				in_scorefxn->set_energy_method_options( emoptions );
			}
		} // Mod tags

		// weights for arbitrary ScoreFunctions should be tampered with only as a consequence of user input--NEVER by default

		// hotspot hash constraint
		if ( scorefxn_tag->hasOption("hs_hash") ) {
			core::Real hotspot_hash( 0.0 ); // APL FIX THIS!  This used to be initialized when the HotspotHashingConstraints were read in.
			core::Real const hs_hash( scorefxn_tag->getOption<core::Real>( "hs_hash", hotspot_hash ) );
			TR<<"setting "<<scorefxn_name<<" backbone_stub_constraint to "<<hs_hash<<'\n';
			in_scorefxn->set_weight( backbone_stub_constraint, hs_hash );
		}

		//fpd should we symmetrize scorefunction?
		bool const scorefxn_symm( scorefxn_tag->getOption<bool>( "symmetric", 0 ) );
		if (scorefxn_symm) {
			in_scorefxn = ScoreFunctionOP( new SymmetricScoreFunction( in_scorefxn ) );
			TR<<"symmetrizing "<<scorefxn_name<<'\n';
		}

		data.add( "scorefxns" , scorefxn_name, in_scorefxn );
	}//end user-defined scorefxns
	TR.flush();
}

DataLoaderOP
ScoreFunctionLoaderCreator::create_loader() const { return new ScoreFunctionLoader; }

std::string
ScoreFunctionLoaderCreator::keyname() const { return "SCOREFXNS"; }


} //namespace parser
} //namespace jd2
} //namespace protocols
