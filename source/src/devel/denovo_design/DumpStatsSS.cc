// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/DumpStatsSS.cc
/// @brief
/// @author TJ Brunette
//
// Unit headers
#include <devel/denovo_design/DumpStatsSS.hh>
#include <devel/denovo_design/DumpStatsSSCreator.hh>
#include <core/io/external/PsiPredInterface.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/parser/BluePrint.fwd.hh>
#include <protocols/parser/BluePrint.hh>
#include <protocols/ss_prediction/SS_predictor.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <utility/sys_util.hh>


#include <utility/io/ozstream.hh>
#include <ObjexxFCL/format.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace devel {
namespace denovo_design {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "devel/denovo_design/movers/DumpStatsSS" );

// XRW TEMP std::string
// XRW TEMP DumpStatsSSCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return DumpStatsSS::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP DumpStatsSSCreator::create_mover() const {
// XRW TEMP  return utility::pointer::make_shared< DumpStatsSS >();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP DumpStatsSS::mover_name()
// XRW TEMP {
// XRW TEMP  return "DumpStatsSS";
// XRW TEMP }

DumpStatsSS::DumpStatsSS():
	protocols::moves::Mover( DumpStatsSS::mover_name() ),
	scorefxn_(/* 0 */)
{}

DumpStatsSS::DumpStatsSS(DumpStatsSS const &rval):
	protocols::moves::Mover( DumpStatsSS::mover_name() ),
	fname_(rval.fname_),
	output_(rval.output_),
	scorefxn_(rval.scorefxn_->clone()),
	psipred_cmd_(rval.psipred_cmd_),
	psipred_interface_(rval.psipred_interface_),
	ss_predictor_(rval.ss_predictor_),
	blueprint_(rval.blueprint_),
	start_time_(rval.start_time_)
{}


DumpStatsSS::~DumpStatsSS() = default;

void DumpStatsSS::apply( core::pose::Pose & pose ) {
	using namespace ObjexxFCL::format;
	//Keep file open. & dump to it.
	core::Real current_time = clock();
	core::Real time = (current_time-start_time_)/core::Real(CLOCKS_PER_SEC);
	core::Real fa_score = 0;
	if ( scorefxn_ ) {
		fa_score = scorefxn_->score(pose);
	}
	std::string wanted_ss;
	if ( blueprint_ ) {
		wanted_ss = blueprint_->secstruct();
	} else {
		core::scoring::dssp::Dssp dssp( pose );
		wanted_ss = dssp.get_dssp_secstruct();
	}
	std::string sequence;
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		if ( pose.residue( i ).is_protein() ) sequence += pose.residue( i ).name1();
	}
	core::Real svm_prob = compute_svm_prob(sequence,wanted_ss);
	core::Real psipred_prob = compute_psipred_prob(pose,wanted_ss);
	*output_ << F(8,3,fa_score) << "  " << F(8,3,svm_prob) << "  " << F(8,3,psipred_prob) << "  " << time << "  " << sequence << std::endl;
}

core::Real DumpStatsSS::compute_svm_prob(std::string sequence, std::string wanted_ss){
	runtime_assert( ss_predictor_ != nullptr );
	runtime_assert( sequence.size() == wanted_ss.size() );
	utility::vector1< utility::vector1< core::Real > > ss_pred( ss_predictor_->predict_ss( sequence ) );
	utility::vector1< core::Real > probabilities;
	for ( core::Size i=1; i<=wanted_ss.size(); ++i ) {
		probabilities.push_back( protocols::ss_prediction::get_prob( wanted_ss[i-1], ss_pred[i] ) );
	}
	return compute_boltz_sum( probabilities );
}

core::Real
DumpStatsSS::compute_psipred_prob(core::pose::Pose & pose, std::string wanted_ss)
{
	runtime_assert( psipred_interface_ != nullptr );
	core::io::external::PsiPredResult const & psipred_result =
		psipred_interface_->run_psipred( pose, wanted_ss );
	return compute_boltz_sum( generate_prob( psipred_result, wanted_ss ) );
}

/// @brief computes the weighted boltzmann sum of the passed vector
core::Real
DumpStatsSS::compute_boltz_sum( utility::vector1< core::Real > const & probabilities ) const
{
	core::Real temp(0.6); //hard coded from SSPredictionFilter.cc
	core::Real sum( 0.0 );
	for ( core::Size i=1; i<=probabilities.size(); ++i ) {
		sum += exp( -probabilities[i] / temp );
	}
	return sum/probabilities.size();
}

void DumpStatsSS::set_scorefxn( core::scoring::ScoreFunctionOP scorefxn) {
	scorefxn_ = scorefxn;
}

void
DumpStatsSS::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	fname_ = tag->getOption<std::string>( "fname", "SS_stats.txt" );
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}
	psipred_cmd_ = tag->getOption<std:: string > ("cmd", "");
	psipred_interface_ = utility::pointer::make_shared< core::io::external::PsiPredInterface >( psipred_cmd_ );
	ss_predictor_ = utility::pointer::make_shared< protocols::ss_prediction::SS_predictor >( "HLE" );
	std::string blueprint_file = tag->getOption< std::string >( "blueprint", "" );
	if ( blueprint_file != "" ) {
		TR << "Dssp-derived secondary structure will be overridden by user specified blueprint file." << std::endl;
		blueprint_ = utility::pointer::make_shared< protocols::parser::BluePrint >( blueprint_file );
		if ( ! blueprint_ ) {
			utility_exit_with_message("There was an error getting the blueprint file   loaded.");
		}
	}

	start_time_ = clock();
	output_ = new utility::io::ozstream(fname_);
	*output_ << "score svmProb psiProb time sequence" << std::endl;
}

// XRW TEMP std::string
// XRW TEMP DumpStatsSS::get_name() const {
// XRW TEMP  return DumpStatsSS::mover_name();
// XRW TEMP }

std::string DumpStatsSS::get_name() const {
	return mover_name();
}

std::string DumpStatsSS::mover_name() {
	return "DumpStatsSS";
}

void DumpStatsSS::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "fname", xs_string, "Secondary structure file to which to dump statistics", "SS_stats.txt" );
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute( "cmd", xs_string, "Command to run psipred" )
		+ XMLSchemaAttribute( "blueprint", xs_string, "Blueprint file that could be used to override DSSP secondary structure." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string DumpStatsSSCreator::keyname() const {
	return DumpStatsSS::mover_name();
}

protocols::moves::MoverOP
DumpStatsSSCreator::create_mover() const {
	return utility::pointer::make_shared< DumpStatsSS >();
}

void DumpStatsSSCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DumpStatsSS::provide_xml_schema( xsd );
}

}//devel
}//denovo_design

