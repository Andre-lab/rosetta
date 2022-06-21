// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Zhe Zhang

#include <devel/replica_docking/CaIrmsdFilter.hh>
#include <devel/replica_docking/CaIrmsdFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/docking/metrics.hh>


#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

#include <core/scoring/ScoreFunction.hh> // AUTO IWYU For ScoreFunction
#include <basic/Tracer.hh> // AUTO IWYU For Tracer, Tracer::TracerProxy
#include <basic/options/option.hh> // AUTO IWYU For OptionCollection, option

// Project Headers


static basic::Tracer TR( "devel.replica_docking.CaIrmsdFilter" );

namespace devel {
namespace replica_docking {



void CaIrmsdFilter::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( in::file::native );
}


CaIrmsdFilter::CaIrmsdFilter() :
	Filter( "Ca_Irms" ),
	lower_threshold_( 0.0 ),
	upper_threshold_(9999)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	scorefxn_ = core::scoring::get_score_function();
	// scorefxn_->show(TR.Info);
	movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
	TR << "End constructer"<<std::endl;

}

CaIrmsdFilter::CaIrmsdFilter( core::scoring::ScoreFunctionOP sfxn, core::Size const rb_jump,core::Real const lower_threshold, core::Real const upper_threshold ) :
	Filter( "Ca_Irms" ),
	lower_threshold_( lower_threshold ),
	upper_threshold_(upper_threshold)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose() );
		core::import_pose::pose_from_file( *native_pose, option[ in::file::native ], core::import_pose::PDB_file);
		native_pose_ = native_pose;
	} else {
		utility_exit_with_message("need to specify native pdb to calculate Ca-Irms");
	}

	if ( !sfxn ) {
		scorefxn_ = core::scoring::get_score_function();
	} else {
		scorefxn_ = sfxn->clone();
	}
	TR.Info <<"CaIrmsdEvaluator: "<<"score" << std::endl;
	// scorefxn_->show(TR.Info);
	movable_jumps_.push_back( rb_jump );
	TR << "End constructer"<<std::endl;

}

CaIrmsdFilter::~CaIrmsdFilter() = default;

protocols::filters::FilterOP
CaIrmsdFilter::clone() const{
	return utility::pointer::make_shared< CaIrmsdFilter >( *this );
}

protocols::filters::FilterOP
CaIrmsdFilter::fresh_instance() const{
	return utility::pointer::make_shared< CaIrmsdFilter >();
}

void
CaIrmsdFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
) {

	std::string const scorefxn_name(
		protocols::rosetta_scripts::get_score_function_name(tag) );
	scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( scorefxn_name );
	// //  scorefxn_ = new core::scoring::ScoreFunction( *(data.get< core::scoring::ScoreFunction * >( "scorefxns", scorefxn_name )) );

	// // scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );

	lower_threshold_ = tag->getOption<core::Real>( "threshold", 0.0 );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", 9999.0);
	jump( tag->getOption< core::Size >( "jump", 1 ));

}

bool
CaIrmsdFilter::apply( core::pose::Pose const & pose ) const {

	if ( native_pose_ == nullptr ) {
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			core::pose::PoseOP native_pose( utility::pointer::make_shared< core::pose::Pose >() );
			core::import_pose::pose_from_file(
				*native_pose,
				basic::options::option[ basic::options::OptionKeys::in::file::native ](),
				core::import_pose::PDB_file
			);
			native_pose_ = native_pose;
		} else {
			utility_exit_with_message("need to specify native pdb with the -in:file:native commandline option to calculate Irms.");
		}
	}

	core::Real const Irms( compute( pose ) );

	TR<<"Ca_Irms is "<<Irms<<". ";
	if ( Irms >= lower_threshold_ && Irms <= upper_threshold_ ) {
		TR<<"passing." <<std::endl;
		return true;
	} else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
CaIrmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const Irms( compute( pose ));
	out<<"Ca_Irms= "<< Irms<<'\n';
}

core::Real
CaIrmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const Irms( compute( pose ));
	return( Irms );
}

void
CaIrmsdFilter::jump( core::Size const jump_id )
{
	movable_jumps_.push_back( jump_id );
}


core::Real
CaIrmsdFilter::compute( core::pose::Pose const & pose ) const {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( native_pose_ == nullptr ) {
		if ( option[ in::file::native ].user() ) {
			core::pose::PoseOP native_pose( utility::pointer::make_shared< core::pose::Pose >() );
			core::import_pose::pose_from_file( *native_pose, option[ in::file::native ](), core::import_pose::PDB_file);
			native_pose_ = native_pose;
		} else {
			utility_exit_with_message("Need to specify native pdb with -in:file:native to calculate Irms.");
		}
	}

	TR<<"compute Irms"<< std::endl;
	core::Real irms = protocols::docking::calc_CA_Irmsd( pose, *native_pose_, scorefxn_, movable_jumps_ );
	return( irms );
}

std::string CaIrmsdFilter::name() const {
	return class_name();
}

std::string CaIrmsdFilter::class_name() {
	return "Ca_Irms";
}

void CaIrmsdFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::rosetta_scripts::attributes_for_get_score_function_name( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Threshold below which the interaction score filter fails", "-30" )
		+ XMLSchemaAttribute::attribute_w_default( "upper_threshold", xsct_real, "Threshold above which the interaction score filter fails", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "jump", xsct_positive_integer, "Jump across which the interface is defined, numbered sequentially from 1", "1" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string CaIrmsdFilterCreator::keyname() const {
	return CaIrmsdFilter::class_name();
}

protocols::filters::FilterOP
CaIrmsdFilterCreator::create_filter() const {
	return utility::pointer::make_shared< CaIrmsdFilter >();
}

void CaIrmsdFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CaIrmsdFilter::provide_xml_schema( xsd );
}


}
}
