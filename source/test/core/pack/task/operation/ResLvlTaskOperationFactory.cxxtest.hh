// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/OperateOnResidueSubset.cxxtest.hh
/// @brief  test suite for core::select::residue_selector::OperateOnResidueSubset
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <core/pack/task/operation/ResLvlTaskOperationFactory.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>
#include <core/pack/task/operation/ResLvlTaskOperationCreator.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/XMLSchemaValidation.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

#include <core/init_util.hh> // AUTO IWYU For core_init
#include <iostream> // AUTO IWYU For endl, basic_ostream, cout, ostream

using namespace core::pack::task::operation;
using namespace utility::tag;

class DummyRLTOCreator : public ResLvlTaskOperationCreator {
public:
	ResLvlTaskOperationOP create_res_level_task_operation() const override { return ResLvlTaskOperationOP( 0 ); }
	std::string keyname() const override { return "DummyResLvlTaskOperation"; }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override
	{
		res_lvl_task_op_schema_empty( xsd, keyname() );
	}
	bool skip_citation_info_in_schema() const override { return true; }

};

class DummyRLTO2Creator : public ResLvlTaskOperationCreator {
public:
	ResLvlTaskOperationOP create_res_level_task_operation() const override { return ResLvlTaskOperationOP( 0 ); }
	std::string keyname() const override { return "DummyResLvlTaskOperation2"; }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override
	{
		XMLSchemaComplexType ct;

		// this is the error: the ComplexType name needs to come from
		// the complex_type_name_for_res_lvl_task_op function
		ct.name( keyname() );

		xsd.add_top_level_element( ct );
	}
	bool skip_citation_info_in_schema() const override { return true; }

};

class ResLvlTaskOperationFactoryTests : public CxxTest::TestSuite {
private:
	bool dummy_initialized_;
public:

	// ctor is called exactly once per unit test suite execution (and once per execution of core.test, for that matter)
	// register the dummy creator here to avoid it happening multiple times.
	ResLvlTaskOperationFactoryTests() : dummy_initialized_( false ) {}

	void setUp() {
		if ( !dummy_initialized_ ) {
			dummy_initialized_ = true;
			ResLvlTaskOperationFactory::get_instance()->factory_register( utility::pointer::make_shared< DummyRLTOCreator >());
		}
		core_init();
	}

	void test_res_lvl_task_op_factory_xsd_creation() {
		XMLSchemaDefinition xsd;
		ResLvlTaskOperationFactory::get_instance()->define_res_lvl_task_op_xml_schema( xsd );
		TS_ASSERT( xsd.has_top_level_element( complex_type_name_for_res_lvl_task_op( "DummyResLvlTaskOperation" )));
		//std::cout << "Res lvl task operation factory XSD:\n" << xsd.full_definition() << std::endl;
	}

	void test_all_res_lvl_task_ops_define_valid_xsds() {
		XMLSchemaDefinition xsd;
		ResLvlTaskOperationFactory::get_instance()->define_res_lvl_task_op_xml_schema( xsd );

		XMLSchemaModelGroupOP rosetta_scripts_seq( new XMLSchemaModelGroup( xsmgt_sequence ));

		XMLSchemaComplexTypeOP rosetta_scripts_type( new XMLSchemaComplexType );
		rosetta_scripts_type->set_model_group( rosetta_scripts_seq );

		XMLSchemaElementOP rosetta_scripts_element( new XMLSchemaElement );
		rosetta_scripts_element->name( "ROSETTASCRIPTS" );
		rosetta_scripts_element->element_type_def( rosetta_scripts_type );

		xsd.add_top_level_element( *rosetta_scripts_element );

		//std::cout << "xsd for res lvl task operations:\n" << xsd.full_definition() << std::endl;
		XMLValidationOutput output = test_if_schema_is_valid( xsd.full_definition() );
		TS_ASSERT( output.valid() );
		if ( ! output.valid() ) {
			std::cout << output.error_messages() << std::endl;
		}
	}

	void test_res_lvl_task_op_factory_w_bad_xsd() {
		ResLvlTaskOperationFactory::get_instance()->factory_register( utility::pointer::make_shared< DummyRLTO2Creator >());
		try {
			XMLSchemaDefinition xsd;
			ResLvlTaskOperationFactory::get_instance()->define_res_lvl_task_op_xml_schema( xsd );
			TS_ASSERT( false );
		} catch (utility::excn::Exception const & e ) {
			std::string expected_message =
				"Could not generate an XML Schema for ResLvlTaskOperations from ResLvlTaskOperationFactory; offending class"
				" must call core::pack::task::operation::complex_type_name_for_res_lvl_task_op when defining its XML Schema"; // "\ndefine_xml_schema_group: failed to detect a complex type of name \"" + complex_type_name_for_res_lvl_task_op( "DummyResLvlTaskOperation2" ) + "\" for \"DummyResLvlTaskOperation2\"\n";

			//std::cout << e.msg() << "\n" << expected_message << std::endl;
			TS_ASSERT( e.msg().find(expected_message) != std::string::npos );
		}
	}

};
