// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/apps/pilot/vmullig/bundle_gui/mainwindow.cpp
/// @brief  Main window class for the bundle GUI project.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include "mainwindow.h"
#include "ui_mainwindow.h"

//Qt includes:
#include <QScrollArea>
#include <QDebug>
#include <QVBoxLayout>
#include <QFileDialog>
#include <QDebug>

//Rosetta numeric headers
#include <numeric/conversions.hh>

//Rosetta core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/LayerSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/BooleanValuedParameter.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/import_pose/import_pose.hh>

//Rosetta protocols headers
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/helical_bundle/PerturbBundle.hh>
#include <protocols/helical_bundle/MakeBundleHelix.hh>
#include <protocols/helical_bundle/parameters/OmegaBundleParameter.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>

//Rosetta utility headers
#include <utility/pointer/memory.hh>

//Rosetta ui headers
#include <ui/ui_protocols/helical_bundle/HelicalBundleDialogueWidget.h>
#include <ui/ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.h>
#include <ui/ui_protocols/helical_bundle/HelixOptionWidget.h>

//C++ headers
#include <sstream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
	ui(new Ui::MainWindow),
	pose_(new core::pose::Pose),
	nonparametric_pose_(nullptr), //Null unless a pose is loaded.
	nonparametric_pose_offset_(numeric::xyzVector<core::Real>(0,0,0)),
	nonparametric_pose_rotation_( numeric::xyzTransform< core::Real >::identity() ),
	sfxn_(nullptr),
	layer_selector_( new core::select::residue_selector::LayerSelector ),
	residue_selector_( layer_selector_ ),
	need_to_update_pose_(false),
	need_to_rebuild_pose_(false),
	timer_(new QTimer),
	pose_draw_widget_( new ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget(this) ),
	perturb_bundle_( utility::pointer::make_shared< protocols::helical_bundle::PerturbBundle >() )
{

    ui->setupUi(this);
	connect( timer_, SIGNAL(timeout()), this, SLOT( on_timer_refresh() ) );

	ui->tab_1->deleteLater();

	update_layer_selector_settings();

	ui->radiobutton_colourby_selection->setAutoExclusive(true);
	ui->radiobutton_colourby_totalscore->setAutoExclusive(true);

	ui->verticalLayout_opengl->addWidget(pose_draw_widget_);
	pose_draw_widget_->set_pose(pose_);
	pose_draw_widget_->set_color_mode( ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget::ColorMode::selection );
	pose_draw_widget_->set_residue_selector(residue_selector_);
	pose_draw_widget_->show();

	ui->presets_comboBox->addItem(QString("Select preset"));
    ui->presets_comboBox->addItem(QString("3-helix bundle"));
    ui->presets_comboBox->addItem(QString("Antiparallel beta-barrel"));

	connect( pose_draw_widget_, SIGNAL(nonparametric_geometry_has_moved()), this, SLOT(on_nonparametric_geometry_moved()) );

	timer_->start(50);

}

MainWindow::~MainWindow()
{
    delete ui;
}

/////////////PRIVATE FUNCTIONS/////////////////////////////


/// @brief Use the BundleGridSampler to update the pose.
void
MainWindow::update_pose() {
	using namespace protocols::helical_bundle;

	PerturbBundle & pertbundle( *perturb_bundle_ );
	BundleParametrizationCalculatorOP pertbundle_calculator( pertbundle.default_calculator() );

	bool const helix_count_changed( pertbundle.n_helices() != static_cast<core::Size>(ui->tabWidget->count()) );

	if( helix_count_changed ) {
		pertbundle_calculator->set_use_degrees(true);
		pertbundle.reset_helices();
	}

	core::Size nhelices(0);

	for(core::Size i(0), imax(ui->tabWidget->count()); i<imax; ++i) {
		QScrollArea* curscrollarea(  dynamic_cast< QScrollArea* >( ui->tabWidget->widget(i) ) );
		if(curscrollarea == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget_ptr(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
		if(curwidget_ptr == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget& curwidget( *curwidget_ptr );

		++nhelices;
		debug_assert(curwidget.helix_index() == nhelices); //Should be true

		protocols::helical_bundle::BundleParametrizationCalculator & curhelix_calculator( *(helix_count_changed ? pertbundle.add_helix(nhelices) : pertbundle.individual_helix_calculator( i+1 )) );
		for( core::Size iparam(1); iparam<=static_cast<core::Size>( BPC_last_parameter_to_be_sampled ); ++iparam ) {
			core::conformation::parametric::RealValuedParameterOP curparam_op( curhelix_calculator.real_parameter( iparam ) );
			if( curparam_op == nullptr ) continue;
			core::conformation::parametric::RealValuedParameter & curparam( *curparam_op );
			if( iparam == static_cast< core::Size >( BPC_omega0 ) && curwidget.omega0_pitchcopy_helix() != 0  ) {
#ifdef NDEBUG
				protocols::helical_bundle::parameters::OmegaBundleParameterOP curparam_omega_op( utility::pointer::static_pointer_cast< protocols::helical_bundle::parameters::OmegaBundleParameter >( curparam_op ) );
#else
				protocols::helical_bundle::parameters::OmegaBundleParameterOP curparam_omega_op( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::OmegaBundleParameter >( curparam_op ) );
				debug_assert( curparam_omega_op != nullptr );
#endif
				protocols::helical_bundle::parameters::OmegaBundleParameter & curparam_omega( *curparam_omega_op );
				curparam_omega.set_copies_pitch(true);
				curparam_omega.set_copy_from_parameters_index( curwidget.omega0_pitchcopy_helix() );
				continue;
			}

			ui::ui_protocols::helical_bundle::HelixOptionWidget & widget( *(curwidget.get_subwidget_nonconst(iparam) ) );
			if( widget.widget_mode() == ui::ui_protocols::helical_bundle::HOW_set_value ) {
				curparam.reset_copying_settings();
				curparam.set_value( widget.value() );
			} else {
				curparam.set_copy_from_parameters_index( widget.helix_to_copy() );
			}
		}
	}
	if(nhelices > 0) {
		pertbundle.apply(*pose_);
		//qDebug() << "Updating a helix"; //DELETE ME
	}

	if( sfxn_ != nullptr ) (*sfxn_)(*pose_);

	pose_draw_widget_->update_pose_draw();
}

void
MainWindow::rebuild_pose_from_scratch() {
	pose_->clear();

	if( nonparametric_pose_ != nullptr ) {
		(*pose_) = (*nonparametric_pose_);
		on_nonparametric_geometry_moved(); //Translate the geometry appropriately.
	}

	protocols::helical_bundle::MakeBundle mkbundle;
	mkbundle.set_use_degrees(true);
	mkbundle.set_reset_pose(pose_->total_residue() == 0);

	core::Size nhelices(0);

	for(core::Size i(0), imax(ui->tabWidget->count()); i<=imax; ++i) {
		QScrollArea* curscrollarea(  dynamic_cast< QScrollArea* >( ui->tabWidget->widget(i) ) );
		if(curscrollarea == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget_ptr(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
		if(curwidget_ptr == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget &curwidget( *curwidget_ptr );

		++nhelices;
		debug_assert(curwidget.helix_index() == nhelices); //Should be true.
		mkbundle.add_helix();

		protocols::helical_bundle::MakeBundleHelix &helix( *(mkbundle.helix(nhelices)) );

		if( curwidget.custom_params_file_selected() ) {
			if( curwidget.custom_params_file() == "A_form_RNA_helix" ) {
				utility::vector1< std::string> URAvec(1, "URA");
				helix.set_residue_name(URAvec);
			} else if( curwidget.custom_params_file() == "14_helix" ) {
				utility::vector1< std::string> B3Avec(1, "B3A");
				helix.set_residue_name(B3Avec);
			} else if ( curwidget.custom_params_file() == "L_alpha_helix" ) {
				utility::vector1< std::string> DALAvec(1, "DALA");
				helix.set_residue_name(DALAvec);
			} else if ( curwidget.custom_params_file() == "a3b" ) {
				utility::vector1< std::string> a3bvec(3, "ALA");
				a3bvec.emplace_back("B3A");
				helix.set_residue_name(a3bvec);
			} else {
				utility::vector1<std::string> ALAvec(1, "ALA");
				helix.set_residue_name(ALAvec);
			}
		} else {
			utility::vector1<std::string> ALAvec(1, "ALA");
			helix.set_residue_name(ALAvec);
		}

		helix.set_helix_length(curwidget.helix_length());
		helix.set_minor_helix_params_from_file(curwidget.params_file()); //NEED TO SET UP CACHING FOR THIS!  CURRENTLY READS FROM DISK WHENEVER THE POSE IS REBUILT FROM SCRATCH!

		protocols::helical_bundle::BundleParametrizationCalculatorOP calculator_op( helix.calculator_op() );
		debug_assert( calculator_op != nullptr );
		protocols::helical_bundle::BundleParametrizationCalculator calculator( *calculator_op );

		for( core::Size iparam(1); iparam<=static_cast< core::Size >(protocols::helical_bundle::BPC_last_parameter_to_be_sampled); ++iparam ) {
			core::conformation::parametric::RealValuedParameterOP curparam_op( calculator.real_parameter(iparam) );
			debug_assert(curparam_op != nullptr);
			core::conformation::parametric::RealValuedParameter curparam( *curparam_op );

			if( iparam == static_cast<core::Size>( protocols::helical_bundle::BPC_omega0 ) && curwidget.omega0_pitchcopy_helix() != 0 ) {
#ifdef NDEBUG
				protocols::helical_bundle::parameters::OmegaBundleParameterOP curparam_omega_op( utility::pointer::static_pointer_cast< protocols::helical_bundle::parameters::OmegaBundleParameter >( curparam_op ) );
#else
				protocols::helical_bundle::parameters::OmegaBundleParameterOP curparam_omega_op( utility::pointer::dynamic_pointer_cast< protocols::helical_bundle::parameters::OmegaBundleParameter >( curparam_op ) );
				debug_assert( curparam_omega_op != nullptr );
#endif
				protocols::helical_bundle::parameters::OmegaBundleParameter & curparam_omega( *curparam_omega_op );
				curparam_omega.set_copies_pitch(true);
				curparam_omega.set_copy_from_parameters_index( curwidget.omega0_pitchcopy_helix() );
				continue;
			}

			ui::ui_protocols::helical_bundle::HelixOptionWidget & widget( *(curwidget.get_subwidget_nonconst(iparam) ) );
			if( widget.widget_mode() == ui::ui_protocols::helical_bundle::HOW_set_value ) {
				curparam.set_value( widget.value() );
			} else {
				curparam.set_copy_from_parameters_index( widget.helix_to_copy() );
			}
		}
		debug_assert(helix.calculator_op()->boolean_parameter( protocols::helical_bundle::BPC_invert_helix ) != nullptr );
		helix.calculator_op()->boolean_parameter( protocols::helical_bundle::BPC_invert_helix )->set_value( curwidget.invert_helix() );
	}
	if(nhelices > 0) {
			mkbundle.apply(*pose_);
			update_pose(); //Needed for omega0 pitch copying.  Calls pose_draw_widget_->update_pose_draw().
	} else {
		pose_draw_widget_->update_pose_draw();
	}
}

/// @brief Update the indices of helices in UI elements.
void
MainWindow::update_helix_indices() {
	core::Size nhelices(0);
	for(core::Size i(0), imax(ui->tabWidget->count()); i<imax; ++i) {
		QScrollArea* curscrollarea(  dynamic_cast< QScrollArea* >( ui->tabWidget->widget(i) ) );
		if(curscrollarea == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
		if(curwidget == nullptr) continue;
		++nhelices;
		curwidget->set_helix_index(nhelices);

		std::stringstream curname;
		curname << "Helix " << nhelices;
		ui->tabWidget->setTabText( i, QString( curname.str().c_str() ) );
	}
}

/// @brief Given a helix to remove, update the dependencies of helices that depend on that helix.
void
MainWindow::update_later_helix_dependencies(
	core::Size const helix_to_remove
) {
	core::Size const helix_count( ui->tabWidget->count() );
	if(helix_to_remove == helix_count) return;
	for(core::Size i(helix_to_remove+1); i<=helix_count; ++i) {
		QScrollArea* curscrollarea(  dynamic_cast< QScrollArea* >( ui->tabWidget->widget(i) ) );
		if(curscrollarea == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
		if(curwidget == nullptr) continue;
		curwidget->reset_controls_dependent_on_other_helix( helix_to_remove );
		curwidget->update_control_dependencies_given_helix_to_be_removed( helix_to_remove );
	}
}

/// @brief Update the layer selector's settings from the UI.
void
MainWindow::update_layer_selector_settings() {
	layer_selector_->set_layers( ui->checkbox_select_core->isChecked() , ui->checkbox_select_boundary->isChecked(), ui->checkbox_select_surface->isChecked() );
}

/// @brief Switch the residue selector to a composite selector that does not select the imported geometry,
/// or just use the layer selector, depending on whether we have imported geometry.
void
MainWindow::update_residue_selector() {
	if( nonparametric_pose_ == nullptr ) {
		residue_selector_ = layer_selector_;
	} else {
		core::select::residue_selector::AndResidueSelectorOP andsel( new core::select::residue_selector::AndResidueSelector );
		core::select::residue_selector::NotResidueSelectorOP notsel( new core::select::residue_selector::NotResidueSelector );
		core::select::residue_selector::ResidueIndexSelectorOP indexsel( new core::select::residue_selector::ResidueIndexSelector );
		for(core::Size i(1), imax( nonparametric_pose_->total_residue() ); i<=imax; ++i) indexsel->append_index(i);
		notsel->set_residue_selector(indexsel);
		andsel->add_residue_selector(notsel);
		andsel->add_residue_selector(layer_selector_);
		residue_selector_ = andsel;
	}

	pose_draw_widget_->set_residue_selector(residue_selector_);
}

/// @brief Given a centre-of-mass vector, calculate the centre of mass of the nonparametric
/// geometry in the pose.
/// @details Does nothing (returns 0, 0, 0) if no nonparametric geometry is loaded.
void
MainWindow::compute_com_vect(
	numeric::xyzVector< core::Real > & comvect
) const {
	comvect = numeric::xyzVector< core::Real >(0, 0, 0);
	if( nonparametric_pose_ == nullptr ) return;

	core::pose::Pose const & pose( *pose_ ); //Get a reference to the pose for speed

	core::Size counter(0);
	for( core::Size ir(1), irmax(nonparametric_pose_->total_residue()); ir<=irmax; ++ir ) {
		for( core::Size ia(1), iamax( pose.residue_type(ir).natoms() ); ia<=iamax; ++ia ) {
			if( pose.residue_type(ir).atom_is_hydrogen(ia) ) continue;
			comvect += pose.xyz( core::id::AtomID(ia,ir) );
			++counter;
		}
	}
	if(counter > 0) {
		comvect /= static_cast<core::Real>(counter);
	}
}



/////////////PRIVATE Q_SLOTS///////////////////////////////

/// @brief What do do when the add_helix button is clicked.
void MainWindow::on_add_helix_button_clicked()
{
	ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* new_widget( new ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget );
	QScrollArea* new_scrollarea( new QScrollArea);
	new_scrollarea->setLayout( new QVBoxLayout );
	new_scrollarea->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	new_scrollarea->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
	new_scrollarea->setWidget(new_widget);
	new_scrollarea->layout()->addWidget( new_widget );
	new_scrollarea->setWidgetResizable(true);

	new_widget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
	ui->tabWidget->addTab( new_scrollarea, QString("Helix") );
	//ui->verticalLayout_3->addWidget(new_widget);
	new_scrollarea->show();
	new_widget->show();
	connect(new_widget, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()));
	connect(new_widget, SIGNAL(control_has_changed_requiring_pose_rebuild()), this, SLOT(on_control_changed_requiring_pose_rebuild()));
	connect( new_widget, SIGNAL( helix_has_been_removed(core::Size const) ), this, SLOT( on_helix_removed(core::Size const) ) );
	ui->tabWidget->update();
	update_helix_indices();
	rebuild_pose_from_scratch();
}

/// @brief What do do when a remove helix button is clicked.
void MainWindow::on_helix_removed(core::Size const helix_to_remove) {
	update_later_helix_dependencies( helix_to_remove );
	ui->tabWidget->removeTab( static_cast<int>( helix_to_remove - 1 ) );
	ui->tabWidget->update();
	update_helix_indices();
	rebuild_pose_from_scratch();
}

/// @brief What do to when a control reports that the user has changed a value.  (Update the pose).
void
MainWindow::on_control_changed() {
	need_to_update_pose_ = true;
}

/// @brief What do do when a control reports that the user has changed a value, AND a full rebuild of the pose is required.
void
MainWindow::on_control_changed_requiring_pose_rebuild() {
	need_to_rebuild_pose_ = true;
}

/// @brief Periodically, check whether we need to update the pose and the display, and do so if we need to.
void MainWindow::on_timer_refresh() {
	if( need_to_rebuild_pose_ ) {
		need_to_rebuild_pose_ = false;
		need_to_update_pose_ = false;
		rebuild_pose_from_scratch();
	} else if(need_to_update_pose_) {
		need_to_update_pose_ = false;
		need_to_rebuild_pose_ = false;
		update_pose();
		//pose_->dump_pdb("temp.pdb"); //DELETE ME -- for testing only.
	}
}

/// @brief The "export PDB" button is clicked.  Dump a PDB file.
void MainWindow::on_pushButton_clicked()
{
	QString const filename( QFileDialog::getSaveFileName(this, "Save PDB file", "", "PDB file (*.pdb)") );
	if( filename.isEmpty() ) return;

	std::string pdb_filename( filename.toStdString() );

	/*if( !pdb_filename.substr(pdb_filename.length()-5, 4).compare(".pdb") ) {
		pdb_filename += ".pdb";
	}*/
	pose_->dump_pdb( pdb_filename );
}

/// @brief The user has clicked the "select core" checkbox.  Update the layer selector.
void MainWindow::on_checkbox_select_core_clicked()
{
	update_layer_selector_settings();
	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has clicked the "select boundary" checkbox.  Update the layer selector.
void MainWindow::on_checkbox_select_boundary_clicked()
{
	update_layer_selector_settings();
	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has clicked the "select surface" checkbox.  Update the layer selector.
void MainWindow::on_checkbox_select_surface_clicked()
{
	update_layer_selector_settings();
	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has clicked the run packer button.
void MainWindow::on_button_run_packer_clicked()
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace core::chemical;

	core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
	core::pack::task::TaskFactoryOP task( new core::pack::task::TaskFactory );

	if( !ui->checkbox_allow_design->isChecked() ) {
		RestrictToRepackingOP repack_only( new RestrictToRepacking );
		task->push_back(repack_only);
	} else {
		utility::vector1< bool > aa_list(20, true);
		aa_list[aa_ala] = ui->checkbox_ala->isChecked();
		aa_list[aa_cys] = ui->checkbox_cys->isChecked();
		aa_list[aa_gly] = ui->checkbox_gly->isChecked();
		aa_list[aa_pro] = ui->checkbox_pro->isChecked();
		aa_list[aa_his] = ui->checkbox_his->isChecked();
		aa_list[aa_met] = ui->checkbox_met->isChecked();

		aa_list[aa_val] = aa_list[aa_ile] = aa_list[aa_leu] = ui->checkbox_aliphatic->isChecked();
		aa_list[aa_trp] = aa_list[aa_tyr] = aa_list[aa_phe] = ui->checkbox_aromatic->isChecked();
		aa_list[aa_lys] = aa_list[aa_arg] = aa_list[aa_asp] = aa_list[aa_glu] = ui->checkbox_charged->isChecked();
		aa_list[aa_gln] = aa_list[aa_asn] = aa_list[aa_thr] = aa_list[aa_ser] = ui->checkbox_polar->isChecked();

		RestrictAbsentCanonicalAASOP restrict_canonical( new RestrictAbsentCanonicalAAS( 0, aa_list ) );
		task->push_back(restrict_canonical);
	}

	PreventRepackingRLTOP norepack( new PreventRepackingRLT );
	OperateOnResidueSubsetOP op_subset( new OperateOnResidueSubset(norepack, residue_selector_, true) );
	task->push_back(op_subset);

	protocols::minimization_packing::PackRotamersMover packer;
	packer.score_function( sfxn );
	packer.task_factory( task );

	packer.apply( *pose_ );
	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has checked or unchecked the allow design checkbox.
void MainWindow::on_checkbox_allow_design_clicked()
{
	bool const checked( ui->checkbox_allow_design->isChecked() );
	ui->checkbox_ala->setEnabled( checked );
	ui->checkbox_met->setEnabled( checked );
	ui->checkbox_aliphatic->setEnabled( checked );
	ui->checkbox_gly->setEnabled( checked );
	ui->checkbox_pro->setEnabled( checked );
	ui->checkbox_aromatic->setEnabled( checked );
	ui->checkbox_charged->setEnabled( checked );
	ui->checkbox_polar->setEnabled( checked );
	ui->checkbox_his->setEnabled( checked );
	ui->checkbox_cys->setEnabled( checked );
}

/// @brief The user has clicked on the "run minimizer" button.
void MainWindow::on_button_run_minimizer_clicked()
{
	core::select::residue_selector::ResidueSubset selection( residue_selector_->apply(*pose_) );
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb(false);
	movemap->set_chi(false);
	movemap->set_jump(false);
	for(core::Size i(1), imax(selection.size()); i<=imax; ++i) {
		movemap->set_chi( i, selection[i] );
	}

	protocols::minimization_packing::MinMover minmove( movemap, core::scoring::get_score_function(), "lbfgs_armijo_nonmonotone", 0.0000001, false );
	minmove.apply(*pose_);
	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has selected that the pose be coloured by what's selected.  Update accordingly.
void MainWindow::on_radiobutton_colourby_selection_clicked()
{
	//ui->radiobutton_colourby_totalscore->setChecked(false);
	sfxn_ = nullptr;
	pose_draw_widget_->set_color_mode( ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget::ColorMode::selection );
	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has selected that the pose be coloured by total score (per residue).  Update accordingly.
void MainWindow::on_radiobutton_colourby_totalscore_clicked()
{
	//ui->radiobutton_colourby_selection->setChecked(false);
	sfxn_ = core::scoring::get_score_function();
	(*sfxn_)(*pose_); //Score the pose.
	pose_draw_widget_->set_color_mode( ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget::ColorMode::total_energy );
	pose_draw_widget_->update_pose_draw();
}

/// @brief When the user clicks the "clear imported pose" button, clear the imported pose.
void MainWindow::on_clearimported_button_clicked()
{
	nonparametric_pose_.reset();
	nonparametric_pose_ = nullptr;
	nonparametric_pose_offset_.x( 0.0 );
	nonparametric_pose_offset_.y( 0.0 );
	nonparametric_pose_offset_.z( 0.0 );

	pose_draw_widget_->reset_nonparametric_translation_vect();
	pose_draw_widget_->reset_nonparametric_rotation_xform();

	update_residue_selector();

	rebuild_pose_from_scratch();
}

/// @brief When the user clicks the "import pose" button, import a pose.
void MainWindow::on_importPDB_button_clicked()
{
	QString const filename( QFileDialog::getOpenFileName(this, "Load PDB file", "", "PDB file (*.pdb)") );
	if( filename.isEmpty() ) return;

	std::string pdb_filename( filename.toStdString() );

	nonparametric_pose_ = utility::pointer::make_shared< core::pose::Pose >();
	core::import_pose::pose_from_file( *nonparametric_pose_, pdb_filename, false, core::import_pose::PDB_file );

	nonparametric_pose_offset_.x( 0.0 );
	nonparametric_pose_offset_.y( 0.0 );
	nonparametric_pose_offset_.z( 0.0 );

	pose_draw_widget_->reset_nonparametric_translation_vect();
	pose_draw_widget_->reset_nonparametric_rotation_xform();

	update_residue_selector();

	rebuild_pose_from_scratch();
}

/// @brief When the user specifies that the mouse should be used for dragging the
/// nonparametric geometry around.
void MainWindow::on_tool_drag_nonparametric_clicked()
{
	pose_draw_widget_->set_extended_drag_mode( ui::ui_protocols::helical_bundle::HBPDOGLW_drag_nonparametric_geometry );
	pose_draw_widget_->determine_xyz_unit_vectors();
}

/// @brief When the user specifies that the mouse should be used for rotating the
/// viewing direction around.
void MainWindow::on_tool_rotate_zoom_view_clicked()
{
	pose_draw_widget_->set_drag_mode( ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget::DragMode::rotate_and_zoom_viewport );
	pose_draw_widget_->set_extended_drag_mode( ui::ui_protocols::helical_bundle::HBPDOGLW_base_class_mode );
}

/// @brief When the user specifies that the mouse should be used for rotating the
/// nonparametric geometry around.
void MainWindow::on_tool_rotate_nonparametric_clicked()
{
	pose_draw_widget_->set_extended_drag_mode( ui::ui_protocols::helical_bundle::HBPDOGLW_rotate_nonparametric_geometry );
	pose_draw_widget_->determine_xyz_unit_vectors();
	numeric::xyzVector< core::Real > comvect(0,0,0);
	compute_com_vect( comvect );
	pose_draw_widget_->set_nonparametric_COM_vect( comvect );
}


/// @brief When the user indicates that the nonparametric geometry has moved.
void
MainWindow::on_nonparametric_geometry_moved() {
	if(nonparametric_pose_ == nullptr) return; //Do nothing if there's no nonparametric geometry.

	//Get Pose references for speed:
	core::pose::Pose const & nppose( *nonparametric_pose_ );
	core::pose::Pose & pose( *pose_ );

	nonparametric_pose_offset_ = pose_draw_widget_->get_nonparametric_translation_vect(); //Get the translation vector.
	nonparametric_pose_rotation_ = pose_draw_widget_->get_nonparametric_rotation_xform(); //Get the rotation transform.

	for(core::Size ir(1), irmax(nppose.total_residue()); ir<=irmax; ++ir) { //Loop through all residues
		for(core::Size ia(1), iamax( nppose.residue(ir).natoms() ); ia<=iamax; ++ia ) { //Loop through all atoms
			core::id::AtomID const curat( ia, ir );
			pose.set_xyz( curat, nonparametric_pose_rotation_*nppose.xyz( curat ) + nonparametric_pose_offset_ );
		}
	}
	pose.update_residue_neighbors();

	if( sfxn_ != nullptr ) {
		(*sfxn_)(pose);
	}

	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has selected a preset in the drop-down list.
void MainWindow::on_presets_comboBox_currentIndexChanged(const QString &arg1)
{
	if( arg1 == "Select preset" ) { return; } //Do nothing if a preset wasn't picked.

	//Clear current geometry:
	core::Size const tabcount( ui->tabWidget->count() );
	if( tabcount > 0 ) {
		for(core::Size i(1); i<=tabcount; ++i) {
			QScrollArea* curscrollarea( dynamic_cast< QScrollArea*>( ui->tabWidget->widget(0) ) );
			ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget_ptr(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
			curwidget_ptr->on_remove_button_clicked();
		}
	}

	//Apply the presets:
	if( arg1 == "3-helix bundle" ) {
		//Add three helices:
		for(core::Size i(0); i<3; ++i) {
			on_add_helix_button_clicked();
			QScrollArea* curscrollarea( dynamic_cast< QScrollArea*>( ui->tabWidget->widget(i) ) );
			ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget_ptr(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
			curwidget_ptr->set_real_subwidget( protocols::helical_bundle::BPC_delta_omega0, static_cast< core::Real >(i * 120) );
			if(i > 0) { curwidget_ptr->set_everything_except_delta_omega0_copies_helix1(); }
			else {
				curwidget_ptr->set_real_subwidget( protocols::helical_bundle::BPC_r0, 6.5 );
				curwidget_ptr->set_real_subwidget( protocols::helical_bundle::BPC_omega0, 3.0 );
			}
			curwidget_ptr->set_invert( i % 2 == 1 );
		}
	}
	else if( arg1 == "Antiparallel beta-barrel" ) {
		//Add eight strands:
		for(core::Size i(0); i<8; ++i) {
			on_add_helix_button_clicked();
			QScrollArea* curscrollarea( dynamic_cast< QScrollArea*>( ui->tabWidget->widget(i) ) );
			ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget_ptr(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
			curwidget_ptr->set_real_subwidget( protocols::helical_bundle::BPC_delta_omega0, static_cast<core::Real>( i * 45 ) );
			curwidget_ptr->set_helix_length(7);
			curwidget_ptr->set_to_beta_strand();
			if(i > 0) { curwidget_ptr->set_everything_except_delta_omega0_copies_helix1(); }
			else {
				curwidget_ptr->set_real_subwidget( protocols::helical_bundle::BPC_r0, 6.5 );
				curwidget_ptr->set_real_subwidget( protocols::helical_bundle::BPC_omega0, 16.0 );
				curwidget_ptr->set_real_subwidget( protocols::helical_bundle::BPC_delta_omega1, 81.0 );
			}
			curwidget_ptr->set_invert( i % 2 == 1 );
		}
	}


	ui->presets_comboBox->setCurrentIndex(0);
	ui->presets_comboBox->update();
}
