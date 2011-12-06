// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/canonical_sampling/PDBTrajectoryRecorder.cc
///
/// @brief
/// @author


// Unit header or inline function header
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorderCreator.hh>

// Other project headers or inline function headers
#include <core/io/raw_data/ScoreStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>  // required for Windows build
#include <protocols/canonical_sampling/TemperingBase.hh>
#include <protocols/jd2/ScoreMap.hh>
#include <protocols/jd2/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <utility/tag/Tag.hh>

// External library headers

// C++ headers
#include <iomanip>

// Operating system headers

// Forward declarations


namespace protocols {
namespace canonical_sampling {

std::string
PDBTrajectoryRecorderCreator::keyname() const {
	return PDBTrajectoryRecorderCreator::mover_name();
}

protocols::moves::MoverOP
PDBTrajectoryRecorderCreator::create_mover() const {
	return new PDBTrajectoryRecorder;
}

std::string
PDBTrajectoryRecorderCreator::mover_name() {
	return "PDBTrajectoryRecorder";
}

PDBTrajectoryRecorder::PDBTrajectoryRecorder()
{
	file_name("traj.pdb");
}

PDBTrajectoryRecorder::~PDBTrajectoryRecorder() {}

PDBTrajectoryRecorder::PDBTrajectoryRecorder(
	PDBTrajectoryRecorder const & other
) : protocols::canonical_sampling::TrajectoryRecorder( other ) {}

PDBTrajectoryRecorder&
PDBTrajectoryRecorder::operator=( PDBTrajectoryRecorder const & /* other */ )
{
	// assignment not allowed
	runtime_assert(false);
	return *this;
}

protocols::moves::MoverOP
PDBTrajectoryRecorder::clone() const
{
	return new protocols::canonical_sampling::PDBTrajectoryRecorder( *this );
}

protocols::moves::MoverOP
PDBTrajectoryRecorder::fresh_instance() const
{
	return new PDBTrajectoryRecorder;
}

std::string
PDBTrajectoryRecorder::get_name() const
{
	return "PDBTrajectoryRecorder";
}

void
PDBTrajectoryRecorder::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap& data /* data */,
	protocols::filters::Filters_map const& filter /* filters */,
	protocols::moves::Movers_map const& mover /* movers */,
	core::pose::Pose const& pose/* pose */
)
{
	Parent::parse_my_tag( tag, data, filter, mover, pose );
}

void
PDBTrajectoryRecorder::reset(
	protocols::moves::MonteCarlo const & mc,
	protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover //= 0
)
{
	Parent::reset(mc, metropolis_hastings_mover);
	write_model(mc.last_accepted_pose(), metropolis_hastings_mover);
}

void
PDBTrajectoryRecorder::write_model(
	core::pose::Pose const & pose,
	protocols::canonical_sampling::MetropolisHastingsMoverCAP metropolis_hastings_mover //= 0
)
{
	if (trajectory_stream_.filename() == "") {

		std::string filename( metropolis_hastings_mover ? metropolis_hastings_mover->output_file_name(file_name(), cumulate_replicas()) : file_name() );

		trajectory_stream_.open( filename );
	}

	core::Size replica = protocols::jd2::current_replica();

	TemperingBaseCAP tempering = 0;
	if (metropolis_hastings_mover) {
		tempering = dynamic_cast< TemperingBase const * >( metropolis_hastings_mover->tempering()() );
	}

	std::map < std::string, core::Real > score_map;
	std::map < std::string, std::string > string_map;
	protocols::jd2::ScoreMap::score_map_from_scored_pose(score_map, pose);
	core::io::raw_data::ScoreStruct score_struct;

	trajectory_stream_ << "MODEL     " << std::setw(4) << model_count() << std::endl;
	trajectory_stream_ << "REMARK  99 Trial: " << step_count() << std::endl;
	if (replica) trajectory_stream_ << "REMARK  99 Replica: " << replica << std::endl;
	if (tempering) trajectory_stream_ << "REMARK  99 Temperature: " << metropolis_hastings_mover->monte_carlo()->temperature() << std::endl;
	trajectory_stream_ << "REMARK  99 Score: " << pose.energies().total_energy() << std::endl;
	trajectory_stream_ << "REMARK  99 ";
	score_struct.print_header(trajectory_stream_, score_map, string_map, false);
	trajectory_stream_ << "REMARK  99 ";
	score_struct.print_scores(trajectory_stream_, score_map, string_map);

	pose.dump_pdb(trajectory_stream_);
	trajectory_stream_ << "ENDMDL" << std::endl;
	trajectory_stream_.flush();
}

void
PDBTrajectoryRecorder::finalize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover
)
{
	Parent::finalize_simulation( pose, metropolis_hastings_mover );
	trajectory_stream_.close();
}


} // namespace canonical_sampling
} // namespace protocols
