// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author  Ingemar Andre

#ifndef INCLUDED_protocols_moves_symmetry_SymRotamerTrialsMover_hh
#define INCLUDED_protocols_moves_symmetry_SymRotamerTrialsMover_hh

// Unit headers
#include <protocols/moves/symmetry/SymRotamerTrialsMover.fwd.hh>
#include <protocols/moves/RotamerTrialsMover.fwd.hh>
#include <protocols/moves/RotamerTrialsMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

namespace protocols {
namespace moves {
namespace symmetry {

class SymRotamerTrialsMover : public RotamerTrialsMover {
public:

typedef core::conformation::symmetry::SymmetricConformation SymmetricConformation;
typedef core::conformation::symmetry::SymmetryInfo SymmetryInfo;

public:

	// default constructor
	SymRotamerTrialsMover();

	///@brief constructor with PackerTask. use a PackerTask ONLY for fixed-sequence work.
	/// WARNING TO ANY DESIGNER WHO PASSES IN A TASK: YOUR DESIGN STEPS WILL BE UNDONE
	/// AS THIS TASK CONCEIVES OF THE INPUT SEQUENCE THAT CORRESPONDS TO THE ORIGINAL SEQUENCE
	SymRotamerTrialsMover(
		ScoreFunctionCOP scorefxn_in,
		PackerTask & task_in
	);

	///@brief constructor with TaskFactory
	SymRotamerTrialsMover(
		ScoreFunctionCOP scorefxn_in,
		TaskFactoryCOP factory_in
	);

	~SymRotamerTrialsMover();

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	void
	make_symmetric_task(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskOP task
	);

};

class SymEnergyCutRotamerTrialsMover : public SymRotamerTrialsMover {
public:

typedef core::conformation::symmetry::SymmetricConformation SymmetricConformation;
typedef core::conformation::symmetry::SymmetryInfo SymmetryInfo;


public:

	// default constructor
	SymEnergyCutRotamerTrialsMover();

	// constructor with arguments
	SymEnergyCutRotamerTrialsMover(
		ScoreFunctionCOP scorefxn_in,
		PackerTask & task_in,
		MonteCarloOP mc_in,
		core::Real energycut_in
	);

	// constructor with arguments
	SymEnergyCutRotamerTrialsMover(
		ScoreFunctionCOP scorefxn_in,
		TaskFactoryCOP factory_in,
		MonteCarloOP mc_in,
		core::Real energycut_in
	);

	~SymEnergyCutRotamerTrialsMover();

	void
	make_symmetric_task(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskOP task
	);

public:

	/// @brief apply this mover to a pose
	virtual
	void
	apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

protected:

	/// @brief selects a subset of residues to repack based on the per
	/// residue energies of the last accepted pose in the MC object.
	void
	setup_energycut_task(
		core::pose::Pose const & pose,
		MonteCarlo const & mc,
		core::pack::task::PackerTask & task_in
	) const;

	MonteCarloOP
	mc();

private:

	// data
	MonteCarloOP mc_;
	core::Real energycut_;

};

} // symmetry
} // moves
} // protocols

#endif

