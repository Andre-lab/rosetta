	// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
	// vi: set ts=2 noet:
	//
	// This file is part of the Rosetta software suite and is made available under license.
	// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
	// (C) 199x-2009 Rosetta Commons participating institutions and developers.
	// For more information, see http://www.rosettacommons.org/.

	/// @brief initialization protocols for symmetrical docking
	/// @author Ingemar Andre


#ifndef INCLUDED_protocols_symmetric_docking_SymDockingInitialPerturbation_hh
#define INCLUDED_protocols_symmetric_docking_SymDockingInitialPerturbation_hh

// Package headers
// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
//#include <core/scoring/symmetry/SymmetricScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/RigidBodyMover.hh>

// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/kinematics/Jump.fwd.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

//utility
#include <utility/vector1.hh>

//Auto Headers
#include <core/types.hh>
#include <core/conformation/symmetry/SymDof.hh>


namespace protocols {
namespace symmetric_docking {

	class SymDockingInitialPerturbation; // fwd declaration
	typedef utility::pointer::owning_ptr< SymDockingInitialPerturbation > SymDockingInitialPerturbationOP;
	typedef utility::pointer::owning_ptr< SymDockingInitialPerturbation const > SymDockingInitialPerturbationCOP;

	class SymDockingInitialPerturbation : public moves::Mover
	{
	public:

	typedef core::conformation::symmetry::SymmetricConformation SymmetricConformation;
	typedef core::conformation::symmetry::SymmetryInfo SymmetryInfo;

	public:

		// default constructor
		SymDockingInitialPerturbation();

		// constructor with arguments
		SymDockingInitialPerturbation(
			bool const slide_in
		);

		// destructor
		~SymDockingInitialPerturbation();

		// protocol functions
		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const;

	private:
		/// do slide into context?
		bool slide_;
	};

	/// @brief Contrary to the name, slides things apart first, then together.
	class SymDockingSlideIntoContact : public moves::Mover
	{
	public:

		// default constructor
		SymDockingSlideIntoContact();

		// constructor with arguments
		SymDockingSlideIntoContact(
			std::map< Size, core::conformation::symmetry::SymDof > dofs
		);
		// destructor
		~SymDockingSlideIntoContact();

	// protocol functions
		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const;

	private:
		// dof to allow
		std::map< Size, core::conformation::symmetry::SymDof > dofs_;
		// scorefxn to use
		core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn_;
	};

	class FaSymDockingSlideTogether; // fwd declaration
	typedef utility::pointer::owning_ptr< FaSymDockingSlideTogether > FaSymDockingSlideTogetherOP;
	typedef utility::pointer::owning_ptr< FaSymDockingSlideTogether const > FaSymDockingSlideTogetherCOP;

	/// @brief Slides docking partners together by monitoring fa_rep.
	/// @details
	///		If partners are already touching, no change is made.
	///		Separation will be 1A or less after calling this function.
	class FaSymDockingSlideTogether : public moves::Mover
	{
	public:
		FaSymDockingSlideTogether(
			std::map< Size, core::conformation::symmetry::SymDof > dofs
		);

	 ~FaSymDockingSlideTogether();

		virtual void apply( core::pose::Pose & pose );
		virtual std::string get_name() const;

	private:
		 // dof to allow
		std::map< Size, core::conformation::symmetry::SymDof > dofs_;
		// scorefxn to use
		core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn_;
		core::Real tolerance_; ///< how accurate do you want to be?

};

	class SymmetrySlider; // fwd declaration
	typedef utility::pointer::owning_ptr< SymmetrySlider > SymmetrySliderOP;
	typedef utility::pointer::owning_ptr< SymmetrySlider const > SymmetrySliderCOP;

	class SymmetrySlider : public utility::pointer::ReferenceCount
	{
		public:
			SymmetrySlider(core::pose::Pose & pose);
			SymmetrySlider(
				core::pose::Pose & pose,
				core::conformation::symmetry::SlideCriteriaType score_criteria,
				std::string SlideCriteriaVal = "AUTOMATIC"
			);
			~SymmetrySlider();

			void setup( core::pose::Pose & pose);

			virtual void select_jump() = 0;
	//		virtual void select_dir( core::pose::Pose & pose ) = 0;
			// get functions
			core::Size get_current_jump();
			core::Real get_slide_threshold();
			std::map< core::Size, bool > get_allow_slide_jump_map() const;
			// set functions
			void set_current_jump( core::Size jump_nbr );
			void set_slide_criteria(std::string SlideCriteria );
			void set_slide_threshold( std::string threshold );

			core::Real slide_score( core::pose::Pose & pose );
			void slide_away( core::pose::Pose & pose );
			core::Real step_size();
			void disallow_current_slide();
			bool allowed_current_slide();

			bool continue_slide( core::pose::Pose & pose );
			bool very_far_away( core::pose::Pose & pose );
			bool finished();
			// this is where it happens...
			void slide(core::pose::Pose & pose);
			void apply(core::pose::Pose & pose);

			//fpd a helper functions to get the right direction for slide moves
			core::Real rg( core::pose::Pose const & pose );
			
			//fpd will the mover compress or expand the system?
			bool dofmover_compresses( core::pose::Pose & pose, protocols::moves::RigidBodyDofTransMover & dofmover );

		private:
			std::map< Size, core::conformation::symmetry::SymDof > dofs_;
			std::map< core::Size, bool > AllowSlideJumpMap_;
			std::map< core::Size, core::kinematics::Jump > InitialJumps_;
			std::map< core::Size, bool > InvertJump_;
			core::conformation::symmetry::SlideCriteriaType SlideCriteriaType_;
			core::scoring::symmetry::SymmetricScoreFunctionOP scorefxn_;
			std::string SlideThreshold_;
			bool reset_slide_;
			core::Size total_num_slides_;
			core::Size current_jump_;
	};

	// Derive your own slide class by specifying what jumps are selected in the
	// function select_jump()

	class SequentialSymmetrySlider; // fwd declaration
	typedef utility::pointer::owning_ptr< SequentialSymmetrySlider > SequentialSymmetrySliderOP;
	typedef utility::pointer::owning_ptr< SequentialSymmetrySlider const > SequentialSymmetrySliderCOP;

	class SequentialSymmetrySlider : public SymmetrySlider
	{
		public:
		SequentialSymmetrySlider( core::pose::Pose & pose );
		SequentialSymmetrySlider( SymmetrySlider const & Slider );
		SequentialSymmetrySlider(
			core::pose::Pose & pose,
			core::conformation::symmetry::SlideCriteriaType score_criteria,
			std::string SlideCriteriaVal = "AUTOMATIC"
		);
		void init();
		void select_jump();
		std::vector<core::Size> slide_order_;
//		void select_dir( core::pose::Pose & pose );

	};

	// Derive your own slide class by specifying what jumps are selected in the
	// function select_jump()

	class OrderedSequentialSymmetrySlider; // fwd declaration
	typedef utility::pointer::owning_ptr< OrderedSequentialSymmetrySlider > OrderedSequentialSymmetrySliderOP;
	typedef utility::pointer::owning_ptr< OrderedSequentialSymmetrySlider const > OrderedSequentialSymmetrySliderCOP;

	class OrderedSequentialSymmetrySlider : public SymmetrySlider
	{
		public:
		OrderedSequentialSymmetrySlider( core::pose::Pose & pose, std::vector<core::Size> slide_order );
		OrderedSequentialSymmetrySlider( SymmetrySlider const & Slider, std::vector<core::Size> slide_order );
		OrderedSequentialSymmetrySlider(
			core::pose::Pose & pose,
			core::conformation::symmetry::SlideCriteriaType score_criteria,
			std::string SlideCriteriaVal = "AUTOMATIC",
			std::vector<core::Size> slide_order = std::vector<core::Size>()
		);

		void select_jump();
//		void select_dir( core::pose::Pose & pose );
		private:
			std::vector<core::Size> slide_order_;

	};


	// Derive your own slide class by specifying what jumps are selected in the
	// function select_jump()
	class RandomSymmetrySlider; // fwd declaration
	typedef utility::pointer::owning_ptr< RandomSymmetrySlider > RandomSymmetrySliderOP;
	typedef utility::pointer::owning_ptr< RandomSymmetrySlider const > RandomSymmetrySliderCOP;

	class RandomSymmetrySlider : public SymmetrySlider
	{
		public:
		RandomSymmetrySlider( core::pose::Pose & pose );
		RandomSymmetrySlider( SymmetrySlider const & Slider );
		RandomSymmetrySlider(
			core::pose::Pose & pose,
			core::conformation::symmetry::SlideCriteriaType score_criteria,
			std::string SlideCriteriaVal = "AUTOMATIC"
		);

		void select_jump();
//		void select_dir( core::pose::Pose & pose );

	};


} // docking
} // protocols

#endif
