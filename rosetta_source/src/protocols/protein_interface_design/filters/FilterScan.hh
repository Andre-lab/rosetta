// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/FilterScanFilter.hh
/// @brief Scans a task factory, mutates to each designable residue and evaluates a filter. Mutations that pass the filter are output in a resfile format
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_filters_FilterScanFilter_hh
#define INCLUDED_protocols_protein_interface_design_filters_FilterScanFilter_hh


// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/protein_interface_design/filters/FilterScan.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/vector1.hh>

// Unit headers

namespace protocols {
namespace protein_interface_design{
namespace filters {

class FilterScanFilter : public protocols::filters::Filter
{
private:
	typedef protocols::filters::Filter parent;
public:
	/// @brief default ctor
	FilterScanFilter();
	virtual ~FilterScanFilter();

	///@brief Constructor with a single target residue
	virtual bool apply( core::pose::Pose const & pose ) const;
	virtual void report( std::ostream & out, core::pose::Pose const & pose ) const;
	virtual core::Real report_sm( core::pose::Pose const & pose ) const;
	virtual protocols::filters::FilterOP clone() const;
	virtual protocols::filters::FilterOP fresh_instance() const;
	/// Undefined, commenting out to fix PyRosetta build  core::Real compute( core::pose::Pose const & pose ) const;
	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	void resfile_name( std::string const resfile_name);
	std::string resfile_name() const;
	protocols::filters::FilterOP triage_filter() const;
	void triage_filter( protocols::filters::FilterOP filter );

	protocols::filters::FilterOP filter() const;
	void filter( protocols::filters::FilterOP filter );
	std::string resfile_general_property() const;
	void resfile_general_property( std::string const );
	void relax_mover( protocols::moves::MoverOP  mover );
	protocols::moves::MoverOP relax_mover() const;
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	void parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	bool unbound() const;
	void unbound( bool const u );
	core::Size jump() const;
	void jump( core::Size const j );
	bool delta() const;
	void delta( bool const d );
	bool report_all() const;
	void report_all( bool const ra );
	void dump_pdb( bool const d );
	bool dump_pdb() const;
private:
	core::pack::task::TaskFactoryOP task_factory_;
	protocols::filters::FilterOP triage_filter_;//dflt null; mutations that are categorically rejected, no matter what
	protocols::filters::FilterOP filter_;//dflt null; a filter to use for its report functionality to probe the pose's state
	std::string resfile_name_;
	std::string resfile_general_property_; //dflt nataa; what to write in the resfile above the 'start' line
	protocols::moves::MoverOP relax_mover_; //dflt nullmover; what to do after mutation
	core::scoring::ScoreFunctionOP scorefxn_; // which scorefxn to use during packing
	bool delta_; // dflt false; compute as delta? If true, all values are reported relative to the baseline input pose's filter evaluation.
	bool unbound_;
	bool report_all_;
	core::Size jump_;
	void unbind( core::pose::Pose & ) const; //utility function for unbinding the pose
	bool dump_pdb_; // dflt false; dump a pdb for each substitution (with extensions signifying the substitution).
};

} // filters
} //protein_interface_design
} // protocols

#endif //INCLUDED_protocols_Filters_FilterScanFilter_HH_

