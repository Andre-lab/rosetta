// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/Tracer.hh
/// @brief  Tracer IO system
/// @author Sergey Lyskov
///

#ifndef INCLUDED_basic_Tracer_hh
#define INCLUDED_basic_Tracer_hh

#include <basic/Tracer.fwd.hh>
// AUTO-REMOVED #include <utility/stream_util.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <sstream>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/vector1_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <ostream>


namespace basic {


extern std::ostream *final_channel;

/// @brief
/// Priority levels for T() and Tracer object, modeled on the log4j project and its offspring.
/// Priorities in Tracer are still ints so users can pass other arbitrary integer values (for now).
enum TracerPriority {
	t_fatal   = 0,   //< The FATAL level designates very severe error events that will presumably lead the application to abort.
	t_error   = 100, //< The ERROR level designates error events that might still allow the application to continue running.
	t_warning = 200, //< The WARN level designates potentially harmful situations.
	t_info    = 300, //< The INFO level designates informational messages that highlight the progress of the application at coarse-grained level.
	t_debug   = 400, //< The DEBUG level designates fine-grained informational events that are most useful to debug an application.
	t_trace   = 500  //< The TRACE level designates finer-grained informational events than the DEBUG level.
};




/// @brief Base class for Tracer, TracerProxy and UTracer objects.
template <class CharT, class Traits = std::char_traits<CharT> >
class basic_otstream : public std::basic_ostream<CharT, Traits>, public utility::pointer::ReferenceCount
{
protected: /// Inner class declaration
	/// @brief Wrapper class for std::stringbuf
	template <class _CharT, class _Traits = std::char_traits<_CharT> >
	class basic_tstringbuf : public std::basic_stringbuf<_CharT, _Traits> {
	public:
		basic_tstringbuf(basic_otstream *ot) : otsream_(ot) {}
		virtual ~basic_tstringbuf() {}

	protected:
		virtual int sync() {
			otsream_->t_flush( this->str() ); //std::basic_stringbuf<CharT, Traits>::str() );
			//std::basic_stringbuf<CharT, Traits>::str("");
			this->str("");
			return 0;
		}
	private:
		basic_otstream *otsream_;
	};


public:
	basic_otstream() : std::basic_ostream<CharT, Traits> ( new basic_tstringbuf<CharT, Traits> (this) ) {}
	virtual ~basic_otstream() { delete this->rdbuf(); }


	/// @brief Return true if inner string buffer is empty.
	bool is_flushed() const {
			basic_tstringbuf<char> * buf = dynamic_cast< basic_tstringbuf<char> * >( this->rdbuf() );
			return buf->str().size() == 0;
	}


protected:
	/// @brief notification that flush function was called and inner buffer should be outputed.
	virtual void t_flush(std::string const &) { assert("basic_otstream::t_flush"); };

private:
	basic_otstream(basic_otstream const & );


	/// Data members
	/// @brief inner string buffer
	//std::basic_stringbuf<CharT, Traits> * tstringbuf_;
};


typedef basic_otstream<char> otstream;

typedef utility::pointer::owning_ptr< otstream > otstreamOP;


/// @brief data structure to store all system level options for Tracer system.
struct TracerOptions
{
	/// @brief system priority level
	int level;

	/// @brief should channel name be printed during the IO?
	bool print_channel_name;

	/// @brief list of muted channels
	utility::vector1<std::string> muted;

	/// @brief list of unmuted channels
	utility::vector1<std::string> unmuted;

	#ifdef EXPERIMENTAL_TRACER_FEATURES
		/// @brief channel is muted for all but error-level
		utility::vector1<std::string> muted_warning;

		/// @brief channel is muted for all but warning and error level
		utility::vector1<std::string> muted_info;

		/// @brief and so on...
		utility::vector1<std::string> muted_debug;
		utility::vector1<std::string> muted_trace;

		/// @brief channel is unmuted for error level
		utility::vector1<std::string> unmuted_error;

		/// @brief channel is unmuted for error and warning level
		utility::vector1<std::string> unmuted_warning;

		/// @brief and so on...
		utility::vector1<std::string> unmuted_info;
		utility::vector1<std::string> unmuted_debug;
	#endif // EXPERIMENTAL_TRACER_FEATURES
};




/// @brief Class for handling user debug/warnings/errors.
///  Use instance of this class instead of 'std::cout' for all your regular io.
///  Channel argument must be related to the location of the source file. For example if you create
///  Tracer object in src/basic/scoring/myfile.cc,
///    then channel must be something like 'src.basic.scoring.myfile'
class Tracer :  public otstream
{

public:

	/// @brief Create Tracer object with given channel and priority
	Tracer(std::string const & channel="", TracerPriority priority=t_info, bool muted_by_default=false);

	virtual ~Tracer();

	/// @brief re-init using data from another tracer object.
	void init(Tracer const & tr);

	/// @brief flush tracer buffer and flush buffers of all
	///        sub-channels ie: Fatal, Error, Warning, Info, Debug, Trace
	void flush_all_channels();

	/// @brief set ios hook for all tracer io operation.
	/// @param monitoring_channels_list is space separated list of channels.
	static void set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list);

	static void set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list, bool raw);

	static std::string const AllChannels;

	/// @brief Is this tracer currently visible?.
	bool visible() const;

	/// @brief is this tracer visible, if it used the given priority value?
	bool visible( int priority ) const;

	/// @brief get/set tracer priority level.
	int priority() { return priority_; }
	Tracer & operator ()(int priority);
	void priority(int priority);

	std::string const& channel() { return channel_; }

	/// @brief get/set tracer options - global options for Tracer IO.
	static TracerOptions & tracer_options() { return tracer_options_; }

	/// @brief global super mute flag that allow to mute all io no matter what.
	static bool super_mute() { return super_mute_; }
	static void super_mute(bool f) { super_mute_ = f; }

	static void flush_all_tracers();

public: /// Inner Classes
	 /// @brief Small inner class acting as a proxy to an object that hold it.
	class TracerProxy : public otstream // std::ostringstream //
	{
	public:
		TracerProxy(Tracer &tracer, int priority, std::string const &channel) :
			tracer_(tracer), priority_(priority), channel_(channel), visible_(true), visibility_calculated_(false) {}

		virtual ~TracerProxy();

		bool visible() const;

	protected:
		virtual void t_flush(std::string const &);

	private:
		Tracer & tracer_;
		int priority_;

		/// @brief We need to copy channel name here so we can generate appropriate 'warning' message
		/// in destructor, where tracer_ object is no longer valid.
		std::string channel_;

		/// @brief is channel visible?
		mutable bool visible_;

		/// @brief is channel visibility already calculated?
		mutable bool visibility_calculated_;
	};

	/// @brief channels with predefined priority levels.
	TracerProxy Fatal, Error, Warning, Info, Debug, Trace;


protected:
	/// @brief overload member function.
	virtual void t_flush(std::string const &);


private: /// Functions
	/// @brief copy constructor.
	Tracer(Tracer const & tr);


	/// @brief return true if channel is inside vector, some logic apply.
	static bool in(utility::vector1<std::string> const &, std::string const channel, bool strict);

	template <class out_stream>
	void prepend_channel_name( out_stream& sout, std::string const &str);

	/// @brief calcualte visibility of the current object depending of the channel name and priority.
	void calculate_visibility(void) const;
	static void calculate_visibility(std::string const &channel, int priority, bool &visible, bool &muted, bool muted_by_default);

	#ifdef EXPERIMENTAL_TRACER_FEATURES
		static void experimental_calculate_visibility(std::string const &channel, int priority, bool &visible, int &mute_level, bool muted_by_default);
	#endif // EXPERIMENTAL_TRACER_FEATURES


private: /// Data members

	/// @brief channel name
	std::string channel_;

	/// @brief channel priority level
	int priority_;

	/// @brief is channel visible?
	mutable bool visible_;

	/// @brief is channel muted ?
	mutable bool muted_;

	/// @brief above which level is channel muted ?
	mutable int mute_level_;

	/// @brief is channel muted by default?
	bool muted_by_default_;

	/// @brief should channel name be printed during the io?
	//mutable bool print_channel_name_;

	/// @brief is current printing position a begining of the line?
	bool begining_of_the_line_;

	/// @brief is channel visibility already calculated?
	mutable bool visibility_calculated_;

	/// @brief system priority level
	//mutable int level_;

	/// static data members
	/// @brief link to Tracer like object where all output for selecting channels should go.
	static otstreamOP ios_hook_;

	/// @brief should the ios_hook_ the raw output?
	static bool ios_hook_raw_;

	/// @brief list of channels for which outout should be redirected.
	static utility::vector1<std::string> monitoring_list_;


	/// @brief global option collection for Tracer IO.
	static TracerOptions tracer_options_;


	/// @brief global super mute flag that allow to mute all io no matter what.
	static bool super_mute_;

	/// @which Mpi rank is this process
	static int mpi_rank_;

	/// @brief static collection of all Tracer objects
	static std::vector< Tracer * > & all_tracers();

	/// @brief T is special function for assign tracer property on the static object.
	friend Tracer & T(std::string const &, TracerPriority);
};


/// @brief T is special function for assign tracer property on the static object.
Tracer & T(std::string const & channel, TracerPriority priority=t_info);

/// @brief Predefined Error tracer.
inline Tracer & Error(TracerPriority priority=t_error) { return T("Error", priority); }

/// @brief Predefined Warning tracer.
inline Tracer & Warning(TracerPriority priority=t_warning) { return T("Warning", priority); }




} // namespace basic

#endif // INCLUDED_basic_tracer_hh

