// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief unresizable vector whose size is known at compile time,
/// which may be allocated on the stack, and which indexes from 1.

#ifndef INCLUDED_utility_fixedsizearray1_hh
#define INCLUDED_utility_fixedsizearray1_hh

// Unit headers
#include <utility/fixedsizearray1.fwd.hh>

// C++ headers
#include <iterator>
#include <cassert>

// #ifdef USEBOOSTSERIALIZE
// #include ???
// #endif

namespace utility {

/// Requirements:
/// S must be a positive integer
/// T must be interpretable as 0

template< typename T, platform::Size S >
class fixedsizearray1iterator {
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef ptrdiff_t                         difference_type;
	typedef T                                      value_type;
	typedef T *                                       pointer;
	typedef T &                                     reference;
	typedef ptrdiff_t                                distance;

public:
	fixedsizearray1iterator( T * array_ptr, T * position ) :
		array_ptr_( array_ptr ),
		position_( position )
	{}

	fixedsizearray1iterator( fixedsizearray1iterator< T, S > const & other ) :
		array_ptr_( other.array_ptr_ ),
		position_( other.position_ )
	{}

	fixedsizearray1iterator< T, S >
	operator = ( fixedsizearray1iterator< T, S > const & rhs ) {
		array_ptr_ = rhs.array_ptr_;
		position_ = rhs.position_;
        return *this;
	}

	T & operator * () {
		assert( valid() );
		return * position_;
	}

	fixedsizearray1iterator< T, S >
	operator ++ () {
		assert( valid() );
		position_ += 1;
		return *this;
	}

	/// @brief random access iterator jump by d
	fixedsizearray1iterator< T, S >
	operator + ( distance d ) {
		assert( valid() );
		return fixedsizearray1iterator< T, S >( array_ptr_, position_ + d );
	}

	/// @brief random access iterator jump by -d
	fixedsizearray1iterator< T, S >
	operator - ( distance d ) {
		assert( valid() );
		return fixedsizearray1iterator< T, S >( array_ptr_, position_ - d );
	}

	distance
	operator - ( fixedsizearray1const_iterator< T, S > const & other ) const;

	/// @brief random access increment
	fixedsizearray1iterator< T, S > const &
	operator += ( distance d ) {
		assert( valid() );
		position_ += d;
		return this;
	}

	/// @brief random access decrement
	fixedsizearray1iterator< T, S > const &
	operator -= ( distance d ) {
		assert( valid() );
		position_ -= d;
		return this;
	}


	friend
	inline
	bool
	operator <( fixedsizearray1iterator< T, S > const & a,
	            fixedsizearray1iterator< T, S > const & b )
	{
		return ( a.position_ < b.position_ );
	}

	bool
	operator == ( fixedsizearray1iterator< T, S > const & rhs ) {
		return position_ == rhs.position_ && array_ptr_ == rhs.array_ptr_;
	}

	bool
	operator != ( fixedsizearray1iterator< T, S > const & rhs ) {
		return position_ != rhs.position_ || array_ptr_ != rhs.array_ptr_;
	}

	/// @brief < comparison
	bool
	operator < ( fixedsizearray1iterator< T, S > const & rhs ) {
		return position_ < rhs.position_;
	}

	/// @brief <= comparison
	bool
	operator <= ( fixedsizearray1iterator< T, S > const & rhs ) {
		return position_ <= rhs.position_;
	}

	/// @brief > comparison
	bool
	operator > ( fixedsizearray1iterator< T, S > const & rhs ) {
		return position_ > rhs.position_;
	}

	/// @brief >= comparison
	bool
	operator >= ( fixedsizearray1iterator< T, S > const & rhs ) {
		return position_ >= rhs.position_;
	}

	friend class fixedsizearray1const_iterator< T, S >;

protected:

	bool
	valid() const {
		return array_ptr_ + S >  position_ && position_ >= array_ptr_;
	}

private:
	T * array_ptr_;
	T * position_;
};

template< typename T, platform::Size S >
class fixedsizearray1const_iterator {
public:
	typedef std::random_access_iterator_tag iterator_category;
	typedef ptrdiff_t                         difference_type;
	typedef T                                      value_type;
	typedef T *                                       pointer;
	typedef T &                                     reference;
	typedef ptrdiff_t                                distance;

public:
	fixedsizearray1const_iterator( T const * array_ptr, T const * position ) :
		array_ptr_( array_ptr ),
		position_( position )
	{}

	fixedsizearray1const_iterator( fixedsizearray1const_iterator< T, S > const & other ) :
		array_ptr_( other.array_ptr_ ),
		position_( other.position_ )
	{}

	fixedsizearray1const_iterator( fixedsizearray1iterator< T, S > const & other ) :
		array_ptr_( other.array_ptr_ ),
		position_( other.position_ )
	{}


	fixedsizearray1const_iterator< T, S >
	operator = ( fixedsizearray1const_iterator< T, S > const & rhs ) {
		array_ptr_ = rhs.array_ptr_;
		position_ = rhs.position_;
	}

	T const & operator * () {
		assert( valid() );
		return * position_;
	}

	fixedsizearray1const_iterator< T, S >
	operator ++ () {
		assert( valid() );
		position_ += 1;
		return *this;
	}

	/// @brief random access iterator jump by d
	fixedsizearray1const_iterator< T, S >
	operator + ( distance d ) {
		assert( valid() );
		return fixedsizearray1iterator< T, S >( array_ptr_, position_ + d );
	}

	/// @brief random access iterator jump by -d
	fixedsizearray1const_iterator< T, S >
	operator - ( distance d ) {
		assert( valid() );
		return fixedsizearray1iterator< T, S >( array_ptr_, position_ - d );
	}


	friend
	inline
	bool
	operator <( fixedsizearray1const_iterator< T, S > const & a,
	            fixedsizearray1const_iterator< T, S > const & b )
	{
		return ( a.position_ < b.position_ );
	}

	//// @brief distance between two iterators
	distance
	operator - ( fixedsizearray1const_iterator< T, S > const & other ) const {
		return position_ - other.position_;
	}

	/// @brief random access increment
	fixedsizearray1const_iterator< T, S > const &
	operator += ( distance d ) {
		assert( valid() );
		position_ += d;
		return this;
	}

	/// @brief random access decrement
	fixedsizearray1const_iterator< T, S > const &
	operator -= ( distance d ) {
		assert( valid() );
		position_ -= d;
		return this;
	}

	bool
	operator == ( fixedsizearray1const_iterator< T, S > const & rhs ) {
		return position_ == rhs.position_ && array_ptr_ == rhs.array_ptr_;
	}

	bool
	operator != ( fixedsizearray1const_iterator< T, S > const & rhs ) {
		return position_ != rhs.position_ || array_ptr_ != rhs.array_ptr_;
	}

	/// @brief < comparison
	bool
	operator < ( fixedsizearray1const_iterator< T, S > const & rhs ) {
		return position_ < rhs.position_;
	}

	/// @brief <= comparison
	bool
	operator <= ( fixedsizearray1const_iterator< T, S > const & rhs ) {
		return position_ <= rhs.position_;
	}

	/// @brief > comparison
	bool
	operator > ( fixedsizearray1const_iterator< T, S > const & rhs ) {
		return position_ > rhs.position_;
	}

	/// @brief >= comparison
	bool
	operator >= ( fixedsizearray1const_iterator< T, S > const & rhs ) {
		return position_ >= rhs.position_;
	}

protected:

	bool
	valid() const {
		return array_ptr_ + S >  position_ && position_ >= array_ptr_;
	}


private:
	T const * array_ptr_;
	T const * position_;
};

//template < typename T, platform::Size S >
//fixedsizearray1iterator< T, S >::distance
//fixedsizearray1iterator< T, S >::operator - (
//	fixedsizearray1const_iterator< T, S > const & other
//) const {
//	return position_ - other.position_;
//}


template < typename T, platform::Size S >
class fixedsizearray1
{
public:
	typedef platform::Size Size;
	typedef platform::SSize SSize;

	typedef T value_type;
	typedef fixedsizearray1iterator< T, S > iterator;
	typedef fixedsizearray1const_iterator< T, S > const_iterator;

	//typedef  typename T &        reference;
	//typedef  typename T const &  const_reference;

public:
	/// Constructors and the assigmnet operator

	fixedsizearray1() {
		for ( Size ii = 0; ii < S; ++ii ) array_[ ii ] = value_type ( 0 );
	}

	fixedsizearray1( value_type def ) {
		for ( Size ii = 0; ii < S; ++ii ) array_[ ii ] = def;
	}

	fixedsizearray1( fixedsizearray1< T, S > const & source ) {
		for ( Size ii = 0; ii < S; ++ii ) array_[ ii ] = source.array_[ ii ];
	}

	fixedsizearray1< T, S > const &
	operator = ( fixedsizearray1< T, S > const & rhs ) {
		for ( Size ii = 0; ii < S; ++ii ) array_[ ii ] = rhs.array_[ ii ];
		return *this;
	}

#ifdef USEBOOSTSERIALIZE
public:
	template<class Archive> 
	void 
	serialize(Archive & ar, const unsigned int /*version*/){
		for( Size ii = 0; ii < S; ++ii ) ar & array_[ii];
	}
#endif


public:
	/// Mutators and accessors

	value_type &
	operator [] ( Size index ) {
		assert( range( index ) );
		return array_[ index - 1 ];
	}

	value_type const &
	operator [] ( Size index ) const  {
		assert( range( index ) );
		return array_[ index - 1 ];
	}


	Size
	size() const {
		return S;
	}

public:
	/// Iterators
	iterator
	begin() {
		return fixedsizearray1iterator< T, S >( array_, array_ );
	}

	iterator
	end() {
		return fixedsizearray1iterator< T, S >( array_, array_ + S );
	}

	const_iterator
	begin() const {
		return fixedsizearray1const_iterator< T, S >( array_, array_ );
	}

	const_iterator
	end() const {
		return fixedsizearray1const_iterator< T, S >( array_, array_ + S);
	}

protected:

	bool
	range( Size index ) const {
		return index > (Size) 0 && index <= (Size) S;
	}

private:
	/// Data
	T array_[ S ];
};


} // namespace utility

#endif
