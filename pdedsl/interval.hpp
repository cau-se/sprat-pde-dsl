/*
 * Copyright 2014-2015 Arne Johanson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef INTERVAL_HPP_
#define INTERVAL_HPP_


#include <boost/iterator/iterator_facade.hpp>
#include "config.hpp"



struct Interval {
	real a;
	real b;

	constexpr Interval() : a(0.0), b(0.0) {}
	constexpr Interval(real inA, real inB) : a(inA), b(inB) {}
	constexpr real length() { return b-a; }
};


template<typename IntType = int>
struct IntRange {
	IntType a;
	IntType b; // one behind the last element

	constexpr IntRange() : a(0), b(0) {}
	constexpr IntRange(IntType inA, IntType inB) : a(inA), b(inB) {}
	constexpr IntRange(IntType inB) : a(0), b(inB) {}
	constexpr IntType length() { return b-a; }

	// ******************* Iterator
	class const_iterator : public boost::iterator_facade<
	const_iterator
	, IntType
	, boost::random_access_traversal_tag
	, IntType> {
	public:
		const_iterator() : _index(0) {}
		explicit const_iterator(IntType i) : _index(i) {}

	private:
		friend class boost::iterator_core_access;

		void increment() { ++_index; }
		void decrement() { --_index; }
		void advance(std::ptrdiff_t n) { _index+=n; }
		std::ptrdiff_t distance_to(const_iterator const& other) const { return (std::ptrdiff_t)other._index - (std::ptrdiff_t)_index; }

		bool equal(const_iterator const& other) const {	return _index == other._index; }

		IntType dereference() const { return _index; }

		IntType _index;
	};
	typedef const_iterator iterator;

	iterator begin() const { return iterator(a); }
	iterator end() const { return iterator(b); }
	const_iterator cbegin() const { return const_iterator(a); }
	const_iterator cend() const { return const_iterator(b); }
};

typedef IntRange<index_t> IndexRange;
typedef IntRange<uint> UIntRange;



#endif /* INTERVAL_HPP_ */

