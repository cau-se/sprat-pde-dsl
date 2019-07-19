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

#ifndef UTIL_HPP_
#define UTIL_HPP_

#include <iostream>
//#include <iomanip>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <boost/math/special_functions/next.hpp>
#include "config.hpp"


template <typename T>
inline void printVector(std::vector<T> const& vec) {
	//std::cout << std::setprecision(precision);
	for(uint i=0; i<vec.size(); ++i) {
		std::cout << vec[i];
		if(i!=vec.size()-1) {
			std::cout << ", ";
		}
	}
	std::cout << std::endl;
}

// cf. http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
static inline bool almostEqual(real a, real b) {
	// Check if the numbers are really close -- needed when comparing numbers near zero.
	if(fabs(a - b) <= ARITHMETIC_MAX_ABS_DIFF) {
		return true;
	}

	//if (uA.Negative() != uB.Negative())
	//	return false;

	if(fabs(boost::math::float_distance(a, b)) <= ARITHMETIC_MAX_ULP_DIFF) {
		return true;
	}

	return false;
}

static inline real max(real a, real b) {
	return (a>=b ? a : b);
}

static inline real min(real a, real b) {
	return (a<b ? a : b);
}

static inline real max(real a, real b, real c) {
	return (a<b ? (b<c ? c : b) : (a<c ? c : a));
}

static inline real min(real a, real b, real c) {
	return (a>b ? (b>c ? c : b) : (a>c ? c : a));
}


#endif /* UTIL_HPP_ */
