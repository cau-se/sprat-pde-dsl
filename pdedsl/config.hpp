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

#ifndef CONFIG_HPP_
#define CONFIG_HPP_


#include <limits>
#include <cassert>
#ifdef SPRAT_BUILD_WITH_MPI
#include <mpi.h>
#endif

typedef unsigned int uint;
typedef uint64_t     index_t;
typedef int64_t      signed_index_t;
typedef uint16_t     short_index_t;
typedef int16_t      signed_short_index_t;
typedef double       real;

#ifdef SPRAT_BUILD_WITH_MPI
#define SPRAT_MPI_INDEX_T           MPI_UINT64_T
#define SPRAT_MPI_SHORT_INDEX_T     MPI_UINT16_T
#define SPRAT_MPI_REAL              MPI_DOUBLE
#endif //#ifdef SPRAT_BUILD_WITH_MPI

#define INDEX_T_NO_NEIGHBOR           (std::numeric_limits<index_t>::max())
#define INDEX_T_NO_LOCAL_NEIGHBOR     (std::numeric_limits<index_t>::max()-1)
#define INDEX_T_MAX_VALID             (std::numeric_limits<index_t>::max()-2)
//#define INDEX_T_NO_NEIGHBOR           UINT_MAX
//#define INDEX_T_NO_LOCAL_NEIGHBOR     (UINT_MAX-1)
//#define INDEX_T_MAX_VALID             (UINT_MAX-2)

#define runtimeAssert(COND) \
	assert(COND)


#define ARITHMETIC_MAX_ABS_DIFF		1e-14
#define ARITHMETIC_MAX_ULP_DIFF		4


#endif /* CONFIG_HPP_ */
