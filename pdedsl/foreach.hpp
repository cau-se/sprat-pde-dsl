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

#ifndef FOREACH_HPP_
#define FOREACH_HPP_

#include <boost/preprocessor/cat.hpp>
#ifdef SPRAT_BUILD_WITH_OPENMP
#include <omp.h>
#endif


#define FOREACH_ID(x) BOOST_PP_CAT(x, __LINE__) // Boost.Foreach: problems with MSVC?

#define foreach(VAR, CONTAINER, BODY) \
	{ \
		auto const& FOREACH_ID(fe_container) = CONTAINER; \
		const auto FOREACH_ID(fe_container_cend) = FOREACH_ID(fe_container).cend(); \
		for(auto FOREACH_ID(fe_it) = FOREACH_ID(fe_container).cbegin(); FOREACH_ID(fe_it) != FOREACH_ID(fe_container_cend); ++FOREACH_ID(fe_it)) { \
			VAR = *FOREACH_ID(fe_it); \
			BODY \
		} \
	}


#ifdef SPRAT_BUILD_WITH_OPENMP
#define FOREACH_STRINGIFY(x) #x
#define FOREACH_OMP_HELPER(x) FOREACH_STRINGIFY(omp parallel x)


#define foreach_omp_index(ITERVAR, NELEMENTS, OMPEXPR, BODY) \
	_Pragma(FOREACH_OMP_HELPER(OMPEXPR)) \
	{ \
		const size_t FOREACH_ID(fe_nElements) = NELEMENTS; \
		const unsigned int FOREACH_ID(fe_thread_count) = omp_get_num_threads(); \
		const unsigned int FOREACH_ID(fe_thread_id) = omp_get_thread_num(); \
		const size_t FOREACH_ID(fe_chunk_size) = FOREACH_ID(fe_nElements) / FOREACH_ID(fe_thread_count); \
		const auto FOREACH_ID(fe_index_begin) = FOREACH_ID(fe_thread_id) * FOREACH_ID(fe_chunk_size); \
		const auto FOREACH_ID(fe_index_end) = ((FOREACH_ID(fe_thread_id) == FOREACH_ID(fe_thread_count) - 1 || FOREACH_ID(fe_nElements) < FOREACH_ID(fe_thread_count)) ? FOREACH_ID(fe_nElements) : FOREACH_ID(fe_index_begin) + FOREACH_ID(fe_chunk_size)); \
		if(FOREACH_ID(fe_nElements) >= FOREACH_ID(fe_thread_count) || FOREACH_ID(fe_thread_id) == 0) { \
			for(size_t ITERVAR = FOREACH_ID(fe_index_begin); ITERVAR < FOREACH_ID(fe_index_end); ++ITERVAR) { \
				BODY \
			} \
		} \
	}


#define foreach_omp_iterator(VAR, BEGIN, END, OMPEXPR, BODY) \
	_Pragma(FOREACH_OMP_HELPER(OMPEXPR)) \
	{ \
		const auto FOREACH_ID(fe_container_begin) = BEGIN; \
		const auto FOREACH_ID(fe_container_end) = END; \
		const size_t FOREACH_ID(fe_container_size) = FOREACH_ID(fe_container_end) - FOREACH_ID(fe_container_begin); \
		const unsigned int FOREACH_ID(fe_thread_count) = omp_get_num_threads(); \
		const unsigned int FOREACH_ID(fe_thread_id) = omp_get_thread_num(); \
		const size_t FOREACH_ID(fe_chunk_size) = FOREACH_ID(fe_container_size) / FOREACH_ID(fe_thread_count); \
		const auto FOREACH_ID(fe_it_begin) = FOREACH_ID(fe_container_begin) + FOREACH_ID(fe_thread_id) * FOREACH_ID(fe_chunk_size); \
		const auto FOREACH_ID(fe_it_end) = ((FOREACH_ID(fe_thread_id) == FOREACH_ID(fe_thread_count) - 1 || FOREACH_ID(fe_container_size) < FOREACH_ID(fe_thread_count)) ? FOREACH_ID(fe_container_end) : FOREACH_ID(fe_it_begin) + FOREACH_ID(fe_chunk_size)); \
		if(FOREACH_ID(fe_container_size) >= FOREACH_ID(fe_thread_count) || FOREACH_ID(fe_thread_id) == 0) { \
			for(auto FOREACH_ID(fe_it) = FOREACH_ID(fe_it_begin); FOREACH_ID(fe_it) != FOREACH_ID(fe_it_end); ++FOREACH_ID(fe_it)) { \
				VAR = *FOREACH_ID(fe_it); \
				BODY \
			} \
		} \
	}


#define foreach_omp(VAR, CONTAINER, OMPEXPR, BODY) \
	_Pragma(FOREACH_OMP_HELPER(OMPEXPR)) \
	{ \
		auto const& FOREACH_ID(fe_container) = CONTAINER; \
		const auto FOREACH_ID(fe_container_cbegin) = FOREACH_ID(fe_container).cbegin(); \
		const auto FOREACH_ID(fe_container_cend) = FOREACH_ID(fe_container).cend(); \
		const size_t FOREACH_ID(fe_container_size) = FOREACH_ID(fe_container_cend) - FOREACH_ID(fe_container_cbegin); \
		const unsigned int FOREACH_ID(fe_thread_count) = omp_get_num_threads(); \
		const unsigned int FOREACH_ID(fe_thread_id) = omp_get_thread_num(); \
		const size_t FOREACH_ID(fe_chunk_size) = FOREACH_ID(fe_container_size) / FOREACH_ID(fe_thread_count); \
		const auto FOREACH_ID(fe_it_begin) = FOREACH_ID(fe_container_cbegin) + FOREACH_ID(fe_thread_id) * FOREACH_ID(fe_chunk_size); \
		const auto FOREACH_ID(fe_it_end) = ((FOREACH_ID(fe_thread_id) == FOREACH_ID(fe_thread_count) - 1 || FOREACH_ID(fe_container_size) < FOREACH_ID(fe_thread_count)) ? FOREACH_ID(fe_container_cend) : FOREACH_ID(fe_it_begin) + FOREACH_ID(fe_chunk_size)); \
		if(FOREACH_ID(fe_container_size) >= FOREACH_ID(fe_thread_count) || FOREACH_ID(fe_thread_id) == 0) { \
			for(auto FOREACH_ID(fe_it) = FOREACH_ID(fe_it_begin); FOREACH_ID(fe_it) != FOREACH_ID(fe_it_end); ++FOREACH_ID(fe_it)) { \
				VAR = *FOREACH_ID(fe_it); \
				BODY \
			} \
		} \
	}


#define foreachElementIndependently(VAR, MESH, OMPEXPR, BODY) \
	_Pragma(FOREACH_OMP_HELPER(OMPEXPR)) \
	{ \
		const unsigned int FOREACH_ID(fe_thread_count) = omp_get_num_threads(); \
		const unsigned int FOREACH_ID(fe_thread_id) = omp_get_thread_num(); \
		auto const& FOREACH_ID(fe_mesh) = MESH; \
		auto const& FOREACH_ID(fe_decomposition) = FOREACH_ID(fe_mesh).getIndependentElementDecomposition(); \
		const auto FOREACH_ID(fe_decomposition_cend) = FOREACH_ID(fe_decomposition).cend(); \
		for(auto FOREACH_ID(fe_decomp_it) = FOREACH_ID(fe_decomposition).cbegin(); FOREACH_ID(fe_decomp_it) != FOREACH_ID(fe_decomposition_cend); ++FOREACH_ID(fe_decomp_it)) { \
			auto const& FOREACH_ID(fe_container) = *FOREACH_ID(fe_decomp_it); \
			const auto FOREACH_ID(fe_container_cbegin) = FOREACH_ID(fe_container).cbegin(); \
			const auto FOREACH_ID(fe_container_cend) = FOREACH_ID(fe_container).cend(); \
			const size_t FOREACH_ID(fe_container_size) = FOREACH_ID(fe_container).size(); \
			const size_t FOREACH_ID(fe_chunk_size) = FOREACH_ID(fe_container_size) / FOREACH_ID(fe_thread_count); \
			const auto FOREACH_ID(fe_it_begin) = FOREACH_ID(fe_container_cbegin) + FOREACH_ID(fe_thread_id) * FOREACH_ID(fe_chunk_size); \
			const auto FOREACH_ID(fe_it_end) = ((FOREACH_ID(fe_thread_id) == FOREACH_ID(fe_thread_count) - 1 || FOREACH_ID(fe_container_size) < FOREACH_ID(fe_thread_count)) ? FOREACH_ID(fe_container_cend) : FOREACH_ID(fe_it_begin) + FOREACH_ID(fe_chunk_size)); \
			if(FOREACH_ID(fe_container_size) >= FOREACH_ID(fe_thread_count) || FOREACH_ID(fe_thread_id) == 0) { \
				for(auto FOREACH_ID(fe_it) = FOREACH_ID(fe_it_begin); FOREACH_ID(fe_it) != FOREACH_ID(fe_it_end); ++FOREACH_ID(fe_it)) { \
					VAR = FOREACH_ID(fe_mesh).element(*FOREACH_ID(fe_it)); \
					BODY \
				} \
			} \
			_Pragma("omp barrier") \
		} \
	}

#else // #ifdef SPRAT_BUILD_WITH_OPENMP

#define foreach_omp(VAR, CONTAINER, OMPEXPR, BODY) \
	foreach(VAR, CONTAINER, BODY)

#define foreach_omp_index(ITERVAR, NELEMENTS, OMPEXPR, BODY) \
	{ \
		const size_t FOREACH_ID(fe_nElements) = NELEMENTS; \
		for(size_t ITERVAR=0; ITERVAR<FOREACH_ID(fe_nElements); ++ITERVAR) { \
			BODY \
		} \
	}

#define foreach_omp_iterator(VAR, BEGIN, END, OMPEXPR, BODY) \
	{ \
		const auto FOREACH_ID(fe_it_begin) = BEGIN; \
		const auto FOREACH_ID(fe_it_end) = END; \
		for(auto FOREACH_ID(fe_it) = FOREACH_ID(fe_it_begin); FOREACH_ID(fe_it) != FOREACH_ID(fe_it_end); ++FOREACH_ID(fe_it)) { \
			VAR = *FOREACH_ID(fe_it); \
			BODY \
		} \
	}

#define foreachElementIndependently(VAR, MESH, OMPEXPR, BODY) \
	foreach(VAR, Elements(MESH), BODY)

#endif // #ifdef SPRAT_BUILD_WITH_OPENMP




#endif /* FOREACH_HPP_ */
