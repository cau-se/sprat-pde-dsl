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

#ifndef NUMA_HPP_
#define NUMA_HPP_


#include <iostream>
#include <cstdlib>
#include <cstddef>
#include <functional>
#include "foreach.hpp"

template <typename T>
inline T* numa_new(size_t nInstances) {
#ifndef SPRAT_NO_NUMA
	T* ptr = static_cast<T*>(malloc(nInstances*sizeof(T)));
#else
	T* ptr = new T[nInstances];
#endif
	if(ptr == 0) {
		std::cout << "numa_new(): malloc failed!" << std::endl;
		exit(1);
	}
	return ptr;
}

template <typename T>
inline void numa_delete(T* ptr) {
#ifndef SPRAT_NO_NUMA
	if(ptr != 0) {
		free(ptr);
	}
#else
	delete[] ptr;
#endif
}

template <typename T>
inline T* numa_new_omp_init(size_t nInstances, T const& initValue) {
	T* ptr = numa_new<T>(nInstances);
	foreach_omp_iterator(T& v, ptr, ptr+nInstances, shared(ptr), {
		v = initValue;
	})
	return ptr;
}

template <typename T>
class numa_vector {
private:
	T* values;
	size_t nValues;

public:
	typedef T                 value_type;
	typedef value_type&       reference;
	typedef value_type const& const_reference;
	typedef value_type*       pointer;
	typedef value_type const* const_pointer;
	typedef pointer           iterator;
	typedef const_pointer     const_iterator;
	typedef ptrdiff_t         difference_type;
	typedef size_t            size_type;

	numa_vector() : values(0), nValues(0) {}
	numa_vector(size_type n) :
		values(numa_new<T>(n)),
		nValues(n)
	{}
	numa_vector(size_type n, const_reference initValue) :
		values(numa_new_omp_init<T>(n, initValue)),
		nValues(n)
	{}
	numa_vector(numa_vector<T> const& other) :
		values(numa_new<T>(other.size())),
		nValues(other.size())
	{
		foreach_omp_iterator(reference v, values, values+nValues, , {
				v = other[&v-values];
		})
	}

	~numa_vector() {
		numa_delete<T>(values);
	}

	void touch(const_reference initValue) {
		foreach_omp_iterator(reference v, values, values+nValues, shared(initValue), {
			v = initValue;
		})
	}
	void touch(std::function<void (pointer, size_type)> touchFunction) {
		if(touchFunction) {
			touchFunction(values, nValues);
		}
	}




	const_reference operator[](size_t i) const {
		return values[i];
	}
	reference operator[](size_t i) {
		return values[i];
	}

	numa_vector<T>& operator=(numa_vector<T> const& other) {
		if(nValues != other.size()) {
			resize(other.size());
		}
		foreach_omp_iterator(reference v, values, values+nValues, , {
			v = other[&v-values];
		})

		return *this;
	}

	size_type size() const {
		return nValues;
	}

	void resize(size_type n) {
		numa_delete<T>(values);
		values = numa_new<T>(n);
		nValues = n;
	}

	void resize(size_type n, const_reference initValue) {
		numa_delete<T>(values);
		values = numa_new_omp_init<T>(n, initValue);
		nValues = n;
	}

	bool empty() const {
		return (nValues == 0);
	}

	void clear() {
		numa_delete<T>(values);
		nValues = 0;
	}



	iterator begin() noexcept {
		return values;
	}
	const_iterator begin() const noexcept {
		return values;
	}
	iterator end() noexcept {
		return values + nValues;
	}
	const_iterator end() const noexcept {
		return values + nValues;
	}

	const_iterator cbegin() const noexcept {
		return values;
	}
	const_iterator cend() const noexcept {
		return values + nValues;
	}
};




#endif /* NUMA_HPP_ */
