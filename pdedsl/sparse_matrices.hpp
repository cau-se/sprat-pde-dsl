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

#ifndef SPARSE_MATRICES_HPP_
#define SPARSE_MATRICES_HPP_

#include <vector>
#include <forward_list>
#include <utility>
#include <cfloat>
#include "config.hpp"
#include "util.hpp"
#include "numa.hpp"



class BareMatrix {
private:
	index_t _nRows;
	index_t _nCols;
	real * _entries; // row-major

	void freeMemory() {
		numa_delete(_entries);
	}

public:
	BareMatrix() : _nRows(0), _nCols(0), _entries(0) {}
	BareMatrix(index_t n) : _nRows(0), _nCols(0), _entries(0) {
		reset(n, n);
	}
	BareMatrix(index_t nRows, index_t nCols) : _nRows(0), _nCols(0), _entries(0) {
		reset(nRows, nCols);
	}

	~BareMatrix() {
		freeMemory();
	}

	void reset(index_t nRows, index_t nCols) {
		freeMemory();
		_nRows = nRows;
		_nCols = nCols;
		_entries = numa_new<real>(nEntries());
	}

	index_t nEntries() const {
		return _nRows * _nCols;
	}

	index_t nRows() const {
		return _nRows;
	}

	index_t nCols() const {
		return _nCols;
	}

	real get(index_t i, index_t j) const {
		return _entries[i*_nCols + j];
	}

	void set(index_t i, index_t j, real v) {
		_entries[i*_nCols + j] = v;
	}

	template <typename DestVecT, typename SrcVecT>
	DestVecT & multiplyAdd(DestVecT & dest, SrcVecT const& src) const {
		foreach_omp_index(i, _nRows, shared(dest, src), {
			real dotProduct = 0.0;
			const index_t rowOffset = i*_nCols;
			for(index_t j=0; j<_nCols; ++j) {
				dotProduct += _entries[rowOffset + j] * src[j];
			}
			dest[i] += dotProduct;
		})
		return dest;
	}
	template <typename DestVecT, typename SrcVecT>
	DestVecT & multiplySubtract(DestVecT & dest, SrcVecT const& src) const {
		foreach_omp_index(i, _nRows, shared(dest, src), {
			real dotProduct = 0.0;
			const index_t rowOffset = i*_nCols;
			for(index_t j=0; j<_nCols; ++j) {
				dotProduct += _entries[rowOffset + j] * src[j];
			}
			dest[i] -= dotProduct;
		})
		return dest;
	}
};

class LILMatrix {
private:
	std::vector< std::forward_list< std::pair<index_t, real> > > _rows;
	index_t _nEntries;

public:
	LILMatrix() : _rows(0), _nEntries(0) {}
	LILMatrix(index_t nRows) : _rows(nRows), _nEntries(0) {}

	void reset(index_t nRows) {
		_rows.clear();
		_rows.resize(nRows);
		_nEntries = 0;
	}

	void addElement(index_t row, index_t col, real value) {
		auto& list = _rows[row];

		if(list.empty()) {
			list.push_front(std::pair<index_t, real>(col, value));
			++_nEntries;
		}
		else {
			for(auto& e : list) {
				if(e.first > col) {
					break;
				}
				if(e.first == col) {
					e.second += value;
					return;
				}
			}

			// Column does not exist; find appropriate insertion point.
			if(list.front().first > col) {
				list.push_front(std::pair<index_t, real>(col, value));
				++_nEntries;
			}
			else {
				auto itPrev = list.begin();
				auto it = itPrev;
				while(it != list.end()) {
					if((*it).first > col) {
						list.insert_after(itPrev, std::pair<index_t, real>(col, value));
						++_nEntries;
						return;
					}
					itPrev = it;
					++it;
				}
				list.insert_after(itPrev, std::pair<index_t, real>(col, value));
				++_nEntries;
			}
		}
	}

	std::forward_list< std::pair<index_t, real> > const& getRow(index_t row) const {
			return _rows[row];
	}
	std::forward_list< std::pair<index_t, real> >& getRow(index_t row) {
		return _rows[row];
	}

	index_t nEntries() const {
		return _nEntries;
	}
	index_t nRows() const {
		return _rows.size();
	}

	void print() const {
		for(index_t i=0; i<_rows.size(); ++i) {
			std::cout << "Row " << i << ": ";
			for(auto const& e : _rows[i]) {
				std::cout << "(" << e.first << ", " << e.second << ") ";
			}
			std::cout << std::endl;
		}
	}
};


class BareCSRMatrix {
private:
	index_t _nRows;
	index_t _nNonZeros;

	real * _entries;
	index_t * _columns;

	index_t * _rowStart; // n+1 elements (the begin of the (n+1)-th row is also stored)

	void freeMemory() {
		numa_delete(_rowStart);
		numa_delete(_columns);
		numa_delete(_entries);
	}

public:
	BareCSRMatrix() :
		_nRows(0),
		_nNonZeros(0),
		_entries(0),
		_columns(0),
		_rowStart(0)
	{}

	BareCSRMatrix(LILMatrix const& lilMatrix) :
		_nRows(0),
		_nNonZeros(0),
		_entries(0),
		_columns(0),
		_rowStart(0)
	{
		reset(lilMatrix);
	}

	~BareCSRMatrix() {
		freeMemory();
	}

	void reset(LILMatrix const& lilMatrix) {
		freeMemory();

		_nRows = lilMatrix.nRows();
		_nNonZeros = lilMatrix.nEntries();

		_entries = numa_new<real>(_nNonZeros);
		_columns = numa_new<index_t>(_nNonZeros);
		_rowStart = numa_new<index_t>(_nRows+1);

		index_t element = 0;
		for(index_t row=0; row<_nRows; ++row) {
			_rowStart[row] = element;

			auto const& list = lilMatrix.getRow(row);
			for(auto const& e : list) {
				_entries[element] = e.second;
				_columns[element] = e.first;
				++element;
			}
		}
		_rowStart[_nRows] = element;
	}



	index_t nEntries() const {
		return _nNonZeros;
	}
	index_t nRows() const {
		return _nRows;
	}

	template <typename DestVecT, typename SrcVecT>
	DestVecT & multiplyAdd(DestVecT & dest, SrcVecT const& src) const {
		foreach_omp_index(row, _nRows, shared(dest, src), {
			real rowResult = 0.0;
			index_t endIndex = _rowStart[row+1];
			for(index_t itemIndex=_rowStart[row]; itemIndex<endIndex; ++itemIndex) {
				rowResult += _entries[itemIndex] * src[_columns[itemIndex]];
			}
			dest[row] += rowResult;
		})
		return dest;
	}
	template <typename DestVecT, typename SrcVecT>
	DestVecT & multiplySubtract(DestVecT & dest, SrcVecT const& src) const {
		foreach_omp_index(row, _nRows, shared(dest, src), {
			real rowResult = 0.0;
			index_t endIndex = _rowStart[row+1];
			for(index_t itemIndex=_rowStart[row]; itemIndex<endIndex; ++itemIndex) {
				rowResult += _entries[itemIndex] * src[_columns[itemIndex]];
			}
			dest[row] -= rowResult;
		})
		return dest;
	}

	std::pair<index_t, real> findMaximum(index_t row, real * vector, index_t offset, index_t stride) const {
		real max_value = -1.0*DBL_MAX;
		index_t max_index = 0;

		for(index_t index=_rowStart[row]; index<_rowStart[row+1]; ++index) {
			const index_t column = _columns[index];
			const real value = _entries[index] * vector[stride*column + offset];
			if(value > max_value) {
				max_value = value;
				max_index = column;
			}
		}

		return std::pair<index_t, real>(max_index, max_value);
	}

	std::pair<index_t, real> findMinimum(index_t row, real * vector, index_t offset, index_t stride) const {
		real min_value = DBL_MAX;
		index_t min_index = 0;

		for(index_t index=_rowStart[row]; index<_rowStart[row+1]; ++index) {
			const index_t column = _columns[index];
			const real value = _entries[index] * vector[stride*column + offset];
			if(value < min_value) {
				min_value = value;
				min_index = column;
			}
		}

		return std::pair<index_t, real>(min_index, min_value);
	}

	void print() const {
		for(index_t i=0; i<_nRows; ++i) {
			std::cout << "Row " << i << " from " << _rowStart[i] << " to " << (_rowStart[i+1]-1) << ": ";
			for(index_t j=_rowStart[i]; j<_rowStart[i+1]; ++j) {
				std::cout << "(" << _columns[j] << ", " << _entries[j] << ") ";
			}
			std::cout << std::endl;
		}
	}
};


#endif /* SPARSE_MATRICES_HPP_ */
