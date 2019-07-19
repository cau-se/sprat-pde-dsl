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

#ifndef PVECTOR_HPP_
#define PVECTOR_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector>
#ifdef SPRAT_BUILD_WITH_MPI
#include <mpi.h>
#endif
#include "config.hpp"
#include "util.hpp"
#include "parallel.hpp"
#include "numa.hpp"




class BareDistributedVectorBase {
protected:
	real * _values;
	index_t _nLocalValues;
	index_t _nValues;

	ParallelExecutionEnvironment const* _pEE;
	MPCommunicationRegistry const* _comRegistry;
#ifdef SPRAT_BUILD_WITH_MPI
	std::vector<MPI_Request> _sendRequests;
	std::vector<MPI_Request> _recvRequests;
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	int _tag;

	//BareDistributedVectorBase(BareDistributedVectorBase const& other) :
	//	_values(other._values),
	//	_nValues(other._nValues) {}
	BareDistributedVectorBase(real * values, int tag, ParallelExecutionEnvironment const* pEE, MPCommunicationRegistry const* comRegistry) :
		_values(values),
		_nLocalValues(0),
		_nValues(0),
		_pEE(pEE),
		_comRegistry(comRegistry),
#ifdef SPRAT_BUILD_WITH_MPI
		_sendRequests(0),
		_recvRequests(0),
#endif //#ifdef SPRAT_BUILD_WITH_MPI
		//_singleRequest(),
		_tag(tag) {
#ifdef SPRAT_BUILD_WITH_MPI
		if(pEE) {
			_sendRequests.reserve(2*pEE->nComputeProcesses());
			_recvRequests.reserve(2*pEE->nComputeProcesses());
		}
#endif //#ifdef SPRAT_BUILD_WITH_MPI
		if(comRegistry) {
			_nLocalValues = comRegistry->nLocalDoF;
			_nValues = comRegistry->nDoF;
		}
	}

public:

	void touch(real initValue = 0.0) {
		foreach_omp_iterator(real& v, _values, _values+_nValues, shared(initValue), {
			v = initValue;
		})
	}
	void touch(std::function<void (real *, index_t)> touchFunction) {
		if(touchFunction) {
			touchFunction(_values, _nValues);
		}
	}

	void touchLocal(real initValue = 0.0) {
		foreach_omp_iterator(real& v, _values, _values+_nLocalValues, shared(initValue), {
			v = initValue;
		})
	}
	void touchLocal(std::function<void (real *, index_t)> touchFunction) {
		if(touchFunction) {
			touchFunction(_values, _nLocalValues);
		}
	}

	void printLocal() const {
		for(index_t i=0; i<_nLocalValues; ++i) {
			std::cout << std::setprecision(6) << std::scientific << _values[i] << std::endl;
		}
		std::cout << std::endl;
	}
	void print() const {
		for(index_t i=0; i<_nValues; ++i) {
			std::cout << std::setprecision(6) << std::scientific << _values[i] << std::endl;
		}
		std::cout << std::endl;
	}

	index_t size() const {
		return _nValues;
	}
	index_t localSize() const {
		return _nLocalValues;
	}

	real * getValues() const {
		return _values;
	}

	real& operator[](index_t index) const {
		return _values[index];
	}
	real get(index_t index) const {
		return _values[index];
	}
	real getLocal(index_t index) const {
		return _values[index];
	}


#ifdef SPRAT_BUILD_WITH_MPI
	void collectFromMasterAsync() {
		_recvRequests.push_back(_pEE->recvAsync(_values, _nLocalValues, 0, _tag));
	}
	void distributeToMasterAsync() {
		_sendRequests.push_back(_pEE->sendAsync(_values, _nLocalValues, 0, _tag));
	}
#endif //#ifdef SPRAT_BUILD_WITH_MPI

	void finishPendingSends() {
#ifdef SPRAT_BUILD_WITH_MPI
		_pEE->waitFor(_sendRequests);
		_sendRequests.clear();
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}
	void finishPendingRecvs() {
#ifdef SPRAT_BUILD_WITH_MPI
		_pEE->waitFor(_recvRequests);
		_recvRequests.clear();
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}
	void finishPendingCommunication() {
		finishPendingSends();
		finishPendingRecvs();
	}


	real norm() const {
		real localResult = 0.0;
		foreach_omp_index(i, _nLocalValues, reduction(+:localResult),  {
			localResult += _values[i] * _values[i];
		})

		return sqrt(_pEE->globalSum(localResult));
	}

	template <typename VecT>
	real dotProduct(VecT const& other) const {
		real localResult = 0.0;
		foreach_omp_index(i, _nLocalValues, shared(other) reduction(+:localResult), {
			localResult += _values[i] * other.get(i);
		})

		return _pEE->globalSum(localResult);
	}

	real average() const {
		real result = 0.0;
		for(index_t i=0; i<_nLocalValues; ++i) {
			result += _values[i];
		}
		return result / (real)_nLocalValues;
	}

	void setTag(int tag) {
		_tag=tag;
	}
	int getTag() const {
		return _tag;
	}
	//int getNeighborhoodTag(uint i) const {
	//	return _baseTag + 1 + _neighborhoods[i].neighborhoodID;
	//}

	void exchangeDataAsync() {
#ifdef SPRAT_BUILD_WITH_MPI
		finishPendingCommunication();

		for(uint i=0; i<_comRegistry->dofRecvFromNeighbor.size(); ++i) {
			_recvRequests.push_back(_pEE->recvAsyncMPIType(_values,
					_comRegistry->dofRecvFromNeighborMPIType[i],
					_comRegistry->dofRecvFromNeighbor[i].neighborID,
					getTag()));
		}
		for(uint i=0; i<_comRegistry->dofSendToNeighbor.size(); ++i) {
			_sendRequests.push_back(_pEE->sendAsyncMPIType(_values,
					_comRegistry->dofSendToNeighborMPIType[i],
					_comRegistry->dofSendToNeighbor[i].neighborID,
					getTag()));
		}
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}

	void exchangeDataBlocked() {
		exchangeDataAsync();
		finishPendingCommunication();
	}


	void setTo(real value) const {
		foreach_omp_index(i, _nValues, shared(value), {
			_values[i] = value;
		})
	}
	template <typename VecT>
	void setTo(VecT const& other) const {
		foreach_omp_index(i, _nLocalValues, shared(other), {
			_values[i] = other.getLocal(i);
		})
	}

//	using Vector::setTo;
//	void setTo(ParallelVector const& other) const {
//		for(uint i=0; i<_nValues; ++i) {
//			_values[i] = other.get(i);
//		}
//	}
//	void addScaled(real lambda, ParallelVector const& other) {
//		for(uint i=0; i<_nValues; ++i) {
//			_values[i] += lambda*other.get(i);
//		}
//	}
};


class BareDistributedVector : public BareDistributedVectorBase {
private:
	friend class BareDistributedVectorView;

public:
	BareDistributedVector() : BareDistributedVectorBase(0, 0, 0, 0) {}
	BareDistributedVector(std::string const& name, ParallelExecutionEnvironment& pEE, MPCommunicationRegistry& comRegistry) :
		BareDistributedVectorBase(numa_new<real>(comRegistry.nDoF),
				pEE.getTagFromName(name),
				&pEE,
				&comRegistry) {}
	~BareDistributedVector() {
		numa_delete(_values);
	}

	void reset(std::string name, ParallelExecutionEnvironment& pEE, MPCommunicationRegistry& comRegistry) {
		numa_delete(_values);
		_values = numa_new<real>(comRegistry.nDoF);
		_nValues = comRegistry.nDoF;
		_nLocalValues = comRegistry.nLocalDoF;
		_tag = pEE.getTagFromName(name);
		_pEE = &pEE;
		_comRegistry = &comRegistry;

#ifdef SPRAT_BUILD_WITH_MPI
		_sendRequests.resize(0);
		_sendRequests.reserve(2*pEE.nComputeProcesses());
		_recvRequests.resize(0);
		_recvRequests.reserve(2*pEE.nComputeProcesses());
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}

	template <typename VecT>
	BareDistributedVector const& operator=(VecT const& other) {
		foreach_omp_index(i, _nLocalValues, shared(other), {
			_values[i] = other.getLocal(i);
		})
		return *this;
	}
};

class BareDistributedVectorView : public BareDistributedVectorBase {
private:

public:
	BareDistributedVectorView() : BareDistributedVectorBase(0, 0, 0, 0) {}
	BareDistributedVectorView(BareDistributedVector const& other) :
		BareDistributedVectorBase(other.getValues(),
				other.getTag(),
				other._pEE,
				other._comRegistry) {}

	void reset(BareDistributedVector const& other) {
		_values = other.getValues();
		_nValues = other.size();
		_nLocalValues = other.localSize();
		_tag = other.getTag();
		_pEE = other._pEE;
		_comRegistry = other._comRegistry;

#ifdef SPRAT_BUILD_WITH_MPI
		_sendRequests.resize(0);
		_recvRequests.resize(0);

		if(_pEE) {
			_sendRequests.reserve(2*_pEE->nComputeProcesses());
			_recvRequests.reserve(2*_pEE->nComputeProcesses());
		}
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}

	template <typename VecT>
	BareDistributedVectorView const& operator=(VecT const& other) {
		foreach_omp_index(i, _nLocalValues, shared(other), {
			_values[i] = other.getLocal(i);
		})
		return *this;
	}
};




class BareVectorBase {
private:
	BareVectorBase() {}
protected:
	real * _values;
	index_t _nValues;

	BareVectorBase(BareVectorBase const& other) : _values(other._values), _nValues(other._nValues) {}
	BareVectorBase(real nValues, real * values) : _values(values), _nValues(nValues) {}

public:

	void touch(real initValue = 0.0) {
		foreach_omp_iterator(real& v, _values, _values+_nValues, shared(initValue), {
			v = initValue;
		})
	}
	void touch(std::function<void (real *, index_t)> touchFunction) {
		if(touchFunction) {
			touchFunction(_values, _nValues);
		}
	}

	void printToFile(std::string prefix, uint id, std::string suffix = ".txt") const {
		std::ofstream myfile;
		const int outputPrecision = 16;
		myfile.open((prefix + std::to_string(id) + suffix).c_str());
		for(index_t i=0; i<_nValues; ++i) {
			myfile << std::setprecision(outputPrecision) << std::scientific << _values[i] << std::endl;
		}
		myfile.close();
	}

	void writeBinaryToFile(std::string prefix, uint id, std::string suffix = ".bin") const {
		std::ofstream myfile;
		myfile.open((prefix + std::to_string(id) + suffix).c_str(), std::ios_base::out | std::ios_base::binary);
		for(index_t i=0; i<_nValues; ++i) {
			char * const location = (char *)&_values[i];
			myfile.write(location, sizeof(real));
		}
		myfile.close();
	}

	void readBinaryFromFile(std::string prefix, uint id, std::string suffix = ".bin") {
		std::ifstream myfile;
		myfile.open((prefix + std::to_string(id) + suffix).c_str(), std::ios_base::in | std::ios_base::binary);
		for(index_t i=0; i<_nValues; ++i) {
			char * const location = (char *)&_values[i];
			myfile.read(location, sizeof(real));
		}
		myfile.close();
	}
	void print() const {
		for(index_t i=0; i<_nValues; ++i) {
			std::cout << std::setprecision(6) << std::scientific << _values[i] << std::endl;
		}
		std::cout << std::endl;
	}

	index_t size() const {
		return _nValues;
	}
	index_t localSize() const {
		return _nValues;
	}

	real * getValues() const {
		return _values;
	}
	real * getLocalValues() const {
		return _values;
	}

	real& operator[](index_t index) const {
		return _values[index];
	}
	real get(index_t index) const {
		return _values[index];
	}
	real getLocal(index_t index) const {
		return _values[index];
	}

#ifdef SPRAT_BUILD_WITH_MPI
	std::vector<MPI_Request> collectFromSlavesAsync(std::string const& name, ParallelExecutionEnvironment const& pEE, std::vector<index_t> const& dofStartIndexForProcess) const {
		int tag = pEE.getTagFromName(name);
		std::vector<MPI_Request> requests(pEE.nComputeProcesses());

		for(int i=0; i<pEE.nComputeProcesses(); ++i) {
			requests[i] = pEE.recvAsync(_values+dofStartIndexForProcess[i], dofStartIndexForProcess[i+1]-dofStartIndexForProcess[i], i+1, tag);
		}

		return requests;
	}
	std::vector<MPI_Request> distributeToSlavesAsync(std::string const& name, ParallelExecutionEnvironment const& pEE, std::vector<index_t> const& dofStartIndexForProcess) const {
		int tag = pEE.getTagFromName(name);
		std::vector<MPI_Request> requests(pEE.nComputeProcesses());

		for(int i=0; i<pEE.nComputeProcesses(); ++i) {
			requests[i] = pEE.sendAsync(_values+dofStartIndexForProcess[i], dofStartIndexForProcess[i+1]-dofStartIndexForProcess[i], i+1, tag);
		}

		return requests;
	}
#endif // SPRAT_BUILD_WITH_MPI

	real dotProduct(BareDistributedVector const& other) const {
		return other.dotProduct(*this);
	}
	real dotProduct(BareDistributedVectorView const& other) const {
		return other.dotProduct(*this);
	}

	template <typename VecT>
	real dotProduct(VecT const& other) const {
		real result = 0.0;
		foreach_omp_index(i, _nValues, shared(other) reduction(+:result), {
			result += _values[i] * other.getLocal(i);
		})
		return result;
	}
	template <typename VecT>
	real dotProductPartial(VecT const& other, index_t nEntries) const {
		real result = 0.0;
		foreach_omp_index(i, nEntries, shared(other,nEntries) reduction(+:result), {
			result += _values[i] * other.getLocal(i);
		})
		return result;
	}


	real norm() const {
		real result = 0.0;
		foreach_omp_index(i, _nValues, reduction(+:result), {
			result += _values[i] * _values[i];
		})
		return sqrt(result);
	}

	real componentSum() const {
		real result = 0.0;
		foreach_omp_index(i, _nValues, reduction(+:result), {
			result += _values[i];
		})
		return result;
	}

	real average() const {
		real result = 0.0;
		for(index_t i=0; i<_nValues; ++i) {
			result += _values[i];
		}
		return result / (real)_nValues;
	}



	void setTo(real value) const {
		foreach_omp_index(i, _nValues, shared(value), {
			_values[i] = value;
		})
	}
	template <typename VecT>
	void setTo(VecT const& other) const {
		foreach_omp_index(i, _nValues, shared(other), {
			_values[i] = other.getLocal(i);
		})
	}

//	void scale(real lambda) const {
//		for(uint i=0; i<_nValues; ++i) {
//			_values[i] *= lambda;
//		}
//	}
//
//	template <typename VecT>
//	void addScaled(real lambda, VecT const& other) const {
//		for(uint i=0; i<_nValues; ++i) {
//			_values[i] += lambda*other.getLocal(i);
//		}
//	}
//
//	template <typename VecT>
//	void divideBy(VecT const& other) const {
//		for(uint i=0; i<_nValues; ++i) {
//			_values[i] /= other.getLocal(i);
//		}
//	}
};

class BareVector : public BareVectorBase {
public:
	BareVector() : BareVectorBase(0, 0) {}
	explicit BareVector(index_t nValues) : BareVectorBase(nValues, numa_new<real>(nValues)) {}
	explicit BareVector(BareVector const& other) : BareVectorBase(other.size(), numa_new<real>(other.size())) {
		foreach_omp_index(i, other.size(), , {
			_values[i] = other._values[i];
		})
	}

	~BareVector() {
		//std::cout << "~BareVector: Freeing " << _values << std::endl;
		numa_delete(_values);
	}

	void reset(index_t nValues) {
		numa_delete(_values);
		_values = numa_new<real>(nValues);
		//std::cout << "BareVector::reset: Allocated " << _values << std::endl;
		_nValues = nValues;
	}

	template <typename VecT>
	BareVector const& operator=(VecT const& other) {
		foreach_omp_index(i, _nValues, shared(other), {
			_values[i] = other.getLocal(i);
		})
		return *this;
	}
};

class BareVectorView : public BareVectorBase {
public:
	BareVectorView() : BareVectorBase(0, 0) {}
	BareVectorView(index_t nValues, real * values) : BareVectorBase(nValues, values) {}
	explicit BareVectorView(BareVector const& other) : BareVectorBase(other.size(), other.getLocalValues()) {}
	explicit BareVectorView(BareVectorView const& other) : BareVectorBase(other.size(), other.getLocalValues()) {}

	void reset(index_t nValues, real * values) {
		_values = values;
		_nValues = nValues;
	}

	void reset(BareVector const& other) {
		_values = other.getLocalValues();
		_nValues = other.size();
	}

	template <typename VecT>
	BareVectorView const& operator=(VecT const& other) {
		foreach_omp_index(i, _nValues, shared(other), {
			_values[i] = other.getLocal(i);
		})
		return *this;
	}
};



#endif //PVECTOR_HPP_
