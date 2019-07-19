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

#ifndef PARALLEL_HPP_
#define PARALLEL_HPP_



#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#ifdef SPRAT_BUILD_WITH_MPI
#include <mpi.h>
#endif //#ifdef SPRAT_BUILD_WITH_MPI
#include "config.hpp"
#include "util.hpp"

#ifdef SPRAT_BUILD_WITH_MPI
#define MPI_CALL(OPERATION, ERRORVAR) \
	if((ERRORVAR = OPERATION) != MPI_SUCCESS) { \
		ParallelExecutionEnvironment::handleError(ERRORVAR); \
	}
#endif //#ifdef SPRAT_BUILD_WITH_MPI



class ParallelExecutionEnvironment {
#ifdef SPRAT_BUILD_WITH_MPI
private:
	int _nProcesses;
	int _ID;

	int _tagStride;
	int _currentTag;

	std::vector<std::string> _tagDictionary;

	MPI_Comm _workerComm;

	bool _isSetUp;

public:
	ParallelExecutionEnvironment() :
		_nProcesses(0),
		_ID(0),
		_tagStride(1),
		_currentTag(0),
		_tagDictionary(0),
		_isSetUp(false) {}
	~ParallelExecutionEnvironment() {
		if(_isSetUp) {
			finish();
		}
	}

	static void handleError(int error) {
		char errorMessage[MPI_MAX_ERROR_STRING];
		int messageLength;

		MPI_Error_string(error, errorMessage, &messageLength);
		std::cout << "MPI error: " << errorMessage << std::endl;

		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	void init(int * argc, char *** argv) {
		int mpiErrorID;
		MPI_CALL(MPI_Init(argc, argv), mpiErrorID);

		MPI_CALL(MPI_Comm_size(MPI_COMM_WORLD, &_nProcesses), mpiErrorID);
		MPI_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &_ID), mpiErrorID);

		MPI_Comm oldWorldComm;
		MPI_CALL(MPI_Comm_dup(MPI_COMM_WORLD, &oldWorldComm), mpiErrorID);
		MPI_CALL(MPI_Comm_split(oldWorldComm, (isMaster() ? MPI_UNDEFINED : 1), _ID, &_workerComm), mpiErrorID);

		if(_nProcesses<=1) {
			std::cout << "ERROR: Cannot run IO-Master without compute slaves!" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		_isSetUp = true;
	}

	void finish() {
		int mpiErrorID;
		/*if(!isMaster()) {
			MPI_CALL(MPI_Comm_free(&_workerComm), mpiErrorID);
		}*/
		MPI_CALL(MPI_Finalize(), mpiErrorID);
		_isSetUp = false;
	}

	bool isMaster() const {
		return (_ID == 0);
	}
	bool isSetUp() const { return _isSetUp; }
	int processID() const { return _ID; }
	int computeProcessID() const { return _ID-1; }
	int nProcesses() const { return _nProcesses; }
	int nComputeProcesses() const { return _nProcesses-1; }


	int getTagFromName(std::string const& name) const {
		auto result = std::find(_tagDictionary.begin(), _tagDictionary.end(), name);
		if(result == _tagDictionary.end()) {
			std::cout << "WARNING: '" << name << "' is not a registered name!" << std::endl;
		}
		return result-_tagDictionary.begin();
	}
	void registerName(std::string const& name) {
		if(std::find(_tagDictionary.begin(), _tagDictionary.end(), name) == _tagDictionary.end()) {
			_tagDictionary.push_back(name);
		}
		else {
			std::cout << "WARNING: '" << name << "' has already be registered as a name!" << std::endl;
		}
	}




	// ********************** MPI-Abstraktionen

	void send(real const * data, int nItems, int destID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Send(const_cast<real *>(data), nItems, SPRAT_MPI_REAL, destID, tag, MPI_COMM_WORLD), mpiErrorID);
	}
	void send(index_t const * data, int nItems, int destID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Send(const_cast<index_t *>(data), nItems, SPRAT_MPI_INDEX_T, destID, tag, MPI_COMM_WORLD), mpiErrorID);
	}
	void sendShortIndices(short_index_t const * data, int nItems, int destID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Send(const_cast<short_index_t *>(data), nItems, SPRAT_MPI_SHORT_INDEX_T, destID, tag, MPI_COMM_WORLD), mpiErrorID);
	}
	MPI_Request sendAsync(real const * data, int nItems, int destID=0, int tag=0) const {
		MPI_Request request;
		int mpiErrorID;
		MPI_CALL(MPI_Isend(const_cast<real *>(data), nItems, SPRAT_MPI_REAL, destID, tag, MPI_COMM_WORLD, &request), mpiErrorID);

		return request;
	}
	MPI_Request sendAsyncMPIType(real const * data, MPI_Datatype type, int destID=0, int tag=0) const {
		MPI_Request request;
		int mpiErrorID;
		MPI_CALL(MPI_Isend(const_cast<real *>(data), 1, type, destID, tag, MPI_COMM_WORLD, &request), mpiErrorID);

		return request;
	}
	void sendReal(real datum, int destID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Send(&datum, 1, SPRAT_MPI_REAL, destID, tag, MPI_COMM_WORLD), mpiErrorID);
	}
	void sendUInt(uint datum, int destID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Send(&datum, 1, MPI_UNSIGNED, destID, tag, MPI_COMM_WORLD), mpiErrorID);
	}
	void sendInt(int datum, int destID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Send(&datum, 1, MPI_INT, destID, tag, MPI_COMM_WORLD), mpiErrorID);
	}
	void sendShortIndex(short_index_t datum, int destID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Send(&datum, 1, SPRAT_MPI_SHORT_INDEX_T, destID, tag, MPI_COMM_WORLD), mpiErrorID);
	}
	void sendIndex(index_t datum, int destID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Send(&datum, 1, SPRAT_MPI_INDEX_T, destID, tag, MPI_COMM_WORLD), mpiErrorID);
	}
	void sendVectorContent(std::vector<real> const& vec, int destID=0, int tag=0) const {
		send(&vec[0], vec.size(), destID, tag);
	}
	void sendVectorContent(std::vector<index_t> const& vec, int destID=0, int tag=0) const {
		send(&vec[0], vec.size(), destID, tag);
	}
	void sendShortIndexVectorContent(std::vector<short_index_t> const& vec, int destID=0, int tag=0) const {
		sendShortIndices(&vec[0], vec.size(), destID, tag);
	}

	void recv(real * data, int nItems, int sourceID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Recv(data, nItems, SPRAT_MPI_REAL, sourceID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE), mpiErrorID);
	}
	void recv(index_t * data, int nItems, int sourceID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Recv(data, nItems, SPRAT_MPI_INDEX_T, sourceID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE), mpiErrorID);
	}
	void recvShortIndices(short_index_t * data, int nItems, int sourceID=0, int tag=0) const {
		int mpiErrorID;
		MPI_CALL(MPI_Recv(data, nItems, SPRAT_MPI_SHORT_INDEX_T, sourceID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE), mpiErrorID);
	}
	MPI_Request recvAsync(real * data, int nItems, int sourceID=0, int tag=0) const {
		MPI_Request request;
		int mpiErrorID;
		MPI_CALL(MPI_Irecv(data, nItems, SPRAT_MPI_REAL, sourceID, tag, MPI_COMM_WORLD, &request), mpiErrorID);

		return request;
	}
	MPI_Request recvAsyncMPIType(real * data, MPI_Datatype type, int sourceID=0, int tag=0) const {
		MPI_Request request;
		int mpiErrorID;
		MPI_CALL(MPI_Irecv(data, 1, type, sourceID, tag, MPI_COMM_WORLD, &request), mpiErrorID);

		return request;
	}
	real recvReal(int sourceID=0, int tag=0) const {
		real result;
		int mpiErrorID;
		MPI_CALL(MPI_Recv(&result, 1, SPRAT_MPI_REAL, sourceID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE), mpiErrorID);
		return result;
	}
	uint recvUInt(int sourceID=0, int tag=0) const {
		uint result;
		int mpiErrorID;
		MPI_CALL(MPI_Recv(&result, 1, MPI_UNSIGNED, sourceID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE), mpiErrorID);
		return result;
	}
	int recvInt(int sourceID=0, int tag=0) const {
		int result;
		int mpiErrorID;
		MPI_CALL(MPI_Recv(&result, 1, MPI_INT, sourceID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE), mpiErrorID);
		return result;
	}
	short_index_t recvShortIndex(int sourceID=0, int tag=0) const {
		short_index_t result;
		int mpiErrorID;
		MPI_CALL(MPI_Recv(&result, 1, SPRAT_MPI_SHORT_INDEX_T, sourceID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE), mpiErrorID);
		return result;
	}
	index_t recvIndex(int sourceID=0, int tag=0) const {
		index_t result;
		int mpiErrorID;
		MPI_CALL(MPI_Recv(&result, 1, SPRAT_MPI_INDEX_T, sourceID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE), mpiErrorID);
		return result;
	}
	void recvVectorContent(std::vector<real>& vec, int sourceID=0, int tag=0) const {
		recv(&vec[0], vec.size(), sourceID, tag);
	}
	void recvVectorContent(std::vector<index_t>& vec, int sourceID=0, int tag=0) const {
		recv(&vec[0], vec.size(), sourceID, tag);
	}
	void recvShortIndexVectorContent(std::vector<short_index_t>& vec, int sourceID=0, int tag=0) const {
		recvShortIndices(&vec[0], vec.size(), sourceID, tag);
	}

	void waitFor(MPI_Request request) const {
		int mpiErrorID;
		MPI_CALL(MPI_Wait(&request, MPI_STATUS_IGNORE), mpiErrorID);
	}
	void waitFor(std::vector<MPI_Request>& requests) const {
		int mpiErrorID;
		//MPI_CALL(MPI_Waitall(requests.size(), &requests[0], MPI_STATUS_IGNORE), mpiErrorID);
		for(uint i=0; i<requests.size(); ++i) {
			MPI_CALL(MPI_Wait(&requests[i], MPI_STATUS_IGNORE), mpiErrorID);
		}
	}

	void barrier() const {
		int mpiErrorID;
		MPI_CALL(MPI_Barrier(MPI_COMM_WORLD), mpiErrorID);
	}
	void computeNodeBarrier() const {
		int mpiErrorID;
		MPI_CALL(MPI_Barrier(_workerComm), mpiErrorID);
	}


	real globalSum(real localSummand) const {
		int mpiErrorID;
		real globalResult;
		MPI_CALL(MPI_Allreduce(&localSummand, &globalResult, 1, SPRAT_MPI_REAL, MPI_SUM, _workerComm), mpiErrorID);
		return globalResult;
	}
	index_t globalSum(index_t localSummand) const {
		int mpiErrorID;
		index_t globalResult;
		MPI_CALL(MPI_Allreduce(&localSummand, &globalResult, 1, SPRAT_MPI_INDEX_T, MPI_SUM, _workerComm), mpiErrorID);
		return globalResult;
	}
#else //#ifdef SPRAT_BUILD_WITH_MPI

public:
	int getTagFromName(std::string const& name) const {
		return 0;
	}
	void registerName(std::string const& name) const {}

	real globalSum(real localSummand) const {
		return localSummand;
	}
	index_t globalSum(index_t localSummand) const {
		return localSummand;
	}

#endif //#ifdef SPRAT_BUILD_WITH_MPI
};


#ifdef SPRAT_BUILD_WITH_MPI

class IndexVector : public std::vector<index_t> {
public:
	IndexVector() : std::vector<index_t>() {}
	explicit IndexVector(index_t size) : std::vector<index_t>(size) {}

	index_t inverse(index_t mappedValue) const {
		return std::find(begin(), end(), mappedValue) - begin();
	}

	//void send(ParallelExecutionEnvironment const& pEE, int toID) const {

	//}
	//void recv(ParallelExecutionEnvironment const& pEE, int fromID=0) {

	//}
};

class MPNeighborDoFIndexVector : public std::vector<index_t> {
public:
	int neighborID;
	//index_t startIndex;

	MPNeighborDoFIndexVector() : std::vector<index_t>(), neighborID(0) {}
	MPNeighborDoFIndexVector(index_t size, int neighborID_in) :
			std::vector<index_t>(size),
			neighborID(neighborID_in) {}


	void invertWith(IndexVector const& dictionary) {
		for(index_t i=0; i<size(); ++i) {
			this->operator[](i) = dictionary.inverse(this->operator[](i));
		}
	}


	void send(ParallelExecutionEnvironment const& pEE) const {
		//pEE.sendIndex(startIndex, neighborID);
		pEE.sendIndex(size(), neighborID);
		pEE.sendVectorContent(*this, neighborID);
	}
	void recv(ParallelExecutionEnvironment const& pEE, int fromID) {
		neighborID = fromID;
		//startIndex = pEE.recvIndex(neighborID);
		resize(pEE.recvIndex(neighborID));
		pEE.recvVectorContent(*this, neighborID);
	}

	MPI_Datatype createMPIType() const {
		int mpiErrorID;
		MPI_Datatype newType;
		std::vector<int> blockLengths(0); // WARNING: This should be std::vector<index_t> which is currently not supported by MPI.
		std::vector<int> displacements(0);

		displacements.push_back(this->operator[](0));
		for(index_t i=1; i<this->size(); ++i) {
			if(this->operator[](i) != this->operator[](i-1)+1) {
				// Start a new block
				blockLengths.push_back(this->operator[](i-1)-displacements.back()+1);
				displacements.push_back(this->operator[](i));
			}
		}
		blockLengths.push_back(this->back()-displacements.back()+1);

		//std::cout << "Number of Blocks for neighbor " << neighborID << ": " << blockLengths.size() << std::endl;
		//for(uint i=0; i<blockLengths.size(); ++i) {
		//	std::cout << "Block " << i << " starts at " << displacements[i] << " and is " << blockLengths[i] << " entries long." << std::endl;
		//}

		MPI_CALL(MPI_Type_indexed(blockLengths.size(),
				&blockLengths[0],
				&displacements[0],
            SPRAT_MPI_REAL,
            &newType), mpiErrorID);
		MPI_CALL(MPI_Type_commit(&newType), mpiErrorID);

		return newType;
	}
};

#endif //#ifdef SPRAT_BUILD_WITH_MPI


struct MPCommunicationRegistry {
	index_t nLocalDoF;
	index_t nDoF;
#ifdef SPRAT_BUILD_WITH_MPI
	IndexVector globalDofIndexVector;
	std::vector<MPNeighborDoFIndexVector> dofRecvFromNeighbor;
	std::vector<MPNeighborDoFIndexVector> dofSendToNeighbor;

	std::vector<MPNeighborDoFIndexVector> elementDoFRecvFromNeighbor;
	std::vector<MPNeighborDoFIndexVector> elementDoFSendToNeighbor;

	std::vector<MPI_Datatype> dofRecvFromNeighborMPIType;
	std::vector<MPI_Datatype> dofSendToNeighborMPIType;

	std::vector<MPI_Datatype> elementDoFRecvFromNeighborMPIType;
	std::vector<MPI_Datatype> elementDoFSendToNeighborMPIType;

	void createMPITypes() {
		dofRecvFromNeighborMPIType.resize(dofRecvFromNeighbor.size());
		for(uint i=0; i<dofRecvFromNeighborMPIType.size(); ++i) {
			dofRecvFromNeighborMPIType[i] = dofRecvFromNeighbor[i].createMPIType();
		}

		dofSendToNeighborMPIType.resize(dofSendToNeighbor.size());
		for(uint i=0; i<dofSendToNeighborMPIType.size(); ++i) {
			dofSendToNeighborMPIType[i] = dofSendToNeighbor[i].createMPIType();
		}


		elementDoFRecvFromNeighborMPIType.resize(elementDoFRecvFromNeighbor.size());
		for(uint i=0; i<elementDoFRecvFromNeighborMPIType.size(); ++i) {
			elementDoFRecvFromNeighborMPIType[i] = elementDoFRecvFromNeighbor[i].createMPIType();
		}

		elementDoFSendToNeighborMPIType.resize(elementDoFSendToNeighbor.size());
		for(uint i=0; i<elementDoFSendToNeighborMPIType.size(); ++i) {
			elementDoFSendToNeighborMPIType[i] = elementDoFSendToNeighbor[i].createMPIType();
		}
	}
#endif //#ifdef SPRAT_BUILD_WITH_MPI
};


#endif /* PARALLEL_HPP_ */
