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

#ifndef ELEMENT_ARRAYTYPES_HPP_
#define ELEMENT_ARRAYTYPES_HPP_

class ElementVectorArray {
protected:
	real * _values;
	FEMMeshT const* _femMesh;

	friend class DistributedElementVectorArray;
public:
	ElementVectorArray() : _values(0), _femMesh(0) {}
	ElementVectorArray(FEMMeshT const& femMesh, real initialValue = 0.0) {
		reset(femMesh, initialValue);
	}
	~ElementVectorArray() {
		numa_delete(_values);
	}

	void reset(FEMMeshT const& femMesh, real initialValue = 0.0) {
		_values = numa_new<real>(femMesh.nElements() * FEMMeshT::nDoFPerElement());
		_femMesh = &femMesh;
		touch(initialValue);
	}

	void touch(real value = 0.0) {
		foreachElementIndependently(auto tau, *_femMesh, shared(value) , {
			real * const basePtr = _values+((index_t)tau)*nDoFPerElement();
			for(uint i=0; i<FEMMeshT::nDoFPerElement(); ++i) {
				basePtr[i] = value;
			}
		})
	}

	void printToFile(std::string prefix, uint id, std::string suffix = ".txt") const {
		std::ofstream myfile;
		const int outputPrecision = 16;
		myfile.open((prefix + std::to_string(id) + suffix).c_str());
		for(index_t i=0; i<_femMesh->nElements() * nDoFPerElement(); ++i) {
			myfile << std::setprecision(outputPrecision) << std::scientific << _values[i] << std::endl;
		}
		myfile.close();
	}

	void writeBinaryToFile(std::string prefix, uint id, std::string suffix = ".bin") const {
		std::ofstream myfile;
		myfile.open((prefix + std::to_string(id) + suffix).c_str(), std::ios_base::out | std::ios_base::binary);
		for(index_t i=0; i<_femMesh->nElements() * nDoFPerElement(); ++i) {
			char * const location = (char *)&_values[i];
			myfile.write(location, sizeof(real));
		}
		myfile.close();
	}

	void readBinaryFromFile(std::string prefix, uint id, std::string suffix = ".bin") {
		std::ifstream myfile;
		myfile.open((prefix + std::to_string(id) + suffix).c_str(), std::ios_base::in | std::ios_base::binary);
		for(index_t i=0; i<_femMesh->nElements() * nDoFPerElement(); ++i) {
			char * const location = (char *)&_values[i];
			myfile.read(location, sizeof(real));
		}
		myfile.close();
	}

	void swapContent(ElementVectorArray& other) {
		real * myOldValues = _values;
		_values = other._values;
		other._values = myOldValues;
	}

	ElementVectorView operator[](ElementT const& element) const {
		return ElementVectorView(element, _values+((index_t)element)*nDoFPerElement());
	}
};


class DistributedElementVectorArray : public ElementVectorArray {
private:
	ParallelExecutionEnvironment const& _pEE;
	MPCommunicationRegistry const& _comRegistry;
#ifdef SPRAT_BUILD_WITH_MPI
	std::vector<MPI_Request> _sendRequests;
	std::vector<MPI_Request> _recvRequests;
#endif // #ifdef SPRAT_BUILD_WITH_MPI
	int _tag;

public:
	DistributedElementVectorArray(FEMMeshT const& femMesh, std::string const& name, ParallelExecutionEnvironment const& pEE, MPCommunicationRegistry const& comRegistry, real initialValue = 0.0) :
		ElementVectorArray(femMesh, initialValue),
		_pEE(pEE),
		_comRegistry(comRegistry),
#ifdef SPRAT_BUILD_WITH_MPI
		_sendRequests(0),
		_recvRequests(0),
#endif //#ifdef SPRAT_BUILD_WITH_MPI
		_tag(pEE.getTagFromName(name))
	{
#ifdef SPRAT_BUILD_WITH_MPI
		_sendRequests.reserve(2*pEE.nComputeProcesses());
		_recvRequests.reserve(2*pEE.nComputeProcesses());
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}

#ifdef SPRAT_BUILD_WITH_MPI
	void swapContent(ElementVectorArray& other) {
		finishPendingCommunication();

		real * myOldValues = _values;
		_values = other._values;
		other._values = myOldValues;
	}
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	void swapContent(DistributedElementVectorArray& other) {
		finishPendingCommunication();
		other.finishPendingCommunication();

		real * myOldValues = _values;
		_values = other._values;
		other._values = myOldValues;
	}

	void setTag(int tag) {
		_tag=tag;
	}
	int getTag() const {
		return _tag;
	}

	void exchangeDataAsync() {
#ifdef SPRAT_BUILD_WITH_MPI
		//std::cout << "exchangeDataAsync()" << " in Prozess " << _pEE.processID() << std::endl;
		finishPendingCommunication();

		for(uint i=0; i<_comRegistry.elementDoFRecvFromNeighbor.size(); ++i) {
			//std::cout << "Prozess " << _pEE.processID() << " empfängt von " << _comRegistry.elementDoFRecvFromNeighbor[i].neighborID << std::endl;
			_recvRequests.push_back(_pEE.recvAsyncMPIType(_values,
					_comRegistry.elementDoFRecvFromNeighborMPIType[i],
					_comRegistry.elementDoFRecvFromNeighbor[i].neighborID,
					getTag()));
		}
		for(uint i=0; i<_comRegistry.elementDoFSendToNeighbor.size(); ++i) {
			//std::cout << "Prozess " << _pEE.processID() << " sendet an " << _comRegistry.elementDoFSendToNeighbor[i].neighborID << std::endl;
			_sendRequests.push_back(_pEE.sendAsyncMPIType(_values,
					_comRegistry.elementDoFSendToNeighborMPIType[i],
					_comRegistry.elementDoFSendToNeighbor[i].neighborID,
					getTag()));
		}
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}

	void exchangeDataBlocked() {
		exchangeDataAsync();
		finishPendingCommunication();
	}

	void finishPendingSends() {
#ifdef SPRAT_BUILD_WITH_MPI
		//std::cout << "finishPendingSends()" << " in Prozess " << _pEE.processID() << std::endl;
		_pEE.waitFor(_sendRequests);
		_sendRequests.clear();
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}
	void finishPendingRecvs() {
#ifdef SPRAT_BUILD_WITH_MPI
		//std::cout << "finishPendingRecvs()" << " in Prozess " << _pEE.processID() << std::endl;
		_pEE.waitFor(_recvRequests);
		_recvRequests.clear();
#endif //#ifdef SPRAT_BUILD_WITH_MPI
	}
	void finishPendingCommunication() {
		//std::cout << "finishPendingCommunication()" << " in Prozess " << _pEE.processID() << std::endl;
		finishPendingSends();
		finishPendingRecvs();
	}
};





class LeightweightElementMatrixArray {
protected:
	real * _referenceValues;
	FEMMeshT const * _femMesh;

public:
	LeightweightElementMatrixArray() : _referenceValues(0) {} // Leave everything uninitialized for performance reasons.
	LeightweightElementMatrixArray(FEMMeshT const& femMesh, ElementMatrix const& reference)  : _referenceValues(0) {
		reset(femMesh, reference);
	}
	~LeightweightElementMatrixArray() {
		numa_delete(_referenceValues);
	}

	void reset(FEMMeshT const& femMesh, ElementMatrix const& reference) {
		_femMesh = &femMesh;

		numa_delete(_referenceValues);
		_referenceValues = numa_new<real>(nDoFPerElement() * nDoFPerElement());
		for(uint i=0; i<nDoFPerElement() * nDoFPerElement(); ++i) {
			_referenceValues[i] = reference.getValues()[i];
		}
	}

	//LightweightScaledElementMatrix operator[](ElementT const& element) const {
	//	return LightweightScaledElementMatrix(element, ElementScalingFunction<ScalingFunction>(element), _referenceValues);
	//}

//	boost::proto::result_of::make_expr<boost::proto::tag::terminal, ElementMVDomain, BareLightweightScaledElementMatrix>::type
//	operator[](ElementT const& element) const {
//		return boost::proto::make_expr<boost::proto::tag::terminal, ElementMVDomain>(
//				BareLightweightScaledElementMatrix(element, _scalingFunction(element), _referenceValues)
//		);
//	}
};



class BorderElementMatrixArray {
private:
	numa_vector<real> _matrixValues; // nDoFPerElement()^2-matrices in row-major order contiguously aligned in memory
	real _zeroMatrix[nDoFPerElement()*nDoFPerElement()];
	FEMMeshT const * _femMesh;

	void init_zeroMatrix() {
		for(uint i=0; i<nDoFPerElement()*nDoFPerElement(); ++i) {
			_zeroMatrix[i] = 0.0;
		}
	}

public:
	BorderElementMatrixArray() : _matrixValues(0), _femMesh(0) {
		init_zeroMatrix();
	}
	BorderElementMatrixArray(FEMMeshT const& femMesh)  : _matrixValues(0), _femMesh(&femMesh) {
		init_zeroMatrix();
		reset(femMesh);
	}

	void reset(FEMMeshT const& femMesh) {
		_matrixValues.resize(femMesh.nElementsWithDomainBorder()* nDoFPerElement()*nDoFPerElement(), 0.0);
		_femMesh = &femMesh;
		touch();
	}

	void touch(real value = 0.0) {
		foreachElementIndependently(auto tau, (*_femMesh), shared(value), {
			if(tau.hasDomainBorder()) {
				real * const basePtr = &_matrixValues[tau.iAmTheNThElementWithDomainBorder() * nDoFPerElement()*nDoFPerElement()];
				for(uint i=0; i<nDoFPerElement()*nDoFPerElement(); ++i) {
					basePtr[i] = value;
				}
			}
		})
	}

	ElementMatrixView operator[](ElementT const& element) {
		if(element.hasDomainBorder()) {
			return ElementMatrixView(element, &_matrixValues[element.iAmTheNThElementWithDomainBorder() * nDoFPerElement()*nDoFPerElement()]);
		}
		return ElementMatrixView(element, _zeroMatrix);
	}

	void print(std::string name) {
		std::cout << "This is BorderElementMatrixArray " << name << std::endl;
		const index_t count = _matrixValues.size()/(nDoFPerElement()*nDoFPerElement());
		for(index_t k=0; k<count; ++k) {
			for(uint i=0; i<nDoFPerElement(); ++i) {
				for(uint j=0; j<nDoFPerElement(); ++j) {
					printf("%12.4e", _matrixValues[k*nDoFPerElement()*nDoFPerElement() + i*nDoFPerElement() + j]);
					if(j!=nDoFPerElement()-1) {
						std::cout << ", ";
					}
				}
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
	}

//	boost::proto::result_of::make_expr<boost::proto::tag::terminal, ElementMVDomain, BareLightweightScaledElementMatrix>::type
//	operator[](ElementT const& element) const {
//		return boost::proto::make_expr<boost::proto::tag::terminal, ElementMVDomain>(
//				BareLightweightScaledElementMatrix(element, _scalingFunction(element), _referenceValues)
//		);
//	}
};

#endif //ELEMENT_ARRAYTYPES_HPP_
