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

#ifndef MESH_RECT_HPP_
#define MESH_RECT_HPP_

#include <vector>
#include <utility>
#include <algorithm>
#include <cstdlib>
#include <limits>
#include <forward_list>
#include <boost/iterator/iterator_facade.hpp>
#include "config.hpp"
#include "util.hpp"
#include "interval.hpp"
#include "numa.hpp"
#ifdef SPRAT_BUILD_WITH_MPI
#include "parallel.hpp"
#endif


struct RectMeshDimension {
	std::string name;

	uint nElements;
	uint nNodes;

	//real totalIntervalLength;

	std::vector<real> nodes;
	std::vector<real> elementDiameter;

	real minElementDiameter;


	RectMeshDimension() {}

	RectMeshDimension(std::string name_in, Interval x_in, uint nElements_in, bool logarithmic = false) :
		name(name_in),
		nElements(nElements_in),
		nNodes(nElements_in+1),
		nodes(nElements_in+1),
		elementDiameter(nElements_in)
	{
		if(!logarithmic) {
			for(auto i : UIntRange(0, nNodes)) {
				nodes[i] = x_in.a + ((real)i/(real)nElements) * (x_in.b - x_in.a);
			}
		}
		else {
			const real L_a = log10(x_in.a);
			const real L_b = log10(x_in.b);

			for(auto i : UIntRange(0, nNodes)) {
				nodes[i] = pow(10.0, L_a + ((real)i/(real)nElements) * (L_b - L_a));
			}
		}

		for(auto i : UIntRange(0, nElements)) {
			elementDiameter[i] = nodes[i+1] - nodes[i];
		}
		minElementDiameter = *std::min_element(elementDiameter.cbegin(), elementDiameter.cend());
	}

#ifdef SPRAT_BUILD_WITH_MPI
	void send(ParallelExecutionEnvironment const& pEE, int toID) const {
		pEE.sendUInt(nElements, toID);
		pEE.sendUInt(nNodes, toID);
		pEE.sendReal(minElementDiameter, toID);
		pEE.sendVectorContent(nodes, toID);
		pEE.sendVectorContent(elementDiameter, toID);
	}
	void recv(ParallelExecutionEnvironment const& pEE, int fromID=0) {
		name = "";
		nElements = pEE.recvUInt(fromID);
		nNodes = pEE.recvUInt(fromID);
		minElementDiameter = pEE.recvReal(fromID);
		nodes.resize(nNodes);
		pEE.recvVectorContent(nodes);
		elementDiameter.resize(nElements);
		pEE.recvVectorContent(elementDiameter);
	}
#endif

	void print(std::string name="") const {
		if(name.length()>0) {
			std::cout << "This is RectMeshDimension " << name << "." << std::endl;
		}
		else {
			std::cout << "This is a RectMeshDimension." << std::endl;
		}

		std::cout << "I consist of " << nElements << " Elements and " << nNodes << " Nodes." << std::endl;
		std::cout << "My Interval is [" << nodes[0] << ", " << nodes[nNodes-1] << "] and my shortest element diameter is " << minElementDiameter << "." << std::endl;
		std::cout << "My nodes: " << std::endl;
		printVector(nodes);
		std::cout << "My element diameters: " << std::endl;
		printVector(elementDiameter);
		std::cout << std::endl;
	}
};




template <unsigned int Dim>
class FEMMeshRectP1 {
public:
	constexpr static index_t noNeighbor       = INDEX_T_NO_NEIGHBOR;
	constexpr static index_t noLocalNeighbor  = INDEX_T_NO_LOCAL_NEIGHBOR;
	constexpr static index_t maxValidNeighbor = INDEX_T_MAX_VALID;

	constexpr static unsigned int nDimensions() {
		return Dim;
	}
	constexpr static unsigned int nDoFPerElement() {
		return twoToThePowerOfN(Dim);
	}
	constexpr static unsigned int nHypersufacesPerElement() {
		return 2*Dim;
	}
private:
	struct ElementData {
		short_index_t multiIndex[Dim];
		index_t dofIndices[nDoFPerElement()];
		index_t neighborIndices[nHypersufacesPerElement()];

#ifdef SPRAT_BUILD_WITH_MPI
		void send(ParallelExecutionEnvironment const& pEE, int toID) const {
			for(uint i=0; i<Dim; ++i) {
				pEE.sendShortIndex(multiIndex[i], toID);
			}
			for(uint i=0; i<nDoFPerElement(); ++i) {
				pEE.sendIndex(dofIndices[i], toID);
			}
			for(uint i=0; i<nHypersufacesPerElement(); ++i) {
				pEE.sendIndex(neighborIndices[i], toID);
			}
		}
		void recv(ParallelExecutionEnvironment const& pEE, int fromID=0) {
			for(uint i=0; i<Dim; ++i) {
				multiIndex[i] = pEE.recvShortIndex(fromID);
			}
			for(uint i=0; i<nDoFPerElement(); ++i) {
				dofIndices[i] = pEE.recvIndex(fromID);
			}
			for(uint i=0; i<nHypersufacesPerElement(); ++i) {
				neighborIndices[i] = pEE.recvIndex(fromID);
			}
		}
#endif
		void print(index_t i) const {
			std::cout << "Element " << i << " mit Multiindex ( ";
			for(uint i=0; i<Dim; ++i) {
				std::cout << multiIndex[i] << " ";
			}
			std::cout << ") und Freiheitsgraden ";
			for(uint i=0; i<nDoFPerElement(); ++i) {
				std::cout << dofIndices[i] << " ";
			}
			std::cout << std::endl;
		}
	};
	struct DoFData {
		short_index_t multiIndex[Dim];

#ifdef SPRAT_BUILD_WITH_MPI
		void send(ParallelExecutionEnvironment const& pEE, int toID) const {
			for(uint i=0; i<Dim; ++i) {
				pEE.sendShortIndex(multiIndex[i], toID);
			}
		}
		void recv(ParallelExecutionEnvironment const& pEE, int fromID=0) {
			for(uint i=0; i<Dim; ++i) {
				multiIndex[i] = pEE.recvShortIndex(fromID);
			}
		}
#endif
		void print(index_t i) const {
			std::cout << "DoF " << i << " mit Multiindex ( ";
			for(uint i=0; i<Dim; ++i) {
				std::cout << multiIndex[i] << " ";
			}
			std::cout << ")" << std::endl;
		}
	};

	numa_vector<ElementData> _elementData;
	std::vector<bool> _elementHasDomainBorder;
	std::vector<index_t> _elementsWithDomainBorder;
	numa_vector<DoFData> _dofData;
	index_t _nLocalElements;
	index_t _nLocalDoF;

#ifdef SPRAT_BUILD_WITH_OPENMP
	std::vector<std::forward_list<index_t>> _elementsOfDoF;
	std::vector<numa_vector<index_t>> _independentElementDecomposition;
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP

	constexpr static unsigned int twoToThePowerOfN(unsigned int n) {
		return (n>0 ? 2*twoToThePowerOfN(n-1) : 1);
	}

	index_t multiIndexToElementIndexInFullyPopulatedMesh(short_index_t multiIndex[Dim]) {
		index_t result = multiIndex[0];
		for(uint i=1; i<Dim; ++i) {
			result *= dimensions[i].nElements;
			result += multiIndex[i];
		}
		return result;
	}

	void findElementsWithDomainBorder() {
		_elementsWithDomainBorder.resize(0);
		_elementHasDomainBorder.resize(_elementData.size());
		for(index_t i=0; i<_elementData.size(); ++i) {
			_elementHasDomainBorder[i] = false;
			for(uint k=0; k<nHypersufacesPerElement(); ++k) {
				if(_elementData[i].neighborIndices[k] == noNeighbor) {
					_elementHasDomainBorder[i] = true;
					_elementsWithDomainBorder.push_back(i);
					break;
				}
			}
		}
	}

	void createElementsOfDoF(index_t numDoF, std::vector<ElementData> const& elemData) {
#ifdef SPRAT_BUILD_WITH_OPENMP
		_elementsOfDoF.resize(numDoF);

		for(index_t i=0; i<numDoF; ++i) {
			_elementsOfDoF[i].clear();
		}

		for(index_t i=0; i<elemData.size(); ++i) {
			for(uint k=0; k<nDoFPerElement(); ++k) {
				_elementsOfDoF[elemData[i].dofIndices[k]].push_front(i);
			}
		}
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP
	}

	void createMaximalIndependentElementSets(std::vector<ElementData> const& elemData) {
#ifdef SPRAT_BUILD_WITH_OPENMP
		static constexpr uint noColor = std::numeric_limits<uint>::max();
		//static constexpr uint noColor = UINT_MAX; // For Intel compiler bug
		std::vector<uint> elementColor(elemData.size(), noColor);
		index_t nUncolored = elemData.size();
		uint nColors = 0;

		for(uint currentColor=0; nUncolored>0; ++currentColor) {
			for(index_t i=0; i<elemData.size(); ++i) {
				if(elementColor[i] != noColor) {
					continue;
				}

				bool hasNeighbourOfCurrentColor = false;
				for(uint k=0; k<nDoFPerElement() && !hasNeighbourOfCurrentColor; ++k) {
					const index_t dof = elemData[i].dofIndices[k];
					for(auto it=_elementsOfDoF[dof].cbegin(); it!=_elementsOfDoF[dof].cend(); ++it) {
						if(elementColor[*it] == currentColor) {
							hasNeighbourOfCurrentColor = true;
							break;
						}
					}
				}
				if(!hasNeighbourOfCurrentColor) {
					elementColor[i] = currentColor;
					nUncolored -= 1;
				}
			}

			nColors += 1;
		}

		std::vector<std::vector<index_t>> independentElementDecompositionTemp(nColors);
		for(uint color=0; color<nColors; ++color) {
			independentElementDecompositionTemp[color].clear();
		}
		for(index_t i=0; i<elemData.size(); ++i) {
			independentElementDecompositionTemp[elementColor[i]].push_back(i);
		}

		_independentElementDecomposition.resize(nColors);
		for(uint color=0; color<nColors; ++color) {
			_independentElementDecomposition[color].resize(independentElementDecompositionTemp[color].size());
			foreach_omp_index(i, _independentElementDecomposition[color].size(), shared(color,independentElementDecompositionTemp), {
					_independentElementDecomposition[color][i] = independentElementDecompositionTemp[color][i];
			})
		}


		// Debug output
		//for(index_t i=0; i<nElements(); ++i) {
		//	std::cout << "Element " << i << " hat Farbe " << elementColor[i] << std::endl;
		//}
		//std::cout << std::endl;
		//for(uint color=0; color<nColors; ++color) {
		//	std::cout << "Farbe " << color << " haben die Elemente:";
		//	for(index_t i=0; i<_independentElementDecomposition[color].size(); ++i) {
		//		std::cout << " " << _independentElementDecomposition[color][i];
		//	}
		//	std::cout << std::endl;
		//}
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP
	}

public:
	RectMeshDimension dimensions[Dim];

	FEMMeshRectP1() {}

	FEMMeshRectP1(std::vector<RectMeshDimension> const& dimensions_in) :
		FEMMeshRectP1(&dimensions_in[0]) {}

	FEMMeshRectP1(RectMeshDimension const * const dimensions_in) {
		index_t totalElements = 1;
		index_t totalDoF = 1;
		for(uint i=0; i<Dim; ++i) {
			dimensions[i] = dimensions_in[i];
			totalElements *= dimensions[i].nElements;
			totalDoF *= dimensions[i].nNodes;
		}

		std::vector<ElementData> elementDataTemp(totalElements);
		_nLocalElements = totalElements;
		_dofData.resize(totalDoF);
		_nLocalDoF = totalDoF;
		//std::cout << "FEMMeshRectP1::FEMMeshRectP1: totalDoF="<<totalDoF<<"_dofData.size()="<<_dofData.size()<<std::endl;

		for(index_t index=0; index<totalElements; ++index) {
			for(uint k=0; k<Dim; ++k) {
				// *** multiIndex
				index_t divideAway = index;
				for(uint j=Dim-1; j>k; --j) {
					divideAway /= dimensions[j].nElements;
				}
				elementDataTemp[index].multiIndex[k] = divideAway % dimensions[k].nElements;
			}
			for(uint k=0; k<Dim; ++k) {
				// *** neighborIndices
				if(elementDataTemp[index].multiIndex[k] == 0) { // find left neighbor
					elementDataTemp[index].neighborIndices[2*k  ] = noNeighbor;
				}
				else {
					short_index_t leftNeighborIndex[Dim];
					for(uint j=0; j<Dim; ++j) {
						leftNeighborIndex[j] = elementDataTemp[index].multiIndex[j];
					}
					leftNeighborIndex[k] -= 1;
					elementDataTemp[index].neighborIndices[2*k  ] = multiIndexToElementIndexInFullyPopulatedMesh(leftNeighborIndex);
				}

				if(elementDataTemp[index].multiIndex[k] == dimensions[k].nElements-1) { // find right neighbor
					elementDataTemp[index].neighborIndices[2*k+1] = noNeighbor;
				}
				else {
					short_index_t rightNeighborIndex[Dim];
					for(uint j=0; j<Dim; ++j) {
						rightNeighborIndex[j] = elementDataTemp[index].multiIndex[j];
					}
					rightNeighborIndex[k] += 1;
					elementDataTemp[index].neighborIndices[2*k+1] = multiIndexToElementIndexInFullyPopulatedMesh(rightNeighborIndex);
				}
			}
			for(uint k=0; k<nDoFPerElement(); ++k) {
				// *** dofIndices
				uint dofPosition[Dim];
				uint divideAway = k;
				for(uint j=0; j<Dim; ++j) {
					dofPosition[Dim-1-j] = divideAway%2;
					divideAway /= 2;
				}
				index_t dofIndex = 0;
				for(uint j=0; j<Dim; ++j) {
					dofIndex *= dimensions[j].nNodes;
					dofIndex += elementDataTemp[index].multiIndex[j] + dofPosition[j];
				}
				elementDataTemp[index].dofIndices[k] = dofIndex;
			}
		}


		createElementsOfDoF(totalDoF, elementDataTemp);
		createMaximalIndependentElementSets(elementDataTemp);

		_elementData.resize(elementDataTemp.size());
		foreachElementIndependently(auto tau, *this, shared(elementDataTemp), {
				_elementData[tau] = elementDataTemp[tau];
		})
		elementDataTemp.clear();

		findElementsWithDomainBorder();

		//std::cout << "There are " << _elementData.size() << " elements in the master mesh." << std::endl;


		foreach_omp_index(index, totalDoF, , {
			for(uint k=0; k<Dim; ++k) {
				index_t divideAway = index;
				for(uint j=Dim-1; j>k; --j) {
					divideAway /= dimensions[j].nNodes;
				}
				_dofData[index].multiIndex[k] = divideAway % dimensions[k].nNodes;
			}
		})
	}

#ifdef SPRAT_BUILD_WITH_MPI
	void receiveStrideDecomposition(ParallelExecutionEnvironment const& pEE, MPCommunicationRegistry& comRegistry) {
		/*
		 * Empfange:
		 * - Dim RectMeshDimensions
		 * - Anzahl an Elementen
		 * - Anzahl an lokalen Elementen
		 * - Die Elemente (ElementData + den globalen Index des Elements)  -- zuerst die eigenen Elemente!
		 *   Es ist garantiert, dass die globalen Indizes der empfangen Elemente aufsteigend sind (gilt natürlich nur jeweils für eigene und Nachbarelemente)
		 * - Anzahl an DoF
		 * - Anzahl an lokalen DoF
		 * - Die DoF (DoFData + den globalen Index des DoF + ID dem der Knoten gehört) -- zuerst die eigenen DoF!
		 *   Es ist garantiert, dass die DoF der Nachbarn in sich auch jeweils zusammenhängen (daher der Sort beim senden!)
		 *
		 */

		//sleep(pEE.processID());

		for(uint k=0; k<Dim; ++k) {
			dimensions[k].recv(pEE); // Receive Dimensions
			//dimensions[k].print();
		}

		// Receive Elements
		std::vector<ElementData> elementDataTemp(pEE.recvIndex());
		//_elementData.resize(pEE.recvIndex());
		IndexVector globalElementIDs(elementDataTemp.size());
		_nLocalElements = pEE.recvIndex();
		for(index_t i=0; i<elementDataTemp.size(); ++i) {
			elementDataTemp[i].recv(pEE);
			//elementDataTemp[i].print(i);
			globalElementIDs[i] = pEE.recvIndex();
		}

		// Receive local DoF
		std::vector<DoFData> dofDataTemp(pEE.recvIndex());
		//_dofData.resize(pEE.recvIndex());
		comRegistry.nDoF = dofDataTemp.size();
		comRegistry.nLocalDoF = _nLocalDoF = pEE.recvIndex();
		comRegistry.globalDofIndexVector.resize(dofDataTemp.size());
		for(index_t i=0; i<_nLocalDoF; ++i) {
			dofDataTemp[i].recv(pEE);
			//dofDataTemp[i].print(i);
			comRegistry.globalDofIndexVector[i] = pEE.recvIndex();
			pEE.recvInt();
		}
		// Receive ghost DoF
		comRegistry.dofRecvFromNeighbor.resize(0);
		uint lastNeighbor = 0;
		for(index_t i=_nLocalDoF; i<dofDataTemp.size(); ++i) {
			dofDataTemp[i].recv(pEE);
			//dofDataTemp[i].print(i);
			comRegistry.globalDofIndexVector[i] = pEE.recvIndex();
			uint recvNeighbor = pEE.recvInt();
			if(lastNeighbor != recvNeighbor) {
				//std::cout << "Empfange Geisterpkt von " << recvNeighbor << std::endl;
				lastNeighbor = recvNeighbor;
				comRegistry.dofRecvFromNeighbor.push_back(MPNeighborDoFIndexVector(0, recvNeighbor));
			}
			comRegistry.dofRecvFromNeighbor.back().push_back(comRegistry.globalDofIndexVector[i]);
			//std::cout << i << " -> " << globalDofIndexVector[i] << std::endl;
		}


		// In den Element-Daten müssen die Indizes der DoF auf lokale angepasst werden..
		for(index_t i=0; i<elementDataTemp.size(); ++i) {
			for(uint k=0; k<nDoFPerElement(); ++k) {
				elementDataTemp[i].dofIndices[k] = comRegistry.globalDofIndexVector.inverse(elementDataTemp[i].dofIndices[k]);
			}
			//elementDataTemp[i].print(i);
		}


		//std::cout << "This is node " << pEE.processID() << ". I received " << elementDataTemp.size() << " Elements, "
		//		<< _nLocalDoF << " " << dofDataTemp.size() << std::endl;


		/*
		 * Für jeden Nachbarn: Ich schicke dem Nachbarn meine dofRecvFromNeighbor (das sind globale Indizes), die er mit regelmäßig zusenden soll
		 * Er schickt mir die dofSendToNeighbor, die diejenigen globalen Indized beinhalten, die ich ihm schicken soll.
		 */
		comRegistry.dofSendToNeighbor.resize(comRegistry.dofRecvFromNeighbor.size());
		for(uint i=0; i<pEE.nComputeProcesses(); ++i) { // Knoten i ist mit Senden dran
			for(uint k=0; k<comRegistry.dofRecvFromNeighbor.size(); ++k) {
				if(pEE.processID() == i+1) { // Wenn ich dran bin verschicke ich alle.
					comRegistry.dofRecvFromNeighbor[k].send(pEE);
				}
				else if(comRegistry.dofRecvFromNeighbor[k].neighborID == i+1) { // Wenn ich nicht dran bin empfange ich ggf. von dem der dran ist.
					//dofSendToNeighbor[k].neighborID = i+1;
					comRegistry.dofSendToNeighbor[k].recv(pEE, i+1);
				}
			}
		}


		//sleep(pEE.nComputeProcesses() + pEE.processID());

		// Wandele die empfangenen dofSendToNeighbor in ein local-to-local mapping um (d.h. dofSendToNeighbor[i] beinhaltet
		// dann diejenigen lokalen Indizes, die ich i in der dadurch vorgegebenen Reihenfolge schicken soll (also als erstes dofSendToNeighbor[i][0]...).
		for(uint i=0; i<comRegistry.dofSendToNeighbor.size(); ++i) {
			comRegistry.dofSendToNeighbor[i].invertWith(comRegistry.globalDofIndexVector);
			//std::cout << "Prozess " << pEE.processID() << " schickt an " << comRegistry.dofSendToNeighbor[i].neighborID << ":" << std::endl;
			//for(uint k=0; k<comRegistry.dofSendToNeighbor[i].size(); ++k) {
			//	std::cout << comRegistry.dofSendToNeighbor[i][k] << std::endl;
			//	assert(comRegistry.dofSendToNeighbor[i][k] < _nLocalDoF);
			//}
		}

		// Tue das gleiche auch für die zu empangnenden DoF; Denn da stehen ja bis jetzt *globale* DoF drin!
		for(uint i=0; i<comRegistry.dofRecvFromNeighbor.size(); ++i) {
			comRegistry.dofRecvFromNeighbor[i].invertWith(comRegistry.globalDofIndexVector);
			//std::cout << "Prozess " << pEE.processID() << " empfängt von " << comRegistry.dofRecvFromNeighbor[i].neighborID << ":" << std::endl;
			//for(uint k=0; k<comRegistry.dofRecvFromNeighbor[i].size(); ++k) {
			//	std::cout << comRegistry.dofRecvFromNeighbor[i][k] << std::endl;
			//	assert(comRegistry.dofRecvFromNeighbor[i][k] >= _nLocalDoF);
			//}
		}

		// Jetzt muss ich noch die elementDoFRecvFrom/SendToNeighbor aus dofRecvFrom/SendToNeighbor kreieren.
		// Die elementDoFRecvFrom/SendToNeighbor beschreiben das gleiche wie die dofRecvFrom/SendToNeighbor
		// nur für Arrays von Elementvektoren.
		// Dafür gehe ich für jeden zu sendenden/empfangenden DoF alle Nachbarelemente durch und notiere mir
		// alle Vorkomnisse dieses DoF nacheinander. Sobald die Nachbarelemente bei Sender und Empfänger die
		// gleiche Sortierung haben, können Sender und Empfänger ohne weitere Kommunikation einen
		// "Austauschplan" für die Element-DoF erzeugen.
		comRegistry.elementDoFRecvFromNeighbor.resize(comRegistry.dofRecvFromNeighbor.size());
		for(uint i=0; i<comRegistry.dofRecvFromNeighbor.size(); ++i) {
			comRegistry.elementDoFRecvFromNeighbor[i].resize(0);
			comRegistry.elementDoFRecvFromNeighbor[i].neighborID = comRegistry.dofRecvFromNeighbor[i].neighborID;
			for(index_t j=0; j<comRegistry.dofRecvFromNeighbor[i].size(); ++j) {
				for(index_t k=_nLocalElements; k<elementDataTemp.size(); ++k) {
					for(uint l=0; l<nDoFPerElement(); ++l) {
						if(elementDataTemp[k].dofIndices[l] == comRegistry.dofRecvFromNeighbor[i][j]) {
							comRegistry.elementDoFRecvFromNeighbor[i].push_back(k*nDoFPerElement() + l);
							//std::cout << "Prozess " << pEE.processID() << " empfängt von " << comRegistry.elementDoFRecvFromNeighbor[i].neighborID << " an Position " << comRegistry.elementDoFRecvFromNeighbor[i].size()-1 << " den lokalen ElementDoF " << k*nDoFPerElement() + l << std::endl;
							//assert(elementDataTemp[k].dofIndices[l] >= _nLocalDoF);
						}
					}
				}
			}
		}

		comRegistry.elementDoFSendToNeighbor.resize(comRegistry.dofSendToNeighbor.size());
		for(uint i=0; i<comRegistry.dofSendToNeighbor.size(); ++i) {
			comRegistry.elementDoFSendToNeighbor[i].resize(0);
			comRegistry.elementDoFSendToNeighbor[i].neighborID = comRegistry.dofSendToNeighbor[i].neighborID;
			for(index_t j=0; j<comRegistry.dofSendToNeighbor[i].size(); ++j) {
				//std::cout << "Prozess " << pEE.processID() << " sendet an " << comRegistry.dofSendToNeighbor[i].neighborID << " an Position " << j << " den lokalen DoF " << comRegistry.dofSendToNeighbor[i][j] << std::endl;
				for(index_t k=_nLocalElements; k<elementDataTemp.size(); ++k) {
					for(uint l=0; l<nDoFPerElement(); ++l) {
						if(elementDataTemp[k].dofIndices[l] == comRegistry.dofSendToNeighbor[i][j]) {
							comRegistry.elementDoFSendToNeighbor[i].push_back(k*nDoFPerElement() + l);
							//std::cout << "Prozess " << pEE.processID() << " sendet an " << comRegistry.elementDoFSendToNeighbor[i].neighborID << " an Position " << comRegistry.elementDoFSendToNeighbor[i].size()-1 << " den lokalen ElementDoF " << k*nDoFPerElement() + l << std::endl;
							//assert(elementDataTemp[k].dofIndices[l] < _nLocalDoF);
						}
					}
				}
			}
		}

		//sleep(pEE.nProcesses() + pEE.processID());

		comRegistry.createMPITypes();

		/*
		 * Abschließend muss ich die lokalen Indizes der Nachbarn der Elemente finden.
		 * Dazu suche ich alle globalen Elementnachbar in globalElementIDs und ersetze sie durch den Index in diesem
		 * Vektor. Falls sie nicht gefunden werden, muss ich sie auf noNeighbour setzen.
		 */
		for(index_t i=0; i<elementDataTemp.size(); ++i) {
			for(uint k=0; k<nHypersufacesPerElement(); ++k) {
				if(elementDataTemp[i].neighborIndices[k] != noNeighbor) {
					index_t localIndex = globalElementIDs.inverse(elementDataTemp[i].neighborIndices[k]);
					elementDataTemp[i].neighborIndices[k] = (localIndex < globalElementIDs.end() - globalElementIDs.begin() ? localIndex : noLocalNeighbor);
				}
			}
		}

		createElementsOfDoF(dofDataTemp.size(), elementDataTemp);
		createMaximalIndependentElementSets(elementDataTemp);

		_elementData.resize(elementDataTemp.size());
		foreachElementIndependently(auto tau, *this, shared(elementDataTemp), {
			_elementData[tau] = elementDataTemp[tau];
		})
		elementDataTemp.clear();
		_dofData.resize(dofDataTemp.size());
		foreach_omp_index(i, dofDataTemp.size(), shared(dofDataTemp), {
			_dofData[i] = dofDataTemp[i];
		})
		dofDataTemp.clear();

		findElementsWithDomainBorder();
	}

	void createAndDistributeStrideDecomposition(ParallelExecutionEnvironment const& pEE, std::vector<index_t>& dofStartIndexForProcess) const {
		/*
		 * Der Master ist dafür verantwortlich, die Elemente/DoF den Slaves so zu schicken, dass
		 * sie lokal richtig geordnet sind. So hat der Master auch gleich eine Zuordnung von
		 * lokalen zu globalen DoF/Elementen für jedes Teilgebiet.
		 */
		//localToGlobalDoF.resize(pEE.nComputeProcesses());

		// Bestimme wie viele DoF jeder Prozess bekommt.
		dofStartIndexForProcess.resize(pEE.nComputeProcesses()+1);
		index_t baseStep = nDoF()/pEE.nComputeProcesses();
		index_t remainder = nDoF()%pEE.nComputeProcesses();
		dofStartIndexForProcess[0] = 0;
		for(uint i=0; i<pEE.nComputeProcesses(); ++i) {
			dofStartIndexForProcess[i+1] = dofStartIndexForProcess[i] + baseStep + (i<remainder ? 1 : 0);
		}


		for(uint i=0; i<pEE.nComputeProcesses(); ++i) {
			for(uint k=0; k<Dim; ++k) {
				dimensions[k].send(pEE, i+1); // Send Dimensions
			}

			// Finde heraus, welche Elemente an Knoten i zu senden sind. Sammele ihre Indizes.
			std::vector<index_t> elementsForThisProcess(0);
			for(index_t j=0; j<_elementData.size(); ++j) {
				for(uint k=0; k<nDoFPerElement(); ++k) {
					// Liegt der Dof im i-ten Gebiet?
					index_t currentDoFIndex = _elementData[j].dofIndices[k];
					if(currentDoFIndex >= dofStartIndexForProcess[i] && currentDoFIndex < dofStartIndexForProcess[i+1]) {
						elementsForThisProcess.push_back(j);
						break;
					}
				}
			}

			// Zerteile die Elemente in i-rein-lokale und i-Nachbarelemente
			std::vector<index_t> localElements(0);
			std::vector<index_t> neighborElements(0);
			for(index_t j=0; j<elementsForThisProcess.size(); ++j) {
				uint nLocalDoF = 0;
				for(uint k=0; k<nDoFPerElement(); ++k) {
					index_t currentDoFIndex = _elementData[elementsForThisProcess[j]].dofIndices[k];
					if(currentDoFIndex >= dofStartIndexForProcess[i] && currentDoFIndex < dofStartIndexForProcess[i+1]) {
						nLocalDoF+=1;
					}
				}
				if(nLocalDoF == nDoFPerElement()) { // rein lokal
					localElements.push_back(elementsForThisProcess[j]);
				}
				else {
					neighborElements.push_back(elementsForThisProcess[j]);
				}
			}

			// Diese Sorts sind eigentlich überflüssig, da durch die Art der Konstruktion die aufsteigende globale Reihenfolge bereits festgelegt ist,.
			std::sort(localElements.begin(), localElements.end());
			std::sort(neighborElements.begin(), neighborElements.end());
			//std::cout << "Prozess " << i+1 << " bekommt " << localElements.size() + neighborElements.size() << " Elemente. Davon sind " << localElements.size() << " rein lokal." << std::endl;

			// Versende die gefundenen Elemente..
			pEE.sendIndex(localElements.size() + neighborElements.size(), i+1); // nElements
			pEE.sendIndex(localElements.size(), i+1); // nLocalElements
			for(index_t j=0; j<localElements.size(); ++j) {
				_elementData[localElements[j]].send(pEE, i+1);
				pEE.sendIndex(localElements[j], i+1); // der globale Index
			}
			for(index_t j=0; j<neighborElements.size(); ++j) {
				_elementData[neighborElements[j]].send(pEE, i+1);
				pEE.sendIndex(neighborElements[j], i+1); // der globale Index
			}

			// Finde heraus, welche Nodes *zusätzlich* an den i-ten Knoten gesendet werden müssen und wem sie gehören.
			// Iteriere dazu über die Elemente, die uns gehören und notiere alle Knoten, die nicht unsere sind.
			// Passe dabei aber auf, keine DoF doppelt hinzuzufügen!
			std::vector<std::pair<index_t, int>> additionalDoF; // Paar aus Knotenindex und Besitzer.
			for(index_t j=0; j<neighborElements.size(); ++j) {
				for(uint k=0; k<nDoFPerElement(); ++k) {
					// Liegt der Dof im i-ten Gebiet?
					index_t currentDoFIndex = _elementData[neighborElements[j]].dofIndices[k];
					if(currentDoFIndex < dofStartIndexForProcess[i] || currentDoFIndex >= dofStartIndexForProcess[i+1]) {
						int owner = 0;
						for(uint l=0; l<dofStartIndexForProcess.size()-1; ++l) {
							if(currentDoFIndex >= dofStartIndexForProcess[l] && currentDoFIndex < dofStartIndexForProcess[l+1]) {
								owner = l+1;
								break;
							}
						}

						if(std::find(additionalDoF.begin(), additionalDoF.end(), std::pair<index_t, int>(currentDoFIndex, owner))==additionalDoF.end()) {
							additionalDoF.push_back(std::pair<index_t, int>(currentDoFIndex, owner));
						}
					}
				}
			}

			// Sortiere die Freiheitsgrade so, dass diejenigen, die zu einem Nachbarn gehören zusammenhängen.
			std::sort(additionalDoF.begin(), additionalDoF.end(),
					[](const std::pair<index_t, int>& lhs, const std::pair<index_t, int>& rhs) -> bool {
				if(lhs.second != rhs.second) {
					return lhs.second < rhs.second;
				}
				return lhs.first < rhs.first; // Dies hier sorgt dafür, dass die DoF für einen Nachbarn in sich noch nach globalem Index geordnetet sind.
			} );

			index_t nLocalDoFForThisProcess = dofStartIndexForProcess[i+1]-dofStartIndexForProcess[i];
			pEE.sendIndex(nLocalDoFForThisProcess + additionalDoF.size(), i+1); // nDoF
			pEE.sendIndex(nLocalDoFForThisProcess, i+1); // nLocalDoF
			//localToGlobalDoF[i].resize(nLocalDoFForThisProcess + additionalDoF.size());
			for(index_t j=0; j<nLocalDoFForThisProcess; ++j) {
				_dofData[dofStartIndexForProcess[i] + j].send(pEE, i+1);
				pEE.sendIndex(dofStartIndexForProcess[i] + j, i+1); // der globale Index
				//localToGlobalDoF[i][j] = dofStartIndexForProcess[i] + j;
				pEE.sendInt(i+1, i+1); // der besitzer
			}
			for(index_t j=0; j<additionalDoF.size(); ++j) {
				_dofData[additionalDoF[j].first].send(pEE, i+1);
				pEE.sendIndex(additionalDoF[j].first, i+1); // der globale Index
				//localToGlobalDoF[i][nLocalDoFForThisProcess + j] = additionalDoF[j].first;
				pEE.sendInt(additionalDoF[j].second, i+1); // der besitzer
			}
		}
	}
#endif // SPRAT_BUILD_WITH_MPI


#ifdef SPRAT_BUILD_WITH_OPENMP
	std::vector<numa_vector<index_t>> const& getIndependentElementDecomposition() const {
		return _independentElementDecomposition;
	}
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP

	index_t nElements() const {
		return _elementData.size();
	}
	index_t nLocalElements() const {
		return _nLocalElements;
	}
	index_t nElementsWithDomainBorder() const {
		return _elementsWithDomainBorder.size();
	}
	index_t nDoF() const {
		return _dofData.size();
	}
	index_t nLocalDoF() const {
		return _nLocalDoF;
	}


	class ElementT {
	public:
		constexpr static unsigned int nDoFPerElement() {
			return FEMMeshRectP1<Dim>::twoToThePowerOfN(Dim);
		}
	private:

		FEMMeshRectP1<Dim> const& _femMesh;
		const index_t _index;

		// WARNING: Initially, end must be vec.size()
		index_t binarySearchInSortedIndexVector(std::vector<index_t> const& vec, index_t begin, index_t end, index_t value) const {
			const index_t mid = (begin+end)/2;
			const index_t mid_value = vec[mid];
			if(mid_value == value) {
				return mid;
			}
			else if(value < mid_value) {
				return binarySearchInSortedIndexVector(vec, begin, mid, value);
			}
			else {
				return binarySearchInSortedIndexVector(vec, mid, end, value);
			}
		}

	public:
		ElementT(FEMMeshRectP1<Dim> const& femMesh, index_t index) : _femMesh(femMesh), _index(index) {}

		operator index_t() const {
			return _index;
		}
		index_t index() const {
			return _index;
		}

		uint indexInDimension(uint dimension) const {
			return _femMesh._elementData[_index].multiIndex[dimension];
		}

		real diamInDimension(uint dimension) const {
			return _femMesh.dimensions[dimension].elementDiameter[indexInDimension(dimension)];
		}

		index_t globalDoFIndex(uint localDoFIndex) const {
			return _femMesh._elementData[_index].dofIndices[localDoFIndex];
		}

		index_t const * globalDoFIndices() const {
			return _femMesh._elementData[_index].dofIndices;
		}

		FEMMeshRectP1<Dim> const& getMesh() const {
			return _femMesh;
		}

		bool hasDomainBorder() const {
			return _femMesh._elementHasDomainBorder[_index];
		}

		index_t iAmTheNThElementWithDomainBorder() const {
			//const index_t result = std::find(_femMesh._elementsWithDomainBorder.begin(), _femMesh._elementsWithDomainBorder.end(), _index) - _femMesh._elementsWithDomainBorder.begin();
			//if(result >= _femMesh._elementsWithDomainBorder.size()) {
			//	std::cout << "IAmTheNThElementWithDomainBorder: This should NOT happen!!!" << std::endl;
			//}

			return binarySearchInSortedIndexVector(_femMesh._elementsWithDomainBorder, 0, _femMesh._elementsWithDomainBorder.size(), _index);
		}
		uint nDomainBorderHypesurfaces() const {
			uint result = 0;
			for(uint i=0; i<nHypersufacesPerElement(); ++i) {
				if(neighborElementIndices()[i] == noNeighbor) {
					result += 1;
				}
			}
			return result;
		}
		index_t neighborElementIndex(short_index_t index) const {
			return _femMesh._elementData[_index].neighborIndices[index];
		}
		index_t const * neighborElementIndices() const {
			return _femMesh._elementData[_index].neighborIndices;
		}

		/*
		 * Terminology:
		 * (Global) DoF -- all DoF accessible by this compute node
		 * LocalDoF     -- DoF genuinely belonging to this compute node
		 * ElementDoF   -- DoF belonging to one element
		 */
		class ElementDoFT {
		private:
			uint _localIndex;
			index_t _globalIndex;

		public:
			ElementDoFT(uint localIndex, index_t globalIndex) : _localIndex(localIndex), _globalIndex(globalIndex) {
				//std::cout << "Constructing ElementDoFT with localIndex " << _localIndex << " and global index " << _globalIndex << std::endl;
			}

			operator uint() const {
				return _localIndex;
			}
			uint index() const {
				return _localIndex;
			}

			index_t globalIndex() const {
				return _globalIndex;
			}

			bool operator==(ElementDoFT const& other) const {
				return _localIndex == other._localIndex;
			}
			bool operator!=(ElementDoFT const& other) const {
				return _localIndex != other._localIndex;
			}
			bool operator<(ElementDoFT const& other) const {
				return _localIndex < other._localIndex;
			}
			bool operator<=(ElementDoFT const& other) const {
				return _localIndex <= other._localIndex;
			}
			bool operator>(ElementDoFT const& other) const {
				return _localIndex > other._localIndex;
			}
			bool operator>=(ElementDoFT const& other) const {
				return _localIndex >= other._localIndex;
			}
		};

		typename FEMMeshRectP1<Dim>::ElementT::ElementDoFT elementDoF(uint localIndex) const {
			return typename FEMMeshRectP1<Dim>::ElementT::ElementDoFT(localIndex, globalDoFIndex(localIndex));
		}

		class IterateElementDoF {
		private:
			index_t const*const _globalDoF;

		public:
			IterateElementDoF(typename FEMMeshRectP1<Dim>::ElementT element) :
				_globalDoF(element.globalDoFIndices()) {}


			class const_iterator : public boost::iterator_facade<
			const_iterator
			, typename FEMMeshRectP1<Dim>::ElementT::ElementDoFT const
			, boost::random_access_traversal_tag
			, typename FEMMeshRectP1<Dim>::ElementT::ElementDoFT const
			, signed_short_index_t> {
			public:
				//const_iterator() : _element(0), _elementDoFIndex(0) {}
				explicit const_iterator(index_t const* globalDoF) : _globalDoF(globalDoF), _localIndex(0) {}
				const_iterator(index_t const* globalDoF, short_index_t localIndex) : _globalDoF(globalDoF), _localIndex(localIndex) {}

			private:
				friend class boost::iterator_core_access;

				void increment() { ++_localIndex; }
				void decrement() { --_localIndex; }
				void advance(signed_short_index_t n) { _localIndex+=n; }
				signed_short_index_t distance_to(const_iterator const& other) const { return other._localIndex - this->_localIndex; }

				bool equal(const_iterator const& other) const {	return this->_localIndex == other._localIndex; }

				typename FEMMeshRectP1<Dim>::ElementT::ElementDoFT const dereference() const {
					return typename FEMMeshRectP1<Dim>::ElementT::ElementDoFT(_localIndex, _globalDoF[_localIndex]);
				}

				index_t const* _globalDoF;
				signed_short_index_t _localIndex;
			};
			typedef const_iterator iterator;

			iterator begin() const { return iterator(_globalDoF); }
			iterator end() const { return iterator(_globalDoF, FEMMeshRectP1<Dim>::nDoFPerElement()); }
			const_iterator cbegin() const { return const_iterator(_globalDoF); }
			const_iterator cend() const { return const_iterator(_globalDoF, FEMMeshRectP1<Dim>::nDoFPerElement()); }
		};

		class HypersurfaceT {
		public:
			constexpr static unsigned int nDoFPerHypersurface() {
				return Dim;
			}
		private:
			typename FEMMeshRectP1<Dim>::ElementT const& _element;
			short_index_t _surfaceIndex;
		public:
			HypersurfaceT(typename FEMMeshRectP1<Dim>::ElementT const& element, uint localIndex) :
				_element(element),
				_surfaceIndex(localIndex)
			{}

			operator short_index_t() const {
				return _surfaceIndex;
			}
			short_index_t index() const {
				return _surfaceIndex;
			}

			uint elementDoFIndex(uint hypersurfaceDoFIndex) const {
				const uint normalDimension = _surfaceIndex/2;
				const uint normalEntry = _surfaceIndex%2;
				const uint allOnes = ~(0U);
				const uint maskLowerPart = ~(allOnes << (Dim-1-normalDimension));
				const uint maskHigherPart = allOnes << (Dim-1-normalDimension);
				const uint result = (normalEntry<<(Dim-1-normalDimension))    |
				                    (hypersurfaceDoFIndex & maskLowerPart)    |
				                    ((hypersurfaceDoFIndex & maskHigherPart)<<1);
				//std::cout << "Für hypersurfaceDoFIndex=" << hypersurfaceDoFIndex << " von Hypersurface " << _surfaceIndex <<
				//		" ist elementDoFIndex=" << result << std::endl;
				return result;
			}

			class HypersurfaceDoFT {
			private:
				uint _hypersurfaceIndex;
				uint _elementIndex;

			public:
				HypersurfaceDoFT(uint hypersurfaceIndex, index_t elementIndex) : _hypersurfaceIndex(hypersurfaceIndex), _elementIndex(elementIndex) {}

				operator uint() const {
					return _hypersurfaceIndex;
				}
				uint index() const {
					return _hypersurfaceIndex;
				}

				uint elementDoFIndex() const {
					return _elementIndex;
				}

				bool operator==(HypersurfaceDoFT const& other) const {
					return _hypersurfaceIndex == other._hypersurfaceIndex;
				}
				bool operator!=(HypersurfaceDoFT const& other) const {
					return _hypersurfaceIndex != other._hypersurfaceIndex;
				}
				bool operator<(HypersurfaceDoFT const& other) const {
					return _hypersurfaceIndex < other._hypersurfaceIndex;
				}
				bool operator<=(HypersurfaceDoFT const& other) const {
					return _hypersurfaceIndex <= other._hypersurfaceIndex;
				}
				bool operator>(HypersurfaceDoFT const& other) const {
					return _hypersurfaceIndex > other._hypersurfaceIndex;
				}
				bool operator>=(HypersurfaceDoFT const& other) const {
					return _hypersurfaceIndex >= other._hypersurfaceIndex;
				}
			};

			typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT hypersurfaceDoF(uint hypersurfaceDoFIndex) const {
				return typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT(hypersurfaceDoFIndex, elementDoFIndex(hypersurfaceDoFIndex));
			}

			real surfaceIntegral(HypersurfaceDoFT const& i, HypersurfaceDoFT const& j, uint dimension) {
				uint i_elem = i.elementDoFIndex();
				uint j_elem = j.elementDoFIndex();

				const uint normalDimension = _surfaceIndex/2;
				if(normalDimension != dimension) {
					return 0.0;
				}

				uint iMultiIndex[Dim];
				uint jMultiIndex[Dim];

				for(uint k=0; k<Dim; ++k) {
					iMultiIndex[Dim-1-k] = i_elem%2;
					i_elem /= 2;
					jMultiIndex[Dim-1-k] = j_elem%2;
					j_elem /= 2;
				}


				const real intValue[2] = {1.0/3.0, 1.0/6.0};

				const bool isLeftSurface = (_surfaceIndex%2 == 0);
				real result = (isLeftSurface ? -1.0 : 1.0);
				for(uint k=0; k<Dim; ++k) {
					if(k != normalDimension) {
						result *= intValue[(iMultiIndex[k]!=jMultiIndex[k] ? 1 : 0)] * _element.diamInDimension(k);
					}
				}
				return result;
			}

			class IterateHypersurfaceDoF {
			private:
				typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT const& _hypersurface;

			public:
				IterateHypersurfaceDoF(typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT const& hypersurface) :
					_hypersurface(hypersurface) {
					//std::cout << "IterateHypersurfaceDoF(hypersurface), hypersurface = " << &hypersurface << std::endl;
				}


				class const_iterator : public boost::iterator_facade<
				const_iterator
				, typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT const
				, boost::random_access_traversal_tag
				, typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT const
				, signed_short_index_t> {
				public:
					//const_iterator() : _element(0), _elementDoFIndex(0) {}
					explicit const_iterator(typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT const& hypersurface) : _hypersurface(hypersurface), _localIndex(0) {
						//std::cout << "IterateHypersurfaceDoF::const_iterator(hypersurface), hypersurface = " << &hypersurface << std::endl;
					}
					const_iterator(typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT const& hypersurface, short_index_t localIndex) : _hypersurface(hypersurface), _localIndex(localIndex) {
						//std::cout << "IterateHypersurfaceDoF::const_iterator(hypersurface, index), hypersurface = " << &hypersurface << std::endl;
					}

				private:
					friend class boost::iterator_core_access;

					void increment() { ++_localIndex; }
					void decrement() { --_localIndex; }
					void advance(signed_short_index_t n) { _localIndex+=n; }
					signed_short_index_t distance_to(const_iterator const& other) const { return other._localIndex - this->_localIndex; }

					bool equal(const_iterator const& other) const {	return this->_localIndex == other._localIndex; }

					typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT const dereference() const {
						return _hypersurface.hypersurfaceDoF(_localIndex);
					}

					typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT const& _hypersurface;
					signed_short_index_t _localIndex;
				};
				typedef const_iterator iterator;

				iterator begin() const { return iterator(_hypersurface); }
				iterator end() const { return iterator(_hypersurface, FEMMeshRectP1<Dim>::ElementT::HypersurfaceT::nDoFPerHypersurface()); }
				const_iterator cbegin() const { return const_iterator(_hypersurface); }
				const_iterator cend() const { return const_iterator(_hypersurface, FEMMeshRectP1<Dim>::ElementT::HypersurfaceT::nDoFPerHypersurface()); }
			};
		};

		typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT hypersurface(uint surfaceIndex) const {
			return typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT(*this, surfaceIndex);
		}

		class IterateElementDomainBoundaryHypersurfaces {
		private:
			typename FEMMeshRectP1<Dim>::ElementT const& _element;
			const uint _nDomainBoundaryHypersurfaces;
		public:
			IterateElementDomainBoundaryHypersurfaces(typename FEMMeshRectP1<Dim>::ElementT const& element) :
				_element(element),
				_nDomainBoundaryHypersurfaces(element.nDomainBorderHypesurfaces()){}

			class const_iterator : public boost::iterator_facade<
			const_iterator
			, typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT const
			, boost::random_access_traversal_tag
			, typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT const
			, int> {
			public:
				const_iterator() : _element(0), _hypersurfaceIndex(0) {}
				//explicit const_iterator(FEMMeshRectP1<Dim> const& femMesh) : _femMesh(femMesh), _elementIndex(0) {}
				explicit const_iterator(typename FEMMeshRectP1<Dim>::ElementT const& element) : _element(&element), _hypersurfaceIndex(0) {}
				const_iterator(typename FEMMeshRectP1<Dim>::ElementT const& element, int hypersurfaceIndex) : _element(&element), _hypersurfaceIndex(hypersurfaceIndex) {}

			private:
				friend class boost::iterator_core_access;

				void increment() { ++_hypersurfaceIndex; }
				void decrement() { --_hypersurfaceIndex; }
				void advance(int n) { _hypersurfaceIndex+=n; }
				int distance_to(const_iterator const& other) const { return  other._hypersurfaceIndex - this->_hypersurfaceIndex; }

				bool equal(const_iterator const& other) const {	return this->_hypersurfaceIndex == other._hypersurfaceIndex; }

				typename FEMMeshRectP1<Dim>::ElementT::HypersurfaceT const dereference() const {
					uint nBorderSurface = 0;
					for(int i=0; i<nHypersufacesPerElement(); ++i) {
						if(_element->neighborElementIndex(i) == noNeighbor) {
							if(nBorderSurface == _hypersurfaceIndex) {
								nBorderSurface = i;
								break;
							}
							nBorderSurface += 1;
						}
					}
					return _element->hypersurface(nBorderSurface);
				}

				typename FEMMeshRectP1<Dim>::ElementT const* _element;
				int _hypersurfaceIndex;
			};
			typedef const_iterator iterator;

			iterator begin() const { return iterator(_element); }
			iterator end() const { return iterator(_element, _nDomainBoundaryHypersurfaces); }
			const_iterator cbegin() const { return const_iterator(_element); }
			const_iterator cend() const { return const_iterator(_element, _nDomainBoundaryHypersurfaces); }
		};
	};

	ElementT const element(index_t index) const {
		return ElementT(*this, index);
	}

	class DoFT {
	private:
		FEMMeshRectP1<Dim> const& _femMesh;
		const index_t _index;

	public:
		DoFT(FEMMeshRectP1<Dim> const& femMesh, index_t index) : _femMesh(femMesh), _index(index) {}

		operator index_t() const {
			return _index;
		}
		index_t index() const {
			return _index;
		}

		uint indexInDimension(uint dimension) const {
			return _femMesh._dofData[_index].multiIndex[dimension];
		}

		real positionInDimension(uint dimension) const {
			return _femMesh.dimensions[dimension].nodes[indexInDimension(dimension)];
		}
	};

	DoFT dof(index_t index) const {
		return DoFT(*this, index);
	}


	/*
	 * Iterations
	 */

	class IterateMeshElements {
	private:
		FEMMeshRectP1<Dim> const& _femMesh;
		const index_t _beginIndex;
		const index_t _endIndex;

	public:
		IterateMeshElements(FEMMeshRectP1<Dim> const& femMesh, index_t beginIndex, index_t endIndex) :
			_femMesh(femMesh), _beginIndex(beginIndex), _endIndex(endIndex) {}


		class const_iterator : public boost::iterator_facade<
		const_iterator
		, typename FEMMeshRectP1<Dim>::ElementT const
		, boost::random_access_traversal_tag
		, typename FEMMeshRectP1<Dim>::ElementT const
		, signed_index_t> {
		public:
			const_iterator() : _femMesh(0), _elementIndex(0) {}
			//explicit const_iterator(FEMMeshRectP1<Dim> const& femMesh) : _femMesh(femMesh), _elementIndex(0) {}
			const_iterator(FEMMeshRectP1<Dim> const& femMesh, index_t elementIndex) : _femMesh(&femMesh), _elementIndex(elementIndex) {}

		private:
			friend class boost::iterator_core_access;

			void increment() { ++_elementIndex; }
			void decrement() { --_elementIndex; }
			void advance(signed_index_t n) { _elementIndex+=n; }
			signed_index_t distance_to(const_iterator const& other) const { return  other._elementIndex - this->_elementIndex; }

			bool equal(const_iterator const& other) const {	return this->_elementIndex == other._elementIndex; }

			typename FEMMeshRectP1<Dim>::ElementT const dereference() const {
				return _femMesh->element(_elementIndex);
			}

			FEMMeshRectP1<Dim> const * _femMesh;
			signed_index_t _elementIndex;
		};
		typedef const_iterator iterator;

		iterator begin() const { return iterator(_femMesh, _beginIndex); }
		iterator end() const { return iterator(_femMesh, _endIndex); }
		const_iterator cbegin() const { return const_iterator(_femMesh, _beginIndex); }
		const_iterator cend() const { return const_iterator(_femMesh, _endIndex); }
	};

	class IterateDoF {
	private:
		FEMMeshRectP1<Dim> const& _femMesh;
		const index_t _beginIndex;
		const index_t _endIndex;
	public:
		IterateDoF(FEMMeshRectP1<Dim> const& femMesh, index_t beginIndex, index_t endIndex) :
			_femMesh(femMesh), _beginIndex(beginIndex), _endIndex(endIndex) {}


		class const_iterator : public boost::iterator_facade<
		const_iterator
		, typename FEMMeshRectP1<Dim>::DoFT const
		, boost::random_access_traversal_tag
		, typename FEMMeshRectP1<Dim>::DoFT const
		, signed_index_t> {
		public:
			//const_iterator() : _element(0), _elementDoFIndex(0) {}
			//explicit const_iterator(Mesh_2DRect const& mesh) : _mesh(mesh), _dofIndex(0) {}
			const_iterator(FEMMeshRectP1<Dim> const& femMesh, index_t dofIndex) : _femMesh(&femMesh), _dofIndex(dofIndex) {}

		private:
			friend class boost::iterator_core_access;

			void increment() { ++_dofIndex; }
			void decrement() { --_dofIndex; }
			void advance(signed_index_t n) { _dofIndex+=n; }
			signed_index_t distance_to(const_iterator const& other) const { return other._dofIndex - this->_dofIndex; }

			bool equal(const_iterator const& other) const {	return this->_dofIndex == other._dofIndex; }

			typename FEMMeshRectP1<Dim>::DoFT const dereference() const {
				return _femMesh->dof(_dofIndex);
			}

			FEMMeshRectP1<Dim> const* _femMesh;
			signed_index_t _dofIndex;
		};
		typedef const_iterator iterator;

		iterator begin() const { return iterator(_femMesh, _beginIndex); }
		iterator end() const { return iterator(_femMesh, _endIndex); }
		const_iterator cbegin() const { return const_iterator(_femMesh, _beginIndex); }
		const_iterator cend() const { return const_iterator(_femMesh, _endIndex); }
	};


	/*
	 * Returns \prod_{k=1}^n  \int_0^1 \tilde{\phi}_i_k (x_k) * \tilde{\phi}_j_k (x_k) dx_k
	 * Assertion: i and j are *local* DoF Indices (\in {0,1,2, ..., nDofPerElement()})
	 */
	static real integrateReferenceElement(uint i, uint j) {
		uint iMultiIndex[Dim];
		uint jMultiIndex[Dim];

		for(uint k=0; k<Dim; ++k) {
			iMultiIndex[Dim-1-k] = i%2;
			i /= 2;
			jMultiIndex[Dim-1-k] = j%2;
			j /= 2;
		}

		const real intValue[4] = {1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/3.0}; // (i_x/y, j_x/y): (0,0), (0,1), (1,0), (1,1)
		real result = 1.0;
		for(uint k=0; k<Dim; ++k) {
			result *= intValue[2*iMultiIndex[k] + jMultiIndex[k]];
		}
		return result;
	}

	/*
	 * Returns \prod_{k=1,k!=l}^n  \int_0^1 \tilde{\phi}_i_k (x_k) * \tilde{\phi}_j_k (x_k) dx_k *
	 *         \int_0^1 (d/dx_l \tilde{\phi}_i_l (x_l)) * \tilde{\phi}_j_l (x_l) dx
	 * Assertion: i and j are *local* DoF Indices (\in {0,1,2, ..., nDofPerElement()})
	 */
	static real integratePartIntDerivativeReferenceElement(uint i, uint j, uint derivativeDimension) {
		uint iMultiIndex[Dim];
		uint jMultiIndex[Dim];

		for(uint k=0; k<Dim; ++k) {
			iMultiIndex[Dim-1-k] = i%2;
			i /= 2;
			jMultiIndex[Dim-1-k] = j%2;
			j /= 2;
		}

		const real intValue[4] = {1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/3.0}; // (i_x/y, j_x/y): (0,0), (0,1), (1,0), (1,1)
		const real intDerivativeXValue[4] = {-0.5, -0.5, 0.5, 0.5}; // (i_x, j_x): (0,0), (0,1), (1,0), (1,1)
		real result = 1.0;
		for(uint k=0; k<Dim; ++k) {
			if(k == derivativeDimension) {
				result *= intDerivativeXValue[2*iMultiIndex[k] + jMultiIndex[k]];
			}
			else {
				result *= intValue[2*iMultiIndex[k] + jMultiIndex[k]];
			}
		}
		return result;
	}
};



template <unsigned int Dim>
static inline typename FEMMeshRectP1<Dim>::IterateMeshElements Elements(FEMMeshRectP1<Dim> const& femMesh) {
	return typename FEMMeshRectP1<Dim>::IterateMeshElements(femMesh, 0, femMesh.nElements());
}
template <unsigned int Dim>
static inline typename FEMMeshRectP1<Dim>::IterateMeshElements LocalElements(FEMMeshRectP1<Dim> const& femMesh) {
	return typename FEMMeshRectP1<Dim>::IterateMeshElements(femMesh, 0, femMesh.nLocalElements());
}
template <unsigned int Dim>
static inline typename FEMMeshRectP1<Dim>::IterateMeshElements NeighborElements(FEMMeshRectP1<Dim> const& femMesh) {
	return typename FEMMeshRectP1<Dim>::IterateMeshElements(femMesh, femMesh.nLocalElements(), femMesh.nElements());
}

template <unsigned int Dim>
static inline UIntRange ElementDoFIndices(FEMMeshRectP1<Dim> const& femMesh) {
	return UIntRange(0, FEMMeshRectP1<Dim>::nDoFPerElement());
}
template <typename ElemT>
static inline typename ElemT::IterateElementDoF ElementDoF(ElemT const& element) {
	return typename ElemT::IterateElementDoF(element);
}
template <typename ElemT>
static inline typename ElemT::IterateElementDomainBoundaryHypersurfaces DomainBoundaryHypersurfaces(ElemT const& element) {
	return typename ElemT::IterateElementDomainBoundaryHypersurfaces(element);
}
template <typename HSurfaceT>
static inline typename HSurfaceT::IterateHypersurfaceDoF HypersurfaceDoF(HSurfaceT const& hypersurface) {
	return typename HSurfaceT::IterateHypersurfaceDoF(hypersurface);
}

template <unsigned int Dim>
static inline typename FEMMeshRectP1<Dim>::IterateDoF DoF(FEMMeshRectP1<Dim> const& femMesh) {
	return typename FEMMeshRectP1<Dim>::IterateDoF(femMesh, 0, femMesh.nDoF());
}
template <unsigned int Dim>
static inline typename FEMMeshRectP1<Dim>::IterateDoF LocalDoF(FEMMeshRectP1<Dim> const& femMesh) {
	return typename FEMMeshRectP1<Dim>::IterateDoF(femMesh, 0, femMesh.nLocalDoF());
}
template <unsigned int Dim>
static inline typename FEMMeshRectP1<Dim>::IterateDoF GhostDoF(FEMMeshRectP1<Dim> const& femMesh) {
	return typename FEMMeshRectP1<Dim>::IterateDoF(femMesh, femMesh.nLocalDoF(), femMesh.nDoF());
}



#endif/* MESH_RECT_HPP_ */
