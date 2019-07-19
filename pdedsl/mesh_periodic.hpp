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

#ifndef MESH_PERIODIC_HPP_
#define MESH_PERIODIC_HPP_

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
#include "mesh_rect.hpp"




template <unsigned int Dim>
class FEMMeshRectP1Periodic {
public:
	constexpr static index_t noNeighbor       = INDEX_T_NO_NEIGHBOR;
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

		void print(index_t j) const {
			std::cout << "Element " << j << " mit Multiindex ( ";
			for(uint i=0; i<Dim; ++i) {
				std::cout << multiIndex[i] << " ";
			}
			std::cout << ") und Freiheitsgraden ";
			for(uint i=0; i<nDoFPerElement(); ++i) {
				std::cout << dofIndices[i] << " ";
			}
//			std::cout << " und Nachbarn ";
//			for(uint i=0; i<nHypersufacesPerElement(); ++i) {
//				std::cout << neighborIndices[i] << " ";
//			}
			std::cout << std::endl;
		}
	};
	struct DoFData {
		short_index_t multiIndex[Dim];

		void print(index_t j) const {
			std::cout << "DoF " << j << " mit Multiindex ( ";
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

	FEMMeshRectP1Periodic() {}

	FEMMeshRectP1Periodic(std::vector<RectMeshDimension> const& dimensions_in) :
		FEMMeshRectP1Periodic(&dimensions_in[0]) {}

	FEMMeshRectP1Periodic(RectMeshDimension const * const dimensions_in) {
		index_t totalElements = 1;
		index_t totalDoF = 1;
		for(uint i=0; i<Dim; ++i) {
			dimensions[i] = dimensions_in[i];
			totalElements *= dimensions[i].nElements;
			totalDoF *= dimensions[i].nNodes-1;
		}

		std::vector<ElementData> elementDataTemp(totalElements);
		_dofData.resize(totalDoF);

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

				// find left neighbor
				short_index_t leftNeighborIndex[Dim];
				for(uint j=0; j<Dim; ++j) {
					leftNeighborIndex[j] = elementDataTemp[index].multiIndex[j];
				}
				if(leftNeighborIndex[k] == 0) {
					leftNeighborIndex[k] = dimensions[k].nElements-1;
				} else {
					leftNeighborIndex[k] -= 1;
				}
				elementDataTemp[index].neighborIndices[2*k  ] = multiIndexToElementIndexInFullyPopulatedMesh(leftNeighborIndex);

				// find right neighbor
				short_index_t rightNeighborIndex[Dim];
				for(uint j=0; j<Dim; ++j) {
					rightNeighborIndex[j] = elementDataTemp[index].multiIndex[j];
				}
				if(rightNeighborIndex[k] == dimensions[k].nElements-1) {
					rightNeighborIndex[k] = 0;
				} else {
					rightNeighborIndex[k] += 1;
				}
				elementDataTemp[index].neighborIndices[2*k+1] = multiIndexToElementIndexInFullyPopulatedMesh(rightNeighborIndex);
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
					dofIndex *= dimensions[j].nNodes-1;
					if(!(elementDataTemp[index].multiIndex[j] == dimensions[j].nElements-1 && dofPosition[j] == 1)) {
						dofIndex += elementDataTemp[index].multiIndex[j] + dofPosition[j];
					}
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
					divideAway /= (dimensions[j].nNodes-1);
				}
				_dofData[index].multiIndex[k] = divideAway % (dimensions[k].nNodes-1);
			}
		})

//		for(uint i=0; i<_elementData.size(); ++i) {
//			_elementData[i].print(i);
//		}
//		for(uint i=0; i<_dofData.size(); ++i) {
//			_dofData[i].print(i);
//		}
	}


#ifdef SPRAT_BUILD_WITH_OPENMP
	std::vector<numa_vector<index_t>> const& getIndependentElementDecomposition() const {
		return _independentElementDecomposition;
	}
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP

	index_t nElements() const {
		return _elementData.size();
	}
	index_t nLocalElements() const {
		return nElements();
	}
	index_t nElementsWithDomainBorder() const {
		return _elementsWithDomainBorder.size();
	}
	index_t nDoF() const {
		return _dofData.size();
	}
	index_t nSpatialDoF() const {
		return _dofData.size() / dimensions[Dim-1].nNodes;
	}
	index_t nLocalDoF() const {
		return nDoF();
	}


	class ElementT {
	public:
		constexpr static unsigned int nDoFPerElement() {
			return FEMMeshRectP1Periodic<Dim>::twoToThePowerOfN(Dim);
		}
	private:

		FEMMeshRectP1Periodic<Dim> const& _femMesh;
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
		ElementT(FEMMeshRectP1Periodic<Dim> const& femMesh, index_t index) : _femMesh(femMesh), _index(index) {}

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

		FEMMeshRectP1Periodic<Dim> const& getMesh() const {
			return _femMesh;
		}


		/*
		 * Returns - \sum_{k=1}^d \int_\tau \frac{\partial}{\partial x_k} \prod_{l} \phi_i_l (x_l)  *  \frac{\partial}{\partial x_k} \prod_{m} \phi_i_m (x_m)  dx
		 * Assertion: i and j are *local* DoF Indices (\in {0,1,2, ..., nDofPerElement()})
		 */
		real integratePartIntLaplace(uint i, uint j) const {
			uint iMultiIndex[Dim];
			uint jMultiIndex[Dim];

			for(uint k=0; k<Dim; ++k) {
				iMultiIndex[Dim-1-k] = i%2;
				i /= 2;
				jMultiIndex[Dim-1-k] = j%2;
				j /= 2;
			}

			const real intValue[4] = {1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/3.0}; // (i_x/y, j_x/y): (0,0), (0,1), (1,0), (1,1)
			const real diffValue[4] = {1.0, -1.0, -1.0, 1.0}; // (i_x/y, j_x/y): (0,0), (0,1), (1,0), (1,1)
			real result = 0.0;
			for(uint k=0; k<Dim; ++k) {
				real dimResult = 1.0;
				for(uint l=0; l<Dim; ++l) {
					if(l != k) {
						dimResult *= diamInDimension(l) * intValue[2*iMultiIndex[l] + jMultiIndex[l]];
					}
					else {
						dimResult *= diffValue[2*iMultiIndex[l] + jMultiIndex[l]] / diamInDimension(l);
					}
				}
				result += dimResult;
			}
			return -1.0 * result;
		}


		bool hasDomainBorder() const {
			return _femMesh._elementHasDomainBorder[_index];
		}

		index_t iAmTheNThElementWithDomainBorder() const {
			//const index_t result = std::find(_femMesh._elementsWithDomainBorder.begin(), _femMesh._elementsWithDomainBorder.end(), _index) - _femMesh._elementsWithDomainBorder.begin();
			//if(result >= _femMesh._elementsWithDomainBorder.size()) {
			//	std::cout << "IAmTheNThElementWithDomainBorder: This should NOT happen!" << std::endl;
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
			uint localIndex() const {
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

		typename FEMMeshRectP1Periodic<Dim>::ElementT::ElementDoFT elementDoF(uint localIndex) const {
			return typename FEMMeshRectP1Periodic<Dim>::ElementT::ElementDoFT(localIndex, globalDoFIndex(localIndex));
		}

		class IterateElementDoF {
		private:
			index_t const*const _globalDoF;

		public:
			IterateElementDoF(typename FEMMeshRectP1Periodic<Dim>::ElementT element) :
				_globalDoF(element.globalDoFIndices()) {}


			class const_iterator : public boost::iterator_facade<
			const_iterator
			, typename FEMMeshRectP1Periodic<Dim>::ElementT::ElementDoFT const
			, boost::random_access_traversal_tag
			, typename FEMMeshRectP1Periodic<Dim>::ElementT::ElementDoFT const
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

				typename FEMMeshRectP1Periodic<Dim>::ElementT::ElementDoFT const dereference() const {
					return typename FEMMeshRectP1Periodic<Dim>::ElementT::ElementDoFT(_localIndex, _globalDoF[_localIndex]);
				}

				index_t const* _globalDoF;
				signed_short_index_t _localIndex;
			};
			typedef const_iterator iterator;

			iterator begin() const { return iterator(_globalDoF); }
			iterator end() const { return iterator(_globalDoF, FEMMeshRectP1Periodic<Dim>::nDoFPerElement()); }
			const_iterator cbegin() const { return const_iterator(_globalDoF); }
			const_iterator cend() const { return const_iterator(_globalDoF, FEMMeshRectP1Periodic<Dim>::nDoFPerElement()); }
		};

		class HypersurfaceT {
		public:
			constexpr static unsigned int nDoFPerHypersurface() {
				return Dim;
			}
		private:
			typename FEMMeshRectP1Periodic<Dim>::ElementT const& _element;
			short_index_t _surfaceIndex;
		public:
			HypersurfaceT(typename FEMMeshRectP1Periodic<Dim>::ElementT const& element, uint localIndex) :
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

			typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT hypersurfaceDoF(uint hypersurfaceDoFIndex) const {
				return typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT(hypersurfaceDoFIndex, elementDoFIndex(hypersurfaceDoFIndex));
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
				typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT const& _hypersurface;

			public:
				IterateHypersurfaceDoF(typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT const& hypersurface) :
					_hypersurface(hypersurface) {
					//std::cout << "IterateHypersurfaceDoF(hypersurface), hypersurface = " << &hypersurface << std::endl;
				}


				class const_iterator : public boost::iterator_facade<
				const_iterator
				, typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT const
				, boost::random_access_traversal_tag
				, typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT const
				, signed_short_index_t> {
				public:
					//const_iterator() : _element(0), _elementDoFIndex(0) {}
					explicit const_iterator(typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT const& hypersurface) : _hypersurface(hypersurface), _localIndex(0) {
						//std::cout << "IterateHypersurfaceDoF::const_iterator(hypersurface), hypersurface = " << &hypersurface << std::endl;
					}
					const_iterator(typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT const& hypersurface, short_index_t localIndex) : _hypersurface(hypersurface), _localIndex(localIndex) {
						//std::cout << "IterateHypersurfaceDoF::const_iterator(hypersurface, index), hypersurface = " << &hypersurface << std::endl;
					}

				private:
					friend class boost::iterator_core_access;

					void increment() { ++_localIndex; }
					void decrement() { --_localIndex; }
					void advance(signed_short_index_t n) { _localIndex+=n; }
					signed_short_index_t distance_to(const_iterator const& other) const { return other._localIndex - this->_localIndex; }

					bool equal(const_iterator const& other) const {	return this->_localIndex == other._localIndex; }

					typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT::HypersurfaceDoFT const dereference() const {
						return _hypersurface.hypersurfaceDoF(_localIndex);
					}

					typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT const& _hypersurface;
					signed_short_index_t _localIndex;
				};
				typedef const_iterator iterator;

				iterator begin() const { return iterator(_hypersurface); }
				iterator end() const { return iterator(_hypersurface, FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT::nDoFPerHypersurface()); }
				const_iterator cbegin() const { return const_iterator(_hypersurface); }
				const_iterator cend() const { return const_iterator(_hypersurface, FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT::nDoFPerHypersurface()); }
			};
		};

		typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT hypersurface(uint surfaceIndex) const {
			return typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT(*this, surfaceIndex);
		}

		class IterateElementDomainBoundaryHypersurfaces {
		private:
			typename FEMMeshRectP1Periodic<Dim>::ElementT const& _element;
			const uint _nDomainBoundaryHypersurfaces;
		public:
			IterateElementDomainBoundaryHypersurfaces(typename FEMMeshRectP1Periodic<Dim>::ElementT const& element) :
				_element(element),
				_nDomainBoundaryHypersurfaces(element.nDomainBorderHypesurfaces()){}

			class const_iterator : public boost::iterator_facade<
			const_iterator
			, typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT const
			, boost::random_access_traversal_tag
			, typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT const
			, int> {
			public:
				const_iterator() : _element(0), _hypersurfaceIndex(0) {}
				//explicit const_iterator(FEMMeshRectP1<Dim> const& femMesh) : _femMesh(femMesh), _elementIndex(0) {}
				explicit const_iterator(typename FEMMeshRectP1Periodic<Dim>::ElementT const& element) : _element(&element), _hypersurfaceIndex(0) {}
				const_iterator(typename FEMMeshRectP1Periodic<Dim>::ElementT const& element, int hypersurfaceIndex) : _element(&element), _hypersurfaceIndex(hypersurfaceIndex) {}

			private:
				friend class boost::iterator_core_access;

				void increment() { ++_hypersurfaceIndex; }
				void decrement() { --_hypersurfaceIndex; }
				void advance(int n) { _hypersurfaceIndex+=n; }
				int distance_to(const_iterator const& other) const { return  other._hypersurfaceIndex - this->_hypersurfaceIndex; }

				bool equal(const_iterator const& other) const {	return this->_hypersurfaceIndex == other._hypersurfaceIndex; }

				typename FEMMeshRectP1Periodic<Dim>::ElementT::HypersurfaceT const dereference() const {
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

				typename FEMMeshRectP1Periodic<Dim>::ElementT const* _element;
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
		FEMMeshRectP1Periodic<Dim> const& _femMesh;
		const index_t _index;

	public:
		DoFT(FEMMeshRectP1Periodic<Dim> const& femMesh, index_t index) : _femMesh(femMesh), _index(index) {}

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

		real spatialDistance(DoFT const& other) const {
			real result = 0.0;

			for(uint i=0; i<Dim-1; ++i) {
				const real otherPos = other.positionInDimension(i);
				const real thisPos = positionInDimension(i);
				const real a = _femMesh.dimensions[i].nodes.front();
				const real b = _femMesh.dimensions[i].nodes.back();

				const real diffNormal = fabs(otherPos - thisPos);
				const real diffWrappedA =
						otherPos - a +
						b - thisPos;
				const real diffWrappedB =
						thisPos - a +
						b - otherPos;

				const real diff = min(diffNormal, diffWrappedA, diffWrappedB);
				result += diff*diff;
			}

			return sqrt(result);
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
		FEMMeshRectP1Periodic<Dim> const& _femMesh;
		const index_t _beginIndex;
		const index_t _endIndex;

	public:
		IterateMeshElements(FEMMeshRectP1Periodic<Dim> const& femMesh, index_t beginIndex, index_t endIndex) :
			_femMesh(femMesh), _beginIndex(beginIndex), _endIndex(endIndex) {}


		class const_iterator : public boost::iterator_facade<
		const_iterator
		, typename FEMMeshRectP1Periodic<Dim>::ElementT const
		, boost::random_access_traversal_tag
		, typename FEMMeshRectP1Periodic<Dim>::ElementT const
		, signed_index_t> {
		public:
			const_iterator() : _femMesh(0), _elementIndex(0) {}
			//explicit const_iterator(FEMMeshRectP1<Dim> const& femMesh) : _femMesh(femMesh), _elementIndex(0) {}
			const_iterator(FEMMeshRectP1Periodic<Dim> const& femMesh, index_t elementIndex) : _femMesh(&femMesh), _elementIndex(elementIndex) {}

		private:
			friend class boost::iterator_core_access;

			void increment() { ++_elementIndex; }
			void decrement() { --_elementIndex; }
			void advance(signed_index_t n) { _elementIndex+=n; }
			signed_index_t distance_to(const_iterator const& other) const { return  other._elementIndex - this->_elementIndex; }

			bool equal(const_iterator const& other) const {	return this->_elementIndex == other._elementIndex; }

			typename FEMMeshRectP1Periodic<Dim>::ElementT const dereference() const {
				return _femMesh->element(_elementIndex);
			}

			FEMMeshRectP1Periodic<Dim> const * _femMesh;
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
		FEMMeshRectP1Periodic<Dim> const& _femMesh;
		const index_t _beginIndex;
		const index_t _endIndex;
	public:
		IterateDoF(FEMMeshRectP1Periodic<Dim> const& femMesh, index_t beginIndex, index_t endIndex) :
			_femMesh(femMesh), _beginIndex(beginIndex), _endIndex(endIndex) {}


		class const_iterator : public boost::iterator_facade<
		const_iterator
		, typename FEMMeshRectP1Periodic<Dim>::DoFT const
		, boost::random_access_traversal_tag
		, typename FEMMeshRectP1Periodic<Dim>::DoFT const
		, signed_index_t> {
		public:
			//const_iterator() : _element(0), _elementDoFIndex(0) {}
			//explicit const_iterator(Mesh_2DRect const& mesh) : _mesh(mesh), _dofIndex(0) {}
			const_iterator(FEMMeshRectP1Periodic<Dim> const& femMesh, index_t dofIndex) : _femMesh(&femMesh), _dofIndex(dofIndex) {}

		private:
			friend class boost::iterator_core_access;

			void increment() { ++_dofIndex; }
			void decrement() { --_dofIndex; }
			void advance(signed_index_t n) { _dofIndex+=n; }
			signed_index_t distance_to(const_iterator const& other) const { return other._dofIndex - this->_dofIndex; }

			bool equal(const_iterator const& other) const {	return this->_dofIndex == other._dofIndex; }

			typename FEMMeshRectP1Periodic<Dim>::DoFT const dereference() const {
				return _femMesh->dof(_dofIndex);
			}

			FEMMeshRectP1Periodic<Dim> const* _femMesh;
			signed_index_t _dofIndex;
		};
		typedef const_iterator iterator;

		iterator begin() const { return iterator(_femMesh, _beginIndex); }
		iterator end() const { return iterator(_femMesh, _endIndex); }
		const_iterator cbegin() const { return const_iterator(_femMesh, _beginIndex); }
		const_iterator cend() const { return const_iterator(_femMesh, _endIndex); }
	};
	class IterateDoFStrided {
	private:
		FEMMeshRectP1Periodic<Dim> const& _femMesh;
		const index_t _beginIndex;
		const index_t _endIndex;
		const index_t _stride;
	public:
		IterateDoFStrided(FEMMeshRectP1Periodic<Dim> const& femMesh, index_t beginIndex, index_t endIndex, index_t stride) :
			_femMesh(femMesh), _beginIndex(beginIndex), _endIndex(endIndex), _stride(stride) {}


		class const_iterator : public boost::iterator_facade<
		const_iterator
		, typename FEMMeshRectP1Periodic<Dim>::DoFT const
		, boost::random_access_traversal_tag
		, typename FEMMeshRectP1Periodic<Dim>::DoFT const
		, signed_index_t> {
		public:
			//const_iterator() : _element(0), _elementDoFIndex(0) {}
			//explicit const_iterator(Mesh_2DRect const& mesh) : _mesh(mesh), _dofIndex(0) {}
			const_iterator(FEMMeshRectP1Periodic<Dim> const& femMesh, index_t dofIndex, index_t stride) : _femMesh(&femMesh), _dofBaseIndex(dofIndex/stride), _stride(stride) {}

		private:
			friend class boost::iterator_core_access;

			void increment() { ++_dofBaseIndex; }
			void decrement() { --_dofBaseIndex; }
			void advance(signed_index_t n) { _dofBaseIndex+=n; }
			signed_index_t distance_to(const_iterator const& other) const { return other._dofBaseIndex - this->_dofBaseIndex; }

			bool equal(const_iterator const& other) const {	return this->_dofBaseIndex == other._dofBaseIndex; }

			typename FEMMeshRectP1Periodic<Dim>::DoFT const dereference() const {
				return _femMesh->dof(_dofBaseIndex*_stride);
			}

			FEMMeshRectP1Periodic<Dim> const* _femMesh;
			signed_index_t _dofBaseIndex;
			signed_index_t _stride;
		};
		typedef const_iterator iterator;

		iterator begin() const { return iterator(_femMesh, _beginIndex, _stride); }
		iterator end() const { return iterator(_femMesh, _endIndex, _stride); }
		const_iterator cbegin() const { return const_iterator(_femMesh, _beginIndex, _stride); }
		const_iterator cend() const { return const_iterator(_femMesh, _endIndex, _stride); }
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
static inline typename FEMMeshRectP1Periodic<Dim>::IterateMeshElements Elements(FEMMeshRectP1Periodic<Dim> const& femMesh) {
	return typename FEMMeshRectP1Periodic<Dim>::IterateMeshElements(femMesh, 0, femMesh.nElements());
}
template <unsigned int Dim>
static inline typename FEMMeshRectP1Periodic<Dim>::IterateMeshElements LocalElements(FEMMeshRectP1Periodic<Dim> const& femMesh) {
	return typename FEMMeshRectP1Periodic<Dim>::IterateMeshElements(femMesh, 0, femMesh.nLocalElements());
}
template <unsigned int Dim>
static inline typename FEMMeshRectP1Periodic<Dim>::IterateMeshElements NeighborElements(FEMMeshRectP1Periodic<Dim> const& femMesh) {
	return typename FEMMeshRectP1Periodic<Dim>::IterateMeshElements(femMesh, femMesh.nLocalElements(), femMesh.nElements());
}

template <unsigned int Dim>
static inline UIntRange ElementDoFIndices(FEMMeshRectP1Periodic<Dim> const& femMesh) {
	return UIntRange(0, FEMMeshRectP1Periodic<Dim>::nDoFPerElement());
}
//template <typename ElemT>
//static inline typename ElemT::IterateElementDoF ElementDoF(ElemT const& element) {
//	return typename ElemT::IterateElementDoF(element);
//}
//template <typename ElemT>
//static inline typename ElemT::IterateElementDomainBoundaryHypersurfaces DomainBoundaryHypersurfaces(ElemT const& element) {
//	return typename ElemT::IterateElementDomainBoundaryHypersurfaces(element);
//}
//template <typename HSurfaceT>
//static inline typename HSurfaceT::IterateHypersurfaceDoF HypersurfaceDoF(HSurfaceT const& hypersurface) {
//	return typename HSurfaceT::IterateHypersurfaceDoF(hypersurface);
//}

template <unsigned int Dim>
static inline typename FEMMeshRectP1Periodic<Dim>::IterateDoF DoF(FEMMeshRectP1Periodic<Dim> const& femMesh) {
	return typename FEMMeshRectP1Periodic<Dim>::IterateDoF(femMesh, 0, femMesh.nDoF());
}
template <unsigned int Dim>
static inline typename FEMMeshRectP1Periodic<Dim>::IterateDoF LocalDoF(FEMMeshRectP1Periodic<Dim> const& femMesh) {
	return typename FEMMeshRectP1Periodic<Dim>::IterateDoF(femMesh, 0, femMesh.nLocalDoF());
}
template <unsigned int Dim>
static inline typename FEMMeshRectP1Periodic<Dim>::IterateDoF GhostDoF(FEMMeshRectP1Periodic<Dim> const& femMesh) {
	return typename FEMMeshRectP1Periodic<Dim>::IterateDoF(femMesh, femMesh.nLocalDoF(), femMesh.nDoF());
}

template <unsigned int Dim>
static inline typename FEMMeshRectP1Periodic<Dim>::IterateDoFStrided BaseLevelDoF(FEMMeshRectP1Periodic<Dim> const& femMesh) {
	return typename FEMMeshRectP1Periodic<Dim>::IterateDoFStrided(femMesh, 0, femMesh.nDoF(), femMesh.dimensions[Dim-1].nNodes);
}




#endif /* MESH_PERIODIC_HPP_ */
