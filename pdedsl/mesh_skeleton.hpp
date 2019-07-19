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

#ifndef MESH_SKELETON_HPP_
#define MESH_SKELETON_HPP_

#include <vector>
#include <cstdlib>
#include <limits>
#include <forward_list>
#include <boost/iterator/iterator_facade.hpp>
#include "config.hpp"
#include "util.hpp"
#include "interval.hpp"
#include "numa.hpp"


/*
 * Steps for creating a new mesh type:
 * 1. Copy this file and implement your own mesh type (instructions below)
 * 2. Include your mesh in "mesh.hpp" (see instructions there)
 */


// Dim is the number of dimensions of the domain
template <unsigned int Dim>
class FEMMeshSkeleton { // Choose your own name
public:
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
	constexpr static unsigned int twoToThePowerOfN(unsigned int n) {
		return (n>0 ? 2*twoToThePowerOfN(n-1) : 1);
	}

	struct ElementData {
		// Data to be stored for each element (e.g., its multi index and its neighbour elements)
	};
	struct DoFData {
		// Data to be stored for each DoF (e.g., its multi index)
	};

	numa_vector<ElementData> _elementData;
	numa_vector<DoFData> _dofData;

#ifdef SPRAT_BUILD_WITH_OPENMP
	// For foreachElementIndependently..
	// .. to which elements does DoF i belong?
	std::vector<std::forward_list<index_t>> _elementsOfDoF;
	// .. which sets of elements are there, whichs DoF can each be manipulated in parallel?
	std::vector<numa_vector<index_t>> _independentElementDecomposition;
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP

	
	void createElementsOfDoF(index_t numDoF, std::vector<ElementData> const& elemData) {
#ifdef SPRAT_BUILD_WITH_OPENMP
		_elementsOfDoF.resize(numDoF);

		for(index_t i=0; i<numDoF; ++i) {
			_elementsOfDoF[i].clear();
		}

		for(index_t i=0; i<elemData.size(); ++i) {
			for(uint k=0; k<nDoFPerElement(); ++k) {
				// Adapt the following to your ElementData:
				//_elementsOfDoF[elemData[i].dofIndices[k]].push_front(i);
			}
		}
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP
	}

	void createMaximalIndependentElementSets(std::vector<ElementData> const& elemData) {
#ifdef SPRAT_BUILD_WITH_OPENMP
		// Populate _independentElementDecomposition:
		/*
		 * Each element of _independentElementDecomposition is a set of 
		 * element indices whichs DoF can be manipulted in parallel 
		 * without interference (because each DoF is only contained in 
		 * exactly one element of the set).
		 */
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP
	}

public:
	FEMMeshSkeleton(/* ... */) {
		// Implement your own constructor
		
		//createElementsOfDoF(...);
		//createMaximalIndependentElementSets(...);
	}


#ifdef SPRAT_BUILD_WITH_OPENMP
	std::vector<numa_vector<index_t>> const& getIndependentElementDecomposition() const {
		return _independentElementDecomposition;
	}
#endif //#ifdef SPRAT_BUILD_WITH_OPENMP

	index_t nElements() const {
		return _elementData.size();
	}
	index_t nDoF() const {
		return _dofData.size();
	}


	/*
	 * Parts of a mesh (elements, DoF).
	 * Each Part can again have subtypes (elements may have local DoF or
	 * hypersurfaces).
	 */
	
	class ElementT {
	public:
		constexpr static unsigned int nDoFPerElement() {
			return FEMMeshSkeleton<Dim>::nDoFPerElement();
		}
	private:
		FEMMeshSkeleton<Dim> const& _femMesh;
		const index_t _index;

	public:
		ElementT(FEMMeshSkeleton<Dim> const& femMesh, index_t index) : 
			_femMesh(femMesh), _index(index)
		{}

		operator index_t() const {
			return _index;
		}
		index_t index() const {
			return _index;
		}

		FEMMeshSkeleton<Dim> const& getMesh() const {
			return _femMesh;
		}
	};

	ElementT const element(index_t index) const {
		return ElementT(*this, index);
	}

	class DoFT {
	private:
		FEMMeshSkeleton<Dim> const& _femMesh;
		const index_t _index;

	public:
		DoFT(FEMMeshSkeleton<Dim> const& femMesh, index_t index) : 
			_femMesh(femMesh), _index(index)
		{}

		operator index_t() const {
			return _index;
		}
		index_t index() const {
			return _index;
		}
	};

	DoFT dof(index_t index) const {
		return DoFT(*this, index);
	}


	/*
	 * These container objects are returned by sets like
	 *     Elements(femMesh)
	 * and implement functions to retrieve iterators (e.g., begin(), 
	 * end()). The iterators in this example are implemented using 
	 * Boost.Iterator.
	 */
	
	class IterateMeshElements {
	private:
		FEMMeshSkeleton<Dim> const& _femMesh;
		const index_t _beginIndex;
		const index_t _endIndex;

	public:
		IterateMeshElements(FEMMeshSkeleton<Dim> const& femMesh, index_t beginIndex, index_t endIndex) :
			_femMesh(femMesh), _beginIndex(beginIndex), _endIndex(endIndex) 
		{}

		class const_iterator : public boost::iterator_facade<
			const_iterator
			, typename FEMMeshSkeleton<Dim>::ElementT const
			, boost::random_access_traversal_tag
			, typename FEMMeshSkeleton<Dim>::ElementT const
			, signed_index_t> {
		public:
			const_iterator() : _femMesh(0), _elementIndex(0) {}
			const_iterator(FEMMeshSkeleton<Dim> const& femMesh, index_t elementIndex) : 
				_femMesh(&femMesh), _elementIndex(elementIndex)
			{}

		private:
			friend class boost::iterator_core_access;

			void increment() { ++_elementIndex; }
			void decrement() { --_elementIndex; }
			void advance(signed_index_t n) { _elementIndex+=n; }
			signed_index_t distance_to(const_iterator const& other) const { 
				return  other._elementIndex - this->_elementIndex; 
			}

			bool equal(const_iterator const& other) const {	
				return this->_elementIndex == other._elementIndex;
			}

			typename FEMMeshSkeleton<Dim>::ElementT const dereference() const {
				return _femMesh->element(_elementIndex);
			}

			FEMMeshSkeleton<Dim> const * _femMesh;
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
		FEMMeshSkeleton<Dim> const& _femMesh;
		const index_t _beginIndex;
		const index_t _endIndex;
	public:
		IterateDoF(FEMMeshSkeleton<Dim> const& femMesh, index_t beginIndex, index_t endIndex) :
			_femMesh(femMesh), _beginIndex(beginIndex), _endIndex(endIndex) 
		{}

		class const_iterator : public boost::iterator_facade<
			const_iterator
			, typename FEMMeshSkeleton<Dim>::DoFT const
			, boost::random_access_traversal_tag
			, typename FEMMeshSkeleton<Dim>::DoFT const
			, signed_index_t> {
		public:
			const_iterator(FEMMeshSkeleton<Dim> const& femMesh, index_t dofIndex) : 
				_femMesh(&femMesh), _dofIndex(dofIndex)
			{}

		private:
			friend class boost::iterator_core_access;

			void increment() { ++_dofIndex; }
			void decrement() { --_dofIndex; }
			void advance(signed_index_t n) { _dofIndex+=n; }
			signed_index_t distance_to(const_iterator const& other) const { 
				return other._dofIndex - this->_dofIndex; 
			}

			bool equal(const_iterator const& other) const {	
				return this->_dofIndex == other._dofIndex; 
			}

			typename FEMMeshSkeleton<Dim>::DoFT const dereference() const {
				return _femMesh->dof(_dofIndex);
			}

			FEMMeshSkeleton<Dim> const* _femMesh;
			signed_index_t _dofIndex;
		};
		typedef const_iterator iterator;

		iterator begin() const { return iterator(_femMesh, _beginIndex); }
		iterator end() const { return iterator(_femMesh, _endIndex); }
		const_iterator cbegin() const { return const_iterator(_femMesh, _beginIndex); }
		const_iterator cend() const { return const_iterator(_femMesh, _endIndex); }
	};
};


/*
 * The sets that can be used with this mesh for iterating upon them.
 */

template <unsigned int Dim>
static inline typename FEMMeshSkeleton<Dim>::IterateMeshElements Elements(FEMMeshSkeleton<Dim> const& femMesh) {
	return typename FEMMeshSkeleton<Dim>::IterateMeshElements(femMesh, 0, femMesh.nElements());
}

template <unsigned int Dim>
static inline typename FEMMeshSkeleton<Dim>::IterateDoF DoF(FEMMeshSkeleton<Dim> const& femMesh) {
	return typename FEMMeshSkeleton<Dim>::IterateDoF(femMesh, 0, femMesh.nDoF());
}


#endif
