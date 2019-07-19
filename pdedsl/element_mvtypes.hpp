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

#ifndef ELEMENT_MVTYPES_HPP_
#define ELEMENT_MVTYPES_HPP_

#define MATRIX_INDEX(i, j, n) ((i)*(n) + (j))



class BareElementVectorBase {
protected:
	real * _values;
	index_t const * _globalDoF;

public:
	BareElementVectorBase() {}
	BareElementVectorBase(real * values) : _values(values) {}
	BareElementVectorBase(real * values, index_t const * globalDoF) : _values(values), _globalDoF(globalDoF) {}

	static constexpr uint size() {
		return ElementT::nDoFPerElement();
	}

	real * getValues() const {
		return _values;
	}
	index_t const * getGlobalDoFIndices() const {
		return _globalDoF;
	}

	index_t globalDoFIndex(uint localIndex) const {
		return _globalDoF[localIndex];
	}

	real& operator[](uint index) const {
		return _values[index];
	}

	real get(uint index) const {
		return _values[index];
	}
};

class BareElementVector : public BareElementVectorBase {
public:
	BareElementVector() : BareElementVectorBase(0) {}
	BareElementVector(ElementT const& element) : BareElementVectorBase(numa_new<real>(ElementT::nDoFPerElement()), element.globalDoFIndices()) {}
	BareElementVector(BareElementVector const& other) : BareElementVectorBase(numa_new<real>(ElementT::nDoFPerElement()), other._globalDoF) {}
	~BareElementVector() {
		numa_delete(_values);
	}

	void reset() {
		numa_delete(_values);
		_values = numa_new<real>(nDoFPerElement());
		_globalDoF = 0;
	}

	void reset(BareElementVector const& other) {
		numa_delete(_values);
		_values = numa_new<real>(nDoFPerElement());
		for(uint i=0; i<nDoFPerElement(); ++i) {
			_values[i] = other.get(i);
		}
		_globalDoF = other._globalDoF;
	}

	template <typename ElementT>
	void reset(ElementT const& element) {
		numa_delete(_values);
		_values = numa_new<real>(nDoFPerElement());
		_globalDoF = element.globalDoFIndices();
	}

	template <typename ElementT>
	void setElement(ElementT const& element) {
		_globalDoF = element.globalDoFIndices();
	}

	BareElementVector const& operator=(BareElementVector const& other) {
		for(uint i=0; i<nDoFPerElement(); ++i) {
			_values[i] = other.get(i);
		}
		return *this;
	}
};

class BareElementVectorView : public BareElementVectorBase {
public:
	BareElementVectorView() : BareElementVectorBase() {}
	BareElementVectorView(ElementT const& element, real * values) : BareElementVectorBase(values, element.globalDoFIndices()) {}

	void reset(ElementT const& element, real * values) {
		_values = values;
		_globalDoF = element.globalDoFIndices();
	}
	void reset(BareElementVectorView const& other) {
		_values = other.getValues();
		_globalDoF = other.getGlobalDoFIndices();
	}
	void reset(BareElementVector const& other) {
		_values = other.getValues();
		_globalDoF = other.getGlobalDoFIndices();
	}

	BareElementVectorView const& operator=(BareElementVectorView const& other) {
		for(index_t i=0; i<nDoFPerElement(); ++i) {
			_values[i] = other.get(i);
		}
		return *this;
	}
};




class BareElementMatrix {
private:
	real * _values; // row-major
	index_t const * _globalDoF;

public:
	BareElementMatrix() :
		_values(0) {}
	BareElementMatrix(ElementT const& element) :
		_values(numa_new<real>(nDoFPerElement() * nDoFPerElement())),
		_globalDoF(element.globalDoFIndices()) {}
	~BareElementMatrix() {
		numa_delete(_values);
	}

	void reset() {
		numa_delete(_values);
		_values = numa_new<real>(nDoFPerElement() * nDoFPerElement());
		_globalDoF = 0;
	}

	template <typename ElementT>
	void reset(ElementT const& element) {
		numa_delete(_values);
		_values = numa_new<real>(nDoFPerElement() * nDoFPerElement());
		_globalDoF = element.globalDoFIndices();
	}

	static constexpr index_t dimension() {
		return nDoFPerElement() * nDoFPerElement();
	}

	index_t globalDoFIndex(uint localIndex) const {
		return _globalDoF[localIndex];
	}

	real * getValues() const {
		return _values;
	}

	template <typename ElementT>
	void setElement(ElementT const& element) {
		_globalDoF = element.globalDoFIndices();
	}

	real& operator()(uint i, uint j) const {
		return _values[MATRIX_INDEX(i, j, nDoFPerElement())];
	}

	real get(uint i, uint j) const {
		return _values[MATRIX_INDEX(i, j, nDoFPerElement())];
	}
};


class BareElementMatrixView {
private:
	real * _values; // row-major
	index_t const * _globalDoF;

public:
	BareElementMatrixView() {}
	BareElementMatrixView(ElementT const& element, real * values) :
		_values(values), _globalDoF(element.globalDoFIndices())  {}

	void reset(ElementT const& element, real * values) {
		_values =  values;
		_globalDoF = element.globalDoFIndices();
	}

	index_t globalDoFIndex(uint localIndex) const {
		return _globalDoF[localIndex];
	}

	real operator()(uint i, uint j) const {
		return _values[MATRIX_INDEX(i, j, nDoFPerElement())];
	}
	real get(uint i, uint j) const {
		return _values[MATRIX_INDEX(i, j, nDoFPerElement())];
	}

	void set(uint i, uint j, real value) const {
		_values[MATRIX_INDEX(i, j, nDoFPerElement())] = value;
	}
	void add(uint i, uint j, real value) const {
		_values[MATRIX_INDEX(i, j, nDoFPerElement())] += value;
	}
};

class BareLightweightScaledElementMatrix {
private:
	real const * _values; // row-major
	real _scaling;
	index_t const * _globalDoF;

public:
	BareLightweightScaledElementMatrix() {}
	BareLightweightScaledElementMatrix(ElementT const& element, real scaling, real * values) :
		_values(values), _scaling(scaling), _globalDoF(element.globalDoFIndices())  {}

	void reset(ElementT const& element, real scaling, real * values) {
		_values =  values;
		_scaling = scaling;
		_globalDoF = element.globalDoFIndices();
	}

	index_t globalDoFIndex(uint localIndex) const {
		return _globalDoF[localIndex];
	}

	real operator()(uint i, uint j) const {
		return _scaling * _values[MATRIX_INDEX(i, j, nDoFPerElement())];
	}

	real get(uint i, uint j) const {
		return _scaling * _values[MATRIX_INDEX(i, j, nDoFPerElement())];
	}
};

#endif //ELEMENT_MVTYPES_HPP_
