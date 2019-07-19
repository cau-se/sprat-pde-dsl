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

#ifndef PROTO_MATRIX_HPP_
#define PROTO_MATRIX_HPP_

template<typename = boost::proto::is_proto_expr>
struct CSRMatrix_ : ProtoCSRMatrix
{
	CSRMatrix_() {}
	CSRMatrix_(LILMatrix const& lilMatrix) {
		boost::proto::value(*this).reset(lilMatrix);
	}

	index_t nEntries() const {
		return boost::proto::value(*this).nEntries();
	}
	index_t nRows() const {
		return boost::proto::value(*this).nRows();
	}

	void print() const {
		boost::proto::value(*this).print();
	}

	template <typename DestVecT, typename SrcVecT>
	DestVecT & multiplyAdd(DestVecT & dest, SrcVecT const& src) const {
		return boost::proto::value(*this).multiplyAdd(dest, src);
	}
	template <typename DestVecT, typename SrcVecT>
	DestVecT & multiplySubtract(DestVecT & dest, SrcVecT const& src) const {
		return boost::proto::value(*this).multiplySubtract(dest, src);
	}
};

template<typename = boost::proto::is_proto_expr>
struct Matrix_ : ProtoMatrix
{
	Matrix_() {}
	Matrix_(index_t n) {
		boost::proto::value(*this).reset(n, n);
	}
	Matrix_(index_t nRows, index_t nCols) {
		boost::proto::value(*this).reset(nRows, nCols);
	}

	index_t nEntries() const {
		return boost::proto::value(*this).nEntries();
	}
	index_t nRows() const {
		return boost::proto::value(*this).nRows();
	}

	index_t nCols() const {
		return boost::proto::value(*this).nCols();
	}

	real get(index_t i, index_t j) const {
		return boost::proto::value(*this).get(i, j);
	}

	void set(index_t i, index_t j, real v) {
		boost::proto::value(*this).set(i, j, v);
	}

	template <typename DestVecT, typename SrcVecT>
	DestVecT & multiplyAdd(DestVecT & dest, SrcVecT const& src) const {
		return boost::proto::value(*this).multiplyAdd(dest, src);
	}
	template <typename DestVecT, typename SrcVecT>
	DestVecT & multiplySubtract(DestVecT & dest, SrcVecT const& src) const {
		return boost::proto::value(*this).multiplySubtract(dest, src);
	}
};


template<typename = boost::proto::is_proto_expr>
struct ElementMatrix_ : ProtoElementMatrix
{
	ElementMatrix_() {
		boost::proto::value(*this).reset();
	}
	ElementMatrix_(ElementT const& element) {
		boost::proto::value(*this).reset(element);
	}

	static constexpr index_t dimension() {
		return nDoFPerElement() * nDoFPerElement();
	}

	real * getValues() const {
		return boost::proto::value(*this).getValues();
	}

	template <typename ElementT>
	void setElement(ElementT const& element) {
		boost::proto::value(*this).setElement(element);
	}

	real& operator()(uint i, uint j) const {
		return boost::proto::value(*this).operator()(i,j);
	}

	real get(uint i, uint j) const {
		return boost::proto::value(*this).get(i, j);
	}
};

template<typename = boost::proto::is_proto_expr>
struct LightweightScaledElementMatrix_ : ProtoLightweightScaledElementMatrix
{
	LightweightScaledElementMatrix_() {}
	LightweightScaledElementMatrix_(ElementT const& element, real scaling, real * values) {
		boost::proto::value(*this).reset(element, scaling, values);
	}

	static constexpr index_t dimension() {
		return nDoFPerElement() * nDoFPerElement();
	}

	real operator()(uint i, uint j) const {
		return boost::proto::value(*this).operator()(i,j);
	}

	real get(uint i, uint j) const {
		return boost::proto::value(*this).get(i, j);
	}
};

template<typename = boost::proto::is_proto_expr>
struct ElementMatrixView_ : ProtoElementMatrixView
{
	ElementMatrixView_() {}
	ElementMatrixView_(ElementT const& element, real * values) {
		boost::proto::value(*this).reset(element, values);
	}

	static constexpr index_t dimension() {
		return nDoFPerElement() * nDoFPerElement();
	}

	real operator()(uint i, uint j) const {
		return boost::proto::value(*this).operator()(i,j);
	}
	real get(uint i, uint j) const {
		return boost::proto::value(*this).get(i, j);
	}

	void set(uint i, uint j, real value) const {
		boost::proto::value(*this).set(i, j, value);
	}
	void add(uint i, uint j, real value) const {
		boost::proto::value(*this).add(i, j, value);
	}
};

#endif //PROTO_MATRIX_HPP_
