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

#ifndef PROTO_CONTEXT_HPP_
#define PROTO_CONTEXT_HPP_

// Contexts
// 1. GlobalVectorsDoAppear
struct GlobalVectorsDoAppear : boost::proto::callable_context<
	  GlobalVectorsDoAppear const // derived context
	, boost::proto::null_context const  // fall-back context
> {
	typedef bool result_type;

	bool operator()(boost::proto::tag::terminal, BareVector const& v) const {
		return true;
	}
	bool operator()(boost::proto::tag::terminal, BareVectorView const& v) const {
		return true;
	}
#ifdef SPRAT_BUILD_WITH_MPI
	bool operator()(boost::proto::tag::terminal, BareDistributedVector const& v) const {
		return true;
	}
	bool operator()(boost::proto::tag::terminal, BareDistributedVectorView const& v) const {
		return true;
	}
#endif
	template <typename Terminal>
	bool operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
		return false;
	}


	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::plus, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) || boost::proto::eval(right, *this));
	}

	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::minus, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) || boost::proto::eval(right, *this));
	}

	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, ElementMatrix const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, LightweightScaledElementMatrix const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, ElementMatrixView const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, Matrix const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, CSRMatrix const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::multiplies, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) || boost::proto::eval(right, *this));
	}

	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::divides, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) || boost::proto::eval(right, *this));
	}

	template <typename Child>
	bool operator()(boost::proto::tag::negate, Child const& child) const {
		return boost::proto::eval(child, *this);
	}
};
// 1b. VectorValuesDoAppear
struct VectorValuesDoAppear : boost::proto::callable_context<
	  VectorValuesDoAppear const // derived context
	, boost::proto::null_context const  // fall-back context
> {
	typedef bool result_type;

	template <typename Terminal>
	bool operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
		return true;
	}


	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::plus, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) || boost::proto::eval(right, *this));
	}

	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::minus, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) || boost::proto::eval(right, *this));
	}

	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, ElementMatrix const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, LightweightScaledElementMatrix const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, ElementMatrixView const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, Matrix const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	real operator()(boost::proto::tag::multiplies, CSRMatrix const& left, Right const& right) const {
		return false;
	}
	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::multiplies, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) || boost::proto::eval(right, *this));
	}

	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::divides, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) || boost::proto::eval(right, *this));
	}

	template <typename Child>
	bool operator()(boost::proto::tag::negate, Child const& child) const {
		return boost::proto::eval(child, *this);
	}
};
// 2. VectorSubscriptWithLocalAndGlobalIndex
struct VectorSubscriptWithLocalAndGlobalIndex : boost::proto::callable_context< VectorSubscriptWithLocalAndGlobalIndex const >
{
private:
	uint const _localIndex;
	index_t const _globalIndex;

public:
	typedef real result_type;

	explicit VectorSubscriptWithLocalAndGlobalIndex(uint localIndex, index_t globalIndex) : _localIndex(localIndex), _globalIndex(globalIndex) {}

	real operator()(boost::proto::tag::terminal, BareVector const& v) const {
		return v.get(_globalIndex);
	}
	real operator()(boost::proto::tag::terminal, BareVectorView const& v) const {
		return v.get(_globalIndex);
	}
#ifdef SPRAT_BUILD_WITH_MPI
	real operator()(boost::proto::tag::terminal, BareDistributedVector const& v) const {
		return v.get(_globalIndex);
	}
	real operator()(boost::proto::tag::terminal, BareDistributedVectorView const& v) const {
		return v.get(_globalIndex);
	}
#endif
	real operator()(boost::proto::tag::terminal, BareElementVectorView const& v) const {
		return v.get(_localIndex);
	}
	real operator()(boost::proto::tag::terminal, BareElementVector const& v) const {
		return v.get(_localIndex);
	}
	real operator()(boost::proto::tag::terminal, BareElementMatrix const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, BareElementMatrixView const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, BareLightweightScaledElementMatrix const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, BareMatrix const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, BareCSRMatrix const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, real v) const {
		return v;
	}
	template <typename Terminal>
	real operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
		return terminal;
	}
};
// 3. VectorSubscript (index can be global or local)
struct VectorSubscriptContext : boost::proto::callable_context< VectorSubscriptContext const >
{
private:
	index_t const _index;

public:
	typedef real result_type;

	explicit VectorSubscriptContext(index_t i) : _index(i) {}

	real operator()(boost::proto::tag::terminal, BareVector const& v) const {
		return v.get(_index);
	}
	real operator()(boost::proto::tag::terminal, BareVectorView const& v) const {
		return v.get(_index);
	}
#ifdef SPRAT_BUILD_WITH_MPI
	real operator()(boost::proto::tag::terminal, BareDistributedVector const& v) const {
		return v.get(_index);
	}
	real operator()(boost::proto::tag::terminal, BareDistributedVectorView const& v) const {
		return v.get(_index);
	}
#endif
	real operator()(boost::proto::tag::terminal, BareElementVectorView const& v) const {
		return v.get(_index);
	}
	real operator()(boost::proto::tag::terminal, BareElementVector const& v) const {
		return v.get(_index);
	}
	real operator()(boost::proto::tag::terminal, BareElementMatrix const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, BareElementMatrixView const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, BareLightweightScaledElementMatrix const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, BareMatrix const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, BareCSRMatrix const& v) const {
		return 0.0;
	}
	real operator()(boost::proto::tag::terminal, real v) const {
		return v;
	}
	template <typename Terminal>
	real operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
		return terminal;
	}
};
// 4. VectorSubscriptWithFixedValueForElementVector (has a global index as well as the value for the element vector)
//struct VectorSubscriptWithFixedValueForElementVector : boost::proto::callable_context< VectorSubscriptWithFixedValueForElementVector const >
//{
//private:
//	index_t const _globalIndex;
//	real const _elementVecValue;
//
//public:
//	typedef real result_type;
//
//	explicit VectorSubscriptWithFixedValueForElementVector(index_t globalIndex, real elementVecValue) : _globalIndex(globalIndex), _elementVecValue(elementVecValue) {}
//
//	real operator()(boost::proto::tag::terminal, BareVector const& v) const {
//		return v.get(_globalIndex);
//	}
//	real operator()(boost::proto::tag::terminal, BareVectorView const& v) const {
//		return v.get(_globalIndex);
//	}
//#ifdef SPRAT_BUILD_WITH_MPI
//	real operator()(boost::proto::tag::terminal, BareDistributedVector const& v) const {
//		return v.get(_globalIndex);
//	}
//	real operator()(boost::proto::tag::terminal, BareDistributedVectorView const& v) const {
//		return v.get(_globalIndex);
//	}
//#endif
//	real operator()(boost::proto::tag::terminal, BareElementVectorView const& v) const {
//		return _elementVecValue;
//	}
//	real operator()(boost::proto::tag::terminal, BareElementVector const& v) const {
//		return _elementVecValue;
//	}
//};
// 5. CountElementVectors
struct CountElementVectors : boost::proto::callable_context< CountElementVectors const >
{
	typedef uint result_type;

	result_type operator()(boost::proto::tag::terminal, BareElementVectorView const& v) const {
		return 1;
	}
	result_type operator()(boost::proto::tag::terminal, BareElementVector const& v) const {
		return 1;
	}
	template <typename Terminal>
	result_type operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
		return 0;
	}


	template <typename Left, typename Right>
	result_type operator()(boost::proto::tag::plus, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) + boost::proto::eval(right, *this));
	}
	template <typename Left, typename Right>
	result_type operator()(boost::proto::tag::minus, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) + boost::proto::eval(right, *this));
	}
	template <typename Left, typename Right>
	result_type operator()(boost::proto::tag::multiplies, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) + boost::proto::eval(right, *this));
	}
	template <typename Left, typename Right>
	result_type operator()(boost::proto::tag::divides, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) + boost::proto::eval(right, *this));
	}

	template <typename Child>
	result_type operator()(boost::proto::tag::negate, Child const& child) const {
		return boost::proto::eval(child, *this);
	}
};

// 6. ExpressionIsJustAFactor
struct ExpressionIsJustAFactor : boost::proto::callable_context<
	  ExpressionIsJustAFactor const // derived context
	, boost::proto::null_context const  // fall-back context
> {
	typedef bool result_type;

	template <typename Terminal>
	bool operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
		return true;
	}


	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::plus, Left const& left, Right const& right) const {
		return false;
	}

	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::minus, Left const& left, Right const& right) const {
		return false;
	}

	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::multiplies, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) && boost::proto::eval(right, *this));
	}

	template <typename Left, typename Right>
	bool operator()(boost::proto::tag::divides, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) && boost::proto::eval(right, *this));
	}

	template <typename Child>
	bool operator()(boost::proto::tag::negate, Child const& child) const {
		return boost::proto::eval(child, *this);
	}
};
//7. Local2GlobalDoFIndex
struct Local2GlobalDoFIndex : boost::proto::callable_context< Local2GlobalDoFIndex const >
{
private:
	uint const _localIndex;

public:
	typedef index_t result_type;

	Local2GlobalDoFIndex(uint localIndex) : _localIndex(localIndex) {}

	index_t operator()(boost::proto::tag::terminal, BareElementVectorView const& v) const {
		return v.globalDoFIndex(_localIndex);
	}
	index_t operator()(boost::proto::tag::terminal, BareElementVector const& v) const {
		return v.globalDoFIndex(_localIndex);
	}
	index_t operator()(boost::proto::tag::terminal, BareElementMatrix const& M) const {
		return M.globalDoFIndex(_localIndex);
	}
	index_t operator()(boost::proto::tag::terminal, BareLightweightScaledElementMatrix const& M) const {
		return M.globalDoFIndex(_localIndex);
	}
	index_t operator()(boost::proto::tag::terminal, BareElementMatrixView const& M) const {
		return M.globalDoFIndex(_localIndex);
	}
	template <typename Terminal>
	index_t operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
		return 0;
	}


	template <typename Left, typename Right>
	index_t operator()(boost::proto::tag::plus, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) + boost::proto::eval(right, *this));
	}
	template <typename Left, typename Right>
	index_t operator()(boost::proto::tag::minus, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) + boost::proto::eval(right, *this));
	}
	template <typename Left, typename Right>
	index_t operator()(boost::proto::tag::multiplies, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) + boost::proto::eval(right, *this));
	}
	template <typename Left, typename Right>
	index_t operator()(boost::proto::tag::divides, Left const& left, Right const& right) const {
		return (boost::proto::eval(left, *this) + boost::proto::eval(right, *this));
	}

	template <typename Child>
	index_t operator()(boost::proto::tag::negate, Child const& child) const {
		return boost::proto::eval(child, *this);
	}
};


// 8. CalculateElemMatrixContribution
struct CalculateElemMatrixContribution : boost::proto::callable_context<CalculateElemMatrixContribution const> {
	typedef real result_type;

	uint i, j;
	index_t j_global;

	CalculateElemMatrixContribution(uint i_in, uint j_in, index_t j_global_in) : i(i_in), j(j_in), j_global(j_global_in) {}

	template <typename Terminal>
	real operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
		return 0.0;
	}

	template <typename Right>
	real operator()(boost::proto::tag::multiplies, ElementMatrix const& left, Right const& right) const
	{
		return boost::proto::value(left).get(i, j) * right[j_global];
	}

	template <typename Right>
	real operator()(boost::proto::tag::multiplies, LightweightScaledElementMatrix const& left, Right const& right) const
	{
		return boost::proto::value(left).get(i, j) * right[j_global];
	}

	template <typename Right>
	real operator()(boost::proto::tag::multiplies, ElementMatrixView const& left, Right const& right) const
	{
		return boost::proto::value(left).get(i, j) * right[j_global];
	}
};

//9. GetSomeElement
// Alternative to Local2GlobalDoFIndex; slower according to measurements
//struct GetSomeElement : boost::proto::callable_context< GetSomeElement const >
//{
//	typedef ElementT const * const result_type;
//
//	result_type operator()(boost::proto::tag::terminal, BareElementVectorView const& v) const {
//		return &(v.getElement());
//	}
//	result_type operator()(boost::proto::tag::terminal, BareElementVector const& v) const {
//		return &(v.getElement());
//	}
//	result_type operator()(boost::proto::tag::terminal, BareElementMatrix const& v) const {
//		return &(v.getElement());
//	}
//	result_type operator()(boost::proto::tag::terminal, BareLightweightScaledElementMatrix const& v) const {
//		return &(v.getElement());
//	}
//	template <typename Terminal>
//	result_type operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
//		return 0;
//	}
//
//
//	template <typename Left, typename Right>
//	result_type operator()(boost::proto::tag::plus, Left const& left, Right const& right) const {
//		result_type result = boost::proto::eval(left, *this);
//		return (result ? result : boost::proto::eval(right, *this));
//	}
//	template <typename Left, typename Right>
//	result_type operator()(boost::proto::tag::minus, Left const& left, Right const& right) const {
//		result_type result = boost::proto::eval(left, *this);
//		return (result ? result : boost::proto::eval(right, *this));
//	}
//	template <typename Left, typename Right>
//	result_type operator()(boost::proto::tag::multiplies, Left const& left, Right const& right) const {
//		result_type result = boost::proto::eval(left, *this);
//		return (result ? result : boost::proto::eval(right, *this));
//	}
//	template <typename Left, typename Right>
//	result_type operator()(boost::proto::tag::divides, Left const& left, Right const& right) const {
//		result_type result = boost::proto::eval(left, *this);
//		return (result ? result : boost::proto::eval(right, *this));
//	}
//
//	template <typename Child>
//	result_type operator()(boost::proto::tag::negate, Child const& child) const {
//		return boost::proto::eval(child, *this);
//	}
//};

// 10. AddMatrixContribution
template <typename VecT>
struct AddMatrixContribution : boost::proto::callable_context<
	AddMatrixContribution<VecT> const
	, boost::proto::null_context const  // fall-back context
> {
	typedef void result_type;

	VecT const& target;
	int sign; // Even: positive, odd: negative


	AddMatrixContribution(VecT const& target_in, int sign_in)
			: target(target_in), sign(sign_in) {}

	template <typename Terminal>
	void operator()(boost::proto::tag::terminal, Terminal const& terminal) const {
	}

	template <typename Left, typename Right>
	void operator()(boost::proto::tag::plus, Left const& left, Right const& right) const {
		boost::proto::eval(left, AddMatrixContribution<VecT>(target, sign));
		boost::proto::eval(right, AddMatrixContribution<VecT>(target, sign));
	}
	template <typename Left, typename Right>
	void operator()(boost::proto::tag::minus, Left const& left, Right const& right) const {
		boost::proto::eval(left, AddMatrixContribution<VecT>(target, sign));
		boost::proto::eval(right, AddMatrixContribution<VecT>(target, sign+1));
	}
	template <typename Right>
	void operator()(boost::proto::tag::multiplies, CSRMatrix const& left, Right const& right) const {
		if(sign%2 == 0) {
			left.multiplyAdd(target, right);
		}
		else {
			left.multiplySubtract(target, right);
		}
	}
	template <typename Right>
	void operator()(boost::proto::tag::multiplies, Matrix const& left, Right const& right) const {
		if(sign%2 == 0) {
			left.multiplyAdd(target, right);
		}
		else {
			left.multiplySubtract(target, right);
		}
	}
	template <typename Left, typename Right>
	void operator()(boost::proto::tag::multiplies, Left const& left, Right const& right) const {
	}
	template <typename Left, typename Right>
	void operator()(boost::proto::tag::divides, Left const& left, Right const& right) const {
	}

	template <typename Child>
	void operator()(boost::proto::tag::negate, Child const& child) const {
		sign += 1;
		boost::proto::eval(child, *this);
	}
};

#endif //PROTO_CONTEXT_HPP_
