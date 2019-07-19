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

#ifndef PROTO_VECTOR_HPP_
#define PROTO_VECTOR_HPP_

// ************************************* Vector Functionality

#define ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(OP, COMPFUNCTION, FAIL, SUCCESS)	\
		bool operator OP(real rhs) const { \
			const index_t numElements = boost::proto::value(*this).localSize(); \
			bool answer = SUCCESS; \
			foreach_omp_index(i, numElements, shared(rhs,answer), { \
				if(COMPFUNCTION(boost::proto::value(*this).get(i), rhs)) { \
					answer = FAIL; \
				} \
			}) \
			return answer; \
		}


#define ARITHMETIC_ELEMENTVECTOR_REALALMOSTCOMPARISON(OP, COMPFUNCTION, FAIL, SUCCESS) \
		bool operator OP(real rhs) const { \
			for(index_t i=0; i<nDoFPerElement(); ++i) { \
				if(COMPFUNCTION(boost::proto::value(*this).get(i), rhs)) { \
					return FAIL; \
				} \
			} \
			return SUCCESS; \
		}



#define ARITHMETIC_VECTOR_REALCOMPARISON(OP) \
		bool operator OP(real rhs) const { \
			const index_t numElements = boost::proto::value(*this).localSize(); \
			bool answer = true; \
			foreach_omp_index(i, numElements, shared(rhs,answer), { \
				if(!(boost::proto::value(*this).get(i) OP rhs)) { \
					answer = false; \
				} \
			}) \
			return answer; \
		}


#define ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(OP) \
		bool operator OP(real rhs) const { \
			for(uint i=0; i<nDoFPerElement(); ++i) { \
				if(!(boost::proto::value(*this).get(i) OP rhs)) { \
					return false; \
				} \
			} \
			return true; \
		}


#define ARITHMETIC_VECTOR_REALASSIGNMENT(VECTYPE, OP) \
		real operator OP(real value) const { \
			const index_t numElements = boost::proto::value(*this).localSize(); \
			foreach_omp_index(i, numElements, shared(value), { \
				boost::proto::value(*this)[i] OP value; \
			}) \
			return value; \
		}



#define ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(VECTYPE, OP) \
		template<typename Expr> \
		VECTYPE const& operator OP(VectorExprWrapper<Expr> const& expr) const { \
			const index_t numElements = boost::proto::value(*this).localSize(); \
			foreach_omp_index(i, numElements, shared(expr), { \
				boost::proto::value(*this)[i] OP expr[i]; \
			}) \
			return *this; \
		}


#define ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(VECTYPE, OP) \
		template<typename Expr> \
		VECTYPE const& operator OP(ElemVectorExprWrapper<Expr> const& expr) const { \
			mvAssert("Error: Only one ElementVector per Vector assignment supported!", 1 == boost::proto::eval(expr, CountElementVectors())); \
			mvAssert("Error: Sums of ElementVector expressions not supported for Vector assignment!", boost::proto::eval(expr, ExpressionIsJustAFactor())); \
			for(uint i=0; i<nDoFPerElement(); ++i) { \
				const index_t globalIndex = boost::proto::eval(expr, Local2GlobalDoFIndex(i)); \
				const real opResult = boost::proto::eval(expr, VectorSubscriptWithLocalAndGlobalIndex(i, globalIndex)); \
				boost::proto::value(*this)[globalIndex] OP opResult; \
			} \
			return *this; \
		}

#define ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(VECTYPE, MATRIXOP) \
		template<typename Expr> \
		VECTYPE const& operator MATRIXOP(ElemMVExprWrapper<Expr> const& expr) const { \
			mvAssert("Error: Sums of ElementMatrix expressions not supported for Vector assignment!", boost::proto::eval(expr, ExpressionIsJustAFactor())); \
			index_t globalDoF[nDoFPerElement()]; \
			for(uint i=0; i<nDoFPerElement(); ++i) { \
				globalDoF[i] = boost::proto::eval(expr, Local2GlobalDoFIndex(i)); \
			} \
			for(uint i=0; i<nDoFPerElement(); ++i) { \
				if(i < boost::proto::value(*this).localSize()) { \
					for(uint j=0; j<nDoFPerElement(); ++j) { \
						const real opResult = boost::proto::eval(expr, CalculateElemMatrixContribution(i, j, globalDoF[j])); \
						boost::proto::value(*this)[globalDoF[i]] MATRIXOP opResult; \
					} \
				} \
			} \
			return *this; \
		}

// (OP,MATRIXSIGN) must be either (+=,0) or (-=,1)
#define ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(VECTYPE, OP, MATRIXSIGN) \
		template<typename Expr> \
		VECTYPE const& operator OP(MVExprWrapper<Expr> const& expr) const { \
			if(boost::proto::eval(expr, VectorValuesDoAppear())) { \
				foreach_omp_index(i, boost::proto::value(*this).localSize(), shared(expr), { \
					boost::proto::value(*this)[i] OP boost::proto::eval(expr, VectorSubscriptContext(i)); \
				}) \
			} \
			boost::proto::eval(expr, AddMatrixContribution<VECTYPE>(*this, MATRIXSIGN)); \
			return *this; \
		}

// OP is =
#define ARITHMETIC_VECTOR_MVEXPRASSIGNMENT_EQUAL(VECTYPE) \
		template<typename Expr> \
		VECTYPE const& operator =(MVExprWrapper<Expr> const& expr) const { \
			if(boost::proto::eval(expr, VectorValuesDoAppear())) { \
				foreach_omp_index(i, boost::proto::value(*this).localSize(), shared(expr), { \
					boost::proto::value(*this)[i] = boost::proto::eval(expr, VectorSubscriptContext(i)); \
				}) \
			} \
			else { \
				foreach_omp_index(i, boost::proto::value(*this).localSize(), shared(expr), { \
					boost::proto::value(*this)[i] = 0.0; \
				}) \
			} \
			boost::proto::eval(expr, AddMatrixContribution<VECTYPE>(*this, 0)); \
			return *this; \
		}



#define ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(VECTYPE, OP) \
		real operator OP(real value) const { \
			for(uint i=0; i<nDoFPerElement(); ++i) { \
				boost::proto::value(*this)[i] OP value; \
			} \
			return value; \
		}

#define ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(VECTYPE, OP) \
		template<typename Expr> \
		VECTYPE const& operator OP(VectorExprWrapper<Expr> const& expr) const { \
			for(uint i=0; i<nDoFPerElement(); ++i) { \
				index_t globalIndex = boost::proto::value(*this).globalDoFIndex(i); \
				boost::proto::value(*this)[i] OP boost::proto::eval(expr, VectorSubscriptWithLocalAndGlobalIndex(i, globalIndex)); \
			} \
			return *this; \
		}

#define ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(VECTYPE, OP) \
		template<typename Expr> \
		VECTYPE const& operator OP(ElemVectorExprWrapper<Expr> const& expr) const { \
			if(boost::proto::eval(expr, GlobalVectorsDoAppear())) { \
				for(uint i=0; i<nDoFPerElement(); ++i) { \
					index_t globalIndex = boost::proto::value(*this).globalDoFIndex(i); \
					boost::proto::value(*this)[i] OP boost::proto::eval(expr, VectorSubscriptWithLocalAndGlobalIndex(i, globalIndex)); \
				} \
			} \
			else { \
				for(uint i=0; i<nDoFPerElement(); ++i) { \
					boost::proto::value(*this)[i] OP boost::proto::eval(expr, VectorSubscriptContext(i)); \
				} \
			} \
			return *this; \
		}

// OP must be either =, += or -= but not *= or /=!
#define ARITHMETIC_ELEMENTVECTOR_ELEMMVEXPRASSIGNMENT(VECTYPE, OP, MATRIXOP) \
		template<typename Expr> \
		VECTYPE const& operator OP(ElemMVExprWrapper<Expr> const& expr) const { \
			index_t globalDoF[nDoFPerElement()]; \
			for(uint i=0; i<nDoFPerElement(); ++i) { \
				globalDoF[i] = boost::proto::value(*this).globalDoFIndex(i); \
			} \
			if(boost::proto::eval(expr, GlobalVectorsDoAppear())) { \
				for(uint i=0; i<nDoFPerElement(); ++i) { \
					boost::proto::value(*this)[i] OP boost::proto::eval(expr, VectorSubscriptWithLocalAndGlobalIndex(i, globalDoF[i])); \
				} \
			} \
			else { \
				for(uint i=0; i<nDoFPerElement(); ++i) { \
					boost::proto::value(*this)[i] OP boost::proto::eval(expr, VectorSubscriptContext(i)); \
				} \
			} \
			for(uint i=0; i<nDoFPerElement(); ++i) { \
				for(uint j=0; j<nDoFPerElement(); ++j) { \
					boost::proto::value(*this)[i] MATRIXOP boost::proto::eval(expr, CalculateElemMatrixContribution(i, j, globalDoF[j])); \
				} \
			} \
			return *this; \
		}


// ************************************** Proto Vector Types

template<typename = boost::proto::is_proto_expr>
struct Vector_ : ProtoVector
{
public:

	Vector_() {}
	explicit Vector_(index_t const n) {
		boost::proto::value(*this).reset(n);
	}

	void touch(real initValue = 0.0) {
		boost::proto::value(*this).touch(initValue);
	}
	void touch(std::function<void (real *, index_t)> touchFunction) {
		boost::proto::value(*this).touch(touchFunction);
	}

	void printToFile(std::string prefix, uint id, std::string suffix = ".txt") const {
		boost::proto::value(*this).printToFile(prefix, id, suffix);
	}

	void writeBinaryToFile(std::string prefix, uint id, std::string suffix = ".bin") const {
		boost::proto::value(*this).writeBinaryToFile(prefix, id, suffix);
	}

	void readBinaryFromFile(std::string prefix, uint id, std::string suffix = ".bin") {
		boost::proto::value(*this).readBinaryFromFile(prefix, id, suffix);
	}

	void print() const { boost::proto::value(*this).print(); }


#ifdef SPRAT_BUILD_WITH_MPI
	std::vector<MPI_Request> collectFromSlavesAsync(std::string const& name, ParallelExecutionEnvironment const& pEE, std::vector<index_t> const& dofStartIndexForProcess) const {
		return boost::proto::value(*this).collectFromSlavesAsync(name, pEE, dofStartIndexForProcess);
	}
	std::vector<MPI_Request> distributeToSlavesAsync(std::string const& name, ParallelExecutionEnvironment const& pEE, std::vector<index_t> const& dofStartIndexForProcess) const {
		return boost::proto::value(*this).distributeToSlavesAsync(name, pEE, dofStartIndexForProcess);
	}
#endif


	void reset(index_t const n) { boost::proto::value(*this).reset(n); }
	index_t size() const { return boost::proto::value(*this).size(); }
	index_t localSize() const { return boost::proto::value(*this).localSize(); }
	//void swap(Vector_& other) { boost::proto::value(*this).swap(boost::proto::value(other)); }
	real * getValues() const { return boost::proto::value(*this).getValues(); }

	real& operator[](index_t const i) const { return boost::proto::value(*this)[i]; }
//	real operator[](index_t const i) const { return boost::proto::value(*this).get(i); }


	Vector_ const& operator=(Vector_ const& other) {
		boost::proto::value(*this).setTo(boost::proto::value(other));
		return *this;
	}
//	Vector_ const& operator=(VectorView const& other) {
//		boost::proto::value(*this).setTo(boost::proto::value(other));
//		return *this;
//	}

//	template<typename Expr>
//	Vector_ const& operator =(MVExprWrapper<Expr> const& expr) const {
//		return *this;
//	}

	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT_EQUAL(Vector_)
	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(Vector_, +=, 0)
	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(Vector_, -=, 1)

	ARITHMETIC_VECTOR_REALASSIGNMENT(Vector_, =)
	ARITHMETIC_VECTOR_REALASSIGNMENT(Vector_, +=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(Vector_, -=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(Vector_, *=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(Vector_, /=)

	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(Vector_, =)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(Vector_, +=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(Vector_, -=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(Vector_, *=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(Vector_, /=)

	//ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(Vector_, =);
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(Vector_, +=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(Vector_, -=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(Vector_, *=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(Vector_, /=)

	ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(Vector_, +=)
	ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(Vector_, -=)


	ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(==, !almostEqual, false, true)
	ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(!=, !almostEqual, true, false)
	ARITHMETIC_VECTOR_REALCOMPARISON(<)
	ARITHMETIC_VECTOR_REALCOMPARISON(>)
	ARITHMETIC_VECTOR_REALCOMPARISON(<=)
	ARITHMETIC_VECTOR_REALCOMPARISON(>=)


//	template<typename Expr>
//	VectorView_ const& operator=(ElemMVExprWrapper<Expr> expr) {
//		return Arithmetic_Vector_ElemMVExprAssignment(*this, expr);
//	}

	real norm() const {
		return boost::proto::value(*this).norm();
	}
	real average() const {
		return boost::proto::value(*this).average();
	}
	template<typename VecType>
	real dotProduct(VecType const& other) const {
		return boost::proto::value(*this).dotProduct(boost::proto::value(other));
	}
	template<typename VecType>
	real dotProductPartial(VecType const& other, index_t nEntries) const {
		return boost::proto::value(*this).dotProductPartial(boost::proto::value(other), nEntries);
	}
};

typedef Vector_<> Vector;



template<typename = boost::proto::is_proto_expr>
struct VectorView_ : ProtoVectorView
{
	VectorView_() {}
	VectorView_(index_t n, real * values) {
		boost::proto::value(*this).reset(n, values);
	}
	VectorView_(Vector const& other) {
		boost::proto::value(*this).reset(other.size(), other.getValues());
	}

	void touch(real initValue = 0.0) {
		boost::proto::value(*this).touch(initValue);
	}
	void touch(std::function<void (real *, index_t)> touchFunction) {
		boost::proto::value(*this).touch(touchFunction);
	}

	void printToFile(std::string prefix, uint id, std::string suffix = ".txt") {
		boost::proto::value(*this).printToFile(prefix, id, suffix);
	}

	void print() const { boost::proto::value(*this).print(); }

#ifdef SPRAT_BUILD_WITH_MPI
	std::vector<MPI_Request> collectFromSlavesAsync(std::string const& name, ParallelExecutionEnvironment const& pEE, std::vector<index_t> const& dofStartIndexForProcess) const {
		return boost::proto::value(*this).collectFromSlavesAsync(name, pEE, dofStartIndexForProcess);
	}
	std::vector<MPI_Request> distributeToSlavesAsync(std::string const& name, ParallelExecutionEnvironment const& pEE, std::vector<index_t> const& dofStartIndexForProcess) const {
		return boost::proto::value(*this).distributeToSlavesAsync(name, pEE, dofStartIndexForProcess);
	}
#endif

	void reset(index_t const n, real * values) { boost::proto::value(*this).reset(n, values); }
	index_t size() const { return boost::proto::value(*this).size(); }
	//void swap(Vector_& other) { boost::proto::value(*this).swap(boost::proto::value(other)); }
	real * getValues() const { return boost::proto::value(*this).getValues(); }

	real& operator[](index_t const i) const { return boost::proto::value(*this)[i]; }
	//real operator[](index_t const i) const { return boost::proto::value(*this).get(i); }

	//VectorView_ const& operator=(VectorView_ const& other) {
	//	boost::proto::value(*this).operator=(boost::proto::value(other));
	//	return *this;
	//}

	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT_EQUAL(VectorView_)
	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(VectorView_, +=, 0)
	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(VectorView_, -=, 1)

	ARITHMETIC_VECTOR_REALASSIGNMENT(VectorView_, =)
	ARITHMETIC_VECTOR_REALASSIGNMENT(VectorView_, +=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(VectorView_, -=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(VectorView_, *=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(VectorView_, /=)

	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(VectorView_, =)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(VectorView_, +=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(VectorView_, -=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(VectorView_, *=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(VectorView_, /=)

	//ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(VectorView_, =)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(VectorView_, +=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(VectorView_, -=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(VectorView_, *=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(VectorView_, /=)

	ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(VectorView_, +=)
	ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(VectorView_, -=)


	ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(==, !almostEqual, false, true)
	ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(!=, !almostEqual, true, false)
	ARITHMETIC_VECTOR_REALCOMPARISON(<)
	ARITHMETIC_VECTOR_REALCOMPARISON(>)
	ARITHMETIC_VECTOR_REALCOMPARISON(<=)
	ARITHMETIC_VECTOR_REALCOMPARISON(>=)

	//template<typename Expr>
	//VectorView_ const& operator=(ElemMVExprWrapper<Expr> expr) {
	//	return Arithmetic_Vector_ElemMVExprAssignment(*this, expr);
	//}

	real norm() const {
		return boost::proto::value(*this).norm();
	}
	real average() const {
		return boost::proto::value(*this).average();
	}
	template<typename VecType>
	real dotProduct(VecType const& other) const {
		return boost::proto::value(*this).dotProduct(boost::proto::value(other));
	}
	template<typename VecType>
	real dotProductPartial(VecType const& other, index_t nEntries) const {
		return boost::proto::value(*this).dotProductPartial(boost::proto::value(other), nEntries);
	}
};

typedef VectorView_<> VectorView;



template<typename = boost::proto::is_proto_expr>
struct DistributedVector_ : ProtoDistributedVector
{
	DistributedVector_() {}
	DistributedVector_(std::string name, ParallelExecutionEnvironment& pEE, MPCommunicationRegistry& comRegistry) {
		boost::proto::value(*this).reset(name, pEE, comRegistry);
	}

	void reset(std::string name, ParallelExecutionEnvironment& pEE, MPCommunicationRegistry& comRegistry) {
		boost::proto::value(*this).reset(name, pEE, comRegistry);
	}

	void touch(real initValue = 0.0) {
		boost::proto::value(*this).touch(initValue);
	}
	void touch(std::function<void (real *, index_t)> touchFunction) {
		boost::proto::value(*this).touch(touchFunction);
	}
	void touchLocal(real initValue = 0.0) {
		boost::proto::value(*this).touchLocal(initValue);
	}
	void touchLocal(std::function<void (real *, index_t)> touchFunction) {
		boost::proto::value(*this).touchLocal(touchFunction);
	}

	VectorView globalAsLocalVector() const {
		return VectorView(boost::proto::value(*this).size(), boost::proto::value(*this).getValues());
	}

	void print() const { boost::proto::value(*this).print(); }

	index_t size() const { return boost::proto::value(*this).size(); }
	index_t localSize() const { return boost::proto::value(*this).localSize(); }
	//void swap(Vector_& other) { boost::proto::value(*this).swap(boost::proto::value(other)); }
	real * getValues() const { return boost::proto::value(*this).getValues(); }
	int getTag() const { return boost::proto::value(*this).getTag(); }

	real& operator[](index_t const i) const { return boost::proto::value(*this)[i]; }
	//real operator[](index_t const i) const { return boost::proto::value(*this).get(i); }


	void collectFromMasterAsync() {
		boost::proto::value(*this).collectFromMasterAsync();
	}
	void distributeToMasterAsync() {
		boost::proto::value(*this).distributeToMasterAsync();
	}
	void finishPendingSends() {
		boost::proto::value(*this).finishPendingSends();
	}
	void finishPendingRecvs() {
		boost::proto::value(*this).finishPendingRecvs();
	}
	void finishPendingCommunication() {
		boost::proto::value(*this).finishPendingCommunication();
	}
	void exchangeDataAsync() {
		boost::proto::value(*this).exchangeDataAsync();
	}
	void exchangeDataBlocked() {
		boost::proto::value(*this).exchangeDataBlocked();
	}

	DistributedVector_ const& operator=(DistributedVector_ const& other) {
		boost::proto::value(*this).setTo(boost::proto::value(other));
		return *this;
	}
	//DistributedVector_ const& operator=(DistributedVectorView const& other) {
	//	boost::proto::value(*this).setTo(boost::proto::value(other));
	//	return *this;
	//}

	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT_EQUAL(DistributedVector_)
	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(DistributedVector_, +=, 0)
	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(DistributedVector_, -=, 1)

	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVector_, =)
	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVector_, +=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVector_, -=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVector_, *=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVector_, /=)

	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVector_, =)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVector_, +=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVector_, -=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVector_, *=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVector_, /=)

	//ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVector_, =)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVector_, +=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVector_, -=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVector_, *=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVector_, /=)

	ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(DistributedVector_, +=)
	ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(DistributedVector_, -=)


	ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(==, !almostEqual, false, true)
	ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(!=, !almostEqual, true, false)
	ARITHMETIC_VECTOR_REALCOMPARISON(<)
	ARITHMETIC_VECTOR_REALCOMPARISON(>)
	ARITHMETIC_VECTOR_REALCOMPARISON(<=)
	ARITHMETIC_VECTOR_REALCOMPARISON(>=)


	//template<typename Expr>
	//VectorView_ const& operator=(ElemMVExprWrapper<Expr> expr) {
	//	return Arithmetic_Vector_ElemMVExprAssignment(*this, expr);
	//}

	real norm() const {
		return boost::proto::value(*this).norm();
	}
	real average() const {
		return boost::proto::value(*this).average();
	}
	template<typename VecType>
	real dotProduct(VecType const& other) const {
		return boost::proto::value(*this).dotProduct(boost::proto::value(other));
	}
};

typedef DistributedVector_<> DistributedVector;



template<typename = boost::proto::is_proto_expr>
struct DistributedVectorView_ : ProtoDistributedVectorView
{
	DistributedVectorView_() {}
	DistributedVectorView_(DistributedVector const& other) {
		boost::proto::value(*this).reset(boost::proto::value(other));
	}

	void touch(real initValue = 0.0) {
		boost::proto::value(*this).touch(initValue);
	}
	void touch(std::function<void (real *, index_t)> touchFunction) {
		boost::proto::value(*this).touch(touchFunction);
	}
	void touchLocal(real initValue = 0.0) {
		boost::proto::value(*this).touchLocal(initValue);
	}
	void touchLocal(std::function<void (real *, index_t)> touchFunction) {
		boost::proto::value(*this).touchLocal(touchFunction);
	}

	VectorView globalAsLocalVector() const {
		return VectorView(boost::proto::value(*this).size(), boost::proto::value(*this).getValues());
	}

	index_t size() const { return boost::proto::value(*this).size(); }
	index_t localSize() const { return boost::proto::value(*this).localSize(); }
	//void swap(Vector_& other) { boost::proto::value(*this).swap(boost::proto::value(other)); }
	real * getValues() const { return boost::proto::value(*this).getValues(); }
	int getTag() const { return boost::proto::value(*this).getTag(); }

	real& operator[](index_t const i) const { return boost::proto::value(*this)[i]; }
	//real operator[](index_t const i) const { return boost::proto::value(*this).get(i); }

	void collectFromMasterAsync() {
		boost::proto::value(*this).collectFromMasterAsync();
	}
	void distributeToMasterAsync() {
		boost::proto::value(*this).distributeToMasterAsync();
	}
	void finishPendingSends() {
		boost::proto::value(*this).finishPendingSends();
	}
	void finishPendingRecvs() {
		boost::proto::value(*this).finishPendingRecvs();
	}
	void finishPendingCommunication() {
		boost::proto::value(*this).finishPendingCommunication();
	}
	void exchangeDataAsync() {
		boost::proto::value(*this).exchangeDataAsync();
	}
	void exchangeDataBlocked() {
		boost::proto::value(*this).exchangeDataBlocked();
	}

	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT_EQUAL(DistributedVectorView_)
	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(DistributedVectorView_, +=, 0)
	ARITHMETIC_VECTOR_MVEXPRASSIGNMENT(DistributedVectorView_, -=, 1)

	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVectorView_, =)
	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVectorView_, +=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVectorView_, -=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVectorView_, *=)
	ARITHMETIC_VECTOR_REALASSIGNMENT(DistributedVectorView_, /=)

	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVectorView_, =)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVectorView_, +=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVectorView_, -=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVectorView_, *=)
	ARITHMETIC_VECTOR_VECEXPRASSIGNMENT(DistributedVectorView_, /=)

	//ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVectorView_, =)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVectorView_, +=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVectorView_, -=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVectorView_, *=)
	ARITHMETIC_VECTOR_ELEMVECEXPRASSIGNMENT(DistributedVectorView_, /=)

	ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(DistributedVectorView_, +=)
	ARITHMETIC_VECTOR_ELEMMVEXPRASSIGNMENT(DistributedVectorView_, -=)


	ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(==, !almostEqual, false, true)
	ARITHMETIC_VECTOR_REALALMOSTCOMPARISON(!=, !almostEqual, true, false)
	ARITHMETIC_VECTOR_REALCOMPARISON(<)
	ARITHMETIC_VECTOR_REALCOMPARISON(>)
	ARITHMETIC_VECTOR_REALCOMPARISON(<=)
	ARITHMETIC_VECTOR_REALCOMPARISON(>=)

	//template<typename Expr>
	//VectorView_ const& operator=(ElemMVExprWrapper<Expr> expr) {
	//	return Arithmetic_Vector_ElemMVExprAssignment(*this, expr);
	//}

	real norm() const {
		return boost::proto::value(*this).norm();
	}
	real average() const {
		return boost::proto::value(*this).average();
	}
	template<typename VecType>
	real dotProduct(VecType const& other) const {
		return boost::proto::value(*this).dotProduct(boost::proto::value(other));
	}
};

typedef DistributedVectorView_<> DistributedVectorView;



template<typename = boost::proto::is_proto_expr>
struct ElementVector_ : ProtoElementVector
{
	ElementVector_() {
		boost::proto::value(*this).reset();
	}
	ElementVector_(ElementT const& element) {
		boost::proto::value(*this).reset(element);
	}
	ElementVector_(ElementVector_ const& other) {
		boost::proto::value(*this).reset(boost::proto::value(other));
	}

	void reset(ElementT const& element) {
		boost::proto::value(*this).reset(element);
	}
	static constexpr uint size() {
		return nDoFPerElement();
	}
	//void swap(Vector_& other) { boost::proto::value(*this).swap(boost::proto::value(other)); }
	real * getValues() const { return boost::proto::value(*this).getValues(); }

	real& operator[](uint const i) const { return boost::proto::value(*this)[i]; }
	real get(uint const i) const { return boost::proto::value(*this).get(i); }
	//real operator[](index_t const i) const { return boost::proto::value(*this).get(i); }

	ElementVector_ const& operator=(ElementVector_ const& other) {
		boost::proto::value(*this) = boost::proto::value(other);
		return *this;
	}

	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVector_, =)
	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVector_, +=)
	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVector_, -=)
	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVector_, *=)
	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVector_, /=)

	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVector_, =)
	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVector_, +=)
	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVector_, -=)
	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVector_, *=)
	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVector_, /=)

	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVector_, =)
	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVector_, +=)
	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVector_, -=)
	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVector_, *=)
	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVector_, /=)

	ARITHMETIC_ELEMENTVECTOR_ELEMMVEXPRASSIGNMENT(ElementVector_, =, +=)
	ARITHMETIC_ELEMENTVECTOR_ELEMMVEXPRASSIGNMENT(ElementVector_, +=, +=)
	ARITHMETIC_ELEMENTVECTOR_ELEMMVEXPRASSIGNMENT(ElementVector_, -=, -=)

	ARITHMETIC_ELEMENTVECTOR_REALALMOSTCOMPARISON(==, !almostEqual, false, true)
	ARITHMETIC_ELEMENTVECTOR_REALALMOSTCOMPARISON(!=, !almostEqual, true, false)
	ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(<)
	ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(>)
	ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(<=)
	ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(>=)

	real norm() const {
		return boost::proto::value(*this).norm();
	}
	template<typename VecType>
	real dotProduct(VecType const& other) const {
		return boost::proto::value(*this).dotProduct(boost::proto::value(other));
	}

	template <typename ElementT>
	void setElement(ElementT const& element) {
		boost::proto::value(*this).setElement(element);
	}
};

typedef ElementVector_<> ElementVector;



template<typename = boost::proto::is_proto_expr>
struct ElementVectorView_ : ProtoElementVectorView
{
	ElementVectorView_() {}
	ElementVectorView_(ElementT const& element, real * values) {
		boost::proto::value(*this).reset(element, values);
	}
	ElementVectorView_(ElementVectorView_ const& other) {
		boost::proto::value(*this).reset(boost::proto::value(other));
	}
	ElementVectorView_(ElementVector_<> const& other) {
		boost::proto::value(*this).reset(boost::proto::value(other));
	}

	void reset(ElementT const& element, real * values) {
		boost::proto::value(*this).reset(element, values);
	}
	static constexpr uint size() {
		return nDoFPerElement();
	}
	//void swap(Vector_& other) { boost::proto::value(*this).swap(boost::proto::value(other)); }
	real * getValues() const { return boost::proto::value(*this).getValues(); }

	real& operator[](uint const i) const { return boost::proto::value(*this)[i]; }
	real get(uint const i) const { return boost::proto::value(*this).get(i); }
	//real operator[](index_t const i) const { return boost::proto::value(*this).get(i); }

	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVectorView_, =)
	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVectorView_, +=)
	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVectorView_, -=)
	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVectorView_, *=)
	ARITHMETIC_ELEMENTVECTOR_REALASSIGNMENT(ElementVectorView_, /=)

	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVectorView_, =)
	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVectorView_, +=)
	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVectorView_, -=)
	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVectorView_, *=)
	ARITHMETIC_ELEMENTVECTOR_VECEXPRASSIGNMENT(ElementVectorView_, /=)

	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVectorView_, =)
	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVectorView_, +=)
	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVectorView_, -=)
	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVectorView_, *=)
	ARITHMETIC_ELEMENTVECTOR_ELEMVECEXPRASSIGNMENT(ElementVectorView_, /=)

	ARITHMETIC_ELEMENTVECTOR_ELEMMVEXPRASSIGNMENT(ElementVectorView_, =, +=)
	ARITHMETIC_ELEMENTVECTOR_ELEMMVEXPRASSIGNMENT(ElementVectorView_, +=, +=)
	ARITHMETIC_ELEMENTVECTOR_ELEMMVEXPRASSIGNMENT(ElementVectorView_, -=, -=)

	ARITHMETIC_ELEMENTVECTOR_REALALMOSTCOMPARISON(==, !almostEqual, false, true)
	ARITHMETIC_ELEMENTVECTOR_REALALMOSTCOMPARISON(!=, !almostEqual, true, false)
	ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(<)
	ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(>)
	ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(<=)
	ARITHMETIC_ELEMENTVECTOR_REALCOMPARISON(>=)

	real norm() const {
		return boost::proto::value(*this).norm();
	}
	template<typename VecType>
	real dotProduct(VecType const& other) const {
		return boost::proto::value(*this).dotProduct(boost::proto::value(other));
	}
};

typedef ElementVectorView_<> ElementVectorView;

#endif //PROTO_VECTOR_HPP_
