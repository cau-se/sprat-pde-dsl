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

#ifndef ARITHMETIC_HPP_
#define ARITHMETIC_HPP_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <boost/proto/proto.hpp>
#ifdef SPRAT_BUILD_WITH_MPI
#include <mpi.h>
#endif
#include "config.hpp"
#include "util.hpp"
#include "pvector.hpp"
#include "sparse_matrices.hpp"
#include "mesh.hpp"
#include "mesh_traits.hpp"



#define mvAssert(MSG, COND) \
	assert((COND) && MSG);



// *************************************** Vectors

// See pvectors.hpp


// *************************************** Element Vectors & Matrices

#include "element_mvtypes.hpp"


// *************************************** Proto Grammars

struct VectorType
	: boost::proto::or_<
	    boost::proto::terminal< BareVector >
	  , boost::proto::terminal< BareVectorView >
	  , boost::proto::terminal< BareDistributedVector >
	  , boost::proto::terminal< BareDistributedVectorView >
> {};

struct ElementVectorType
	: boost::proto::or_<
	    boost::proto::terminal< BareElementVector >
	  , boost::proto::terminal< BareElementVectorView >
> {};

struct ElemetMatrixType
	: boost::proto::or_<
	  	  boost::proto::terminal< BareElementMatrix >
		, boost::proto::terminal< BareLightweightScaledElementMatrix >
		, boost::proto::terminal< BareElementMatrixView >
> {};

struct MatrixType
	: boost::proto::or_<
	  boost::proto::terminal< BareMatrix >
	, boost::proto::terminal< BareCSRMatrix >
//	, ElemetMatrixType
> {};

struct VectorGrammar
	: boost::proto::or_<
	    VectorType
	  , boost::proto::terminal< boost::proto::convertible_to<real> >
	  , boost::proto::plus< VectorGrammar, VectorGrammar >
	  , boost::proto::minus< VectorGrammar, VectorGrammar >
	  , boost::proto::multiplies< VectorGrammar, VectorGrammar >
	  , boost::proto::divides< VectorGrammar, VectorGrammar >
	  , boost::proto::negate< VectorGrammar >
> {};


struct ElementVectorGrammar
	: boost::proto::or_<
	    VectorType
	  , ElementVectorType
	  , boost::proto::terminal< boost::proto::convertible_to<real> >
	  , boost::proto::plus< ElementVectorGrammar, ElementVectorGrammar >
	  , boost::proto::minus< ElementVectorGrammar, ElementVectorGrammar >
	  , boost::proto::multiplies< ElementVectorGrammar, ElementVectorGrammar >
	  , boost::proto::divides< ElementVectorGrammar, ElementVectorGrammar >
	  , boost::proto::negate< ElementVectorGrammar >
> {};

struct ElementMatrixGrammar
	: boost::proto::or_<
	    boost::proto::plus< ElementMatrixGrammar, ElementMatrixGrammar >
	  , boost::proto::minus< ElementMatrixGrammar, ElementMatrixGrammar >
	  , boost::proto::multiplies< ElemetMatrixType, ElementVectorGrammar >
> {};


struct ElementMVGrammar
	: boost::proto::or_<
	  	 boost::proto::or_<
	    	 VectorType
	     , ElementVectorType
	     , boost::proto::terminal< boost::proto::convertible_to<real> >
		 >
	  , boost::proto::plus< ElementMVGrammar, ElementMVGrammar >
	  , boost::proto::minus< ElementMVGrammar, ElementMVGrammar >
	  , boost::proto::multiplies< ElementVectorGrammar, ElementVectorGrammar >
	  , boost::proto::multiplies< ElemetMatrixType, VectorGrammar >
	  , boost::proto::divides< ElementVectorGrammar, ElementVectorGrammar >
	  , boost::proto::negate< ElementMVGrammar >
> {};


struct MVGrammar
	: boost::proto::or_<
		 boost::proto::or_<
			 VectorType
		 , boost::proto::terminal< boost::proto::convertible_to<real> >
//		 , ElementVectorType
		 >
	  , boost::proto::plus< MVGrammar, MVGrammar >
	  , boost::proto::minus< MVGrammar, MVGrammar >
	  , boost::proto::multiplies< MVGrammar, MVGrammar >
	  , boost::proto::multiplies< MatrixType, MVGrammar >
	  , boost::proto::divides< MVGrammar, MVGrammar >
	  , boost::proto::negate< MVGrammar >
> {};


// ******************************** Forward-Declared Expression Wrappers

template<typename Expr>
struct ElemMVExprWrapper;

template<typename Expr>
struct ElemVectorExprWrapper;

template<typename Expr>
struct VectorExprWrapper;

template<typename Expr>
struct MVExprWrapper;


// ******************************** Proto-fying Bare Data Types

typedef VectorExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareVector > > > ProtoVector;
typedef VectorExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareVectorView > > > ProtoVectorView;
typedef VectorExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareDistributedVector > > > ProtoDistributedVector;
typedef VectorExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareDistributedVectorView > > > ProtoDistributedVectorView;

typedef ElemVectorExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareElementVector > > > ProtoElementVector;
typedef ElemVectorExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareElementVectorView > > > ProtoElementVectorView;

typedef ElemMVExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareElementMatrix > > > ProtoElementMatrix;
typedef ElemMVExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareLightweightScaledElementMatrix > > > ProtoLightweightScaledElementMatrix;
typedef ElemMVExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareElementMatrixView > > > ProtoElementMatrixView;

typedef MVExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareCSRMatrix > > > ProtoCSRMatrix;
typedef MVExprWrapper<boost::proto::basic_expr<boost::proto::tag::terminal, boost::proto::term< BareMatrix > > > ProtoMatrix;


// ******************************** Proto Domains

struct MVDomain : boost::proto::domain< boost::proto::generator<MVExprWrapper> , MVGrammar > {};
struct ElementMVDomain : boost::proto::domain< boost::proto::generator<ElemMVExprWrapper> , ElementMVGrammar, MVDomain > {};
struct ElementVectorDomain : boost::proto::domain< boost::proto::generator<ElemVectorExprWrapper> , ElementVectorGrammar, ElementMVDomain > {};
struct VectorDomain : boost::proto::domain< boost::proto::generator<VectorExprWrapper> , VectorGrammar, ElementVectorDomain > {};


// ******************************** Proto Expression Wrappers I

template<typename Expr>
struct ElemMVExprWrapper : boost::proto::extends<Expr, ElemMVExprWrapper<Expr>, ElementMVDomain>
{
	typedef boost::proto::extends<Expr, ElemMVExprWrapper<Expr>, ElementMVDomain>   base_type;

	ElemMVExprWrapper(Expr const& expr = Expr()) : base_type(expr) {}
};

template<typename Expr>
struct MVExprWrapper : boost::proto::extends<Expr, MVExprWrapper<Expr>, MVDomain>
{
	typedef boost::proto::extends<Expr, MVExprWrapper<Expr>, MVDomain>   base_type;

	MVExprWrapper(Expr const& expr = Expr()) : base_type(expr) {}
};



// *************************** Matrix-Proto-Types

#include "proto_matrix.hpp"

typedef CSRMatrix_<> CSRMatrix;
typedef Matrix_<> Matrix;

typedef ElementMatrix_<> ElementMatrix;
typedef LightweightScaledElementMatrix_<> LightweightScaledElementMatrix;
typedef ElementMatrixView_<> ElementMatrixView;




// *************************** Proto Execution Contexts

#include "proto_context.hpp"





// *************************** Proto Expression Wrappers II

template<typename Expr>
struct ElemVectorExprWrapper : boost::proto::extends<Expr, ElemVectorExprWrapper<Expr>, ElementVectorDomain>
{
	typedef boost::proto::extends<Expr, ElemVectorExprWrapper<Expr>, ElementVectorDomain>   base_type;

	ElemVectorExprWrapper(Expr const& expr = Expr()) : base_type(expr) {}
};


template<typename Expr>
struct VectorExprWrapper : boost::proto::extends<Expr, VectorExprWrapper<Expr>, VectorDomain>
{
	typedef boost::proto::extends<Expr, VectorExprWrapper<Expr>, VectorDomain>   base_type;

	VectorExprWrapper(Expr const& expr = Expr()) : base_type(expr) {}

	typedef real result_type;

	real operator[](index_t i) const
	{
		VectorSubscriptContext const ctx(i);
		return boost::proto::eval(*this, ctx);
	}
};



// ************************************** Vector-Proto-Types

#include "proto_vector.hpp"




// ************************************** Proto Element Array Types

#include "element_arraytypes.hpp"





#endif /* ARITHMETIC_HPP_ */
