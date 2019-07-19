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

#include <iostream>

/*
 * This example solves the heat transport equation in 2D
 *
 * d/dt u(t,x,y) = c * (d^2/dx^2 u(t,x,y) + d^2/dy^2 u(t,x,y))
 *               = c * LAPLACE * u(t,x,y)
 *
 * on the domain [0,1]^2 with periodic boundary conditions using
 * the implicit Euler method. The solution of the resulting
 * linear system is approximated with the CG method. For
 * discretization in space we employ P1 finite elements on a
 * rectangular 5x5 grid.
 *
 * Applying the implicit Euler method:
 *     u*(t_{n+1}) = u*(t_{n}) + Delta_t * c * LAPLACE * u*(t_{n+1})
 * <=> (I - LAPLACE * Delta_t * c) u*(t_{n+1}) = u*(t_{n})
 *
 * The degrees of freedom are numbered in the following way
 * (wrapped edges!):
 *
 *    0   5  10  15  20   0
 *
 *    4   9  14  19  24   4
 *
 *    3   8  13  18  23   3
 *
 *    2   7  12  17  22   2
 * ^
 * |  1   6  11  16  21   1
 * y
 * |  0   5  10  15  20   0
 *   --x-->
 */

#define N_DIMENSIONS   2

#include "../pdedsl/sprat_pde_dsl_base.hpp"
typedef FEMMeshRectP1Periodic<N_DIMENSIONS> FEMMeshT;
#include "../pdedsl/sprat_pde_dsl_arithmetic.hpp"

int main(int argc, char ** argv) {
	RectMeshDimension dimensions[N_DIMENSIONS];
	dimensions[0] = RectMeshDimension("x", Interval(0.0, 1.0), 5);
	dimensions[1] = RectMeshDimension("y", Interval(0.0, 1.0), 5);
	FEMMeshT femMesh(dimensions);

	Vector u(femMesh.nDoF());

	// Initial conditions
	foreach_omp(auto i, DoF(femMesh), ,{
		const real x = i.positionInDimension(0);
		const real y = i.positionInDimension(1);
		u[i] = (x<0.5 ? 1.0 : 0.0);
	})
	u.print();

	// Assemble Laplace operator in linked list matrix
	LILMatrix laplaceLIL(femMesh.nDoF());
	foreach(auto tau, Elements(femMesh), {
		foreach(auto i, ElementDoF(tau), {
			foreach(auto j, ElementDoF(tau), {
				laplaceLIL.addElement(i.globalIndex(), j.globalIndex(),
					tau.integratePartIntLaplace(i, j));
			})
		})
	})
	laplaceLIL.print();
	CSRMatrix laplace(laplaceLIL);
	laplace.print();

	const real delta_t = 0.1;
	const real c = 0.1;
	real t = 0.0;

	while(t < 1.0) {
		Vector p(femMesh.nDoF());
		Vector r(femMesh.nDoF());
		Vector a(femMesh.nDoF());

		r = u - (u - laplace * (delta_t * c * u));
		p = r;

		uint cgIteration = 0;
		while(r.norm() > 1.0e-10 * u.norm()  &&  p.norm() > 0.0) {
			a = p - laplace * (delta_t * c * p);

			const real p_dot_a = p.dotProduct(a);
			const real lambda_opt = p.dotProduct(r) / p_dot_a;

			u += lambda_opt * p;
			r -= lambda_opt * a;

			const real mu = r.dotProduct(a) / p_dot_a;
			p = r - mu * p;
			cgIteration++;
		}
		std::cout << "cg iterations: " << cgIteration << std::endl;

		t += delta_t;
		u.print();
	}

	return 0;
}
