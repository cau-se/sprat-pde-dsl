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
#include <cfloat>

/*
 * This example solves the 2D advection equation
 *
 * d/dt u(t,x,y) + d/dx (q_x(t,x,y,u) * u) + d/dy (q_y(t,x,y,u) * u) = 0
 *
 * on the domain [0,1]^2 using a flux-corrected FEM scheme with P1 elements.
 * An explicit Runge-Kutta method is employed for time stepping.
 *
 * For more information about the algorithm see "example2.pdf". Note, that
 * in the referenced document the algorithm is presented with periodic
 * boundary conditions while here we introduce specific boundary operators
 * to handle the in-/outflow.
 */

#define REAL_FIND_MIN_INIT			DBL_MAX
#define REAL_FIND_MAX_INIT			-1.0

#define N_DIMENSIONS   2

#include "../pdedsl/sprat_pde_dsl_base.hpp"
typedef FEMMeshRectP1<N_DIMENSIONS> FEMMeshT;
#include "../pdedsl/sprat_pde_dsl_arithmetic.hpp"


void print2DVector(Vector v, uint n_x, uint n_y) {
	for(uint row=0; row<n_y; ++row) {
		for(uint col=0; col<n_x; ++col) {
			printf("%12.6f", v[col*n_y + row]);
		}
		printf("\n");
	}
	printf("\n");
}

int main(int argc, char ** argv) {
	RectMeshDimension dimensions[N_DIMENSIONS];
	dimensions[0] = RectMeshDimension("x", Interval(0.0, 1.0), 12);
	dimensions[1] = RectMeshDimension("y", Interval(0.0, 1.0), 12);
	FEMMeshT femMesh(dimensions);

	// Solution
	Vector u(femMesh.nDoF());

	// Set initial conditions
	foreach_omp(auto i, DoF(femMesh), ,{
		const real x = i.positionInDimension(0);
		const real y = i.positionInDimension(1);
		u[i] = (x<0.5 ? 1.0 : 0.0);
	})
	print2DVector(u, femMesh.dimensions[0].nNodes, femMesh.dimensions[1].nNodes);


	// Flux vectors
	Vector q_x(femMesh.nDoF());
	Vector q_y(femMesh.nDoF());

	// Vectors for flux corrected transport (FCT) solver
	Vector Delta_u_L(femMesh.nDoF());
	Vector Delta_u_H(femMesh.nDoF());
	Vector P_plus(femMesh.nDoF());
	Vector P_minus(femMesh.nDoF());
	Vector u_L_max(femMesh.nDoF());
	Vector u_L_min(femMesh.nDoF());

	// Vectors for cg solver
	Vector b(femMesh.nDoF());
	Vector p(femMesh.nDoF());
	Vector a(femMesh.nDoF());
	Vector r(femMesh.nDoF());


	// Element vectors for FCT solver
	ElementVectorArray F_A(femMesh);
	ElementVectorArray F_H_0(femMesh);
	ElementVectorArray F_H_1(femMesh);
	ElementVectorArray F_H_2(femMesh);
	ElementVectorArray F_H(femMesh);
	ElementVectorArray F_L(femMesh);


	// Mass matrix M_C
	class : public LeightweightElementMatrixArray {
	public:
		LightweightScaledElementMatrix operator[](ElementT const& element) const {
			return LightweightScaledElementMatrix(element, element.diamInDimension(0) * element.diamInDimension(1), _referenceValues);
		}
	} M_C;

	// Lumped mass matrix
	Vector M_L_assembled(femMesh.nDoF());

	// Internal transport operators Cin_x/y
	class : public LeightweightElementMatrixArray {
	public:
		LightweightScaledElementMatrix operator[](ElementT const& element) const {
			return LightweightScaledElementMatrix(element, element.diamInDimension(1), _referenceValues);
		}
	} Cin_x;
	class : public LeightweightElementMatrixArray {
	public:
		LightweightScaledElementMatrix operator[](ElementT const& element) const {
			return LightweightScaledElementMatrix(element, element.diamInDimension(0), _referenceValues);
		}
	} Cin_y;

	// Border part of transport operators Cbo_x/y
	BorderElementMatrixArray Cbo_x(femMesh);
	BorderElementMatrixArray Cbo_y(femMesh);


	// Reference element matrices
	ElementVector M_L_reference(ElementT(femMesh, 0));
	ElementMatrix M_C_reference(ElementT(femMesh, 0));
	ElementMatrix M_C_sp_reference(ElementT(femMesh, 0));
	ElementMatrix Cin_x_reference(ElementT(femMesh, 0));
	ElementMatrix Cin_y_reference(ElementT(femMesh, 0));
	ElementMatrix Cin_r_reference(ElementT(femMesh, 0));

	// Compute reference element matrices
	foreach(auto i, ElementDoFIndices(femMesh), {
		M_L_reference[i] = 0.0;
		foreach(auto j, ElementDoFIndices(femMesh), {
			M_C_reference(i, j) = femMesh.integrateReferenceElement(i, j);
			M_L_reference[i] += M_C_reference(i, j);
			Cin_x_reference(i, j) = femMesh.integratePartIntDerivativeReferenceElement(i, j, 0);
			Cin_y_reference(i, j) = femMesh.integratePartIntDerivativeReferenceElement(i, j, 1);
		})
	})

	M_C.reset(femMesh, M_C_reference);
	Cin_x.reset(femMesh, Cin_x_reference);
	Cin_y.reset(femMesh, Cin_y_reference);

	// Lump mass matrix
	M_L_assembled = 0.0;
	foreachElementIndependently(auto tau, femMesh, , {
		foreach(auto i_local, ElementDoF(tau), {
			M_L_assembled[i_local.globalIndex()] += tau.diamInDimension(0) * tau.diamInDimension(1) * M_L_reference[i_local];
		})
	})

	real volOmega = 0.0;
	foreach_omp(auto i, DoF(femMesh), reduction(+:volOmega), {
		volOmega += M_L_assembled[i];
	})

	// Compute border part of transport operators
	foreach_omp(auto tau, Elements(femMesh), , {
		foreach(auto edge, DomainBoundaryHypersurfaces(tau), {
			foreach(auto i, HypersurfaceDoF(edge), {
				foreach(auto j, HypersurfaceDoF(edge), {
					Cbo_x[tau].add(i.elementDoFIndex(), j.elementDoFIndex(), -1.0*edge.surfaceIntegral(i, j, 0));
					Cbo_y[tau].add(i.elementDoFIndex(), j.elementDoFIndex(), -1.0*edge.surfaceIntegral(i, j, 1));
				})
			})
		})
	})




	const real Delta_t = 0.1;
	real t = 0.0;
	uint iteration = 0;
	while(t < 1.0) {
		q_x = 0.1 * u;
		q_y = 0.0;
		// Prescribe Neumann boundary conditions on inflow edge
		foreach_omp(auto i, DoF(femMesh), shared(q_x), {
			if(i.indexInDimension(0) == 0) {
				q_x[i] = 0.0;
			}
		})

		// Preparation
		F_H_1.swapContent(F_H_2); // (0 2 1)
		F_H_0.swapContent(F_H_1); // (2 0 1)

		// Step 1.
		Delta_u_L = b = 0.0;

		ElementMatrix D(ElementT(femMesh, 0));
		foreachElementIndependently(auto tau, femMesh, private(D), {
			// a) Update D^\tau
			foreach(auto i, ElementDoF(tau), {
				const auto i_global = i.globalIndex();
				D(i, i) = 0.0;
				foreach(auto j, ElementDoF(tau), {
					if(i!=j) {
						const auto j_global = j.globalIndex();
						const real u_i = u[i_global];
						const real u_j = u[j_global];
						real k_ij = 0.0;
						real k_ji = 0.0;

						if(fabs(u_j - u_i) > 1e-14) {
							const real inverse_u_diff = 1.0 / (u_j - u_i);
							const real v_x_ij = (q_x[j_global] - q_x[i_global]) * inverse_u_diff; // == v_x_ji
							const real v_y_ij = (q_y[j_global] - q_y[i_global]) * inverse_u_diff;

							k_ij = (Cin_x[tau](j, i) + Cbo_x[tau](j, i)) * v_x_ij
								 + (Cin_y[tau](j, i) + Cbo_y[tau](j, i)) * v_y_ij;
							k_ji = (Cin_x[tau](i, j) + Cbo_x[tau](i, j)) * v_x_ij
								 + (Cin_y[tau](i, j) + Cbo_y[tau](i, j)) * v_y_ij;
//							k_ij = Cin_x[tau](j, i) * v_x_ij
//								 + Cin_y[tau](j, i) * v_y_ij;
//							k_ji = Cin_x[tau](i, j) * v_x_ij
//								 + Cin_y[tau](i, j) * v_y_ij;
						}

						D(i, j) = max(0.0, -1.0*k_ij, -1.0*k_ji);
						D(i, i) -= D(i, j);
					}
				})
			})

			// b) Calculate F^H_\tau and F^L_\tau
			F_H_0[tau] = Cin_x[tau]*q_x + Cbo_x[tau]*q_x
			           + Cin_y[tau]*q_y + Cbo_y[tau]*q_y;
//			F_H_0[tau] = Cin_x[tau]*q_x
//					   + Cin_y[tau]*q_y;

			//D.setElement(tau);
			F_L[tau]   = F_H_0[tau] + D*u;

			// Runge-Kutta method
			if(iteration>=2) {
				const real alpha30 = 23.0/12.0;
				const real alpha31 = -4.0/3.0;
				const real alpha32 = 5.0/12.0;
				F_H[tau] = alpha30*F_H_0[tau]
				         + alpha31*F_H_1[tau]
				         + alpha32*F_H_2[tau];
			}
			else if(iteration==1) {
				const real alpha20 = 1.5;
				const real alpha21 = -0.5;
				F_H[tau] = alpha20*F_H_0[tau]
				         + alpha21*F_H_1[tau];
			}
			else if(iteration==0) {
				F_H[tau] = F_H_0[tau];
			}

			// c)
			b         += F_H[tau];
			Delta_u_L += F_L[tau];
		})


		// Step 2.
		Delta_u_L = Delta_t * Delta_u_L/M_L_assembled;


		//Step 3.
		b *= Delta_t;
		Delta_u_H = Delta_u_L;


		// CG
		// r <- b - Ax
		r = b;
		foreachElementIndependently(auto tau, femMesh, , {
			r -= M_C[tau] * Delta_u_H;
		})

		p = r;

		const uint cgIterations = 8;
		for(uint i=0; i<cgIterations && p.norm()>0.0; ++i) {
			a = 0.0;
			foreachElementIndependently(auto tau, femMesh, , {
				a += M_C[tau] * p;
			})

			const real p_dot_a = p.dotProduct(a);
			const real lambda_opt = p.dotProduct(r) / p_dot_a;
			Delta_u_H += lambda_opt * p;
			r -= lambda_opt * a;

			const real mu = r.dotProduct(a) / p_dot_a;
			p = r - mu * p;
		}


		// Step 4. (for the time being, choose u to be u_L)
		u += Delta_u_L;
		VectorView u_L(u);

		// Step 5.
		P_minus = P_plus = 0.0;
		u_L_max = REAL_FIND_MAX_INIT;
		u_L_min = REAL_FIND_MIN_INIT;
		foreachElementIndependently(auto tau, femMesh, , {
			// a)
			F_A[tau] = tau.diamInDimension(0) * tau.diamInDimension(1) * M_L_reference * Delta_u_H
			         - M_C[tau] * Delta_u_H
			         + Delta_t * (F_H[tau] - F_L[tau]);
			// b)
			real u_L_element_max = REAL_FIND_MAX_INIT;
			real u_L_element_min = REAL_FIND_MIN_INIT;
			foreach(auto i, ElementDoF(tau), {
				const auto i_global = i.globalIndex();
				// i.
				u_L_element_max = max(u_L_element_max, u_L[i_global]);
				u_L_element_min = min(u_L_element_min, u_L[i_global]);

				// ii.
				P_plus[i_global]  += max(0.0, F_A[tau][i]);
				P_minus[i_global] += min(0.0, F_A[tau][i]);
			})

			// c)
			foreach(auto i, ElementDoF(tau), {
				const auto i_global = i.globalIndex();
				u_L_max[i_global] = max(u_L_max[i_global], u_L_element_max);
				u_L_min[i_global] = min(u_L_min[i_global], u_L_element_min);
			})
		})

		// Step 6.
		foreach_omp(auto i, DoF(femMesh), , {
			if(P_plus[i] > 0.0) {
				P_plus[i] = min(1.0, (u_L_max[i] - u_L[i]) / P_plus[i]);
			}
			if(P_minus[i] < 0.0) {
				P_minus[i] = min(1.0, (u_L_min[i] - u_L[i]) / P_minus[i]);
			}
		})

		// Step 7.
		// nop; (It already holds: u == u_L)

		// Step 8.
		VectorView Delta_u_A(Delta_u_H);
		Delta_u_A = 0.0;
		foreachElementIndependently(auto tau, femMesh, , {
			// a)
			real lambda_tau = 1.0;
			foreach(auto i, ElementDoF(tau), {
				if(F_A[tau][i] > 0.0) {
					lambda_tau = min(lambda_tau, P_plus[i.globalIndex()]);
				}
				else if(F_A[tau][i] < 0.0) {
					lambda_tau = min(lambda_tau, P_minus[i.globalIndex()]);
				}
			})

			// b)
			Delta_u_A += lambda_tau * F_A[tau];
		})

		// Step 9.
		u += Delta_u_A / M_L_assembled;

		// Step 10: PostProcessLowOrderSolution
		index_t positive_nodes = 0;
		real mass_to_redistribute = 0.0;
		real positive_mass = 0.0;
		foreach_omp(auto i, DoF(femMesh), reduction(+:mass_to_redistribute,positive_mass,positive_nodes), {
			if(u[i]<=0.0) {
				mass_to_redistribute += u[i] * M_L_assembled[i];
				u[i] = 0.0;
			}
			else {
				positive_mass += u[i] * M_L_assembled[i];
				++positive_nodes;
			}
		})

		if(positive_mass >= -1.0*mass_to_redistribute) {
			foreach_omp(auto i, LocalDoF(femMesh), , {
				if(u[i]>0.0) {
					u[i] += (u[i]/positive_mass) * mass_to_redistribute;
				}
			})
		}
		else {
			foreach_omp(auto i, LocalDoF(femMesh), , {
				u[i] += (1.0/volOmega) * mass_to_redistribute;
			})
		}

		iteration += 1;
		t += Delta_t;

		print2DVector(u, femMesh.dimensions[0].nNodes, femMesh.dimensions[1].nNodes);
	}

	return 0;
}
