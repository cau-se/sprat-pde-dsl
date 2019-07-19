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
 * This example demonstrates how to use the Sprat PDE Solver DSL for
 * parallel distributed computations with MPI+OpenMP.
 *
 * To run the example use
 * mpirun -np <number of processes> ./example3
 */

#define N_DIMENSIONS   2

#include "../pdedsl/sprat_pde_dsl_base.hpp"
typedef FEMMeshRectP1<N_DIMENSIONS> FEMMeshT;
#include "../pdedsl/sprat_pde_dsl_arithmetic.hpp"


void print2DVector(Vector v, uint n_x, uint n_y) {
	for(uint row=0; row<n_y; ++row) {
		for(uint col=0; col<n_x; ++col) {
			printf("%14.6e", v[col*n_y + row]);
		}
		printf("\n");
	}
	printf("\n");
}

int main(int argc, char ** argv) {
	ParallelExecutionEnvironment pEE;
	pEE.init(&argc, &argv);

	// Register vector name(s) for data exchange
	pEE.registerName("u");

	const uint maxIterations = 10;

	if(pEE.isMaster()) {
		RectMeshDimension dimensions[N_DIMENSIONS];
		dimensions[0] = RectMeshDimension("x", Interval(0.0, 1.0), 12);
		dimensions[1] = RectMeshDimension("y", Interval(0.0, 1.0), 12);
		FEMMeshT femMesh(dimensions);

		// Slice mesh orthogonal to the first dimension in as many strides as there are compute slaves
		std::vector<index_t> dofStartIndexForProcess;
		femMesh.createAndDistributeStrideDecomposition(pEE, dofStartIndexForProcess);

		// Set and distribute initial data to the slaves
		Vector u(femMesh.nDoF());
		foreach_omp(auto i, DoF(femMesh), shared(u), {
			u[i] = (real)i;
		})
		auto requests = u.distributeToSlavesAsync("u", pEE, dofStartIndexForProcess);
		pEE.waitFor(requests);

		// Do something for every compute iteration (e.g., record something every n-th frame)
		for(uint iteration=0; iteration<maxIterations; ++iteration) {
			// ...
		}

		// Gather and display results
		requests = u.collectFromSlavesAsync("u", pEE, dofStartIndexForProcess);
		pEE.waitFor(requests);

		print2DVector(u, femMesh.dimensions[0].nNodes, femMesh.dimensions[1].nNodes);
	}
	else { // As a compute slave..
		FEMMeshT femMesh;
		MPCommunicationRegistry comRegistry;
		// Receive mesh stride. The mesh contains local degrees of freedom (DoF)
		// as well as ghost DoF (DoF that do not belong to this compute slave but
		// share an element with local DoF).
		femMesh.receiveStrideDecomposition(pEE, comRegistry);


		DistributedVector u("u", pEE, comRegistry);
		u.collectFromMasterAsync();
		u.finishPendingRecvs();

		Vector v(u.size());

		for(uint iteration=0; iteration<maxIterations; ++iteration) {
			// Each data exchange operation updates the values for the ghost DoF
			u.exchangeDataBlocked();

			v = u;

			foreachElementIndependently(auto tau, femMesh, shared(v,u), {
				real elementSum = 0.0;
				foreach(auto i, ElementDoF(tau), {
					elementSum += u[i.globalIndex()];
				})
				foreach(auto i, ElementDoF(tau), {
					v[i.globalIndex()] += elementSum;
				})
			})

			u = v;
		}

		u.distributeToMasterAsync();
		u.finishPendingSends();
	}

	pEE.barrier();

	return 0;
}
