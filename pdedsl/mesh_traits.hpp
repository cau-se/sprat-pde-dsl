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

#ifndef MESH_TRAITS_HPP_
#define MESH_TRAITS_HPP_

#include "mesh.hpp"

// Defines types for selected FEMMesh type

typedef FEMMeshT::ElementT ElementT;
typedef FEMMeshT::DoFT DoFT;

static constexpr uint nDoFPerElement() {
	return FEMMeshT::nDoFPerElement();
}


#endif /* MESH_TRAITS_HPP_ */
